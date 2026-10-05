"""Assemble the viewer runtime from its ES-module sources.

The viewer is authored as ordinary ES modules under ``web/src`` so editors,
linters and ``node --test`` can load them directly. A figure, however, must be
one self-contained HTML file that also works from ``file://``, where module
imports between inline scripts are unavailable. This module performs the small
amount of bundling that requires:

* modules are ordered topologically from the entry point;
* ``import { a, b } from "./x.js";`` lines are removed;
* ``export`` keywords on top-level declarations are dropped;
* everything is wrapped in one strict IIFE.

The rules are deliberately narrow so the transform stays obviously correct:
only named relative imports without aliases are allowed, every imported name
must be exported by its module, every exported top-level name must be unique
across the whole runtime, and ``export default`` / re-exports are rejected.
Lines inside template literals, strings and comments are left alone.
Violations raise ``ViewerBuildError`` at build time rather than producing a
subtly broken figure.
"""

from __future__ import annotations

import base64
import bisect
import functools
import re
from pathlib import Path

WEB_ROOT = Path(__file__).resolve().parent / "web"
SRC_ROOT = WEB_ROOT / "src"
STYLE_ROOT = WEB_ROOT / "styles"
SKY_THUMB_ROOT = WEB_ROOT / "assets" / "sky"
ENTRY_MODULE = "app/main.js"

_IMPORT_RE = re.compile(
    r"^import\s*\{(?P<names>[^}]*)\}\s*from\s*[\"'](?P<path>\.{1,2}/[^\"']+)[\"'];?[ \t]*$",
    re.M,
)
_ANY_IMPORT_RE = re.compile(r"^\s*import\b(?!\s*\()", re.M)
_EXPORT_DECL_RE = re.compile(
    r"^export\s+(?=(?:async\s+)?function\b|class\b|const\b|let\b)", re.M
)
_ANY_EXPORT_RE = re.compile(r"^export\b", re.M)
_BINDING_RE = re.compile(r"[A-Za-z_$][\w$]*")


class ViewerBuildError(RuntimeError):
    """Raised when the viewer sources violate the bundling rules."""


# Lexical scan: where do strings, template literals and comments span lines?
# The transforms below are line-anchored, so this is all they need to know to
# leave such text alone (a code sample in a help string, a shader, ...).
_CODE_STOP_RE = re.compile(r"[`'\"/]")
_NESTED_CODE_STOP_RE = re.compile(r"[`'\"/{}]")  # code inside a template's ${...}
_TEMPLATE_STOP_RE = re.compile(r"`|\\|\$\{")
_QUOTE_STOP_RE = {q: re.compile(r"\\|\n|" + q) for q in "'\""}
_REGEX_AFTER = frozenset("(,=:[!&|?{};+-*%<>~^}")
_REGEX_AFTER_WORDS = frozenset(
    "return typeof instanceof in of new delete void throw case do else yield await".split()
)


def _regex_allowed(source: str, at: int) -> bool:
    """Whether a ``/`` at ``at`` starts a regular expression (not a division)."""

    k = at - 1
    while k >= 0 and source[k] in " \t\r\n":
        k -= 1
    if k < 0 or source[k] in _REGEX_AFTER:
        return True
    end = k + 1
    while k >= 0 and (source[k].isalnum() or source[k] in "_$"):
        k -= 1
    return source[k + 1:end] in _REGEX_AFTER_WORDS


def _regex_end(source: str, i: int) -> int | None:
    """End of a regular-expression literal whose body starts at ``i``."""

    in_class = False
    while i < len(source):
        c = source[i]
        if c == "\\":
            i += 2
            continue
        if c == "\n":
            return None  # not a regex after all: a division
        if in_class:
            in_class = c != "]"
        elif c == "[":
            in_class = True
        elif c == "/":
            return i + 1
        i += 1
    return None


class _Text:
    """Top-level spans of a module's source that are text, not code."""

    def __init__(self, rel: str, source: str):
        spans: list[tuple[int, int]] = []
        # One entry per open template literal: "T" in its text, or the brace
        # depth of the code inside one of its ${...}.
        stack: list = []
        start = i = 0
        n = len(source)

        def unterminated(kind: str, at: int) -> ViewerBuildError:
            line = source.count("\n", 0, at) + 1
            return ViewerBuildError(
                f"{rel}: unterminated {kind} at line {line}; the bundler cannot tell code from text here"
            )

        while i < n:
            if stack and stack[-1] == "T":
                m = _TEMPLATE_STOP_RE.search(source, i)
                if m is None:
                    raise unterminated("template literal", start)
                i = m.end()
                if m.group() == "\\":
                    i += 1
                elif m.group() == "`":
                    stack.pop()
                    if not stack:
                        spans.append((start, i))
                else:
                    stack.append(0)
                continue
            m = (_NESTED_CODE_STOP_RE if stack else _CODE_STOP_RE).search(source, i)
            if m is None:
                break
            j, ch = m.start(), m.group()
            i = j + 1
            if ch == "{":
                stack[-1] += 1
            elif ch == "}":
                if stack[-1]:
                    stack[-1] -= 1
                else:
                    stack.pop()  # back in the template's text
            elif ch == "`":
                if not stack:
                    start = j
                stack.append("T")
            elif ch != "/":  # a quoted string
                while True:
                    q = _QUOTE_STOP_RE[ch].search(source, i)
                    if q is None or q.group() == "\n":
                        raise unterminated("string", j)
                    i = q.end() + (1 if q.group() == "\\" else 0)
                    if q.group() == ch:
                        break
                if not stack and source.find("\n", j, i) >= 0:
                    spans.append((j, i))
            elif source.startswith("/", i):
                k = source.find("\n", i)
                i = n if k < 0 else k
            elif source.startswith("*", i):
                k = source.find("*/", i + 1)
                if k < 0:
                    raise unterminated("comment", j)
                i = k + 2
                if not stack and source.find("\n", j, i) >= 0:
                    spans.append((j, i))
            elif _regex_allowed(source, j):
                i = _regex_end(source, i) or i
        if stack:
            raise unterminated("template literal", start)
        self.spans = spans
        self._starts = [s for s, _ in spans]

    def is_code(self, pos: int) -> bool:
        k = bisect.bisect_right(self._starts, pos) - 1
        return k < 0 or not (self.spans[k][0] < pos < self.spans[k][1])


@functools.lru_cache(maxsize=256)
def _text(rel: str, source: str) -> _Text:
    return _Text(rel, source)


def _imports(rel: str, source: str) -> list[tuple[re.Match, str, list[str]]]:
    """``(match, module, names)`` for each named relative import in code."""

    text = _text(rel, source)
    out = []
    for m in _IMPORT_RE.finditer(source):
        if not text.is_code(m.start()):
            continue
        names = []
        for raw in m.group("names").split(","):
            name = raw.strip()
            if not name:
                continue
            if not _BINDING_RE.fullmatch(name):
                raise ViewerBuildError(
                    f"{rel}: unsupported import binding `{' '.join(name.split())}` from {m.group('path')}; "
                    "aliases (`a as b`) are not supported, import the exported name itself"
                )
            names.append(name)
        out.append((m, _resolve(rel, m.group("path")), names))
    return out


def _read_module(rel_path: str) -> str:
    path = (SRC_ROOT / rel_path).resolve()
    if not path.is_file():
        raise ViewerBuildError(f"Viewer module not found: {rel_path}")
    return path.read_text(encoding="utf-8")


def _resolve(from_rel: str, spec: str) -> str:
    base = (SRC_ROOT / from_rel).parent
    resolved = (base / spec).resolve()
    try:
        return resolved.relative_to(SRC_ROOT.resolve()).as_posix()
    except ValueError as exc:  # pragma: no cover - defensive
        raise ViewerBuildError(f"Import escapes the source tree: {spec} in {from_rel}") from exc


def module_order(entry: str = ENTRY_MODULE) -> list[str]:
    """Return modules reachable from ``entry`` in dependency order."""

    order: list[str] = []
    state: dict[str, str] = {}

    def visit(rel: str, chain: tuple[str, ...]) -> None:
        mark = state.get(rel)
        if mark == "done":
            return
        if mark == "active":
            cycle = " -> ".join(chain + (rel,))
            raise ViewerBuildError(f"Circular viewer import: {cycle}")
        state[rel] = "active"
        source = _read_module(rel)
        for _match, target, _names in _imports(rel, source):
            visit(target, chain + (rel,))
        state[rel] = "done"
        order.append(rel)

    visit(entry, ())
    return order


def _strip_module(rel: str, source: str) -> str:
    text = _text(rel, source)
    for m in _ANY_EXPORT_RE.finditer(source):
        if text.is_code(m.start()) and not _EXPORT_DECL_RE.match(source, m.start()):
            raise ViewerBuildError(
                f"{rel}: only `export function|class|const|let` declarations are supported"
            )
    cuts = [m.span() for m, _target, _names in _imports(rel, source)]
    imported = {start for start, _end in cuts}
    for m in _ANY_IMPORT_RE.finditer(source):
        at = m.end() - len("import")
        if text.is_code(at) and at not in imported:
            end = source.find("\n", at)
            line = source[at: end if end >= 0 else len(source)]
            raise ViewerBuildError(f"{rel}: unsupported import form: {line.strip()}")
    cuts += [m.span() for m in _EXPORT_DECL_RE.finditer(source) if text.is_code(m.start())]
    body, at = [], 0
    for start, end in sorted(cuts):
        body.append(source[at:start])
        at = end
    body.append(source[at:])
    return "".join(body)


_EXPORT_NAME_RE = re.compile(
    r"^export\s+(?:async\s+)?(?:function\*?|class|const|let)\s+([A-Za-z_$][\w$]*)",
    re.M,
)


def _exported_names(rel: str, source: str) -> list[str]:
    text = _text(rel, source)
    return [m.group(1) for m in _EXPORT_NAME_RE.finditer(source) if text.is_code(m.start())]


@functools.lru_cache(maxsize=4)
def _bundle_cached(entry: str, fingerprint: tuple) -> str:
    del fingerprint  # only used as a cache key
    order = module_order(entry)
    owners: dict[str, str] = {}
    exports: dict[str, set[str]] = {}
    parts: list[str] = []
    for rel in order:
        source = _read_module(rel)
        exported = _exported_names(rel, source)
        # Dependencies come first in `order`, so their exports are known.
        for _match, target, names in _imports(rel, source):
            missing = [name for name in names if name not in exports.get(target, ())]
            if missing:
                listed = ", ".join(f"`{name}`" for name in missing)
                raise ViewerBuildError(
                    f"{rel} imports {listed} from {target}, which does not export "
                    + ("it" if len(missing) == 1 else "them")
                )
        exports[rel] = set(exported)
        for name in exported:
            other = owners.get(name)
            if other and other != rel:
                raise ViewerBuildError(
                    f"Exported name `{name}` is declared in both {other} and {rel}"
                )
            owners[name] = rel
        body = _strip_module(rel, source).strip()
        # Each module gets a private scope; only its exports are shared.
        if exported:
            names = ", ".join(exported)
            parts.append(
                f"// ---- {rel} ----\nconst {{ {names} }} = (() => {{\n{body}\nreturn {{ {names} }};\n}})();\n"
            )
        else:
            parts.append(f"// ---- {rel} ----\n(() => {{\n{body}\n}})();\n")
    js = '(() => {\n"use strict";\n' + "\n".join(parts) + "})();\n"
    # The runtime is inlined in a <script> element: these sequences would end
    # it early (or switch the HTML tokenizer into its escaped states).
    lowered = js.lower()
    for bad in ("</script", "<!--"):
        at = lowered.find(bad)
        if at >= 0:
            line = js.count("\n", 0, at) + 1
            raise ViewerBuildError(f"Bundled runtime contains {bad!r} (bundle line {line}); build it from parts")
    return js


def _fingerprint() -> tuple:
    return tuple(
        (p.as_posix(), p.stat().st_mtime_ns)
        for p in sorted(SRC_ROOT.rglob("*.js"))
    )


def bundle_js(entry: str = ENTRY_MODULE) -> str:
    """Return the complete runtime as one classic script."""

    return _bundle_cached(entry, _fingerprint())


def bundle_css() -> str:
    """Return the viewer stylesheet (files concatenated in name order)."""

    files = sorted(STYLE_ROOT.glob("*.css"))
    return "\n".join(f.read_text(encoding="utf-8").strip() + "\n" for f in files)


@functools.lru_cache(maxsize=2)
def _sky_thumbnails_cached(fingerprint: tuple) -> tuple:
    return tuple(
        (name, "data:image/jpeg;base64," + base64.b64encode((SKY_THUMB_ROOT / name).read_bytes()).decode("ascii"))
        for name, _ in fingerprint
    )


def sky_thumbnails() -> dict[str, str]:
    """All-sky thumbnails of the Sky background picker's surveys.

    Keys are file stems (``p_dss2_color``, as ``thumbKey`` in
    ``web/src/sky/catalog.js`` derives them); values are JPEG data URIs.
    Regenerate the files with ``scripts/make_sky_thumbnails.py``.
    """

    files = sorted(SKY_THUMB_ROOT.glob("*.jpg"))
    pairs = _sky_thumbnails_cached(tuple((p.name, p.stat().st_mtime_ns) for p in files))
    return {Path(name).stem: uri for name, uri in pairs}


def template_html() -> str:
    return (WEB_ROOT / "template.html").read_text(encoding="utf-8")
