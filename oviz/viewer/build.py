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
only named relative imports are allowed, every exported top-level name must be
unique across the whole runtime, and ``export default`` / re-exports are
rejected. Violations raise ``ViewerBuildError`` at build time rather than
producing a subtly broken figure.
"""

from __future__ import annotations

import functools
import re
from pathlib import Path

WEB_ROOT = Path(__file__).resolve().parent / "web"
SRC_ROOT = WEB_ROOT / "src"
STYLE_ROOT = WEB_ROOT / "styles"
ENTRY_MODULE = "app/main.js"

_IMPORT_RE = re.compile(
    r"^import\s*\{(?P<names>[^}]*)\}\s*from\s*[\"'](?P<path>\.{1,2}/[^\"']+)[\"'];?[ \t]*$",
    re.M,
)
_ANY_IMPORT_RE = re.compile(r"^\s*import\b(?!\s*\()", re.M)
_EXPORT_DECL_RE = re.compile(
    r"^export\s+(?=(?:async\s+)?function\b|class\b|const\b|let\b)", re.M
)
_BAD_EXPORT_RE = re.compile(r"^export\s+(?:default\b|\{|\*)", re.M)


class ViewerBuildError(RuntimeError):
    """Raised when the viewer sources violate the bundling rules."""


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
        for match in _IMPORT_RE.finditer(source):
            visit(_resolve(rel, match.group("path")), chain + (rel,))
        state[rel] = "done"
        order.append(rel)

    visit(entry, ())
    return order


def _strip_module(rel: str, source: str) -> str:
    if _BAD_EXPORT_RE.search(source):
        raise ViewerBuildError(
            f"{rel}: only `export function|class|const|let` declarations are supported"
        )
    body = _IMPORT_RE.sub("", source)
    leftover = _ANY_IMPORT_RE.search(body)
    if leftover:
        line = body[leftover.start(): body.find("\n", leftover.start())]
        raise ViewerBuildError(f"{rel}: unsupported import form: {line.strip()}")
    return _EXPORT_DECL_RE.sub("", body)


_EXPORT_NAME_RE = re.compile(
    r"^export\s+(?:async\s+)?(?:function\*?|class|const|let)\s+([A-Za-z_$][\w$]*)",
    re.M,
)


@functools.lru_cache(maxsize=4)
def _bundle_cached(entry: str, fingerprint: tuple) -> str:
    del fingerprint  # only used as a cache key
    order = module_order(entry)
    owners: dict[str, str] = {}
    parts: list[str] = []
    for rel in order:
        source = _read_module(rel)
        exported = _EXPORT_NAME_RE.findall(source)
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
    return '(() => {\n"use strict";\n' + "\n".join(parts) + "})();\n"


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


def template_html() -> str:
    return (WEB_ROOT / "template.html").read_text(encoding="utf-8")
