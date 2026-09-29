"""Standalone HTML figures for the Oviz viewer."""

from __future__ import annotations

import hashlib
import html
import json
import re
import tempfile
import webbrowser
from pathlib import Path
from typing import Any

import numpy as np

from .build import bundle_css, bundle_js, template_html
from .bundle import Bundle, encode_base64_lines
from .compile import compile_scene_spec

VIEWER_VERSION = "2.0.0"
_PLACEHOLDER_RE = re.compile(r"__(?:TITLE|VERSION|THEME|MODE|CSS|RUNTIME|MANIFEST|BLOBS)__")

#: Viewer modes (see ``web/src/ui/modes.js``): "focus" shows the figure, its
#: key and one quiet bar; "detailed" keeps every panel and control in view.
VIEWER_MODES = ("focus", "detailed")


def _check_mode(mode: str | None) -> str:
    mode = str(mode or "focus").strip().lower()
    if mode not in VIEWER_MODES:
        raise ValueError(f"viewer mode must be one of {VIEWER_MODES}, got {mode!r}")
    return mode


def _json_default(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        if value.dtype == object:
            return value.tolist()
        return [value.dtype.str, list(value.shape), hashlib.sha1(np.ascontiguousarray(value).tobytes()).hexdigest()]
    if isinstance(value, np.generic):
        return value.item()
    raise TypeError(f"cannot fingerprint {type(value).__name__}")


def _hash_json(digest: Any, value: Any, depth: int) -> None:
    """Feed ``value`` to ``digest`` as JSON, split into pieces ``depth`` levels deep."""

    if depth and isinstance(value, dict):
        digest.update(b"{")
        for k, v in value.items():
            digest.update(json.dumps(k, default=_json_default).encode() + b":")
            _hash_json(digest, v, depth - 1)
            digest.update(b",")
        digest.update(b"}")
    elif depth and isinstance(value, (list, tuple)):
        digest.update(b"[")
        for v in value:
            _hash_json(digest, v, depth - 1)
            digest.update(b",")
        digest.update(b"]")
    else:
        digest.update(json.dumps(value, default=_json_default, separators=(",", ":")).encode())


def _spec_fingerprint(spec: dict[str, Any]) -> str | None:
    """Content hash of the scene spec, or None when it cannot be hashed.

    Any in-place edit (a frame's colour, a point, a volume setting, States
    attached after construction) must trigger a recompile, so the whole
    spec is hashed. Its JSON is fed in pieces (per frame, volume layer,
    member cluster, ...) so a spec of hundreds of MB is never one string;
    hashing a 413 MB spec takes ~3 s against ~13 s to compile it.
    """

    digest = hashlib.sha1()
    try:
        _hash_json(digest, spec, 3)
    except (TypeError, ValueError, RecursionError):
        return None
    return digest.hexdigest()


def _json_for_script(obj: Any) -> str:
    text = json.dumps(obj, separators=(",", ":"), ensure_ascii=False, allow_nan=False)
    # Never let data terminate the surrounding <script> element; \u003c is a
    # valid JSON escape (unlike "<\\!--"), so JSON.parse still accepts it.
    return text.replace("<", "\\u003c")


def render_bundle_html(bundle: Bundle, *, title: str | None = None, theme: str = "dark",
                       mode: str = "focus") -> str:
    mode = _check_mode(mode)
    manifest = dict(bundle.manifest)
    manifest["blobs"] = [b.descriptor() for b in bundle.blobs]
    manifest["viewer"] = {"version": VIEWER_VERSION, "mode": mode}
    blob_html = []
    for b in sorted(bundle.blobs, key=lambda b: (b.priority, b.id)):
        attrs = (
            f'type="application/octet-stream" data-oviz-blob="{html.escape(b.id, quote=True)}"'
        )
        blob_html.append(f"<script {attrs}>\n{encode_base64_lines(b.payload)}\n</script>")
    page_title = title if title is not None else (manifest.get("title") or "Oviz figure")
    replacements = {
        "__TITLE__": html.escape(str(page_title)),
        "__VERSION__": VIEWER_VERSION,
        "__THEME__": html.escape(theme, quote=True),
        "__MODE__": mode,
        "__CSS__": bundle_css(),
        "__MANIFEST__": _json_for_script(manifest),
        "__BLOBS__": "\n".join(blob_html),
        "__RUNTIME__": bundle_js(),
    }
    # One pass: inserted text is never re-scanned for placeholders.
    return _PLACEHOLDER_RE.sub(lambda m: replacements[m.group(0)], template_html())


class OvizFigure:
    """An interactive figure rendered by the Oviz WebGL2 viewer.

    Build one from a scene spec (what :meth:`oviz.Animate3D.make_plot`
    produces) or from an existing :class:`Bundle`, then call
    :meth:`write_html`. ``scene_spec`` stays mutable until the HTML is
    written, so scripts can attach saved States (``fig.scene_spec["states"]``)
    exactly as with the classic figure.
    """

    def __init__(self, scene_spec: dict[str, Any] | None = None, *, bundle: Bundle | None = None,
                 theme: str = "dark", title: str | None = None, mode: str = "focus",
                 **_legacy_options: Any):
        if bundle is None and scene_spec is None:
            raise ValueError("OvizFigure needs a scene_spec or a bundle")
        #: Mode the figure opens in: one of :data:`VIEWER_MODES`.
        self.mode = _check_mode(mode)
        self.scene_spec = scene_spec
        self._bundle = bundle
        self._bundle_key: str | None = None
        self.theme = theme
        self.title = title

    @property
    def bundle(self) -> Bundle:
        """The compiled bundle (recompiled if the scene spec's content changed)."""

        if self.scene_spec is None:
            return self._bundle
        key = _spec_fingerprint(self.scene_spec)
        if self._bundle is None or key is None or key != self._bundle_key:
            self._bundle = compile_scene_spec(self.scene_spec)
            self._bundle_key = key
        return self._bundle

    # Mirror the classic figure's API so callers can switch transparently.
    def to_dict(self) -> dict[str, Any]:
        return self.scene_spec if self.scene_spec is not None else self.bundle.manifest

    def to_html(self, *, mode: str | None = None, **_: Any) -> str:
        return render_bundle_html(self.bundle, title=self.title, theme=self.theme, mode=mode or self.mode)

    def write_html(self, file: str | Path, **kwargs: Any) -> Path:
        path = Path(file)
        path.write_text(self.to_html(**kwargs), encoding="utf-8")
        return path

    def size_report(self) -> dict[str, Any]:
        return self.bundle.size_report()

    def _repr_html_(self) -> str:
        # srcdoc rather than a data: URL: Chromium does not load data: URL
        # documents over 2 MiB, and figures are usually far larger. A srcdoc
        # frame would share the notebook's origin, so the sandbox (without
        # allow-same-origin) keeps the figure, and the Aladin script it
        # loads, isolated the way a data: URL did.
        return (
            f'<iframe srcdoc="{html.escape(self.to_html(), quote=True)}" style="width:100%;height:78vh;min-height:520px;'
            'border:0;border-radius:12px;display:block" loading="eager" '
            'sandbox="allow-scripts allow-downloads allow-popups allow-popups-to-escape-sandbox" '
            'referrerpolicy="no-referrer" allow="fullscreen; clipboard-write"></iframe>'
        )

    def show(self) -> str | None:
        try:
            from IPython.display import HTML, display

            display(HTML(self._repr_html_()))
            return None
        except Exception:
            pass
        with tempfile.NamedTemporaryFile("w", suffix=".html", delete=False, encoding="utf-8") as tmp:
            tmp.write(self.to_html())
        webbrowser.open(Path(tmp.name).resolve().as_uri())
        return tmp.name
