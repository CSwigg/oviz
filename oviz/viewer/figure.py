"""Standalone HTML figures for the Oviz viewer."""

from __future__ import annotations

import base64
import html
import json
import re
import tempfile
import webbrowser
from pathlib import Path
from typing import Any

from .build import bundle_css, bundle_js, template_html
from .bundle import Bundle, encode_base64_lines
from .compile import compile_scene_spec

VIEWER_VERSION = "2.0.0"
_PLACEHOLDER_RE = re.compile(r"__(?:TITLE|VERSION|THEME|SKIN|CSS|RUNTIME|MANIFEST|BLOBS)__")

#: Interface styles the viewer ships with (see ``web/src/ui/skins.js``).
VIEWER_STYLES = ("observatory", "spatial", "island", "cards")


def _check_style(style: str) -> str:
    style = str(style or "observatory").strip().lower()
    if style not in VIEWER_STYLES:
        raise ValueError(f"viewer style must be one of {VIEWER_STYLES}, got {style!r}")
    return style


def _spec_fingerprint(spec: dict[str, Any]) -> str:
    """Cheap change detector: identity of the big members + small ones' JSON."""

    parts = []
    for k in sorted(spec):
        v = spec[k]
        if k in ("frames", "volumes", "sky_panel", "image_planes", "dendrogram"):
            parts.append(f"{k}:{id(v)}:{len(v) if hasattr(v, '__len__') else 0}")
        else:
            try:
                parts.append(f"{k}:{json.dumps(v, sort_keys=True, default=str)}")
            except (TypeError, ValueError):
                parts.append(f"{k}:{id(v)}")
    return str(hash("|".join(parts)))


def _json_for_script(obj: Any) -> str:
    text = json.dumps(obj, separators=(",", ":"), ensure_ascii=False, allow_nan=False)
    # Never let data terminate the surrounding <script> element; \u003c is a
    # valid JSON escape (unlike "<\\!--"), so JSON.parse still accepts it.
    return text.replace("<", "\\u003c")


def render_bundle_html(bundle: Bundle, *, title: str | None = None, theme: str = "dark",
                       style: str = "observatory") -> str:
    style = _check_style(style)
    manifest = dict(bundle.manifest)
    manifest["blobs"] = [b.descriptor() for b in bundle.blobs]
    manifest["viewer"] = {"version": VIEWER_VERSION, "style": style}
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
        "__SKIN__": style,
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
                 theme: str = "dark", title: str | None = None, style: str = "observatory",
                 **_legacy_options: Any):
        if bundle is None and scene_spec is None:
            raise ValueError("OvizFigure needs a scene_spec or a bundle")
        #: Interface style the figure opens in: one of :data:`VIEWER_STYLES`.
        self.style = _check_style(style)
        self.scene_spec = scene_spec
        self._bundle = bundle
        self._bundle_key: str | None = None
        self.theme = theme
        self.title = title

    @property
    def bundle(self) -> Bundle:
        """The compiled bundle (recompiled if the scene spec changed)."""

        if self.scene_spec is None:
            return self._bundle
        key = _spec_fingerprint(self.scene_spec)
        if self._bundle is None or key != self._bundle_key:
            self._bundle = compile_scene_spec(self.scene_spec)
            self._bundle_key = key
        return self._bundle

    # Mirror the classic figure's API so callers can switch transparently.
    def to_dict(self) -> dict[str, Any]:
        return self.scene_spec if self.scene_spec is not None else self.bundle.manifest

    def to_html(self, *, style: str | None = None, **_: Any) -> str:
        return render_bundle_html(self.bundle, title=self.title, theme=self.theme, style=style or self.style)

    def write_html(self, file: str | Path, **kwargs: Any) -> Path:
        path = Path(file)
        path.write_text(self.to_html(**kwargs), encoding="utf-8")
        return path

    def size_report(self) -> dict[str, Any]:
        return self.bundle.size_report()

    def _data_url(self) -> str:
        encoded = base64.b64encode(self.to_html().encode("utf-8")).decode("ascii")
        return f"data:text/html;charset=utf-8;base64,{encoded}"

    def _repr_html_(self) -> str:
        return (
            f'<iframe src="{self._data_url()}" style="width:100%;height:78vh;min-height:520px;'
            'border:0;border-radius:12px;display:block" loading="eager" '
            'referrerpolicy="no-referrer" allow="fullscreen"></iframe>'
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
