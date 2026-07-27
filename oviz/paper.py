"""Authoring helpers for the Oviz paper reading mode.

Builds the ``scene_spec["paper"]`` section programmatically: sections of HTML
blocks, embedded figures (downscaled to data URLs), and anchors that bind
blocks to saved States by name. Name resolution happens in ``to_spec`` so a
typo fails at build time with the list of available states.
"""

from __future__ import annotations

import base64
import html as _html
import io
import mimetypes
from pathlib import Path
from typing import Any

from .threejs_paper import normalize_paper_spec, resolve_paper_state_bindings


def _escape(text: str) -> str:
    return _html.escape(str(text or ""), quote=False)


def _image_to_data_url(
    source: str | Path | bytes,
    *,
    max_width_px: int = 1600,
    jpeg_quality: int = 85,
) -> tuple[str, int]:
    """Return (data_url, output_width). Downscales via Pillow when available."""
    if isinstance(source, (str, Path)):
        path = Path(source)
        raw = path.read_bytes()
        mime = mimetypes.guess_type(str(path))[0] or "image/png"
    else:
        raw = bytes(source)
        mime = "image/png"
        if raw[:3] == b"\xff\xd8\xff":
            mime = "image/jpeg"
    width = int(max_width_px)
    try:
        from PIL import Image

        with Image.open(io.BytesIO(raw)) as image:
            image.load()
            width = image.width
            if image.width > max_width_px:
                height = max(1, round(image.height * max_width_px / image.width))
                image = image.resize((int(max_width_px), int(height)), Image.LANCZOS)
                width = image.width
            buffer = io.BytesIO()
            if image.mode in ("RGBA", "LA", "P") and mime != "image/jpeg":
                image.save(buffer, format="PNG", optimize=True)
                mime = "image/png"
            else:
                if image.mode != "RGB":
                    image = image.convert("RGB")
                image.save(buffer, format="JPEG", quality=int(jpeg_quality), optimize=True)
                mime = "image/jpeg"
            raw = buffer.getvalue()
    except ImportError:
        pass
    encoded = base64.b64encode(raw).decode("ascii")
    return f"data:{mime};base64,{encoded}", width


class Paper:
    """Programmatic builder for the paper spec."""

    def __init__(
        self,
        title: str,
        *,
        authors: list[str] | None = None,
        venue_html: str = "",
        link_url: str = "",
        panel_width_fraction: float = 0.42,
        reading_line: float = 0.35,
        resume_policy: str = "next-anchor",
        math_enabled: bool = True,
        enabled: bool = True,
    ) -> None:
        self.title = str(title)
        self.authors = list(authors or [])
        self.venue_html = str(venue_html)
        self.link_url = str(link_url)
        self.panel_width_fraction = float(panel_width_fraction)
        self.reading_line = float(reading_line)
        self.resume_policy = str(resume_policy)
        self.math_enabled = bool(math_enabled)
        self.enabled = bool(enabled)
        self._sections: list[dict[str, Any]] = []

    # -- structure -----------------------------------------------------
    def add_section(self, title: str, *, level: int = 1, section_id: str | None = None) -> str:
        section_id = section_id or f"section-{len(self._sections) + 1}"
        self._sections.append(
            {
                "id": section_id,
                "level": int(level),
                "title_html": _escape(title),
                "blocks": [],
            }
        )
        return section_id

    def _target_section(self) -> dict[str, Any]:
        if not self._sections:
            self.add_section("", level=1)
        return self._sections[-1]

    @staticmethod
    def _anchor(
        state: str | None,
        label: str | None,
        transition_ms: float | None,
        easing: str | None,
    ) -> dict[str, Any] | None:
        if not state and not label:
            return None
        anchor: dict[str, Any] = {"state": str(state or ""), "label": str(label or "")}
        if transition_ms is not None or easing is not None:
            anchor["transition"] = {
                "duration_ms": float(transition_ms if transition_ms is not None else 1200.0),
                "easing": str(easing or "easeInOutCubic"),
            }
        return anchor

    def add_html(
        self,
        html: str,
        *,
        state: str | None = None,
        label: str | None = None,
        transition_ms: float | None = None,
        easing: str | None = None,
    ) -> None:
        self._target_section()["blocks"].append(
            {
                "type": "html",
                "html": str(html),
                "anchor": self._anchor(state, label, transition_ms, easing),
            }
        )

    def add_paragraph(self, text: str, **anchor_kwargs: Any) -> None:
        self.add_html(f"<p>{_escape(text)}</p>", **anchor_kwargs)

    def add_figure(
        self,
        image: str | Path | bytes,
        caption_html: str = "",
        *,
        max_width_px: int = 1600,
        jpeg_quality: int = 85,
        live: bool = False,
        state: str | None = None,
        label: str | None = None,
        transition_ms: float | None = None,
        easing: str | None = None,
    ) -> None:
        data_url, width = _image_to_data_url(
            image, max_width_px=max_width_px, jpeg_quality=jpeg_quality
        )
        self._target_section()["blocks"].append(
            {
                "type": "figure",
                "image_data_url": data_url,
                "caption_html": str(caption_html),
                "width_px": int(width),
                "live": bool(live),
                "anchor": self._anchor(state, label, transition_ms, easing),
            }
        )

    # -- output ---------------------------------------------------------
    def to_spec(self, *, states: dict[str, Any] | None = None) -> dict[str, Any]:
        spec = normalize_paper_spec(
            {
                "available": True,
                "enabled": self.enabled,
                "title": self.title,
                "authors": self.authors,
                "venue_html": self.venue_html,
                "link_url": self.link_url,
                "panel": {"width_fraction": self.panel_width_fraction},
                "sync": {
                    "mode": "scroll",
                    "reading_line": self.reading_line,
                    "resume_policy": self.resume_policy,
                },
                "math": {"enabled": self.math_enabled},
                "sections": self._sections,
            }
        )
        if states is not None:
            spec = resolve_paper_state_bindings(spec, states)
        return spec
