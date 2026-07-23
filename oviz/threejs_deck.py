"""Schema helpers for Reveal.js-backed Oviz presentation decks."""

from __future__ import annotations

import uuid
from copy import deepcopy
from typing import Any


DECK_SCHEMA_VERSION = 2
DEFAULT_REVEAL_VERSION = "5.2.1"
DECK_WIDTH = 1600.0
DECK_HEIGHT = 900.0
SUPPORTED_BLOCK_TYPES = {"title", "subtitle", "text"}
SUPPORTED_OBJECT_KINDS = {"text", "shape"}
SUPPORTED_SHAPE_TYPES = {
    "rectangle",
    "rounded_rectangle",
    "ellipse",
    "triangle",
    "diamond",
    "line",
    "arrow",
}
SUPPORTED_TEXT_ALIGNMENTS = {"left", "center", "right"}
SUPPORTED_FONT_STYLES = {"normal", "italic"}
SUPPORTED_LIST_STYLES = {"none", "bullet", "number"}
SUPPORTED_BORDER_STYLES = {"solid", "dashed", "dotted"}


def _stable_id(prefix: str) -> str:
    return f"{prefix}-{uuid.uuid4()}"


def _finite_number(value: Any, fallback: float) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return float(fallback)
    if number != number or number in (float("inf"), float("-inf")):
        return float(fallback)
    return number


def _clamp(value: Any, minimum: float, maximum: float, fallback: float) -> float:
    return min(max(_finite_number(value, fallback), minimum), maximum)


def _normalize_color(value: Any, fallback: str = "#ffffff") -> str:
    text = str(value or "").strip()
    if len(text) == 7 and text.startswith("#"):
        try:
            int(text[1:], 16)
        except ValueError:
            pass
        else:
            return text.lower()
    return fallback


def _normalize_font_family(value: Any) -> str:
    text = str(value or "sans").strip()
    return text[:100] or "sans"


def normalize_deck_object(
    value: dict[str, Any] | None,
    *,
    index: int = 0,
    legacy_percent: bool = False,
) -> dict[str, Any]:
    """Normalize a v2 slide object, migrating v1 percentage text blocks."""

    source = value if isinstance(value, dict) else {}
    kind = str(source.get("kind") or "text").strip().lower()
    if kind not in SUPPORTED_OBJECT_KINDS:
        kind = "text"
    block_type = str(source.get("type") or "text").strip().lower()
    if block_type not in SUPPORTED_BLOCK_TYPES:
        block_type = "text"
    shape_type = str(source.get("shape_type") or source.get("shape") or "rectangle").strip().lower()
    if shape_type not in SUPPORTED_SHAPE_TYPES:
        shape_type = "rectangle"

    default_font_size = {"title": 64.0, "subtitle": 34.0, "text": 26.0}[block_type]
    default_weight = {"title": 700, "subtitle": 500, "text": 400}[block_type]
    default_height = {"title": 120.0, "subtitle": 76.0, "text": 110.0}[block_type]
    if kind == "shape":
        default_height = 180.0

    raw_x = _finite_number(source.get("x"), 112.0)
    raw_y = _finite_number(source.get("y"), 126.0 + index * 90.0)
    raw_width = _finite_number(source.get("width"), 1024.0 if kind == "text" else 280.0)
    raw_height = _finite_number(source.get("height"), default_height)
    if legacy_percent:
        raw_x *= DECK_WIDTH / 100.0
        raw_y *= DECK_HEIGHT / 100.0
        raw_width *= DECK_WIDTH / 100.0
        # v1 blocks had no height; a supplied height is nevertheless percentage based.
        if "height" in source:
            raw_height *= DECK_HEIGHT / 100.0

    text_align = str(source.get("align") or "left").strip().lower()
    if text_align not in SUPPORTED_TEXT_ALIGNMENTS:
        text_align = "left"
    font_style = str(source.get("font_style") or "normal").strip().lower()
    if font_style not in SUPPORTED_FONT_STYLES:
        font_style = "normal"
    list_style = str(source.get("list_style") or "none").strip().lower()
    if list_style not in SUPPORTED_LIST_STYLES:
        list_style = "none"
    border_style = str(source.get("border_style") or "solid").strip().lower()
    if border_style not in SUPPORTED_BORDER_STYLES:
        border_style = "solid"

    object_id = str(source.get("id") or "").strip() or _stable_id("object")
    default_name = (
        str(source.get("text") or "").strip().splitlines()[0][:36]
        if kind == "text" and str(source.get("text") or "").strip()
        else shape_type.replace("_", " ").title()
    )
    return {
        "id": object_id,
        "name": str(source.get("name") or default_name or f"Object {index + 1}").strip()[:80]
        or f"Object {index + 1}",
        "kind": kind,
        "type": block_type,
        "shape_type": shape_type,
        "text": str(source.get("text") or ""),
        "x": _clamp(raw_x, 0.0, DECK_WIDTH - 20.0, 112.0),
        "y": _clamp(raw_y, 0.0, DECK_HEIGHT - 20.0, 126.0 + index * 90.0),
        "width": _clamp(raw_width, 20.0, DECK_WIDTH, 1024.0 if kind == "text" else 280.0),
        "height": _clamp(raw_height, 20.0, DECK_HEIGHT, default_height),
        "rotation": _clamp(source.get("rotation"), -360.0, 360.0, 0.0),
        "opacity": _clamp(source.get("opacity"), 0.0, 1.0, 1.0),
        "locked": bool(source.get("locked", False)),
        "group_id": str(source.get("group_id") or "").strip() or None,
        "font_size": _clamp(source.get("font_size"), 10.0, 240.0, default_font_size),
        "font_weight": int(_clamp(source.get("font_weight"), 100.0, 900.0, default_weight)),
        "font_family": _normalize_font_family(source.get("font_family")),
        "font_style": font_style,
        "underline": bool(source.get("underline", False)),
        "strikethrough": bool(source.get("strikethrough", False)),
        "line_height": _clamp(source.get("line_height"), 0.8, 2.4, 1.08),
        "character_spacing": _clamp(source.get("character_spacing"), -5.0, 30.0, 0.0),
        "list_style": list_style,
        "color": _normalize_color(source.get("color"), "#ffffff"),
        "align": text_align,
        "text_outline": bool(source.get("text_outline", False)),
        "text_outline_color": _normalize_color(source.get("text_outline_color"), "#000000"),
        "text_outline_width": _clamp(source.get("text_outline_width"), 0.0, 8.0, 1.0),
        "text_shadow": bool(source.get("text_shadow", False)),
        "fill_color": _normalize_color(source.get("fill_color"), "#2b6cb0"),
        "fill_opacity": _clamp(source.get("fill_opacity"), 0.0, 1.0, 0.7),
        "border_color": _normalize_color(source.get("border_color"), "#ffffff"),
        "border_width": _clamp(source.get("border_width"), 0.0, 24.0, 2.0),
        "border_style": border_style,
        "corner_radius": _clamp(source.get("corner_radius"), 0.0, 200.0, 24.0),
        "shadow": bool(source.get("shadow", False)),
    }


def normalize_deck_block(value: dict[str, Any] | None, *, index: int = 0) -> dict[str, Any]:
    """Backward-compatible alias for callers that still create v1 text blocks."""

    return normalize_deck_object(value, index=index, legacy_percent=True)


def normalize_deck_spec(value: dict[str, Any] | None) -> dict[str, Any]:
    source = value if isinstance(value, dict) else {}
    source_version = int(_finite_number(source.get("schema_version"), 1.0))
    raw_slides = source.get("slides", [])
    if not isinstance(raw_slides, list):
        raw_slides = []
    seen_slides: set[str] = set()
    seen_objects: set[str] = set()
    slides: list[dict[str, Any]] = []
    for slide_index, raw_slide in enumerate(raw_slides):
        slide = raw_slide if isinstance(raw_slide, dict) else {}
        slide_id = str(slide.get("id") or "").strip() or _stable_id("slide")
        while slide_id in seen_slides:
            slide_id = _stable_id("slide")
        seen_slides.add(slide_id)
        uses_legacy_blocks = "objects" not in slide
        raw_objects = slide.get("objects", slide.get("blocks", []))
        if not isinstance(raw_objects, list):
            raw_objects = []
        objects: list[dict[str, Any]] = []
        for object_index, raw_object in enumerate(raw_objects):
            obj = normalize_deck_object(
                raw_object,
                index=object_index,
                legacy_percent=uses_legacy_blocks or source_version < 2,
            )
            while obj["id"] in seen_objects:
                obj["id"] = _stable_id("object")
            seen_objects.add(obj["id"])
            objects.append(obj)
        state_id = slide.get("state_id")
        slides.append(
            {
                "id": slide_id,
                "name": str(slide.get("name") or f"Slide {slide_index + 1}").strip()
                or f"Slide {slide_index + 1}",
                "state_id": None if state_id in (None, "", "original") else str(state_id),
                "objects": objects,
                "notes": str(slide.get("notes") or ""),
            }
        )
    reveal = source.get("reveal") if isinstance(source.get("reveal"), dict) else {}
    guides = source.get("guides") if isinstance(source.get("guides"), dict) else {}
    return {
        "schema_version": DECK_SCHEMA_VERSION,
        "available": bool(source.get("available", True)),
        "enabled": bool(source.get("enabled", bool(slides))),
        "embedded": bool(source.get("embedded", False)),
        "revision": max(int(source.get("revision") or 0), 0),
        "aspect_ratio": "16:9",
        "guides": {
            "smart": bool(guides.get("smart", True)),
            "grid": bool(guides.get("grid", False)),
            "grid_size": _clamp(guides.get("grid_size"), 5.0, 200.0, 20.0),
        },
        "reveal": {
            "version": str(reveal.get("version") or DEFAULT_REVEAL_VERSION),
            "transition": str(reveal.get("transition") or "fade"),
            "background_transition": str(reveal.get("background_transition") or "fade"),
        },
        "slides": deepcopy(slides),
    }
