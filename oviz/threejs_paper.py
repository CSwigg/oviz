"""Schema helpers for the Oviz paper reading mode."""

from __future__ import annotations

import uuid
from copy import deepcopy
from typing import Any

from .threejs_states import deduplicate_state_assets

PAPER_SCHEMA_VERSION = 1
SUPPORTED_BLOCK_TYPES = {"html", "figure"}
SUPPORTED_SYNC_MODES = {"scroll", "manual"}
SUPPORTED_RESUME_POLICIES = {"next-anchor", "manual"}
DEFAULT_PANEL_WIDTH_FRACTION = 0.42
DEFAULT_READING_LINE = 0.35
DEFAULT_KATEX_VERSION = "0.16.11"


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


def _unique_id(raw: Any, prefix: str, seen: set[str]) -> str:
    candidate = str(raw or "").strip() or _stable_id(prefix)
    while candidate in seen:
        candidate = _stable_id(prefix)
    seen.add(candidate)
    return candidate


def normalize_paper_anchor(value: Any, *, seen_ids: set[str]) -> dict[str, Any] | None:
    if not isinstance(value, dict):
        return None
    state_id = value.get("state_id")
    state_name = str(value.get("state") or "").strip()
    if state_id in (None, "", "original"):
        state_id = None if state_id in (None, "") else "original"
    else:
        state_id = str(state_id)
    transition = value.get("transition")
    if isinstance(transition, dict):
        transition = {
            "duration_ms": _clamp(
                transition.get("duration_ms", transition.get("duration")),
                0.0,
                60_000.0,
                1200.0,
            ),
            "easing": str(transition.get("easing") or "easeInOutCubic"),
        }
    else:
        transition = None
    return {
        "id": _unique_id(value.get("id"), "anchor", seen_ids),
        "state": state_name,
        "state_id": state_id,
        "transition": transition,
        "label": str(value.get("label") or "").strip(),
    }


def normalize_paper_block(
    value: Any,
    *,
    index: int,
    seen_block_ids: set[str],
    seen_anchor_ids: set[str],
) -> dict[str, Any]:
    source = value if isinstance(value, dict) else {}
    block_type = str(source.get("type") or "html").strip().lower()
    if block_type not in SUPPORTED_BLOCK_TYPES:
        block_type = "html"
    block: dict[str, Any] = {
        "id": _unique_id(source.get("id"), "block", seen_block_ids),
        "type": block_type,
        "anchor": normalize_paper_anchor(source.get("anchor"), seen_ids=seen_anchor_ids),
    }
    if block_type == "figure":
        block["asset"] = str(source.get("asset") or "").strip() or None
        block["image_data_url"] = (
            str(source.get("image_data_url") or "").strip() or None
        )
        block["caption_html"] = str(source.get("caption_html") or "")
        block["width_px"] = int(_clamp(source.get("width_px"), 64.0, 4096.0, 1400.0))
        block["live"] = bool(source.get("live", False))
    else:
        block["html"] = str(source.get("html") or "")
    return block


def normalize_paper_spec(value: dict[str, Any] | None) -> dict[str, Any]:
    source = value if isinstance(value, dict) else {}
    raw_sections = source.get("sections", [])
    if not isinstance(raw_sections, list):
        raw_sections = []
    seen_sections: set[str] = set()
    seen_blocks: set[str] = set()
    seen_anchors: set[str] = set()
    sections: list[dict[str, Any]] = []
    for section_index, raw_section in enumerate(raw_sections):
        section = raw_section if isinstance(raw_section, dict) else {}
        raw_blocks = section.get("blocks", [])
        if not isinstance(raw_blocks, list):
            raw_blocks = []
        sections.append(
            {
                "id": _unique_id(section.get("id"), "section", seen_sections),
                "level": int(_clamp(section.get("level"), 1.0, 4.0, 1.0)),
                "title_html": str(section.get("title_html") or section.get("title") or ""),
                "blocks": [
                    normalize_paper_block(
                        raw_block,
                        index=block_index,
                        seen_block_ids=seen_blocks,
                        seen_anchor_ids=seen_anchors,
                    )
                    for block_index, raw_block in enumerate(raw_blocks)
                ],
            }
        )

    panel = source.get("panel") if isinstance(source.get("panel"), dict) else {}
    sync = source.get("sync") if isinstance(source.get("sync"), dict) else {}
    math = source.get("math") if isinstance(source.get("math"), dict) else {}
    sync_mode = str(sync.get("mode") or "scroll").strip().lower()
    if sync_mode not in SUPPORTED_SYNC_MODES:
        sync_mode = "scroll"
    resume_policy = str(sync.get("resume_policy") or "next-anchor").strip().lower()
    if resume_policy not in SUPPORTED_RESUME_POLICIES:
        resume_policy = "next-anchor"

    authors = source.get("authors")
    if isinstance(authors, str):
        authors = [authors]
    if not isinstance(authors, list):
        authors = []

    assets = source.get("assets") if isinstance(source.get("assets"), dict) else {}
    normalized = {
        "schema_version": PAPER_SCHEMA_VERSION,
        "available": bool(source.get("available", bool(sections))),
        "enabled": bool(source.get("enabled", bool(sections))),
        "title": str(source.get("title") or ""),
        "authors": [str(author) for author in authors if str(author or "").strip()],
        "venue_html": str(source.get("venue_html") or source.get("venue") or ""),
        "link_url": str(source.get("link_url") or ""),
        "panel": {
            "width_fraction": _clamp(
                panel.get("width_fraction"), 0.25, 0.60, DEFAULT_PANEL_WIDTH_FRACTION
            ),
        },
        "sync": {
            "mode": sync_mode,
            "reading_line": _clamp(sync.get("reading_line"), 0.10, 0.90, DEFAULT_READING_LINE),
            "resume_policy": resume_policy,
        },
        "math": {
            "enabled": bool(math.get("enabled", True)),
            "katex_version": str(math.get("katex_version") or DEFAULT_KATEX_VERSION),
        },
        "sections": sections,
        "assets": {str(key): str(item) for key, item in assets.items()},
    }
    # Content-address any large inline data URLs (figures) exactly like the
    # states system so repeated images are stored once.
    compact_sections, collected_assets = deduplicate_state_assets(
        normalized["sections"], assets=dict(normalized["assets"])
    )
    normalized["sections"] = compact_sections
    normalized["assets"] = collected_assets
    return normalized


def resolve_paper_state_bindings(
    paper_spec: dict[str, Any],
    states_spec: dict[str, Any] | None,
) -> dict[str, Any]:
    """Resolve anchor ``state`` names to state ids against a states spec."""
    resolved = deepcopy(paper_spec)
    items = (states_spec or {}).get("items") or []
    by_name: dict[str, str] = {}
    by_lower: dict[str, str] = {}
    for item in items:
        if not isinstance(item, dict):
            continue
        name = str(item.get("name") or "").strip()
        state_id = str(item.get("id") or "").strip()
        if not name or not state_id:
            continue
        by_name.setdefault(name, state_id)
        by_lower.setdefault(name.lower(), state_id)
    missing: list[str] = []
    for section in resolved.get("sections", []):
        for block in section.get("blocks", []):
            anchor = block.get("anchor")
            if not isinstance(anchor, dict):
                continue
            if anchor.get("state_id"):
                continue
            name = str(anchor.get("state") or "").strip()
            if not name:
                continue
            state_id = by_name.get(name) or by_lower.get(name.lower())
            if state_id:
                anchor["state_id"] = state_id
            else:
                missing.append(name)
    if missing:
        available = ", ".join(sorted({str(item.get("name")) for item in items if isinstance(item, dict)}))
        raise ValueError(
            "Unknown paper anchor state name(s): "
            + ", ".join(sorted(set(missing)))
            + (f". Available states: {available}" if available else ". No states defined.")
        )
    return resolved
