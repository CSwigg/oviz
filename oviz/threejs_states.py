"""Schema helpers for ordered Three.js viewer states."""

from __future__ import annotations

import hashlib
import uuid
from copy import deepcopy
from typing import Any


STATES_SCHEMA_VERSION = 1
DEFAULT_STATE_TRANSITION = {"duration_ms": 1200, "easing": "easeInOutCubic"}
SUPPORTED_STATE_EASINGS = {"linear", "easeInOutCubic", "easeInOutQuad", "easeOutCubic"}


def _stable_id(prefix: str) -> str:
    return f"{prefix}-{uuid.uuid4()}"


def default_states_project_id(scene_spec: dict[str, Any]) -> str:
    """Return a deterministic draft key for repeated exports of one figure."""
    explicit = str(scene_spec.get("states_project_id") or scene_spec.get("project_id") or "").strip()
    if explicit:
        return explicit
    identity = ":".join(
        [
            "oviz-states",
            str(scene_spec.get("title") or "untitled"),
            str(scene_spec.get("width") or ""),
            str(scene_spec.get("height") or ""),
        ]
    )
    return f"project-{uuid.uuid5(uuid.NAMESPACE_URL, identity)}"


def normalize_transition(
    value: dict[str, Any] | None,
    *,
    fallback: dict[str, Any] | None = None,
) -> dict[str, Any]:
    base = fallback or DEFAULT_STATE_TRANSITION
    source = value if isinstance(value, dict) else {}
    try:
        duration_ms = float(source.get("duration_ms", source.get("duration", base["duration_ms"])))
    except (TypeError, ValueError):
        duration_ms = float(base["duration_ms"])
    duration_ms = min(max(duration_ms, 0.0), 60_000.0)
    easing = str(source.get("easing") or base["easing"])
    if easing not in SUPPORTED_STATE_EASINGS:
        easing = str(base["easing"])
    return {
        "duration_ms": int(duration_ms) if duration_ms.is_integer() else duration_ms,
        "easing": easing,
    }


def normalize_states_spec(
    value: dict[str, Any] | None,
    *,
    project_id: str | None = None,
) -> dict[str, Any]:
    source = value if isinstance(value, dict) else {}
    items = source.get("items", source.get("ordered_states", []))
    if not isinstance(items, list):
        items = []
    seen: set[str] = set()
    normalized_items: list[dict[str, Any]] = []
    for index, raw_item in enumerate(items):
        item = raw_item if isinstance(raw_item, dict) else {}
        state_id = str(item.get("id") or "").strip() or _stable_id("state")
        while state_id in seen:
            state_id = _stable_id("state")
        seen.add(state_id)
        transition = item.get("transition")
        normalized_items.append(
            {
                "id": state_id,
                "name": str(item.get("name") or f"State {index + 1}").strip() or f"State {index + 1}",
                "transition": normalize_transition(transition) if isinstance(transition, dict) else None,
                "snapshot": deepcopy(item.get("snapshot", item.get("state", {})))
                if isinstance(item.get("snapshot", item.get("state", {})), dict)
                else {},
                "degraded": bool(item.get("degraded", False)),
            }
        )
    return {
        "schema_version": STATES_SCHEMA_VERSION,
        "project_id": str(source.get("project_id") or project_id or _stable_id("project")),
        "revision": max(int(source.get("revision") or 0), 0),
        "synchronized_revision": max(int(source.get("synchronized_revision") or 0), 0),
        "default_mode": "present" if source.get("default_mode") == "present" else "edit",
        "default_transition": normalize_transition(source.get("default_transition")),
        "items": normalized_items,
        "assets": deepcopy(source.get("assets")) if isinstance(source.get("assets"), dict) else {},
    }


def deduplicate_state_assets(
    value: Any,
    *,
    assets: dict[str, str] | None = None,
    minimum_length: int = 4096,
) -> tuple[Any, dict[str, str]]:
    """Replace repeated large data URLs with content-addressed references."""
    output_assets = assets if assets is not None else {}
    if isinstance(value, str) and value.startswith("data:") and len(value) >= minimum_length:
        digest = hashlib.sha256(value.encode("utf-8")).hexdigest()
        output_assets[digest] = value
        return {"__oviz_asset_ref__": digest}, output_assets
    if isinstance(value, list):
        output = []
        for item in value:
            compact, output_assets = deduplicate_state_assets(
                item, assets=output_assets, minimum_length=minimum_length
            )
            output.append(compact)
        return output, output_assets
    if isinstance(value, dict):
        output_dict = {}
        for key, item in value.items():
            compact, output_assets = deduplicate_state_assets(
                item, assets=output_assets, minimum_length=minimum_length
            )
            output_dict[key] = compact
        return output_dict, output_assets
    return deepcopy(value), output_assets
