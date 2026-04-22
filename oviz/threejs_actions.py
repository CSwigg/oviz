from __future__ import annotations

import copy
from typing import Any


_ACTION_STEP_START_VALUES = {"after_previous", "with_previous"}
_ACTION_CAMERA_TARGET_KINDS = {"group", "trace", "point"}
_ACTION_TIME_DIRECTIONS = {"forward", "backward"}


def _require_mapping(value: Any, *, label: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise ValueError(f"{label} must be a dict.")
    return value


def _coerce_non_negative_int(value: Any, *, label: str, default: int) -> int:
    if value is None:
        return int(default)
    try:
        number = int(round(float(value)))
    except Exception as exc:  # pragma: no cover - defensive
        raise ValueError(f"{label} must be a non-negative integer.") from exc
    if number < 0:
        raise ValueError(f"{label} must be a non-negative integer.")
    return number


def _coerce_positive_float(value: Any, *, label: str, default: float) -> float:
    if value is None:
        return float(default)
    try:
        number = float(value)
    except Exception as exc:  # pragma: no cover - defensive
        raise ValueError(f"{label} must be a positive number.") from exc
    if number <= 0.0:
        raise ValueError(f"{label} must be a positive number.")
    return number


def _coerce_xyz_mapping(value: Any, *, label: str) -> dict[str, float]:
    mapping = _require_mapping(value, label=label)
    coords = {}
    for axis in ("x", "y", "z"):
        if axis not in mapping:
            raise ValueError(f"{label} must include {axis}.")
        try:
            coords[axis] = float(mapping[axis])
        except Exception as exc:  # pragma: no cover - defensive
            raise ValueError(f"{label}.{axis} must be numeric.") from exc
    return coords


def _coerce_xy_mapping(value: Any, *, label: str) -> dict[str, float]:
    mapping = _require_mapping(value, label=label)
    coords = {}
    for axis in ("x", "y"):
        if axis not in mapping:
            raise ValueError(f"{label} must include {axis}.")
        try:
            coords[axis] = float(mapping[axis])
        except Exception as exc:  # pragma: no cover - defensive
            raise ValueError(f"{label}.{axis} must be numeric.") from exc
    return coords


def _normalize_camera_target(
    target_value: Any,
    *,
    group_names: set[str],
    trace_key_by_name: dict[str, str],
    label_prefix: str,
) -> dict[str, Any]:
    target = _require_mapping(target_value, label=f"{label_prefix}.target")
    kind = str(target.get("kind") or "").strip().lower()
    if kind not in _ACTION_CAMERA_TARGET_KINDS:
        raise ValueError(
            f"{label_prefix}.target.kind must be one of "
            f"{sorted(_ACTION_CAMERA_TARGET_KINDS)}."
        )

    if kind == "group":
        group_name = str(target.get("name") or "").strip()
        if not group_name:
            raise ValueError(f"{label_prefix}.target.name is required for group targets.")
        if group_name not in group_names:
            raise ValueError(f"{label_prefix}.target.name references unknown group {group_name!r}.")
        return {"kind": "group", "name": group_name}

    if kind == "trace":
        trace_name = str(target.get("name") or "").strip()
        if not trace_name:
            raise ValueError(f"{label_prefix}.target.name is required for trace targets.")
        trace_key = trace_key_by_name.get(trace_name)
        if not trace_key:
            raise ValueError(f"{label_prefix}.target.name references unknown trace {trace_name!r}.")
        return {"kind": "trace", "name": trace_name, "key": trace_key}

    return {"kind": "point", **_coerce_xyz_mapping(target, label=f"{label_prefix}.target")}


def _normalize_camera_overrides(value: Any, *, label_prefix: str) -> dict[str, Any] | None:
    if value is None:
        return None
    camera = _require_mapping(value, label=f"{label_prefix}.camera")
    normalized: dict[str, Any] = {}
    if "position" in camera:
        normalized["position"] = _coerce_xyz_mapping(
            camera.get("position"),
            label=f"{label_prefix}.camera.position",
        )
    if "target" in camera:
        normalized["target"] = _coerce_xyz_mapping(
            camera.get("target"),
            label=f"{label_prefix}.camera.target",
        )
    if "up" in camera:
        normalized["up"] = _coerce_xyz_mapping(
            camera.get("up"),
            label=f"{label_prefix}.camera.up",
        )
    if "view_offset" in camera:
        normalized["view_offset"] = _coerce_xy_mapping(
            camera.get("view_offset"),
            label=f"{label_prefix}.camera.view_offset",
        )
    if "fov" in camera:
        try:
            normalized["fov"] = float(camera.get("fov"))
        except Exception as exc:  # pragma: no cover - defensive
            raise ValueError(f"{label_prefix}.camera.fov must be numeric.") from exc
        if normalized["fov"] <= 0.0:
            raise ValueError(f"{label_prefix}.camera.fov must be positive.")
    return normalized or None


def _normalize_orbit_config(value: Any, *, label_prefix: str) -> dict[str, Any]:
    if value is None:
        return {"enabled": False, "speed_multiplier": 1.0, "persist": True, "direction": 1.0}
    orbit = _require_mapping(value, label=f"{label_prefix}.orbit")
    speed_multiplier = orbit.get("speed_multiplier", 1.0)
    try:
        speed_multiplier = float(speed_multiplier)
    except Exception as exc:  # pragma: no cover - defensive
        raise ValueError(f"{label_prefix}.orbit.speed_multiplier must be numeric.") from exc
    if speed_multiplier <= 0.0:
        raise ValueError(f"{label_prefix}.orbit.speed_multiplier must be positive.")
    direction = orbit.get("direction", 1.0)
    if isinstance(direction, str):
        normalized_direction = direction.strip().lower()
        if normalized_direction in {"reverse", "clockwise", "cw", "-1", "negative"}:
            direction = -1.0
        elif normalized_direction in {"forward", "counterclockwise", "anticlockwise", "ccw", "1", "positive"}:
            direction = 1.0
        else:
            raise ValueError(
                f"{label_prefix}.orbit.direction must be numeric or one of "
                "['forward', 'reverse', 'clockwise', 'counterclockwise']."
            )
    else:
        try:
            direction = float(direction)
        except Exception as exc:  # pragma: no cover - defensive
            raise ValueError(f"{label_prefix}.orbit.direction must be numeric.") from exc
        direction = -1.0 if direction < 0.0 else 1.0
    return {
        "enabled": bool(orbit.get("enabled")),
        "speed_multiplier": speed_multiplier,
        "persist": bool(orbit.get("persist", True)),
        "direction": direction,
    }


def normalize_threejs_actions(
    actions: Any,
    *,
    group_names: list[str],
    trace_key_by_name: dict[str, str],
    playback_interval_ms: int,
) -> list[dict[str, Any]]:
    if actions in (None, []):
        return []
    if not isinstance(actions, (list, tuple)):
        raise ValueError("actions must be a list of action definitions.")

    group_name_set = {str(name) for name in group_names}
    normalized_actions: list[dict[str, Any]] = []
    seen_action_keys: set[str] = set()
    seen_action_labels: set[str] = set()

    for action_index, raw_action in enumerate(actions):
        action_label = f"actions[{action_index}]"
        action = _require_mapping(raw_action, label=action_label)

        key = str(action.get("key") or "").strip()
        label = str(action.get("label") or "").strip()
        if not key:
            raise ValueError(f"{action_label}.key is required.")
        if not label:
            raise ValueError(f"{action_label}.label is required.")
        if key in seen_action_keys:
            raise ValueError(f"Duplicate action key {key!r}.")
        if label in seen_action_labels:
            raise ValueError(f"Duplicate action label {label!r}.")
        seen_action_keys.add(key)
        seen_action_labels.add(label)

        raw_steps = action.get("steps")
        if not isinstance(raw_steps, (list, tuple)) or not raw_steps:
            raise ValueError(f"{action_label}.steps must be a non-empty list.")

        normalized_steps: list[dict[str, Any]] = []
        for step_index, raw_step in enumerate(raw_steps):
            step_label = f"{action_label}.steps[{step_index}]"
            step = _require_mapping(raw_step, label=step_label)
            step_type = str(step.get("type") or "").strip().lower()
            if step_type not in {"legend_group", "camera", "time"}:
                raise ValueError(
                    f"{step_label}.type must be one of ['camera', 'legend_group', 'time']."
                )

            start_mode = str(step.get("start") or "after_previous").strip().lower()
            if start_mode not in _ACTION_STEP_START_VALUES:
                raise ValueError(
                    f"{step_label}.start must be one of {sorted(_ACTION_STEP_START_VALUES)}."
                )

            normalized_step: dict[str, Any] = {
                "type": step_type,
                "start": start_mode,
                "delay_ms": _coerce_non_negative_int(
                    step.get("delay_ms"),
                    label=f"{step_label}.delay_ms",
                    default=0,
                ),
            }

            if step_type == "legend_group":
                group_name = str(step.get("group") or "").strip()
                if not group_name:
                    raise ValueError(f"{step_label}.group is required.")
                if group_name not in group_name_set:
                    raise ValueError(f"{step_label}.group references unknown group {group_name!r}.")
                normalized_step.update(
                    {
                        "group": group_name,
                        "duration_ms": _coerce_non_negative_int(
                            step.get("duration_ms"),
                            label=f"{step_label}.duration_ms",
                            default=500,
                        ),
                        "easing": str(step.get("easing") or "easeInOutCubic"),
                    }
                )
            elif step_type == "camera":
                normalized_step.update(
                    {
                        "target": _normalize_camera_target(
                            step.get("target"),
                            group_names=group_name_set,
                            trace_key_by_name=trace_key_by_name,
                            label_prefix=step_label,
                        ),
                        "duration_ms": _coerce_non_negative_int(
                            step.get("duration_ms"),
                            label=f"{step_label}.duration_ms",
                            default=1200,
                        ),
                        "easing": str(step.get("easing") or "easeInOutCubic"),
                        "fit_padding": _coerce_positive_float(
                            step.get("fit_padding"),
                            label=f"{step_label}.fit_padding",
                            default=1.15,
                        ),
                        "distance_scale": _coerce_positive_float(
                            step.get("distance_scale"),
                            label=f"{step_label}.distance_scale",
                            default=1.0,
                        ),
                        "orbit": _normalize_orbit_config(
                            step.get("orbit"),
                            label_prefix=step_label,
                        ),
                    }
                )
                if step.get("anchor_target") is not None:
                    normalized_step["anchor_target"] = _normalize_camera_target(
                        step.get("anchor_target"),
                        group_names=group_name_set,
                        trace_key_by_name=trace_key_by_name,
                        label_prefix=f"{step_label}.anchor_target",
                    )
                camera_overrides = _normalize_camera_overrides(
                    step.get("camera"),
                    label_prefix=step_label,
                )
                if camera_overrides:
                    normalized_step["camera"] = camera_overrides
            else:
                direction = str(step.get("direction") or "").strip().lower()
                if direction not in _ACTION_TIME_DIRECTIONS:
                    raise ValueError(
                        f"{step_label}.direction must be one of {sorted(_ACTION_TIME_DIRECTIONS)}."
                    )
                stop_rule_count = sum(
                    key in step and step.get(key) is not None
                    for key in ("stop_at_time_myr", "stop_after_frames", "stop_after_ms")
                )
                if stop_rule_count > 1:
                    raise ValueError(
                        f"{step_label} may define at most one of stop_at_time_myr, "
                        "stop_after_frames, or stop_after_ms."
                    )
                normalized_step.update(
                    {
                        "direction": direction,
                        "interval_ms": _coerce_non_negative_int(
                            step.get("interval_ms"),
                            label=f"{step_label}.interval_ms",
                            default=playback_interval_ms,
                        ),
                    }
                )
                if normalized_step["interval_ms"] <= 0:
                    raise ValueError(f"{step_label}.interval_ms must be a positive integer.")
                if step.get("stop_at_time_myr") is not None:
                    try:
                        normalized_step["stop_at_time_myr"] = float(step.get("stop_at_time_myr"))
                    except Exception as exc:  # pragma: no cover - defensive
                        raise ValueError(f"{step_label}.stop_at_time_myr must be numeric.") from exc
                if step.get("stop_after_frames") is not None:
                    normalized_step["stop_after_frames"] = _coerce_non_negative_int(
                        step.get("stop_after_frames"),
                        label=f"{step_label}.stop_after_frames",
                        default=0,
                    )
                    if normalized_step["stop_after_frames"] <= 0:
                        raise ValueError(f"{step_label}.stop_after_frames must be positive.")
                if step.get("stop_after_ms") is not None:
                    normalized_step["stop_after_ms"] = _coerce_non_negative_int(
                        step.get("stop_after_ms"),
                        label=f"{step_label}.stop_after_ms",
                        default=0,
                    )
                    if normalized_step["stop_after_ms"] <= 0:
                        raise ValueError(f"{step_label}.stop_after_ms must be positive.")

            normalized_steps.append(normalized_step)

        normalized_actions.append(
            {
                "key": key,
                "label": label,
                "description": str(action.get("description") or ""),
                "steps": normalized_steps,
            }
        )

    return copy.deepcopy(normalized_actions)
