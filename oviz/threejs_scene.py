from __future__ import annotations

import copy

import numpy as np


def _minimal_mode_enabled(plot) -> bool:
    initial_state = getattr(plot, "threejs_initial_state", {}) or {}
    return bool(initial_state.get("minimal_mode_enabled"))


def _galactic_simple_config(plot, *, file_to_data_url, coerce_float):
    initial_state = getattr(plot, "threejs_initial_state", {}) or {}
    if not bool(initial_state.get("galactic_simple_mode_enabled")):
        return {"enabled": False}

    image_path = initial_state.get("galactic_plane_image_path")
    image_data_url = file_to_data_url(image_path) if image_path else None
    return {
        "enabled": True,
        "key": "galactic-plane-overlay",
        "image_data_url": image_data_url,
        "size_pc": float(coerce_float(initial_state.get("galactic_plane_size_pc"), 40000.0)),
        "opacity": float(np.clip(coerce_float(initial_state.get("galactic_plane_opacity"), 0.6), 0.0, 1.0)),
        "hide_below_scale_bar_pc": 200.0,
        "fade_start_scale_bar_pc": 350.0,
    }


def build_threejs_scene_spec(
    plot,
    frames,
    *,
    trace_to_plotly_json,
    coerce_range,
    format_time_label,
    coerce_float,
    file_to_data_url,
    catalog_from_frame_spec,
    annotate_point_motion_ranges,
):
    minimal_mode = _minimal_mode_enabled(plot)
    galactic_simple = _galactic_simple_config(
        plot,
        file_to_data_url=file_to_data_url,
        coerce_float=coerce_float,
    )
    if galactic_simple.get("enabled") and frames:
        galactic_simple = copy.deepcopy(galactic_simple)
        earliest_frame_json = None
        for raw_frame in frames:
            try:
                frame_json = raw_frame.to_plotly_json()
            except Exception:
                frame_json = raw_frame if isinstance(raw_frame, dict) else {}
            if earliest_frame_json is None:
                earliest_frame_json = frame_json
        galactic_simple["earliest_sun_center"] = plot._threejs_sun_position(
            earliest_frame_json or {},
            {"x": 0.0, "y": 0.0, "z": 0.0},
        )

    layout_json = plot.figure_layout.to_plotly_json()
    scene_layout = layout_json.get("scene", {})
    x_range = coerce_range(scene_layout.get("xaxis", {}).get("range"), [-1.0, 1.0])
    y_range = coerce_range(scene_layout.get("yaxis", {}).get("range"), [-1.0, 1.0])
    z_range = coerce_range(scene_layout.get("zaxis", {}).get("range"), [-1.0, 1.0])
    center = {
        "x": 0.5 * (x_range[0] + x_range[1]),
        "y": 0.5 * (y_range[0] + y_range[1]),
        "z": 0.5 * (z_range[0] + z_range[1]),
    }
    max_span = max(
        x_range[1] - x_range[0],
        y_range[1] - y_range[0],
        z_range[1] - z_range[0],
        1.0,
    )

    axis_x = scene_layout.get("xaxis", {})
    axis_y = scene_layout.get("yaxis", {})
    axis_z = scene_layout.get("zaxis", {})
    show_axes = not (
        axis_x.get("visible") is False
        and axis_y.get("visible") is False
        and axis_z.get("visible") is False
    )

    volume_layers = plot._build_threejs_volume_layers()
    trace_keys = []
    legend_items = []
    for idx, trace in enumerate(plot.initial_data):
        trace_json = trace_to_plotly_json(trace)
        trace_key = f"trace-{idx}"
        trace_spec = plot._threejs_trace_spec(
            trace_key,
            trace_json,
            minimal_mode=minimal_mode,
            galactic_simple_mode=bool(galactic_simple.get("enabled")),
        )
        if trace_spec is None:
            continue

        trace_keys.append((trace_key, trace_json.get("name")))
        if trace_spec.get("showlegend"):
            legend_items.append(
                {
                    "key": trace_key,
                    "name": trace_spec.get("name") or trace_key,
                    "color": trace_spec.get("legend_color"),
                    "kind": "trace",
                    "has_points": bool(trace_spec.get("points")),
                    "has_segments": bool(trace_spec.get("segments")),
                    "has_labels": bool(trace_spec.get("labels")),
                    "has_n_stars": bool(trace_spec.get("has_n_stars")),
                    "size_by_n_stars_default": bool(trace_spec.get("size_by_n_stars_default")),
                    "default_color": trace_spec.get("legend_color"),
                    "default_opacity": float(trace_spec.get("default_opacity", 1.0)),
                    "default_point_size": trace_spec.get("default_point_size"),
                }
            )

    trace_key_by_name = {
        str(trace_name): trace_key
        for trace_key, trace_name in trace_keys
        if isinstance(trace_name, str) and trace_name
    }

    seen_volume_state_keys = set()
    for layer in volume_layers:
        state_key = str(layer.get("state_key") or layer.get("key"))
        if state_key in seen_volume_state_keys:
            continue
        seen_volume_state_keys.add(state_key)
        legend_items.append(
            {
                "key": state_key,
                "name": str(layer.get("state_name") or layer.get("name") or state_key),
                "color": layer.get("legend_color"),
                "kind": "volume",
            }
        )

    group_visibility = {}
    for group_name, traces_list in plot.trace_grouping_dict.items():
        group_visibility[group_name] = {}
        for trace_key, trace_name in trace_keys:
            visible = plot.get_visibility(trace_name, traces_list)
            if visible == "legendonly":
                group_visibility[group_name][trace_key] = "legendonly"
            else:
                group_visibility[group_name][trace_key] = bool(visible)

    for group_name in group_visibility:
        for layer in volume_layers:
            state_key = str(layer.get("state_key") or layer.get("key"))
            group_visibility[group_name][state_key] = "legendonly"

    frames_by_time = {}
    for time_val, frame in zip(plot.time, frames):
        frames_by_time[round(float(time_val), 12)] = frame.to_plotly_json()

    frame_specs = []
    ordered_times = [float(t) for t in plot._ordered_slider_times()]
    for time_val in ordered_times:
        frame_json = frames_by_time.get(round(float(time_val), 12))
        if frame_json is None:
            continue

        traces = []
        for idx, trace_json in enumerate(frame_json.get("data", [])):
            trace_spec = plot._threejs_trace_spec(
                f"trace-{idx}",
                trace_json,
                minimal_mode=minimal_mode,
                galactic_simple_mode=bool(galactic_simple.get("enabled")),
            )
            if trace_spec is not None:
                traces.append(trace_spec)

        frame_specs.append(
            {
                "name": format_time_label(time_val),
                "time": float(time_val),
                "traces": traces,
                "decorations": plot._threejs_frame_decorations(
                    frame_json=frame_json,
                    time_value=float(time_val),
                    x_range=x_range,
                    y_range=y_range,
                    z_range=z_range,
                    fallback_center=center,
                    volume_layers=volume_layers,
                    galactic_simple_config=galactic_simple,
                ),
            }
        )

    default_group = (
        "Clusters" if "Clusters" in plot.trace_grouping_dict else list(plot.trace_grouping_dict.keys())[0]
    )
    initial_frame_index = 0
    for idx, frame_spec in enumerate(frame_specs):
        if np.isclose(float(frame_spec["time"]), 0.0, atol=1e-9):
            initial_frame_index = idx
            break

    default_sky_catalog = {}
    if getattr(plot, "enable_sky_panel", False) and frame_specs:
        default_sky_catalog = catalog_from_frame_spec(frame_specs[initial_frame_index])

    if minimal_mode:
        default_sky_catalog = {}
    else:
        annotate_point_motion_ranges(frame_specs)

    title_cfg = layout_json.get("title", {})
    if isinstance(title_cfg, dict):
        title_text = title_cfg.get("text") or ""
    else:
        title_text = str(title_cfg)

    compact_layout = {"scene": copy.deepcopy(scene_layout)}
    initial_state = {
        "click_selection_enabled": False,
        "current_group": default_group,
        **copy.deepcopy(getattr(plot, "threejs_initial_state", {}) or {}),
    }
    if minimal_mode:
        for state_key in (
            "current_selection",
            "current_selections",
            "selected_cluster_keys",
            "lasso_selection_mask",
            "widgets",
            "drawers",
            "active_manual_label_id",
            "legend_panel_state",
            "legend_panel_user_sized",
            "galactic_simple_mode_enabled",
            "galactic_plane_image_path",
            "galactic_plane_size_pc",
            "galactic_plane_opacity",
        ):
            initial_state.pop(state_key, None)
        initial_state["click_selection_enabled"] = False
        initial_state["lasso_volume_selection_enabled"] = False
        initial_state["lasso_armed"] = False

    legend_payload = []
    for item in legend_items:
        legend_payload.append(
            {
                "key": item.get("key"),
                "name": item.get("name"),
                "color": item.get("color"),
                "kind": item.get("kind"),
            }
        )

    saved_global_controls = (
        initial_state.get("global_controls") if isinstance(initial_state.get("global_controls"), dict) else {}
    )
    animation_spec = {
        "fade_in_time_default": float(coerce_float(saved_global_controls.get("fade_in_time_myr"), plot.fade_in_time)),
        "fade_in_and_out_default": (
            bool(saved_global_controls.get("fade_in_and_out_enabled"))
            if isinstance(saved_global_controls.get("fade_in_and_out_enabled"), bool)
            else bool(getattr(plot, "fade_in_and_out", False))
        ),
    }
    if not minimal_mode:
        animation_spec["focus_trace_key_default"] = (
            trace_key_by_name.get(str(plot.focus_group)) if plot.focus_group else ""
        )
        animation_spec["focus_options"] = [
            {
                "key": trace_key,
                "name": str(trace_name),
            }
            for trace_key, trace_name in trace_keys
            if isinstance(trace_name, str)
            and trace_name
            and not str(trace_name).endswith(" Track")
            and any(
                trace.get("key") == trace_key and trace.get("points")
                for frame_spec in frame_specs[:1]
                for trace in frame_spec.get("traces", [])
            )
        ]

    image_planes = []
    if galactic_simple.get("enabled") and galactic_simple.get("image_data_url"):
        image_planes.append(
            {
                "key": str(galactic_simple.get("key") or "galactic-plane-overlay"),
                "image_data_url": galactic_simple.get("image_data_url"),
                "width_pc": float(coerce_float(galactic_simple.get("size_pc"), 40000.0)),
                "height_pc": float(coerce_float(galactic_simple.get("size_pc"), 40000.0)),
                "render_order": -20,
                "hide_below_scale_bar_pc": float(
                    coerce_float(galactic_simple.get("hide_below_scale_bar_pc"), 400.0)
                ),
                "fade_start_scale_bar_pc": float(
                    coerce_float(galactic_simple.get("fade_start_scale_bar_pc"), 1000.0)
                ),
            }
        )

    return {
        "renderer": "threejs",
        "export_profile": "minimal" if minimal_mode else "full",
        "galactic_simple": {
            "enabled": bool(galactic_simple.get("enabled")),
            "track_orbit_target_to_sun": False,
        },
        "title": title_text,
        "width": int(layout_json.get("width") or 900),
        "height": int(layout_json.get("height") or 700),
        "center": center,
        "max_span": float(max_span),
        "ranges": {
            "x": x_range,
            "y": y_range,
            "z": z_range,
        },
        "layout": compact_layout if minimal_mode else layout_json,
        "axes": {
            "x": axis_x,
            "y": axis_y,
            "z": axis_z,
        },
        "theme": plot._threejs_theme(layout_json),
        "frames": frame_specs,
        "initial_frame_index": int(initial_frame_index),
        "group_order": list(plot.trace_grouping_dict.keys()),
        "default_group": default_group,
        "group_visibility": group_visibility,
        "legend": {
            "items": legend_payload if minimal_mode else legend_items,
        },
        "show_axes": bool(show_axes),
        "playback_interval_ms": 240,
        "camera_up": {"x": 0.0, "y": 0.0, "z": 1.0},
        "image_planes": image_planes,
        "sky_panel": (
            {"enabled": False}
            if minimal_mode
            else plot._build_threejs_sky_panel_spec(default_sky_catalog)
        ),
        "selection_box": (
            {"enabled": False}
            if minimal_mode
            else copy.deepcopy(getattr(plot, "selection_box_spec", {"enabled": False}))
        ),
        "age_kde": (
            {"enabled": False}
            if minimal_mode
            else plot._build_threejs_age_kde_spec(trace_key_by_name)
        ),
        "cluster_filter": (
            {"enabled": False}
            if minimal_mode
            else plot._build_threejs_cluster_filter_spec(trace_key_by_name)
        ),
        "dendrogram": (
            {"enabled": False}
            if minimal_mode
            else plot._build_threejs_dendrogram_spec(trace_key_by_name)
        ),
        "volumes": {
            "enabled": bool(volume_layers),
            "layers": volume_layers,
            "co_rotation_rate_rad_per_myr": float((plot.vo / plot.ro) * 0.001022) if plot.ro else 0.0,
        },
        "animation": animation_spec,
        "initial_state": initial_state,
        "note": "" if minimal_mode else plot._threejs_note_text(),
    }
