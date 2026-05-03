from __future__ import annotations

import copy

import numpy as np

from .threejs_actions import normalize_threejs_actions
from .threejs_profiles import normalize_threejs_initial_state


def _lite_mode_enabled(plot) -> bool:
    initial_state = normalize_threejs_initial_state(getattr(plot, "threejs_initial_state", {}) or {})
    return bool(
        initial_state.get("lite_mode_enabled")
        or initial_state.get("minimal_mode_enabled")
    )


def _galactic_simple_config(plot, *, file_to_data_url, coerce_float):
    initial_state = normalize_threejs_initial_state(getattr(plot, "threejs_initial_state", {}) or {})
    if not bool(
        initial_state.get("galactic_lite_mode_enabled")
        or initial_state.get("galactic_simple_mode_enabled")
    ):
        return {"enabled": False}

    image_enabled = initial_state.get("galaxy_image")
    if image_enabled is None:
        image_enabled = bool(
            initial_state.get("galaxy_image_path")
            or initial_state.get("galactic_plane_image_path")
        )
    image_path = initial_state.get("galaxy_image_path") or initial_state.get("galactic_plane_image_path")
    image_data_url = file_to_data_url(image_path) if image_enabled and image_path else None
    return {
        "enabled": True,
        "key": "galactic-plane-overlay",
        "track_orbit_target_to_sun": bool(initial_state.get("track_orbit_target_to_sun")),
        "image_data_url": image_data_url,
        "size_pc": float(
            coerce_float(
                initial_state.get("galaxy_image_size_pc", initial_state.get("galactic_plane_size_pc")),
                40000.0,
            )
        ),
        "opacity": float(
            np.clip(
                coerce_float(
                    initial_state.get("galaxy_image_opacity", initial_state.get("galactic_plane_opacity")),
                    0.6,
                ),
                0.0,
                1.0,
            )
        ),
        "hide_below_scale_bar_pc": 200.0,
        "fade_start_scale_bar_pc": 350.0,
    }


def _galaxy_image_config(plot, *, file_to_data_url, coerce_float):
    initial_state = normalize_threejs_initial_state(getattr(plot, "threejs_initial_state", {}) or {})
    image_enabled = initial_state.get("galaxy_image")
    image_path = initial_state.get("galaxy_image_path")
    if image_enabled is None:
        image_enabled = bool(image_path)
    if not image_enabled or not image_path:
        return {"enabled": False}

    image_data_url = file_to_data_url(image_path)
    if not image_data_url:
        return {"enabled": False}

    return {
        "enabled": True,
        "key": "galaxy-image-overlay",
        "image_data_url": image_data_url,
        "size_pc": float(coerce_float(initial_state.get("galaxy_image_size_pc"), 40000.0)),
        "opacity": float(
            np.clip(coerce_float(initial_state.get("galaxy_image_opacity"), 0.6), 0.0, 1.0)
        ),
        "hide_below_scale_bar_pc": float(
            coerce_float(initial_state.get("galaxy_image_hide_below_scale_bar_pc"), 200.0)
        ),
        "fade_start_scale_bar_pc": float(
            coerce_float(initial_state.get("galaxy_image_fade_start_scale_bar_pc"), 350.0)
        ),
        "only_at_t0": bool(initial_state.get("galaxy_image_only_at_t0", True)),
    }


def _image_plane_spec(config, *, coerce_float):
    return {
        "key": str(config.get("key") or "galaxy-image-overlay"),
        "image_data_url": config.get("image_data_url"),
        "width_pc": float(coerce_float(config.get("size_pc"), 40000.0)),
        "height_pc": float(coerce_float(config.get("size_pc"), 40000.0)),
        "render_order": -20,
        "hide_below_scale_bar_pc": float(
            coerce_float(config.get("hide_below_scale_bar_pc"), 400.0)
        ),
        "fade_start_scale_bar_pc": float(
            coerce_float(config.get("fade_start_scale_bar_pc"), 1000.0)
        ),
    }


def _sky_dome_value(raw_config, initial_state, key, default=None):
    if key in raw_config:
        return raw_config.get(key)
    return initial_state.get(f"sky_dome_{key}", default)


def _sky_dome_config(plot, *, file_to_data_url, coerce_float):
    initial_state = normalize_threejs_initial_state(getattr(plot, "threejs_initial_state", {}) or {})
    raw_config = initial_state.get("sky_dome") if isinstance(initial_state.get("sky_dome"), dict) else {}

    image_data_url = _sky_dome_value(raw_config, initial_state, "image_data_url")
    image_path = _sky_dome_value(raw_config, initial_state, "image_path")
    if not image_data_url and image_path:
        image_data_url = file_to_data_url(image_path)

    enabled = _sky_dome_value(raw_config, initial_state, "enabled")
    if enabled is None:
        enabled = bool(image_data_url)
    if not bool(enabled):
        return {"enabled": False}

    projection_metadata = _sky_dome_value(raw_config, initial_state, "projection_metadata", {})
    if isinstance(projection_metadata, dict):
        projection_metadata = copy.deepcopy(projection_metadata)
    else:
        projection_metadata = {}

    projection = str(
        _sky_dome_value(raw_config, initial_state, "projection", "MOL")
        or "MOL"
    )
    background_mode = (
        _sky_dome_value(raw_config, initial_state, "background_mode")
        or _sky_dome_value(raw_config, initial_state, "mode")
        or _sky_dome_value(raw_config, initial_state, "render_mode")
    )
    source_value = _sky_dome_value(raw_config, initial_state, "source")
    if source_value:
        source = str(source_value)
    elif image_data_url:
        source = "local_image"
    elif str(background_mode or "").strip().lower() in {"native_hips", "native-hips", "hips"}:
        source = "hips"
    elif str(background_mode or "").strip().lower() in {"hips2fits", "hips-2-fits"}:
        source = "hips2fits"
    else:
        source = "aladin"

    spec = {
        "enabled": True,
        "source": source,
        "projection": projection,
        "projection_metadata": projection_metadata,
        "force_visible": bool(_sky_dome_value(raw_config, initial_state, "force_visible", False)),
        "radius_pc": float(coerce_float(_sky_dome_value(raw_config, initial_state, "radius_pc"), 40000.0)),
        "opacity": float(
            np.clip(coerce_float(_sky_dome_value(raw_config, initial_state, "opacity"), 0.55), 0.0, 1.0)
        ),
        "flip_x": bool(_sky_dome_value(raw_config, initial_state, "flip_x", False)),
        "flip_y": bool(_sky_dome_value(raw_config, initial_state, "flip_y", False)),
        "full_opacity_scale_bar_pc": float(
            coerce_float(
                _sky_dome_value(raw_config, initial_state, "full_opacity_scale_bar_pc"),
                120.0,
            )
        ),
        "fade_out_scale_bar_pc": float(
            coerce_float(
                _sky_dome_value(raw_config, initial_state, "fade_out_scale_bar_pc"),
                360.0,
            )
        ),
    }
    if background_mode:
        spec["background_mode"] = str(background_mode)
    if image_data_url:
        spec["image_data_url"] = str(image_data_url)
    elif str(source).strip().lower() in {"hips2fits", "hips-2-fits"}:
        hips_survey = _sky_dome_value(raw_config, initial_state, "hips_survey")
        hips_frame = _sky_dome_value(raw_config, initial_state, "hips_frame", "galactic")
        if hips_survey:
            spec["hips_survey"] = str(hips_survey)
        spec["hips_frame"] = str(hips_frame or "galactic")
        spec["hips2fits_service_url"] = str(
            _sky_dome_value(
                raw_config,
                initial_state,
                "hips2fits_service_url",
                "https://alasky.cds.unistra.fr/hips-image-services/hips2fits",
            )
            or "https://alasky.cds.unistra.fr/hips-image-services/hips2fits"
        )
        spec["hips2fits_width"] = int(
            np.clip(
                coerce_float(_sky_dome_value(raw_config, initial_state, "hips2fits_width"), 8192),
                1024,
                12000,
            )
        )
        spec["hips2fits_height"] = int(
            np.clip(
                coerce_float(_sky_dome_value(raw_config, initial_state, "hips2fits_height"), 4096),
                512,
                6000,
            )
        )
        spec["hips2fits_projection"] = str(
            _sky_dome_value(raw_config, initial_state, "hips2fits_projection", "CAR")
            or "CAR"
        )
        spec["hips2fits_coordsys"] = str(
            _sky_dome_value(raw_config, initial_state, "hips2fits_coordsys", "galactic")
            or "galactic"
        )
        spec["hips2fits_format"] = str(
            _sky_dome_value(raw_config, initial_state, "hips2fits_format", "jpg")
            or "jpg"
        ).lstrip(".").lower()
        spec["hips2fits_preview_width"] = int(
            np.clip(
                coerce_float(_sky_dome_value(raw_config, initial_state, "hips2fits_preview_width"), 2048),
                512,
                12000,
            )
        )
        spec["hips2fits_preview_height"] = int(
            np.clip(
                coerce_float(_sky_dome_value(raw_config, initial_state, "hips2fits_preview_height"), 512),
                256,
                6000,
            )
        )
        spec["hips2fits_medium_width"] = int(
            np.clip(
                coerce_float(_sky_dome_value(raw_config, initial_state, "hips2fits_medium_width"), 4096),
                512,
                12000,
            )
        )
        spec["hips2fits_medium_height"] = int(
            np.clip(
                coerce_float(_sky_dome_value(raw_config, initial_state, "hips2fits_medium_height"), 2048),
                256,
                6000,
            )
        )
        spec["hips2fits_center_frame"] = str(
            _sky_dome_value(raw_config, initial_state, "hips2fits_center_frame", "galactic")
            or "galactic"
        )
        spec["hips2fits_l_deg"] = float(
            coerce_float(_sky_dome_value(raw_config, initial_state, "hips2fits_l_deg"), 0.0)
        )
        spec["hips2fits_b_deg"] = float(
            coerce_float(_sky_dome_value(raw_config, initial_state, "hips2fits_b_deg"), 0.0)
        )
        spec["hips2fits_ra_deg"] = float(
            coerce_float(_sky_dome_value(raw_config, initial_state, "hips2fits_ra_deg"), 0.0)
        )
        spec["hips2fits_dec_deg"] = float(
            coerce_float(_sky_dome_value(raw_config, initial_state, "hips2fits_dec_deg"), 0.0)
        )
        spec["hips2fits_fov_deg"] = float(
            np.clip(
                coerce_float(_sky_dome_value(raw_config, initial_state, "hips2fits_fov_deg"), 360.0),
                1.0,
                360.0,
            )
        )
    elif str(source).strip().lower() in {"hips", "native_hips", "native-hips"}:
        hips_base_url = _sky_dome_value(raw_config, initial_state, "hips_base_url")
        hips_survey = _sky_dome_value(raw_config, initial_state, "hips_survey")
        hips_frame = _sky_dome_value(raw_config, initial_state, "hips_frame", "icrs")
        if hips_base_url:
            spec["hips_base_url"] = str(hips_base_url).rstrip("/")
        if hips_survey:
            spec["hips_survey"] = str(hips_survey)
        spec["hips_frame"] = str(hips_frame or "icrs")
        spec["hips_tile_format"] = str(
            _sky_dome_value(raw_config, initial_state, "hips_tile_format", "jpg")
            or "jpg"
        ).lstrip(".").lower()
        spec["hips_allsky_order"] = int(
            np.clip(
                coerce_float(_sky_dome_value(raw_config, initial_state, "hips_allsky_order"), 3),
                0,
                6,
            )
        )
        spec["hips_tile_order"] = int(
            np.clip(
                coerce_float(_sky_dome_value(raw_config, initial_state, "hips_tile_order"), 4),
                0,
                9,
            )
        )
        spec["hips_tile_subdivisions"] = int(
            np.clip(
                coerce_float(_sky_dome_value(raw_config, initial_state, "hips_tile_subdivisions"), 16),
                2,
                64,
            )
        )
        spec["hips_allsky_tile_subdivisions"] = int(
            np.clip(
                coerce_float(_sky_dome_value(raw_config, initial_state, "hips_allsky_tile_subdivisions"), 16),
                3,
                64,
            )
        )
        spec["hips_max_active_tiles"] = int(
            np.clip(
                coerce_float(_sky_dome_value(raw_config, initial_state, "hips_max_active_tiles"), 160),
                12,
                512,
            )
        )
        spec["hips_max_concurrent_tile_loads"] = int(
            np.clip(
                coerce_float(_sky_dome_value(raw_config, initial_state, "hips_max_concurrent_tile_loads"), 8),
                1,
                32,
            )
        )
        spec["hips_startup_preload_tiles"] = int(
            np.clip(
                coerce_float(_sky_dome_value(raw_config, initial_state, "hips_startup_preload_tiles"), 96),
                0,
                256,
            )
        )
        spec["hips_startup_wait_ms"] = float(
            np.clip(
                coerce_float(_sky_dome_value(raw_config, initial_state, "hips_startup_wait_ms"), 900.0),
                0.0,
                3000.0,
            )
        )
        spec["hips_brightness"] = float(
            np.clip(
                coerce_float(_sky_dome_value(raw_config, initial_state, "hips_brightness"), 2.4),
                0.1,
                8.0,
            )
        )
        spec["hips_contrast"] = float(
            np.clip(
                coerce_float(_sky_dome_value(raw_config, initial_state, "hips_contrast"), 1.25),
                0.1,
                4.0,
            )
        )
        spec["hips_gamma"] = float(
            np.clip(
                coerce_float(_sky_dome_value(raw_config, initial_state, "hips_gamma"), 1.35),
                0.2,
                4.0,
            )
        )
        spec["hips_update_interval_ms"] = float(
            np.clip(
                coerce_float(_sky_dome_value(raw_config, initial_state, "hips_update_interval_ms"), 250.0),
                50.0,
                2000.0,
            )
        )
    else:
        spec["capture_width_px"] = int(
            np.clip(
                coerce_float(_sky_dome_value(raw_config, initial_state, "capture_width_px"), 4096),
                512,
                8192,
            )
        )
        spec["capture_height_px"] = int(
            np.clip(
                coerce_float(_sky_dome_value(raw_config, initial_state, "capture_height_px"), 2048),
                256,
                4096,
            )
        )
        spec["capture_format"] = str(
            _sky_dome_value(raw_config, initial_state, "capture_format", "image/jpeg")
            or "image/jpeg"
        )
        spec["capture_quality"] = float(
            np.clip(
                coerce_float(_sky_dome_value(raw_config, initial_state, "capture_quality"), 0.94),
                0.1,
                1.0,
            )
        )
    return spec


def _timeline_enabled(plot, frame_specs) -> bool:
    if len(frame_specs) <= 1:
        return False

    data_collection = getattr(plot, "data_collection", None)
    if data_collection is None or not hasattr(data_collection, "get_all_clusters"):
        return True

    clusters = list(data_collection.get_all_clusters() or [])
    if not clusters:
        return False

    return any(bool(getattr(cluster, "has_time_varying_geometry", True)) for cluster in clusters)


def _compact_motion_payload(motion):
    if not isinstance(motion, dict):
        return None

    compact = {}
    key = str(motion.get("key") or "").strip()
    if key:
        compact["key"] = key
    for numeric_key in ("age_now_myr", "time_myr"):
        value = motion.get(numeric_key)
        try:
            numeric_value = float(value)
        except Exception:
            continue
        if np.isfinite(numeric_value):
            compact[numeric_key] = numeric_value
    return compact or None


def _compact_threejs_frame_payloads(frame_specs):
    """Trim repeated per-frame point metadata while preserving t=0 selection data."""
    for frame_spec in frame_specs:
        keep_selection = np.isclose(float(frame_spec.get("time", np.nan)), 0.0, atol=1e-9)
        for trace in frame_spec.get("traces", []) or []:
            for point in trace.get("points", []) or []:
                if not keep_selection:
                    point.pop("selection", None)
                    point.pop("hovertext", None)
                compact_motion = _compact_motion_payload(point.get("motion"))
                if compact_motion is None:
                    point.pop("motion", None)
                else:
                    point["motion"] = compact_motion


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
    minimal_mode = _lite_mode_enabled(plot)
    galactic_simple = _galactic_simple_config(
        plot,
        file_to_data_url=file_to_data_url,
        coerce_float=coerce_float,
    )
    galaxy_image = (
        {"enabled": False}
        if galactic_simple.get("enabled")
        else _galaxy_image_config(
            plot,
            file_to_data_url=file_to_data_url,
            coerce_float=coerce_float,
        )
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
    normalized_actions = normalize_threejs_actions(
        getattr(plot, "threejs_actions", []) or [],
        group_names=list(plot.trace_grouping_dict.keys()),
        trace_key_by_name=trace_key_by_name,
        playback_interval_ms=240,
    )

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
                    galaxy_image_config=galaxy_image,
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
    sky_dome = _sky_dome_config(
        plot,
        file_to_data_url=file_to_data_url,
        coerce_float=coerce_float,
    )

    normalized_initial_state = normalize_threejs_initial_state(getattr(plot, "threejs_initial_state", {}) or {})
    compact_payload = bool(normalized_initial_state.get("compact_payload_enabled"))

    if minimal_mode:
        default_sky_catalog = {}
    else:
        annotate_point_motion_ranges(frame_specs)

    if compact_payload:
        _compact_threejs_frame_payloads(frame_specs)

    title_cfg = layout_json.get("title", {})
    if isinstance(title_cfg, dict):
        title_text = title_cfg.get("text") or ""
    else:
        title_text = str(title_cfg)

    compact_layout = {"scene": copy.deepcopy(scene_layout)}
    initial_state = {
        "click_selection_enabled": False,
        "lasso_volume_selection_enabled": True,
        "current_group": default_group,
        **copy.deepcopy(normalized_initial_state),
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
            "galactic_lite_mode_enabled",
            "galactic_plane_image_path",
            "galactic_plane_size_pc",
            "galactic_plane_opacity",
            "galaxy_image_path",
            "galaxy_image_size_pc",
            "galaxy_image_opacity",
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
        image_planes.append(_image_plane_spec(galactic_simple, coerce_float=coerce_float))
    if galaxy_image.get("enabled") and galaxy_image.get("image_data_url"):
        image_planes.append(_image_plane_spec(galaxy_image, coerce_float=coerce_float))

    return {
        "renderer": "threejs",
        "export_profile": "lite" if minimal_mode else "full",
        "timeline": {
            "enabled": bool(_timeline_enabled(plot, frame_specs)),
            "frame_count": len(frame_specs),
        },
        "galactic_simple": {
            "enabled": bool(galactic_simple.get("enabled")),
            "track_orbit_target_to_sun": bool(galactic_simple.get("track_orbit_target_to_sun")),
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
        "actions": {
            "enabled": bool(normalized_actions),
            "items": normalized_actions,
        },
        "camera_up": {"x": 0.0, "y": 0.0, "z": 1.0},
        "image_planes": image_planes,
        "sky_panel": (
            {"enabled": False}
            if minimal_mode
            else plot._build_threejs_sky_panel_spec(default_sky_catalog)
        ),
        "sky_dome": {"enabled": False} if minimal_mode else sky_dome,
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
