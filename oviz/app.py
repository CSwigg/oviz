from __future__ import annotations

from typing import Optional
import base64
import copy
import json

from dash import Dash, dcc, html, Input, Output, State, no_update, ctx, Patch
import plotly.graph_objects as go
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
import astropy.units as u


CUSTOMDATA_IDX_AGE_NOW = 0
CUSTOMDATA_IDX_AGE_AT_T = 1
CUSTOMDATA_IDX_L0_DEG = 2
CUSTOMDATA_IDX_B0_DEG = 3
CUSTOMDATA_IDX_DIST0_PC = 4
CUSTOMDATA_IDX_X0 = 5
CUSTOMDATA_IDX_Y0 = 6
CUSTOMDATA_IDX_Z0 = 7
CUSTOMDATA_IDX_CLUSTER_NAME = 8
CUSTOMDATA_IDX_CLUSTER_COLOR = 9

FOOTPRINT_CONE_NAME = "__oviz_sky_footprint_cone__"
FOOTPRINT_RIM_NAME = "__oviz_sky_footprint_rim__"
FOOTPRINT_TRACE_NAMES = {FOOTPRINT_CONE_NAME, FOOTPRINT_RIM_NAME}
MAX_SELECTED_MEMBER_POINTS = 1200


def _is_array_like(value):
    return isinstance(value, (list, tuple, np.ndarray))


def _coerce_plotly_array(value):
    if isinstance(value, dict) and "bdata" in value:
        try:
            decoded = np.frombuffer(base64.b64decode(value["bdata"]), dtype=np.dtype(value.get("dtype", "f8")))
            shape = value.get("shape")
            if shape:
                if isinstance(shape, str):
                    dims = tuple(int(part.strip()) for part in shape.split(",") if part.strip())
                elif _is_array_like(shape):
                    dims = tuple(int(part) for part in shape)
                else:
                    dims = ()
                if dims:
                    decoded = decoded.reshape(dims)
            return decoded
        except Exception:
            return None

    if _is_array_like(value):
        return np.asarray(value)
    return None


def _extract_age_values_from_trace(trace: dict):
    customdata = trace.get("customdata")
    if customdata is None:
        return None

    arr = _coerce_plotly_array(customdata)
    if arr is None:
        return None
    if arr.ndim != 2 or arr.shape[1] < 1:
        return None

    try:
        return arr[:, CUSTOMDATA_IDX_AGE_NOW].astype(float)
    except Exception:
        return None


def _iter_cluster_traces(fig_dict: dict):
    for trace in fig_dict.get("data", []):
        if isinstance(trace, dict) and trace.get("meta", {}).get("trace_kind") == "cluster":
            yield trace

    for frame in fig_dict.get("frames", []):
        for trace in frame.get("data", []):
            if isinstance(trace, dict) and trace.get("meta", {}).get("trace_kind") == "cluster":
                yield trace


def _get_age_bounds(fig_dict: dict):
    all_ages = []
    for trace in _iter_cluster_traces(fig_dict):
        ages = _extract_age_values_from_trace(trace)
        if ages is not None and ages.size:
            all_ages.append(ages)

    if not all_ages:
        return None

    flat = np.concatenate(all_ages)
    if flat.size == 0:
        return None

    return float(np.nanmin(flat)), float(np.nanmax(flat))


def _theme_colors_from_figure(fig_dict: dict):
    scene = fig_dict.get("layout", {}).get("scene", {})
    xaxis = scene.get("xaxis", {}) if isinstance(scene, dict) else {}
    bg = str(xaxis.get("backgroundcolor", "")).lower()

    if "rgb(0" in bg or "rgba(0" in bg or "black" in bg:
        return {
            "text": "#b7b7b7",
            "panel_bg": "rgba(30,30,30,0.75)",
            "rail": "rgba(110,110,110,0.40)",
            "track": "#b7b7b7",
            "handle_border": "#d0d0d0",
            "handle_bg": "#6d6d6d",
            "footprint": "#6ec5ff",
            "panel_solid": "#121212",
        }

    return {
        "text": "#1f1f1f",
        "panel_bg": "rgba(255,255,255,0.82)",
        "rail": "rgba(25,25,25,0.25)",
        "track": "#2f2f2f",
        "handle_border": "#2f2f2f",
        "handle_bg": "#f2f2f2",
        "footprint": "#0b67b1",
        "panel_solid": "#f7f7f7",
    }


def _apply_mask_to_trace(trace: dict, mask: np.ndarray):
    n = len(mask)

    for key in ("x", "y", "z", "text", "hovertext", "customdata"):
        val = trace.get(key)
        arr = _coerce_plotly_array(val)
        if arr is not None and len(arr) == n:
            trace[key] = arr[mask].tolist()

    marker = trace.get("marker")
    if isinstance(marker, dict):
        for mk in ("size", "color", "opacity", "symbol"):
            val = marker.get(mk)
            arr = _coerce_plotly_array(val)
            if arr is not None and len(arr) == n:
                marker[mk] = arr[mask].tolist()


def _filter_cluster_trace_by_age(trace: dict, age_min: float, age_max: float):
    ages = _extract_age_values_from_trace(trace)
    if ages is None:
        return

    mask = (ages >= float(age_min)) & (ages <= float(age_max))
    _apply_mask_to_trace(trace, mask)


def _filter_figure_by_age(fig_dict: dict, age_min: float, age_max: float):
    filtered = copy.deepcopy(fig_dict)

    for trace in filtered.get("data", []):
        if isinstance(trace, dict) and trace.get("meta", {}).get("trace_kind") == "cluster":
            _filter_cluster_trace_by_age(trace, age_min, age_max)

    for frame in filtered.get("frames", []):
        for trace in frame.get("data", []):
            if isinstance(trace, dict) and trace.get("meta", {}).get("trace_kind") == "cluster":
                _filter_cluster_trace_by_age(trace, age_min, age_max)

    return filtered


def _to_float(value):
    try:
        return float(value)
    except Exception:
        return np.nan


def _extract_selected_point(click_data: dict | None):
    if not click_data:
        return None
    points = click_data.get("points")
    if not points:
        return None

    point = points[0]
    custom = point.get("customdata")
    if custom is None:
        return None

    arr = _coerce_plotly_array(custom)
    if arr is None:
        return None
    arr = np.asarray(arr).flatten()
    if arr.size <= CUSTOMDATA_IDX_Z0:
        return None

    age_now = np.nan
    age_at_t = np.nan
    if arr.size > CUSTOMDATA_IDX_AGE_NOW:
        age_now = _to_float(arr[CUSTOMDATA_IDX_AGE_NOW])
    if arr.size > CUSTOMDATA_IDX_AGE_AT_T:
        age_at_t = _to_float(arr[CUSTOMDATA_IDX_AGE_AT_T])

    click_time_myr = np.nan
    if np.isfinite(age_now) and np.isfinite(age_at_t):
        click_time_myr = float(age_at_t - age_now)

    l_deg = _to_float(arr[CUSTOMDATA_IDX_L0_DEG])
    b_deg = _to_float(arr[CUSTOMDATA_IDX_B0_DEG])
    dist_pc = _to_float(arr[CUSTOMDATA_IDX_DIST0_PC])
    x0 = _to_float(arr[CUSTOMDATA_IDX_X0])
    y0 = _to_float(arr[CUSTOMDATA_IDX_Y0])
    z0 = _to_float(arr[CUSTOMDATA_IDX_Z0])

    if not np.isfinite(l_deg) or not np.isfinite(b_deg):
        return None

    cluster_name = None
    if arr.size > CUSTOMDATA_IDX_CLUSTER_NAME:
        cluster_name = str(arr[CUSTOMDATA_IDX_CLUSTER_NAME])

    cluster_color = None
    if arr.size > CUSTOMDATA_IDX_CLUSTER_COLOR:
        cluster_color = str(arr[CUSTOMDATA_IDX_CLUSTER_COLOR])

    return {
        "l_deg": l_deg,
        "b_deg": b_deg,
        "dist_pc": dist_pc,
        "x0": x0,
        "y0": y0,
        "z0": z0,
        "click_time_myr": click_time_myr,
        "cluster_name": cluster_name,
        "cluster_color": cluster_color,
    }


def _cone_mesh_from_axis(axis_xyz: np.ndarray, half_angle_deg: float, n_theta: int = 72, n_steps: int = 24):
    axis_len = float(np.linalg.norm(axis_xyz))
    if (not np.isfinite(axis_len)) or axis_len <= 0:
        return None

    half_angle_rad = np.deg2rad(float(half_angle_deg))
    if (not np.isfinite(half_angle_rad)) or half_angle_rad <= 0:
        return None

    w_hat = axis_xyz / axis_len
    ref = np.array([0.0, 0.0, 1.0], dtype=float)
    if abs(float(np.dot(w_hat, ref))) > 0.95:
        ref = np.array([0.0, 1.0, 0.0], dtype=float)
    u_hat = np.cross(w_hat, ref)
    u_norm = np.linalg.norm(u_hat)
    if u_norm <= 0:
        return None
    u_hat /= u_norm
    v_hat = np.cross(w_hat, u_hat)

    theta = np.linspace(0.0, 2.0 * np.pi, int(n_theta), endpoint=False)
    path = np.linspace(0.0, axis_len, int(n_steps))
    tan_half = float(np.tan(half_angle_rad))

    x_vals = []
    y_vals = []
    z_vals = []

    for s in path:
        ring_radius = s * tan_half
        center = w_hat * s
        for ang in theta:
            offset = ring_radius * ((np.cos(ang) * u_hat) + (np.sin(ang) * v_hat))
            point = center + offset
            x_vals.append(float(point[0]))
            y_vals.append(float(point[1]))
            z_vals.append(float(point[2]))

    i_idx = []
    j_idx = []
    k_idx = []

    def idx(step_i, theta_i):
        return int(step_i) * int(n_theta) + int(theta_i)

    for step_i in range(int(n_steps) - 1):
        for theta_i in range(int(n_theta)):
            theta_next = (theta_i + 1) % int(n_theta)
            p00 = idx(step_i, theta_i)
            p01 = idx(step_i, theta_next)
            p10 = idx(step_i + 1, theta_i)
            p11 = idx(step_i + 1, theta_next)
            i_idx.extend([p00, p00])
            j_idx.extend([p10, p11])
            k_idx.extend([p11, p01])

    rim_radius = axis_len * tan_half
    rim_x = []
    rim_y = []
    rim_z = []
    center = w_hat * axis_len
    for ang in np.linspace(0.0, 2.0 * np.pi, int(n_theta) + 1):
        offset = rim_radius * ((np.cos(ang) * u_hat) + (np.sin(ang) * v_hat))
        point = center + offset
        rim_x.append(float(point[0]))
        rim_y.append(float(point[1]))
        rim_z.append(float(point[2]))

    return {
        "mesh": {
            "x": x_vals,
            "y": y_vals,
            "z": z_vals,
            "i": i_idx,
            "j": j_idx,
            "k": k_idx,
        },
        "rim": {"x": rim_x, "y": rim_y, "z": rim_z},
    }


def _build_footprint_traces(selection: dict, radius_deg: float, theme_colors: dict):
    axis_xyz = np.array([selection["x0"], selection["y0"], selection["z0"]], dtype=float)
    if not np.all(np.isfinite(axis_xyz)):
        return []
    cone = _cone_mesh_from_axis(axis_xyz=axis_xyz, half_angle_deg=float(radius_deg))
    if cone is None:
        return []

    cone_trace = go.Mesh3d(
        x=cone["mesh"]["x"],
        y=cone["mesh"]["y"],
        z=cone["mesh"]["z"],
        i=cone["mesh"]["i"],
        j=cone["mesh"]["j"],
        k=cone["mesh"]["k"],
        color=theme_colors["footprint"],
        opacity=0.18,
        name=FOOTPRINT_CONE_NAME,
        hoverinfo="skip",
        showlegend=False,
    )
    rim_trace = go.Scatter3d(
        x=cone["rim"]["x"],
        y=cone["rim"]["y"],
        z=cone["rim"]["z"],
        mode="lines",
        line=dict(color=theme_colors["footprint"], width=5),
        name=FOOTPRINT_RIM_NAME,
        hoverinfo="skip",
        showlegend=False,
    )
    return [
        cone_trace.to_plotly_json(),
        rim_trace.to_plotly_json(),
    ]


def _strip_footprint_from_data(data: list):
    cleaned = []
    for trace in data:
        trace_name = trace.get("name") if isinstance(trace, dict) else None
        if trace_name not in FOOTPRINT_TRACE_NAMES:
            cleaned.append(trace)
    return cleaned


def _empty_footprint_traces(theme_colors: dict):
    cone_trace = go.Mesh3d(
        x=[],
        y=[],
        z=[],
        i=[],
        j=[],
        k=[],
        color=theme_colors["footprint"],
        opacity=0.18,
        name=FOOTPRINT_CONE_NAME,
        hoverinfo="skip",
        showlegend=False,
    )
    rim_trace = go.Scatter3d(
        x=[],
        y=[],
        z=[],
        mode="lines",
        line=dict(color=theme_colors["footprint"], width=5),
        name=FOOTPRINT_RIM_NAME,
        hoverinfo="skip",
        showlegend=False,
    )
    return [cone_trace.to_plotly_json(), rim_trace.to_plotly_json()]


def _infer_frame_time_from_data(frame: dict):
    for trace in frame.get("data", []):
        if not (isinstance(trace, dict) and trace.get("meta", {}).get("trace_kind") == "cluster"):
            continue
        custom = trace.get("customdata")
        if custom is None:
            continue
        arr = _coerce_plotly_array(custom)
        if arr is None:
            continue
        if arr.ndim != 2 or arr.shape[0] == 0 or arr.shape[1] <= CUSTOMDATA_IDX_AGE_AT_T:
            continue
        age_now = _to_float(arr[0][CUSTOMDATA_IDX_AGE_NOW])
        age_at_t = _to_float(arr[0][CUSTOMDATA_IDX_AGE_AT_T])
        if np.isfinite(age_now) and np.isfinite(age_at_t):
            return float(age_at_t - age_now)
    return np.nan


def _frame_times(fig_dict: dict):
    frames = fig_dict.get("frames", [])
    if not isinstance(frames, list) or not frames:
        return []

    slider_times = _slider_times(fig_dict)
    if len(slider_times) == len(frames):
        return [float(t) if np.isfinite(t) else np.nan for t in slider_times]

    out = []
    for frame in frames:
        t_val = _to_float(frame.get("name")) if isinstance(frame, dict) else np.nan
        if np.isfinite(t_val):
            out.append(float(t_val))
            continue
        out.append(_infer_frame_time_from_data(frame if isinstance(frame, dict) else {}))
    return out


def _overlay_indices(data: list):
    cone_idx = None
    rim_idx = None
    if not isinstance(data, list):
        return cone_idx, rim_idx
    for idx, trace in enumerate(data):
        if not isinstance(trace, dict):
            continue
        trace_name = trace.get("name")
        if trace_name == FOOTPRINT_CONE_NAME:
            cone_idx = idx
        elif trace_name == FOOTPRINT_RIM_NAME:
            rim_idx = idx
    return cone_idx, rim_idx


def _footprint_overlays(selection: dict | None, sky_radius_deg: float, theme_colors: dict, active_time_myr: float | None):
    empty_overlay = _empty_footprint_traces(theme_colors)

    full_overlay = None
    if selection is not None:
        candidate = _build_footprint_traces(selection, radius_deg=sky_radius_deg, theme_colors=theme_colors)
        if candidate:
            full_overlay = candidate

    show_current = False
    if full_overlay is not None:
        t_now = _to_float(active_time_myr)
        show_current = (not np.isfinite(t_now)) or np.isclose(float(t_now), 0.0, atol=1e-9)

    current_overlay = full_overlay if (show_current and full_overlay is not None) else empty_overlay
    t0_overlay = full_overlay if full_overlay is not None else empty_overlay
    return current_overlay, t0_overlay, empty_overlay


def _build_footprint_patch(
    current_figure: dict | None,
    selection: dict | None,
    sky_radius_deg: float,
    theme_colors: dict,
):
    if not isinstance(current_figure, dict):
        return None

    data = current_figure.get("data", [])
    cone_idx, rim_idx = _overlay_indices(data)
    if cone_idx is None or rim_idx is None:
        return None

    active_time_now = _active_slider_time(current_figure)
    current_overlay, t0_overlay, empty_overlay = _footprint_overlays(
        selection=selection,
        sky_radius_deg=sky_radius_deg,
        theme_colors=theme_colors,
        active_time_myr=active_time_now,
    )
    base_cone, base_rim = current_overlay[0], current_overlay[1]

    patch = Patch()
    patch["data"][cone_idx]["x"] = copy.deepcopy(base_cone.get("x", []))
    patch["data"][cone_idx]["y"] = copy.deepcopy(base_cone.get("y", []))
    patch["data"][cone_idx]["z"] = copy.deepcopy(base_cone.get("z", []))
    patch["data"][cone_idx]["i"] = copy.deepcopy(base_cone.get("i", []))
    patch["data"][cone_idx]["j"] = copy.deepcopy(base_cone.get("j", []))
    patch["data"][cone_idx]["k"] = copy.deepcopy(base_cone.get("k", []))
    patch["data"][rim_idx]["x"] = copy.deepcopy(base_rim.get("x", []))
    patch["data"][rim_idx]["y"] = copy.deepcopy(base_rim.get("y", []))
    patch["data"][rim_idx]["z"] = copy.deepcopy(base_rim.get("z", []))

    frames = current_figure.get("frames", [])
    if not isinstance(frames, list):
        return patch

    frame_times = _frame_times(current_figure)
    for fi, frame in enumerate(frames):
        if not isinstance(frame, dict):
            continue
        fdata = frame.get("data", [])
        fcone_idx, frim_idx = _overlay_indices(fdata)
        if fcone_idx is None or frim_idx is None:
            continue

        t_val = frame_times[fi] if fi < len(frame_times) else np.nan
        if np.isfinite(t_val) and np.isclose(float(t_val), 0.0, atol=1e-9):
            ov = t0_overlay
        else:
            ov = empty_overlay
        fcone, frim = ov[0], ov[1]

        patch["frames"][fi]["data"][fcone_idx]["x"] = copy.deepcopy(fcone.get("x", []))
        patch["frames"][fi]["data"][fcone_idx]["y"] = copy.deepcopy(fcone.get("y", []))
        patch["frames"][fi]["data"][fcone_idx]["z"] = copy.deepcopy(fcone.get("z", []))
        patch["frames"][fi]["data"][fcone_idx]["i"] = copy.deepcopy(fcone.get("i", []))
        patch["frames"][fi]["data"][fcone_idx]["j"] = copy.deepcopy(fcone.get("j", []))
        patch["frames"][fi]["data"][fcone_idx]["k"] = copy.deepcopy(fcone.get("k", []))

        patch["frames"][fi]["data"][frim_idx]["x"] = copy.deepcopy(frim.get("x", []))
        patch["frames"][fi]["data"][frim_idx]["y"] = copy.deepcopy(frim.get("y", []))
        patch["frames"][fi]["data"][frim_idx]["z"] = copy.deepcopy(frim.get("z", []))

    return patch


def _apply_footprint_to_figure(
    fig_dict: dict,
    selection: dict | None,
    sky_radius_deg: float,
    theme_colors: dict,
    active_time_myr: float | None = None,
):
    updated = copy.deepcopy(fig_dict)
    updated["data"] = _strip_footprint_from_data(updated.get("data", []))
    current_overlay, t0_overlay, empty_overlay = _footprint_overlays(
        selection=selection,
        sky_radius_deg=sky_radius_deg,
        theme_colors=theme_colors,
        active_time_myr=active_time_myr,
    )
    updated["data"].extend(copy.deepcopy(current_overlay))

    frames = updated.get("frames", [])
    if not isinstance(frames, list) or not frames:
        return updated

    frame_times = _frame_times(updated)
    for idx, frame in enumerate(frames):
        if not isinstance(frame, dict):
            continue
        frame_data = _strip_footprint_from_data(frame.get("data", []))
        t_val = frame_times[idx] if idx < len(frame_times) else np.nan
        if np.isfinite(t_val) and np.isclose(float(t_val), 0.0, atol=1e-9):
            frame_data.extend(copy.deepcopy(t0_overlay))
        else:
            frame_data.extend(copy.deepcopy(empty_overlay))
        frame["data"] = frame_data

    return updated


def _sky_readout_text(selection: dict | None, sky_radius_deg: float, sky_survey: str):
    if selection is None:
        return "Click a cluster member in the 3D panel to load sky imagery."

    cluster_name = selection.get("cluster_name")
    cluster_label = f"Cluster: {cluster_name}\n" if cluster_name else ""

    return (
        f"{cluster_label}"
        f"Selected direction (t=0): l={selection['l_deg']:.4f} deg, b={selection['b_deg']:.4f} deg\n"
        f"Beam radius: {float(sky_radius_deg):.2f} deg\n"
        f"Survey: {sky_survey}"
    )


def _sky_fov_from_radius(radius_deg: float):
    return float(np.clip(float(radius_deg) * 2.4, 1.2, 180.0))


def _build_empty_sky_srcdoc(theme_colors: dict):
    text_color = theme_colors["text"]
    bg = theme_colors["panel_solid"]
    return (
        "<!doctype html><html><head><meta charset='utf-8'></head>"
        f"<body style=\"margin:0;padding:0;background:{bg};color:{text_color};"
        "display:flex;align-items:center;justify-content:center;height:100vh;"
        "font-family:Helvetica,Arial,sans-serif;font-size:14px;\">"
        "Click a cluster member in the 3D panel to open Aladin Lite."
        "</body></html>"
    )


def _slider_times(fig_dict: dict):
    try:
        steps = fig_dict.get("layout", {}).get("sliders", [])[0].get("steps", [])
    except Exception:
        return []

    times = []
    for st in steps:
        label = st.get("label")
        t_val = _to_float(label)
        if np.isfinite(t_val):
            times.append(float(t_val))
            continue

        try:
            args_time = st.get("args", [])[0][0]
            t_val = _to_float(args_time)
            if np.isfinite(t_val):
                times.append(float(t_val))
                continue
        except Exception:
            pass

        times.append(np.nan)
    return times


def _active_slider_time(fig_dict: dict | None):
    if not isinstance(fig_dict, dict):
        return np.nan
    try:
        sliders = fig_dict.get("layout", {}).get("sliders", [])
        if not isinstance(sliders, list) or not sliders:
            return np.nan
        active_idx = int(sliders[0].get("active"))
    except Exception:
        return np.nan

    times = _slider_times(fig_dict)
    if active_idx < 0 or active_idx >= len(times):
        return np.nan
    t_val = _to_float(times[active_idx])
    return float(t_val) if np.isfinite(t_val) else np.nan


def _extract_scene_ranges(fig_dict: dict):
    scene = fig_dict.get("layout", {}).get("scene", {})
    if not isinstance(scene, dict):
        return None

    out = {}
    for ax in ("xaxis", "yaxis", "zaxis"):
        axis_obj = scene.get(ax, {})
        if not isinstance(axis_obj, dict):
            continue
        rng = axis_obj.get("range")
        if not isinstance(rng, (list, tuple)) or len(rng) != 2:
            continue
        lo = _to_float(rng[0])
        hi = _to_float(rng[1])
        if np.isfinite(lo) and np.isfinite(hi) and lo != hi:
            out[ax] = [float(lo), float(hi)]

    return out if out else None


def _compute_global_cluster_ranges(fig_dict: dict, pad_frac: float = 0.03):
    mins = np.array([np.inf, np.inf, np.inf], dtype=float)
    maxs = np.array([-np.inf, -np.inf, -np.inf], dtype=float)

    for trace in _iter_cluster_traces(fig_dict):
        x = np.asarray(trace.get("x", []), dtype=float)
        y = np.asarray(trace.get("y", []), dtype=float)
        z = np.asarray(trace.get("z", []), dtype=float)
        if x.size == 0 or y.size == 0 or z.size == 0:
            continue
        finite = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
        if not np.any(finite):
            continue
        x = x[finite]
        y = y[finite]
        z = z[finite]
        mins = np.minimum(mins, [np.nanmin(x), np.nanmin(y), np.nanmin(z)])
        maxs = np.maximum(maxs, [np.nanmax(x), np.nanmax(y), np.nanmax(z)])

    if not (np.all(np.isfinite(mins)) and np.all(np.isfinite(maxs))):
        return None

    ranges = {}
    for i, ax in enumerate(("xaxis", "yaxis", "zaxis")):
        lo = float(mins[i])
        hi = float(maxs[i])
        width = hi - lo
        if not np.isfinite(width):
            continue
        if width <= 0:
            width = 1.0
            lo -= 0.5
            hi += 0.5
        pad = max(width * float(pad_frac), 1e-6)
        ranges[ax] = [lo - pad, hi + pad]
    return ranges if ranges else None


def _scene_range_lock_from_figure(fig_dict: dict):
    explicit = _extract_scene_ranges(fig_dict)
    if explicit is not None:
        return explicit
    return _compute_global_cluster_ranges(fig_dict)


def _apply_scene_range_lock(fig_dict: dict, range_lock: dict | None):
    if not range_lock:
        return
    scene = fig_dict.setdefault("layout", {}).setdefault("scene", {})
    if not isinstance(scene, dict):
        return
    for ax, rng in range_lock.items():
        axis_obj = scene.setdefault(ax, {})
        if not isinstance(axis_obj, dict):
            axis_obj = {}
            scene[ax] = axis_obj
        axis_obj["range"] = [float(rng[0]), float(rng[1])]
        axis_obj["autorange"] = False


def _drop_scene_camera(fig_dict: dict):
    scene = fig_dict.get("layout", {}).get("scene", {})
    if isinstance(scene, dict):
        scene.pop("camera", None)


def _copy_slider_active(target_fig: dict, source_fig: dict | None):
    if not isinstance(source_fig, dict):
        return
    try:
        active = source_fig.get("layout", {}).get("sliders", [])[0].get("active")
    except Exception:
        active = None
    if active is None:
        return
    try:
        target_fig["layout"]["sliders"][0]["active"] = int(active)
    except Exception:
        pass


def _copy_updatemenus_active(target_fig: dict, source_fig: dict | None):
    if not isinstance(source_fig, dict):
        return
    src_menus = source_fig.get("layout", {}).get("updatemenus", [])
    tgt_menus = target_fig.get("layout", {}).get("updatemenus", [])
    if not (isinstance(src_menus, list) and isinstance(tgt_menus, list)):
        return
    for idx in range(min(len(src_menus), len(tgt_menus))):
        src_active = src_menus[idx].get("active") if isinstance(src_menus[idx], dict) else None
        if src_active is None:
            continue
        if isinstance(tgt_menus[idx], dict):
            tgt_menus[idx]["active"] = src_active


def _copy_trace_visibility(target_fig: dict, source_fig: dict | None):
    if not isinstance(source_fig, dict):
        return
    src_data = source_fig.get("data", [])
    tgt_data = target_fig.get("data", [])
    if not (isinstance(src_data, list) and isinstance(tgt_data, list)):
        return

    src_filtered = [tr for tr in src_data if isinstance(tr, dict) and tr.get("name") not in FOOTPRINT_TRACE_NAMES]
    tgt_filtered = [tr for tr in tgt_data if isinstance(tr, dict) and tr.get("name") not in FOOTPRINT_TRACE_NAMES]

    if len(src_filtered) == len(tgt_filtered):
        pairs = list(zip(tgt_filtered, src_filtered))
    else:
        src_by_key = {}
        for tr in src_filtered:
            key = (tr.get("name"), tr.get("type"), tr.get("legendgroup"))
            src_by_key.setdefault(key, []).append(tr)
        pairs = []
        for tr in tgt_filtered:
            key = (tr.get("name"), tr.get("type"), tr.get("legendgroup"))
            candidates = src_by_key.get(key, [])
            if candidates:
                pairs.append((tr, candidates.pop(0)))

    for tgt, src in pairs:
        if "visible" in src:
            tgt["visible"] = src["visible"]
        else:
            tgt.pop("visible", None)


def _trace_match_key(trace: dict):
    if not isinstance(trace, dict):
        return None
    meta = trace.get("meta", {})
    trace_kind = meta.get("trace_kind") if isinstance(meta, dict) else None
    return (trace.get("name"), trace.get("type"), trace.get("legendgroup"), trace_kind)


def _copy_trace_uids(target_fig: dict, source_fig: dict | None):
    if not isinstance(source_fig, dict):
        return

    src_data = source_fig.get("data", [])
    tgt_data = target_fig.get("data", [])
    if not (isinstance(src_data, list) and isinstance(tgt_data, list)):
        return

    if len(src_data) == len(tgt_data):
        pairs = list(zip(tgt_data, src_data))
    else:
        src_by_key = {}
        for tr in src_data:
            if not isinstance(tr, dict):
                continue
            src_by_key.setdefault(_trace_match_key(tr), []).append(tr)
        pairs = []
        for tr in tgt_data:
            if not isinstance(tr, dict):
                continue
            key = _trace_match_key(tr)
            candidates = src_by_key.get(key, [])
            if candidates:
                pairs.append((tr, candidates.pop(0)))

    for tgt, src in pairs:
        uid = src.get("uid") if isinstance(src, dict) else None
        if uid is not None:
            tgt["uid"] = uid

    for tr in tgt_data:
        if not isinstance(tr, dict):
            continue
        name = tr.get("name")
        if name == FOOTPRINT_CONE_NAME:
            tr["uid"] = "oviz-footprint-cone"
        elif name == FOOTPRINT_RIM_NAME:
            tr["uid"] = "oviz-footprint-rim"

    src_frames = source_fig.get("frames", [])
    tgt_frames = target_fig.get("frames", [])
    if not (isinstance(src_frames, list) and isinstance(tgt_frames, list)):
        return

    for fi in range(min(len(src_frames), len(tgt_frames))):
        src_frame = src_frames[fi] if isinstance(src_frames[fi], dict) else {}
        tgt_frame = tgt_frames[fi] if isinstance(tgt_frames[fi], dict) else {}
        src_fd = src_frame.get("data", [])
        tgt_fd = tgt_frame.get("data", [])
        if not (isinstance(src_fd, list) and isinstance(tgt_fd, list)):
            continue

        if len(src_fd) == len(tgt_fd):
            frame_pairs = list(zip(tgt_fd, src_fd))
        else:
            src_by_key = {}
            for tr in src_fd:
                if not isinstance(tr, dict):
                    continue
                src_by_key.setdefault(_trace_match_key(tr), []).append(tr)
            frame_pairs = []
            for tr in tgt_fd:
                if not isinstance(tr, dict):
                    continue
                key = _trace_match_key(tr)
                candidates = src_by_key.get(key, [])
                if candidates:
                    frame_pairs.append((tr, candidates.pop(0)))

        for tgt, src in frame_pairs:
            uid = src.get("uid") if isinstance(src, dict) else None
            if uid is not None:
                tgt["uid"] = uid

        for tr in tgt_fd:
            if not isinstance(tr, dict):
                continue
            name = tr.get("name")
            if name == FOOTPRINT_CONE_NAME:
                tr["uid"] = "oviz-footprint-cone"
            elif name == FOOTPRINT_RIM_NAME:
                tr["uid"] = "oviz-footprint-rim"


def _ensure_uirevision(fig_dict: dict):
    layout = fig_dict.setdefault("layout", {})
    layout.setdefault("uirevision", "oviz-3d")
    scene = layout.setdefault("scene", {})
    if isinstance(scene, dict):
        scene.setdefault("uirevision", "oviz-3d-scene")


def _resolve_trace_color(trace: dict):
    marker = trace.get("marker", {}) if isinstance(trace, dict) else {}
    color = marker.get("color") if isinstance(marker, dict) else None
    if isinstance(color, str) and color:
        return color
    return "white"


def _collect_cluster_centers(fig_dict: dict):
    centers = {}
    accum = {}
    for trace in fig_dict.get("data", []):
        if not (isinstance(trace, dict) and trace.get("meta", {}).get("trace_kind") == "cluster"):
            continue

        trace_color = _resolve_trace_color(trace)
        custom = trace.get("customdata")
        if custom is None:
            continue

        arr = _coerce_plotly_array(custom)
        if arr is None:
            continue
        if arr.ndim != 2 or arr.shape[1] <= CUSTOMDATA_IDX_B0_DEG:
            continue

        for row in arr:
            l_deg = _to_float(row[CUSTOMDATA_IDX_L0_DEG])
            b_deg = _to_float(row[CUSTOMDATA_IDX_B0_DEG])
            if not np.isfinite(l_deg) or not np.isfinite(b_deg):
                continue

            cluster_name = None
            if len(row) > CUSTOMDATA_IDX_CLUSTER_NAME:
                cluster_name = str(row[CUSTOMDATA_IDX_CLUSTER_NAME])
            if not cluster_name:
                continue

            color = trace_color
            if len(row) > CUSTOMDATA_IDX_CLUSTER_COLOR and isinstance(row[CUSTOMDATA_IDX_CLUSTER_COLOR], str):
                color = str(row[CUSTOMDATA_IDX_CLUSTER_COLOR])
            payload = accum.setdefault(cluster_name, {"l": [], "b": [], "color": color})
            payload["l"].append(float(l_deg))
            payload["b"].append(float(b_deg))
            payload["color"] = color

    for cluster_name, payload in accum.items():
        l_arr = np.asarray(payload["l"], dtype=float)
        b_arr = np.asarray(payload["b"], dtype=float)
        if l_arr.size == 0 or b_arr.size == 0:
            continue
        centers[cluster_name] = {
            "l_deg": float(np.nanmedian(l_arr)),
            "b_deg": float(np.nanmedian(b_arr)),
            "color": payload["color"],
        }
    return centers


def _load_cluster_members(cluster_members_file: Optional[str]):
    if not cluster_members_file:
        return None

    try:
        df = pd.read_csv(cluster_members_file, usecols=["name", "l", "b"])
    except Exception:
        return None

    if df.empty:
        return None

    df = df.rename(columns={"name": "cluster_name", "l": "l_deg", "b": "b_deg"})
    df["cluster_name"] = df["cluster_name"].astype(str)
    df["l_deg"] = pd.to_numeric(df["l_deg"], errors="coerce")
    df["b_deg"] = pd.to_numeric(df["b_deg"], errors="coerce")
    df = df.loc[df["l_deg"].notnull() & df["b_deg"].notnull()].copy()
    return df if not df.empty else None


def _group_members_by_cluster(members_df: Optional[pd.DataFrame]):
    if members_df is None or members_df.empty:
        return None

    grouped = {}
    for cluster_name, grp in members_df.groupby("cluster_name", sort=False):
        l_vals = grp["l_deg"].to_numpy(dtype=np.float64)
        b_vals = grp["b_deg"].to_numpy(dtype=np.float64)
        ra_vals, dec_vals = _gal_to_icrs(l_vals, b_vals)
        grouped[str(cluster_name)] = {
            "l_deg": l_vals,
            "b_deg": b_vals,
            "ra_deg": np.asarray(ra_vals, dtype=np.float64),
            "dec_deg": np.asarray(dec_vals, dtype=np.float64),
        }
    return grouped


def _angular_sep_deg(l1_deg, b1_deg, l2_deg_arr, b2_deg_arr):
    l1 = np.deg2rad(float(l1_deg))
    b1 = np.deg2rad(float(b1_deg))
    l2 = np.deg2rad(np.asarray(l2_deg_arr, dtype=float))
    b2 = np.deg2rad(np.asarray(b2_deg_arr, dtype=float))

    cos_sep = (
        np.sin(b1) * np.sin(b2)
        + np.cos(b1) * np.cos(b2) * np.cos(l1 - l2)
    )
    cos_sep = np.clip(cos_sep, -1.0, 1.0)
    return np.rad2deg(np.arccos(cos_sep))


def _gal_to_icrs(l_deg_arr, b_deg_arr):
    c = SkyCoord(l=np.asarray(l_deg_arr, dtype=float) * u.deg, b=np.asarray(b_deg_arr, dtype=float) * u.deg, frame="galactic")
    icrs = c.icrs
    return icrs.ra.deg, icrs.dec.deg


def _selection_equal(a: dict | None, b: dict | None):
    if a is None and b is None:
        return True
    if (a is None) != (b is None):
        return False

    if a.get("cluster_name") != b.get("cluster_name"):
        return False

    for key in ("l_deg", "b_deg", "x0", "y0", "z0"):
        av = _to_float(a.get(key))
        bv = _to_float(b.get(key))
        if not (np.isfinite(av) and np.isfinite(bv)):
            return False
        if not np.isclose(av, bv, atol=1e-8):
            return False
    return True


def _build_aladin_catalog_payload(selection: dict, sky_radius_deg: float, centers_by_name: dict, members_by_cluster: Optional[dict]):
    if selection is None:
        return []

    l0 = float(selection["l_deg"])
    b0 = float(selection["b_deg"])
    selected_cluster = selection.get("cluster_name")
    center = centers_by_name.get(selected_cluster, {}) if selected_cluster else {}
    cluster_color = center.get("color") or selection.get("cluster_color") or "white"
    cluster_label = str(selected_cluster) if selected_cluster else "Selection"
    points = []
    if members_by_cluster is not None and selected_cluster in members_by_cluster:
        rows = members_by_cluster[selected_cluster]
        l_vals = rows["l_deg"]
        b_vals = rows["b_deg"]
        ra_vals = rows.get("ra_deg")
        dec_vals = rows.get("dec_deg")

        if (
            l_vals is not None
            and b_vals is not None
            and ra_vals is not None
            and dec_vals is not None
            and l_vals.size
            and b_vals.size
            and ra_vals.size
            and dec_vals.size
        ):
            idx = np.arange(l_vals.size, dtype=int)
            if idx.size > MAX_SELECTED_MEMBER_POINTS:
                choose = np.linspace(0, idx.size - 1, int(MAX_SELECTED_MEMBER_POINTS), dtype=int)
                idx = idx[choose]

            points = [
                {
                    "l": float(l_vals[i]),
                    "b": float(b_vals[i]),
                    "ra": float(ra_vals[i]),
                    "dec": float(dec_vals[i]),
                    "label": cluster_label,
                }
                for i in idx
            ]

    if not points:
        ra_one, dec_one = _gal_to_icrs([l0], [b0])
        points = [
            {
                "l": l0,
                "b": b0,
                "ra": float(ra_one[0]),
                "dec": float(dec_one[0]),
                "label": cluster_label,
            }
        ]

    return [
        {
            "name": cluster_label,
            "color": cluster_color,
            "opacity": 1.0,
            "sourceSize": 7,
            "points": points,
        }
    ]


def _build_aladin_srcdoc(
    selection: dict,
    sky_radius_deg: float,
    sky_survey: str,
    sky_frame: str,
    theme_colors: dict,
    catalog_payload: list,
):
    fov_deg = _sky_fov_from_radius(sky_radius_deg)
    l_deg = float(selection["l_deg"])
    b_deg = float(selection["b_deg"])
    ra_sel, dec_sel = _gal_to_icrs([l_deg], [b_deg])
    ra_sel = float(ra_sel[0])
    dec_sel = float(dec_sel[0])
    survey = str(sky_survey)
    coo_frame = "galactic" if sky_frame == "galactic" else "equatorial"
    bg = theme_colors["panel_solid"]
    txt = theme_colors["text"]
    beam_color = theme_colors["footprint"]

    return f"""<!doctype html>
<html>
  <head>
    <meta charset="utf-8" />
    <meta name="viewport" content="width=device-width, initial-scale=1" />
    <style>
      html, body {{ margin: 0; padding: 0; width: 100%; height: 100%; background: {bg}; color: {txt}; overflow:hidden; }}
      #oviz-wrap {{ position: relative; width: 100%; height: 100%; }}
      #aladin-lite-div {{ width: 100%; height: 100%; }}
    </style>
  </head>
  <body>
    <div id="oviz-wrap">
      <div id="aladin-lite-div"></div>
    </div>
    <script src="https://aladin.u-strasbg.fr/AladinLite/api/v3/latest/aladin.js" charset="utf-8"></script>
    <script>
      (function() {{
        const payload = {json.dumps(catalog_payload)};
        if (typeof window.A === "undefined") {{
          console.error("Aladin Lite library failed to load.");
          return;
        }}

        A.init.then(() => {{
          const options = {{
            survey: {json.dumps(survey)},
            fov: {json.dumps(fov_deg)},
            target: {json.dumps(f"{ra_sel:.6f} {dec_sel:.6f}")},
            cooFrame: {json.dumps(coo_frame)},
            showReticle: true,
            showLayersControl: true,
            showGotoControl: true,
            showFrame: true
          }};
          const aladin = A.aladin("#aladin-lite-div", options);
          if (aladin && typeof aladin.setImageSurvey === "function") {{
            aladin.setImageSurvey({json.dumps(survey)});
          }}
          if (aladin && typeof aladin.gotoRaDec === "function") {{
            aladin.gotoRaDec({ra_sel:.8f}, {dec_sel:.8f});
          }}

          const beam = A.graphicOverlay({{ color: {json.dumps(beam_color)}, lineWidth: 2, opacity: 0.95 }});
          aladin.addOverlay(beam);
          beam.add(A.circle({ra_sel:.8f}, {dec_sel:.8f}, {float(sky_radius_deg):.8f}));

          payload.forEach((catDef) => {{
            const cat = A.catalog({{
              name: catDef.name,
              color: catDef.color,
              sourceSize: catDef.sourceSize || 7,
              shape: "circle",
              opacity: catDef.opacity
            }});
            aladin.addCatalog(cat);
            const sources = [];
            (catDef.points || []).forEach((pt) => {{
              sources.push(A.source(Number(pt.ra), Number(pt.dec), {{ popupTitle: pt.label || catDef.name }}));
            }});
            if (sources.length) {{
              cat.addSources(sources);
            }}
          }});
        }}).catch((err) => {{
          console.error("Aladin initialization error:", err);
        }});
      }})();
    </script>
  </body>
</html>
"""


def _build_age_controls(theme_colors: dict, age_bounds):
    controls = []
    if age_bounds is None:
        return controls

    age_min, age_max = age_bounds
    controls.append(
        html.Div(
            [
                html.Div("Age (Myr)", style={"fontWeight": 600, "marginBottom": "2px", "fontSize": "12px"}),
                dcc.RangeSlider(
                    id="oviz-age-range",
                    min=age_min,
                    max=age_max,
                    step=0.1,
                    value=[age_min, age_max],
                    allowCross=False,
                    marks={},
                    vertical=True,
                    verticalHeight=300,
                    tooltip={"placement": "bottom", "always_visible": False},
                    className="oviz-age-slider",
                ),
            ],
            style={
                "position": "absolute",
                "right": "8px",
                "top": "50%",
                "transform": "translateY(-50%)",
                "zIndex": 25,
                "padding": "4px",
                "width": "58px",
                "background": "transparent",
                "display": "flex",
                "flexDirection": "column",
                "alignItems": "center",
                "gap": "2px",
                "color": theme_colors["text"],
            },
        )
    )
    return controls


def create_dash_app(
    figure: go.Figure,
    title: str = "oviz",
    graph_id: str = "oviz-graph",
    enable_age_filter: bool = True,
    enable_sky_panel: bool = False,
    sky_radius_deg: float = 1.0,
    sky_frame: str = "galactic",
    sky_survey: str = "P/DSS2/color",
    cluster_members_file: Optional[str] = None,
) -> Dash:
    """Create a Dash app that renders an existing Plotly figure."""
    if sky_radius_deg <= 0:
        raise ValueError("sky_radius_deg must be > 0.")
    if sky_frame not in ("galactic", "icrs"):
        raise ValueError("sky_frame must be either 'galactic' or 'icrs'.")
    if sky_frame != "galactic":
        raise ValueError("sky_frame='icrs' is not yet supported; use 'galactic'.")

    fig_local = go.Figure(figure)
    fig_local.update_layout(autosize=True)
    if fig_local.layout and fig_local.layout.to_plotly_json().get("scene") is not None:
        fig_local.layout.scene.uirevision = "oviz-3d-scene"
    fig_local.layout.uirevision = "oviz-3d"
    fig_local.update_layout(width=None, height=None)

    app = Dash(__name__)
    app.title = title

    base_figure_dict = fig_local.to_dict()
    scene_range_lock = _scene_range_lock_from_figure(base_figure_dict)
    _apply_scene_range_lock(base_figure_dict, scene_range_lock)
    _ensure_uirevision(base_figure_dict)
    age_bounds = _get_age_bounds(base_figure_dict)
    theme_colors = _theme_colors_from_figure(base_figure_dict)

    if enable_sky_panel:
        base_figure_dict = _apply_footprint_to_figure(
            fig_dict=base_figure_dict,
            selection=None,
            sky_radius_deg=sky_radius_deg,
            theme_colors=theme_colors,
            active_time_myr=0.0,
        )
        _apply_scene_range_lock(base_figure_dict, scene_range_lock)
        _ensure_uirevision(base_figure_dict)

    fig_local = go.Figure(base_figure_dict)
    centers_by_name = _collect_cluster_centers(base_figure_dict)
    members_df = _load_cluster_members(cluster_members_file)
    members_by_cluster = _group_members_by_cluster(members_df)

    slider_css = f"""
    html, body, #react-entry-point, #_dash-app-content {{
        margin: 0;
        padding: 0;
        width: 100%;
        height: 100%;
        overflow: hidden;
        background: {theme_colors['panel_solid']};
    }}
    .oviz-age-slider .rc-slider-rail {{
        background-color: {theme_colors['rail']};
        width: 4px;
    }}
    .oviz-age-slider .rc-slider-track {{
        background-color: {theme_colors['track']};
        width: 4px;
    }}
    .oviz-age-slider .rc-slider-handle {{
        border: 2px solid {theme_colors['handle_border']};
        background-color: {theme_colors['handle_bg']};
        box-shadow: none;
        width: 14px;
        height: 14px;
        margin-left: -5px;
    }}
    .oviz-age-slider .rc-slider-handle:hover,
    .oviz-age-slider .rc-slider-handle:active,
    .oviz-age-slider .rc-slider-handle:focus {{
        border-color: {theme_colors['handle_border']};
        box-shadow: none;
    }}
    .oviz-age-slider .rc-slider-mark,
    .oviz-age-slider .rc-slider-dot {{
        display: none;
    }}
    """

    app.index_string = f"""
    <!DOCTYPE html>
    <html>
        <head>
            {{%metas%}}
            <title>{{%title%}}</title>
            {{%favicon%}}
            {{%css%}}
            <style>{slider_css}</style>
        </head>
        <body>
            {{%app_entry%}}
            <footer>
                {{%config%}}
                {{%scripts%}}
                {{%renderer%}}
            </footer>
        </body>
    </html>
    """

    controls = []
    age_filter_enabled = bool(enable_age_filter and age_bounds is not None)
    if age_filter_enabled:
        controls = _build_age_controls(theme_colors, age_bounds)

    sky_container_normal_style = {
        "position": "absolute",
        "top": "10px",
        "right": "10px",
        "width": "min(33vw, 420px)",
        "height": "min(33vw, 420px)",
        "aspectRatio": "1 / 1",
        "zIndex": 40,
        "overflow": "hidden",
        "background": theme_colors["panel_solid"],
        "border": "1px solid rgba(140,140,140,0.5)",
        "borderRadius": "8px",
        "boxShadow": "0 8px 24px rgba(0,0,0,0.30)",
    }
    sky_container_fullscreen_style = {
        "position": "absolute",
        "top": "0",
        "left": "0",
        "width": "100vw",
        "height": "100vh",
        "zIndex": 80,
        "overflow": "hidden",
        "background": theme_colors["panel_solid"],
        "border": "0",
        "borderRadius": "0",
        "boxShadow": "none",
    }
    sky_show_btn_base_style = {
        "position": "absolute",
        "top": "10px",
        "right": "10px",
        "zIndex": 45,
        "padding": "6px 10px",
        "fontSize": "11px",
        "fontFamily": "Menlo,Monaco,Consolas,monospace",
        "color": theme_colors["text"],
        "background": theme_colors["panel_bg"],
        "border": "1px solid rgba(140,140,140,0.7)",
        "borderRadius": "6px",
        "cursor": "pointer",
        "display": "none",
    }
    sky_show_full_btn_base_style = {
        "position": "absolute",
        "top": "44px",
        "right": "10px",
        "zIndex": 45,
        "padding": "6px 10px",
        "fontSize": "11px",
        "fontFamily": "Menlo,Monaco,Consolas,monospace",
        "color": theme_colors["text"],
        "background": theme_colors["panel_bg"],
        "border": "1px solid rgba(140,140,140,0.7)",
        "borderRadius": "6px",
        "cursor": "pointer",
        "display": "none",
    }
    sky_drag_handle_base_style = {
        "position": "absolute",
        "top": "0",
        "left": "0",
        "right": "100px",
        "height": "24px",
        "zIndex": 31,
        "cursor": "move",
        "background": "rgba(120,120,120,0.14)",
        "borderBottom": "1px solid rgba(140,140,140,0.35)",
        "display": "block",
        "userSelect": "none",
        "touchAction": "none",
    }
    sky_resize_nw_style = {
        "position": "absolute",
        "top": "0",
        "left": "0",
        "width": "14px",
        "height": "14px",
        "zIndex": 32,
        "cursor": "nwse-resize",
        "background": "rgba(120,120,120,0.10)",
        "display": "block",
        "touchAction": "none",
    }
    sky_resize_ne_style = {
        "position": "absolute",
        "top": "0",
        "right": "0",
        "width": "14px",
        "height": "14px",
        "zIndex": 32,
        "cursor": "nesw-resize",
        "background": "rgba(120,120,120,0.10)",
        "display": "block",
        "touchAction": "none",
    }
    sky_resize_sw_style = {
        "position": "absolute",
        "bottom": "0",
        "left": "0",
        "width": "14px",
        "height": "14px",
        "zIndex": 32,
        "cursor": "nesw-resize",
        "background": "rgba(120,120,120,0.10)",
        "display": "block",
        "touchAction": "none",
    }
    sky_resize_se_style = {
        "position": "absolute",
        "bottom": "0",
        "right": "0",
        "width": "14px",
        "height": "14px",
        "zIndex": 32,
        "cursor": "nwse-resize",
        "background": "rgba(120,120,120,0.10)",
        "display": "block",
        "touchAction": "none",
    }

    graph_children = [
        dcc.Graph(
            id=graph_id,
            figure=fig_local,
            style={"position": "absolute", "inset": "0", "height": "100%", "width": "100%"},
            config={"displaylogo": False, "responsive": True},
        ),
        *controls,
    ]

    if enable_sky_panel:
        graph_children.append(
            html.Div(
                [
                    html.Div(
                        id="oviz-sky-drag-handle",
                        title="Drag sky panel",
                        style=sky_drag_handle_base_style,
                    ),
                    html.Iframe(
                        id="oviz-sky-frame",
                        srcDoc=_build_empty_sky_srcdoc(theme_colors),
                        style={
                            "width": "100%",
                            "height": "100%",
                            "border": "0",
                            "display": "block",
                            "background": theme_colors["panel_solid"],
                        },
                    ),
                    html.Pre(
                        id="oviz-sky-readout",
                        children=_sky_readout_text(None, sky_radius_deg=sky_radius_deg, sky_survey=sky_survey),
                        style={
                            "display": "none",
                        },
                    ),
                    html.Button(
                        "Full",
                        id="oviz-sky-size-btn",
                        n_clicks=0,
                        title="Toggle fullscreen sky panel",
                        style={
                            "position": "absolute",
                            "top": "6px",
                            "right": "52px",
                            "zIndex": 30,
                            "padding": "2px 6px",
                            "fontSize": "10px",
                            "fontFamily": "Menlo,Monaco,Consolas,monospace",
                            "color": theme_colors["text"],
                            "background": theme_colors["panel_bg"],
                            "border": "1px solid rgba(140,140,140,0.6)",
                            "borderRadius": "4px",
                            "cursor": "pointer",
                        },
                    ),
                    html.Button(
                        "Hide",
                        id="oviz-sky-hide-btn",
                        n_clicks=0,
                        title="Hide sky panel",
                        style={
                            "position": "absolute",
                            "top": "6px",
                            "right": "6px",
                            "zIndex": 30,
                            "padding": "2px 6px",
                            "fontSize": "10px",
                            "fontFamily": "Menlo,Monaco,Consolas,monospace",
                            "color": theme_colors["text"],
                            "background": theme_colors["panel_bg"],
                            "border": "1px solid rgba(140,140,140,0.6)",
                            "borderRadius": "4px",
                            "cursor": "pointer",
                        },
                    ),
                    html.Div(
                        id="oviz-sky-resize-nw",
                        className="oviz-sky-resize-handle",
                        **{"data-dir": "nw"},
                        style=sky_resize_nw_style,
                    ),
                    html.Div(
                        id="oviz-sky-resize-ne",
                        className="oviz-sky-resize-handle",
                        **{"data-dir": "ne"},
                        style=sky_resize_ne_style,
                    ),
                    html.Div(
                        id="oviz-sky-resize-sw",
                        className="oviz-sky-resize-handle",
                        **{"data-dir": "sw"},
                        style=sky_resize_sw_style,
                    ),
                    html.Div(
                        id="oviz-sky-resize-se",
                        className="oviz-sky-resize-handle",
                        **{"data-dir": "se"},
                        style=sky_resize_se_style,
                    ),
                ],
                id="oviz-sky-container",
                style=sky_container_normal_style,
            )
        )
        graph_children.append(
            html.Button(
                "Sky",
                id="oviz-sky-show-btn",
                n_clicks=0,
                title="Show sky panel",
                style=sky_show_btn_base_style,
            )
        )
        graph_children.append(
            html.Button(
                "Sky Full",
                id="oviz-sky-show-full-btn",
                n_clicks=0,
                title="Show fullscreen sky panel",
                style=sky_show_full_btn_base_style,
            )
        )

    layout_children = [
        dcc.Store(id="oviz-base-figure", data=base_figure_dict),
        dcc.Store(id="oviz-selection-store", data=None),
        dcc.Store(id="oviz-age-range-store", data=list(age_bounds) if age_filter_enabled else None),
    ]
    if enable_sky_panel:
        layout_children.append(dcc.Store(id="oviz-sky-mode-store", data="normal"))
    layout_children.append(
        html.Div(
            graph_children,
            style={
                "position": "relative",
                "width": "100%",
                "height": "100%",
                "overflow": "hidden",
                "background": theme_colors["panel_solid"],
            },
        )
    )

    app.layout = html.Div(
        layout_children,
        style={
            "width": "100vw",
            "height": "100vh",
            "margin": "0",
            "padding": "0",
            "overflow": "hidden",
            "background": theme_colors["panel_solid"],
        },
    )

    if enable_sky_panel:
        @app.callback(
            Output("oviz-sky-mode-store", "data"),
            Input("oviz-sky-hide-btn", "n_clicks"),
            Input("oviz-sky-size-btn", "n_clicks"),
            Input("oviz-sky-show-btn", "n_clicks"),
            Input("oviz-sky-show-full-btn", "n_clicks"),
            State("oviz-sky-mode-store", "data"),
            prevent_initial_call=True,
        )
        def _toggle_sky_panel_mode(n_hide, n_size, n_show, n_show_full, current_mode):
            triggered = ctx.triggered_id
            mode_now = str(current_mode) if current_mode in {"normal", "fullscreen", "hidden"} else "normal"
            if triggered == "oviz-sky-hide-btn":
                return "hidden"
            if triggered == "oviz-sky-size-btn":
                return "normal" if mode_now == "fullscreen" else "fullscreen"
            if triggered == "oviz-sky-show-btn":
                return "normal"
            if triggered == "oviz-sky-show-full-btn":
                return "fullscreen"
            return mode_now

        @app.callback(
            Output("oviz-sky-container", "style"),
            Output("oviz-sky-show-btn", "style"),
            Output("oviz-sky-show-full-btn", "style"),
            Output("oviz-sky-size-btn", "children"),
            Output("oviz-sky-drag-handle", "style"),
            Input("oviz-sky-mode-store", "data"),
            prevent_initial_call=False,
        )
        def _sync_sky_panel_visibility_styles(mode):
            mode_now = str(mode) if mode in {"normal", "fullscreen", "hidden"} else "normal"
            panel_style = dict(sky_container_normal_style)
            show_btn_style = dict(sky_show_btn_base_style)
            show_full_btn_style = dict(sky_show_full_btn_base_style)
            size_label = "Full"
            drag_handle_style = dict(sky_drag_handle_base_style)

            if mode_now == "fullscreen":
                panel_style = dict(sky_container_fullscreen_style)
                size_label = "Restore"
                drag_handle_style["display"] = "none"
            elif mode_now == "hidden":
                panel_style["display"] = "none"
                show_btn_style["display"] = "block"
                show_full_btn_style["display"] = "block"
                drag_handle_style["display"] = "none"

            return (
                panel_style,
                show_btn_style,
                show_full_btn_style,
                size_label,
                drag_handle_style,
            )

    def _is_t0_click(selection: dict | None):
        if not isinstance(selection, dict):
            return False
        click_time = _to_float(selection.get("click_time_myr"))
        return np.isfinite(click_time) and np.isclose(click_time, 0.0, atol=1e-9)

    if age_filter_enabled and enable_sky_panel:
        @app.callback(
            Output(graph_id, "figure"),
            Output("oviz-sky-frame", "srcDoc"),
            Output("oviz-sky-readout", "children"),
            Output("oviz-selection-store", "data"),
            Output("oviz-age-range-store", "data"),
            Input("oviz-age-range", "value"),
            Input(graph_id, "clickData"),
            State("oviz-base-figure", "data"),
            State(graph_id, "figure"),
            State("oviz-selection-store", "data"),
            State("oviz-age-range-store", "data"),
            prevent_initial_call=False,
        )
        def _update_age_and_sky(
            age_range,
            click_data,
            base_figure,
            current_figure,
            selection_store,
            age_range_store,
        ):
            age_now = None
            if age_range and len(age_range) == 2:
                age_now = [float(age_range[0]), float(age_range[1])]

            if age_now is None:
                if isinstance(age_range_store, (list, tuple)) and len(age_range_store) == 2:
                    age_now = [float(age_range_store[0]), float(age_range_store[1])]
                else:
                    age_now = [float(age_bounds[0]), float(age_bounds[1])]

            age_changed = True
            if isinstance(age_range_store, (list, tuple)) and len(age_range_store) == 2:
                age_changed = not (
                    np.isclose(float(age_range_store[0]), age_now[0], atol=1e-12)
                    and np.isclose(float(age_range_store[1]), age_now[1], atol=1e-12)
                )

            fig_base = _filter_figure_by_age(base_figure, age_now[0], age_now[1])
            _apply_scene_range_lock(fig_base, scene_range_lock)
            _ensure_uirevision(fig_base)

            selection_prev = selection_store if isinstance(selection_store, dict) else None
            selection = selection_prev
            click_selection = _extract_selected_point(click_data)
            if _is_t0_click(click_selection):
                selection = click_selection

            selection_changed = not _selection_equal(selection_prev, selection)
            figure_changed = age_changed or selection_changed

            if figure_changed:
                fig_out = None
                if (not age_changed) and isinstance(current_figure, dict):
                    fig_out = _build_footprint_patch(
                        current_figure=current_figure,
                        selection=selection,
                        sky_radius_deg=sky_radius_deg,
                        theme_colors=theme_colors,
                    )

                if fig_out is None:
                    active_time_now = _active_slider_time(current_figure)
                    fig_out = _apply_footprint_to_figure(
                        fig_dict=fig_base,
                        selection=selection,
                        sky_radius_deg=sky_radius_deg,
                        theme_colors=theme_colors,
                        active_time_myr=active_time_now,
                    )
                    _copy_trace_uids(fig_out, current_figure)
                    _copy_trace_visibility(fig_out, current_figure)
                    _copy_slider_active(fig_out, current_figure)
                    _copy_updatemenus_active(fig_out, current_figure)
                    _drop_scene_camera(fig_out)
                    _apply_scene_range_lock(fig_out, scene_range_lock)
                    _ensure_uirevision(fig_out)
            else:
                fig_out = no_update

            if selection_changed:
                if selection is None:
                    sky_src = _build_empty_sky_srcdoc(theme_colors)
                    sky_text = _sky_readout_text(None, sky_radius_deg=sky_radius_deg, sky_survey=sky_survey)
                else:
                    payload = _build_aladin_catalog_payload(
                        selection=selection,
                        sky_radius_deg=sky_radius_deg,
                        centers_by_name=centers_by_name,
                        members_by_cluster=members_by_cluster,
                    )
                    sky_src = _build_aladin_srcdoc(
                        selection=selection,
                        sky_radius_deg=sky_radius_deg,
                        sky_survey=sky_survey,
                        sky_frame=sky_frame,
                        theme_colors=theme_colors,
                        catalog_payload=payload,
                    )
                    sky_text = _sky_readout_text(selection, sky_radius_deg=sky_radius_deg, sky_survey=sky_survey)
            else:
                sky_src = no_update
                sky_text = no_update

            selection_out = selection if selection_changed else no_update
            age_store_out = age_now if age_changed else no_update

            return (fig_out, sky_src, sky_text, selection_out, age_store_out)

    elif age_filter_enabled:
        @app.callback(
            Output(graph_id, "figure"),
            Input("oviz-age-range", "value"),
            State("oviz-base-figure", "data"),
            State(graph_id, "figure"),
            prevent_initial_call=False,
        )
        def _update_age_filtered_figure(age_range, base_figure, current_figure):
            if not age_range or len(age_range) != 2:
                fig_out = copy.deepcopy(base_figure)
            else:
                fig_out = _filter_figure_by_age(base_figure, age_range[0], age_range[1])
            _copy_trace_uids(fig_out, current_figure)
            _copy_trace_visibility(fig_out, current_figure)
            _copy_slider_active(fig_out, current_figure)
            _copy_updatemenus_active(fig_out, current_figure)
            _drop_scene_camera(fig_out)
            _apply_scene_range_lock(fig_out, scene_range_lock)
            _ensure_uirevision(fig_out)
            return fig_out

    elif enable_sky_panel:
        @app.callback(
            Output(graph_id, "figure"),
            Output("oviz-sky-frame", "srcDoc"),
            Output("oviz-sky-readout", "children"),
            Output("oviz-selection-store", "data"),
            Input(graph_id, "clickData"),
            State("oviz-base-figure", "data"),
            State(graph_id, "figure"),
            State("oviz-selection-store", "data"),
            prevent_initial_call=False,
        )
        def _update_sky_only(click_data, base_figure, current_figure, selection_store):
            selection_prev = selection_store if isinstance(selection_store, dict) else None
            selection = selection_prev
            click_selection = _extract_selected_point(click_data)
            if _is_t0_click(click_selection):
                selection = click_selection

            selection_changed = not _selection_equal(selection_prev, selection)
            figure_changed = selection_changed

            if not figure_changed:
                return (no_update, no_update, no_update, no_update)

            fig_out = _build_footprint_patch(
                current_figure=current_figure,
                selection=selection,
                sky_radius_deg=sky_radius_deg,
                theme_colors=theme_colors,
            )
            if fig_out is None:
                active_time_now = _active_slider_time(current_figure)
                fig_out = _apply_footprint_to_figure(
                    fig_dict=base_figure,
                    selection=selection,
                    sky_radius_deg=sky_radius_deg,
                    theme_colors=theme_colors,
                    active_time_myr=active_time_now,
                )
                _copy_trace_uids(fig_out, current_figure)
                _copy_trace_visibility(fig_out, current_figure)
                _copy_slider_active(fig_out, current_figure)
                _copy_updatemenus_active(fig_out, current_figure)
                _drop_scene_camera(fig_out)
                _apply_scene_range_lock(fig_out, scene_range_lock)
                _ensure_uirevision(fig_out)

            if selection is None:
                sky_src = _build_empty_sky_srcdoc(theme_colors)
                sky_text = _sky_readout_text(None, sky_radius_deg=sky_radius_deg, sky_survey=sky_survey)
            elif selection_changed:
                payload = _build_aladin_catalog_payload(
                    selection=selection,
                    sky_radius_deg=sky_radius_deg,
                    centers_by_name=centers_by_name,
                    members_by_cluster=members_by_cluster,
                )
                sky_src = _build_aladin_srcdoc(
                    selection=selection,
                    sky_radius_deg=sky_radius_deg,
                    sky_survey=sky_survey,
                    sky_frame=sky_frame,
                    theme_colors=theme_colors,
                    catalog_payload=payload,
                )
                sky_text = _sky_readout_text(selection, sky_radius_deg=sky_radius_deg, sky_survey=sky_survey)
            else:
                sky_src = no_update
                sky_text = no_update

            return (fig_out, sky_src, sky_text, selection if selection_changed else no_update)

    return app


def run_dash_app(
    figure: go.Figure,
    host: str = "127.0.0.1",
    port: int = 8050,
    debug: bool = False,
    title: str = "oviz",
    **create_dash_kwargs,
) -> Dash:
    """Run the figure as a standard Dash app and return the app instance."""
    app = create_dash_app(figure=figure, title=title, **create_dash_kwargs)
    app.run(host=host, port=port, debug=debug)
    return app


def run_dash_app_in_notebook(
    figure: go.Figure,
    mode: str = "inline",
    host: str = "127.0.0.1",
    port: int = 8050,
    debug: bool = False,
    title: str = "oviz",
    height: int = 750,
    **create_dash_kwargs,
) -> Dash:
    """Run a Dash app from a notebook."""
    app = create_dash_app(figure=figure, title=title, **create_dash_kwargs)

    app.run(
        host=host,
        port=port,
        debug=debug,
        jupyter_mode=mode,
        jupyter_height=height,
    )
    return app


def launch_from_animate3d(
    animate3d,
    *,
    time,
    make_plot_kwargs: Optional[dict] = None,
    dash_kwargs: Optional[dict] = None,
    notebook: bool = True,
    mode: str = "inline",
    host: str = "127.0.0.1",
    port: int = 8050,
    debug: bool = False,
    title: str = "oviz",
    height: int = 750,
):
    """Build a figure via Animate3D.make_plot and launch a Dash app."""
    kwargs = dict(make_plot_kwargs or {})
    kwargs.setdefault("show", False)
    kwargs.setdefault("save_name", None)
    app_kwargs = dict(dash_kwargs or {})

    fig = animate3d.make_plot(time=time, **kwargs)
    if notebook:
        return run_dash_app_in_notebook(
            figure=fig,
            mode=mode,
            host=host,
            port=port,
            debug=debug,
            title=title,
            height=height,
            **app_kwargs,
        )

    return run_dash_app(
        figure=fig,
        host=host,
        port=port,
        debug=debug,
        title=title,
        **app_kwargs,
    )
