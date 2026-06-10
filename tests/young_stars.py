#!/usr/bin/env python3
"""Render a young-cluster member-star version of the main figure."""

from __future__ import annotations

import argparse
import math
import os
import re
import shutil
import sys
from pathlib import Path

TESTS_DIR = Path(__file__).resolve().parent
if str(TESTS_DIR) not in sys.path:
    sys.path.insert(0, str(TESTS_DIR))

from main_figure import (
    DESKTOP_ROOT,
    GALACTIC_PLANE_IMAGE_PATH,
    HOME_DIR,
    MAIN_FIGURE_DUST_MAX_RESOLUTION,
    MAIN_FIGURE_DUST_MAX_RESOLUTION_CAP,
    MAIN_FIGURE_DUST_SAMPLES,
    NOTEBOOK_PATH,
    REPO_ROOT,
    build_mccallum_ne_volume_source_block,
    convert_notebook_to_script_source,
    patch_edenhofer_volume_integer_setting,
    run_script_source,
)


CHRONOS_CLUSTER_RESULTS_PATH = (
    HOME_DIR
    / "Desktop"
    / "astro_research"
    / "supernovae_map"
    / "outputs"
    / "chronos"
    / "dual_model_refit_with_masses"
    / "cluster_results.csv"
)
YOUNG_STARS_MEMBERS_PATH = HOME_DIR / "Downloads" / "members-2.csv"
DEFAULT_OUTPUT_HTML = Path(__file__).resolve().with_suffix(".html")
WEBSITE_OUTPUT_HTML = (
    HOME_DIR
    / "Desktop"
    / "astro_research"
    / "cam_website"
    / "interactive_figures"
    / "young_stars.html"
)
YOUNG_STARS_MAX_AGE_MYR = 40.0
YOUNG_STARS_RED_MAX_AGE_MYR = 20.0
YOUNG_STARS_XY_BOUNDS_PC = (-1000.0, 1000.0)
YOUNG_STARS_Z_BOUNDS_PC = (-400.0, 400.0)
YOUNG_STARS_POINT_SIZE_REFERENCE_SPAN_PC = 20000.0
YOUNG_STARS_XYZ_RANGES = (
    YOUNG_STARS_XY_BOUNDS_PC,
    YOUNG_STARS_XY_BOUNDS_PC,
    YOUNG_STARS_Z_BOUNDS_PC,
)
YOUNG_STARS_VOLUME_CLIP_BOUNDS = {
    "x": [float(YOUNG_STARS_XY_BOUNDS_PC[0]), float(YOUNG_STARS_XY_BOUNDS_PC[1])],
    "y": [float(YOUNG_STARS_XY_BOUNDS_PC[0]), float(YOUNG_STARS_XY_BOUNDS_PC[1])],
    "z": [float(YOUNG_STARS_Z_BOUNDS_PC[0]), float(YOUNG_STARS_Z_BOUNDS_PC[1])],
}
YOUNG_STARS_XY_SPAN_PC = YOUNG_STARS_XY_BOUNDS_PC[1] - YOUNG_STARS_XY_BOUNDS_PC[0]
YOUNG_STARS_Z_SPAN_PC = YOUNG_STARS_Z_BOUNDS_PC[1] - YOUNG_STARS_Z_BOUNDS_PC[0]
YOUNG_STARS_ASPECT_Z = YOUNG_STARS_Z_SPAN_PC / YOUNG_STARS_XY_SPAN_PC


def build_young_stars_axis_style_source_block() -> str:
    return """
for _axis_name, _axis_title in (
    ("xaxis", "x (pc)"),
    ("yaxis", "y (pc)"),
    ("zaxis", "z (pc)"),
):
    _axis_layout = plot_3d.figure_layout_dict["scene"][_axis_name]
    _axis_layout.update(
        title=_axis_title,
        visible=True,
        showline=True,
        showgrid=False,
        zeroline=False,
        linecolor="#777777",
        linewidth=1.0,
        tickfont={"size": 14, "color": "#777777", "family": "Helvetica"},
        title_font={"size": 18, "color": "#777777", "family": "Helvetica"},
        dtick=200,
        nticks=6,
    )
""".strip()


def build_young_stars_scene_postprocess_source_block() -> str:
    return f"""
import math as _young_stars_math

YOUNG_STARS_SCENE_BOUNDS = {{
    "x": {tuple(float(v) for v in YOUNG_STARS_XY_BOUNDS_PC)!r},
    "y": {tuple(float(v) for v in YOUNG_STARS_XY_BOUNDS_PC)!r},
    "z": {tuple(float(v) for v in YOUNG_STARS_Z_BOUNDS_PC)!r},
}}
YOUNG_STARS_SCENE_MAX_SPAN_PC = {max(YOUNG_STARS_XY_SPAN_PC, YOUNG_STARS_Z_SPAN_PC)!r}
YOUNG_STARS_SCENE_ASPECT_Z = {float(YOUNG_STARS_ASPECT_Z)!r}
YOUNG_STARS_POINT_SIZE_REFERENCE_SPAN_PC = {float(YOUNG_STARS_POINT_SIZE_REFERENCE_SPAN_PC)!r}


def _young_stars_point_inside_scene_bounds(point):
    try:
        coords = [float(point[axis]) for axis in ("x", "y", "z")]
    except Exception:
        return False
    if not all(_young_stars_math.isfinite(value) for value in coords):
        return False
    return all(
        YOUNG_STARS_SCENE_BOUNDS[axis][0] <= coords[idx] <= YOUNG_STARS_SCENE_BOUNDS[axis][1]
        for idx, axis in enumerate(("x", "y", "z"))
    )


def _young_stars_segment_inside_scene_bounds(segment):
    if not isinstance(segment, (list, tuple)) or len(segment) < 6:
        return False
    try:
        start = {{"x": float(segment[0]), "y": float(segment[1]), "z": float(segment[2])}}
        end = {{"x": float(segment[3]), "y": float(segment[4]), "z": float(segment[5])}}
    except Exception:
        return False
    return (
        _young_stars_point_inside_scene_bounds(start)
        and _young_stars_point_inside_scene_bounds(end)
    )


def _young_stars_reduced_focus_point(point):
    try:
        x = float(point["x"])
        y = float(point["y"])
        z = float(point["z"])
    except Exception:
        return None
    if not all(_young_stars_math.isfinite(value) for value in (x, y, z)):
        return None
    if not _young_stars_point_inside_scene_bounds({{"x": x, "y": y, "z": z}}):
        return None

    focus_point = {{"x": x, "y": y, "z": z}}
    selection = point.get("selection") if isinstance(point, dict) else None
    if isinstance(selection, dict):
        focus_selection = {{}}
        for key in ("cluster_name", "trace_name"):
            value = str(selection.get(key, "")).strip()
            if value:
                focus_selection[key] = value
        for key in ("x0", "y0", "z0", "ra_deg", "dec_deg", "age_now_myr", "n_stars"):
            value = selection.get(key)
            try:
                numeric = float(value)
            except Exception:
                continue
            if _young_stars_math.isfinite(numeric):
                focus_selection[key] = numeric
        if focus_selection:
            focus_point["selection"] = focus_selection
    point_n_stars = point.get("n_stars") if isinstance(point, dict) else None
    try:
        point_n_stars = float(point_n_stars)
    except Exception:
        point_n_stars = float("nan")
    if _young_stars_math.isfinite(point_n_stars):
        focus_point["n_stars"] = point_n_stars
    return focus_point


def _finalize_young_stars_scene_spec(scene_spec):
    present_frames = [
        frame
        for frame in scene_spec.get("frames", []) or []
        if _young_stars_math.isclose(float(frame.get("time", float("nan"))), 0.0, abs_tol=1e-9)
    ]
    if present_frames:
        scene_spec["frames"] = present_frames[:1]
        scene_spec["initial_frame_index"] = 0
        timeline = scene_spec.setdefault("timeline", {{}})
        timeline["enabled"] = False
        timeline["frame_count"] = 1

    scene_spec["ranges"] = {{axis: list(bounds) for axis, bounds in YOUNG_STARS_SCENE_BOUNDS.items()}}
    scene_spec["center"] = {{"x": 0.0, "y": 0.0, "z": 0.0}}
    scene_spec["max_span"] = float(YOUNG_STARS_SCENE_MAX_SPAN_PC)
    scene_spec["point_size_reference_span_pc"] = float(YOUNG_STARS_POINT_SIZE_REFERENCE_SPAN_PC)

    layout_scene = scene_spec.setdefault("layout", {{}}).setdefault("scene", {{}})
    for axis_name, axis in (("xaxis", "x"), ("yaxis", "y"), ("zaxis", "z")):
        axis_layout = layout_scene.setdefault(axis_name, {{}})
        axis_layout.update(
            range=list(YOUNG_STARS_SCENE_BOUNDS[axis]),
            visible=True,
            showline=True,
            showgrid=False,
            zeroline=False,
            linecolor="#777777",
            linewidth=1.0,
        )
    layout_scene["aspectratio"] = {{"x": 1.0, "y": 1.0, "z": YOUNG_STARS_SCENE_ASPECT_Z}}

    for axis, axis_spec in (scene_spec.get("axes") or {{}}).items():
        axis_spec.update(
            range=list(YOUNG_STARS_SCENE_BOUNDS.get(axis, (-500.0, 500.0))),
            visible=True,
            showline=True,
            showgrid=False,
            zeroline=False,
            linecolor="#777777",
            linewidth=1.0,
        )

    initial_state = scene_spec.setdefault("initial_state", {{}})
    initial_state["clean_box_axes"] = True
    initial_state["galaxy_image"] = False
    for _galaxy_key in (
        "galaxy_image_path",
        "galaxy_image_size_pc",
        "galaxy_image_opacity",
        "galaxy_image_hide_below_scale_bar_pc",
        "galaxy_image_fade_start_scale_bar_pc",
    ):
        initial_state.pop(_galaxy_key, None)
    initial_state["initial_zoom_anchor"] = {{"x": 0.0, "y": 0.0, "z": 0.0}}
    initial_state["camera"] = {{
        "position": {{"x": 140.0, "y": -1310.0, "z": 940.0}},
        "target": {{"x": 0.0, "y": 0.0, "z": 0.0}},
        "up": {{"x": 0.0, "y": 0.0, "z": 1.0}},
    }}
    global_controls = initial_state.setdefault("global_controls", {{}})
    global_controls["camera_fov"] = 55.0
    global_controls["point_size_scale"] = 0.75
    global_controls["point_glow_strength"] = 0.0
    global_controls["fade_in_time_myr"] = 8.0

    age_kde = scene_spec.get("age_kde")
    if isinstance(age_kde, dict):
        age_kde["traces"] = [
            trace
            for trace in age_kde.get("traces", []) or []
            if str(trace.get("trace_name", "")).strip() != "Sun"
        ]
        age_kde["cluster_points"] = [
            point
            for point in age_kde.get("cluster_points", []) or []
            if str(point.get("trace_name", "")).strip() != "Sun"
        ]

    for frame in scene_spec.get("frames", []) or []:
        for trace in frame.get("traces", []) or []:
            if trace.get("points"):
                for point in trace["points"]:
                    point.pop("motion", None)
                trace["focus_points"] = [
                    focus_point
                    for focus_point in (
                        _young_stars_reduced_focus_point(point) for point in trace["points"]
                    )
                    if focus_point is not None
                ]
                trace["points"] = [
                    point for point in trace["points"] if _young_stars_point_inside_scene_bounds(point)
                ]
            elif trace.get("focus_points"):
                trace["focus_points"] = [
                    point
                    for point in trace["focus_points"]
                    if _young_stars_point_inside_scene_bounds(point)
                ]
            if trace.get("labels"):
                trace["labels"] = [
                    label for label in trace["labels"] if _young_stars_point_inside_scene_bounds(label)
                ]
            if trace.get("segments"):
                trace["segments"] = [
                    segment for segment in trace["segments"] if _young_stars_segment_inside_scene_bounds(segment)
                ]
""".strip()


def patch_edenhofer_volume_clip_xyz_bounds(source: str) -> str:
    edenhofer_block_pattern = (
        r'(?ms)(edenhofer_volume\s*=\s*\{.*?'
        r'"name"\s*:\s*"Edenhofer\+2024 Dust".*?'
        r'^\s*\})'
    )
    block_match = re.search(edenhofer_block_pattern, source)
    if not block_match:
        raise RuntimeError("Could not find the Edenhofer dust volume block.")

    block = block_match.group(1)
    clip_text = f'"clip_bounds": {YOUNG_STARS_VOLUME_CLIP_BOUNDS!r},'
    updated_block, replaced = re.subn(
        r'(?m)^(\s*)"clip_bounds"\s*:\s*\{.*?\}\s*,\s*$',
        rf"\1{clip_text}",
        block,
        count=1,
    )
    if replaced == 0:
        insert_pattern = r'(?m)^(\s*)"hdu"\s*:\s*["\'][^"\']+["\']\s*,\s*$'

        def insert_clip_bounds(match):
            indent = match.group(1)
            return f"{match.group(0)}\n{indent}{clip_text}"

        updated_block, replaced = re.subn(
            insert_pattern,
            insert_clip_bounds,
            block,
            count=1,
        )
    if replaced != 1:
        raise RuntimeError("Could not update the Edenhofer dust clip_bounds.")

    return source[:block_match.start(1)] + updated_block + source[block_match.end(1):]


def build_cropped_mccallum_ne_volume_source_block() -> str:
    block = build_mccallum_ne_volume_source_block()
    block, replaced = re.subn(
        r'"clip_bounds"\s*:\s*\{[^}]*\}',
        f'"clip_bounds": {YOUNG_STARS_VOLUME_CLIP_BOUNDS!r}',
        block,
        count=1,
    )
    if replaced != 1:
        raise RuntimeError("Could not crop the McCallum electron-density clip_bounds.")
    return block


def build_young_star_layer_source_block() -> str:
    return f"""
import os
from pathlib import Path as _OvizPath

YOUNG_STARS_MEMBERS_PATH = _OvizPath(
    os.environ.get("OVIZ_YOUNG_STARS_MEMBERS_CSV", {str(YOUNG_STARS_MEMBERS_PATH)!r})
).expanduser()
CHRONOS_CLUSTER_RESULTS_PATH = _OvizPath(
    os.environ.get("OVIZ_CHRONOS_CLUSTER_RESULTS_CSV", {str(CHRONOS_CLUSTER_RESULTS_PATH)!r})
).expanduser()
YOUNG_STARS_MAX_AGE_MYR = {float(YOUNG_STARS_MAX_AGE_MYR)!r}
YOUNG_STARS_RED_MAX_AGE_MYR = {float(YOUNG_STARS_RED_MAX_AGE_MYR)!r}
YOUNG_STARS_XY_BOUNDS_PC = {tuple(float(v) for v in YOUNG_STARS_XY_BOUNDS_PC)!r}
YOUNG_STARS_Z_BOUNDS_PC = {tuple(float(v) for v in YOUNG_STARS_Z_BOUNDS_PC)!r}


def _young_stars_bool_series(series):
    return series.astype(str).str.strip().str.lower().isin(["true", "1", "t", "yes", "y"])


def _load_young_cluster_member_stars():
    if not CHRONOS_CLUSTER_RESULTS_PATH.exists():
        raise FileNotFoundError(f"Missing Chronos cluster results: {{CHRONOS_CLUSTER_RESULTS_PATH}}")
    if not YOUNG_STARS_MEMBERS_PATH.exists():
        raise FileNotFoundError(f"Missing cluster member catalog: {{YOUNG_STARS_MEMBERS_PATH}}")

    chronos_ages = pd.read_csv(
        CHRONOS_CLUSTER_RESULTS_PATH,
        usecols=[
            "name",
            "parsec_age_lo",
            "parsec_age_mode",
            "parsec_age_hi",
            "parsec_status",
        ],
    ).rename(
        columns={{
            "parsec_age_lo": "chronos_age_lo_myr",
            "parsec_age_mode": "chronos_age_myr",
            "parsec_age_hi": "chronos_age_hi_myr",
        }}
    )
    for _age_col in ["chronos_age_lo_myr", "chronos_age_myr", "chronos_age_hi_myr"]:
        chronos_ages[_age_col] = pd.to_numeric(chronos_ages[_age_col], errors="coerce")
    chronos_ages = chronos_ages.loc[
        chronos_ages["parsec_status"].eq("success")
        & chronos_ages["chronos_age_myr"].notnull()
    ].copy()

    cluster_positions = df_hunt_full[["name", "x", "y", "z", "class_50", "cst"]].copy()
    for _cluster_col in ["x", "y", "z", "class_50", "cst"]:
        cluster_positions[_cluster_col] = pd.to_numeric(cluster_positions[_cluster_col], errors="coerce")
    young_clusters = cluster_positions.merge(
        chronos_ages,
        on="name",
        how="inner",
        validate="m:1",
    )
    young_clusters = young_clusters.loc[
        (young_clusters["chronos_age_myr"] < YOUNG_STARS_MAX_AGE_MYR)
        & (young_clusters["x"].between(*YOUNG_STARS_XY_BOUNDS_PC))
        & (young_clusters["y"].between(*YOUNG_STARS_XY_BOUNDS_PC))
        & (young_clusters["z"].between(*YOUNG_STARS_Z_BOUNDS_PC))
        & (young_clusters["class_50"] > 0.5)
        & (young_clusters["cst"] > 5.0)
    ].copy()
    if young_clusters.empty:
        raise RuntimeError("No young Chronos clusters passed the age, distance, and quality cuts.")

    young_cluster_names = set(young_clusters["name"].astype(str))
    age_by_cluster = young_clusters.set_index("name")["chronos_age_myr"].to_dict()
    age_lo_by_cluster = young_clusters.set_index("name")["chronos_age_lo_myr"].to_dict()
    age_hi_by_cluster = young_clusters.set_index("name")["chronos_age_hi_myr"].to_dict()

    member_chunks = []
    member_usecols = [
        "name",
        "source_id",
        "within_r_t",
        "probability",
        "x",
        "y",
        "z",
        "phot_g_mean_mag",
    ]
    for member_chunk in pd.read_csv(
        YOUNG_STARS_MEMBERS_PATH,
        usecols=member_usecols,
        chunksize=200_000,
    ):
        member_chunk["name"] = member_chunk["name"].astype(str)
        member_chunk = member_chunk.loc[member_chunk["name"].isin(young_cluster_names)].copy()
        if member_chunk.empty:
            continue
        member_chunk = member_chunk.loc[_young_stars_bool_series(member_chunk["within_r_t"])].copy()
        if member_chunk.empty:
            continue
        for _member_col in ["x", "y", "z", "probability", "phot_g_mean_mag"]:
            member_chunk[_member_col] = pd.to_numeric(member_chunk[_member_col], errors="coerce")
        member_chunk = member_chunk.loc[
            np.isfinite(member_chunk[["x", "y", "z"]].to_numpy(dtype=float)).all(axis=1)
            & member_chunk["x"].between(*YOUNG_STARS_XY_BOUNDS_PC)
            & member_chunk["y"].between(*YOUNG_STARS_XY_BOUNDS_PC)
            & member_chunk["z"].between(*YOUNG_STARS_Z_BOUNDS_PC)
        ].copy()
        if not member_chunk.empty:
            member_chunks.append(member_chunk)

    if not member_chunks:
        raise RuntimeError("No member stars matched the selected young Chronos clusters.")

    members = pd.concat(member_chunks, ignore_index=True)
    members = members.sort_values(["source_id", "probability"], ascending=[True, False])
    members = members.drop_duplicates(subset=["source_id"], keep="first").reset_index(drop=True)
    members["cluster_name"] = members["name"].astype(str)
    members["n_stars"] = (
        members.groupby("cluster_name")["source_id"].transform("count").astype(int)
    )
    members["age_myr"] = members["cluster_name"].map(age_by_cluster)
    members["chronos_age_lo_myr"] = members["cluster_name"].map(age_lo_by_cluster)
    members["chronos_age_hi_myr"] = members["cluster_name"].map(age_hi_by_cluster)
    members["source_id"] = members["source_id"].astype(str)
    members["name"] = members["cluster_name"]

    print(
        "Loaded "
        f"{{len(members)}} member stars from {{len(young_clusters)}} Chronos clusters "
        f"younger than {{YOUNG_STARS_MAX_AGE_MYR:g}} Myr in the local "
        "2 x 2 kpc, |z| <= 400 pc box."
    )
    return members


class YoungStarLayer(Layer):
    def set_age_based_sizes(self, fade_in_time=None, fade_in_and_out=None, fade_in_and_disp=None, disp_time=None):
        if self.df_int is None or self.df_int.empty:
            return

        actual_fade_in_time = fade_in_time if fade_in_time is not None else self.fade_in_time
        if actual_fade_in_time is None:
            actual_fade_in_time = 5.0
        actual_fade_in_time = max(float(actual_fade_in_time), 1e-6)

        age_at_t = (
            pd.to_numeric(self.df_int["age_myr"], errors="coerce")
            + pd.to_numeric(self.df_int["time"], errors="coerce")
        )
        progress = np.clip((age_at_t + actual_fade_in_time) / actual_fade_in_time, 0.0, 1.0)
        # Match the main figure's soft birth transition while supporting many
        # member stars with the same cluster name.
        smooth_progress = 1.0 / (1.0 + np.exp(-(progress - 0.5) / 0.18))
        smooth_progress = np.where(age_at_t >= 0.0, 1.0, smooth_progress)
        self.df_int["size"] = self.min_size + (self.max_size - self.min_size) * smooth_progress
        self.sizes_set = True


def _build_young_star_trace(members, *, layer_name, color, opacity):
    if members.empty:
        print(f"Skipping empty {{layer_name}} trace.")
        return None
    return YoungStarLayer(
        members.copy(),
        layer_name=layer_name,
        min_size=0.0,
        max_size=2.0,
        color=color,
        opacity=opacity,
        marker_style="circle",
        show_tracks=False,
        size_by_n_stars=False,
        assume_stationary=True,
    )


young_star_members = _load_young_cluster_member_stars()
_young_star_ages = pd.to_numeric(young_star_members["age_myr"], errors="coerce")
young_star_members_0_20 = young_star_members.loc[
    _young_star_ages < YOUNG_STARS_RED_MAX_AGE_MYR
].copy()
young_star_members_20_40 = young_star_members.loc[
    (_young_star_ages >= YOUNG_STARS_RED_MAX_AGE_MYR)
    & (_young_star_ages < YOUNG_STARS_MAX_AGE_MYR)
].copy()
print(
    "Age-bin split: "
    f"{{len(young_star_members_0_20)}} stars younger than {{YOUNG_STARS_RED_MAX_AGE_MYR:g}} Myr; "
    f"{{len(young_star_members_20_40)}} stars from {{YOUNG_STARS_RED_MAX_AGE_MYR:g}}-{{YOUNG_STARS_MAX_AGE_MYR:g}} Myr."
)
young_star_trace_0_20 = _build_young_star_trace(
    young_star_members_0_20,
    layer_name="Young Cluster Stars (0-20 Myr)",
    color="red",
    opacity=0.9,
)
young_star_trace_20_40 = _build_young_star_trace(
    young_star_members_20_40,
    layer_name="Young Cluster Stars (20-40 Myr)",
    color="#00ffff",
    opacity=0.9,
)
young_star_traces = [
    trace for trace in (young_star_trace_0_20, young_star_trace_20_40) if trace is not None
]
""".strip()


def build_quintana_association_star_layer_source_block() -> str:
    return """
QUINTANA_ASSOCIATION_MEMBERS_CSV = os.environ.get("OVIZ_QUINTANA_ASSOCIATION_MEMBERS_CSV", "").strip()
QUINTANA_CYAN_LAYER_NAME = "Quintana Association Stars"


def _quintana_first_column(columns, *candidates):
    column_by_lower = {{str(column).strip().lower(): column for column in columns}}
    for candidate in candidates:
        column = column_by_lower.get(str(candidate).strip().lower())
        if column is not None:
            return column
    return None


def _load_quintana_association_stars():
    if not QUINTANA_ASSOCIATION_MEMBERS_CSV:
        print(
            "Quintana association-member table is not configured; "
            "skipping the cyan association-star layer instead of plotting non-members."
        )
        return None

    catalog_path = _OvizPath(QUINTANA_ASSOCIATION_MEMBERS_CSV).expanduser()

    if not catalog_path.exists():
        raise FileNotFoundError(f"Missing Quintana/Q25 star catalog: {{catalog_path}}")

    stars = pd.read_csv(catalog_path)
    x_col = _quintana_first_column(stars.columns, "x", "X", "x_pc", "X_pc")
    y_col = _quintana_first_column(stars.columns, "y", "Y", "y_pc", "Y_pc")
    z_col = _quintana_first_column(stars.columns, "z", "Z", "z_pc", "Z_pc")
    if x_col is None or y_col is None or z_col is None:
        raise RuntimeError(
            "Quintana/Q25 catalog must contain heliocentric Cartesian x, y, z columns."
        )

    assoc_col = _quintana_first_column(
        stars.columns,
        "association",
        "association_name",
        "ob_association",
        "obassociation",
        "assoc",
        "group",
        "cluster",
        "name",
        "Name",
    )
    source_id_col = _quintana_first_column(
        stars.columns,
        "source_id",
        "sourceid",
        "gaia",
        "Gaia",
        "gaia_dr3",
        "GaiaDR3",
        "gaiaedr3",
        "GaiaEDR3",
    )

    out = pd.DataFrame({{
        "x": pd.to_numeric(stars[x_col], errors="coerce"),
        "y": pd.to_numeric(stars[y_col], errors="coerce"),
        "z": pd.to_numeric(stars[z_col], errors="coerce"),
        "age_myr": np.nan,
    }})
    if assoc_col is not None:
        out["name"] = stars[assoc_col].astype(str).str.strip()
        out.loc[out["name"].eq("") | out["name"].str.lower().eq("nan"), "name"] = (
            "Quintana association"
        )
    else:
        out["name"] = QUINTANA_CYAN_LAYER_NAME
    if source_id_col is not None:
        out["source_id"] = stars[source_id_col].astype(str).str.strip()

    out = out.loc[
        np.isfinite(out[["x", "y", "z"]].to_numpy(dtype=float)).all(axis=1)
        & out["x"].between(*YOUNG_STARS_XY_BOUNDS_PC)
        & out["y"].between(*YOUNG_STARS_XY_BOUNDS_PC)
        & out["z"].between(*YOUNG_STARS_Z_BOUNDS_PC)
    ].copy()
    if source_id_col is not None:
        out = out.sort_values("source_id").drop_duplicates("source_id", keep="first")
    out = out.reset_index(drop=True)

    if out.empty:
        print(
            "Loaded 0 Quintana/Q25 association-member stars in the local "
            "2 x 2 kpc, |z| <= 400 pc box; skipping the cyan layer."
        )
        return None

    print(f"Loaded {{len(out)}} stars for Quintana/Q25 association cyan layer from {{catalog_path}}.")
    return out


class StaticAssociationStarLayer(Layer):
    def set_age_based_sizes(self, fade_in_time=None, fade_in_and_out=None, fade_in_and_disp=None, disp_time=None):
        if self.df_int is None or self.df_int.empty:
            return
        self.df_int["size"] = float(self.max_size)
        self.sizes_set = True


quintana_association_stars = _load_quintana_association_stars()
if quintana_association_stars is not None:
    quintana_association_trace = StaticAssociationStarLayer(
        quintana_association_stars,
        layer_name=QUINTANA_CYAN_LAYER_NAME,
        min_size=2.0,
        max_size=2.0,
        color="#00ffff",
        opacity=0.65,
        marker_style="circle",
        show_tracks=False,
        size_by_n_stars=False,
        assume_stationary=True,
        default_age_myr=float("nan"),
    )
else:
    quintana_association_trace = None
""".strip()


def patch_young_stars_script_source(
    source: str,
    output_html: Path,
    *,
    compact_payload: bool = True,
) -> str:
    source = source.replace("/Users/cam", str(HOME_DIR))
    source = source.replace("/Users/cam/Desktop", str(DESKTOP_ROOT))
    source = source.replace(
        "from oviz import Trace, TraceCollection, Animate3D",
        "from oviz import Trace, TraceCollection, Animate3D, Layer",
    )

    source, replaced_time_grid = re.subn(
        r"(?m)^time_int\s*=\s*np\.round\(np\.arange\(0,\s*-66,\s*-1\),\s*1\)\s*$",
        f"time_int = np.array([0.0, -{float(YOUNG_STARS_MAX_AGE_MYR)!r}])",
        source,
        count=1,
    )
    if replaced_time_grid != 1:
        raise RuntimeError("Could not replace the young-star time grid with a present-day-only frame.")

    static_sun_block = """
class StaticPointLayer(Layer):
    def set_age_based_sizes(self, fade_in_time=None, fade_in_and_out=None, fade_in_and_disp=None, disp_time=None):
        if self.df_int is None or self.df_int.empty:
            return
        self.df_int["size"] = float(self.max_size)
        self.sizes_set = True


sun_trace = StaticPointLayer(
    sun[["x", "y", "z"]].copy(),
    layer_name="Sun",
    min_size=5,
    max_size=5,
    color="yellow",
    opacity=1,
    marker_style="circle",
    show_tracks=False,
    size_by_n_stars=False,
    assume_stationary=True,
    default_age_myr=float("nan"),
)
""".strip()
    young_star_block = build_young_star_layer_source_block()
    quintana_association_block = build_quintana_association_star_layer_source_block()
    source, inserted_young_stars = re.subn(
        r"(?m)^(sun_trace\s*=\s*Trace\(sun,.*\)\s*)$",
        static_sun_block + "\n\n" + young_star_block + "\n\n" + quintana_association_block,
        source,
        count=1,
    )
    if inserted_young_stars != 1:
        raise RuntimeError("Could not inject the young-cluster member-star layer.")

    traces_block = """
_young_stars_trace_list = [
    sun_trace,
    *young_star_traces,
]
_young_stars_group_names = ["Sun", *[trace.layer_name for trace in young_star_traces]]
if quintana_association_trace is not None:
    _young_stars_trace_list.append(quintana_association_trace)
    _young_stars_group_names.append(quintana_association_trace.layer_name)

trace_groupings = {
    "Young Stars": _young_stars_group_names,
}

traces = TraceCollection(_young_stars_trace_list)

### NEW THREEJS
""".lstrip()
    source, replaced_traces = re.subn(
        r"(?ms)^traces\s*=\s*TraceCollection\(\[.*?^### NEW THREEJS\s*$",
        traces_block,
        source,
        count=1,
    )
    if replaced_traces != 1:
        raise RuntimeError("Could not replace the main trace collection with the young-star layer.")

    source, replaced_xyz_ranges = re.subn(
        r"(?m)^xyz_widths\s*=\s*\([^)]+\)\s*$",
        f"xyz_widths = (1000, 1000, 400)\nxyz_ranges = {YOUNG_STARS_XYZ_RANGES!r}",
        source,
        count=1,
    )
    if replaced_xyz_ranges != 1:
        raise RuntimeError("Could not set the young-star local xyz ranges.")

    source, inserted_xyz_ranges = re.subn(
        r"(?m)^(\s*)xyz_widths\s*=\s*xyz_widths\s*,\s*$",
        r"\1xyz_widths = xyz_widths,\n\1xyz_ranges = xyz_ranges,",
        source,
        count=1,
    )
    if inserted_xyz_ranges != 1:
        raise RuntimeError("Could not pass the young-star xyz ranges to Animate3D.")

    axis_style_block = build_young_stars_axis_style_source_block()
    source, inserted_axis_style = re.subn(
        r"(?m)^(save_name\s*=)",
        axis_style_block + "\n\n" + r"\1",
        source,
        count=1,
    )
    if inserted_axis_style != 1:
        raise RuntimeError("Could not inject the young-star box axis style.")

    source, replaced_age_kde = re.subn(
        r"(?m)^(\s*)show_age_kde_inset\s*=\s*(?:False|True)\s*,\s*$",
        r"\1show_age_kde_inset=True,",
        source,
        count=1,
    )
    if replaced_age_kde != 1:
        raise RuntimeError("Could not enable the age KDE inset.")

    source = patch_edenhofer_volume_integer_setting(
        source,
        "max_resolution",
        MAIN_FIGURE_DUST_MAX_RESOLUTION,
    )
    source = patch_edenhofer_volume_integer_setting(
        source,
        "max_resolution_cap",
        MAIN_FIGURE_DUST_MAX_RESOLUTION_CAP,
        insert_after="max_resolution",
    )
    source = patch_edenhofer_volume_integer_setting(
        source,
        "samples",
        MAIN_FIGURE_DUST_SAMPLES,
    )
    source = patch_edenhofer_volume_clip_xyz_bounds(source)

    output_html_str = str(output_html)
    source, replaced_output = re.subn(
        r"(?m)^save_name\s*=\s*['\"][^'\"]+['\"]",
        f"save_name = {output_html_str!r}",
        source,
        count=1,
    )
    if replaced_output != 1:
        raise RuntimeError("Could not rewrite save_name in converted notebook script.")

    if "mccallum_ne_volumes = _build_mccallum_ne_volumes_for_main_figure()" not in source:
        ne_block = build_cropped_mccallum_ne_volume_source_block()
        source, inserted_ne = re.subn(
            r"(?m)^(fig3d\s*=\s*plot_3d\.make_plot\()",
            ne_block + "\n\n" + r"\1",
            source,
            count=1,
        )
        if inserted_ne != 1:
            raise RuntimeError("Could not inject McCallum electron-density volume helper block.")

    source, replaced_volume_list = re.subn(
        r"volumes\s*=\s*\[\s*edenhofer_volume\s*,\s*mccallum_ne\s*\]",
        "volumes=[edenhofer_volume, *mccallum_ne_volumes]",
        source,
        count=1,
    )
    if replaced_volume_list != 1:
        raise RuntimeError("Could not replace the notebook volume list.")

    source, replaced_galactic_mode = re.subn(
        r"(?m)^(\s*)galactic_mode\s*=\s*True\s*,\s*$",
        r"\1galactic_mode=False,",
        source,
        count=1,
    )
    if replaced_galactic_mode != 1:
        raise RuntimeError("Could not disable the galactic-mode 10 kpc range override.")

    source, replaced_gc_line = re.subn(
        r"(?m)^(\s*)show_gc_line\s*=\s*True\s*,?\s*$",
        r"\1show_gc_line=False,",
        source,
        count=1,
    )
    if replaced_gc_line != 1:
        raise RuntimeError("Could not disable the galactic radius guide line.")

    source, replaced_make_plot_save = re.subn(
        r"(?m)^(\s*)save_name\s*=\s*save_name\s*,\s*$",
        r"\1save_name=None,",
        source,
        count=1,
    )
    if replaced_make_plot_save != 1:
        raise RuntimeError("Could not defer writing until after scene post-processing.")

    source, edenhofer_time_patch = re.subn(
        r'("colormap"\s*:\s*greys_cmap,\s*# or just "ice"\s*\n)(\s*})',
        r'\1    "supports_show_all_times": True,\n    "co_rotate_with_frame": True,\n\2',
        source,
        count=1,
    )
    if edenhofer_time_patch != 1:
        raise RuntimeError("Could not patch Edenhofer volume frame behavior.")

    if "threejs_initial_state=" not in source:
        initial_state_bits = [
            "'current_group': 'Young Stars'",
            "'click_selection_enabled': False",
            f"'compact_payload_enabled': {bool(compact_payload)!r}",
            "'clean_box_axes': True",
            "'scene_float_precision': 2",
            "'active_volume_key': 'volume-0'",
            "'legend_state': {'volume-0': True, **({'mccallum-ne': False} if mccallum_ne_volumes else {})}",
            "'volume_state_by_key': {'volume-0': {'visible': True}, **({'mccallum-ne': {'visible': False}} if mccallum_ne_volumes else {})}",
            "'galaxy_image': False",
            f"'galaxy_image_path': {str(GALACTIC_PLANE_IMAGE_PATH)!r}",
            "'galaxy_image_size_pc': 40000.0",
            "'galaxy_image_opacity': 0.35",
            "'galaxy_image_hide_below_scale_bar_pc': 420.0",
            "'galaxy_image_fade_start_scale_bar_pc': 700.0",
            "'sky_dome_enabled': True",
            "'sky_dome_background_mode': 'live_aladin'",
            "'sky_dome_source': 'aladin'",
            "'sky_dome_projection': 'TAN'",
            "'sky_dome_capture_width_px': 4096",
            "'sky_dome_capture_height_px': 2048",
            "'sky_dome_capture_format': 'image/jpeg'",
            "'sky_dome_capture_quality': 0.94",
            "'sky_dome_radius_pc': 40000.0",
            "'sky_dome_opacity': 1.0",
            "'sky_dome_force_visible': False",
            "'sky_dome_full_opacity_scale_bar_pc': 120.0",
            "'sky_dome_fade_out_scale_bar_pc': 360.0",
            "'sky_layers': ["
            "{'key': 'P/DSS2/color', 'label': 'DSS2 Color', 'survey': 'P/DSS2/color', 'opacity': 1.0, 'visible': True}, "
            "{'key': 'P/Mellinger/color', 'label': 'Mellinger Color', 'survey': 'P/Mellinger/color', 'opacity': 1.0, 'visible': False}, "
            "{'key': 'P/PLANCK/R2/HFI/color', 'label': 'Planck Dust Emission Color', 'survey': 'P/PLANCK/R2/HFI/color', 'opacity': 1.0, 'visible': True}"
            "]",
            "'active_sky_layer_key': 'P/DSS2/color'",
            "'initial_zoom_anchor': {'x': 0.0, 'y': 0.0, 'z': 0.0}",
            "'camera': {'position': {'x': 140.0, 'y': -1310.0, 'z': 940.0}, 'target': {'x': 0.0, 'y': 0.0, 'z': 0.0}, 'up': {'x': 0.0, 'y': 0.0, 'z': 1.0}}",
            "'global_controls': {'camera_fov': 55.0, 'point_size_scale': 0.75, 'point_glow_strength': 0.0, 'fade_in_time_myr': 8.0}",
        ]
        source, injected_initial_state = re.subn(
            r'(renderer\s*=\s*"threejs",\s*\n)',
            '\\1    threejs_initial_state={'
            + ", ".join(initial_state_bits)
            + '},\n',
            source,
            count=1,
        )
        if injected_initial_state != 1:
            raise RuntimeError("Could not inject threejs_initial_state into make_plot call.")

    scene_postprocess_block = build_young_stars_scene_postprocess_source_block()
    source, inserted_scene_postprocess = re.subn(
        r"(?m)^(mccallum_ne_volumes\s*=\s*_build_mccallum_ne_volumes_for_main_figure\(\)\s*)$",
        r"\1\n\n" + scene_postprocess_block,
        source,
        count=1,
    )
    if inserted_scene_postprocess != 1:
        raise RuntimeError("Could not inject the young-star scene post-processor.")

    source, inserted_scene_write = re.subn(
        r"(?ms)^(fig3d\s*=\s*plot_3d\.make_plot\(.*?^\)\s*)$",
        r"\1\n_finalize_young_stars_scene_spec(fig3d.scene_spec)\nfig3d.write_html(save_name)\n",
        source,
        count=1,
    )
    if inserted_scene_write != 1:
        raise RuntimeError("Could not write the young-star scene after post-processing.")

    return source


def run_young_stars_figure(
    output_html: Path = DEFAULT_OUTPUT_HTML,
    *,
    compact_payload: bool = True,
    website_output_html: Path | None = None,
) -> Path:
    output_html.parent.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl")
    os.environ.setdefault("XDG_CACHE_HOME", "/tmp")
    os.environ.setdefault("MPLBACKEND", "Agg")

    if not NOTEBOOK_PATH.exists():
        raise FileNotFoundError(f"Missing notebook: {NOTEBOOK_PATH}")
    if not REPO_ROOT.exists():
        raise FileNotFoundError(f"Missing repo root: {REPO_ROOT}")

    notebook_source = convert_notebook_to_script_source(NOTEBOOK_PATH)
    patched_source = patch_young_stars_script_source(
        notebook_source,
        output_html=output_html,
        compact_payload=compact_payload,
    )

    run_script_source(patched_source)
    if website_output_html is not None:
        website_output_html = Path(website_output_html)
        if output_html.resolve() != website_output_html.resolve():
            website_output_html.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(output_html, website_output_html)
    return output_html


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-html",
        type=Path,
        default=DEFAULT_OUTPUT_HTML,
        help="HTML output path for the rendered figure.",
    )
    parser.add_argument(
        "--full-payload",
        action="store_true",
        help="Keep repeated per-frame hover, selection, and motion metadata in the HTML.",
    )
    parser.add_argument(
        "--website-copy",
        action="store_true",
        help=f"Also copy the generated HTML to {WEBSITE_OUTPUT_HTML}.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    output_html = run_young_stars_figure(
        output_html=args.output_html,
        compact_payload=not bool(args.full_payload),
        website_output_html=WEBSITE_OUTPUT_HTML if args.website_copy else None,
    )
    print(f"Wrote {output_html}")
    if args.website_copy:
        print(f"Copied {output_html} to {WEBSITE_OUTPUT_HTML}")


if __name__ == "__main__":
    main()
