#!/usr/bin/env python3
"""Render a Taurus/Sigma-focused oviz figure."""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path


HOME_DIR = Path.home()
REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT_CSV = HOME_DIR / "Downloads" / "taurus_core_sigma_age-Feb-2025-v2.csv"
DEFAULT_OUTPUT_HTML = Path(__file__).resolve().with_suffix(".html")
HUNT_CLUSTER_CATALOG_PATH = (
    HOME_DIR
    / "Desktop"
    / "astro_research"
    / "supernovae_map_work"
    / "clusters"
    / "vels_output"
    / "2026-02-04"
    / "cluster_velocities_jan2026.csv"
)
CHRONOS_AGES_PATH = HOME_DIR / "Downloads" / "hunt_sample_chronos_ages_multiprocessing_feb_2026.csv"
CLUSTER_SAMPLE_DATA_PATH = HOME_DIR / "Downloads" / "cluster_sample_data.csv"
EDENHOFER_FITS_PATH = HOME_DIR / "Downloads" / "mean_and_std_xyz-2.fits"
GALACTIC_PLANE_IMAGE_PATH = HOME_DIR / "Downloads" / "Top-down_view_of_the_Milky_Way.jpg"
TAURUS_TRACE_NAME = "Taurus/Sigma Bulk Clusters"
TAURUS_STARS_TRACE_NAME = "Taurus/Sigma Stars (t=0)"
TAURUS_LABEL_TRACE_NAME = "Taurus/Sigma Labels (t=0)"
TAURUS_STAR_COLOR = "#56ff96"

FAMILY_DEFS = (
    ("alphaPer", "Alpha Persei Family", "violet"),
    ("cr135", "Cr 135 Family", "orange"),
    ("m6", "M6 Family", "cyan"),
)
BY_EYE_FAMILY_DEFS = (
    (
        "Lacerta Family",
        [
            "CWNU_1243",
            "HSC_661",
            "UPK_109",
            "HSC_705",
            "CWNU_96",
            "Theia_420",
            "UPK_166",
            "UPK_168",
            "Teutsch_39",
            "OCSN_32",
            "Theia_100",
            "CWNU_311",
        ],
        "#F0E68C",
    ),
    (
        "Trumpler 3 Family",
        [
            "Trumpler_3",
            "Theia_850",
            "CWNU_525",
            "FSR_0686",
            "Theia_1722",
            "UPK_325",
            "UPK_307",
            "FSR_0732",
            "NGC_1960",
            "HSC_1350",
            "HSC_1341",
            "HXHWL_18",
            "HXWHB_8",
            "NGC_1502",
            "CWNU_409",
            "CWNU_364",
            "NGC_1444",
            "UBC_51",
            "Berkeley_14A",
            "CWNU_518",
            "CWNU_205",
            "HSC_1308",
        ],
        "#00CED1",
    ),
    (
        "Proto Orion Family",
        [
            "HSC_1340",
            "NGC_2232",
            "CWNU_1057",
            "CWNU_1111",
            "CWNU_1052",
            "Theia_71",
            "ZHBJZ_1",
        ],
        "#DCDCDC",
    ),
)


def configure_runtime_environment() -> None:
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl")
    os.environ.setdefault("XDG_CACHE_HOME", "/tmp")
    os.environ.setdefault("MPLBACKEND", "Agg")


def require_existing_path(path: Path, label: str) -> Path:
    path = Path(path).expanduser()
    if not path.exists():
        raise FileNotFoundError(f"Missing {label}: {path}")
    return path


def coerce_numeric_columns(df, columns):
    import pandas as pd

    for column in columns:
        if column in df.columns:
            df[column] = pd.to_numeric(df[column], errors="coerce")
    return df


def load_cluster_context():
    import numpy as np
    import pandas as pd

    require_existing_path(HUNT_CLUSTER_CATALOG_PATH, "Hunt cluster catalog")
    require_existing_path(CHRONOS_AGES_PATH, "Chronos ages catalog")
    require_existing_path(CLUSTER_SAMPLE_DATA_PATH, "cluster family catalog")

    df_hunt_full = pd.read_csv(HUNT_CLUSTER_CATALOG_PATH)
    ages_chronos = pd.read_csv(CHRONOS_AGES_PATH, usecols=["name", "age_chronos_mode"])
    df_hunt_full = pd.merge(df_hunt_full, ages_chronos, on="name", how="left")
    df_hunt_full["age_myr"] = df_hunt_full["age_chronos_mode"]

    df_hunt_full = df_hunt_full.drop(
        columns=["x", "y", "z", "U", "V", "W", "U_err", "V_err", "W_err"],
        errors="ignore",
    ).rename(
        columns={
            "U_2026": "U",
            "V_2026": "V",
            "W_2026": "W",
            "U_err_2026": "U_err",
            "V_err_2026": "V_err",
            "W_err_2026": "W_err",
            "x_new": "x",
            "y_new": "y",
            "z_new": "z",
        }
    )

    df_hunt_full = coerce_numeric_columns(
        df_hunt_full,
        [
            "x",
            "y",
            "z",
            "U",
            "V",
            "W",
            "U_err",
            "V_err",
            "W_err",
            "age_myr",
            "n_stars",
            "n_rvs_2026",
            "class_50",
            "cst",
        ],
    )

    df_hunt_good = df_hunt_full.loc[
        (df_hunt_full["U_err"] < 10)
        & (df_hunt_full["V_err"] < 10)
        & (df_hunt_full["W_err"] < 10)
        & (df_hunt_full["U"].notnull())
        & (df_hunt_full["V"].notnull())
        & (df_hunt_full["W"].notnull())
        & (df_hunt_full["age_myr"].notnull())
        & (df_hunt_full["x"].between(-2000, 2000))
        & (df_hunt_full["y"].between(-2000, 2000))
        & (df_hunt_full["z"].between(-300, 300))
        & (df_hunt_full["n_rvs_2026"] >= 3)
        & (df_hunt_full["class_50"] > 0.5)
        & (df_hunt_full["cst"] > 5)
    ].copy()

    df_hunt_60 = df_hunt_good.loc[df_hunt_good["age_myr"] < 60].copy()
    df_hunt_young = df_hunt_good.loc[df_hunt_good["age_myr"] < 15].copy()

    c_s24 = pd.read_csv(CLUSTER_SAMPLE_DATA_PATH, usecols=["name", "family"])
    family_frames = []
    for family_key, label, color in FAMILY_DEFS:
        frame = df_hunt_60.loc[df_hunt_60["name"].isin(c_s24.loc[c_s24["family"] == family_key, "name"])].copy()
        family_frames.append((label, frame, color))

    for label, names, color in BY_EYE_FAMILY_DEFS:
        frame = df_hunt_good.loc[df_hunt_good["name"].isin(names)].copy()
        family_frames.append((label, frame, color))

    finite_columns = ["x", "y", "z", "U", "V", "W", "age_myr"]
    df_hunt_60 = df_hunt_60.replace([np.inf, -np.inf], np.nan).dropna(subset=finite_columns)
    df_hunt_young = df_hunt_young.replace([np.inf, -np.inf], np.nan).dropna(subset=finite_columns)
    family_frames = [
        (label, frame.replace([np.inf, -np.inf], np.nan).dropna(subset=finite_columns), color)
        for label, frame, color in family_frames
    ]
    return df_hunt_60, df_hunt_young, family_frames


def load_taurus_bulk_clusters(csv_path: Path):
    import numpy as np
    import pandas as pd

    csv_path = require_existing_path(csv_path, "Taurus/Sigma CSV")
    taurus = pd.read_csv(csv_path).rename(columns={"X": "x", "Y": "y", "Z": "z"})
    taurus = coerce_numeric_columns(taurus, ["x", "y", "z", "U", "V", "W", "age_myr"])
    taurus = taurus.replace([np.inf, -np.inf], np.nan).dropna(subset=["x", "y", "z", "age_myr"]).copy()
    if "name" not in taurus.columns:
        taurus["name"] = "Taurus/Sigma"
    taurus["name"] = taurus["name"].fillna("Taurus/Sigma").astype(str)
    taurus["has_velocity"] = taurus[["U", "V", "W"]].notna().all(axis=1)

    source_column = "source_id" if "source_id" in taurus.columns else "name"
    bulk = (
        taurus.groupby("name", as_index=False)
        .agg(
            x=("x", "median"),
            y=("y", "median"),
            z=("z", "median"),
            U=("U", "median"),
            V=("V", "median"),
            W=("W", "median"),
            age_myr=("age_myr", "median"),
            n_stars=(source_column, "size"),
            n_velocity_stars=("has_velocity", "sum"),
        )
        .replace([np.inf, -np.inf], np.nan)
    )
    required_columns = ["x", "y", "z", "U", "V", "W", "age_myr"]
    missing_velocity = bulk.loc[bulk[["U", "V", "W"]].isna().any(axis=1), "name"].tolist()
    if missing_velocity:
        print(
            "Skipping Taurus/Sigma bulk clusters without complete median U,V,W: "
            + ", ".join(missing_velocity)
        )
    bulk = bulk.dropna(subset=required_columns).copy()
    if bulk.empty:
        raise ValueError("No Taurus/Sigma bulk clusters have complete phase-space coordinates.")
    return bulk


def load_taurus_star_sample(csv_path: Path):
    import numpy as np
    import pandas as pd

    csv_path = require_existing_path(csv_path, "Taurus/Sigma CSV")
    taurus = pd.read_csv(csv_path).rename(columns={"X": "x", "Y": "y", "Z": "z"})
    taurus = coerce_numeric_columns(taurus, ["x", "y", "z", "age_myr"])
    taurus = taurus.replace([np.inf, -np.inf], np.nan).dropna(subset=["x", "y", "z"]).copy()
    if "name" not in taurus.columns:
        taurus["name"] = "Taurus/Sigma"
    taurus["name"] = taurus["name"].fillna("Taurus/Sigma").astype(str)
    return taurus


def build_taurus_star_customdata(taurus_stars):
    import numpy as np

    x = taurus_stars["x"].to_numpy(dtype=float)
    y = taurus_stars["y"].to_numpy(dtype=float)
    z = taurus_stars["z"].to_numpy(dtype=float)
    dist = np.sqrt(x**2 + y**2 + z**2)
    with np.errstate(invalid="ignore", divide="ignore"):
        l_deg = np.mod(np.rad2deg(np.arctan2(y, x)), 360.0)
        b_deg = np.rad2deg(np.arcsin(np.clip(z / np.where(dist > 0, dist, np.nan), -1.0, 1.0)))

    if "source_id" in taurus_stars.columns:
        star_ids = taurus_stars["source_id"].astype(str).to_numpy(dtype=object)
    else:
        star_ids = np.arange(len(taurus_stars)).astype(str)
    cluster_names = np.array([f"Taurus/Sigma star {star_id}" for star_id in star_ids], dtype=object)
    trace_colors = np.repeat(TAURUS_STAR_COLOR, len(taurus_stars)).astype(object)
    n_stars = np.ones(len(taurus_stars), dtype=float)
    no_motion_age = np.full(len(taurus_stars), np.nan, dtype=float)

    return np.column_stack(
        (
            no_motion_age,
            no_motion_age,
            l_deg,
            b_deg,
            dist,
            x,
            y,
            z,
            cluster_names,
            trace_colors,
            n_stars,
        )
    )


def build_taurus_static_traces(taurus_bulk, taurus_stars):
    import plotly.graph_objects as go

    star_trace = go.Scatter3d(
        x=taurus_stars["x"].to_numpy(dtype=float),
        y=taurus_stars["y"].to_numpy(dtype=float),
        z=taurus_stars["z"].to_numpy(dtype=float),
        mode="markers",
        marker=dict(
            size=0.35,
            color=TAURUS_STAR_COLOR,
            opacity=0.75,
            symbol="circle",
            line=dict(width=0),
        ),
        customdata=build_taurus_star_customdata(taurus_stars),
        name=TAURUS_STARS_TRACE_NAME,
        hoverinfo="skip",
        meta={"trace_kind": "stars", "static": True},
    )

    label_trace = go.Scatter3d(
        x=taurus_bulk["x"].to_numpy(dtype=float),
        y=taurus_bulk["y"].to_numpy(dtype=float),
        z=(taurus_bulk["z"] + 7.0).to_numpy(dtype=float),
        mode="text",
        text=taurus_bulk["name"].astype(str).to_numpy(dtype=object),
        textposition="top center",
        textfont=dict(color="#dcffe8", size=18, family="Helvetica"),
        name=TAURUS_LABEL_TRACE_NAME,
        hoverinfo="skip",
        meta={
            "screen_stable_text": True,
            "screen_px": 18.0,
            "trace_kind": "labels",
        },
    )
    return [star_trace, label_trace]


def default_trace_off(scene_spec: dict, trace_name: str) -> None:
    legend_items = (scene_spec.get("legend") or {}).get("items") or []
    trace_key = next(
        (
            item.get("key")
            for item in legend_items
            if item.get("kind") == "trace" and item.get("name") == trace_name
        ),
        None,
    )
    if not trace_key:
        return

    scene_spec.setdefault("initial_state", {}).setdefault("legend_state", {})[str(trace_key)] = False
    for group_defaults in (scene_spec.get("group_visibility") or {}).values():
        if isinstance(group_defaults, dict) and str(trace_key) in group_defaults:
            group_defaults[str(trace_key)] = "legendonly"


def main_style_initial_state(*, compact_payload: bool, taurus_center: dict[str, float]):
    camera_target = {
        "x": 0.0,
        "y": 0.0,
        "z": 0.0,
    }
    camera_position = {
        "x": camera_target["x"] + 280.0,
        "y": camera_target["y"] - 470.0,
        "z": camera_target["z"] + 280.0,
    }
    return {
        "current_group": "Taurus Sigma",
        "click_selection_enabled": False,
        "compact_payload_enabled": bool(compact_payload),
        "active_volume_key": "volume-0",
        "legend_state": {"volume-0": True},
        "volume_state_by_key": {"volume-0": {"visible": True}},
        "galaxy_image": True,
        "galaxy_image_path": str(GALACTIC_PLANE_IMAGE_PATH),
        "galaxy_image_size_pc": 40000.0,
        "galaxy_image_opacity": 0.35,
        "galaxy_image_hide_below_scale_bar_pc": 420.0,
        "galaxy_image_fade_start_scale_bar_pc": 700.0,
        "sky_dome_enabled": True,
        "sky_dome_background_mode": "live_aladin",
        "sky_dome_source": "aladin",
        "sky_dome_projection": "TAN",
        "sky_dome_capture_width_px": 4096,
        "sky_dome_capture_height_px": 2048,
        "sky_dome_capture_format": "image/jpeg",
        "sky_dome_capture_quality": 0.94,
        "sky_dome_radius_pc": 40000.0,
        "sky_dome_opacity": 1.0,
        "sky_dome_force_visible": False,
        "sky_dome_full_opacity_scale_bar_pc": 120.0,
        "sky_dome_fade_out_scale_bar_pc": 360.0,
        "sky_layers": [
            {
                "key": "P/Mellinger/color",
                "label": "Mellinger Color",
                "survey": "P/Mellinger/color",
                "opacity": 1.0,
                "visible": False,
            },
            {
                "key": "P/PLANCK/R2/HFI/color",
                "label": "Planck Dust Emission Color",
                "survey": "P/PLANCK/R2/HFI/color",
                "opacity": 1.0,
                "visible": True,
            },
        ],
        "active_sky_layer_key": "P/PLANCK/R2/HFI/color",
        "initial_zoom_anchor": camera_target,
        "camera": {
            "position": camera_position,
            "target": camera_target,
            "up": {"x": 0.0, "y": 0.0, "z": 1.0},
        },
        "global_controls": {
            "camera_fov": 55.0,
            "point_size_scale": 0.5,
            "point_glow_strength": 0.75,
            "fade_in_time_myr": 8.0,
        },
    }


def render_figure(input_csv: Path, output_html: Path, *, compact_payload: bool = True) -> Path:
    configure_runtime_environment()

    if str(REPO_ROOT) not in sys.path:
        sys.path.insert(0, str(REPO_ROOT))

    import numpy as np
    import pandas as pd

    from oviz import Animate3D, Trace, TraceCollection

    output_html = Path(output_html)
    output_html.parent.mkdir(parents=True, exist_ok=True)

    df_hunt_60, df_hunt_young, family_frames = load_cluster_context()
    taurus = load_taurus_bulk_clusters(input_csv)
    taurus_stars = load_taurus_star_sample(input_csv)
    taurus_center = {
        "x": float(np.nanmedian(taurus["x"].to_numpy(dtype=float))),
        "y": float(np.nanmedian(taurus["y"].to_numpy(dtype=float))),
        "z": float(np.nanmedian(taurus["z"].to_numpy(dtype=float))),
    }

    sun = pd.DataFrame(
        {
            "name": ["Sun"],
            "family": ["Sun"],
            "age_myr": [8000.0],
            "U": [0.0],
            "V": [0.0],
            "W": [0.0],
            "x": [0.0],
            "y": [0.0],
            "z": [27.0],
            "n_stars": [1],
        }
    )

    taurus_trace = Trace(
        taurus[["x", "y", "z", "U", "V", "W", "name", "age_myr", "n_stars"]],
        data_name=TAURUS_TRACE_NAME,
        min_size=8.0,
        max_size=8.0,
        color="#56ff96",
        opacity=0.25,
        marker_style="circle",
        show_tracks=False,
    )
    traces = TraceCollection(
        [
            Trace(sun, data_name="Sun", min_size=5, max_size=5, color="yellow"),
            taurus_trace,
            Trace(
                df_hunt_60,
                data_name="Clusters (< 60 Myr)",
                min_size=0,
                max_size=7,
                color="grey",
                opacity=0.34,
            ),
            Trace(
                df_hunt_young,
                data_name="Young Clusters (< 15 Myr)",
                min_size=1.25,
                max_size=7,
                color="red",
                opacity=0.86,
            ),
            *[
                Trace(
                    family_df,
                    data_name=label,
                    min_size=0,
                    max_size=7,
                    color=color,
                    opacity=0.95,
                )
                for label, family_df, color in family_frames
                if not family_df.empty
            ],
        ]
    )

    family_names = [
        label
        for label, family_df, _color in family_frames
        if not family_df.empty
    ]
    trace_groupings = {
        "Taurus Sigma": [
            "Sun",
            TAURUS_TRACE_NAME,
        ],
        "Taurus + Context": [
            "Sun",
            TAURUS_TRACE_NAME,
            "Clusters (< 60 Myr)",
            "Young Clusters (< 15 Myr)",
            *family_names,
        ],
        "Cluster Families": ["Sun", TAURUS_TRACE_NAME, *family_names],
        "Young Clusters": ["Sun", TAURUS_TRACE_NAME, "Young Clusters (< 15 Myr)"],
        "<60 Myr Cluster Sample": ["Sun", TAURUS_TRACE_NAME, "Clusters (< 60 Myr)"],
    }

    plot_3d = Animate3D(
        data_collection=traces,
        xyz_widths=(2000, 2000, 400),
        figure_theme="dark",
        trace_grouping_dict=trace_groupings,
    )

    edenhofer_volume = {
        "name": "Edenhofer+2024 Dust",
        "path": str(require_existing_path(EDENHOFER_FITS_PATH, "Edenhofer FITS map")),
        "hdu": "MEAN",
        "max_resolution": 640,
        "max_resolution_cap": 640,
        "opacity": 1.0,
        "samples": 200,
        "alpha_coef": 200,
        "vmin": 0,
        "vmax": 0.0385,
        "colormap": "Greys",
        "supports_show_all_times": True,
        "co_rotate_with_frame": True,
    }

    figure = plot_3d.make_plot(
        time=np.round(np.arange(0, -66, -1), 1),
        show=False,
        save_name=None,
        static_traces=build_taurus_static_traces(taurus, taurus_stars),
        static_traces_times=[[0.0], [0.0]],
        static_traces_legendonly=False,
        focus_group=None,
        fade_in_time=8,
        fade_in_and_out=False,
        show_age_kde_inset=True,
        age_kde_bandwidth_myr=2,
        include_spiral_arms=False,
        galactic_mode=True,
        show_galactic_guides=False,
        camera_zoom_factor=5,
        show_gc_line=True,
        show_milky_way_model=False,
        renderer="threejs",
        enable_sky_panel=True,
        volumes=[edenhofer_volume],
        threejs_initial_state=main_style_initial_state(
            compact_payload=compact_payload,
            taurus_center=taurus_center,
        ),
    )
    default_trace_off(figure.scene_spec, TAURUS_LABEL_TRACE_NAME)
    figure.write_html(str(output_html))
    return output_html


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input-csv",
        type=Path,
        default=DEFAULT_INPUT_CSV,
        help="Taurus/Sigma source CSV.",
    )
    parser.add_argument(
        "--output-html",
        type=Path,
        default=DEFAULT_OUTPUT_HTML,
        help="HTML output path.",
    )
    parser.add_argument(
        "--full-payload",
        action="store_true",
        help="Disable compact scene payload storage.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    output_html = render_figure(
        input_csv=args.input_csv,
        output_html=args.output_html,
        compact_payload=not bool(args.full_payload),
    )
    print(f"Wrote {output_html}")


if __name__ == "__main__":
    main()
