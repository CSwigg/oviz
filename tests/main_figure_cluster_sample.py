#!/usr/bin/env python3
"""Render a main-figure variant with the Figure 3 cluster sample and dust only."""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl")
os.environ.setdefault("XDG_CACHE_HOME", "/tmp")
os.environ.setdefault("MPLBACKEND", "Agg")

import numpy as np
import pandas as pd


HOME_DIR = Path.home()
REPO_ROOT = Path(__file__).resolve().parents[1]
FULL_CLUSTER_SAMPLE_PATH = (
    HOME_DIR
    / "Desktop"
    / "astro_research"
    / "supernovae_map"
    / "paper"
    / "analysis"
    / "Figure_3_full_cluster_sample_clusters.csv"
)
SWIGGUM_FAMILY_SAMPLE_PATH = (
    HOME_DIR
    / "Desktop"
    / "astro_research"
    / "supernovae_map"
    / "paper"
    / "analysis"
    / "Figure_3_full_cluster_sample_swiggum_family_clusters.csv"
)
VARIANT_INPUT_PATH = (
    HOME_DIR
    / "Desktop"
    / "astro_research"
    / "supernovae_map"
    / "visualizer"
    / "plots"
    / "variants"
    / "galpy_almeida2024_radius_total_uniform_sphere"
    / "cluster_input_variant.csv.gz"
)
DUST_FITS_PATH = HOME_DIR / "Downloads" / "mean_and_std_xyz-2.fits"
GALACTIC_PLANE_IMAGE_PATH = HOME_DIR / "Downloads" / "Top-down_view_of_the_Milky_Way.jpg"
DEFAULT_OUTPUT_HTML = Path(__file__).resolve().with_suffix(".html")
DEFAULT_LOOKBACK_MYR = 105
MAIN_FIGURE_VOLUME_Z_CLIP_BOUNDS = (-400.0, 400.0)
MAIN_FIGURE_DUST_MAX_RESOLUTION = 512
MAIN_FIGURE_DUST_MAX_RESOLUTION_CAP = 512
MAIN_FIGURE_DUST_SAMPLES = 200
REQUIRED_TRACE_COLUMNS = ("x", "y", "z", "U", "V", "W", "age_myr")


def _prepare_import_path() -> None:
    repo_root = str(REPO_ROOT)
    if repo_root not in sys.path:
        sys.path.insert(0, repo_root)


def _require_file(path: Path, label: str) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Missing {label}: {path}")


def _load_cluster_sample(path: Path, label: str) -> pd.DataFrame:
    _require_file(path, label)
    df = pd.read_csv(path)
    missing_columns = [col for col in ("name", *REQUIRED_TRACE_COLUMNS) if col not in df.columns]
    if missing_columns:
        raise ValueError(f"{label} is missing required columns: {missing_columns}")

    df = df.copy()
    df["name"] = df["name"].astype(str)
    for column in REQUIRED_TRACE_COLUMNS:
        df[column] = pd.to_numeric(df[column], errors="coerce")
    if "n_stars" in df.columns:
        df["n_stars"] = pd.to_numeric(df["n_stars"], errors="coerce").fillna(1)
    else:
        df["n_stars"] = 1

    valid_mask = df[["name", *REQUIRED_TRACE_COLUMNS]].notna().all(axis=1)
    dropped = int((~valid_mask).sum())
    if dropped:
        print(f"Dropping {dropped} {label} row(s) with missing orbit inputs.")
    df = df.loc[valid_mask].reset_index(drop=True)
    if df.empty:
        raise ValueError(f"{label} has no rows with complete orbit inputs.")
    return df


def _build_sun_trace_df() -> pd.DataFrame:
    return pd.DataFrame(
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


def _build_dust_volume() -> dict:
    import matplotlib.pyplot as plt

    _require_file(DUST_FITS_PATH, "Edenhofer dust FITS cube")
    greys_cmap = plt.get_cmap("Greys")
    return {
        "name": "Edenhofer+2024 Dust",
        "path": str(DUST_FITS_PATH),
        "hdu": "MEAN",
        "clip_bounds": {"z": list(MAIN_FIGURE_VOLUME_Z_CLIP_BOUNDS)},
        "max_resolution": MAIN_FIGURE_DUST_MAX_RESOLUTION,
        "max_resolution_cap": MAIN_FIGURE_DUST_MAX_RESOLUTION_CAP,
        "opacity": 1.0,
        "samples": MAIN_FIGURE_DUST_SAMPLES,
        "alpha_coef": 200,
        "vmin": 0,
        "vmax": 0.0385,
        "colormap": greys_cmap,
        "supports_show_all_times": True,
        "co_rotate_with_frame": True,
    }


def run_cluster_sample_figure(
    output_html: Path = DEFAULT_OUTPUT_HTML,
    *,
    full_cluster_sample_path: Path = FULL_CLUSTER_SAMPLE_PATH,
    swiggum_family_sample_path: Path = SWIGGUM_FAMILY_SAMPLE_PATH,
    lookback_myr: int = DEFAULT_LOOKBACK_MYR,
    theme_key: str | None = None,
    compact_payload: bool = True,
) -> Path:
    _prepare_import_path()
    from oviz import Animate3D, Trace, TraceCollection

    output_html = Path(output_html)
    output_html.parent.mkdir(parents=True, exist_ok=True)

    full_sample = _load_cluster_sample(full_cluster_sample_path, "full cluster sample")
    swiggum_family = _load_cluster_sample(swiggum_family_sample_path, "Swiggum-family sample")
    _require_file(VARIANT_INPUT_PATH, "broader variant input table")

    print(f"Full cluster sample: {len(full_sample)} plotted rows from {full_cluster_sample_path}")
    print(f"Swiggum-family overlay: {len(swiggum_family)} plotted rows from {swiggum_family_sample_path}")
    print(f"Broader variant source: {VARIANT_INPUT_PATH}")

    sun_trace = Trace(
        _build_sun_trace_df(),
        data_name="Sun",
        min_size=4.0,
        max_size=8.0,
        color="yellow",
        opacity=1.0,
        marker_style="circle",
        show_tracks=False,
        size_by_n_stars=False,
    )
    full_sample_trace = Trace(
        full_sample,
        data_name="Full Cluster Sample",
        min_size=0.0,
        max_size=5.5,
        color="#9ba1aa",
        opacity=0.48,
        marker_style="circle",
        show_tracks=False,
        size_by_n_stars=False,
    )
    swiggum_family_trace = Trace(
        swiggum_family,
        data_name="Swiggum Family",
        min_size=1.5,
        max_size=8.0,
        color="#2f80ed",
        opacity=1.0,
        marker_style="circle",
        show_tracks=False,
        size_by_n_stars=False,
    )

    traces = TraceCollection(
        [
            sun_trace,
            full_sample_trace,
            swiggum_family_trace,
        ]
    )
    trace_groupings = {
        "Clusters": ["Sun", "Full Cluster Sample", "Swiggum Family"],
    }

    time_int = np.round(np.arange(0, -(int(lookback_myr) + 1), -1), 1)
    plot_3d = Animate3D(
        data_collection=traces,
        xyz_widths=(2000, 2000, 400),
        figure_theme="dark",
        trace_grouping_dict=trace_groupings,
    )

    initial_state = {
        "current_group": "Clusters",
        "click_selection_enabled": False,
        "compact_payload_enabled": bool(compact_payload),
        "scene_float_precision": 2,
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
                "visible": True,
            },
            {
                "key": "P/PLANCK/R2/HFI/color",
                "label": "Planck Dust Emission Color",
                "survey": "P/PLANCK/R2/HFI/color",
                "opacity": 1.0,
                "visible": True,
            },
        ],
        "active_sky_layer_key": "P/Mellinger/color",
    }
    if theme_key:
        initial_state["global_controls"] = {"theme_key": theme_key}

    plot_3d.make_plot(
        time=time_int,
        show=False,
        save_name=str(output_html),
        static_traces=None,
        static_traces_times=None,
        static_traces_legendonly=True,
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
        volumes=[_build_dust_volume()],
        threejs_initial_state=initial_state,
    )
    return output_html


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-html",
        type=Path,
        default=DEFAULT_OUTPUT_HTML,
        help="HTML output path for the rendered cluster-sample figure.",
    )
    parser.add_argument(
        "--full-cluster-sample",
        type=Path,
        default=FULL_CLUSTER_SAMPLE_PATH,
        help="CSV for the full cluster sample plotted in grey.",
    )
    parser.add_argument(
        "--swiggum-family-sample",
        type=Path,
        default=SWIGGUM_FAMILY_SAMPLE_PATH,
        help="CSV for the Swiggum-family overlay plotted in blue.",
    )
    parser.add_argument(
        "--lookback-myr",
        type=int,
        default=DEFAULT_LOOKBACK_MYR,
        help="Lookback interval for orbit integration.",
    )
    parser.add_argument(
        "--theme-key",
        default=None,
        help="Optional Three.js color theme preset to force into the figure initial state.",
    )
    parser.add_argument(
        "--full-payload",
        action="store_true",
        help="Keep repeated per-frame hover, selection, and motion metadata in the HTML.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    output_html = run_cluster_sample_figure(
        output_html=args.output_html,
        full_cluster_sample_path=args.full_cluster_sample,
        swiggum_family_sample_path=args.swiggum_family_sample,
        lookback_myr=args.lookback_myr,
        theme_key=args.theme_key,
        compact_payload=not bool(args.full_payload),
    )
    print(f"Wrote {output_html}")


if __name__ == "__main__":
    main()
