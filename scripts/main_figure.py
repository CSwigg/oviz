#!/usr/bin/env python3
"""Run the external main_figure notebook as a repeatable smoke test."""

from __future__ import annotations

import argparse
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path


REPO_ROOT = Path("/Users/cam/Desktop/oviz")
NOTEBOOK_PATH = Path("/Users/cam/Desktop/astro_research/radcliffe/oviz_notebooks/main_figure.ipynb")
SUPERNOVAE_ROOT = Path("/Users/cam/Desktop/astro_research/supernovae_map")
SUPERNOVAE_CATALOG_PATH = SUPERNOVAE_ROOT / "paper" / "solar_encounter_catalog_current.csv.gz"
DEFAULT_OUTPUT_HTML = Path("/tmp/main_figure_uncodixified.html")
DEFAULT_UI_DESIGN_KEY = "uncodixified"


def convert_notebook_to_script(notebook_path: Path, workdir: Path) -> Path:
    subprocess.run(
        [
            "jupyter",
            "nbconvert",
            "--to",
            "script",
            str(notebook_path),
            "--output",
            "main_figure_nb.py",
            "--output-dir",
            str(workdir),
        ],
        check=True,
        cwd=str(workdir),
    )
    candidates = sorted(workdir.glob("main_figure_nb.py*"))
    if not candidates:
        raise FileNotFoundError("nbconvert did not produce a Python script.")
    return candidates[0]


def build_supernova_volume_source_block() -> str:
    return f"""
import os
os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl")
os.environ.setdefault("XDG_CACHE_HOME", "/tmp")
os.environ.setdefault("MPLBACKEND", "Agg")

from pathlib import Path
import sys

SUPERNOVAE_MAP_ROOT = Path({str(SUPERNOVAE_ROOT)!r})
SUPERNOVAE_CATALOG_PATH = Path({str(SUPERNOVAE_CATALOG_PATH)!r})
if str(SUPERNOVAE_MAP_ROOT) not in sys.path:
    sys.path.insert(0, str(SUPERNOVAE_MAP_ROOT))

from mapper import OrbitConfig
from mapper.orbits import ClusterOrbiter


def _sn_build_edges(half_width_pc, voxel_size_pc):
    n_bins = int(np.ceil((2.0 * float(half_width_pc)) / float(voxel_size_pc)))
    n_bins = max(n_bins, 1)
    half_extent = 0.5 * n_bins * float(voxel_size_pc)
    return np.linspace(-half_extent, half_extent, n_bins + 1, dtype=float)


def _sn_nearest_values(sorted_times, query_times):
    idx = np.searchsorted(sorted_times, query_times, side="left")
    idx = np.clip(idx, 1, len(sorted_times) - 1)
    left = sorted_times[idx - 1]
    right = sorted_times[idx]
    choose_left = (query_times - left) <= (right - query_times)
    nearest_idx = np.where(choose_left, idx - 1, idx)
    return sorted_times[nearest_idx]


def _sn_rotating_local_to_fixed_galactocentric(x_rot_pc, y_rot_pc, z_rot_pc, time_myr, *, orbit_cfg):
    r_sun_pc = float(orbit_cfg.ro_kpc) * 1000.0
    omega_rot_per_myr = (float(orbit_cfg.vo_kms) / float(orbit_cfg.ro_kpc)) / 10.0
    time_code_units = np.asarray(time_myr, dtype=float) * 0.01022

    local_x = r_sun_pc - np.asarray(x_rot_pc, dtype=float)
    local_y = np.asarray(y_rot_pc, dtype=float)
    r_pc = np.hypot(local_x, local_y)
    phi = np.arctan2(local_y, local_x)
    theta = phi + omega_rot_per_myr * time_code_units - (0.5 * np.pi)

    x_gc_pc = r_pc * np.sin(theta)
    y_gc_pc = r_pc * np.cos(theta)
    z_gc_pc = np.asarray(z_rot_pc, dtype=float)
    return x_gc_pc, y_gc_pc, z_gc_pc


def _sn_build_reference_orbit_df(time_grid, *, orbit_cfg):
    rf = pd.DataFrame(
        {{
            "name": ["rf"],
            "x": [0.0],
            "y": [0.0],
            "z": [0.0],
            "U": [-11.1],
            "V": [-12.24],
            "W": [-7.25],
            "x_err": [0.0],
            "y_err": [0.0],
            "z_err": [0.0],
            "U_err": [0.0],
            "V_err": [0.0],
            "W_err": [0.0],
        }}
    )
    rf_df = ClusterOrbiter(rf, config=orbit_cfg).run_mc_integration(time_grid, n_samples=1).copy()
    if rf_df.empty:
        raise ValueError("Reference orbit integration returned no rows.")
    rf_df = rf_df.loc[rf_df["sample_id"] == 0, ["time_myr", "x_gc_pc", "y_gc_pc", "z_gc_pc"]].copy()
    rf_df = rf_df.rename(
        columns={{
            "x_gc_pc": "rf_x_gc_pc",
            "y_gc_pc": "rf_y_gc_pc",
            "z_gc_pc": "rf_z_gc_pc",
        }}
    )
    if not np.isclose(rf_df["time_myr"].to_numpy(dtype=float), 0.0, atol=1e-9).any():
        x0_gc_pc, y0_gc_pc, z0_gc_pc = _sn_rotating_local_to_fixed_galactocentric(
            np.array([0.0]),
            np.array([0.0]),
            np.array([0.0]),
            np.array([0.0]),
            orbit_cfg=orbit_cfg,
        )
        rf_df = pd.concat(
            [
                rf_df,
                pd.DataFrame(
                    {{
                        "time_myr": [0.0],
                        "rf_x_gc_pc": [float(x0_gc_pc[0])],
                        "rf_y_gc_pc": [float(y0_gc_pc[0])],
                        "rf_z_gc_pc": [float(z0_gc_pc[0])],
                    }}
                ),
            ],
            ignore_index=True,
        )
    return rf_df.sort_values("time_myr").reset_index(drop=True)


def _sn_assign_supernova_positions_in_trace_frame(sne_df, *, rf_df, orbit_cfg):
    df = sne_df.copy()
    time_grid = np.sort(rf_df["time_myr"].unique())
    if len(time_grid) == 0:
        raise ValueError("Reference orbit lookup is empty.")

    df["trace_time_myr"] = _sn_nearest_values(time_grid, df["time_of_death_myr"].to_numpy(dtype=float))
    joined = df[["trace_time_myr"]].merge(
        rf_df.rename(columns={{"time_myr": "trace_time_myr"}}),
        on="trace_time_myr",
        how="left",
        sort=False,
    )
    x_gc_pc, y_gc_pc, z_gc_pc = _sn_rotating_local_to_fixed_galactocentric(
        df["x_pc"].to_numpy(dtype=float),
        df["y_pc"].to_numpy(dtype=float),
        df["z_pc"].to_numpy(dtype=float),
        df["trace_time_myr"].to_numpy(dtype=float),
        orbit_cfg=orbit_cfg,
    )
    df["x_trace_frame_pc"] = x_gc_pc - joined["rf_x_gc_pc"].to_numpy(dtype=float)
    df["y_trace_frame_pc"] = y_gc_pc - joined["rf_y_gc_pc"].to_numpy(dtype=float)
    df["z_trace_frame_pc"] = z_gc_pc - joined["rf_z_gc_pc"].to_numpy(dtype=float)
    return df


def _sn_build_supernova_volume_layers(
    sne_df,
    *,
    time_grid,
    x_edges,
    y_edges,
    z_edges,
    gaussian_sigma_vox=0.8,
    time_window_half_width_myr=5.0,
):
    work = sne_df.copy()
    work = work[np.isfinite(work["time_of_death_myr"])]
    work = work[
        np.isfinite(work["x_trace_frame_pc"])
        & np.isfinite(work["y_trace_frame_pc"])
        & np.isfinite(work["z_trace_frame_pc"])
    ]
    work = work.loc[
        (work["time_of_death_myr"] <= float(np.max(time_grid)) + float(time_window_half_width_myr))
        & (work["time_of_death_myr"] >= float(np.min(time_grid)) - float(time_window_half_width_myr))
    ].copy()

    try:
        from scipy.ndimage import gaussian_filter
    except ImportError:
        gaussian_filter = None

    smoothed_by_time = {{}}
    positive_values = []

    for time_value in time_grid:
        frame_events = work.loc[
            np.abs(work["time_of_death_myr"].to_numpy(dtype=float) - float(time_value))
            <= float(time_window_half_width_myr) + 1e-9
        ]
        if frame_events.empty:
            cube_zyx = np.zeros(
                (len(z_edges) - 1, len(y_edges) - 1, len(x_edges) - 1),
                dtype=np.float32,
            )
        else:
            sample = frame_events[["z_trace_frame_pc", "y_trace_frame_pc", "x_trace_frame_pc"]].to_numpy(dtype=float)
            cube_zyx, _ = np.histogramdd(sample, bins=(z_edges, y_edges, x_edges))
            cube_zyx = cube_zyx.astype(np.float32, copy=False)

        if gaussian_filter is not None and float(gaussian_sigma_vox) > 0:
            cube_zyx = gaussian_filter(cube_zyx, sigma=float(gaussian_sigma_vox), mode="constant").astype(
                np.float32,
                copy=False,
            )

        smoothed_by_time[float(time_value)] = cube_zyx
        positive = cube_zyx[cube_zyx > 0]
        if positive.size:
            positive_values.append(positive)

    if positive_values:
        all_positive = np.concatenate(positive_values)
        data_max = float(np.nanmax(all_positive))
        default_vmin = float(np.nanquantile(all_positive, 0.70))
        default_vmax = float(np.nanquantile(all_positive, 0.995))
    else:
        data_max = 1.0
        default_vmin = 0.05
        default_vmax = 1.0

    if not data_max > 0:
        data_max = 1.0
    if not default_vmax > default_vmin:
        default_vmin = 0.05 * data_max
        default_vmax = data_max

    bounds = {{
        "x": [float(x_edges[0]), float(x_edges[-1])],
        "y": [float(y_edges[0]), float(y_edges[-1])],
        "z": [float(z_edges[0]), float(z_edges[-1])],
    }}

    volumes = []
    for time_index, time_value in enumerate(time_grid):
        volumes.append(
            {{
                "key": f"supernova-density-{{time_index:03d}}",
                "state_key": "supernova-density",
                "state_name": "Supernova Density",
                "name": f"Supernova Density | {{abs(float(time_value)):.0f}} Myr",
                "time_myr": float(time_value),
                "data": smoothed_by_time[float(time_value)],
                "bounds": bounds,
                "data_range": [0.0, float(data_max)],
                "vmin": float(default_vmin),
                "vmax": float(default_vmax),
                "opacity": 0.24,
                "alpha_coef": 95.0,
                "gradient_step": 0.01,
                "samples": 240,
                "colormap": "ice",
                "interpolation": True,
                "visible": False,
            }}
        )
    return volumes


def _build_supernova_volumes_for_main_figure(time_values):
    if not SUPERNOVAE_CATALOG_PATH.exists():
        raise FileNotFoundError(f"Missing supernova catalog: {{SUPERNOVAE_CATALOG_PATH}}")

    time_grid = np.asarray(time_values, dtype=float)
    orbit_cfg = OrbitConfig(backend="galpy", spiral_model="none")
    sne_df = pd.read_csv(
        SUPERNOVAE_CATALOG_PATH,
        usecols=["time_of_death_myr", "x_pc", "y_pc", "z_pc"],
    )
    rf_df = _sn_build_reference_orbit_df(time_grid, orbit_cfg=orbit_cfg)
    sne_trace_df = _sn_assign_supernova_positions_in_trace_frame(
        sne_df,
        rf_df=rf_df,
        orbit_cfg=orbit_cfg,
    )
    x_edges = _sn_build_edges(2000.0, 50.0)
    y_edges = _sn_build_edges(2000.0, 50.0)
    z_edges = _sn_build_edges(400.0, 50.0)
    return _sn_build_supernova_volume_layers(
        sne_trace_df,
        time_grid=time_grid,
        x_edges=x_edges,
        y_edges=y_edges,
        z_edges=z_edges,
        gaussian_sigma_vox=0.8,
        time_window_half_width_myr=5.0,
    )


supernova_volumes = _build_supernova_volumes_for_main_figure(time_int)
""".strip()


def patch_script_source(source: str, output_html: Path, ui_design_key: str | None, theme_key: str | None) -> str:
    output_html_str = str(output_html)
    source, replaced = re.subn(
        r"(?m)^save_name\s*=\s*['\"][^'\"]+['\"]",
        f"save_name = {output_html_str!r}",
        source,
        count=1,
    )
    if replaced != 1:
        raise RuntimeError("Could not rewrite save_name in converted notebook script.")

    control_bits = []
    if theme_key:
        control_bits.append(f"'theme_key': {theme_key!r}")

    make_plot_injections = []
    if ui_design_key:
        make_plot_injections.append(f"threejs_ui_design={ui_design_key!r}")

    if "threejs_initial_state=" not in source:
        initial_state_bits = [
            "'current_group': 'Clusters'",
            "'click_selection_enabled': False",
            "'active_volume_key': 'supernova-density'",
            "'legend_state': {'volume-0': False, 'supernova-density': True}",
            "'volume_state_by_key': {'volume-0': {'visible': False}, 'supernova-density': {'visible': True}}",
        ]
        if control_bits:
            initial_state_bits.append("'global_controls': {" + ", ".join(control_bits) + "}")
        injection_lines = []
        if make_plot_injections:
            injection_lines.extend(f"    {entry},\n" for entry in make_plot_injections)
        injection_lines.append(
            '    threejs_initial_state={'
            + ", ".join(initial_state_bits)
            + '},\n'
        )
        source, injected = re.subn(
            r'(renderer\s*=\s*"threejs",\s*\n)',
            (
                '\\1'
                + "".join(injection_lines)
            ),
            source,
            count=1,
        )
        if injected != 1:
            raise RuntimeError("Could not inject threejs_initial_state into make_plot call.")
    elif make_plot_injections and "threejs_ui_design=" not in source:
        source, injected = re.subn(
            r'(renderer\s*=\s*"threejs",\s*\n)',
            '\\1' + "".join(f"    {entry},\n" for entry in make_plot_injections),
            source,
            count=1,
        )
        if injected != 1:
            raise RuntimeError("Could not inject threejs_ui_design into make_plot call.")

    if "df_hunt_good" in source and "Full Cluster Catalog" not in source:
        catalog_trace_block = """
full_catalog_trace = Trace(
    df_hunt_good,
    data_name='Full Cluster Catalog',
    min_size=0.0,
    max_size=4.0,
    color='#8f959d',
    opacity=0.16,
    marker_style='circle',
    show_tracks=False,
    size_by_n_stars=False,
)
""".strip()

        source, inserted_trace = re.subn(
            r"(?m)^(sun_trace\s*=\s*Trace\(sun,.*\)\s*)$",
            r"\1\n\n" + catalog_trace_block,
            source,
            count=1,
        )
        if inserted_trace != 1:
            raise RuntimeError("Could not inject Full Cluster Catalog trace definition.")

        source, inserted_collection = re.subn(
            r"(traces\s*=\s*TraceCollection\(\[\s*\n\s*sun_trace,\s*\n)",
            r"\1    full_catalog_trace,\n",
            source,
            count=1,
        )
        if inserted_collection != 1:
            raise RuntimeError("Could not add Full Cluster Catalog to the TraceCollection.")

        grouping_patch = """
_all_group_names = ['Sun', 'Full Cluster Catalog']
for _group_names in trace_groupings.values():
    for _name in _group_names:
        if _name not in _all_group_names and _name != 'Full Cluster Catalog':
            _all_group_names.append(_name)
trace_groupings = {
    'All': _all_group_names,
    **{
        key: ['Sun'] + [
            name for name in value
            if name not in {'Sun', 'Full Cluster Catalog'}
        ]
        for key, value in trace_groupings.items()
    },
}
""".strip()
        source, inserted_groupings = re.subn(
            r"(?m)^(plot_3d\s*=\s*Animate3D\()",
            grouping_patch + "\n\n" + r"\1",
            source,
            count=1,
        )
        if inserted_groupings != 1:
            raise RuntimeError("Could not patch trace groupings with Full Cluster Catalog.")

    if "supernova_volumes = _build_supernova_volumes_for_main_figure(time_int)" not in source:
        supernova_block = build_supernova_volume_source_block()
        source, inserted_supernova = re.subn(
            r"(?m)^(fig3d\s*=\s*plot_3d\.make_plot\()",
            supernova_block + "\n\n" + r"\1",
            source,
            count=1,
        )
        if inserted_supernova != 1:
            raise RuntimeError("Could not inject supernova volume helper block.")

    source, replaced_supernova_volumes = re.subn(
        r"volumes\s*=\s*\[\s*edenhofer_volume\s*,\s*mccallum_ne\s*\]",
        "volumes=[edenhofer_volume, *supernova_volumes]",
        source,
        count=1,
    )
    if replaced_supernova_volumes != 1:
        raise RuntimeError("Could not replace notebook volume list with supernova volumes.")

    return source


def run_script_source(source: str) -> None:
    sys.path.insert(0, str(REPO_ROOT))
    globals_dict = {
        "__name__": "__main__",
        "__file__": str(NOTEBOOK_PATH.with_suffix(".py")),
    }
    exec(compile(source, str(NOTEBOOK_PATH.with_suffix(".py")), "exec"), globals_dict)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-html",
        type=Path,
        default=DEFAULT_OUTPUT_HTML,
        help="HTML output path for the rendered figure.",
    )
    parser.add_argument(
        "--ui-design-key",
        default=DEFAULT_UI_DESIGN_KEY,
        help="Three.js UI design preset to force into the figure initial state.",
    )
    parser.add_argument(
        "--theme-key",
        default=None,
        help="Optional Three.js color theme preset to force into the figure initial state.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.output_html.parent.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl")
    os.environ.setdefault("XDG_CACHE_HOME", "/tmp")
    os.environ.setdefault("MPLBACKEND", "Agg")

    with tempfile.TemporaryDirectory(prefix="oviz_main_figure_") as tmp_dir:
        tmp_path = Path(tmp_dir)
        script_path = convert_notebook_to_script(NOTEBOOK_PATH, tmp_path)
        patched_source = patch_script_source(
            script_path.read_text(encoding="utf-8"),
            output_html=args.output_html,
            ui_design_key=args.ui_design_key,
            theme_key=args.theme_key,
        )

    run_script_source(patched_source)
    print(f"Wrote {args.output_html}")


if __name__ == "__main__":
    main()
