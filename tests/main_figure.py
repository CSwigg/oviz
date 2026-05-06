#!/usr/bin/env python3
"""Run the external main_figure notebook as a repeatable smoke test."""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import sys
from pathlib import Path


HOME_DIR = Path.home()
REPO_ROOT = Path(__file__).resolve().parents[1]
NOTEBOOK_PATH = HOME_DIR / "Desktop" / "astro_research" / "radcliffe" / "oviz_notebooks" / "main_figure.ipynb"
SUPERNOVAE_ROOT = HOME_DIR / "Desktop" / "astro_research" / "supernovae_map"
SUPERNOVAE_CATALOG_PATH = SUPERNOVAE_ROOT / "paper" / "solar_encounter_catalog_current.csv.gz"
MCCALLUM_NE_GRID_PATH = HOME_DIR / "Downloads" / "ne_grid.fits"
MIST_PARSEC_AGES_PATH = (
    SUPERNOVAE_ROOT
    / "outputs"
    / "map"
    / "mist_parsec_local1kpc_sfh_compare_200myr"
    / "matched_local1kpc_age_le_200_mist_parsec_clusters.csv"
)
DESKTOP_ROOT = HOME_DIR / "Desktop"
DEFAULT_OUTPUT_HTML = Path(__file__).resolve().with_suffix(".html")
WEBSITE_OUTPUT_HTML = (
    HOME_DIR
    / "Desktop"
    / "astro_research"
    / "cam_website"
    / "interactive_figures"
    / "main_figure.html"
)
GALACTIC_PLANE_IMAGE_PATH = HOME_DIR / "Downloads" / "Top-down_view_of_the_Milky_Way.jpg"
MAIN_FIGURE_DUST_MAX_RESOLUTION = 640
MAIN_FIGURE_DUST_MAX_RESOLUTION_CAP = 640
MAIN_FIGURE_DUST_SAMPLES = 100


def _sanitize_notebook_cell_source(source: str) -> str:
    lines = source.splitlines()
    for line in lines:
        stripped = line.lstrip()
        if not stripped:
            continue
        if stripped.startswith("%%"):
            return "\n".join(f"# {cell_line}" if cell_line else "" for cell_line in lines).strip()
        break

    sanitized_lines = []
    for line in lines:
        stripped = line.lstrip()
        if stripped.startswith("%") or stripped.startswith("!"):
            sanitized_lines.append(f"# {line}")
            continue
        sanitized_lines.append(line)

    return "\n".join(sanitized_lines).strip()


def convert_notebook_to_script_source(notebook_path: Path) -> str:
    notebook = json.loads(notebook_path.read_text(encoding="utf-8"))
    cells = notebook.get("cells", [])
    code_blocks = []

    for cell in cells:
        if cell.get("cell_type") != "code":
            continue
        cell_source = "".join(cell.get("source", []))
        sanitized = _sanitize_notebook_cell_source(cell_source)
        if sanitized:
            code_blocks.append(sanitized)

    if not code_blocks:
        raise RuntimeError(f"No code cells found in notebook: {notebook_path}")

    return "\n\n".join(code_blocks) + "\n"


def patch_edenhofer_volume_integer_setting(
    source: str,
    setting_name: str,
    value: int,
    *,
    insert_after: str | None = None,
) -> str:
    edenhofer_block_pattern = (
        r'(?ms)(edenhofer_volume\s*=\s*\{.*?'
        r'"name"\s*:\s*"Edenhofer\+2024 Dust".*?'
        r'^\s*\})'
    )
    block_match = re.search(edenhofer_block_pattern, source)
    if not block_match:
        raise RuntimeError("Could not find the Edenhofer dust volume block.")

    block = block_match.group(1)
    updated_block, replaced = re.subn(
        rf'(?m)^(\s*"{re.escape(setting_name)}"\s*:\s*)\d+(\s*,)',
        rf"\g<1>{int(value)}\2",
        block,
        count=1,
    )
    if replaced == 0 and insert_after is not None:
        insert_pattern = rf'(?m)^(\s*)"{re.escape(insert_after)}"\s*:\s*\d+\s*,\s*$'

        def insert_setting(match):
            indent = match.group(1)
            return f'{match.group(0)}\n{indent}"{setting_name}": {int(value)},'

        updated_block, replaced = re.subn(
            insert_pattern,
            insert_setting,
            block,
            count=1,
        )
    if replaced != 1:
        raise RuntimeError(f"Could not update the Edenhofer dust {setting_name}.")

    return source[:block_match.start(1)] + updated_block + source[block_match.end(1):]


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
        print(f"Skipping supernova volumes; missing catalog: {{SUPERNOVAE_CATALOG_PATH}}")
        return []

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


def build_mccallum_ne_volume_source_block() -> str:
    return f"""
import os
from pathlib import Path as _OvizPath

MCCALLUM_NE_GRID_PATH = _OvizPath(
    os.environ.get("OVIZ_MCCALLUM_NE_FITS", {str(MCCALLUM_NE_GRID_PATH)!r})
).expanduser()


def _build_mccallum_ne_volumes_for_main_figure():
    if not MCCALLUM_NE_GRID_PATH.exists():
        print(f"Skipping McCallum+2025 electron-density volume; missing FITS cube: {{MCCALLUM_NE_GRID_PATH}}")
        return []

    return [
        {{
            "key": "mccallum-ne",
            "state_key": "mccallum-ne",
            "state_name": "McCallum+2025 electron density",
            "name": "McCallum+2025 Electron Density",
            "path": str(MCCALLUM_NE_GRID_PATH),
            "hdu": "PRIMARY",
            "max_resolution": 384,
            "max_resolution_cap": 384,
            "opacity": 0.12,
            "samples": 100,
            "alpha_coef": 200,
            "gradient_step": 0.006,
            "stretch": "log10",
            "default_vmin_quantile": 0.90,
            "default_vmax_quantile": 0.9995,
            "colormap": "inferno",
            "unit_label": "cm^-3",
            "visible": False,
            "supports_show_all_times": True,
            "co_rotate_with_frame": True,
            "show_all_times": False,
        }}
    ]


mccallum_ne_volumes = _build_mccallum_ne_volumes_for_main_figure()
""".strip()


def build_mist_age_source_block() -> str:
    return f"""
mist_parsec_ages = pd.read_csv({str(MIST_PARSEC_AGES_PATH)!r})
mist_parsec_ages = mist_parsec_ages[
    [
        'name',
        'parsec_age_myr',
        'parsec_initial_mass_msun',
        'mist_age_myr',
        'mist_initial_mass_msun',
    ]
].copy()
for _age_col in ['parsec_age_myr', 'parsec_initial_mass_msun', 'mist_age_myr', 'mist_initial_mass_msun']:
    mist_parsec_ages[_age_col] = pd.to_numeric(mist_parsec_ages[_age_col], errors='coerce')
df_hunt_full = df_hunt_full.merge(
    mist_parsec_ages,
    on='name',
    how='left',
    validate='m:1',
)
df_hunt_full['age_parsec_myr'] = df_hunt_full['parsec_age_myr']
df_hunt_full['age_myr'] = df_hunt_full['mist_age_myr']
df_hunt_full['age_source'] = np.where(df_hunt_full['age_myr'].notnull(), 'MIST', pd.NA)
print(
    'Using MIST ages from {str(MIST_PARSEC_AGES_PATH)}: '
    f"{{df_hunt_full['age_myr'].notnull().sum()}}/{{len(df_hunt_full)}} clusters matched."
)
""".strip()


def patch_script_source(
    source: str,
    output_html: Path,
    theme_key: str | None,
    *,
    minimal_mode: bool = False,
    galactic_simple: bool = False,
    mist_ages: bool = False,
    compact_payload: bool = True,
) -> str:
    source = source.replace("/Users/cam", str(HOME_DIR))
    source = source.replace("/Users/cam/Desktop", str(DESKTOP_ROOT))

    if mist_ages:
        mist_block = build_mist_age_source_block()
        source, replaced_mist_ages = re.subn(
            r"(?m)^df_hunt_full\['age_myr'\]\s*=\s*df_hunt_full\['age_chronos_mode'\]\s*$",
            "df_hunt_full['age_myr'] = df_hunt_full['age_chronos_mode']\n\n" + mist_block,
            source,
            count=1,
        )
        if replaced_mist_ages != 1:
            raise RuntimeError("Could not inject MIST cluster ages into converted notebook script.")

    if galactic_simple:
        source, replaced_time_grid = re.subn(
            r"(?m)^time_int\s*=\s*np\.round\(np\.arange\(0,\s*-66,\s*-1\),\s*1\)\s*$",
            "time_int = np.round(np.arange(0, -66, -1), 1)",
            source,
            count=1,
        )
        if replaced_time_grid != 1:
            raise RuntimeError("Could not restore the galactic simple time grid to 1 Myr steps.")

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
    if galactic_simple:
        control_bits.extend([
            "'camera_fov': 60.0",
            "'point_size_scale': 1.0",
            "'point_glow_strength': 0.60",
            "'fade_in_time_myr': 7.0",
        ])

    if "threejs_initial_state=" not in source:
        initial_state_bits = [
            "'current_group': 'Clusters'",
            "'click_selection_enabled': False",
            f"'compact_payload_enabled': {bool(compact_payload)!r}",
            "'active_volume_key': ('supernova-density' if supernova_volumes else 'volume-0')",
            "'legend_state': ({'volume-0': False, **({'mccallum-ne': False} if mccallum_ne_volumes else {}), 'supernova-density': True} if supernova_volumes else {'volume-0': True, **({'mccallum-ne': False} if mccallum_ne_volumes else {})})",
            "'volume_state_by_key': ({'volume-0': {'visible': False}, **({'mccallum-ne': {'visible': False}} if mccallum_ne_volumes else {}), 'supernova-density': {'visible': True}} if supernova_volumes else {'volume-0': {'visible': True}, **({'mccallum-ne': {'visible': False}} if mccallum_ne_volumes else {})})",
            "'galaxy_image': True",
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
            "{'key': 'P/Mellinger/color', 'label': 'Mellinger Color', 'survey': 'P/Mellinger/color', 'opacity': 1.0, 'visible': True}, "
            "{'key': 'P/PLANCK/R2/HFI/color', 'label': 'Planck Dust Emission Color', 'survey': 'P/PLANCK/R2/HFI/color', 'opacity': 1.0, 'visible': True}"
            "]",
            "'active_sky_layer_key': 'P/Mellinger/color'",
        ]
        if minimal_mode:
            initial_state_bits.append("'minimal_mode_enabled': True")
        if galactic_simple:
            initial_state_bits.extend([
                "'galactic_simple_mode_enabled': True",
                f"'galactic_plane_image_path': {str(GALACTIC_PLANE_IMAGE_PATH)!r}",
                "'galactic_plane_size_pc': 40000.0",
                "'galactic_plane_opacity': 0.6",
                "'initial_zoom_anchor': {'x': 0.0, 'y': 0.0, 'z': 0.0}",
                "'camera': {'position': {'x': 3700.0, 'y': -6550.0, 'z': 4700.0}, 'target': {'x': 3000.0, 'y': 0.0, 'z': 0.0}, 'up': {'x': 0.0, 'y': 0.0, 'z': 1.0}}",
            ])
        if control_bits:
            initial_state_bits.append("'global_controls': {" + ", ".join(control_bits) + "}")
        injection_lines = [
            '    threejs_initial_state={'
            + ", ".join(initial_state_bits)
            + '},\n'
        ]
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

    if minimal_mode:
        source, replaced_young_min_size = re.subn(
            r"(young_trace\s*=\s*Trace\([^\n]*?data_name\s*=\s*'Clusters\s*\(<\s*15\s*Myr\)'[^\n]*?min_size\s*=\s*)([^,]+)",
            r"\g<1>0.0",
            source,
            count=1,
        )
        if replaced_young_min_size != 1:
            raise RuntimeError("Could not adjust young cluster minimum size for minimal mode.")

    if galactic_simple:
        source, recolored_full_sample = re.subn(
            r"(full_sample_trace\s*=\s*Trace\([^\n]*?data_name\s*=\s*'Clusters\s*\(<\s*60\s*Myr\)'[^\n]*?color\s*=\s*)(['\"][^'\"]+['\"])",
            r"\g<1>'#58e1ff'",
            source,
            count=1,
        )
        if recolored_full_sample != 1:
            raise RuntimeError("Could not recolor the <60 Myr cluster trace for galactic simple mode.")

        source, adjusted_full_sample_opacity = re.subn(
            r"(full_sample_trace\s*=\s*Trace\([^\n]*?data_name\s*=\s*'Clusters\s*\(<\s*60\s*Myr\)'[^\n]*?opacity\s*=\s*)([^,\)]+)",
            r"\g<1>0.8",
            source,
            count=1,
        )
        if adjusted_full_sample_opacity != 1:
            raise RuntimeError("Could not adjust the <60 Myr cluster opacity for galactic simple mode.")

        source, removed_young_from_collection = re.subn(
            r"(?m)^(\s*)young_trace,\s*$",
            "",
            source,
            count=1,
        )
        if removed_young_from_collection != 1:
            raise RuntimeError("Could not remove the <15 Myr trace from the galactic simple TraceCollection.")

        source, removed_young_from_grouping = re.subn(
            r'"Clusters"\s*:\s*\[\s*\'Sun\'\s*,\s*\'Clusters\s*\(<\s*60\s*Myr\)\'\s*,\s*\'Clusters\s*\(<\s*15\s*Myr\)\'\s*\]',
            "\"Clusters\": ['Sun', 'Clusters (< 60 Myr)']",
            source,
            count=1,
        )
        if removed_young_from_grouping != 1:
            raise RuntimeError("Could not remove the <15 Myr cluster trace from galactic simple groupings.")

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

    if "mccallum_ne_volumes = _build_mccallum_ne_volumes_for_main_figure()" not in source:
        ne_block = build_mccallum_ne_volume_source_block()
        source, inserted_ne = re.subn(
            r"(?m)^(supernova_volumes\s*=\s*_build_supernova_volumes_for_main_figure\(time_int\)\s*)$",
            r"\1\n\n" + ne_block,
            source,
            count=1,
        )
        if inserted_ne != 1:
            raise RuntimeError("Could not inject McCallum electron-density volume helper block.")

    source, replaced_supernova_volumes = re.subn(
        r"volumes\s*=\s*\[\s*edenhofer_volume\s*,\s*mccallum_ne\s*\]",
        "volumes=[edenhofer_volume, *mccallum_ne_volumes, *supernova_volumes]",
        source,
        count=1,
    )
    if replaced_supernova_volumes != 1:
        raise RuntimeError("Could not replace notebook volume list with supernova volumes.")

    source, edenhofer_time_patch = re.subn(
        r'("colormap"\s*:\s*greys_cmap,\s*# or just "ice"\s*\n)(\s*})',
        r'\1    "supports_show_all_times": True,\n    "co_rotate_with_frame": True,\n\2',
        source,
        count=1,
    )
    if edenhofer_time_patch != 1:
        raise RuntimeError("Could not patch Edenhofer volume frame behavior.")

    return source


def run_script_source(source: str) -> None:
    sys.path.insert(0, str(REPO_ROOT))
    globals_dict = {
        "__name__": "__main__",
        "__file__": str(NOTEBOOK_PATH.with_suffix(".py")),
    }
    exec(compile(source, str(NOTEBOOK_PATH.with_suffix(".py")), "exec"), globals_dict)


def run_main_figure(
    output_html: Path = DEFAULT_OUTPUT_HTML,
    *,
    theme_key: str | None = None,
    minimal_mode: bool = False,
    galactic_simple: bool = False,
    mist_ages: bool = False,
    compact_payload: bool = True,
    website_output_html: Path | None = WEBSITE_OUTPUT_HTML,
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
    patched_source = patch_script_source(
        notebook_source,
        output_html=output_html,
        theme_key=theme_key,
        minimal_mode=minimal_mode,
        galactic_simple=galactic_simple,
        mist_ages=mist_ages,
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
        "--theme-key",
        default=None,
        help="Optional Three.js color theme preset to force into the figure initial state.",
    )
    parser.add_argument(
        "--minimal-mode",
        action="store_true",
        help="Render the figure in minimal presentation mode.",
    )
    parser.add_argument(
        "--mist-ages",
        action="store_true",
        help="Use the matched MIST isochrone ages instead of the default notebook ages.",
    )
    parser.add_argument(
        "--full-payload",
        action="store_true",
        help="Keep repeated per-frame hover, selection, and motion metadata in the HTML.",
    )
    parser.add_argument(
        "--no-website-copy",
        action="store_true",
        help="Do not copy the generated HTML into the cam_website interactive_figures directory.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    output_html = run_main_figure(
        output_html=args.output_html,
        theme_key=args.theme_key,
        minimal_mode=bool(args.minimal_mode),
        mist_ages=bool(args.mist_ages),
        compact_payload=not bool(args.full_payload),
        website_output_html=None if args.no_website_copy else WEBSITE_OUTPUT_HTML,
    )
    print(f"Wrote {output_html}")
    if not args.no_website_copy:
        print(f"Copied {output_html} to {WEBSITE_OUTPUT_HTML}")


if __name__ == "__main__":
    main()
