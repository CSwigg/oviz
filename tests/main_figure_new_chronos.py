#!/usr/bin/env python3
"""Run the external main_figure notebook with newer Chronos cluster ages."""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import sys
import tempfile
from pathlib import Path


HOME_DIR = Path.home()
REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
from oviz.threejs_runtime_ar import THREEJS_AR_QUICKLOOK_SERVICE_WORKER_JS
NOTEBOOK_PATH = HOME_DIR / "Desktop" / "astro_research" / "radcliffe" / "oviz_notebooks" / "main_figure.ipynb"
SUPERNOVAE_ROOT = HOME_DIR / "Desktop" / "astro_research" / "supernovae_map"
SUPERNOVAE_CATALOG_PATH = SUPERNOVAE_ROOT / "paper" / "solar_encounter_catalog_current.csv.gz"
MCCALLUM_NE_GRID_PATH = HOME_DIR / "Downloads" / "ne_grid.fits"
VERGELY_DUST_PATH = HOME_DIR / "Downloads" / "vergely_3D_Dust.fits"
ONEILL_LOCAL_BUBBLE_PATH = HOME_DIR / "Downloads" / "ONeill2024_LocalBubble_Shell_xyz.fits"
GAO_IVS_PEAKS_PATH = HOME_DIR / "Downloads" / "peaks_meanmap_xyz.csv"
GAO_IVS_BOUNDARY_IN_PATH = HOME_DIR / "Downloads" / "npix128boundary_mean_sigma10.csv"
GAO_IVS_BOUNDARY_OUT_PATH = HOME_DIR / "Downloads" / "npix128boundaryout_mean_sigma10.csv"
MIST_PARSEC_AGES_PATH = (
    SUPERNOVAE_ROOT
    / "outputs"
    / "map"
    / "mist_parsec_local1kpc_sfh_compare_200myr"
    / "matched_local1kpc_age_le_200_mist_parsec_clusters.csv"
)
CHRONOS_FASRC_RUN_ROOT = (
    HOME_DIR
    / "Desktop"
    / "astro_research"
    / "chronos_fasrc"
    / "runs"
    / "current"
    / "chronos"
    / "parsec_ybc_huntlt200_xy2kpc_46w_500b_5000s_flatav0to5_12gyr_linearage_192shards"
)
CHRONOS_CLUSTER_RESULTS_PATH = CHRONOS_FASRC_RUN_ROOT / "cluster_results.csv"
CHRONOS_CLUSTER_AGE_STATISTIC = "mode"
DEFAULT_CHRONOS_CLUSTER_MODEL = "parsec"
JUN6_CHRONOS_MODE_AGES_PATH = (
    SUPERNOVAE_ROOT
    / "outputs"
    / "paper_variants"
    / "chronos_parsec_ybc_flatav_posterior_draws_20260605"
    / "parsec"
    / "parsec_chronos_mode_ages.csv"
)
JUN6_CLUSTER_VELOCITIES_PATH = (
    SUPERNOVAE_ROOT
    / "outputs"
    / "velocities"
    / "bulkfit_nomwm_chronos_alt_v2"
    / "cluster_velocities_bulkfit_nomwm.csv"
)
JUN6_MASS_VELOCITY_CATALOG_PATH = (
    SUPERNOVAE_ROOT
    / "outputs"
    / "paper_variants"
    / "chronos_parsec_ybc_flatav_posterior_draws_20260605"
    / "parsec"
    / "parsec_chronos_mf_mass_velocity_catalog_nomwm_bulkfit.csv"
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
MAIN_FIGURE_VOLUME_Z_CLIP_BOUNDS = (-400.0, 400.0)
MAIN_FIGURE_CLUSTER_XY_BOUNDS_PC = (-4000.0, 4000.0)
MAIN_FIGURE_LOOKBACK_MYR = 120
MAIN_FIGURE_TIMESTEP_MYR = 1
MAIN_FIGURE_CLUSTER_YOUNG_AGE_MYR = 15
MAIN_FIGURE_CLUSTER_BLUE_MAX_AGE_MYR = 60
MAIN_FIGURE_CLUSTER_GREY_MAX_AGE_MYR = 150
MAIN_FIGURE_CLUSTER_GREY_COLOR = "#9ca3af"
MAIN_FIGURE_SPIRAL_ARM_TRACE_NAMES = (
    "Spiral Arm: Perseus",
    "Spiral Arm: Local",
    "Spiral Arm: Sagittarius",
    "Spiral Arm: Scutum",
)
MAIN_FIGURE_DUST_MAX_RESOLUTION = 512
MAIN_FIGURE_DUST_MAX_RESOLUTION_CAP = 512
MAIN_FIGURE_DUST_SAMPLES = 200
MAIN_FIGURE_DUST_OPACITY = 0.38
MAIN_FIGURE_DUST_ALPHA_COEF = 105.0
MAIN_FIGURE_MCCALLUM_MAX_RESOLUTION = 512
MAIN_FIGURE_MCCALLUM_MAX_RESOLUTION_CAP = 512
MAIN_FIGURE_MCCALLUM_SAMPLES = 200
MAIN_FIGURE_VERGELY_MAX_RESOLUTION = 512
MAIN_FIGURE_VERGELY_MAX_RESOLUTION_CAP = 512
MAIN_FIGURE_VERGELY_SKY_OVERLAY_MAX_RESOLUTION = 256
MAIN_FIGURE_VERGELY_SAMPLES = 220
MAIN_FIGURE_VERGELY_OPACITY = 0.46
MAIN_FIGURE_VERGELY_ALPHA_COEF = 200.0
MAIN_FIGURE_VERGELY_VMIN = 0.0
MAIN_FIGURE_VERGELY_VMAX = 0.02
MAIN_FIGURE_VERGELY_DEFAULT_VMIN_QUANTILE = 0.70
MAIN_FIGURE_VERGELY_DEFAULT_VMAX_QUANTILE = 0.995
MOBILE_SAFE_TIMESTEP_MYR = 1
MOBILE_SAFE_DUST_MAX_RESOLUTION = 64
MOBILE_SAFE_DUST_MAX_RESOLUTION_CAP = 64
MOBILE_SAFE_DUST_SAMPLES = 40
MOBILE_SAFE_DESKTOP_DUST_MAX_RESOLUTION = MAIN_FIGURE_DUST_MAX_RESOLUTION
MOBILE_SAFE_DESKTOP_DUST_MAX_RESOLUTION_CAP = MAIN_FIGURE_DUST_MAX_RESOLUTION_CAP
MOBILE_SAFE_DESKTOP_DUST_SAMPLES = MAIN_FIGURE_DUST_SAMPLES
MOBILE_SAFE_VERGELY_MAX_RESOLUTION = 64
MOBILE_SAFE_VERGELY_MAX_RESOLUTION_CAP = 64
MOBILE_SAFE_VERGELY_SKY_OVERLAY_MAX_RESOLUTION = 48
MOBILE_SAFE_VERGELY_SAMPLES = 48
MOBILE_SAFE_SKY_DOME_CAPTURE_WIDTH_PX = 640
MOBILE_SAFE_SKY_DOME_CAPTURE_HEIGHT_PX = 320
MOBILE_SAFE_SKY_DOME_CAPTURE_QUALITY = 0.70
MOBILE_SAFE_FULL_CATALOG_MAX_POINTS = 500
MOBILE_SAFE_BACKGROUND_CLUSTER_MAX_POINTS = 500
MOBILE_SAFE_BLUE_CLUSTER_MAX_POINTS = 650
MOBILE_SAFE_GALAXY_IMAGE_MAX_PX = 1024
MOBILE_SAFE_GALAXY_IMAGE_QUALITY = 58
AR_QUICKLOOK_SERVICE_WORKER_NAME = "oviz-ar-quicklook-sw.js"


def write_ar_quicklook_service_worker(output_dir: Path) -> Path:
    worker_path = Path(output_dir) / AR_QUICKLOOK_SERVICE_WORKER_NAME
    worker_path.write_text(THREEJS_AR_QUICKLOOK_SERVICE_WORKER_JS + "\n", encoding="utf-8")
    return worker_path


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


def mobile_safe_galaxy_image_path(source_path: Path = GALACTIC_PLANE_IMAGE_PATH) -> Path:
    """Return a smaller copy of the Milky Way image for upload-safe HTML exports."""
    source_path = Path(source_path).expanduser()
    if not source_path.exists():
        return source_path

    target_path = (
        Path(tempfile.gettempdir())
        / f"oviz_mobile_safe_galaxy_{MOBILE_SAFE_GALAXY_IMAGE_MAX_PX}px_q{MOBILE_SAFE_GALAXY_IMAGE_QUALITY}.jpg"
    )
    if target_path.exists() and target_path.stat().st_mtime >= source_path.stat().st_mtime:
        return target_path

    try:
        from PIL import Image, ImageOps
    except Exception:
        return source_path

    with Image.open(source_path) as image:
        image = ImageOps.exif_transpose(image).convert("RGB")
        resampling = getattr(getattr(Image, "Resampling", Image), "LANCZOS")
        image.thumbnail(
            (MOBILE_SAFE_GALAXY_IMAGE_MAX_PX, MOBILE_SAFE_GALAXY_IMAGE_MAX_PX),
            resampling,
        )
        image.save(
            target_path,
            format="JPEG",
            quality=MOBILE_SAFE_GALAXY_IMAGE_QUALITY,
            optimize=True,
            progressive=True,
        )
    return target_path


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


def patch_edenhofer_volume_numeric_setting(
    source: str,
    setting_name: str,
    value: float,
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
    numeric_literal = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"
    updated_block, replaced = re.subn(
        rf'(?m)^(\s*"{re.escape(setting_name)}"\s*:\s*){numeric_literal}(\s*,)',
        rf"\g<1>{float(value)!r}\2",
        block,
        count=1,
    )
    if replaced != 1:
        raise RuntimeError(f"Could not update the Edenhofer dust {setting_name}.")

    return source[:block_match.start(1)] + updated_block + source[block_match.end(1):]


def patch_edenhofer_volume_clip_bounds(source: str, z_bounds: tuple[float, float]) -> str:
    edenhofer_block_pattern = (
        r'(?ms)(edenhofer_volume\s*=\s*\{.*?'
        r'"name"\s*:\s*"Edenhofer\+2024 Dust".*?'
        r'^\s*\})'
    )
    block_match = re.search(edenhofer_block_pattern, source)
    if not block_match:
        raise RuntimeError("Could not find the Edenhofer dust volume block.")

    block = block_match.group(1)
    clip_text = f'"clip_bounds": {{"z": [{float(z_bounds[0])}, {float(z_bounds[1])}]}},'
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


def patch_base_cluster_filter(source: str) -> str:
    replacement = f"""
df_hunt_good = df_hunt_full.loc[
    (df_hunt_full['U_err'] < 10) &
    (df_hunt_full['V_err'] < 10) &
    (df_hunt_full['W_err'] < 10) &
    (df_hunt_full['U'].notnull()) &
    (df_hunt_full['V'].notnull()) &
    (df_hunt_full['W'].notnull()) &
    (df_hunt_full['x'].notnull()) &
    (df_hunt_full['y'].notnull()) &
    (df_hunt_full['z'].notnull()) &
    (df_hunt_full['age_myr'].notnull()) &
    (df_hunt_full['x'].between({MAIN_FIGURE_CLUSTER_XY_BOUNDS_PC[0]!r}, {MAIN_FIGURE_CLUSTER_XY_BOUNDS_PC[1]!r})) &
    (df_hunt_full['y'].between({MAIN_FIGURE_CLUSTER_XY_BOUNDS_PC[0]!r}, {MAIN_FIGURE_CLUSTER_XY_BOUNDS_PC[1]!r})) &
    (df_hunt_full['n_rvs_2026'] >= 3)
]
""".strip()
    source, replaced = re.subn(
        r"(?ms)^df_hunt_good\s*=\s*df_hunt_full\.loc\[\s*\n.*?^\s*\]\s*$",
        replacement,
        source,
        count=1,
    )
    if replaced != 1:
        raise RuntimeError("Could not rewrite the notebook base cluster filter.")
    return source


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
            "clip_bounds": {{"z": [{float(MAIN_FIGURE_VOLUME_Z_CLIP_BOUNDS[0])}, {float(MAIN_FIGURE_VOLUME_Z_CLIP_BOUNDS[1])}]}},
            "max_resolution": {int(MAIN_FIGURE_MCCALLUM_MAX_RESOLUTION)},
            "max_resolution_cap": {int(MAIN_FIGURE_MCCALLUM_MAX_RESOLUTION_CAP)},
            "opacity": 0.12,
            "samples": {int(MAIN_FIGURE_MCCALLUM_SAMPLES)},
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


def build_vergely_dust_volume_source_block(*, mobile_safe_mode: bool = False) -> str:
    max_resolution = MOBILE_SAFE_VERGELY_MAX_RESOLUTION if mobile_safe_mode else MAIN_FIGURE_VERGELY_MAX_RESOLUTION
    max_resolution_cap = (
        MOBILE_SAFE_VERGELY_MAX_RESOLUTION_CAP
        if mobile_safe_mode
        else MAIN_FIGURE_VERGELY_MAX_RESOLUTION_CAP
    )
    sky_overlay_max_resolution = (
        MOBILE_SAFE_VERGELY_SKY_OVERLAY_MAX_RESOLUTION
        if mobile_safe_mode
        else MAIN_FIGURE_VERGELY_SKY_OVERLAY_MAX_RESOLUTION
    )
    samples = MOBILE_SAFE_VERGELY_SAMPLES if mobile_safe_mode else MAIN_FIGURE_VERGELY_SAMPLES
    return f"""
import os
from pathlib import Path as _OvizPath

VERGELY_DUST_PATH = _OvizPath(
    os.environ.get("OVIZ_VERGELY_DUST_FITS", {str(VERGELY_DUST_PATH)!r})
).expanduser()


def _build_vergely_dust_volumes_for_main_figure():
    if not VERGELY_DUST_PATH.exists():
        print(f"Skipping Vergely 3D dust volume; missing FITS cube: {{VERGELY_DUST_PATH}}")
        return []

    return [
        {{
            "key": "vergely-dust",
            "state_key": "vergely-dust",
            "state_name": "Vergely 3D Dust",
            "name": "Vergely 3D Dust",
            "path": str(VERGELY_DUST_PATH),
            "hdu": "PRIMARY",
            "clip_bounds": {{"z": [{float(MAIN_FIGURE_VOLUME_Z_CLIP_BOUNDS[0])}, {float(MAIN_FIGURE_VOLUME_Z_CLIP_BOUNDS[1])}]}},
            "max_resolution": {int(max_resolution)},
            "max_resolution_cap": {int(max_resolution_cap)},
            "sky_overlay_max_resolution": {int(sky_overlay_max_resolution)},
            "data_encoding": "png_atlas_uint8",
            "opacity": {float(MAIN_FIGURE_VERGELY_OPACITY)!r},
            "samples": {int(samples)},
            "alpha_coef": {float(MAIN_FIGURE_VERGELY_ALPHA_COEF)!r},
            "gradient_step": 0.006,
            "stretch": "asinh",
            "vmin": {float(MAIN_FIGURE_VERGELY_VMIN)!r},
            "vmax": {float(MAIN_FIGURE_VERGELY_VMAX)!r},
            "default_vmin_quantile": {float(MAIN_FIGURE_VERGELY_DEFAULT_VMIN_QUANTILE)!r},
            "default_vmax_quantile": {float(MAIN_FIGURE_VERGELY_DEFAULT_VMAX_QUANTILE)!r},
            "colormap": "magma",
            "unit_label": "mag pc^-1",
            "visible": True,
            "supports_show_all_times": True,
            "co_rotate_with_frame": True,
            "show_all_times": False,
        }}
    ]


vergely_dust_volumes = _build_vergely_dust_volumes_for_main_figure()
""".strip()


def build_adaptive_edenhofer_volume_source_block() -> str:
    return f"""
edenhofer_volume_desktop = dict(edenhofer_volume)
edenhofer_volume_desktop.update({{
    "key": "edenhofer-dust-desktop",
    "state_key": "edenhofer-dust-desktop",
    "state_name": "Edenhofer+2024 Dust",
    "base_state_name": "Edenhofer+2024 Dust",
    "name": "Edenhofer+2024 Dust (desktop)",
    "variant_group": "edenhofer-dust-resolution",
    "variant_label": "Desktop high-res",
    "variant_order": 0,
    "max_resolution": {int(MOBILE_SAFE_DESKTOP_DUST_MAX_RESOLUTION)},
    "max_resolution_cap": {int(MOBILE_SAFE_DESKTOP_DUST_MAX_RESOLUTION_CAP)},
    "sky_overlay_max_resolution": {int(MOBILE_SAFE_DESKTOP_DUST_MAX_RESOLUTION)},
    "samples": {int(MOBILE_SAFE_DESKTOP_DUST_SAMPLES)},
    "visible": True,
}})
edenhofer_volume_mobile = dict(edenhofer_volume)
edenhofer_volume_mobile.update({{
    "key": "edenhofer-dust-mobile",
    "state_key": "edenhofer-dust-mobile",
    "state_name": "Edenhofer+2024 Dust",
    "base_state_name": "Edenhofer+2024 Dust",
    "name": "Edenhofer+2024 Dust (mobile)",
    "variant_group": "edenhofer-dust-resolution",
    "variant_label": "Mobile low-res",
    "variant_order": 1,
    "max_resolution": {int(MOBILE_SAFE_DUST_MAX_RESOLUTION)},
    "max_resolution_cap": {int(MOBILE_SAFE_DUST_MAX_RESOLUTION_CAP)},
    "sky_overlay_max_resolution": {int(MOBILE_SAFE_DUST_MAX_RESOLUTION)},
    "data_encoding": "png_atlas_uint8",
    "samples": {int(MOBILE_SAFE_DUST_SAMPLES)},
    "visible": False,
}})
edenhofer_volumes = [edenhofer_volume_desktop, edenhofer_volume_mobile]
""".strip()


def build_local_shell_volume_source_block() -> str:
    return f"""
import os
from pathlib import Path as _OvizPath

import numpy as _oviz_np
import pandas as _oviz_pd
from matplotlib.colors import LinearSegmentedColormap as _OvizLinearSegmentedColormap

ONEILL_LOCAL_BUBBLE_PATH = _OvizPath(
    os.environ.get("OVIZ_ONEILL_LOCAL_BUBBLE_FITS", {str(ONEILL_LOCAL_BUBBLE_PATH)!r})
).expanduser()
GAO_IVS_PEAKS_PATH = _OvizPath(
    os.environ.get("OVIZ_GAO_IVS_PEAKS_CSV", {str(GAO_IVS_PEAKS_PATH)!r})
).expanduser()
GAO_IVS_BOUNDARY_IN_PATH = _OvizPath(
    os.environ.get("OVIZ_GAO_IVS_BOUNDARY_IN_CSV", {str(GAO_IVS_BOUNDARY_IN_PATH)!r})
).expanduser()
GAO_IVS_BOUNDARY_OUT_PATH = _OvizPath(
    os.environ.get("OVIZ_GAO_IVS_BOUNDARY_OUT_CSV", {str(GAO_IVS_BOUNDARY_OUT_PATH)!r})
).expanduser()

_DODGER_BLUE_VOLUME_CMAP = _OvizLinearSegmentedColormap.from_list(
    "dodger_blue",
    [
        (0.0, (0.0, 0.0, 0.0, 0.0)),
        (0.18, (0.1176470588, 0.5647058824, 1.0, 0.18)),
        (1.0, (0.1176470588, 0.5647058824, 1.0, 1.0)),
    ],
)


def _build_gao_ivs_volume_data(grid_size=192):
    required_paths = [
        GAO_IVS_PEAKS_PATH,
        GAO_IVS_BOUNDARY_IN_PATH,
        GAO_IVS_BOUNDARY_OUT_PATH,
    ]
    missing_paths = [path for path in required_paths if not path.exists()]
    if missing_paths:
        print(
            "Skipping Gao+2024 IVS volume; missing model table(s): "
            + ", ".join(str(path) for path in missing_paths)
        )
        return None

    peaks = _oviz_pd.read_csv(
        GAO_IVS_PEAKS_PATH,
        usecols=["index", "x[pc]", "y[pc]", "z[pc]", "nH[cm-3]"],
    )
    boundary_in = _oviz_pd.read_csv(
        GAO_IVS_BOUNDARY_IN_PATH,
        usecols=["index", "x[pc]", "y[pc]", "z[pc]"],
    )
    boundary_out = _oviz_pd.read_csv(
        GAO_IVS_BOUNDARY_OUT_PATH,
        usecols=["index", "x[pc]", "y[pc]", "z[pc]"],
    )

    peaks["index"] = peaks["index"].astype(int)
    boundary_in["index"] = boundary_in["index"].astype(int)
    boundary_out["index"] = boundary_out["index"].astype(int)
    shell = (
        boundary_in.rename(
            columns={{"x[pc]": "xin", "y[pc]": "yin", "z[pc]": "zin"}}
        )
        .merge(
            boundary_out.rename(
                columns={{"x[pc]": "xout", "y[pc]": "yout", "z[pc]": "zout"}}
            ),
            on="index",
            how="inner",
            validate="1:1",
        )
        .merge(
            peaks.rename(
                columns={{"x[pc]": "xpeak", "y[pc]": "ypeak", "z[pc]": "zpeak"}}
            )[["index", "xpeak", "ypeak", "zpeak", "nH[cm-3]"]],
            on="index",
            how="left",
            validate="1:1",
        )
    )
    if shell.empty:
        print("Skipping Gao+2024 IVS volume; no matched boundary/peak rows.")
        return None

    for col in ["xin", "yin", "zin", "xout", "yout", "zout", "xpeak", "ypeak", "zpeak", "nH[cm-3]"]:
        shell[col] = _oviz_pd.to_numeric(shell[col], errors="coerce")
    shell = shell.dropna(subset=["xin", "yin", "zin", "xout", "yout", "zout"])
    if shell.empty:
        print("Skipping Gao+2024 IVS volume; no finite boundary coordinates.")
        return None

    fractions = _oviz_np.linspace(0.0, 1.0, 11, dtype=_oviz_np.float32)
    xin = shell["xin"].to_numpy(dtype=_oviz_np.float32)
    yin = shell["yin"].to_numpy(dtype=_oviz_np.float32)
    zin = shell["zin"].to_numpy(dtype=_oviz_np.float32)
    xout = shell["xout"].to_numpy(dtype=_oviz_np.float32)
    yout = shell["yout"].to_numpy(dtype=_oviz_np.float32)
    zout = shell["zout"].to_numpy(dtype=_oviz_np.float32)
    x_boundary = (xin[:, None] + (xout - xin)[:, None] * fractions[None, :]).ravel()
    y_boundary = (yin[:, None] + (yout - yin)[:, None] * fractions[None, :]).ravel()
    z_boundary = (zin[:, None] + (zout - zin)[:, None] * fractions[None, :]).ravel()

    density = shell["nH[cm-3]"].to_numpy(dtype=_oviz_np.float32)
    if _oviz_np.isfinite(density).any():
        density = _oviz_np.nan_to_num(density, nan=float(_oviz_np.nanmedian(density)), posinf=0.0, neginf=0.0)
        base_weights = _oviz_np.log1p(_oviz_np.clip(density, 0.0, None))
        boundary_weights = 0.65 * _oviz_np.repeat(base_weights, len(fractions))
    else:
        base_weights = None
        boundary_weights = None

    xpeak = shell["xpeak"].to_numpy(dtype=_oviz_np.float32)
    ypeak = shell["ypeak"].to_numpy(dtype=_oviz_np.float32)
    zpeak = shell["zpeak"].to_numpy(dtype=_oviz_np.float32)
    peak_mask = _oviz_np.isfinite(xpeak) & _oviz_np.isfinite(ypeak) & _oviz_np.isfinite(zpeak)
    if peak_mask.any():
        x_points = _oviz_np.concatenate([x_boundary, xpeak[peak_mask]])
        y_points = _oviz_np.concatenate([y_boundary, ypeak[peak_mask]])
        z_points = _oviz_np.concatenate([z_boundary, zpeak[peak_mask]])
        if base_weights is None:
            weights = None
        else:
            weights = _oviz_np.concatenate([boundary_weights, 2.5 * base_weights[peak_mask]])
    else:
        x_points = x_boundary
        y_points = y_boundary
        z_points = z_boundary
        weights = boundary_weights

    margin_pc = 18.0
    x_bounds = [float(_oviz_np.nanmin(x_points) - margin_pc), float(_oviz_np.nanmax(x_points) + margin_pc)]
    y_bounds = [float(_oviz_np.nanmin(y_points) - margin_pc), float(_oviz_np.nanmax(y_points) + margin_pc)]
    z_bounds = [float(_oviz_np.nanmin(z_points) - margin_pc), float(_oviz_np.nanmax(z_points) + margin_pc)]
    x_edges = _oviz_np.linspace(x_bounds[0], x_bounds[1], int(grid_size) + 1, dtype=_oviz_np.float32)
    y_edges = _oviz_np.linspace(y_bounds[0], y_bounds[1], int(grid_size) + 1, dtype=_oviz_np.float32)
    z_edges = _oviz_np.linspace(z_bounds[0], z_bounds[1], int(grid_size) + 1, dtype=_oviz_np.float32)

    cube_zyx, _ = _oviz_np.histogramdd(
        _oviz_np.column_stack((z_points, y_points, x_points)),
        bins=(z_edges, y_edges, x_edges),
        weights=weights,
    )
    cube_zyx = cube_zyx.astype(_oviz_np.float32, copy=False)

    try:
        from scipy.ndimage import gaussian_filter as _oviz_gaussian_filter
    except ImportError:
        _oviz_gaussian_filter = None
    if _oviz_gaussian_filter is not None:
        cube_zyx = _oviz_gaussian_filter(
            cube_zyx,
            sigma=0.55,
            mode="constant",
        ).astype(_oviz_np.float32, copy=False)
    # Compress the dynamic range before compact uint8 embedding so faint IVS
    # filaments survive instead of being quantized away by the densest voxels.
    cube_zyx = _oviz_np.log1p(_oviz_np.clip(cube_zyx, 0.0, None)).astype(
        _oviz_np.float32,
        copy=False,
    )

    positive = cube_zyx[cube_zyx > 0]
    if not positive.size:
        print("Skipping Gao+2024 IVS volume; voxelization produced an empty cube.")
        return None

    return {{
        "data": cube_zyx,
        "bounds": {{
            "x": x_bounds,
            "y": y_bounds,
            "z": z_bounds,
        }},
    }}


def _build_local_shell_volumes_for_main_figure():
    volumes = []

    if ONEILL_LOCAL_BUBBLE_PATH.exists():
        volumes.append(
            {{
                "key": "oneill-local-bubble-shell",
                "state_key": "oneill-local-bubble-shell",
                "state_name": "O'Neill+2024 Local Bubble",
                "name": "O'Neill+2024 Local Bubble Shell",
                "path": str(ONEILL_LOCAL_BUBBLE_PATH),
                "hdu": "SHELL",
                "time_myr": 0.0,
                "max_resolution": 128,
                "max_resolution_cap": 128,
                "sky_overlay_max_resolution": 128,
                "opacity": 0.22,
                "samples": 120,
                "alpha_coef": 105,
                "gradient_step": 0.006,
                "stretch": "asinh",
                "default_vmin_quantile": 0.86,
                "default_vmax_quantile": 0.999,
                "colormap": _DODGER_BLUE_VOLUME_CMAP,
                "unit_label": "shell density",
                "visible": False,
                "only_at_t0": True,
                "supports_show_all_times": False,
                "show_all_times": False,
            }}
        )
    else:
        print(f"Skipping O'Neill+2024 Local Bubble volume; missing FITS cube: {{ONEILL_LOCAL_BUBBLE_PATH}}")

    ivs_volume_data = _build_gao_ivs_volume_data(grid_size=192)
    if ivs_volume_data is not None:
        volumes.append(
            {{
                "key": "gao-ivs-shell",
                "state_key": "gao-ivs-shell",
                "state_name": "Gao+2024 IRAS Vela Shell",
                "name": "Gao+2024 IRAS Vela Shell",
                "time_myr": 0.0,
                "data": ivs_volume_data["data"],
                "bounds": ivs_volume_data["bounds"],
                "apply_center_offset": True,
                "opacity": 0.26,
                "samples": 192,
                "alpha_coef": 125,
                "gradient_step": 0.006,
                "stretch": "asinh",
                "default_vmin_quantile": 0.45,
                "default_vmax_quantile": 0.995,
                "colormap": _DODGER_BLUE_VOLUME_CMAP,
                "unit_label": "weighted shell density",
                "visible": False,
                "only_at_t0": True,
                "supports_show_all_times": False,
                "show_all_times": False,
            }}
        )

    return volumes


local_shell_volumes = _build_local_shell_volumes_for_main_figure()
optional_static_volume_state = {{
    str(_volume.get("state_key") or _volume.get("key")): {{"visible": False}}
    for _volume in [*mccallum_ne_volumes, *local_shell_volumes]
}}
optional_static_legend_state = {{
    _state_key: False
    for _state_key in optional_static_volume_state
}}
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


def build_jun6_catalog_source_block() -> str:
    velocity_cols = [
        "x_2026",
        "y_2026",
        "z_2026",
        "x_err_2026",
        "y_err_2026",
        "z_err_2026",
        "U_2026",
        "V_2026",
        "W_2026",
        "U_err_2026",
        "V_err_2026",
        "W_err_2026",
        "n_rvs_2026",
        "velocity_fit_status",
        "velocity_fit_quality",
        "velocity_fit_method",
    ]
    return f"""
df_hunt_full = pd.read_csv({str(JUN6_MASS_VELOCITY_CATALOG_PATH)!r})
df_hunt_full['hunt_age_myr'] = pd.to_numeric(
    df_hunt_full['age_myr'] if 'age_myr' in df_hunt_full.columns else pd.Series(np.nan, index=df_hunt_full.index),
    errors='coerce',
)
df_hunt_full['hunt_initial_mass_msun'] = pd.to_numeric(
    df_hunt_full['mass_all_previous'] if 'mass_all_previous' in df_hunt_full.columns else (
        df_hunt_full['mass_all'] if 'mass_all' in df_hunt_full.columns else pd.Series(np.nan, index=df_hunt_full.index)
    ),
    errors='coerce',
)
df_hunt_full['chronos_initial_mass_msun'] = pd.to_numeric(
    df_hunt_full['mass_all'] if 'mass_all' in df_hunt_full.columns else pd.Series(np.nan, index=df_hunt_full.index),
    errors='coerce',
)
df_hunt_full['display_initial_mass_msun'] = df_hunt_full['chronos_initial_mass_msun'].combine_first(
    df_hunt_full['hunt_initial_mass_msun']
)
df_hunt_full['mass_all'] = df_hunt_full['display_initial_mass_msun']
df_hunt_full['mass_source'] = np.where(
    df_hunt_full['chronos_initial_mass_msun'].notnull(),
    'Chronos',
    np.where(df_hunt_full['hunt_initial_mass_msun'].notnull(), 'Hunt', pd.NA),
)
jun6_velocity_cols = ['name', *{velocity_cols!r}]
jun6_velocity_source = pd.read_csv(
    {str(JUN6_CLUSTER_VELOCITIES_PATH)!r},
    usecols=lambda _col: _col in jun6_velocity_cols,
)
if jun6_velocity_source['name'].duplicated().any():
    raise RuntimeError('Jun 6 velocity source has duplicate cluster names.')
df_hunt_full = df_hunt_full.drop(
    columns=[_col for _col in jun6_velocity_cols if _col != 'name' and _col in df_hunt_full.columns],
    errors='ignore',
)
df_hunt_full = df_hunt_full.merge(
    jun6_velocity_source,
    on='name',
    how='left',
    validate='1:1',
)
ages_chronos = pd.read_csv({str(JUN6_CHRONOS_MODE_AGES_PATH)!r})
if ages_chronos['name'].duplicated().any():
    raise RuntimeError('Jun 6 Chronos mode age source has duplicate cluster names.')
df_hunt_full = df_hunt_full.drop(
    columns=['age_chronos_lo', 'age_chronos_mode', 'age_chronos_hi'],
    errors='ignore',
)
df_hunt_full = df_hunt_full.merge(
    ages_chronos[['name', 'age_chronos_lo', 'age_chronos_mode', 'age_chronos_hi']],
    on='name',
    how='left',
    validate='1:1',
)
for _jun6_col in [
    'age_chronos_lo',
    'age_chronos_mode',
    'age_chronos_hi',
    'x_2026',
    'y_2026',
    'z_2026',
    'U_2026',
    'V_2026',
    'W_2026',
    'U_err_2026',
    'V_err_2026',
    'W_err_2026',
    'n_rvs_2026',
]:
    if _jun6_col in df_hunt_full.columns:
        df_hunt_full[_jun6_col] = pd.to_numeric(df_hunt_full[_jun6_col], errors='coerce')
df_hunt_full['age_myr'] = df_hunt_full['age_chronos_mode'].combine_first(df_hunt_full['hunt_age_myr'])
df_hunt_full['age_source'] = np.where(
    df_hunt_full['age_chronos_mode'].notnull(),
    'Chronos',
    np.where(df_hunt_full['hunt_age_myr'].notnull(), 'Hunt', pd.NA),
)
df_hunt_full = df_hunt_full.drop(
    columns=['x', 'y', 'z', 'U', 'V', 'W', 'U_err', 'V_err', 'W_err'],
    errors='ignore',
)
df_hunt_full = df_hunt_full.rename(columns={{
    'U_2026': 'U',
    'V_2026': 'V',
    'W_2026': 'W',
    'U_err_2026': 'U_err',
    'V_err_2026': 'V_err',
    'W_err_2026': 'W_err',
    'x_2026': 'x',
    'y_2026': 'y',
    'z_2026': 'z',
}})
print(
    'Using Jun 6 cluster catalog: '
    f"{{df_hunt_full['age_myr'].notnull().sum()}} clusters with Chronos mode ages, "
    f"{{df_hunt_full[['x', 'y', 'z', 'U', 'V', 'W']].notnull().all(axis=1).sum()}} with complete positions/velocities."
)
""".strip()


def build_chronos_cluster_sample_source_block(
    chronos_results_path: Path,
    chronos_model: str,
) -> str:
    model_key = chronos_model.lower()
    if not re.fullmatch(r"[a-z0-9_]+", model_key):
        raise ValueError(f"Invalid Chronos model name: {chronos_model!r}")

    model_display = model_key.upper()
    age_lo_col = f"{model_key}_age_lo"
    age_value_col = f"{model_key}_age_{CHRONOS_CLUSTER_AGE_STATISTIC}"
    age_hi_col = f"{model_key}_age_hi"
    mass_col = f"{model_key}_mass_cluster_imf_corrected"
    status_col = f"{model_key}_status"

    return f"""
chronos_cluster_columns = set(pd.read_csv(
    {str(chronos_results_path)!r},
    nrows=0,
).columns)
chronos_model_key = {model_key!r}
chronos_model_display = {model_display!r}
chronos_age_value_col = {age_value_col!r}
chronos_required_cols = [
    'name',
    {age_lo_col!r},
    chronos_age_value_col,
    {age_hi_col!r},
    {mass_col!r},
    {status_col!r},
]
chronos_missing_cols = [col for col in chronos_required_cols if col not in chronos_cluster_columns]
if chronos_missing_cols:
    raise RuntimeError(
        f"Could not find {{chronos_missing_cols!r}} in {str(chronos_results_path)!r}."
    )
chronos_cluster_ages = pd.read_csv(
    {str(chronos_results_path)!r},
    usecols=chronos_required_cols,
)
chronos_cluster_ages = chronos_cluster_ages.rename(
    columns={{
        {age_lo_col!r}: 'chronos_age_lo_myr',
        chronos_age_value_col: 'chronos_age_myr',
        {age_hi_col!r}: 'chronos_age_hi_myr',
        {mass_col!r}: 'chronos_initial_mass_msun',
        {status_col!r}: 'chronos_status',
    }}
)
for _chronos_col in [
    'chronos_age_lo_myr',
    'chronos_age_myr',
    'chronos_age_hi_myr',
    'chronos_initial_mass_msun',
]:
    chronos_cluster_ages[_chronos_col] = pd.to_numeric(
        chronos_cluster_ages[_chronos_col],
        errors='coerce',
    )
chronos_cluster_ages = chronos_cluster_ages.loc[
    chronos_cluster_ages['chronos_status'].eq('success')
    & chronos_cluster_ages['chronos_age_myr'].notnull()
].copy()
df_hunt_chronos_full = df_hunt_full.drop(
    columns=[
        'chronos_age_lo_myr',
        'chronos_age_myr',
        'chronos_age_hi_myr',
        'chronos_initial_mass_msun',
        'chronos_status',
    ],
    errors='ignore',
).merge(
    chronos_cluster_ages,
    on='name',
    how='left',
    validate='m:1',
)
if 'hunt_age_myr' not in df_hunt_chronos_full.columns:
    df_hunt_chronos_full['hunt_age_myr'] = pd.to_numeric(
        df_hunt_chronos_full['age_myr'] if 'age_myr' in df_hunt_chronos_full.columns else pd.Series(np.nan, index=df_hunt_chronos_full.index),
        errors='coerce',
    )
if 'hunt_initial_mass_msun' not in df_hunt_chronos_full.columns:
    df_hunt_chronos_full['hunt_initial_mass_msun'] = pd.to_numeric(
        df_hunt_chronos_full['mass_all_previous'] if 'mass_all_previous' in df_hunt_chronos_full.columns else (
            df_hunt_chronos_full['mass_all'] if 'mass_all' in df_hunt_chronos_full.columns else pd.Series(np.nan, index=df_hunt_chronos_full.index)
        ),
        errors='coerce',
    )
df_hunt_chronos_full['display_age_myr'] = df_hunt_chronos_full['chronos_age_myr'].combine_first(
    df_hunt_chronos_full['hunt_age_myr']
)
df_hunt_chronos_full['chronos_initial_mass_msun'] = df_hunt_chronos_full['chronos_initial_mass_msun'].combine_first(
    df_hunt_chronos_full['hunt_initial_mass_msun']
)
df_hunt_chronos_full['display_initial_mass_msun'] = df_hunt_chronos_full['chronos_initial_mass_msun']
df_hunt_chronos_full['mass_all'] = df_hunt_chronos_full['display_initial_mass_msun']
main_figure_cluster_xy_min_pc = {MAIN_FIGURE_CLUSTER_XY_BOUNDS_PC[0]!r}
main_figure_cluster_xy_max_pc = {MAIN_FIGURE_CLUSTER_XY_BOUNDS_PC[1]!r}
main_figure_cluster_young_age_myr = {MAIN_FIGURE_CLUSTER_YOUNG_AGE_MYR!r}
main_figure_cluster_blue_max_age_myr = {MAIN_FIGURE_CLUSTER_BLUE_MAX_AGE_MYR!r}
main_figure_cluster_grey_max_age_myr = {MAIN_FIGURE_CLUSTER_GREY_MAX_AGE_MYR!r}
df_hunt_chronos_sample = df_hunt_chronos_full.loc[
    (df_hunt_chronos_full['U_err'] < 10) &
    (df_hunt_chronos_full['V_err'] < 10) &
    (df_hunt_chronos_full['W_err'] < 10) &
    (df_hunt_chronos_full['U'].notnull()) &
    (df_hunt_chronos_full['V'].notnull()) &
    (df_hunt_chronos_full['W'].notnull()) &
    (df_hunt_chronos_full['x'].notnull()) &
    (df_hunt_chronos_full['y'].notnull()) &
    (df_hunt_chronos_full['z'].notnull()) &
    (df_hunt_chronos_full['x'].between(main_figure_cluster_xy_min_pc, main_figure_cluster_xy_max_pc)) &
    (df_hunt_chronos_full['y'].between(main_figure_cluster_xy_min_pc, main_figure_cluster_xy_max_pc)) &
    (df_hunt_chronos_full['display_age_myr'].notnull()) &
    (df_hunt_chronos_full['n_rvs_2026'] >= 3)
].copy()
df_hunt_chronos_sample['age_myr'] = df_hunt_chronos_sample['display_age_myr']
df_hunt_chronos_sample['age_source'] = np.where(
    df_hunt_chronos_sample['chronos_age_myr'].notnull(),
    chronos_model_display,
    'Hunt',
)
df_hunt_60_chronos = df_hunt_chronos_sample.loc[
    df_hunt_chronos_sample['age_myr'] < main_figure_cluster_blue_max_age_myr
].copy()
df_hunt_0_to_150_chronos = df_hunt_chronos_sample.loc[
    df_hunt_chronos_sample['age_myr'] < main_figure_cluster_grey_max_age_myr
].copy()
df_hunt_young_chronos = df_hunt_chronos_sample.loc[
    df_hunt_chronos_sample['age_myr'] < main_figure_cluster_young_age_myr
].copy()
print(
    f"Using {{chronos_model_display}} ages from {str(chronos_results_path)} for displayed cluster samples: "
    f"{{len(df_hunt_60_chronos)}} clusters <{{main_figure_cluster_blue_max_age_myr:g}} Myr, "
    f"{{len(df_hunt_0_to_150_chronos)}} clusters 0-{{main_figure_cluster_grey_max_age_myr:g}} Myr, "
    f"{{len(df_hunt_young_chronos)}} clusters <{{main_figure_cluster_young_age_myr:g}} Myr. "
    f"Displayed ages use {{chronos_age_value_col}} with Hunt fallbacks when Chronos is missing. "
    f"x/y within [{{main_figure_cluster_xy_min_pc:g}}, {{main_figure_cluster_xy_max_pc:g}}] pc; no z cut applied."
)
""".strip()


def build_chronos_mode_age_sample_source_block(
    chronos_mode_ages_path: Path,
    *,
    model_display: str = "PARSEC",
) -> str:
    return f"""
chronos_mode_age_source = pd.read_csv({str(chronos_mode_ages_path)!r})
chronos_required_cols = ['name', 'age_chronos_lo', 'age_chronos_mode', 'age_chronos_hi']
chronos_missing_cols = [
    _col for _col in chronos_required_cols
    if _col not in chronos_mode_age_source.columns
]
if chronos_missing_cols:
    raise RuntimeError(
        f"Could not find {{chronos_missing_cols!r}} in {str(chronos_mode_ages_path)!r}."
    )
if chronos_mode_age_source['name'].duplicated().any():
    raise RuntimeError('Chronos mode age source has duplicate cluster names.')
chronos_mode_age_source = chronos_mode_age_source[chronos_required_cols].rename(
    columns={{
        'age_chronos_lo': 'chronos_age_lo_myr',
        'age_chronos_mode': 'chronos_age_myr',
        'age_chronos_hi': 'chronos_age_hi_myr',
    }}
)
for _chronos_col in [
    'chronos_age_lo_myr',
    'chronos_age_myr',
    'chronos_age_hi_myr',
]:
    chronos_mode_age_source[_chronos_col] = pd.to_numeric(
        chronos_mode_age_source[_chronos_col],
        errors='coerce',
    )
chronos_mode_age_source = chronos_mode_age_source.loc[
    chronos_mode_age_source['chronos_age_myr'].notnull()
].copy()
df_hunt_chronos_full = df_hunt_full.drop(
    columns=['chronos_age_lo_myr', 'chronos_age_myr', 'chronos_age_hi_myr'],
    errors='ignore',
).merge(
    chronos_mode_age_source,
    on='name',
    how='left',
    validate='1:1',
)
if 'hunt_age_myr' not in df_hunt_chronos_full.columns:
    df_hunt_chronos_full['hunt_age_myr'] = pd.to_numeric(
        df_hunt_chronos_full['age_myr'] if 'age_myr' in df_hunt_chronos_full.columns else pd.Series(np.nan, index=df_hunt_chronos_full.index),
        errors='coerce',
    )
if 'hunt_initial_mass_msun' not in df_hunt_chronos_full.columns:
    df_hunt_chronos_full['hunt_initial_mass_msun'] = pd.to_numeric(
        df_hunt_chronos_full['mass_all_previous'] if 'mass_all_previous' in df_hunt_chronos_full.columns else (
            df_hunt_chronos_full['mass_all'] if 'mass_all' in df_hunt_chronos_full.columns else pd.Series(np.nan, index=df_hunt_chronos_full.index)
        ),
        errors='coerce',
    )
if 'mass_all' in df_hunt_chronos_full.columns:
    df_hunt_chronos_full['chronos_initial_mass_msun'] = pd.to_numeric(
        df_hunt_chronos_full['mass_all'],
        errors='coerce',
    )
else:
    df_hunt_chronos_full['chronos_initial_mass_msun'] = np.nan
df_hunt_chronos_full['display_age_myr'] = df_hunt_chronos_full['chronos_age_myr'].combine_first(
    df_hunt_chronos_full['hunt_age_myr']
)
df_hunt_chronos_full['display_initial_mass_msun'] = df_hunt_chronos_full['chronos_initial_mass_msun'].combine_first(
    df_hunt_chronos_full['hunt_initial_mass_msun']
)
df_hunt_chronos_full['mass_all'] = df_hunt_chronos_full['display_initial_mass_msun']
main_figure_cluster_xy_min_pc = {MAIN_FIGURE_CLUSTER_XY_BOUNDS_PC[0]!r}
main_figure_cluster_xy_max_pc = {MAIN_FIGURE_CLUSTER_XY_BOUNDS_PC[1]!r}
main_figure_cluster_young_age_myr = {MAIN_FIGURE_CLUSTER_YOUNG_AGE_MYR!r}
main_figure_cluster_blue_max_age_myr = {MAIN_FIGURE_CLUSTER_BLUE_MAX_AGE_MYR!r}
main_figure_cluster_grey_max_age_myr = {MAIN_FIGURE_CLUSTER_GREY_MAX_AGE_MYR!r}
df_hunt_chronos_sample = df_hunt_chronos_full.loc[
    (df_hunt_chronos_full['U_err'] < 10) &
    (df_hunt_chronos_full['V_err'] < 10) &
    (df_hunt_chronos_full['W_err'] < 10) &
    (df_hunt_chronos_full['U'].notnull()) &
    (df_hunt_chronos_full['V'].notnull()) &
    (df_hunt_chronos_full['W'].notnull()) &
    (df_hunt_chronos_full['x'].notnull()) &
    (df_hunt_chronos_full['y'].notnull()) &
    (df_hunt_chronos_full['z'].notnull()) &
    (df_hunt_chronos_full['x'].between(main_figure_cluster_xy_min_pc, main_figure_cluster_xy_max_pc)) &
    (df_hunt_chronos_full['y'].between(main_figure_cluster_xy_min_pc, main_figure_cluster_xy_max_pc)) &
    (df_hunt_chronos_full['display_age_myr'].notnull()) &
    (df_hunt_chronos_full['n_rvs_2026'] >= 3)
].copy()
df_hunt_chronos_sample['age_myr'] = df_hunt_chronos_sample['display_age_myr']
df_hunt_chronos_sample['age_source'] = np.where(
    df_hunt_chronos_sample['chronos_age_myr'].notnull(),
    {model_display!r},
    'Hunt',
)
df_hunt_60_chronos = df_hunt_chronos_sample.loc[
    df_hunt_chronos_sample['age_myr'] < main_figure_cluster_blue_max_age_myr
].copy()
df_hunt_0_to_150_chronos = df_hunt_chronos_sample.loc[
    df_hunt_chronos_sample['age_myr'] < main_figure_cluster_grey_max_age_myr
].copy()
df_hunt_young_chronos = df_hunt_chronos_sample.loc[
    df_hunt_chronos_sample['age_myr'] < main_figure_cluster_young_age_myr
].copy()
print(
    'Using Chronos mode ages from {str(chronos_mode_ages_path)} for displayed cluster samples: '
    f"{{len(df_hunt_60_chronos)}} clusters <{{main_figure_cluster_blue_max_age_myr:g}} Myr, "
    f"{{len(df_hunt_0_to_150_chronos)}} clusters 0-{{main_figure_cluster_grey_max_age_myr:g}} Myr, "
    f"{{len(df_hunt_young_chronos)}} clusters <{{main_figure_cluster_young_age_myr:g}} Myr. "
    f"Chronos ages use Hunt fallbacks when missing. "
    f"x/y within [{{main_figure_cluster_xy_min_pc:g}}, {{main_figure_cluster_xy_max_pc:g}}] pc; no z cut applied."
)
""".strip()


def patch_script_source(
    source: str,
    output_html: Path,
    theme_key: str | None,
    *,
    minimal_mode: bool = False,
    mobile_mode: bool = False,
    galactic_simple: bool = False,
    mist_ages: bool = False,
    compact_payload: bool = True,
    mobile_safe_mode: bool = False,
    chronos_results_path: Path = CHRONOS_CLUSTER_RESULTS_PATH,
    chronos_model: str = DEFAULT_CHRONOS_CLUSTER_MODEL,
    include_spiral_arms: bool = False,
    jun6_catalog: bool = False,
) -> str:
    source = source.replace("/Users/cam", str(HOME_DIR))
    source = source.replace("/Users/cam/Desktop", str(DESKTOP_ROOT))

    if jun6_catalog:
        jun6_catalog_block = build_jun6_catalog_source_block()
        source, replaced_jun6_catalog = re.subn(
            (
                r"(?ms)^df_hunt_full\s*=\s*pd\.read_csv\("
                r"'/Users/(?:cam|swiggumc)/Desktop/astro_research/supernovae_map_work/clusters/vels_output/2026-02-04/cluster_velocities_jan2026\.csv'\)"
                r".*?^df_hunt_full\s*=\s*df_hunt_full\.rename\(columns=\{.*?^\}\)\s*$"
            ),
            jun6_catalog_block,
            source,
            count=1,
        )
        if replaced_jun6_catalog != 1:
            raise RuntimeError("Could not replace the notebook cluster catalog block with the Jun 6 catalog block.")

    source = patch_base_cluster_filter(source)

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

    if "df_hunt_60_chronos" not in source:
        if jun6_catalog and Path(chronos_results_path).resolve() == CHRONOS_CLUSTER_RESULTS_PATH.resolve():
            chronos_cluster_block = build_chronos_mode_age_sample_source_block(
                JUN6_CHRONOS_MODE_AGES_PATH,
                model_display=chronos_model.upper(),
            )
        else:
            chronos_cluster_block = build_chronos_cluster_sample_source_block(
                chronos_results_path=chronos_results_path,
                chronos_model=chronos_model,
            )
        source, inserted_chronos_clusters = re.subn(
            r"(?m)^(df_hunt_old\s*=\s*df_hunt_good\.loc\[df_hunt_good\['age_myr'\]\.between\(30,\s*60\)\]\s*)$",
            r"\1\n\n" + chronos_cluster_block,
            source,
            count=1,
        )
        if inserted_chronos_clusters != 1:
            raise RuntimeError("Could not inject Chronos cluster samples into converted notebook script.")

        source, replaced_young_chronos_sample = re.subn(
            r"young_trace\s*=\s*Trace\(df_hunt_young,",
            "young_trace = Trace(df_hunt_young_chronos,",
            source,
            count=1,
        )
        if replaced_young_chronos_sample != 1:
            raise RuntimeError("Could not point the <15 Myr cluster trace at the Chronos sample.")

        full_chronos_source = (
            "df_hunt_60_chronos.sample("
            f"n=min(len(df_hunt_60_chronos), {MOBILE_SAFE_BLUE_CLUSTER_MAX_POINTS}), "
            "random_state=29)"
            if mobile_safe_mode
            else "df_hunt_60_chronos"
        )
        source, replaced_full_chronos_sample = re.subn(
            r"full_sample_trace\s*=\s*Trace\(df_hunt_60,",
            f"full_sample_trace = Trace({full_chronos_source},",
            source,
            count=1,
        )
        if replaced_full_chronos_sample != 1:
            raise RuntimeError("Could not point the <60 Myr cluster trace at the Chronos sample.")

        old_cluster_source = (
            "df_hunt_0_to_150_chronos.sample("
            f"n=min(len(df_hunt_0_to_150_chronos), {MOBILE_SAFE_BACKGROUND_CLUSTER_MAX_POINTS}), "
            "random_state=23)"
            if mobile_safe_mode
            else "df_hunt_0_to_150_chronos"
        )
        old_cluster_trace_source = (
            "old_sample_trace = Trace("
            f"{old_cluster_source}, "
            "data_name = 'Clusters (0-150 Myr)', "
            "min_size = 0, max_size = 7, "
            f"color = {MAIN_FIGURE_CLUSTER_GREY_COLOR!r}, "
            "opacity = .55, marker_style = 'circle', show_tracks = False, size_by_n_stars=False)"
        )
        source, inserted_old_cluster_trace = re.subn(
            r"(?m)^(full_sample_trace\s*=\s*Trace\(.*)$",
            r"\1\n" + old_cluster_trace_source,
            source,
            count=1,
        )
        if inserted_old_cluster_trace != 1:
            raise RuntimeError("Could not add the 0-150 Myr cluster trace.")

        source, inserted_old_cluster_collection = re.subn(
            r"(?m)^(\s*)full_sample_trace,\s*$",
            r"\1old_sample_trace,\n\1full_sample_trace,",
            source,
            count=1,
        )
        if inserted_old_cluster_collection != 1:
            raise RuntimeError("Could not add the 0-150 Myr trace to the TraceCollection.")

        source, inserted_old_cluster_grouping = re.subn(
            r'"Clusters"\s*:\s*\[\s*\'Sun\'\s*,\s*\'Clusters\s*\(<\s*60\s*Myr\)\'\s*,\s*\'Clusters\s*\(<\s*15\s*Myr\)\'\s*\]',
            "\"Clusters\": ['Sun', 'Clusters (0-150 Myr)', 'Clusters (< 60 Myr)', 'Clusters (< 15 Myr)']",
            source,
            count=1,
        )
        if inserted_old_cluster_grouping != 1:
            raise RuntimeError("Could not add the 0-150 Myr trace to the cluster grouping.")

        source, recolored_full_chronos_sample = re.subn(
            r"(full_sample_trace\s*=\s*Trace\([^\n]*?data_name\s*=\s*'Clusters\s*\(<\s*60\s*Myr\)'[^\n]*?color\s*=\s*)(['\"][^'\"]+['\"])",
            r"\g<1>'#2f80ff'",
            source,
            count=1,
        )
        if recolored_full_chronos_sample != 1:
            raise RuntimeError("Could not recolor the <60 Myr Chronos cluster trace.")

    time_step_myr = MOBILE_SAFE_TIMESTEP_MYR if mobile_safe_mode else MAIN_FIGURE_TIMESTEP_MYR
    source, replaced_time_grid = re.subn(
        r"(?m)^time_int\s*=\s*np\.round\(np\.arange\(0,\s*-66,\s*-1\),\s*1\)\s*$",
        f"time_int = np.round(np.arange(0, -{MAIN_FIGURE_LOOKBACK_MYR + 1}, -{time_step_myr}), 1)",
        source,
        count=1,
    )
    if replaced_time_grid != 1:
        raise RuntimeError("Could not rewrite the main figure time grid.")

    if galactic_simple:
        source, replaced_time_grid = re.subn(
            r"(?m)^time_int\s*=\s*np\.round\(np\.arange\(0,\s*-\d+,\s*-\d+\),\s*1\)\s*$",
            f"time_int = np.round(np.arange(0, -{MAIN_FIGURE_LOOKBACK_MYR + 1}, -{time_step_myr}), 1)",
            source,
            count=1,
        )
        if replaced_time_grid != 1:
            raise RuntimeError("Could not restore the galactic simple time grid to 1 Myr steps.")

    source = patch_edenhofer_volume_integer_setting(
        source,
        "max_resolution",
        MOBILE_SAFE_DUST_MAX_RESOLUTION if mobile_safe_mode else MAIN_FIGURE_DUST_MAX_RESOLUTION,
    )
    source = patch_edenhofer_volume_integer_setting(
        source,
        "max_resolution_cap",
        MOBILE_SAFE_DUST_MAX_RESOLUTION_CAP if mobile_safe_mode else MAIN_FIGURE_DUST_MAX_RESOLUTION_CAP,
        insert_after="max_resolution",
    )
    source = patch_edenhofer_volume_integer_setting(
        source,
        "samples",
        MOBILE_SAFE_DUST_SAMPLES if mobile_safe_mode else MAIN_FIGURE_DUST_SAMPLES,
    )
    source = patch_edenhofer_volume_numeric_setting(
        source,
        "opacity",
        MAIN_FIGURE_DUST_OPACITY,
    )
    source = patch_edenhofer_volume_numeric_setting(
        source,
        "alpha_coef",
        MAIN_FIGURE_DUST_ALPHA_COEF,
    )
    source = patch_edenhofer_volume_clip_bounds(
        source,
        MAIN_FIGURE_VOLUME_Z_CLIP_BOUNDS,
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
        galaxy_image_path = (
            mobile_safe_galaxy_image_path()
            if mobile_safe_mode
            else GALACTIC_PLANE_IMAGE_PATH
        )
        if mobile_safe_mode:
            volume_initial_state_bits = [
                "'active_volume_key': 'edenhofer-dust-desktop'",
                "'mobile_active_volume_key': 'edenhofer-dust-mobile'",
                "'mobile_defer_volumes': False",
                "'legend_state': {'edenhofer-dust-desktop': True, 'edenhofer-dust-mobile': False}",
                "'volume_state_by_key': {'edenhofer-dust-desktop': {'visible': True}, 'edenhofer-dust-mobile': {'visible': False}}",
            ]
        else:
            volume_initial_state_bits = [
                "'active_volume_key': ('supernova-density' if supernova_volumes else 'volume-0')",
                "'legend_state': ({'volume-0': False, 'supernova-density': True} if supernova_volumes else {'volume-0': True})",
                "'volume_state_by_key': ({'volume-0': {'visible': False}, 'supernova-density': {'visible': True}} if supernova_volumes else {'volume-0': {'visible': True}})",
            ]
        initial_state_bits = [
            "'current_group': 'Clusters'",
            "'click_selection_enabled': False",
            f"'compact_payload_enabled': {bool(compact_payload)!r}",
            f"'compact_widget_payload_enabled': {bool(mobile_safe_mode)!r}",
            f"'scene_float_precision': {1 if mobile_safe_mode else 2}",
            *volume_initial_state_bits,
            "'galaxy_image': True",
            f"'galaxy_image_path': {str(galaxy_image_path)!r}",
            "'galaxy_image_size_pc': 40000.0",
            "'galaxy_image_opacity': 0.35",
            "'galaxy_image_hide_below_scale_bar_pc': 420.0",
            "'galaxy_image_fade_start_scale_bar_pc': 700.0",
            "'sky_dome_enabled': True",
            "'sky_dome_background_mode': 'live_aladin'",
            "'sky_dome_source': 'aladin'",
            "'sky_dome_projection': 'TAN'",
            f"'sky_dome_capture_width_px': {MOBILE_SAFE_SKY_DOME_CAPTURE_WIDTH_PX if mobile_safe_mode else 4096}",
            f"'sky_dome_capture_height_px': {MOBILE_SAFE_SKY_DOME_CAPTURE_HEIGHT_PX if mobile_safe_mode else 2048}",
            "'sky_dome_capture_format': 'image/jpeg'",
            f"'sky_dome_capture_quality': {MOBILE_SAFE_SKY_DOME_CAPTURE_QUALITY if mobile_safe_mode else 0.94}",
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
        if mobile_mode:
            initial_state_bits.append("'mobile_mode_enabled': True")
            initial_state_bits.append("'legend_open': False")
        if galactic_simple:
            initial_state_bits.extend([
                "'galactic_simple_mode_enabled': True",
                f"'galactic_plane_image_path': {str(galaxy_image_path)!r}",
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
        full_catalog_source = (
            "df_hunt_good.sample("
            f"n=min(len(df_hunt_good), {MOBILE_SAFE_FULL_CATALOG_MAX_POINTS}), "
            "random_state=17)"
            if mobile_safe_mode
            else "df_hunt_good"
        )
        catalog_trace_block = f"""
full_catalog_trace = Trace(
    {full_catalog_source},
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

    if include_spiral_arms:
        source, replaced_spiral_setting = re.subn(
            r"(?m)^(\s*)include_spiral_arms\s*=\s*False\s*,\s*$",
            r"\1include_spiral_arms=True,",
            source,
            count=1,
        )
        if replaced_spiral_setting != 1:
            raise RuntimeError("Could not enable spiral arms in the main figure make_plot call.")

        spiral_grouping_block = f"""
spiral_arm_trace_names = {list(MAIN_FIGURE_SPIRAL_ARM_TRACE_NAMES)!r}
for _spiral_trace_name in spiral_arm_trace_names:
    if _spiral_trace_name not in trace_groupings.get('All', []):
        trace_groupings.setdefault('All', []).append(_spiral_trace_name)
    if _spiral_trace_name not in trace_groupings.get('Clusters', []):
        trace_groupings.setdefault('Clusters', ['Sun']).append(_spiral_trace_name)
trace_groupings['Spiral Arms'] = ['Sun', *spiral_arm_trace_names]
""".strip()
        source, inserted_spiral_group = re.subn(
            r"(?m)^(plot_3d\s*=\s*Animate3D\()",
            spiral_grouping_block + "\n\n" + r"\1",
            source,
            count=1,
        )
        if inserted_spiral_group != 1:
            raise RuntimeError("Could not add spiral arms to the trace groupings.")

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
            r'"Clusters"\s*:\s*\[\s*\'Sun\'\s*,\s*(?:\'Clusters\s*\(0-150\s*Myr\)\'\s*,\s*)?\'Clusters\s*\(<\s*60\s*Myr\)\'\s*,\s*\'Clusters\s*\(<\s*15\s*Myr\)\'\s*\]',
            "\"Clusters\": ['Sun', 'Clusters (0-150 Myr)', 'Clusters (< 60 Myr)']",
            source,
            count=1,
        )
        if removed_young_from_grouping != 1:
            raise RuntimeError("Could not remove the <15 Myr cluster trace from galactic simple groupings.")

    if mobile_safe_mode:
        if "supernova_volumes = []" not in source:
            source, inserted_supernova_stub = re.subn(
                r"(?m)^(fig3d\s*=\s*plot_3d\.make_plot\()",
                "supernova_volumes = []\n\n" + r"\1",
                source,
                count=1,
            )
            if inserted_supernova_stub != 1:
                raise RuntimeError("Could not inject mobile-safe supernova volume stub.")
    elif "supernova_volumes = _build_supernova_volumes_for_main_figure(time_int)" not in source:
        supernova_block = build_supernova_volume_source_block()
        source, inserted_supernova = re.subn(
            r"(?m)^(fig3d\s*=\s*plot_3d\.make_plot\()",
            supernova_block + "\n\n" + r"\1",
            source,
            count=1,
        )
        if inserted_supernova != 1:
            raise RuntimeError("Could not inject supernova volume helper block.")

    if mobile_safe_mode and "edenhofer_volumes = [edenhofer_volume_desktop, edenhofer_volume_mobile]" not in source:
        edenhofer_variant_block = build_adaptive_edenhofer_volume_source_block()
        source, inserted_edenhofer_variants = re.subn(
            r"(?m)^(fig3d\s*=\s*plot_3d\.make_plot\()",
            edenhofer_variant_block + "\n\n" + r"\1",
            source,
            count=1,
        )
        if inserted_edenhofer_variants != 1:
            raise RuntimeError("Could not inject adaptive Edenhofer volume helper block.")

    volume_list_source = (
        "volumes=[*edenhofer_volumes, *supernova_volumes]"
        if mobile_safe_mode
        else "volumes=[edenhofer_volume, *supernova_volumes]"
    )
    source, replaced_supernova_volumes = re.subn(
        r"volumes\s*=\s*\[\s*edenhofer_volume\s*,\s*mccallum_ne\s*\]",
        volume_list_source,
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
    mobile_mode: bool = False,
    galactic_simple: bool = False,
    mist_ages: bool = False,
    compact_payload: bool = True,
    mobile_safe_mode: bool = False,
    chronos_results_path: Path = CHRONOS_CLUSTER_RESULTS_PATH,
    chronos_model: str = DEFAULT_CHRONOS_CLUSTER_MODEL,
    include_spiral_arms: bool = False,
    jun6_catalog: bool = False,
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
        mobile_mode=mobile_mode,
        galactic_simple=galactic_simple,
        mist_ages=mist_ages,
        compact_payload=compact_payload,
        mobile_safe_mode=mobile_safe_mode,
        chronos_results_path=chronos_results_path,
        chronos_model=chronos_model,
        include_spiral_arms=include_spiral_arms,
        jun6_catalog=jun6_catalog,
    )

    run_script_source(patched_source)
    write_ar_quicklook_service_worker(output_html.parent)
    if website_output_html is not None:
        website_output_html = Path(website_output_html)
        if output_html.resolve() != website_output_html.resolve():
            website_output_html.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(output_html, website_output_html)
        write_ar_quicklook_service_worker(website_output_html.parent)
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
        "--mobile-mode",
        action="store_true",
        help="Render the figure with mobile-mode metadata enabled.",
    )
    parser.add_argument(
        "--mobile-safe-mode",
        action="store_true",
        help=(
            "Render a smaller upload/iPhone-safe payload with coarser time steps, "
            "lower-resolution volumes, lower sky capture resolution, and no "
            "embedded supernova density volume."
        ),
    )
    parser.add_argument(
        "--mist-ages",
        action="store_true",
        help="Use the matched MIST isochrone ages instead of the default notebook ages.",
    )
    parser.add_argument(
        "--chronos-results-path",
        type=Path,
        default=CHRONOS_CLUSTER_RESULTS_PATH,
        help="Chronos cluster_results.csv path used for displayed cluster ages.",
    )
    parser.add_argument(
        "--chronos-model",
        choices=("parsec", "mist"),
        default=DEFAULT_CHRONOS_CLUSTER_MODEL,
        help="Model prefix to use from the Chronos cluster_results.csv.",
    )
    parser.add_argument(
        "--include-spiral-arms",
        action="store_true",
        help="Include the rotating spiral arm model traces in the galactic-mode figure.",
    )
    parser.add_argument(
        "--jun6-catalog",
        action="store_true",
        help="Use the Jun 6 Chronos mode ages and bulkfit no-MWM cluster velocity catalog.",
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
        mobile_mode=bool(args.mobile_mode),
        mist_ages=bool(args.mist_ages),
        compact_payload=not bool(args.full_payload),
        mobile_safe_mode=bool(args.mobile_safe_mode),
        chronos_results_path=args.chronos_results_path,
        chronos_model=args.chronos_model,
        include_spiral_arms=bool(args.include_spiral_arms),
        jun6_catalog=bool(args.jun6_catalog),
        website_output_html=None if args.no_website_copy else WEBSITE_OUTPUT_HTML,
    )
    print(f"Wrote {output_html}")
    if not args.no_website_copy:
        print(f"Copied {output_html} to {WEBSITE_OUTPUT_HTML}")


if __name__ == "__main__":
    main()
