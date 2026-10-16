#!/usr/bin/env python3
"""The background figure of Cameren Swiggum's website (cswigg.github.io/cam_website).

The October 1 science (``tests/main_figure_oct1.py``: the July 25 clusters in
the Khalil et al. 2025 spiral potential, read from its stored classic source)
trimmed to what a page background needs: the Sun, the clusters younger than
60 and 15 Myr, the Galactic guides and the Edenhofer et al. (2024) dust. Each
first-author paper whose results live in this space adds layers that start
hidden, and a saved view ("State") that shows them:

``swiggum-2021``  Swiggum et al. 2021, ApJ 917, 21. The Briceño-1 and
                  OBP-near stars (present day) and each star's motion relative
                  to the pair over 3 Myr, which points away from a common centre.
``swiggum-2022``  Swiggum et al. 2022, A&A 664, L13. The Radcliffe Wave model
                  (Alves et al. 2020), the Local Arm masers and arm fit (Reid
                  et al. 2019), with the clusters younger than 15 Myr.
``swiggum-2024``  Swiggum et al. 2024, Nature 631, 49. The Alpha Persei,
                  Collinder 135 and Messier 6 families (already in the October
                  1 science), 30 Myr ago.
``swiggum-2025``  Swiggum et al. 2025, A&A 699, L5. The clusters of the
                  Pleiades, Coma Berenices and Sirius branches (the paper's CDS
                  table), integrated in the same potential, with the families.
``swiggum-2026``  Swiggum et al. 2026, arXiv:2608.20307. The reconstructed
                  core-collapse supernova density (10 Myr trailing window): the
                  events and recipe of the published supernova map.

The figure is written for an iframe under the page's text: the Oviz viewer's
embed mode (no title or search; the key lists the layers in view), opening
in Focus mode, anchored to the LSR. The page passes ``?padding=`` and calls
``setPadding`` / ``goTo`` through the viewer's message bridge.

    python tests/website_figure.py --output-html OUT.html
"""

from __future__ import annotations

import argparse
import copy
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import main_figure_oct1 as oct1  # noqa: E402
from figure_provenance import hash_provenance_paths, short_sha256  # noqa: E402
from oviz import Animate3D, Trace, TraceCollection  # noqa: E402
from oviz.spiral_models import khalil2025_potential  # noqa: E402
from oviz.threejs_embed import scene_spec_jsonable  # noqa: E402
from oviz.viewer.compile import compile_scene_spec  # noqa: E402
from oviz.viewer.figure import OvizFigure  # noqa: E402
from oviz.viewer.upgrade import read_legacy_scene_spec  # noqa: E402

RESEARCH = Path("~/Desktop/astro_research").expanduser()
ORION_STARS = RESEARCH / "orion/orion_with_xyz.csv"
RADCLIFFE_WAVE_MODEL = RESEARCH / "orion/orion_populations/Best_Fit_Wave_Model.csv"
REID_MASERS = RESEARCH / "radcliffe/reid_masers.csv"
REID_LOCAL_ARM = RESEARCH / "radcliffe/reid_local_fit.csv"
S25_CLUSTERS = RESEARCH / "supernovae_map/data/clusters/swiggum2025_clusters.dat"
SUPERNOVA_FIGURE = RESEARCH / "supernovae_map/oviz_supernovae/oviz_supernovae_july29_ice_kde.html"
OUTPUT_HTML = HERE / "website_figure.html"

# What the October 1 science keeps (by legend name, so its trace keys may change).
KEEP_TRACE_NAMES = (
    "Sun", "Clusters (< 60 Myr)", "Clusters (< 15 Myr)",
    "Alpha Persei Family", "Cr 135 Family", "M6 Family",
)
KEEP_VOLUMES = ("volume-0",)  # Edenhofer+2024 dust
GROUP = "Clusters"
TIME_MYR = np.round(np.arange(0, -61, -1), 1)  # as tests/main_figure_oct1.py's science
SUN_ROW = dict(name="Sun", age_myr=8000.0, U=0.0, V=0.0, W=0.0, x=0.0, y=0.0, z=27.0, n_stars=1)
SOLAR_MOTION_KMS = np.array([11.1, 12.24, 7.25])  # Schoenrich et al. (2010): the Orion table adds it
KMS_TO_PC_PER_MYR = 1.0227121650537077

# Swiggum et al. (2021): Table 1 labels of the two groups that expand.
ORION_GROUPS = ((1, "Briceño-1", "#ffb454"), (3, "OBP-near", "#7cc8ff"))
ORION_ARROW_MYR = 4.0
ORION_RV_ERROR_MAX_KMS = 0.3  # precise (mostly APOGEE-2) RVs: the line-of-sight spread then matches the tangential
# Swiggum et al. (2025): the three branches (CDS table group_name).
BRANCHES = (
    ("Pleiades", "Pleiades branch", "#9aa8ff"),
    ("Coma", "Coma Berenices branch", "#ffcf8a"),
    ("Sirius", "Sirius branch", "#9be3a6"),
)
S25_COLUMNS = (  # CDS byte ranges (1-based, inclusive) → 0-based half-open
    ("name", 1, 17), ("age_myr", 19, 26), ("x", 28, 35), ("y", 37, 44), ("z", 46, 53),
    ("U", 55, 61), ("V", 63, 69), ("W", 71, 77), ("group_name", 123, 130),
)

KEYS = {
    "orion_stars": "paper-2021-stars",
    "orion_motion": "paper-2021-motion",
    "rw_model": "paper-2022-radcliffe-wave",
    "masers": "paper-2022-masers",
    "arm_fit": "paper-2022-local-arm",
    "supernovae": "supernova-density",
}


# ---------------------------------------------------------------- the base

def trace_keys_by_name(spec: dict) -> dict[str, str]:
    return {item["name"]: item["key"] for item in spec["legend"]["items"] if item.get("kind", "trace") == "trace"}


def simplify(spec: dict) -> dict:
    """Keep the October 1 science a page background needs; drop the rest (and its bytes)."""
    names = trace_keys_by_name(spec)
    missing = [n for n in KEEP_TRACE_NAMES if n not in names]
    if missing:
        raise RuntimeError(f"The October 1 source lacks {missing}")
    keep_legend = {names[n] for n in KEEP_TRACE_NAMES} | set(KEEP_VOLUMES)
    # Guides (rings, labels, quadrants) have no legend entry: keep every trace not in the legend.
    legend_keys = {item["key"] for item in spec["legend"]["items"]}
    keep = lambda key: key in keep_legend or key not in legend_keys  # noqa: E731
    guides = {t["key"] for t in spec["frames"][0]["traces"] if t["key"] not in legend_keys}
    arms = {t["key"] for t in spec["frames"][0]["traces"] if str(t.get("name", "")).startswith("Spiral arms")}
    keep_trace = lambda key: keep(key) and key not in arms  # noqa: E731
    for frame in spec["frames"]:
        frame["traces"] = [t for t in frame["traces"] if keep_trace(t["key"])]
        frame["decorations"] = [d for d in frame.get("decorations") or []
                                if d.get("kind") != "volume_layer" or d.get("key") in KEEP_VOLUMES]
    spec["legend"]["items"] = [i for i in spec["legend"]["items"] if i["key"] in keep_legend]
    spec["volumes"]["layers"] = [v for v in spec["volumes"]["layers"] if v["key"] in KEEP_VOLUMES]
    base = spec["group_visibility"][GROUP]
    # One group holds every layer; False would leave a layer out of it (and its key).
    spec["group_visibility"] = {GROUP: {k: ("legendonly" if v is False else v) for k, v in base.items() if keep_trace(k)}}
    spec["group_order"] = [GROUP]
    spec["default_group"] = GROUP
    init = spec["initial_state"]
    init["current_group"] = GROUP
    init["legend_state"] = {k: v for k, v in init.get("legend_state", {}).items() if k in KEEP_VOLUMES}
    init["volume_state_by_key"] = {k: v for k, v in init.get("volume_state_by_key", {}).items() if k in KEEP_VOLUMES}
    # A background keeps to the figure: no analysis widgets.
    for widget in ("age_kde", "cluster_filter", "dendrogram"):
        if isinstance(spec.get(widget), dict):
            spec[widget] = {"enabled": False}
    spec["actions"] = {"enabled": False, "items": []}
    members = (spec.get("sky_panel") or {}).get("members_by_cluster")
    if isinstance(members, dict):
        kept_names = set(KEEP_TRACE_NAMES)
        spec["sky_panel"]["members_by_cluster"] = {k: v for k, v in members.items() if k in kept_names}
    spec["title"] = "Star clusters near the Sun"
    return {"guides": sorted(guides), "dropped_arms": sorted(arms)}


def present_index(spec: dict) -> int:
    t0 = int(spec["initial_frame_index"])
    if float(spec["frames"][t0]["time"]) != 0.0:
        raise RuntimeError("The initial frame is not t = 0.")
    return t0


def add_layer(spec: dict, entry: dict, legend: dict, *, frames: list[dict] | None = None) -> None:
    """Add a trace (one dict per frame, or ``entry`` in the present-day frame only), hidden."""
    key = entry["key"] if frames is None else frames[0]["key"]
    if frames is None:
        spec["frames"][present_index(spec)]["traces"].append(entry)
    else:
        for frame, e in zip(spec["frames"], frames):
            frame["traces"].append(e)
    spec["legend"]["items"].insert(
        max(i for i, item in enumerate(spec["legend"]["items"]) if item.get("kind", "trace") == "trace") + 1,
        {**legend, "key": key})
    for visibility in spec["group_visibility"].values():
        visibility[key] = "legendonly"


def points_entry(key: str, name: str, color: str, size: float, points: list[dict], *, opacity: float = 0.95,
                 symbol: str = "circle") -> tuple[dict, dict]:
    entry = {"key": key, "name": name, "showlegend": True, "size_by_n_stars_default": False, "points": points,
             "default_point_size": size, "has_n_stars": False, "default_opacity": opacity, "legend_color": color,
             "default_color": color, "default_symbol": symbol}
    legend = {"name": name, "color": color, "kind": "trace", "has_points": True, "has_segments": False,
              "has_labels": False, "has_n_stars": False, "size_by_n_stars_default": False, "default_color": color,
              "default_opacity": opacity, "default_point_size": size}
    return entry, legend


def segments_entry(key: str, name: str, color: str, segments: np.ndarray, *, width: float = 1.5,
                   opacity: float = 0.9) -> tuple[dict, dict]:
    entry = {"key": key, "name": name, "showlegend": True, "size_by_n_stars_default": False,
             "segments": np.round(segments, 2).tolist(), "line": {"color": color, "width": width, "dash": "solid"},
             "opacity": opacity, "legend_color": color, "default_opacity": opacity}
    legend = {"name": name, "color": color, "kind": "trace", "has_points": False, "has_segments": True,
              "has_labels": False, "has_n_stars": False, "size_by_n_stars_default": False, "default_color": color,
              "default_opacity": opacity}
    return entry, legend


def star_point(x, y, z, label: str, trace_name: str, color: str, extra: str = "") -> dict:
    """A present-day object: Galactic (l, b, d) from its heliocentric position (the scene's frame at t = 0)."""
    d = float(np.sqrt(x * x + y * y + z * z))
    l = float(np.degrees(np.arctan2(y, x)) % 360.0)
    b = float(np.degrees(np.arcsin(z / d))) if d > 0 else 0.0
    x, y, z = (round(float(c), 2) for c in (x, y, z))
    return {
        "x": x, "y": y, "z": z,
        "hovertext": f'<b style="font-size:16px;">{label}</b><br>{trace_name}{extra}<br>d = {d:.0f} pc',
        "selection": {"l_deg": l, "b_deg": b, "dist_pc": d, "x0": x, "y0": y, "z0": z, "click_time_myr": 0.0,
                      "cluster_name": label, "cluster_color": color, "trace_name": trace_name, "name_all": ""},
    }


# ---------------------------------------------------------------- 2021: Orion

def add_orion_layer(spec: dict, path: Path = ORION_STARS) -> dict:
    stars = pd.read_csv(path)
    labels = [lab for lab, _, _ in ORION_GROUPS]
    pair = stars[stars["label"].isin(labels)].copy()
    # The table's velocities include the solar motion; the paper's frame is the Sun's.
    vel = pair[["vx", "vy", "vz"]].to_numpy() - SOLAR_MOTION_KMS
    good = np.isfinite(vel).all(axis=1) & (pair["rv_error"].to_numpy() <= ORION_RV_ERROR_MAX_KMS)
    mean_v = np.median(vel[good], axis=0)
    pos = pair[["x", "y", "z"]].to_numpy()
    counts = {}
    for lab, name, color in ORION_GROUPS:
        sel = (pair["label"] == lab).to_numpy()
        points = [star_point(*pos[i], f"Gaia DR3 {pair['source_id'].iat[i]}", f"{name} (Swiggum+ 2021)", color)
                  for i in np.flatnonzero(sel)]
        entry, legend = points_entry(f"{KEYS['orion_stars']}-{lab}", f"{name} stars", color, 0.6, points, opacity=0.9)
        add_layer(spec, entry, legend)
        counts[name] = int(sel.sum())
    rel = vel[good] - mean_v
    start = pos[good]
    seg = np.hstack([start, start + rel * ORION_ARROW_MYR * KMS_TO_PC_PER_MYR])
    entry, legend = segments_entry(KEYS["orion_motion"], f"Motion relative to the pair ({ORION_ARROW_MYR:g} Myr)",
                                   "#f4f1e8", seg, width=1.4, opacity=0.9)
    add_layer(spec, entry, legend)
    centre = start.mean(axis=0)
    radial = np.einsum("ij,ij->i", rel, (start - centre) / np.linalg.norm(start - centre, axis=1, keepdims=True))
    return {"table": path.name, "table_sha256": short_sha256(path), "stars": counts,
            "stars_with_3d_velocity": int(good.sum()), "rv_error_max_kms": ORION_RV_ERROR_MAX_KMS,
            "pair_velocity_kms": np.round(mean_v, 3).tolist(), "arrow_myr": ORION_ARROW_MYR,
            "median_radial_velocity_kms": round(float(np.median(radial)), 3),
            "centre_pc": np.round(centre, 1).tolist(),
            "velocities": "table vx, vy, vz minus the Schoenrich et al. (2010) solar motion it includes"}


# ---------------------------------------------------------------- 2022: the Radcliffe Wave

def add_radcliffe_layer(spec: dict) -> dict:
    from astropy import units as u
    from astropy.coordinates import ICRS, Galactic

    wave = pd.read_csv(RADCLIFFE_WAVE_MODEL)[["x", "y", "z"]].to_numpy()
    entry, legend = segments_entry(KEYS["rw_model"], "Radcliffe Wave (Alves+ 2020)", "#ff6b6b",
                                   np.hstack([wave[:-1], wave[1:]]), width=2.4, opacity=0.95)
    add_layer(spec, entry, legend)
    masers = pd.read_csv(REID_MASERS)
    local = masers[masers["Arm"].astype(str).str.strip() == "Loc"]
    c = ICRS(local["ra"].to_numpy() * u.deg, local["dec"].to_numpy() * u.deg,
             (1000.0 / local["plx"].to_numpy()) * u.pc).transform_to(Galactic())
    c.representation_type = "cartesian"
    xyz = np.stack([c.u.to_value(u.pc), c.v.to_value(u.pc), c.w.to_value(u.pc)], axis=-1)
    color = "#ffffff"
    points = [star_point(*xyz[i], str(local["Name"].iat[i]).strip(), "Local Arm maser (Reid+ 2019)", color,
                         f"<br>parallax = {local['plx'].iat[i]:g} mas") for i in range(len(local))]
    entry, legend = points_entry(KEYS["masers"], "Local Arm masers (Reid+ 2019)", color, 14.0, points, symbol="diamond")
    add_layer(spec, entry, legend)
    arm = pd.read_csv(REID_LOCAL_ARM)[["x_kpc", "y_kpc", "z_kpc"]].to_numpy() * 1000.0
    entry, legend = segments_entry(KEYS["arm_fit"], "Local Arm (Reid+ 2019 fit)", "#e9eef5",
                                   np.hstack([arm[:-1], arm[1:]]), width=1.2, opacity=0.6)
    add_layer(spec, entry, legend)
    return {"radcliffe_wave_model": RADCLIFFE_WAVE_MODEL.name, "radcliffe_wave_model_sha256": short_sha256(RADCLIFFE_WAVE_MODEL),
            "masers": REID_MASERS.name, "masers_sha256": short_sha256(REID_MASERS), "local_arm_masers": int(len(local)),
            "local_arm_fit": REID_LOCAL_ARM.name, "local_arm_fit_sha256": short_sha256(REID_LOCAL_ARM),
            "frame": "present-day heliocentric Galactic positions (the scene's frame at t = 0); shown at t = 0 only"}


# ---------------------------------------------------------------- 2025: moving-group branches

def read_s25_table(path: Path = S25_CLUSTERS) -> pd.DataFrame:
    table = pd.read_fwf(path, colspecs=[(a - 1, b) for _, a, b in S25_COLUMNS], names=[n for n, _, _ in S25_COLUMNS])
    table["name"] = table["name"].astype(str).str.strip()
    return table


def integrated_traces(spec: dict, rows: list[tuple[str, str, pd.DataFrame]]) -> tuple[dict, list[dict]]:
    """Integrate ``(key, name, table)`` traces exactly as the October 1 science did; check the frames agree."""
    traces = [Trace(pd.DataFrame([SUN_ROW]), data_name="Sun", min_size=5, max_size=5, color="yellow", opacity=1,
                    marker_style="circle", show_tracks=False, size_by_n_stars=False)]
    for _, name, table in rows:
        traces.append(Trace(table, data_name=name, min_size=0.0, max_size=7, color=table.attrs["color"], opacity=1,
                            marker_style="circle", show_tracks=False, size_by_n_stars=False))
    scene = Animate3D(potential=khalil2025_potential(), data_collection=TraceCollection(traces),
                      xyz_widths=(2000, 2000, 400), figure_theme="dark")
    fig = scene.make_plot(time=TIME_MYR, show=False, save_name=None, fade_in_time=8, fade_in_and_out=False,
                          galactic_mode=True, show_galactic_guides=False, renderer="threejs", viewer="classic")
    made = scene_spec_jsonable(fig.to_dict())
    times = [float(f["time"]) for f in made["frames"]]
    if times != [float(f["time"]) for f in spec["frames"]]:
        raise RuntimeError("The integrated frames do not match the October 1 frames.")
    # The same Sun in both: same frame conventions and potential.
    sun = lambda s: np.array([[p["x"], p["y"], p["z"]] for f in s["frames"] for t in f["traces"]  # noqa: E731
                              if t.get("name") == "Sun" for p in t["points"]])
    gap = float(np.abs(sun(made) - sun(spec)).max())
    if gap > 0.05:
        raise RuntimeError(f"The Sun moves differently from the October 1 science ({gap:.3f} pc).")
    legend = {i["name"]: i for i in made["legend"]["items"]}
    for key, name, _ in rows:
        frames = []
        for frame in made["frames"]:
            entry = copy.deepcopy(next(t for t in frame["traces"] if t.get("name") == name))
            entry["key"] = key
            frames.append(entry)
        add_layer(spec, None, legend[name], frames=frames)
    return {"sun_max_difference_pc": gap}, made["legend"]["items"]


def add_moving_group_layer(spec: dict, path: Path = S25_CLUSTERS) -> dict:
    table = read_s25_table(path)
    rows = []
    for group, name, color in BRANCHES:
        part = table[table["group_name"] == group].reset_index(drop=True)
        part.attrs["color"] = color
        rows.append((f"paper-2025-{group.lower()}", name, part))
    check, _ = integrated_traces(spec, rows)
    return {"table": path.name, "table_sha256": short_sha256(path),
            "branches": {name: len(part) for _, name, part in rows}, **check,
            "orbits": "integrated like the October 1 clusters (MWPotential2014 + Khalil et al. 2025 spirals)"}


# ---------------------------------------------------------------- 2026: supernovae

def add_supernova_layer(spec: dict, path: Path = SUPERNOVA_FIGURE) -> dict:
    """The published map's event KDE: its events and recipe, rebuilt by the viewer per time sample."""
    source = read_legacy_scene_spec(path)
    layer = copy.deepcopy(next(v for v in source["volumes"]["layers"] if v.get("key") == KEYS["supernovae"]))
    src_key = layer["procedural_kde_source_key"] = layer.get("procedural_kde_source_key") or "ccsne-latest"
    kde = copy.deepcopy(source["procedural_kde_sources"][src_key])
    # Hidden at first, so the viewer builds the density when it is shown;
    # the lasso statistics are not offered here.
    for field in ("initial_volume_b64", "lasso_stats_events_b64"):
        kde.pop(field, None)
    spec.setdefault("procedural_kde_sources", {})[src_key] = kde
    layer["name"] = layer["state_name"] = "Supernovae (last 10 Myr)"
    layer["visible"] = False
    spec["volumes"]["layers"].append(layer)
    legend = next(i for i in source["legend"]["items"] if i.get("key") == KEYS["supernovae"])
    spec["legend"]["items"].append({**legend, "name": layer["name"]})
    for visibility in spec["group_visibility"].values():
        visibility[KEYS["supernovae"]] = "legendonly"
    init = spec["initial_state"]
    init["legend_state"][KEYS["supernovae"]] = False
    # The published map's look (asinh shows the diffuse past; see the supernova map figure).
    init["volume_state_by_key"][KEYS["supernovae"]] = {**SUPERNOVA_STYLE, "visible": False}
    return {"figure": path.name, "figure_sha256": short_sha256(path), "events": int(kde.get("event_count", 0)),
            "window_myr": kde.get("rolling_window_myr"), "voxel_pc": kde.get("voxel_size_pc"),
            "imf_realizations": kde.get("imf_realization_count")}


SUPERNOVA_STYLE = {"vmin": 0.0238, "vmax": 5.0, "opacity": 1.0, "stretch": "asinh", "colormap": "ice",
                   "kdeTimeWindowMyr": 10.0}


# ---------------------------------------------------------------- saved views

OPEN_ASPECT = 0.65  # width / height of the part of the window the page leaves uncovered (desktop)


def frame_of(spec: dict, time_myr: float) -> int:
    return int(np.argmin([abs(float(f["time"]) - time_myr) for f in spec["frames"]]))


def positions(spec: dict, keys: list[str], time_myr: float) -> np.ndarray:
    """Points of ``keys`` at ``time_myr`` (present-only layers at t = 0)."""
    out = []
    for fi in (frame_of(spec, time_myr), present_index(spec)):
        for t in spec["frames"][fi]["traces"]:
            if t["key"] in keys:
                out += [[p["x"], p["y"], p["z"]] for p in t.get("points") or []]
                out += [seg[:3] for seg in t.get("segments") or []] + [seg[3:] for seg in t.get("segments") or []]
        if out:
            break
    return np.asarray(out, float)


def fit_view(spec: dict, keys: list[str], time_myr: float, toward: tuple[float, float, float], *,
             fov: float = 45.0, margin: float = 1.25, quantile: float = 0.9, box: bool = False,
             distance: float | None = None) -> dict:
    """A camera on ``keys``' points (their median, or the middle of their extent with ``box``),
    from direction ``toward`` (target → eye), fitting them."""
    pts = positions(spec, keys, time_myr)
    target = np.mean(np.percentile(pts, [2, 98], axis=0), axis=0) if box else np.median(pts, axis=0)
    radius = float(np.quantile(np.linalg.norm(pts - target, axis=1), quantile))
    # The page leaves the figure a part about OPEN_ASPECT wide per unit height.
    half = min(np.radians(fov) / 2, np.arctan(np.tan(np.radians(fov) / 2) * OPEN_ASPECT))
    d = distance or radius / np.tan(half) * margin
    eye = target + np.asarray(toward, float) / np.linalg.norm(toward) * d
    return {"camera": {"position": dict(zip("xyz", map(float, eye))), "target": dict(zip("xyz", map(float, target)))},
            "fov": fov}


# Each paper's view: the layers shown (everything else is hidden), the dust,
# the time, the camera (fitted to the ``fit`` layers from direction ``toward``)
# and any style changes for the view.
PAPER_VIEWS = (
    dict(id="swiggum-2021", name="Swiggum et al. (2021)", show=("orion",), dust=0.55, time=0.0,
         fit="orion", toward=(0.95, 0.35, 0.45), margin=1.6, duration=2400),
    dict(id="swiggum-2022", name="Swiggum et al. (2022)", show=("radcliffe", "Clusters (< 15 Myr)", "Sun"),
         dust=0.75, time=0.0, fit="rw", toward=(0.0, -0.02, 1.0), margin=1.1, quantile=1.0, duration=2400,
         styles={"Clusters (< 15 Myr)": {"color": "#8fd0ff"}}),
    dict(id="swiggum-2024", name="Swiggum et al. (2024)", show=("families", "Sun"), dust=None, time=-30.0,
         fit="families", toward=(0.1, -1.0, 0.75), margin=1.5, duration=3600, orbit=True),
    dict(id="swiggum-2025", name="Swiggum et al. (2025)", show=("branches", "families", "Sun"), dust=None,
         time=-40.0, fit="branches", toward=(0.1, -1.0, 1.1), margin=1.05, quantile=0.98, box=True, duration=3600,
         styles={"branches": {"sizeScale": 2.2}}),
    dict(id="swiggum-2026", name="Swiggum et al. (2026)", show=("supernovae", "Sun"), dust=1.0, time=0.0,
         fit=None, view={"camera": {"position": {"x": 0.0, "y": -4300.0, "z": 4700.0}, "target": {"x": 0.0, "y": 0.0, "z": 0.0}},
                         "fov": 45.0}, duration=2400),
)


def layer_sets(spec: dict) -> dict[str, list[str]]:
    names = trace_keys_by_name(spec)
    keys = [i["key"] for i in spec["legend"]["items"] if i.get("kind", "trace") == "trace"]
    return {
        "orion": [k for k in keys if k.startswith(KEYS["orion_stars"]) or k == KEYS["orion_motion"]],
        "radcliffe": [KEYS["rw_model"], KEYS["masers"], KEYS["arm_fit"]],
        "rw": [KEYS["rw_model"]],
        "families": [names[n] for n in ("Alpha Persei Family", "Cr 135 Family", "M6 Family")],
        "branches": [k for k in keys if k.startswith("paper-2025-")],
        **{n: [names[n]] for n in names},
    }


def paper_states(spec: dict) -> dict:
    sets = layer_sets(spec)
    traces = [i["key"] for i in spec["legend"]["items"] if i.get("kind", "trace") == "trace"]
    volumes = [i["key"] for i in spec["legend"]["items"] if i.get("kind") == "volume"]
    dust = spec["initial_state"]["volume_state_by_key"].get("volume-0", {})
    items = []
    for view in PAPER_VIEWS:
        shown = {k for name in view["show"] if name != "supernovae" for k in sets[name]}
        legend_state = {k: k in shown for k in traces}
        sne = "supernovae" in view["show"]
        legend_state.update({"volume-0": view["dust"] is not None, KEYS["supernovae"]: sne})
        volume_state = {"volume-0": {"visible": view["dust"] is not None,
                                     "opacity": float(dust.get("opacity", 1.0)) * (view["dust"] or 1.0)},
                        KEYS["supernovae"]: {**SUPERNOVA_STYLE, "visible": sne}}
        camera = view.get("view") or fit_view(spec, sets[view["fit"]], view["time"], view["toward"],
                                              margin=view.get("margin", 1.25), quantile=view.get("quantile", 0.9),
                                              box=view.get("box", False))
        styles = {k: style for name, style in (view.get("styles") or {}).items() for k in sets[name]}
        controls = {"camera_fov": camera["fov"], "camera_auto_orbit_enabled": bool(view.get("orbit")),
                    "camera_auto_orbit_speed": 0.35}
        items.append({
            "id": view["id"], "name": view["name"], "camera": "follow",
            "transition": {"duration_ms": view["duration"], "easing": "easeInOutCubic"},
            "state": {"current_group": GROUP, "legend_state": legend_state, "volume_state_by_key": volume_state,
                      "trace_style_state": styles, "global_controls": controls, "camera": camera["camera"],
                      "current_frame_value": frame_of(spec, view["time"])},
        })
    unknown = set(volumes) - {"volume-0", KEYS["supernovae"]}
    if unknown:
        raise RuntimeError(f"Volumes without a saved-view setting: {sorted(unknown)}")
    return {"schema_version": 1, "project_id": "cam-website-background", "default_mode": "edit",
            "present_only": False, "autosave_drafts": False, "preserve_camera_on_navigation": False,
            "default_transition": {"duration_ms": 1600, "easing": "easeInOutCubic"}, "items": items, "assets": {}}


# ---------------------------------------------------------------- build

def build_spec(oct1_source: Path = oct1.SOURCE_HTML) -> dict:
    spec = hash_provenance_paths(read_legacy_scene_spec(Path(oct1_source).expanduser()))
    provenance = spec.setdefault("provenance", {})
    provenance["website_base"] = {"source": Path(oct1_source).name, **simplify(spec)}
    provenance["swiggum_2021"] = add_orion_layer(spec)
    provenance["swiggum_2022"] = add_radcliffe_layer(spec)
    provenance["swiggum_2025"] = add_moving_group_layer(spec)
    provenance["swiggum_2026"] = add_supernova_layer(spec)
    spec["states"] = paper_states(spec)
    return spec


def build_figure(output_html: Path = OUTPUT_HTML, *, oct1_source: Path = oct1.SOURCE_HTML) -> dict:
    spec = build_spec(oct1_source)
    bundle = compile_scene_spec(spec)
    out = OvizFigure(bundle=bundle, mode="focus", camera_anchor="lsr", embed=True,
                     title="Star clusters near the Sun").write_html(output_html)
    print(f"Wrote {out} ({out.stat().st_size / 2**20:.1f} MiB)")
    return {"output_bytes": out.stat().st_size, **bundle.size_report(),
            "provenance": {k: spec["provenance"][k] for k in spec["provenance"] if k.startswith("swiggum")}}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output-html", type=Path, default=OUTPUT_HTML)
    parser.add_argument("--oct1-source", type=Path, default=oct1.SOURCE_HTML,
                        help="the October 1 classic source (tests/main_figure_oct1.py writes it)")
    args = parser.parse_args()
    report = build_figure(args.output_html, oct1_source=args.oct1_source)
    print(json.dumps({k: report[k] for k in ("output_bytes", "compile_s", "write_s", "provenance") if k in report},
                     indent=1, default=str))


if __name__ == "__main__":
    main()
