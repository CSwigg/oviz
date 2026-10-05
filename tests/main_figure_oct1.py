#!/usr/bin/env python3
"""The October 1 figure: the July 25 clusters moving in a spiral-arm potential.

Science (slow; ``--rebuild-source``, or when the source is missing): the July
25 source scene (``tests/main_figure_july25.py``: the same catalogs, Chronos
ages, velocities, Sky members and volumes) with the cluster orbits integrated
in MWPotential2014 plus the Khalil et al. (2025) spiral modes
(``oviz.spiral_models.khalil2025_potential``). Two traces start hidden: the
Khalil et al. (2025) arms and the Castro-Ginard et al. (2021) arms, each
turning at its pattern speeds through the timeline. It is written as the
classic source figure ``tests/main_figure_oct1_source.html``.

Viewer (fast, always): that source in the Oviz viewer, opening in focus mode
with the camera anchored to the LSR, as ``tests/main_figure_oct1.html``. It is
published as ``oviz_figures/oviz_oct1.html``.

The viewer step also adds two pulsar traces, both hidden by default
(``--no-pulsars`` leaves them out): the ATNF v2.8.1 pulsars with a measured
proper motion (the web export ``PULSAR_TABLE``) within 3 kpc, with
characteristic age tau_c < 1 Myr and 1-5 Myr. Each moves back along its
plane-of-sky motion relative to local circular rotation (radial velocity
unknown, taken as zero) to t = -tau_c, and sits at that birth site, drifting
with the disk, before it.

    python tests/main_figure_oct1.py                   # viewer from the stored source
    python tests/main_figure_oct1.py --rebuild-source  # re-run the science first

``--baseline OUT`` rebuilds the unchanged July 25 science (MWPotential2014, no
arms) as a classic figure, to check the pipeline against
``tests/main_figure_july25.html``.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import main_figure_july25 as july  # noqa: E402
import main_figure_new_chronos as source_runner  # noqa: E402
from oviz.spiral_models import CASTRO_GINARD2021_ARMS, KHALIL2025_ARMS, provenance  # noqa: E402
from oviz.threejs_embed import scene_spec_jsonable  # noqa: E402
from oviz.threejs_figure import ThreeJSFigure  # noqa: E402
from oviz.viewer.compile import compile_scene_spec  # noqa: E402
from oviz.viewer.figure import OvizFigure  # noqa: E402
from oviz.viewer.upgrade import read_legacy_scene_spec, upgrade_html  # noqa: E402

SOURCE_HTML = HERE / "main_figure_oct1_source.html"
OUTPUT_HTML = HERE / "main_figure_oct1.html"
ARM_MODELS = (KHALIL2025_ARMS, CASTRO_GINARD2021_ARMS)
# ATNF v2.8.1 web query: Name, PMRA, PMDec, Gl, Gb, P0, P1, Dist, Age for every pulsar with exist(PMRA).
PULSAR_TABLE = Path("~/Desktop/astro_research/cfa/sfh_nifty/data/pulsars/psrcat_v2.8.1_pm_only.txt").expanduser()
PULSAR_MAX_DIST_KPC = 3.0
PULSAR_TRACES = (  # key, legend name, tau_c range [Myr], colour, point size
    ("trace-pulsars-young", "Pulsars with PM (τc < 1 Myr)", (0.0, 1.0), "#5cff8d", 8.0),
    ("trace-pulsars-older", "Pulsars with PM (τc 1–5 Myr)", (1.0, 5.0), "#c9d1dc", 6.0),
)
V_SUN_KMS = (11.1, 12.24, 7.25)  # Schoenrich et al. (2010), as the Sun's track in the scene
K_PM = 4.740470463533348  # km/s per (mas/yr x kpc)
KMS_TO_PC_PER_MYR = 1.0227121650537077


def _replace_once(pattern: str, replacement: str, source: str, what: str) -> str:
    source, count = re.subn(pattern, replacement, source, count=1)
    if count != 1:
        raise RuntimeError(f"Could not {what} in the converted notebook script.")
    return source


def patch_source(source: str, *, spiral: bool = True) -> str:
    """The July 25 notebook script, run in the spiral potential with the arm traces."""
    source = _replace_once(r"(?m)^save_name\s*=\s*['\"][^'\"]+['\"]", "save_name = None", source,
                           "skip the notebook's own HTML write")
    if not spiral:
        return source
    source = _replace_once(
        r"(?m)^plot_3d\s*=\s*Animate3D\(",
        "from oviz.spiral_models import CASTRO_GINARD2021_ARMS, KHALIL2025_ARMS, khalil2025_potential\n"
        "plot_3d = Animate3D(\n    potential=khalil2025_potential(),",
        source, "integrate the orbits in the spiral potential",
    )
    return _replace_once(
        r"(?m)^(\s*)include_spiral_arms\s*=\s*False\s*,\s*$",
        r"\1include_spiral_arms=False,\n\1spiral_arm_models=(KHALIL2025_ARMS, CASTRO_GINARD2021_ARMS),",
        source, "add the spiral-arm traces",
    )


def build_source_scene(*, spiral: bool = True, mode_ages_path: Path | None = None) -> dict:
    """Run the July 25 science (in the spiral potential unless ``spiral=False``)."""
    if mode_ages_path is not None:
        # A local copy of the June 6 Chronos mode-age table, for when the
        # configured file is unavailable (e.g. offloaded by iCloud).
        source_runner.JUN6_CHRONOS_MODE_AGES_PATH = Path(mode_ages_path).expanduser().resolve()
    kwargs = july.source_main_figure_kwargs(
        july.DEFAULT_CLUSTER_VELOCITIES_PATH,
        july.DEFAULT_RATZENBOECK_SCOCEN_CATALOG_PATH,
        july.DEFAULT_RATZENBOECK_SCOCEN_METADATA_PATH,
    )
    source = source_runner.build_patched_main_figure_source(
        SOURCE_HTML.with_name("main_figure_oct1_unused.html"), **kwargs,
    )
    namespace = source_runner.run_script_source(patch_source(source, spiral=spiral))
    scene = july.build_state_only_scene(
        scene_spec_jsonable(namespace["fig3d"].to_dict()),
        cluster_velocities_path=july.DEFAULT_CLUSTER_VELOCITIES_PATH,
        ratzenboeck_scocen_catalog_path=july.DEFAULT_RATZENBOECK_SCOCEN_CATALOG_PATH,
        ratzenboeck_scocen_metadata_path=july.DEFAULT_RATZENBOECK_SCOCEN_METADATA_PATH,
        ratzenboeck_scocen_members_path=july.DEFAULT_RATZENBOECK_SCOCEN_MEMBERS_PATH,
    )
    if spiral:
        scene.setdefault("provenance", {})["orbit_potential"] = provenance(ARM_MODELS)
    return scene


def write_classic(scene: dict, path: Path) -> Path:
    path = Path(path).expanduser().resolve()
    if path in {july.DEFAULT_OUTPUT_HTML.resolve(), july.SOURCE_HTML.resolve()}:
        raise ValueError("Refusing to overwrite a July artifact.")
    html = ThreeJSFigure(scene, compress_scene_spec=True).to_html(compress_scene_spec=True)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(html, encoding="utf-8")
    return path


def read_pulsar_table(path: Path):
    """The ATNF web export (``;``-separated value/error/reference columns, ``*`` = no value)."""
    import numpy as np
    import pandas as pd

    rows = []
    for line in Path(path).read_text().splitlines():
        f = line.split(";")
        if f[0].strip().isdigit():
            rows.append(dict(name=f[1], pmra=f[3], pmdec=f[6], gl=f[9], gb=f[10], dist=f[17], age_yr=f[18]))
    psr = pd.DataFrame(rows).replace("*", np.nan)
    for col in psr.columns.drop("name"):
        psr[col] = pd.to_numeric(psr[col], errors="coerce")
    psr["age_myr"] = psr["age_yr"] / 1e6
    lr, br = np.deg2rad(psr["gl"]), np.deg2rad(psr["gb"])
    psr["x_pc"] = 1000.0 * psr["dist"] * np.cos(br) * np.cos(lr)
    psr["y_pc"] = 1000.0 * psr["dist"] * np.cos(br) * np.sin(lr)
    psr["z_pc"] = 1000.0 * psr["dist"] * np.sin(br)
    return psr


def pulsar_tracks(psr, times, ro_kpc: float, vo_kms: float):
    """Positions (N, F, 3) [pc] at the frame times [Myr] in the scene's frame rotating with the LSR.

    Velocity: proper motion at the catalogue distance plus the Sun's motion, minus local circular
    rotation (MWPotential2014 Oort constants), with the radial part dropped (the pulsar's radial
    velocity is unknown). In the rotating frame circular orbits also drift by 2A x along +Y.
    Back to t = -tau_c the pulsar moves; before it, it sits at the birth site drifting with the disk.
    """
    import numpy as np
    from astropy import units as u
    from astropy.coordinates import SkyCoord
    from galpy.potential import MWPotential2014, dvcircdR

    omega0 = vo_kms / ro_kpc
    dvdr = float(dvcircdR(MWPotential2014, 1.0)) * vo_kms / ro_kpc
    oort_a, oort_b = 0.5 * (omega0 - dvdr), -0.5 * (omega0 + dvdr)
    icrs = SkyCoord(l=psr["gl"].to_numpy() * u.deg, b=psr["gb"].to_numpy() * u.deg, frame="galactic").icrs
    pm = SkyCoord(ra=icrs.ra, dec=icrs.dec, pm_ra_cosdec=psr["pmra"].to_numpy() * u.mas / u.yr,
                  pm_dec=psr["pmdec"].to_numpy() * u.mas / u.yr, frame="icrs").galactic
    mul, mub = pm.pm_l_cosb.to_value(u.mas / u.yr), pm.pm_b.to_value(u.mas / u.yr)
    lr, br = np.deg2rad(psr["gl"].to_numpy()), np.deg2rad(psr["gb"].to_numpy())
    rhat = np.stack([np.cos(br) * np.cos(lr), np.cos(br) * np.sin(lr), np.sin(br)], axis=-1)
    e_l = np.stack([-np.sin(lr), np.cos(lr), np.zeros_like(lr)], axis=-1)
    e_b = np.stack([-np.sin(br) * np.cos(lr), -np.sin(br) * np.sin(lr), np.cos(br)], axis=-1)
    r0 = psr[["x_pc", "y_pc", "z_pc"]].to_numpy()
    v_circ = np.stack([(oort_a - oort_b) * r0[:, 1], (oort_a + oort_b) * r0[:, 0], 0.0 * r0[:, 0]], axis=-1) / 1000.0
    v = K_PM * psr["dist"].to_numpy()[:, None] * (mul[:, None] * e_l + mub[:, None] * e_b) + V_SUN_KMS - v_circ
    v -= (v * rhat).sum(-1, keepdims=True) * rhat
    drift = lambda r: np.stack([0.0 * r[..., 0], 2.0 * oort_a * r[..., 0] / 1000.0, 0.0 * r[..., 0]], axis=-1)
    tau = psr["age_myr"].to_numpy()
    birth = r0 - (v + drift(r0)) * tau[:, None] * KMS_TO_PC_PER_MYR
    t = np.asarray(times, float)[None, :, None]
    moving = r0[:, None, :] + (v + drift(r0))[:, None, :] * t * KMS_TO_PC_PER_MYR
    resting = birth[:, None, :] + drift(birth)[:, None, :] * (t + tau[:, None, None]) * KMS_TO_PC_PER_MYR
    pos = np.where(t >= -tau[:, None, None], moving, resting)
    return pos, np.linalg.norm(v, axis=-1), {"oort_A_kms_kpc": oort_a, "oort_B_kms_kpc": oort_b}


def add_pulsar_traces(spec: dict, table_path: Path = PULSAR_TABLE) -> dict[str, int]:
    """Add the proper-motion pulsar traces to every frame of ``spec``, hidden in every group."""
    from astropy import units as u
    from astropy.coordinates import SkyCoord

    psr = read_pulsar_table(table_path)
    psr = psr[(psr["age_myr"] < max(hi for _, _, (_, hi), _, _ in PULSAR_TRACES)) & (psr["dist"] <= PULSAR_MAX_DIST_KPC)
              & psr[["pmra", "pmdec", "gl", "gb"]].notna().all(axis=1)].sort_values("age_myr").reset_index(drop=True)
    radec = SkyCoord(l=psr["gl"].to_numpy() * u.deg, b=psr["gb"].to_numpy() * u.deg, frame="galactic").icrs
    frames = spec["frames"]
    t0 = int(spec["initial_frame_index"])
    if float(frames[t0]["time"]) != 0.0:
        raise RuntimeError("The initial frame is not t = 0.")
    orbit = (spec.get("provenance") or {}).get("orbit_potential") or {}
    pos, v_sky, oort = pulsar_tracks(psr, [float(f["time"]) for f in frames],
                                     float(orbit.get("ro_kpc", 8.122)), float(orbit.get("vo_kms", 236.0)))
    legend = spec["legend"]["items"]
    insert_at = max(i for i, item in enumerate(legend) if item.get("kind", "trace") == "trace") + 1
    counts = {}
    for key, name, (lo, hi), color, size in PULSAR_TRACES:
        rows = psr.index[(psr["age_myr"] >= lo) & (psr["age_myr"] < hi)]
        color_by = {"mode": "age", "label": "τc (Myr)", "cmin": lo, "cmax": hi, "colormap": "turbo",
                    "legend_color": color, "default_color_mode": "fixed"}
        for fi, frame in enumerate(frames):
            points = []
            for i in rows:
                r = psr.loc[i]
                label = f"PSR {r['name']}"
                x, y, z = (round(float(c), 2) for c in pos[i, fi])
                point = {"x": x, "y": y, "z": z, "color_scalar": float(r["age_myr"]), "n_stars": 1.0,
                         "motion": {"key": label, "age_now_myr": float(r["age_myr"]), "time_myr": float(frame["time"])}}
                if fi == t0:
                    point["hovertext"] = (
                        f'<b style="font-size:16px;">{label}</b><br>{name}<br>τc = {r["age_myr"]:.3g} Myr'
                        f'<br>d = {r["dist"]:.2f} kpc<br>μα* = {r["pmra"]:g}, μδ = {r["pmdec"]:g} mas/yr'
                        f"<br>sky-plane speed relative to Galactic rotation = {v_sky[i]:.0f} km/s"
                        "<br>track: proper motion only (radial velocity unknown), back to t = −τc"
                        f"<br>(x,y,z) = ({x:.1f}, {y:.1f}, {z:.1f})")
                    point["selection"] = {
                        "l_deg": float(r["gl"]), "b_deg": float(r["gb"]), "dist_pc": float(r["dist"]) * 1000.0,
                        "x0": x, "y0": y, "z0": z, "age_now_myr": float(r["age_myr"]),
                        "age_at_t_myr": float(r["age_myr"]), "click_time_myr": 0.0,
                        "ra_deg": float(radec.ra.deg[i]), "dec_deg": float(radec.dec.deg[i]), "cluster_name": label,
                        "cluster_color": color, "n_stars": 1.0, "name_all": "", "trace_name": name,
                    }
                points.append(point)
            frame["traces"].append({
                "key": key, "name": name, "showlegend": True, "size_by_n_stars_default": False, "points": points,
                "default_point_size": size, "has_n_stars": False, "default_opacity": 0.95, "color_by": color_by,
                "legend_color": color, "default_color": color, "default_symbol": "diamond",
                "default_color_scalar_kind": "age",
            })
        legend.insert(insert_at, {
            "key": key, "name": name, "color": color, "kind": "trace", "has_points": True, "has_segments": False,
            "has_labels": False, "has_n_stars": False, "size_by_n_stars_default": False, "default_color": color,
            "default_opacity": 0.95, "default_point_size": size, "color_by": color_by,
        })
        insert_at += 1
        for visibility in spec["group_visibility"].values():
            visibility[key] = "legendonly"
        counts[name] = len(rows)
    spec.setdefault("provenance", {})["pulsars"] = {
        "table": str(table_path),
        "query": "ATNF v2.8.1 web query, condition exist(PMRA): Name, PMRA, PMDec, Gl, Gb, P0, P1, Dist, Age",
        "selection": f"tau_c < {max(hi for _, _, (_, hi), _, _ in PULSAR_TRACES):g} Myr, catalogue distance "
                     f"<= {PULSAR_MAX_DIST_KPC:g} kpc",
        "tracks": "proper motion at the catalogue distance relative to local circular rotation, radial velocity "
                  "taken as zero, linear back to t = -tau_c; birth site drifting with the disk before that",
        "galactic_constants": {**oort, "v_sun_kms": list(V_SUN_KMS)},
        "traces": counts,
    }
    return counts


def build_figure(source_html: Path = SOURCE_HTML, output_html: Path = OUTPUT_HTML, *,
                 rebuild_source: bool = False, mode_ages_path: Path | None = None,
                 pulsars: Path | None = PULSAR_TABLE) -> dict:
    """Write the October 1 figure, re-running the science first when asked or needed."""
    source_html = Path(source_html).expanduser().resolve()
    if rebuild_source or not source_html.exists():
        write_classic(build_source_scene(mode_ages_path=mode_ages_path), source_html)
        print(f"Wrote {source_html} ({source_html.stat().st_size / 2**20:.1f} MiB)")
    if pulsars is None:
        return upgrade_html(source_html, output_html, verbose=True, mode="focus", camera_anchor="lsr")
    # upgrade_html, with the pulsar traces added to the scene before it is compiled.
    spec = read_legacy_scene_spec(source_html)
    print("Pulsar traces:", add_pulsar_traces(spec, Path(pulsars).expanduser()))
    bundle = compile_scene_spec(spec)
    out = OvizFigure(bundle=bundle, mode="focus", camera_anchor="lsr").write_html(output_html)
    print(f"Wrote {out} ({out.stat().st_size / 2**20:.1f} MiB)")
    return {"output_bytes": out.stat().st_size, **bundle.size_report()}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--rebuild-source", action="store_true", help="re-run the science pipeline first")
    parser.add_argument("--source-html", type=Path, default=SOURCE_HTML)
    parser.add_argument("--output-html", type=Path, default=OUTPUT_HTML)
    parser.add_argument("--mode-ages", type=Path, default=None,
                        help="local copy of the June 6 Chronos mode-age table, if the configured one is unavailable")
    parser.add_argument("--baseline", type=Path, default=None, metavar="OUT",
                        help="write the unchanged July 25 science (no arms, MWPotential2014) here instead")
    parser.add_argument("--pulsars", type=Path, default=PULSAR_TABLE,
                        help="ATNF proper-motion export for the two hidden pulsar traces")
    parser.add_argument("--no-pulsars", action="store_true", help="leave the pulsar traces out")
    args = parser.parse_args()
    if args.baseline is not None:
        out = write_classic(build_source_scene(spiral=False, mode_ages_path=args.mode_ages), args.baseline)
        print(f"Wrote {out}")
        return
    report = build_figure(args.source_html, args.output_html, rebuild_source=args.rebuild_source,
                          mode_ages_path=args.mode_ages, pulsars=None if args.no_pulsars else args.pulsars)
    print(json.dumps({k: report[k] for k in ("output_bytes", "compile_s", "write_s") if k in report}))


if __name__ == "__main__":
    main()
