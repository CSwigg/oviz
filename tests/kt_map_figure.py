#!/usr/bin/env python3
"""A kinetic tomography (KT) map in the Oviz viewer.

The DIB KT map (``cubes_full_mean_std.h5``: diffuse-interstellar-band
density, its line-of-sight velocity and the residual after subtracting
Galactic rotation, on a 23.5 pc grid out to 6 kpc) is drawn as one volume,
opening coloured by velocity (red receding, blue approaching) with opacity
following density; readers switch it to the residual or to density. Two
animated flow layers run along the sight lines: the full velocity (shown)
and the residual (in the key, hidden). Young clusters from the Figure 3
sample give scale when their catalogue is on this machine.

    python tests/kt_map_figure.py                       # ~/Downloads/cubes_full_mean_std.h5
    python tests/kt_map_figure.py --kt other_kt.h5 --output kt.html
"""

from __future__ import annotations

import argparse
import os
import sys
import time
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl")
os.environ.setdefault("MPLBACKEND", "Agg")

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from oviz import Scene3D, Trace, TraceCollection  # noqa: E402
from oviz.kt import read_kt_map  # noqa: E402

DEFAULT_KT = Path.home() / "Downloads" / "cubes_full_mean_std.h5"
CLUSTERS = (Path.home() / "Desktop" / "astro_research" / "supernovae_map" / "paper" / "analysis"
            / "Figure_3_full_cluster_sample_clusters.csv")
DEFAULT_OUTPUT = Path(__file__).resolve().with_suffix(".html")


def cluster_traces(max_age_myr: float = 30.0) -> list[Trace]:
    """Young clusters for scale, or a lone marker at the Sun's place when the catalogue is absent."""
    if CLUSTERS.exists():
        df = pd.read_csv(CLUSTERS)
        young = df[df["age_myr"] < max_age_myr].copy()
        return [Trace(young[["name", "x", "y", "z", "U", "V", "W", "age_myr", "n_stars"]],
                      data_name=f"Clusters (< {max_age_myr:.0f} Myr)", color="#ffd166", size_by_n_stars=True)]
    sun = pd.DataFrame({"name": ["Sun"], "x": [0.0], "y": [0.0], "z": [0.0], "U": [0.0], "V": [0.0], "W": [0.0], "age_myr": [4600.0]})
    return [Trace(sun, data_name="Sun", color="#ffd166")]


def build(kt_path: Path, output: Path, *, max_resolution: int = 256, max_lines: int = 1500) -> Path:
    t0 = time.perf_counter()
    kt = read_kt_map(kt_path, name="DIB KT map")
    volume = kt.volume(max_resolution=max_resolution)
    flows = [kt.flows(max_lines=max_lines)]
    if "residual" in kt.fields:
        flows.append(kt.flows(field="residual", max_lines=max_lines, visible=False))
    print(f"KT map {kt.density.shape[::-1]} read and prepared in {time.perf_counter() - t0:.1f} s "
          f"({' + '.join(str(f['summary']['lines']) for f in flows)} flow lines)")
    scene = Scene3D(TraceCollection(cluster_traces()), figure_theme="dark")
    figure = scene.make_plot(
        time=np.arange(0, -31, -1.0),
        galactic_mode=True,
        enable_sky_panel=True,
        volumes=[volume],
        flows=flows,
        fade_in_time=0,
    )
    # Open on the whole map, tilted, rather than straight down from 30 kpc.
    figure.scene_spec.setdefault("initial_state", {})["camera"] = {"position": [6000.0, -6800.0, 9300.0], "target": [0.0, 0.0, 0.0]}
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.write_html(str(output))
    print(f"wrote {output} ({output.stat().st_size / 2**20:.1f} MiB) in {time.perf_counter() - t0:.1f} s")
    return output


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--kt", type=Path, default=DEFAULT_KT, help="KT map HDF5 file")
    ap.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    ap.add_argument("--max-resolution", type=int, default=256)
    ap.add_argument("--max-lines", type=int, default=1500)
    args = ap.parse_args(argv)
    if not args.kt.exists():
        print(f"no KT map at {args.kt}; pass --kt", file=sys.stderr)
        return 1
    build(args.kt.expanduser(), args.output, max_resolution=args.max_resolution, max_lines=args.max_lines)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
