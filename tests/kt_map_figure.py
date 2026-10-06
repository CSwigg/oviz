#!/usr/bin/env python3
"""The October 1 figure with a kinetic tomography (KT) map.

Everything in ``tests/main_figure_oct1.py`` (the July 25 clusters in the
Khalil et al. 2025 spiral potential, the dust and Hα volumes, Sky members,
the hidden arm and pulsar traces, the opening view) plus the DIB KT map
(``cubes_full_mean_std.h5``: diffuse-interstellar-band density, its
line-of-sight velocity and the residual after subtracting Galactic rotation,
on a 23.5 pc grid out to 6 kpc). The map is one volume that opens coloured
by velocity (red receding, blue approaching), with opacity following
density; readers switch it to the residual or to density. Two animated flow
layers run along the sight lines: the full velocity (shown) and the residual
(in the key, hidden). Published as ``oviz_figures/oviz_kt_map.html``.

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

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parent
for path in (REPO_ROOT, HERE):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

import main_figure_oct1 as oct1  # noqa: E402
from figure_provenance import short_sha256  # noqa: E402
from oviz.kt import KTMap, read_kt_map  # noqa: E402
from oviz.viz import _build_threejs_inline_volume_layer_spec  # noqa: E402

DEFAULT_KT = Path.home() / "Downloads" / "cubes_full_mean_std.h5"
DEFAULT_OUTPUT = Path(__file__).resolve().with_suffix(".html")
KT_KEY = "kt-dib"


def add_kt_map(spec: dict, kt: KTMap, *, key: str = KT_KEY, max_resolution: int = 256,
               max_lines: int = 1500, source: Path | None = None) -> None:
    """Add ``kt`` to a scene spec as a shown volume (in every group's key) and two flow layers."""
    layers = spec["volumes"]["layers"]
    layer = _build_threejs_inline_volume_layer_spec(kt.volume(key=key, max_resolution=max_resolution), index=len(layers))
    layers.append(layer)
    spec["legend"]["items"].append({"key": key, "name": layer["name"], "color": layer["legend_color"], "kind": "volume"})
    for visibility in spec["group_visibility"].values():
        visibility[key] = "legendonly"
    spec.setdefault("initial_state", {}).setdefault("legend_state", {})[key] = True
    flows = [kt.flows(max_lines=max_lines)]
    if "residual" in kt.fields:
        flows.append(kt.flows(field="residual", max_lines=max_lines, visible=False))
    spec["flows"] = [*(spec.get("flows") or []), *flows]
    spec.setdefault("provenance", {})["kt_map"] = {
        **({"file": source.name, "file_sha256": short_sha256(source)} if source else {}),
        "density": "mean/ew_cube", "velocity": kt.attrs.get("velocity_dataset"),
        "residual": kt.attrs.get("residual_dataset"), "residual_model": kt.residual_note,
        "flow_lines": {f["key"]: f["summary"]["lines"] for f in flows},
    }


def build(kt_path: Path, output: Path, *, max_resolution: int = 256, max_lines: int = 1500) -> Path:
    t0 = time.perf_counter()
    kt = read_kt_map(kt_path, name="DIB KT map")
    report = oct1.build_figure(
        output_html=output,
        extend_spec=lambda spec: add_kt_map(spec, kt, max_resolution=max_resolution, max_lines=max_lines, source=kt_path),
    )
    print(f"KT figure in {time.perf_counter() - t0:.1f} s ({report['output_bytes'] / 2**20:.1f} MiB)")
    return output


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--kt", type=Path, default=DEFAULT_KT, help="KT map HDF5 file")
    ap.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    ap.add_argument("--max-resolution", type=int, default=256)
    ap.add_argument("--max-lines", type=int, default=1500)
    args = ap.parse_args(argv)
    if not args.kt.expanduser().exists():
        print(f"no KT map at {args.kt}; pass --kt", file=sys.stderr)
        return 1
    build(args.kt.expanduser(), args.output, max_resolution=args.max_resolution, max_lines=args.max_lines)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
