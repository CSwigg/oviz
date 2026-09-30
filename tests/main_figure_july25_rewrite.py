"""Rebuild the July 25 figure in the Oviz viewer.

The science is unchanged. The scene spec embedded in the classic
``main_figure_july25.html`` (written by ``main_figure_july25.py``) is compiled
into the viewer's bundle, saved States included, and written next to it as
``main_figure_july25_rewrite.html``.

The camera opens anchored to the Local Standard of Rest. The figure's frame is
LSR-centred: fixed Galactic axes that translate with the LSR's orbit, so the
LSR is the origin at every time. The camera orbits and zooms about the LSR
and stays with it through time, while the Sun and the clusters move.
Readers can anchor it to the Sun or to any cluster instead (Display settings,
or Follow in the inspector).

    python tests/main_figure_july25_rewrite.py
"""

from __future__ import annotations

import argparse
from pathlib import Path

from oviz.viewer.upgrade import upgrade_html

HERE = Path(__file__).resolve().parent
SOURCE_HTML = HERE / "main_figure_july25.html"
OUTPUT_HTML = HERE / "main_figure_july25_rewrite.html"


def build_figure(source_html: Path = SOURCE_HTML, output_html: Path = OUTPUT_HTML, *, mode: str = "focus") -> dict:
    """Write the July 25 figure in the Oviz viewer, anchored to the LSR."""

    source_html = Path(source_html).expanduser().resolve()
    output_html = Path(output_html).expanduser().resolve()
    if source_html == output_html:
        raise ValueError("Refusing to overwrite the classic source figure.")
    return upgrade_html(source_html, output_html, verbose=True, mode=mode, camera_anchor="lsr")


def main() -> None:
    parser = argparse.ArgumentParser(description="Rebuild the July 25 figure in the Oviz viewer.")
    parser.add_argument("--source-html", type=Path, default=SOURCE_HTML)
    parser.add_argument("--output-html", type=Path, default=OUTPUT_HTML)
    parser.add_argument("--mode", choices=("focus", "detailed"), default="focus")
    args = parser.parse_args()
    report = build_figure(args.source_html, args.output_html, mode=args.mode)
    print(f"Wrote {args.output_html} ({report['output_bytes'] / 2**20:.1f} MiB)")


if __name__ == "__main__":
    main()
