#!/usr/bin/env python3
"""Render the main figure in minimal galactic presentation mode."""

from __future__ import annotations

import argparse
from pathlib import Path

from main_figure import run_main_figure


DEFAULT_OUTPUT_HTML = Path("/tmp/main_figure_simple_galactic.html")


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
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    output_html = run_main_figure(
        output_html=args.output_html,
        theme_key=args.theme_key,
        minimal_mode=True,
        galactic_simple=True,
    )
    print(f"Wrote {output_html}")


if __name__ == "__main__":
    main()
