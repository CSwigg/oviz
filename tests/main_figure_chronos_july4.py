#!/usr/bin/env python3
"""Render the main figure with July 4 PARSEC Chronos ages and updated velocities."""

from __future__ import annotations

import argparse
from pathlib import Path

from main_figure_new_chronos import run_main_figure


JULY4_CHRONOS_RESULTS_PATH = Path(
    "/Users/swiggumc/Desktop/astro_research/chronos_fasrc/runs/current/chronos/"
    "parsec_allhunt_46w_500b_5000s_dustav_12gyr_linearage_192shards/cluster_results.csv"
)
DEFAULT_OUTPUT_HTML = Path(__file__).resolve().with_suffix(".html")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-html",
        type=Path,
        default=DEFAULT_OUTPUT_HTML,
        help="HTML output path for the rendered July 4 Chronos figure.",
    )
    parser.add_argument(
        "--include-spiral-arms",
        action="store_true",
        help="Include the rotating spiral arm model traces.",
    )
    parser.add_argument(
        "--full-resolution",
        action="store_true",
        help="Render the original full-resolution payload instead of the upload/iPhone-safe payload.",
    )
    parser.add_argument(
        "--force-mobile-mode",
        action="store_true",
        help="Force mobile controls on every device. By default, mobile controls are auto-detected at runtime.",
    )
    parser.add_argument(
        "--full-payload",
        action="store_true",
        help="Keep repeated per-frame hover, selection, and motion metadata in the HTML.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    output_html = run_main_figure(
        output_html=args.output_html,
        mobile_mode=bool(args.force_mobile_mode),
        compact_payload=not bool(args.full_payload),
        mobile_safe_mode=not bool(args.full_resolution),
        chronos_results_path=JULY4_CHRONOS_RESULTS_PATH,
        chronos_model="parsec",
        include_spiral_arms=bool(args.include_spiral_arms),
        jun6_catalog=True,
        website_output_html=None,
    )
    print(f"Wrote {output_html}")


if __name__ == "__main__":
    main()
