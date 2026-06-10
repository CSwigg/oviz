#!/usr/bin/env python3
"""Render the main figure with Oviz mobile-mode metadata enabled."""

from __future__ import annotations

import argparse
from pathlib import Path

from main_figure import run_main_figure


DEFAULT_OUTPUT_HTML = Path(__file__).resolve().with_suffix(".html")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-html",
        type=Path,
        default=DEFAULT_OUTPUT_HTML,
        help="HTML output path for the rendered mobile figure.",
    )
    parser.add_argument(
        "--theme-key",
        default=None,
        help="Optional Three.js color theme preset to force into the figure initial state.",
    )
    parser.add_argument(
        "--minimal-mode",
        action="store_true",
        help="Also render with the existing minimal presentation profile.",
    )
    parser.add_argument(
        "--full-payload",
        action="store_true",
        help="Keep repeated per-frame hover, selection, and motion metadata in the HTML.",
    )
    parser.add_argument(
        "--website-output-html",
        type=Path,
        default=None,
        help="Optional website copy target. By default, the mobile script only writes its output HTML.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    output_html = run_main_figure(
        output_html=args.output_html,
        theme_key=args.theme_key,
        minimal_mode=bool(args.minimal_mode),
        mobile_mode=True,
        compact_payload=not bool(args.full_payload),
        website_output_html=args.website_output_html,
    )
    print(f"Wrote {output_html}")
    if args.website_output_html is not None:
        print(f"Copied {output_html} to {args.website_output_html}")


if __name__ == "__main__":
    main()
