#!/usr/bin/env python3
"""Render a smaller self-contained main figure for sharing and mobile Safari."""

from __future__ import annotations

import argparse
from pathlib import Path

from main_figure import SHARE_OUTPUT_HTML, run_main_figure


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-html",
        type=Path,
        default=SHARE_OUTPUT_HTML,
        help="HTML output path for the share/mobile figure.",
    )
    parser.add_argument(
        "--theme-key",
        default=None,
        help="Optional Three.js color theme preset to force into the figure initial state.",
    )
    parser.add_argument(
        "--website-output-html",
        type=Path,
        default=None,
        help="Optional website copy target. By default, the share script only writes its output HTML.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    output_html = run_main_figure(
        output_html=args.output_html,
        theme_key=args.theme_key,
        mobile_mode=True,
        compact_payload=True,
        share_mode=True,
        website_output_html=args.website_output_html,
    )
    print(f"Wrote {output_html}")
    if args.website_output_html is not None:
        print(f"Copied {output_html} to {args.website_output_html}")


if __name__ == "__main__":
    main()
