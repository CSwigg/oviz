#!/usr/bin/env python3
"""Regenerate a July 21 Oviz figure with State-only Presentation mode."""

from __future__ import annotations

import argparse
import base64
import gzip
import json
import re
from copy import deepcopy
from pathlib import Path

from oviz.threejs_figure import ThreeJSFigure


SOURCE_HTML = Path(__file__).with_name("main_figure_chronos_july4.html")
DEFAULT_OUTPUT_HTML = Path(__file__).with_suffix(".html")


def read_embedded_scene_spec(path: Path) -> dict:
    html = Path(path).read_text(encoding="utf-8")
    payload_id = None
    for payload_arg in re.findall(r"readOvizSceneSpecPayload\((.*?)\)", html):
        try:
            candidate = json.loads(payload_arg)
        except json.JSONDecodeError:
            continue
        if isinstance(candidate, str) and candidate:
            payload_id = candidate
            break
    if not payload_id:
        raise ValueError(f"Could not locate the compressed scene payload in {path}.")
    chunks = re.findall(
        (
            r"<script\b(?=[^>]*type=[\"']application/octet-stream[\"'])"
            r"(?=[^>]*data-oviz-payload-id=[\"']"
            + re.escape(payload_id)
            + r"[\"'])[^>]*data-oviz-payload-index=[\"'](\d+)[\"'][^>]*>(.*?)</script>"
        ),
        html,
        re.S,
    )
    if not chunks:
        raise ValueError(f"Could not locate payload chunks for {payload_id}.")
    encoded = "".join(
        re.sub(r"\s+", "", chunk)
        for _, chunk in sorted(chunks, key=lambda item: int(item[0]))
    )
    return json.loads(gzip.decompress(base64.b64decode(encoded)))


def build_state_only_scene(scene_spec: dict) -> dict:
    scene = deepcopy(scene_spec)
    scene["title"] = ""
    scene["deck"] = {
        "schema_version": 2,
        "available": False,
        "enabled": False,
        "embedded": False,
        "revision": 0,
        "aspect_ratio": "16:9",
        "guides": {"smart": True, "grid": False, "grid_size": 20},
        "slides": [],
    }
    return scene


def build_figure(source_html: Path, output_html: Path) -> Path:
    source_html = Path(source_html).expanduser().resolve()
    output_html = Path(output_html).expanduser().resolve()
    if source_html == output_html:
        raise ValueError("The July 21 figure must be a sibling copy, not the source artifact.")
    scene = build_state_only_scene(read_embedded_scene_spec(source_html))
    html = ThreeJSFigure(scene, compress_scene_spec=True).to_html(compress_scene_spec=True)
    output_html.parent.mkdir(parents=True, exist_ok=True)
    output_html.write_text(html, encoding="utf-8")
    return output_html


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-html", type=Path, default=SOURCE_HTML)
    parser.add_argument("--output-html", type=Path, default=DEFAULT_OUTPUT_HTML)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    output = build_figure(args.source_html, args.output_html)
    print(f"Wrote {output}")


if __name__ == "__main__":
    main()
