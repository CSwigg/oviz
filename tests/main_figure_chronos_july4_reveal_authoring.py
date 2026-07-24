#!/usr/bin/env python3
"""Create a Reveal-authoring copy of the compact July 4 Oviz figure."""

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
REVEAL_VERSION = "5.2.1"


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


def build_authoring_scene(scene_spec: dict) -> dict:
    scene = deepcopy(scene_spec)
    scene["title"] = "Young Clusters — Reveal Authoring"

    lookback_snapshot = deepcopy(scene.get("initial_state") or {})
    lookback_snapshot["current_frame_index"] = 115
    lookback_snapshot["current_frame_value"] = 115.0
    lookback_snapshot["playback_direction"] = 0
    lookback_snapshot["playback_speed"] = 1.0

    states = deepcopy(scene.get("states") or {})
    states["default_mode"] = "edit"
    states["embedded"] = True
    states["revision"] = max(int(states.get("revision") or 0), 1)
    states["synchronized_revision"] = states["revision"]
    states["items"] = [
        {
            "id": "deck-lookback-five-myr",
            "name": "Five million years ago",
            "transition": {"duration_ms": 1800, "easing": "easeInOutCubic"},
            "snapshot": lookback_snapshot,
            "degraded": False,
        }
    ]
    scene["states"] = states
    scene["deck"] = {
        "schema_version": 2,
        "enabled": True,
        "embedded": True,
        "revision": 1,
        "aspect_ratio": "16:9",
        "guides": {"smart": True, "grid": False, "grid_size": 20},
        "reveal": {
            "version": REVEAL_VERSION,
            "transition": "fade",
            "background_transition": "fade",
        },
        "slides": [
            {
                "id": "deck-opening",
                "name": "Opening",
                "state_id": None,
                "objects": [
                    {
                        "id": "deck-opening-kicker",
                        "kind": "text",
                        "type": "subtitle",
                        "text": "OVIZ",
                        "x": 112,
                        "y": 522,
                        "width": 672,
                        "height": 68,
                        "font_size": 24,
                        "font_weight": 700,
                        "color": "#8fdcff",
                        "align": "left",
                    },
                    {
                        "id": "deck-opening-title",
                        "kind": "text",
                        "type": "title",
                        "text": "Young Clusters in the Local Milky Way",
                        "x": 112,
                        "y": 576,
                        "width": 1072,
                        "height": 134,
                        "font_size": 66,
                        "font_weight": 700,
                        "color": "#ffffff",
                        "align": "left",
                    },
                    {
                        "id": "deck-opening-subtitle",
                        "kind": "text",
                        "type": "subtitle",
                        "text": "An interactive presentation authored directly in Oviz",
                        "x": 112,
                        "y": 720,
                        "width": 992,
                        "height": 82,
                        "font_size": 29,
                        "font_weight": 400,
                        "color": "#e6eef8",
                        "align": "left",
                    },
                ],
                "notes": "Opening title slide.",
            },
            {
                "id": "deck-lookback",
                "name": "Look back in time",
                "state_id": "deck-lookback-five-myr",
                "objects": [
                    {
                        "id": "deck-lookback-title",
                        "kind": "text",
                        "type": "title",
                        "text": "Following the clusters backward",
                        "x": 112,
                        "y": 90,
                        "width": 992,
                        "height": 112,
                        "font_size": 54,
                        "font_weight": 700,
                        "color": "#ffffff",
                        "align": "left",
                    },
                    {
                        "id": "deck-lookback-text",
                        "kind": "text",
                        "type": "text",
                        "text": "The saved Oviz State controls the camera, traces, selections, and time while Reveal.js controls the slide.",
                        "x": 112,
                        "y": 198,
                        "width": 800,
                        "height": 132,
                        "font_size": 25,
                        "font_weight": 400,
                        "color": "#e6eef8",
                        "align": "left",
                    },
                ],
                "notes": "Demonstrates a slide linked to a nonzero-time State.",
            },
        ],
    }
    return scene


def add_reveal_preload(html: str) -> str:
    preload = f"""
    <link data-oviz-deck-reveal="true" rel="stylesheet" href="https://cdn.jsdelivr.net/npm/reveal.js@{REVEAL_VERSION}/dist/reveal.css" />
    <link rel="preload" as="script" href="https://cdn.jsdelivr.net/npm/reveal.js@{REVEAL_VERSION}/dist/reveal.js" />
"""
    closing_head = "  </head>" if "  </head>" in html else "</head>"
    return html.replace(closing_head, preload + closing_head, 1)


def build_authoring_figure(source_html: Path, output_html: Path) -> Path:
    source_html = Path(source_html).expanduser().resolve()
    output_html = Path(output_html).expanduser().resolve()
    if source_html == output_html:
        raise ValueError("The Reveal-authoring figure must be a copy, not the canonical source.")
    scene = build_authoring_scene(read_embedded_scene_spec(source_html))
    html = ThreeJSFigure(scene, compress_scene_spec=True).to_html(compress_scene_spec=True)
    output_html.parent.mkdir(parents=True, exist_ok=True)
    output_html.write_text(add_reveal_preload(html), encoding="utf-8")
    return output_html


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-html", type=Path, default=SOURCE_HTML)
    parser.add_argument("--output-html", type=Path, default=DEFAULT_OUTPUT_HTML)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    output = build_authoring_figure(args.source_html, args.output_html)
    print(f"Wrote {output}")


if __name__ == "__main__":
    main()
