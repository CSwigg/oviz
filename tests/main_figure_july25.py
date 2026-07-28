#!/usr/bin/env python3
"""Regenerate the July 25 non-paper Oviz figure with current cluster velocities.

The source scene is rebuilt with the July 4 Chronos ages and a selectable
cluster-velocity CSV, then converted to the clean July presentation while
keeping the editable Slides authoring option available. The embedded Three.js
runtime carries the volumetric quality and FPS upgrades developed for the
July 25 artifact.
"""

from __future__ import annotations

import argparse
import base64
import gzip
import json
import re
import sys
import tempfile
from copy import deepcopy
from pathlib import Path

from oviz.threejs_figure import ThreeJSFigure


TESTS_DIR = Path(__file__).resolve().parent
if str(TESTS_DIR) not in sys.path:
    sys.path.insert(0, str(TESTS_DIR))

from main_figure_chronos_july4 import (  # noqa: E402
    CLUSTER_MEMBERS_PATH,
    JULY4_CHRONOS_RESULTS_PATH,
)
from main_figure_new_chronos import run_main_figure as build_source_main_figure  # noqa: E402


SOURCE_HTML = Path(__file__).with_name("main_figure_chronos_july4.html")
DEFAULT_OUTPUT_HTML = Path(__file__).with_suffix(".html")
DEFAULT_CLUSTER_VELOCITIES_PATH = Path(
    "/Users/swiggumc/Desktop/astro_research/cfa/velocity analysis/outputs/release/"
    "cluster_velocities_sdssv_covariance_audited.csv"
)
EDENHOFER_GALACTIC_LIGHTING = {
    "lighting_mode": "standard",
    "galactic_center": [8122.0, 0.0, 0.0],
    "galactic_light_intensity": 1.35,
    "galactic_ambient": 0.22,
    "galactic_extinction": 2.4,
    "galactic_scattering": 0.55,
    "galactic_anisotropy": 0.45,
    "galactic_warmth": 0.72,
}


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


def build_state_only_scene(
    scene_spec: dict,
    cluster_velocities_path: Path | None = None,
) -> dict:
    scene = deepcopy(scene_spec)
    scene["title"] = ""
    scene["deck"] = {
        "schema_version": 2,
        "available": True,
        "enabled": False,
        "embedded": False,
        "revision": 0,
        "aspect_ratio": "16:9",
        "guides": {"smart": True, "grid": False, "grid_size": 20},
        "slides": [],
    }
    initial_volume_state = (
        scene.setdefault("initial_state", {}).setdefault("volume_state_by_key", {})
    )
    scene["initial_state"].setdefault("global_controls", {}).update({
        "size_points_by_stars_enabled": True,
        "fade_opacity_by_birth_time_enabled": True,
    })
    for layer in (scene.get("volumes") or {}).get("layers") or []:
        if "edenhofer" not in str(layer.get("name") or "").lower():
            continue
        defaults = layer.setdefault("default_controls", {})
        defaults.update(EDENHOFER_GALACTIC_LIGHTING)
        state_key = str(layer.get("state_key") or layer.get("key") or "")
        if state_key:
            state = initial_volume_state.setdefault(state_key, {})
            state.update({
                "lightingMode": "standard",
                "galacticCenter": [8122.0, 0.0, 0.0],
                "galacticLightIntensity": 1.35,
                "galacticAmbient": 0.22,
                "galacticExtinction": 2.4,
                "galacticScattering": 0.55,
                "galacticAnisotropy": 0.45,
                "galacticWarmth": 0.72,
            })
    if cluster_velocities_path is not None:
        velocity_path = Path(cluster_velocities_path).expanduser().resolve()
        scene.setdefault("provenance", {})["cluster_velocities"] = {
            "filename": velocity_path.name,
            "path": str(velocity_path),
            "release": "SDSS-V covariance audited",
        }
    return scene


def build_figure(
    source_html: Path,
    output_html: Path,
    cluster_velocities_path: Path | None = None,
) -> Path:
    source_html = Path(source_html).expanduser().resolve()
    output_html = Path(output_html).expanduser().resolve()
    if source_html == output_html:
        raise ValueError("The July 25 figure must be a sibling copy, not the source artifact.")
    scene = build_state_only_scene(
        read_embedded_scene_spec(source_html),
        cluster_velocities_path=cluster_velocities_path,
    )
    html = ThreeJSFigure(scene, compress_scene_spec=True).to_html(compress_scene_spec=True)
    output_html.parent.mkdir(parents=True, exist_ok=True)
    output_html.write_text(html, encoding="utf-8")
    return output_html


def build_velocity_source_figure(
    cluster_velocities_path: Path,
    output_html: Path,
) -> Path:
    cluster_velocities_path = Path(cluster_velocities_path).expanduser().resolve()
    if not cluster_velocities_path.exists():
        raise FileNotFoundError(
            f"Missing cluster velocity catalog: {cluster_velocities_path}"
        )
    return build_source_main_figure(
        output_html=output_html,
        mobile_mode=False,
        compact_payload=True,
        mobile_safe_mode=True,
        chronos_results_path=JULY4_CHRONOS_RESULTS_PATH,
        chronos_model="parsec",
        include_spiral_arms=False,
        jun6_catalog=True,
        cluster_velocities_path=cluster_velocities_path,
        include_background_cluster_trace=False,
        cluster_members_file=CLUSTER_MEMBERS_PATH,
        show_cluster_members_in_sky=True,
        website_output_html=None,
    )


def build_figure_from_velocity_catalog(
    cluster_velocities_path: Path,
    output_html: Path,
) -> Path:
    with tempfile.TemporaryDirectory(prefix="oviz-july25-source-") as temp_dir:
        source_html = Path(temp_dir) / "main_figure_july25_source.html"
        build_velocity_source_figure(cluster_velocities_path, source_html)
        return build_figure(
            source_html,
            output_html,
            cluster_velocities_path=cluster_velocities_path,
        )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source-html",
        type=Path,
        default=None,
        help=(
            "Optional prebuilt source scene. If omitted, rebuild the source "
            "from the selected cluster velocity catalog."
        ),
    )
    parser.add_argument(
        "--cluster-velocities-path",
        type=Path,
        default=DEFAULT_CLUSTER_VELOCITIES_PATH,
        help="Cluster velocity CSV used when rebuilding the source scene.",
    )
    parser.add_argument("--output-html", type=Path, default=DEFAULT_OUTPUT_HTML)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.source_html is not None:
        output = build_figure(
            args.source_html,
            args.output_html,
            cluster_velocities_path=args.cluster_velocities_path,
        )
    else:
        output = build_figure_from_velocity_catalog(
            args.cluster_velocities_path,
            args.output_html,
        )
    print(f"Wrote {output}")


if __name__ == "__main__":
    main()
