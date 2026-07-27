#!/usr/bin/env python3
"""Build the paper-mode demo: arXiv 2406.06510 as a read-along Oviz figure.

The scene is the July figure (cluster catalog + Edenhofer dust). The paper
panel on the right carries Swiggum et al. 2024 ("Most nearby young star
clusters formed in three massive complexes"); scroll anchors drive States
that fly the camera through the story: present-day overview, orbit
traceback, the three-family convergence, each family up close, the dust
context, and the all-sky view.
"""

from __future__ import annotations

import argparse
import base64
import gzip
import json
import re
from copy import deepcopy
from pathlib import Path

from oviz.paper import Paper
from oviz.threejs_figure import ThreeJSFigure


SOURCE_HTML = Path(__file__).with_name("main_figure_chronos_july4.html")
DEFAULT_OUTPUT_HTML = Path(__file__).with_suffix(".html")
CONTENT_DIR = Path(__file__).with_name("paper_2406_06510")
PANEL_WIDTH_FRACTION = 0.42
# Push scene subjects toward the visible left region while the panel is open.
VIEW_OFFSET = {"x": 0.21, "y": 0.0}

# Traces from newer chronology work that play no role in arXiv 2406.06510;
# the paper demo keeps only the Sun, the two age-sliced cluster samples, the
# three families, the grid helpers, and the Edenhofer dust volume.
PAPER_IRRELEVANT_TRACES = {
    "Full Cluster Catalog",
    "Lacerta Family",
    "Trumpler 3 Family",
    "Proto Orion Family",
    "Radcliffe Wave Clusters",
    "Split Selection Clusters",
    "Cepheus Spur Clusters",
    "Vela-Sagittarius Clusters",
}


def trace_is_paper_irrelevant(name: object) -> bool:
    base = str(name or "")
    if base in PAPER_IRRELEVANT_TRACES:
        return True
    return any(base == f"{removed} Track" for removed in PAPER_IRRELEVANT_TRACES)


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


def frame_index_for_time(frames: list[dict], time_myr: float) -> int:
    best_index = 0
    best_delta = float("inf")
    for index, frame in enumerate(frames):
        delta = abs(float(frame.get("time", 0.0)) - time_myr)
        if delta < best_delta:
            best_delta = delta
            best_index = index
    return best_index


def trace_centroid(frame: dict, trace_name: str) -> tuple[float, float, float] | None:
    for trace in frame.get("traces", []):
        if str(trace.get("name")) == trace_name:
            points = trace.get("points") or []
            if not points:
                return None
            count = float(len(points))
            return (
                sum(float(p["x"]) for p in points) / count,
                sum(float(p["y"]) for p in points) / count,
                sum(float(p["z"]) for p in points) / count,
            )
    return None


def build_states(scene: dict) -> dict:
    frames = scene["frames"]

    def camera(position, target=(0.0, 0.0, 0.0)):
        return {
            "position": {"x": position[0], "y": position[1], "z": position[2]},
            "target": {"x": target[0], "y": target[1], "z": target[2]},
            "up": {"x": 0.0, "y": 0.0, "z": 1.0},
            "view_offset": dict(VIEW_OFFSET),
        }

    def snapshot(*, time_myr, cam, group="Clusters", globals_extra=None):
        # Snapshots stay minimal and only carry keys that
        # captureRuntimeState() also emits: the post-transition fidelity
        # check iterates the target's keys, so any build-pipeline flag
        # copied from initial_state would fail it and abort the transition.
        frame_index = frame_index_for_time(frames, time_myr)
        controls = {"camera_view_mode": "free"}
        if globals_extra:
            controls.update(globals_extra)
        return {
            "current_group": group,
            "current_frame_index": frame_index,
            "current_frame_value": float(frame_index),
            "playback_state": {"direction": 0, "interval_ms": 1000},
            "camera": cam,
            "zen_mode_enabled": False,
            "global_controls": controls,
        }

    def family_camera(frame_time, family_name, back=420.0, lift=230.0):
        frame = frames[frame_index_for_time(frames, frame_time)]
        centroid = trace_centroid(frame, family_name) or (0.0, 0.0, 0.0)
        cx, cy, cz = centroid
        length = max((cx * cx + cy * cy) ** 0.5, 1.0)
        away = (cx / length, cy / length)
        position = (cx + away[0] * back, cy + away[1] * back, cz + lift)
        return camera(position, centroid)

    overview_camera = camera((60.0, -1750.0, 1080.0))
    dust_camera = camera((720.0, -1150.0, 330.0), (-120.0, 60.0, 30.0))
    converge_camera = camera((160.0, -900.0, 620.0), (10.0, -30.0, 20.0))

    items = [
        {
            "id": "paper-overview",
            "name": "Overview today",
            "transition": {"duration_ms": 1600, "easing": "easeInOutCubic"},
            "snapshot": snapshot(time_myr=0.0, cam=overview_camera, group="Clusters"),
        },
        {
            "id": "paper-traceback",
            "name": "Orbit traceback",
            "transition": {"duration_ms": 2400, "easing": "easeInOutCubic"},
            "snapshot": snapshot(time_myr=-30.0, cam=overview_camera, group="Clusters"),
        },
        {
            "id": "paper-converge",
            "name": "Three families converge",
            "transition": {"duration_ms": 2200, "easing": "easeInOutCubic"},
            "snapshot": snapshot(
                time_myr=-40.0, cam=converge_camera, group="Cluster Families"
            ),
        },
        {
            "id": "paper-cr135",
            "name": "Collinder 135 family",
            "transition": {"duration_ms": 1900, "easing": "easeInOutCubic"},
            "snapshot": snapshot(
                time_myr=-35.0,
                cam=family_camera(-35.0, "Cr 135 Family"),
                group="Cluster Families",
            ),
        },
        {
            "id": "paper-m6",
            "name": "Messier 6 family",
            "transition": {"duration_ms": 1900, "easing": "easeInOutCubic"},
            "snapshot": snapshot(
                time_myr=-35.0,
                cam=family_camera(-35.0, "M6 Family"),
                group="Cluster Families",
            ),
        },
        {
            "id": "paper-alphaper",
            "name": "Alpha Persei family",
            "transition": {"duration_ms": 1900, "easing": "easeInOutCubic"},
            "snapshot": snapshot(
                time_myr=-35.0,
                cam=family_camera(-35.0, "Alpha Persei Family"),
                group="Cluster Families",
            ),
        },
        {
            "id": "paper-dust",
            "name": "Dust context",
            "transition": {"duration_ms": 2000, "easing": "easeInOutCubic"},
            "snapshot": snapshot(time_myr=0.0, cam=dust_camera, group="Dust Structures"),
        },
        {
            "id": "paper-sky",
            "name": "All-sky view",
            "transition": {"duration_ms": 2400, "easing": "easeInOutCubic"},
            "snapshot": snapshot(
                time_myr=0.0,
                cam={
                    "position": {"x": 0.004, "y": -0.007, "z": 0.0},
                    "target": {"x": 0.0, "y": 0.0, "z": 0.0},
                    "up": {"x": 0.0, "y": 0.0, "z": 1.0},
                    "view_offset": dict(VIEW_OFFSET),
                },
                group="Clusters",
                globals_extra={
                    "camera_view_mode": "earth",
                    "sky_member_display_mode": "stars",
                },
            ),
        },
        {
            "id": "paper-return",
            "name": "Overview return",
            "transition": {"duration_ms": 2000, "easing": "easeInOutCubic"},
            "snapshot": snapshot(time_myr=0.0, cam=overview_camera, group="Clusters"),
        },
    ]
    return {
        "schema_version": 1,
        "project_id": "paper-2406-06510-demo",
        "revision": 1,
        "synchronized_revision": 1,
        "default_mode": "edit",
        "default_transition": {"duration_ms": 1600, "easing": "easeInOutCubic"},
        "items": items,
        "assets": {},
    }


PARAGRAPH_ANCHORS = [
    # (keyword regex over plain text, state name, label)
    (r"seventy years|past seventy", "Overview today", "Present-day cluster census"),
    (r"integrate the clusters|Going back 30", "Orbit traceback", "Orbits traced back in time"),
    (r"HDBSCAN|membership list", "Three families converge", "Families identified at birth"),
    (r"supernovae|progenitors", "Dust context", "Feedback shaped the local ISM"),
    (r"Collinder 135|Cr135 family|Cr 135 family", "Collinder 135 family", "The Cr 135 family"),
    (r"M6 family|Messier 6", "Messier 6 family", "The M 6 family"),
    (r"Alpha Persei .{0,10}family|𝛼Per family", "Alpha Persei family", "The α Per family"),
    (r"Local Bubble|GSH ?238", "Dust context", "The Local Bubble and GSH 238"),
]

FIGURE_PLACEMENTS = {
    # figure label -> (after paragraph matching regex, state name, live)
    "Figure 1": (r"HDBSCAN|membership list", "Three families converge", True),
    "Table 1": (r"supernovae|progenitors", None, False),
    "Figure 2": (r"Local Bubble|GSH ?238", "Dust context", True),
    "Figure 3": (r"Alpha Persei .{0,10}family|𝛼Per family", "All-sky view", True),
}


def plain_text(html: str) -> str:
    import html as html_mod

    return html_mod.unescape(re.sub(r"<[^>]+>", "", html))


def build_paper(states: dict) -> dict:
    content = json.loads((CONTENT_DIR / "content.json").read_text(encoding="utf-8"))
    figures_by_label = {fig["label"]: fig for fig in content["figures"]}

    paper = Paper(
        content["title"],
        authors=[plain_text(content["authors_html"])],
        venue_html=(
            'Published in <i>Nature</i> · '
            f'<a href="{content["arxiv_url"]}" target="_blank" rel="noopener">arXiv:2406.06510</a>'
        ),
        link_url=content["arxiv_url"],
        panel_width_fraction=PANEL_WIDTH_FRACTION,
        reading_line=0.35,
        math_enabled=False,
        enabled=True,
    )

    paper.add_section("Abstract", level=1)
    paper.add_html(
        f"<p>{content['abstract_html']}</p>",
        state="Overview today",
        label="Present-day overview",
        transition_ms=1600,
    )

    paper.add_section("Main", level=1)
    used_anchors: set[str] = set()
    pending_figures = dict(FIGURE_PLACEMENTS)
    for paragraph in content["main_paragraphs"]:
        text = plain_text(paragraph)
        anchor_state = None
        anchor_label = None
        for pattern, state_name, label in PARAGRAPH_ANCHORS:
            if state_name in used_anchors:
                continue
            if re.search(pattern, text, re.IGNORECASE):
                anchor_state = state_name
                anchor_label = label
                used_anchors.add(state_name)
                break
        paper.add_html(
            f"<p>{paragraph}</p>",
            state=anchor_state,
            label=anchor_label,
        )
        for figure_label, (pattern, figure_state, live) in list(pending_figures.items()):
            if re.search(pattern, text, re.IGNORECASE):
                figure = figures_by_label.get(figure_label)
                if figure and figure.get("file"):
                    paper.add_figure(
                        CONTENT_DIR / figure["file"],
                        caption_html=figure["caption_html"],
                        max_width_px=1400,
                        jpeg_quality=82,
                        live=bool(live),
                        state=figure_state,
                        label=f"{figure_label} view" if figure_state else None,
                    )
                pending_figures.pop(figure_label)

    paper.add_section("Methods", level=1)
    for block in content["methods_blocks"]:
        if "heading" in block:
            paper.add_html(f"<h4>{block['heading']}</h4>")
        else:
            paper.add_html(f"<p>{block['paragraph']}</p>")
    for label in ("Extended Data Figure 1", "Extended Data Figure 2"):
        figure = figures_by_label.get(label)
        if figure and figure.get("file"):
            paper.add_figure(
                CONTENT_DIR / figure["file"],
                caption_html=figure["caption_html"],
                max_width_px=1200,
                jpeg_quality=82,
            )

    paper.add_section("References", level=1)
    items = "".join(f"<li>{reference}</li>" for reference in content["references"])
    paper.add_html(f'<ol class="oviz-paper-references">{items}</ol>')

    final = paper.add_section("", level=1)
    paper.add_html(
        "<p><i>End of paper — scroll targets return the scene to the overview.</i></p>",
        state="Overview return",
        label="Back to the overview",
    )
    del final
    return paper.to_spec(states=states)


def strip_paper_irrelevant_traces(scene: dict) -> None:
    removed_keys: set[str] = set()
    for frame in scene.get("frames") or []:
        traces = frame.get("traces") or []
        for trace in traces:
            if trace_is_paper_irrelevant(trace.get("name")) and trace.get("key"):
                removed_keys.add(str(trace["key"]))
        frame["traces"] = [
            trace for trace in traces if not trace_is_paper_irrelevant(trace.get("name"))
        ]
    for visibility in (scene.get("group_visibility") or {}).values():
        if isinstance(visibility, dict):
            for key in removed_keys:
                visibility.pop(key, None)
    legend = scene.get("legend") or {}
    if isinstance(legend.get("items"), list):
        legend["items"] = [
            item for item in legend["items"]
            if not trace_is_paper_irrelevant(item.get("name"))
            and str(item.get("key")) not in removed_keys
        ]
    animation = scene.get("animation") or {}
    if isinstance(animation.get("focus_options"), list):
        animation["focus_options"] = [
            option for option in animation["focus_options"]
            if not trace_is_paper_irrelevant(option.get("name"))
            and str(option.get("key")) not in removed_keys
        ]
    members = (scene.get("sky_panel") or {}).get("members_by_cluster")
    if isinstance(members, dict):
        for name in list(members.keys()):
            if trace_is_paper_irrelevant(name):
                members.pop(name, None)


def build_state_only_scene(scene_spec: dict) -> dict:
    scene = deepcopy(scene_spec)
    strip_paper_irrelevant_traces(scene)
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
    initial_volume_state = (
        scene.setdefault("initial_state", {}).setdefault("volume_state_by_key", {})
    )
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
    return scene


def build_figure(source_html: Path, output_html: Path) -> Path:
    source_html = Path(source_html).expanduser().resolve()
    output_html = Path(output_html).expanduser().resolve()
    if source_html == output_html:
        raise ValueError("The paper demo must be a sibling copy, not the source artifact.")
    scene = build_state_only_scene(read_embedded_scene_spec(source_html))
    scene["states"] = build_states(scene)
    scene["paper"] = build_paper(scene["states"])
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
