#!/usr/bin/env python3
"""Build the July 29 Oviz figure with clusters younger than 120 Myr.

The source scene comes from ``main_figure_july29_scene.py`` — a plain,
checked-in Python script (materialized once from the old authoring
notebook, now the authoritative source). Executing it yields the
``ThreeJSFigure`` in memory; the July 25 presentation cleanup then runs on
the scene dict directly. No notebook and no HTML artifact are involved.
"""

from __future__ import annotations

import argparse
import os
import re
import sys
from pathlib import Path

TESTS_DIR = Path(__file__).resolve().parent
REPO_ROOT = TESTS_DIR.parent
for _path in (str(REPO_ROOT), str(TESTS_DIR)):
    if _path not in sys.path:
        sys.path.insert(0, _path)

from oviz.threejs_embed import scene_spec_jsonable  # noqa: E402
from oviz.threejs_figure import ThreeJSFigure  # noqa: E402

import main_figure_july25 as july25  # noqa: E402
import main_figure_new_chronos as source_runner  # noqa: E402


DEFAULT_OUTPUT_HTML = Path(__file__).with_suffix(".html")
SCENE_SCRIPT_PATH = Path(__file__).with_name("main_figure_july29_scene.py")
SCENE_SCRIPT_VELOCITIES_MARKER = "SCENE_SCRIPT_CLUSTER_VELOCITIES"
JULY29_LOOKBACK_MYR = 120
JULY29_CLUSTER_MAX_AGE_MYR = 120
SOURCE_CLUSTER_LABEL = "Clusters (< 60 Myr)"
JULY29_CLUSTER_LABEL = "Clusters (< 120 Myr)"


def _replace_cluster_label(value):
    """Replace the primary sample label everywhere without changing stable keys."""
    if isinstance(value, dict):
        return {
            key: _replace_cluster_label(item)
            for key, item in value.items()
        }
    if isinstance(value, list):
        return [_replace_cluster_label(item) for item in value]
    if isinstance(value, tuple):
        return tuple(_replace_cluster_label(item) for item in value)
    if value == SOURCE_CLUSTER_LABEL:
        return JULY29_CLUSTER_LABEL
    return value


def regenerate_scene_script(
    cluster_velocities_path: Path = july25.DEFAULT_CLUSTER_VELOCITIES_PATH,
    ratzenboeck_scocen_catalog_path: Path = july25.DEFAULT_RATZENBOECK_SCOCEN_CATALOG_PATH,
    ratzenboeck_scocen_metadata_path: Path = july25.DEFAULT_RATZENBOECK_SCOCEN_METADATA_PATH,
    ratzenboeck_scocen_members_path: Path = july25.DEFAULT_RATZENBOECK_SCOCEN_MEMBERS_PATH,
) -> Path:
    """Materialize the July 29 source-scene build as a plain Python script.

    One-time escape hatch from the old notebook pipeline: converts and
    patches the notebook exactly as the runner used to, nulls ``save_name``
    so nothing is written, and saves the result as the checked-in scene
    script. After that the script itself is the source of truth — edit it
    directly rather than regenerating.
    """
    (
        cluster_velocities_path,
        ratzenboeck_scocen_catalog_path,
        ratzenboeck_scocen_metadata_path,
        _members_path,
    ) = july25.resolve_source_input_paths(
        cluster_velocities_path,
        ratzenboeck_scocen_catalog_path,
        ratzenboeck_scocen_metadata_path,
        ratzenboeck_scocen_members_path,
    )
    previous_age_limit = source_runner.MAIN_FIGURE_CLUSTER_BLUE_MAX_AGE_MYR
    source_runner.MAIN_FIGURE_CLUSTER_BLUE_MAX_AGE_MYR = (
        JULY29_CLUSTER_MAX_AGE_MYR
    )
    try:
        source = source_runner.build_patched_main_figure_source(
            output_html=SCENE_SCRIPT_PATH.with_suffix(".unused.html"),
            **july25.source_main_figure_kwargs(
                cluster_velocities_path,
                ratzenboeck_scocen_catalog_path,
                ratzenboeck_scocen_metadata_path,
                lookback_myr=JULY29_LOOKBACK_MYR,
            ),
        )
    finally:
        source_runner.MAIN_FIGURE_CLUSTER_BLUE_MAX_AGE_MYR = previous_age_limit
    source, nulled_save_name = re.subn(
        r"(?m)^save_name\s*=\s*['\"][^'\"]+['\"]",
        "save_name = None",
        source,
        count=1,
    )
    if nulled_save_name != 1:
        raise RuntimeError("Could not null save_name in the scene script.")
    header = (
        '"""July 29 source-scene script (plain Python, no notebook).\n\n'
        "Running this module builds the July 29 source ThreeJSFigure as\n"
        "``fig3d`` without writing any file. It was materialized from the\n"
        "retired notebook pipeline and is now the authoritative source —\n"
        "edit it directly.\n"
        '"""\n\n'
        f"{SCENE_SCRIPT_VELOCITIES_MARKER} = "
        f"{str(cluster_velocities_path)!r}\n\n"
    )
    SCENE_SCRIPT_PATH.write_text(header + source, encoding="utf-8")
    return SCENE_SCRIPT_PATH


def _scene_script_velocities(script_text: str) -> str | None:
    match = re.search(
        SCENE_SCRIPT_VELOCITIES_MARKER + r"\s*=\s*'([^']*)'",
        script_text,
    )
    return match.group(1) if match else None


def build_source_scene(
    cluster_velocities_path: Path,
    ratzenboeck_scocen_catalog_path: Path = july25.DEFAULT_RATZENBOECK_SCOCEN_CATALOG_PATH,
    ratzenboeck_scocen_metadata_path: Path = july25.DEFAULT_RATZENBOECK_SCOCEN_METADATA_PATH,
    ratzenboeck_scocen_members_path: Path = july25.DEFAULT_RATZENBOECK_SCOCEN_MEMBERS_PATH,
) -> dict:
    """Execute the scene script and return the source scene as a dict."""
    cluster_velocities_path = Path(cluster_velocities_path).expanduser().resolve()
    script_is_stale = (
        not SCENE_SCRIPT_PATH.exists()
        or _scene_script_velocities(
            SCENE_SCRIPT_PATH.read_text(encoding="utf-8")
        ) != str(cluster_velocities_path)
    )
    if script_is_stale:
        regenerate_scene_script(
            cluster_velocities_path=cluster_velocities_path,
            ratzenboeck_scocen_catalog_path=ratzenboeck_scocen_catalog_path,
            ratzenboeck_scocen_metadata_path=ratzenboeck_scocen_metadata_path,
            ratzenboeck_scocen_members_path=ratzenboeck_scocen_members_path,
        )
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl")
    os.environ.setdefault("XDG_CACHE_HOME", "/tmp")
    os.environ.setdefault("MPLBACKEND", "Agg")
    namespace = {
        "__name__": "__main__",
        "__file__": str(SCENE_SCRIPT_PATH),
    }
    exec(
        compile(
            SCENE_SCRIPT_PATH.read_text(encoding="utf-8"),
            str(SCENE_SCRIPT_PATH),
            "exec",
        ),
        namespace,
    )
    figure = namespace.get("fig3d")
    if figure is None or not hasattr(figure, "to_dict"):
        raise RuntimeError(
            "The scene script did not produce a ThreeJSFigure named fig3d."
        )
    return scene_spec_jsonable(figure.to_dict())


def build_scene(
    cluster_velocities_path: Path = july25.DEFAULT_CLUSTER_VELOCITIES_PATH,
    ratzenboeck_scocen_catalog_path: Path = july25.DEFAULT_RATZENBOECK_SCOCEN_CATALOG_PATH,
    ratzenboeck_scocen_metadata_path: Path = july25.DEFAULT_RATZENBOECK_SCOCEN_METADATA_PATH,
    ratzenboeck_scocen_members_path: Path = july25.DEFAULT_RATZENBOECK_SCOCEN_MEMBERS_PATH,
) -> dict:
    """Return the finished July 29 presentation scene as a Python dict."""
    source_scene = build_source_scene(
        cluster_velocities_path,
        ratzenboeck_scocen_catalog_path=ratzenboeck_scocen_catalog_path,
        ratzenboeck_scocen_metadata_path=ratzenboeck_scocen_metadata_path,
        ratzenboeck_scocen_members_path=ratzenboeck_scocen_members_path,
    )
    scene = july25.build_state_only_scene(
        source_scene,
        cluster_velocities_path=cluster_velocities_path,
        ratzenboeck_scocen_catalog_path=ratzenboeck_scocen_catalog_path,
        ratzenboeck_scocen_metadata_path=ratzenboeck_scocen_metadata_path,
        ratzenboeck_scocen_members_path=ratzenboeck_scocen_members_path,
        # The July 29 artifact keeps the full catalog: downstream consumers
        # (the ccSN reconstruction figure) select their sample from trace-1.
        drop_full_cluster_catalog=False,
    )
    return _replace_cluster_label(scene)


def build_figure_from_velocity_catalog(
    cluster_velocities_path: Path,
    output_html: Path,
    *,
    ratzenboeck_scocen_catalog_path: Path = july25.DEFAULT_RATZENBOECK_SCOCEN_CATALOG_PATH,
    ratzenboeck_scocen_metadata_path: Path = july25.DEFAULT_RATZENBOECK_SCOCEN_METADATA_PATH,
    ratzenboeck_scocen_members_path: Path = july25.DEFAULT_RATZENBOECK_SCOCEN_MEMBERS_PATH,
) -> Path:
    output_html = Path(output_html).expanduser().resolve()
    scene = build_scene(
        cluster_velocities_path=cluster_velocities_path,
        ratzenboeck_scocen_catalog_path=ratzenboeck_scocen_catalog_path,
        ratzenboeck_scocen_metadata_path=ratzenboeck_scocen_metadata_path,
        ratzenboeck_scocen_members_path=ratzenboeck_scocen_members_path,
    )
    html = ThreeJSFigure(scene, compress_scene_spec=True).to_html(
        compress_scene_spec=True
    )
    output_html.parent.mkdir(parents=True, exist_ok=True)
    output_html.write_text(html, encoding="utf-8")
    return output_html


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--cluster-velocities-path",
        type=Path,
        default=july25.DEFAULT_CLUSTER_VELOCITIES_PATH,
    )
    parser.add_argument("--output-html", type=Path, default=DEFAULT_OUTPUT_HTML)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    output = build_figure_from_velocity_catalog(
        args.cluster_velocities_path,
        args.output_html,
    )
    print(f"Wrote {output}")


if __name__ == "__main__":
    main()
