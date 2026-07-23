import importlib.util
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest


SCRIPT_PATH = Path(__file__).with_name("main_figure_chronos_july4_reveal_authoring.py")
ARTIFACT_PATH = SCRIPT_PATH.with_suffix(".html")


def _load_module():
    spec = importlib.util.spec_from_file_location("july_reveal_authoring", SCRIPT_PATH)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_build_authoring_scene_adds_linked_slides_without_mutating_source():
    module = _load_module()
    source = {
        "title": "Canonical",
        "initial_state": {"current_group": "Clusters"},
        "states": {"project_id": "project", "items": []},
    }

    scene = module.build_authoring_scene(source)

    assert "deck" not in source
    assert scene["title"] == "Young Clusters — Reveal Authoring"
    assert [slide["name"] for slide in scene["deck"]["slides"]] == [
        "Opening",
        "Look back in time",
    ]
    assert scene["deck"]["slides"][0]["state_id"] is None
    assert scene["deck"]["slides"][1]["state_id"] == "deck-lookback-five-myr"
    assert scene["deck"]["schema_version"] == 2
    assert scene["deck"]["slides"][0]["objects"][0]["kind"] == "text"
    assert scene["states"]["items"][0]["snapshot"]["current_frame_value"] == 115.0


def test_build_authoring_figure_refuses_to_overwrite_source(tmp_path):
    module = _load_module()
    source = tmp_path / "scene.html"
    source.write_text("<html></html>", encoding="utf-8")

    with pytest.raises(ValueError, match="must be a copy"):
        module.build_authoring_figure(source, source)


def test_add_reveal_preload_pins_reveal_assets():
    module = _load_module()
    html = module.add_reveal_preload("<html><head></head><body></body></html>")

    assert f"reveal.js@{module.REVEAL_VERSION}/dist/reveal.css" in html
    assert f'as="script" href="https://cdn.jsdelivr.net/npm/reveal.js@{module.REVEAL_VERSION}/dist/reveal.js"' in html
    assert html.index("reveal.css") < html.index("</head>")


@pytest.mark.skipif(not ARTIFACT_PATH.exists(), reason="Reveal-authoring artifact has not been generated")
def test_reveal_authoring_artifact_embeds_editable_deck_and_compact_scene(tmp_path):
    module = _load_module()
    html = ARTIFACT_PATH.read_text(encoding="utf-8")
    scene = module.read_embedded_scene_spec(ARTIFACT_PATH)

    assert ARTIFACT_PATH.stat().st_size < 100 * 1024 * 1024
    assert scene["deck"]["enabled"] is True
    assert scene["deck"]["embedded"] is True
    assert scene["deck"]["schema_version"] == 2
    assert len(scene["deck"]["slides"]) == 2
    assert len(scene["deck"]["slides"][0]["objects"]) == 3
    assert scene["deck"]["slides"][1]["state_id"] == "deck-lookback-five-myr"
    assert scene["states"]["items"][0]["id"] == "deck-lookback-five-myr"
    assert html.count("oviz-deck-authoring-layer") > 1
    assert "function initializeOvizDeck()" in html
    assert "new RevealConstructor(ovizDeckRevealRootEl" in html
    assert "viewer.deck = {" in html
    assert "deck-present-start" in html

    if shutil.which("node"):
        scripts = re.findall(r"<script>(.*?)</script>", html, re.S)
        viewer_script = next(script for script in scripts if "/*__SCENE_SPEC_START__*/" in script)
        script_path = tmp_path / "reveal-authoring.js"
        script_path.write_text(viewer_script, encoding="utf-8")
        subprocess.run(["node", "--check", str(script_path)], check=True)
