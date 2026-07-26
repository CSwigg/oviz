import importlib.util
import re
import shutil
import subprocess
import sys
from pathlib import Path


SCRIPT_PATH = Path(__file__).with_name("main_figure_july25.py")
ARTIFACT_PATH = SCRIPT_PATH.with_suffix(".html")
JULY21_ARTIFACT_PATH = SCRIPT_PATH.with_name("main_figure_july21.html")


def _load_module():
    spec = importlib.util.spec_from_file_location("main_figure_july25", SCRIPT_PATH)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_build_state_only_scene_disables_and_removes_slides():
    module = _load_module()
    source = {
        "title": "Source",
        "deck": {"available": True, "slides": [{"id": "old"}]},
        "volumes": {
            "layers": [{
                "key": "volume-0",
                "state_key": "volume-0",
                "name": "Edenhofer+2024 Dust",
                "default_controls": {},
            }],
        },
        "initial_state": {"volume_state_by_key": {"volume-0": {"visible": True}}},
    }

    scene = module.build_state_only_scene(source)

    assert source["deck"]["slides"] == [{"id": "old"}]
    assert scene["deck"]["available"] is False
    assert scene["deck"]["enabled"] is False
    assert scene["deck"]["slides"] == []
    defaults = scene["volumes"]["layers"][0]["default_controls"]
    volume_state = scene["initial_state"]["volume_state_by_key"]["volume-0"]
    assert defaults["lighting_mode"] == "standard"
    assert defaults["galactic_center"] == [8122.0, 0.0, 0.0]
    assert volume_state["lightingMode"] == "standard"
    assert volume_state["galacticExtinction"] == 2.4


def test_july25_artifact_keeps_presentation_and_adds_runtime_upgrades(tmp_path):
    module = _load_module()
    html = ARTIFACT_PATH.read_text(encoding="utf-8")
    scene = module.read_embedded_scene_spec(ARTIFACT_PATH)

    assert ARTIFACT_PATH.stat().st_size < 100 * 1024 * 1024
    if JULY21_ARTIFACT_PATH.exists():
        assert ARTIFACT_PATH.stat().st_size <= JULY21_ARTIFACT_PATH.stat().st_size
    assert scene["deck"]["available"] is False
    assert scene["deck"]["enabled"] is False
    assert scene["deck"]["slides"] == []
    assert scene["title"] == ""
    edenhofer = next(
        layer
        for layer in scene["volumes"]["layers"]
        if "edenhofer" in str(layer.get("name") or "").lower()
    )
    edenhofer_state = scene["initial_state"]["volume_state_by_key"][edenhofer["state_key"]]
    assert edenhofer["default_controls"]["lighting_mode"] == "standard"
    assert edenhofer_state["lightingMode"] == "standard"
    assert edenhofer_state["galacticCenter"] == [8122.0, 0.0, 0.0]

    # Presentation behavior carried over from the July 21 artifact.
    assert "const experimentalVolumeLightingControlsEnabled = false;" in html
    assert ") ? ovizDeckPrevious() : ovizStatesPresentationPrevious();" in html
    assert ") ? ovizDeckNext() : ovizStatesPresentationNext();" in html
    assert "function ovizStatesPresentationNext()" in html
    assert "function ovizStatesPresentationPrevious()" in html

    # July 25 runtime upgrades: occupancy-grid empty-space skipping, IGN
    # jitter, linear LUT filtering, dithered output, capped pixel ratio.
    assert "function volumeOccupancyFor(" in html
    assert "uniform sampler3D occupancyTexture;" in html
    assert "float jitter = fract(52.9829189 * fract(dot(gl_FragCoord.xy" in html
    assert "renderer.setPixelRatio(Math.min(window.devicePixelRatio || 1, 2));" in html
    assert html.count("const VOLUME_FRAGMENT_SHADER = `") == 1

    if shutil.which("node"):
        scripts = re.findall(r"<script>(.*?)</script>", html, re.S)
        viewer_script = next(script for script in scripts if "/*__SCENE_SPEC_START__*/" in script)
        script_path = tmp_path / "main-figure-july25.js"
        script_path.write_text(viewer_script, encoding="utf-8")
        subprocess.run(["node", "--check", str(script_path)], check=True)
