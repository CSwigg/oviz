import importlib.util
import re
import shutil
import subprocess
import sys
from pathlib import Path


SCRIPT_PATH = Path(__file__).with_name("main_figure_july21.py")
ARTIFACT_PATH = SCRIPT_PATH.with_suffix(".html")


def _load_module():
    spec = importlib.util.spec_from_file_location("main_figure_july21", SCRIPT_PATH)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_build_state_only_scene_disables_and_removes_slides():
    module = _load_module()
    source = {"title": "Source", "deck": {"available": True, "slides": [{"id": "old"}]}}

    scene = module.build_state_only_scene(source)

    assert source["deck"]["slides"] == [{"id": "old"}]
    assert scene["deck"]["available"] is False
    assert scene["deck"]["enabled"] is False
    assert scene["deck"]["slides"] == []


def test_july21_artifact_uses_state_only_presentation(tmp_path):
    module = _load_module()
    html = ARTIFACT_PATH.read_text(encoding="utf-8")
    scene = module.read_embedded_scene_spec(ARTIFACT_PATH)

    assert ARTIFACT_PATH.stat().st_size < 100 * 1024 * 1024
    assert scene["deck"]["available"] is False
    assert scene["deck"]["enabled"] is False
    assert scene["deck"]["slides"] == []
    assert scene["title"] == ""
    assert "Young Clusters — July 21" not in html
    assert ">Slides ▸</button>" not in html
    assert 'class="oviz-three-deck-editor' not in html
    assert ") ? ovizDeckPrevious() : ovizStatesPresentationPrevious();" in html
    assert ") ? ovizDeckNext() : ovizStatesPresentationNext();" in html
    assert "function ovizStatesPresentationNext()" in html
    assert "function ovizStatesPresentationPrevious()" in html
    assert "function ovizQueuePresentationStateNavigation(direction)" in html
    assert "await ovizWaitForStatesControllerReady();" in html
    assert "activeTransition.promise" in html
    assert "ovizResetPresentationStateNavigationQueue(" in html

    if shutil.which("node"):
        scripts = re.findall(r"<script>(.*?)</script>", html, re.S)
        viewer_script = next(script for script in scripts if "/*__SCENE_SPEC_START__*/" in script)
        script_path = tmp_path / "main-figure-july21.js"
        script_path.write_text(viewer_script, encoding="utf-8")
        subprocess.run(["node", "--check", str(script_path)], check=True)
