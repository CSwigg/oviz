import importlib.util
import re
import shutil
import subprocess
import sys
from pathlib import Path


SCRIPT_PATH = Path(__file__).with_name("main_figure_july25.py")
ARTIFACT_PATH = SCRIPT_PATH.with_suffix(".html")


def _load_module():
    spec = importlib.util.spec_from_file_location("main_figure_july25", SCRIPT_PATH)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_build_scene_exposes_empty_slide_authoring_and_display_defaults():
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
    assert scene["deck"]["available"] is True
    assert scene["deck"]["enabled"] is False
    assert scene["deck"]["slides"] == []
    assert scene["initial_state"]["global_controls"] == {
        "size_points_by_stars_enabled": True,
        "fade_opacity_by_birth_time_enabled": True,
    }
    defaults = scene["volumes"]["layers"][0]["default_controls"]
    volume_state = scene["initial_state"]["volume_state_by_key"]["volume-0"]
    assert defaults["lighting_mode"] == "standard"
    assert defaults["galactic_center"] == [8122.0, 0.0, 0.0]
    assert volume_state["lightingMode"] == "standard"
    assert volume_state["galacticExtinction"] == 2.4


def test_default_velocity_catalog_is_the_audited_sdssv_release():
    module = _load_module()

    assert module.DEFAULT_CLUSTER_VELOCITIES_PATH.name == (
        "cluster_velocities_sdssv_covariance_audited.csv"
    )
    assert "velocity analysis/outputs/release" in str(
        module.DEFAULT_CLUSTER_VELOCITIES_PATH
    )


def test_state_only_scene_records_velocity_provenance(tmp_path):
    module = _load_module()
    velocity_path = tmp_path / "cluster_velocities_sdssv_covariance_audited.csv"
    source = {"volumes": {"layers": []}, "initial_state": {}}

    scene = module.build_state_only_scene(source, velocity_path)

    provenance = scene["provenance"]["cluster_velocities"]
    assert provenance["filename"] == velocity_path.name
    assert provenance["path"] == str(velocity_path.resolve())
    assert provenance["release"] == "SDSS-V covariance audited"


def test_velocity_source_build_forwards_the_selected_catalog(tmp_path, monkeypatch):
    module = _load_module()
    velocity_path = tmp_path / "audited.csv"
    output_path = tmp_path / "source.html"
    velocity_path.write_text("name,U_2026,V_2026,W_2026\nA,1,2,3\n", encoding="utf-8")
    captured = {}

    def fake_build_source_main_figure(**kwargs):
        captured.update(kwargs)
        output_path.write_text("<!doctype html>", encoding="utf-8")
        return output_path

    monkeypatch.setattr(
        module,
        "build_source_main_figure",
        fake_build_source_main_figure,
    )

    result = module.build_velocity_source_figure(velocity_path, output_path)

    assert result == output_path
    assert captured["cluster_velocities_path"] == velocity_path.resolve()
    assert captured["chronos_results_path"] == module.JULY4_CHRONOS_RESULTS_PATH
    assert captured["jun6_catalog"] is True
    assert captured["include_background_cluster_trace"] is False
    assert captured["show_cluster_members_in_sky"] is True


def test_july25_artifact_keeps_presentation_and_adds_runtime_upgrades(tmp_path):
    module = _load_module()
    html = ARTIFACT_PATH.read_text(encoding="utf-8")
    scene = module.read_embedded_scene_spec(ARTIFACT_PATH)

    assert ARTIFACT_PATH.stat().st_size < 100 * 1024 * 1024
    # Updated trajectories can have a slightly different compression ratio
    # than the July 21 payload; enforce a direct compact-artifact bound rather
    # than requiring unrelated scientific data to gzip to the same size.
    assert ARTIFACT_PATH.stat().st_size < 40 * 1024 * 1024
    assert scene["deck"]["available"] is True
    assert scene["deck"]["enabled"] is False
    assert scene["deck"]["slides"] == []
    assert scene["initial_state"]["global_controls"]["size_points_by_stars_enabled"] is True
    assert scene["initial_state"]["global_controls"]["fade_opacity_by_birth_time_enabled"] is True
    assert ">Slides ▸</button>" in html
    assert scene["title"] == ""
    velocity_provenance = scene["provenance"]["cluster_velocities"]
    assert velocity_provenance["filename"] == (
        "cluster_velocities_sdssv_covariance_audited.csv"
    )
    assert velocity_provenance["release"] == "SDSS-V covariance audited"
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
