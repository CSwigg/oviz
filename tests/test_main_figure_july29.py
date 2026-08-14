import importlib.util
import sys
from pathlib import Path


SCRIPT_PATH = Path(__file__).with_name("main_figure_july29.py")
ARTIFACT_PATH = SCRIPT_PATH.with_suffix(".html")


def _load_module():
    spec = importlib.util.spec_from_file_location("main_figure_july29", SCRIPT_PATH)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_july29_relabels_only_the_primary_cluster_sample():
    module = _load_module()
    source = {
        "title": module.SOURCE_CLUSTER_LABEL,
        "frames": [{"traces": [{"name": module.SOURCE_CLUSTER_LABEL}]}],
        "other": "Clusters (< 15 Myr)",
    }

    updated = module._replace_cluster_label(source)

    assert source["title"] == module.SOURCE_CLUSTER_LABEL
    assert updated["title"] == module.JULY29_CLUSTER_LABEL
    assert updated["frames"][0]["traces"][0]["name"] == (
        module.JULY29_CLUSTER_LABEL
    )
    assert updated["other"] == "Clusters (< 15 Myr)"


def test_scene_script_regeneration_bakes_120_myr_limits(monkeypatch, tmp_path):
    module = _load_module()
    velocity_path = tmp_path / "audited.csv"
    catalog_path = tmp_path / "ratzi.csv"
    metadata_path = tmp_path / "metadata.fits"
    members_path = tmp_path / "members.fits"
    for path in (velocity_path, catalog_path, metadata_path, members_path):
        path.write_text("placeholder", encoding="utf-8")
    script_path = tmp_path / "scene_script.py"
    monkeypatch.setattr(module, "SCENE_SCRIPT_PATH", script_path)

    captured = {}
    original_age_limit = module.source_runner.MAIN_FIGURE_CLUSTER_BLUE_MAX_AGE_MYR

    def fake_build_patched_main_figure_source(**kwargs):
        captured.update(kwargs)
        captured["active_age_limit"] = (
            module.source_runner.MAIN_FIGURE_CLUSTER_BLUE_MAX_AGE_MYR
        )
        return "save_name = '/tmp/unused.html'\nfig3d = None\n"

    monkeypatch.setattr(
        module.source_runner,
        "build_patched_main_figure_source",
        fake_build_patched_main_figure_source,
    )

    result = module.regenerate_scene_script(
        cluster_velocities_path=velocity_path,
        ratzenboeck_scocen_catalog_path=catalog_path,
        ratzenboeck_scocen_metadata_path=metadata_path,
        ratzenboeck_scocen_members_path=members_path,
    )

    assert result == script_path
    assert captured["lookback_myr"] == 120
    assert captured["active_age_limit"] == 120
    assert module.source_runner.MAIN_FIGURE_CLUSTER_BLUE_MAX_AGE_MYR == (
        original_age_limit
    )
    script_text = script_path.read_text(encoding="utf-8")
    assert "save_name = None" in script_text
    assert module._scene_script_velocities(script_text) == str(
        velocity_path.resolve()
    )


def test_build_source_scene_executes_scene_script(monkeypatch, tmp_path):
    module = _load_module()
    velocity_path = tmp_path / "audited.csv"
    velocity_path.write_text("placeholder", encoding="utf-8")
    script_path = tmp_path / "scene_script.py"
    script_path.write_text(
        f"{module.SCENE_SCRIPT_VELOCITIES_MARKER} = "
        f"{str(velocity_path.resolve())!r}\n"
        "class _Figure:\n"
        "    def to_dict(self):\n"
        "        return {'frames': ['from-script']}\n"
        "fig3d = _Figure()\n",
        encoding="utf-8",
    )
    monkeypatch.setattr(module, "SCENE_SCRIPT_PATH", script_path)
    monkeypatch.setattr(module, "scene_spec_jsonable", lambda scene: scene)

    def fail_regenerate(**kwargs):
        raise AssertionError("A matching scene script must not be regenerated")

    monkeypatch.setattr(module, "regenerate_scene_script", fail_regenerate)

    scene = module.build_source_scene(velocity_path)

    assert scene == {"frames": ["from-script"]}


def test_build_scene_never_touches_html(monkeypatch):
    """The July 29 build must stay a pure scene-dict pipeline."""
    module = _load_module()

    monkeypatch.setattr(
        module,
        "build_source_scene",
        lambda *args, **kwargs: {"frames": [], "title": module.SOURCE_CLUSTER_LABEL},
    )
    monkeypatch.setattr(
        module.july25,
        "build_state_only_scene",
        lambda scene, **kwargs: dict(scene),
    )

    def fail_read(*args, **kwargs):
        raise AssertionError("build_scene must not decode HTML")

    monkeypatch.setattr(module.july25, "read_embedded_scene_spec", fail_read)

    scene = module.build_scene(cluster_velocities_path=Path("unused.csv"))

    assert scene["title"] == module.JULY29_CLUSTER_LABEL


def test_july29_artifact_has_120_myr_sample_and_timeline():
    module = _load_module()
    scene = module.july25.read_embedded_scene_spec(ARTIFACT_PATH)

    assert ARTIFACT_PATH.stat().st_size < 100 * 1024 * 1024
    assert len(scene["frames"]) == 121
    assert scene["frames"][0]["time"] == -120.0
    assert scene["frames"][-1]["time"] == 0.0

    present_traces = {
        trace["name"]: trace
        for trace in scene["frames"][-1]["traces"]
    }
    assert module.JULY29_CLUSTER_LABEL in present_traces
    assert module.SOURCE_CLUSTER_LABEL not in present_traces
    assert len(present_traces[module.JULY29_CLUSTER_LABEL]["points"]) > 1163
    assert "Clusters (0-150 Myr)" not in present_traces

    # The ccSN reconstruction figure selects its cluster sample from the
    # full catalog trace, so the July 29 scene must keep it at trace-1.
    full_catalog = present_traces["Full Cluster Catalog"]
    assert full_catalog["key"] == "trace-1"
    assert len(full_catalog["points"]) > 3000

    volumes = {
        layer.get("state_key") or layer["key"]: layer
        for layer in scene["volumes"]["layers"]
    }
    states = scene["initial_state"]["volume_state_by_key"]
    assert states["volume-0"]["visible"] is True
    assert states["mccallum-ne"]["visible"] is False
    assert states["vergely-dust"]["visible"] is False
    assert volumes["vergely-dust"]["default_controls"]["colormap"] == "Greys"
