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
    assert scene["sky_dome"]["layer_groups"][0] == {
        "key": "all",
        "label": "All",
        "surveys": [],
    }
    assert scene["sky_dome"]["layer_groups"][1]["label"] == "All-Sky"
    assert scene["sky_dome"]["layer_groups"][1]["surveys"] == [
        "P/Mellinger/color",
        "P/DSS2/color",
        "P/PLANCK/R2/HFI/color",
    ]
    assert [group["label"] for group in scene["sky_dome"]["layer_groups"]] == [
        "All",
        "All-Sky",
        "Optical Surveys",
        "Infrared & UV",
    ]
    assert len(scene["initial_state"]["sky_layers"]) == len(
        module.SKY_BACKGROUND_LAYERS
    )
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
    assert module.DEFAULT_RATZENBOECK_SCOCEN_CATALOG_PATH.name == "ratzi_2022.csv"
    assert module.DEFAULT_RATZENBOECK_SCOCEN_METADATA_PATH.name == (
        "ScoCen_SigMA2_groups_averages_34+KERR3_CenFar-to-Cham_37clsuters.fits"
    )
    assert module.DEFAULT_RATZENBOECK_SCOCEN_MEMBERS_PATH.name == (
        "ScoCen_Sept2022.fits"
    )


def test_ratzenboeck_scocen_member_catalog_keeps_all_gaia_sources():
    module = _load_module()
    grouped = module.load_ratzenboeck_scocen_sky_members(
        module.DEFAULT_RATZENBOECK_SCOCEN_MEMBERS_PATH,
        module.DEFAULT_RATZENBOECK_SCOCEN_METADATA_PATH,
    )

    members = [member for cluster_members in grouped.values() for member in cluster_members]
    assert len(grouped) == 37
    assert len(members) == module.RATZENBOECK_SCOCEN_MEMBER_COUNT
    assert len({member["source_id"] for member in members}) == len(members)
    assert all(member["is_cluster_member"] is True for member in members)


def test_state_only_scene_records_velocity_provenance(tmp_path):
    module = _load_module()
    velocity_path = tmp_path / "cluster_velocities_sdssv_covariance_audited.csv"
    source = {"volumes": {"layers": []}, "initial_state": {}}

    scene = module.build_state_only_scene(source, velocity_path)

    provenance = scene["provenance"]["cluster_velocities"]
    assert provenance["filename"] == velocity_path.name
    assert provenance["path"] == str(velocity_path.resolve())
    assert provenance["release"] == "SDSS-V covariance audited"


def test_state_only_scene_keeps_optional_present_day_volumes_hidden_and_greyscale():
    module = _load_module()
    source = {
        "volumes": {
            "layers": [
                {
                    "key": "volume-0",
                    "state_key": "volume-0",
                    "name": "Edenhofer+2024 Dust",
                    "default_controls": {},
                },
                {
                    "key": "mccallum-ne",
                    "state_key": "mccallum-ne",
                    "name": "McCallum Hα",
                    "default_controls": {},
                },
                {
                    "key": "vergely-dust",
                    "state_key": "vergely-dust",
                    "name": "Vergely 3D Dust",
                    "default_controls": {"colormap": "magma"},
                },
            ],
        },
        "initial_state": {
            "volume_state_by_key": {
                "volume-0": {"visible": False},
                "mccallum-ne": {"visible": True, "showAllTimes": True},
                "vergely-dust": {"visible": True, "showAllTimes": True},
            },
        },
    }

    scene = module.build_state_only_scene(source)
    layers = {
        layer["state_key"]: layer
        for layer in scene["volumes"]["layers"]
    }
    states = scene["initial_state"]["volume_state_by_key"]

    assert states["volume-0"]["visible"] is True
    assert states["mccallum-ne"] == {"visible": False, "showAllTimes": False}
    assert states["vergely-dust"] == {
        "visible": False,
        "showAllTimes": False,
        "colormap": "Greys",
    }
    for key in ("mccallum-ne", "vergely-dust"):
        assert layers[key]["only_at_t0"] is True
        assert layers[key]["supports_show_all_times"] is False
    assert layers["vergely-dust"]["default_controls"]["colormap"] == "Greys"


def test_state_only_scene_records_ratzenboeck_scocen_provenance(tmp_path):
    module = _load_module()
    catalog_path = tmp_path / "ratzi_2022.csv"
    metadata_path = tmp_path / "ratzenboeck_sigma_cluster_metadata.fits"
    source = {"volumes": {"layers": []}, "initial_state": {}}

    scene = module.build_state_only_scene(
        source,
        ratzenboeck_scocen_catalog_path=catalog_path,
        ratzenboeck_scocen_metadata_path=metadata_path,
    )

    provenance = scene["provenance"]["ratzenboeck_scocen"]
    assert provenance["catalog_filename"] == catalog_path.name
    assert provenance["metadata_filename"] == metadata_path.name
    assert provenance["n_clusters"] == 37
    assert "A&A 677, A59" in provenance["membership_reference"]
    assert "A&A 678, A71" in provenance["age_reference"]


def test_velocity_source_build_forwards_the_selected_catalog(tmp_path, monkeypatch):
    module = _load_module()
    velocity_path = tmp_path / "audited.csv"
    ratzenboeck_path = tmp_path / "ratzi_2022.csv"
    ratzenboeck_metadata_path = tmp_path / "ratzenboeck_metadata.fits"
    ratzenboeck_members_path = tmp_path / "ratzenboeck_members.fits"
    output_path = tmp_path / "source.html"
    velocity_path.write_text("name,U_2026,V_2026,W_2026\nA,1,2,3\n", encoding="utf-8")
    ratzenboeck_path.write_text(
        "name,age_myr,x,y,z,U,V,W\nA,1,1,2,3,4,5,6\n",
        encoding="utf-8",
    )
    ratzenboeck_metadata_path.write_text(
        "sigma_id,sigma_region,ratzenboeck_group,n_sigma_sources\n1,US,A,10\n",
        encoding="utf-8",
    )
    ratzenboeck_members_path.write_text("placeholder", encoding="utf-8")
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

    result = module.build_velocity_source_figure(
        velocity_path,
        output_path,
        ratzenboeck_scocen_catalog_path=ratzenboeck_path,
        ratzenboeck_scocen_metadata_path=ratzenboeck_metadata_path,
        ratzenboeck_scocen_members_path=ratzenboeck_members_path,
    )

    assert result == output_path
    assert captured["cluster_velocities_path"] == velocity_path.resolve()
    assert captured["chronos_results_path"] == module.JULY4_CHRONOS_RESULTS_PATH
    assert captured["jun6_catalog"] is True
    assert captured["mobile_safe_mode"] is False
    assert captured["lookback_myr"] == 60
    assert captured["include_background_cluster_trace"] is False
    assert captured["show_cluster_members_in_sky"] is True
    assert captured["include_mccallum_halpha_volume"] is True
    assert captured["include_vergely_dust_volume"] is True
    assert captured["ratzenboeck_scocen_catalog_path"] == ratzenboeck_path.resolve()
    assert captured["ratzenboeck_scocen_metadata_path"] == (
        ratzenboeck_metadata_path.resolve()
    )


def test_july25_artifact_keeps_presentation_and_adds_runtime_upgrades(tmp_path):
    module = _load_module()
    html = ARTIFACT_PATH.read_text(encoding="utf-8")
    scene = module.read_embedded_scene_spec(ARTIFACT_PATH)

    assert scene["deck"]["available"] is True
    assert scene["deck"]["enabled"] is False
    assert scene["deck"]["slides"] == []
    assert len(scene["frames"]) == 61
    assert scene["frames"][0]["time"] == -60.0
    assert scene["frames"][-1]["time"] == 0.0
    assert scene["initial_state"]["global_controls"]["size_points_by_stars_enabled"] is True
    assert scene["initial_state"]["global_controls"]["fade_opacity_by_birth_time_enabled"] is True
    sky_layers = scene["initial_state"]["sky_layers"]
    assert [layer["survey"] for layer in sky_layers[:2]] == [
        "P/Mellinger/color",
        "P/DSS2/color",
    ]
    assert len(sky_layers) == len(module.SKY_BACKGROUND_LAYERS)
    assert {layer["survey"] for layer in sky_layers} == {
        layer["survey"] for layer in module.SKY_BACKGROUND_LAYERS
    }
    assert [layer["survey"] for layer in sky_layers if layer["visible"]] == [
        "P/Mellinger/color",
        "P/PLANCK/R2/HFI/color",
    ]
    assert scene["sky_dome"]["layer_groups"] == module.SKY_BACKGROUND_GROUPS
    assert 'data-startup-ready="false"' in html
    assert 'type: "oviz-aladin-sky-layers-ready"' in html
    assert 'type: "oviz-sky-layer-preload"' in html
    assert "const expectedCount = skyLayersForCurrentGroup().length;" in html
    assert "waitForInitialSkyLayers({ timeoutMs: 30000.0 }).then(" in html
    assert "await waitForInitialSkyLayers({ timeoutMs: 30000.0 })" not in html
    assert ">Slides ▸</button>" in html
    assert scene["title"] == ""
    velocity_provenance = scene["provenance"]["cluster_velocities"]
    assert velocity_provenance["filename"] == (
        "cluster_velocities_sdssv_covariance_audited.csv"
    )
    assert velocity_provenance["release"] == "SDSS-V covariance audited"
    present_day_traces = {
        trace["name"]: trace
        for trace in scene["frames"][-1]["traces"]
    }
    assert len(present_day_traces["Full Cluster Catalog"]["points"]) == 3944
    assert "Clusters (0-150 Myr)" not in present_day_traces
    assert len(present_day_traces["Clusters (< 60 Myr)"]["points"]) == 1163
    assert len(present_day_traces["Clusters (< 15 Myr)"]["points"]) == 540
    assert len(
        present_day_traces[module.RATZENBOECK_SCOCEN_TRACE_NAME]["points"]
    ) == 37
    ratzenboeck_trace_key = present_day_traces[
        module.RATZENBOECK_SCOCEN_TRACE_NAME
    ]["key"]
    assert scene["initial_state"]["current_group"] == "Clusters"
    assert not scene["group_visibility"]["Clusters"][ratzenboeck_trace_key]
    assert module.RATZENBOECK_SCOCEN_TRACE_NAME not in scene["group_order"]
    assert (
        module.RATZENBOECK_SCOCEN_TRACE_NAME
        not in scene["group_visibility"]
    )
    assert all(
        visibility[ratzenboeck_trace_key] == (group_name == "All")
        for group_name, visibility in scene["group_visibility"].items()
    )
    ratzenboeck_track_name = f"{module.RATZENBOECK_SCOCEN_TRACE_NAME} Track"
    assert all(
        trace.get("name") != ratzenboeck_track_name
        for frame in scene["frames"]
        for trace in frame["traces"]
    )
    assert all(
        item.get("name") != ratzenboeck_track_name
        for item in scene["legend"]["items"]
    )
    ratzenboeck_provenance = scene["provenance"]["ratzenboeck_scocen"]
    assert ratzenboeck_provenance["n_clusters"] == 37
    assert ratzenboeck_provenance["n_member_stars"] == 13_103
    assert ratzenboeck_provenance["member_identifier"] == "Gaia source_id"
    sigma_members = module.load_ratzenboeck_scocen_sky_members(
        module.DEFAULT_RATZENBOECK_SCOCEN_MEMBERS_PATH,
        module.DEFAULT_RATZENBOECK_SCOCEN_METADATA_PATH,
    )
    embedded_members = scene["sky_panel"]["members_by_cluster"]
    assert sum(len(embedded_members[name]) for name in sigma_members) == 13_103
    assert scene["group_visibility"]["All"][ratzenboeck_trace_key]
    edenhofer = next(
        layer
        for layer in scene["volumes"]["layers"]
        if "edenhofer" in str(layer.get("name") or "").lower()
    )
    edenhofer_state = scene["initial_state"]["volume_state_by_key"][edenhofer["state_key"]]
    assert edenhofer["default_controls"]["lighting_mode"] == "standard"
    assert edenhofer_state["lightingMode"] == "standard"
    assert edenhofer_state["galacticCenter"] == [8122.0, 0.0, 0.0]
    optional_volumes = {
        layer["state_key"]: layer
        for layer in scene["volumes"]["layers"]
        if layer.get("state_key") in {"mccallum-ne", "vergely-dust"}
    }
    assert set(optional_volumes) == {"mccallum-ne", "vergely-dust"}
    assert optional_volumes["mccallum-ne"]["name"] == "McCallum Hα"
    assert optional_volumes["vergely-dust"]["default_controls"]["colormap"] == "Greys"
    for state_key, layer in optional_volumes.items():
        assert layer["only_at_t0"] is True
        assert layer["supports_show_all_times"] is False
        optional_state = scene["initial_state"]["volume_state_by_key"][state_key]
        assert optional_state["visible"] is False
        assert optional_state["showAllTimes"] is False

    # Presentation behavior carried over from the July 21 artifact.
    assert "const experimentalVolumeLightingControlsEnabled = false;" in html
    assert "function ovizNavigatePresentation(direction)" in html
    assert "ovizStatesPresentationPrevious()" in html
    assert "ovizStatesPresentationNext()" in html
    assert "async function ovizApplyPresentOnlyInitialState()" in html
    assert "await ovizApplyPresentOnlyInitialState();" in html
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
