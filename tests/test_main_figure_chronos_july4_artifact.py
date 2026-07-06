import base64
import gzip
import json
import re
import shutil
import subprocess
import textwrap
from pathlib import Path

import pytest


ARTIFACT_HTML = Path(__file__).with_name("main_figure_chronos_july4.html")
MAX_HTML_SIZE_BYTES = 100 * 1024 * 1024
MAX_RAW_SCENE_SIZE_BYTES = 220 * 1024 * 1024
EXPECTED_FRAME_TIMES = [float(value) for value in range(-120, 1)]
MAX_BACKGROUND_POINTS = 500
MAX_BLUE_CLUSTER_POINTS = 650
MAX_MOBILE_VOLUME_AXIS_PIXELS = 64
MIN_DESKTOP_VOLUME_AXIS_PIXELS = 256


def _read_scene_spec(path: Path):
    html = path.read_text(encoding="utf-8")
    payload_id = None
    for payload_arg in re.findall(r"readOvizSceneSpecPayload\((.*?)\)", html):
        try:
            candidate = json.loads(payload_arg)
        except json.JSONDecodeError:
            continue
        if isinstance(candidate, str) and candidate:
            payload_id = candidate
            break
    assert payload_id is not None, "missing compressed scene payload marker"
    chunk_matches = re.findall(
        (
            r"<script\b(?=[^>]*type=[\"']application/octet-stream[\"'])"
            r"(?=[^>]*data-oviz-payload-id=[\"']"
            + re.escape(payload_id)
            + r"[\"'])[^>]*data-oviz-payload-index=[\"'](\d+)[\"'][^>]*>(.*?)</script>"
        ),
        html,
        re.S,
    )
    assert chunk_matches, "missing compressed scene payload chunks"
    encoded = "".join(
        re.sub(r"\s+", "", chunk_text)
        for _, chunk_text in sorted(chunk_matches, key=lambda item: int(item[0]))
    )
    raw_scene = gzip.decompress(base64.b64decode(encoded))
    return html, raw_scene, json.loads(raw_scene)


def _extract_mobile_detection_js(html: str) -> str:
    start = html.index("function ovizQueryFlagValue")
    end = html.index("      const minimalModeEnabled", start)
    return textwrap.dedent(html[start:end])


@pytest.mark.skipif(not ARTIFACT_HTML.exists(), reason="main_figure_chronos_july4.html has not been generated")
def test_main_figure_chronos_july4_artifact_is_mobile_safe():
    html, raw_scene, scene_spec = _read_scene_spec(ARTIFACT_HTML)

    assert ARTIFACT_HTML.stat().st_size <= MAX_HTML_SIZE_BYTES
    assert len(raw_scene) <= MAX_RAW_SCENE_SIZE_BYTES
    assert len(json.dumps(scene_spec["image_planes"], separators=(",", ":"))) <= 512 * 1024

    assert scene_spec["mobile"]["enabled"] is False
    initial_state = scene_spec["initial_state"]
    assert initial_state["compact_payload_enabled"] is True
    assert initial_state["compact_widget_payload_enabled"] is True
    assert "mobile_mode_enabled" not in initial_state
    assert not any(
        isinstance(trace.get("color_by"), dict) and "colormap_options" in trace["color_by"]
        for frame in scene_spec["frames"]
        for trace in frame["traces"]
    )
    assert any(
        isinstance(item.get("color_by"), dict) and item["color_by"].get("colormap_options")
        for item in scene_spec["legend"]["items"]
    )

    assert [frame["time"] for frame in scene_spec["frames"]] == EXPECTED_FRAME_TIMES
    assert scene_spec["timeline"]["frame_count"] == len(EXPECTED_FRAME_TIMES)

    first_frame_traces = {
        trace["name"]: trace
        for trace in scene_spec["frames"][0]["traces"]
    }
    assert len(first_frame_traces["Full Cluster Catalog"]["points"]) <= MAX_BACKGROUND_POINTS
    assert len(first_frame_traces["Clusters (0-150 Myr)"]["points"]) <= MAX_BACKGROUND_POINTS
    assert len(first_frame_traces["Clusters (< 60 Myr)"]["points"]) <= MAX_BLUE_CLUSTER_POINTS
    assert len(first_frame_traces["Clusters (< 60 Myr)"]["points"]) >= 600
    assert len(first_frame_traces["Clusters (< 15 Myr)"]["points"]) >= 400

    sky_panel = scene_spec["sky_panel"]
    assert sky_panel["enabled"] is True
    assert len(sky_panel.get("members_by_cluster", {})) > 1000
    assert scene_spec["sky_dome"]["enabled"] is True
    assert scene_spec["sky_dome"]["source"] == "aladin"
    assert scene_spec["sky_dome"]["background_mode"] == "live_aladin"

    for widget_key in ("age_kde", "cluster_filter", "dendrogram"):
        assert scene_spec[widget_key]["enabled"] is False

    volumes = scene_spec["volumes"]["layers"]
    assert {volume["state_key"] for volume in volumes} == {
        "edenhofer-dust-desktop",
        "edenhofer-dust-mobile",
    }
    assert "vergely" not in json.dumps(scene_spec["volumes"], separators=(",", ":")).lower()
    volume_by_key = {volume["state_key"]: volume for volume in volumes}
    desktop_volume = volume_by_key["edenhofer-dust-desktop"]
    mobile_volume = volume_by_key["edenhofer-dust-mobile"]
    assert desktop_volume["variant_group"] == "edenhofer-dust-resolution"
    assert mobile_volume["variant_group"] == desktop_volume["variant_group"]
    assert desktop_volume["base_state_name"] == "Edenhofer+2024 Dust"
    assert mobile_volume["base_state_name"] == desktop_volume["base_state_name"]
    assert desktop_volume["variant_order"] < mobile_volume["variant_order"]
    assert desktop_volume["data_encoding"] == "uint8"
    assert desktop_volume.get("data_atlas_tiles") is None
    assert mobile_volume["data_encoding"] == "png_atlas_uint8"
    assert mobile_volume.get("data_atlas_tiles")
    desktop_shape = desktop_volume["shape"]
    mobile_shape = mobile_volume["shape"]
    assert max(int(desktop_shape["x"]), int(desktop_shape["y"])) >= MIN_DESKTOP_VOLUME_AXIS_PIXELS
    assert max(int(mobile_shape["x"]), int(mobile_shape["y"]), int(mobile_shape["z"])) <= MAX_MOBILE_VOLUME_AXIS_PIXELS
    assert int(desktop_shape["x"]) > int(mobile_shape["x"])
    assert desktop_volume["visible"] is True
    assert mobile_volume["visible"] is False

    assert initial_state["active_volume_key"] == "edenhofer-dust-desktop"
    assert initial_state["mobile_active_volume_key"] == "edenhofer-dust-mobile"
    assert initial_state["mobile_defer_volumes"] is False
    assert initial_state["volume_state_by_key"]["edenhofer-dust-desktop"]["visible"] is True
    assert initial_state["volume_state_by_key"]["edenhofer-dust-mobile"]["visible"] is False

    assert 'data-mobile="false"' in html
    assert "ovizRuntimeLooksMobile" in html
    assert "ovizMobileModeOverride" in html
    assert 'ovizQueryFlagValue("desktop", "ovizDesktop", "desktopMode")' in html
    assert 'ovizQueryFlagValue("mobile", "ovizMobile", "mobileMode")' in html
    assert "/iPhone|iPod/i.test(userAgent)" in html
    assert "/Android/i.test(userAgent) && /Mobile/i.test(userAgent)" in html
    assert 'window.matchMedia("(hover: none) and (pointer: coarse)")' in html
    assert "const narrowViewport = Math.min(viewportWidth, viewportHeight) <= 760;" in html
    assert "return Boolean(isiPhoneLike || isAndroidPhone || (coarsePointer && narrowViewport));" in html
    assert "mobileVolumesDeferred" in html
    assert "mobile_active_volume_key" in html
    assert "setExclusiveVolumeVariantSelection(activeVolumeKey)" in html
    assert "oviz-three-mobile-sky-view" in html
    assert "oviz-three-mobile-lasso" in html
    assert "oviz-three-mobile-ar" in html
    assert 'class="oviz-three-mobile-ar" type="button" title="Open AR export options"' in html
    assert 'class="oviz-three-mobile-ar" type="button" title="Open AR export options" aria-haspopup="dialog" aria-expanded="false" data-has-selection="false" disabled' not in html
    assert "oviz-three-mobile-legend" in html
    assert "function collectOvizArSnapshot" in html
    assert "function ovizArPresentFrameIndex" in html
    assert "const OVIZ_AR_SKY_TEXTURE_HIGH = { width: 4096, height: 2048 };" in html
    assert "const OVIZ_AR_SKY_TEXTURE_LOW = { width: 2048, height: 1024 };" in html
    assert "function buildOvizArSkyDomeScene" in html
    assert "skyDomeHips2FitsUrl(width, height)" in html
    assert "function loadOvizUSDZExporter" in html
    assert "function normalizeOvizArMaterialForUSDZ" in html
    assert 'const textureSlots = ["map", "normalMap", "aoMap", "roughnessMap", "metalnessMap", "emissiveMap"];' in html
    assert "No selection; exporting a capped t=0 Myr scene snapshot." in html
    assert "USDZExporter" in html
    assert 'rel="ar"' in html


@pytest.mark.skipif(not ARTIFACT_HTML.exists(), reason="main_figure_chronos_july4.html has not been generated")
@pytest.mark.skipif(shutil.which("node") is None, reason="node is not available")
def test_main_figure_chronos_july4_artifact_executes_mobile_detection():
    html, _raw_scene, _scene_spec = _read_scene_spec(ARTIFACT_HTML)
    mobile_detection_js = _extract_mobile_detection_js(html)
    script = f"""
    {mobile_detection_js}

    function evaluateCase(options) {{
      const width = options.width ?? 1280;
      const height = options.height ?? 720;
      globalThis.window = {{
        location: {{ search: options.search || "" }},
        navigator: {{
          userAgent: options.userAgent || "",
          maxTouchPoints: options.maxTouchPoints || 0,
        }},
        innerWidth: width,
        innerHeight: height,
        screen: {{
          width,
          height,
          availWidth: width,
          availHeight: height,
        }},
        matchMedia: () => ({{ matches: Boolean(options.coarsePointer) }}),
      }};
      const override = ovizMobileModeOverride();
      return override === null ? ovizRuntimeLooksMobile() : override;
    }}

    const results = {{
      desktop: evaluateCase({{ userAgent: "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7)" }}),
      iphone: evaluateCase({{ userAgent: "Mozilla/5.0 (iPhone; CPU iPhone OS 17_0 like Mac OS X)" }}),
      androidPhone: evaluateCase({{ userAgent: "Mozilla/5.0 (Linux; Android 14; Pixel 8) Mobile" }}),
      coarseNarrow: evaluateCase({{ coarsePointer: true, maxTouchPoints: 5, width: 390, height: 844 }}),
      coarseWide: evaluateCase({{ coarsePointer: true, maxTouchPoints: 5, width: 1024, height: 900 }}),
      mobileOverride: evaluateCase({{ search: "?mobile=1" }}),
      desktopOverride: evaluateCase({{ search: "?desktop=1", userAgent: "Mozilla/5.0 (iPhone)" }}),
      mobileFalseOverride: evaluateCase({{ search: "?mobile=0", userAgent: "Mozilla/5.0 (iPhone)" }}),
    }};
    process.stdout.write(JSON.stringify(results));
    """
    result = subprocess.run(
        ["node"],
        input=script,
        text=True,
        capture_output=True,
        check=True,
    )
    assert json.loads(result.stdout) == {
        "desktop": False,
        "iphone": True,
        "androidPhone": True,
        "coarseNarrow": True,
        "coarseWide": False,
        "mobileOverride": True,
        "desktopOverride": False,
        "mobileFalseOverride": False,
    }
