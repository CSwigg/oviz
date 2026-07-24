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
MAX_IMAGE_PLANES_SIZE_BYTES = 5 * 1024 * 1024
EXPECTED_FRAME_TIMES = [float(value) for value in range(-120, 1)]
MAX_BACKGROUND_POINTS = 500
MAX_BLUE_CLUSTER_POINTS = 650
MIN_DESKTOP_VOLUME_AXIS_PIXELS = 256
EXPECTED_FILTERED_RETAINED_BOUNDARY_COMPONENTS = 6624


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
    image_planes_size = len(json.dumps(scene_spec["image_planes"], separators=(",", ":")))
    assert 512 * 1024 < image_planes_size <= MAX_IMAGE_PLANES_SIZE_BYTES

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

    cluster_visibility = scene_spec["group_visibility"]["Clusters"]
    boundary_visible_points = [
        sum(
            len(trace.get("points", []))
            for trace in scene_spec["frames"][frame_index]["traces"]
            if cluster_visibility.get(trace["key"]) is True
        )
        for frame_index in (119, 120)
    ]
    # Retained filtering keeps the three existing glow/core/marker components
    # only for traces visible at either endpoint. Ordinary rendering is unchanged.
    assert sum(boundary_visible_points) * 3 == EXPECTED_FILTERED_RETAINED_BOUNDARY_COMPONENTS

    first_frame_traces = {
        trace["name"]: trace
        for trace in scene_spec["frames"][0]["traces"]
    }
    assert len(first_frame_traces["Full Cluster Catalog"]["points"]) <= MAX_BACKGROUND_POINTS
    assert "Clusters (0-150 Myr)" not in first_frame_traces
    assert len(first_frame_traces["Clusters (< 60 Myr)"]["points"]) <= MAX_BLUE_CLUSTER_POINTS
    assert len(first_frame_traces["Clusters (< 60 Myr)"]["points"]) >= 600
    assert len(first_frame_traces["Clusters (< 15 Myr)"]["points"]) >= 400

    sky_panel = scene_spec["sky_panel"]
    assert sky_panel["enabled"] is True
    assert sky_panel["show_cluster_members_in_sky"] is True
    assert sky_panel["member_point_size_denominator"] == 20
    assert sky_panel["member_min_screen_size_px"] == 2.5
    assert len(sky_panel.get("members_by_cluster", {})) > 1000
    explicit_members = [
        member
        for members in sky_panel["members_by_cluster"].values()
        for member in members
        if member.get("is_cluster_member") is True
    ]
    assert len(explicit_members) >= 2500
    assert scene_spec["sky_dome"]["enabled"] is True
    assert scene_spec["sky_dome"]["source"] == "aladin"
    assert scene_spec["sky_dome"]["background_mode"] == "live_aladin"
    assert "function addSkyMemberStars(group, catalog, options = {})" in html
    assert "function skyMemberScaleForPoint(basePointScale, position)" in html
    assert "function normalizeClusterCatalogKey(value)" in html
    assert "new THREE.Points(memberGeometry, glowMaterial)" in html
    assert "new THREE.Points(memberGeometry, coreMaterial)" in html
    assert "root.dataset.skyMemberDrawObjectCount = String(" in html
    assert "function animateSkyMemberReveal(targetProgress" in html
    assert 'starGlowTextureFor("sky_member_halo")' in html
    assert 'starCoreTextureFor("sky_member_core")' in html
    assert "skyMemberBulkOpacityEntries" in html
    assert '"sky_member_glow"' in html
    assert '"sky_member_core"' in html
    assert 'root.dataset.skyClusterBulkPointCount = String(' in html
    assert 'root.dataset.skyMemberStarCount = String(renderedSkyMemberStarCount)' in html

    for widget_key in ("age_kde", "cluster_filter", "dendrogram"):
        assert scene_spec[widget_key]["enabled"] is False

    volumes = scene_spec["volumes"]["layers"]
    assert len(volumes) == 1
    assert "vergely" not in json.dumps(scene_spec["volumes"], separators=(",", ":")).lower()
    edenhofer_volume = volumes[0]
    assert edenhofer_volume["name"] == "Edenhofer+2024 Dust"
    assert edenhofer_volume["data_encoding"] == "uint8"
    assert edenhofer_volume.get("data_atlas_tiles") is None
    assert max(int(value) for value in edenhofer_volume["shape"].values()) >= MIN_DESKTOP_VOLUME_AXIS_PIXELS
    assert edenhofer_volume["visible"] is True
    assert edenhofer_volume["default_controls"]["opacity"] == 1.0
    assert edenhofer_volume["ar_proxy"]["method"] == "block_max"
    assert max(int(value) for value in edenhofer_volume["ar_proxy"]["shape"].values()) <= 64
    assert edenhofer_volume["ar_proxy"]["data_b64"]

    assert initial_state["active_volume_key"] == edenhofer_volume["state_key"]
    assert not initial_state.get("mobile_active_volume_key")
    assert scene_spec["ar"] == {"enabled": False}
    assert initial_state["mobile_defer_volumes"] is False
    assert initial_state["volume_state_by_key"][edenhofer_volume["state_key"]] == {
        "visible": True,
        "opacity": 1.0,
        "stretch": "asinh",
        "vmax": 0.07,
    }
    planck_layer = next(
        layer for layer in initial_state["sky_layers"]
        if layer["survey"] == "P/PLANCK/R2/HFI/color"
    )
    assert planck_layer["stretch"] == "asinh"
    assert planck_layer["cut_max"] == 0.07

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
    assert "lockNorthUp: skyDomeBackgroundOnly" in html
    assert "window.OvizSkyBackgroundBridge" in html
    assert "function setMilkyWayModelOpacityScale(value)" in html
    assert "function galacticReferenceMotionVisible()" in html
    assert "function galacticReferenceTimeOpacity()" in html
    assert "function ovizBuildTransitionPhases(" in html
    assert "phaseMinimumDurationMs = 800.0" in html
    assert "actionHeldTraceOpacityByKey = null" in html
    assert "for (let index = firstIndex; index <= lastIndex; index += 1)" in html
    assert "preserveCamera: true" in html
    assert "Math.round(value * 256.0) / 256.0" in html
    assert "renderer.domElement.style.opacity = String(canvasOpacity)" not in html
    assert "const stableFrameIndex = clampFrameIndex(targetIndex)" in html
    assert "updateTimelineMotionOpacity()" in html
    assert "renderInterpolatedFrameValue(targetIndex, { updateWidgets: false })" not in html
    assert "ovizTimeOpacityScale" in html
    assert "setExclusiveVolumeVariantSelection(activeVolumeKey)" in html

    galaxy_opacity_scale_by_time = {}
    for frame in scene_spec["frames"]:
        if float(frame["time"]) not in {0.0, -1.0, -5.0}:
            continue
        decoration = next(
            item for item in frame.get("decorations", [])
            if item.get("key") == "galaxy-image-overlay"
        )
        galaxy_opacity_scale_by_time[float(frame["time"])] = float(decoration["opacity_scale"])
    assert galaxy_opacity_scale_by_time[0.0] == 1.0
    assert galaxy_opacity_scale_by_time[-1.0] == 0.0
    assert galaxy_opacity_scale_by_time[-5.0] == 0.0
    assert "oviz-three-mobile-sky-view" in html
    assert "oviz-three-mobile-lasso" in html
    assert 'data-ar-enabled="false"' in html
    assert '<button class="oviz-three-mobile-ar"' not in html
    assert "oviz-three-mobile-legend" in html
    assert "function collectOvizArSnapshot" not in html
    assert "function ovizArPresentFrameIndex" not in html
    assert "function buildOvizArSkyDomeScene" not in html
    assert "class OvizUSDZExporter" not in html
    assert 'rel="ar"' not in html


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
