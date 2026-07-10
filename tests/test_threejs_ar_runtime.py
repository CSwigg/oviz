import json
import shutil
import subprocess
import sys
import textwrap
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from oviz.threejs_runtime_ar import THREEJS_AR_RUNTIME_JS


@pytest.mark.skipif(shutil.which("node") is None, reason="node is not available")
def test_ar_snapshot_uses_present_day_selection_and_sky_directions():
    script = f"""
    (async () => {{
    const sceneSpec = {{ initial_frame_index: 0 }};
    const frameSpecs = [
      {{
        time: -1,
        traces: [{{
          name: "Trace A",
          points: [{{
            x: 99,
            y: 0,
            z: 0,
            n_stars: 4,
            selection: {{
              cluster_name: "A",
              trace_name: "Trace A",
              cluster_color: "#ff0000",
              ra_deg: 11,
              dec_deg: 22,
              l_deg: 33,
              b_deg: 44,
            }},
          }}],
        }}],
      }},
      {{
        time: 0,
        traces: [{{
          name: "Trace A",
          points: [{{
            x: 5,
            y: 6,
            z: 7,
            n_stars: 9,
            selection: {{
              cluster_name: "A",
              trace_name: "Trace A",
              cluster_color: "#00ff00",
              ra_deg: 90,
              dec_deg: 0,
              l_deg: 0,
              b_deg: 90,
            }},
          }}],
        }}],
      }},
    ];
    let currentSelection = null;
    let currentSelections = [{{ cluster_name: "A", trace_name: "Trace A" }}];
    let selectedClusterKeys = new Set(["a"]);
    let currentLassoSelectionMask = null;
    let minimalModeEnabled = false;
    let ovizTestCoordsys = "galactic";
    let activeVolumeKey = "volume-0";
    const volumeLayers = [{{
      key: "volume-0",
      state_key: "volume-0",
      name: "Test Volume",
      state_name: "Test Volume",
      shape: {{ x: 2, y: 2, z: 1 }},
      bounds: {{ x: [0, 10], y: [0, 10], z: [0, 10] }},
      data_encoding: "uint8",
      ar_proxy: {{
        data_b64: "AP+AQA==",
        data_encoding: "uint8",
        shape: {{ x: 2, y: 2, z: 1 }},
        method: "block_max",
      }},
      colormap_options: [{{ name: "test", lut_b64: "" }}],
    }}];
    const volumeStateByKey = {{
      "volume-0": {{
        visible: true,
        vmin: 0,
        vmax: 1,
        opacity: 0.4,
        stretch: "linear",
        colormap: "test",
      }},
    }};
    const legendState = {{ "volume-0": true }};
    const volumeScalarDataCache = new Map();
    const volumeScalarDataPendingCache = new Map();
    const plotGroup = {{ position: {{ x: 0, y: 0, z: 0 }} }};
    const THREE = {{
      Vector3: class {{
        constructor(x = 0, y = 0, z = 0) {{
          this.x = x;
          this.y = y;
          this.z = z;
        }}
      }},
      Color: class {{
        constructor(value) {{
          this.value = value;
        }}
        getHex() {{
          return this.value === 0x000000 ? 0 : 1;
        }}
      }},
    }};
    const mobileArButtonEl = {{
      disabled: null,
      dataset: {{}},
      attrs: {{}},
      title: "",
      setAttribute(key, value) {{
        this.attrs[key] = value;
      }},
    }};
    const root = {{ appendChild: () => {{}} }};
    let scalarReadCount = 0;

    function normalizeMemberKey(value) {{
      return String(value || "").trim().toLowerCase().replace(/\\s+/g, "_");
    }}
    function selectionIdentityKeyFor(selection) {{
      return selection && selection.cluster_name ? String(selection.cluster_name) : "";
    }}
    function normalizedSelectionKeyFor(selection) {{
      return normalizeMemberKey(selectionIdentityKeyFor(selection));
    }}
    function selectionMetadataForKey(_key) {{
      return null;
    }}
    function selectionForPoint(point, _trace) {{
      return point.selection || null;
    }}
    function pointBaseColorForTrace(point, trace) {{
      return point.color || trace.default_color || "#ffffff";
    }}
    function skyDomeHips2FitsCoordsys() {{
      return ovizTestCoordsys;
    }}
    function skyDomeHips2FitsUrl(width, height) {{
      return `https://example.test/hips2fits?hips=P%2FPLANCK&width=${{width}}&height=${{height}}&min_cut=1`;
    }}
    function activeSkyLayer() {{
      return {{ survey: "P/Mellinger/color", visible: true, stretch: "linear", cutMin: "", cutMax: "" }};
    }}
    function normalizeSkyLongitude(value) {{
      let lon = Number(value) || 0;
      lon %= 360;
      if (lon < 0) lon += 360;
      return lon;
    }}
    function icrsDegFromGalacticDeg(lDeg, bDeg) {{
      return {{ ra: Number(lDeg), dec: Number(bDeg) }};
    }}
    function galacticDegFromIcrsDeg(raDeg, decDeg) {{
      return {{ l: Number(raDeg), b: Number(decDeg) }};
    }}
    function hasActiveLassoSelectionMask() {{
      return Boolean(currentLassoSelectionMask);
    }}
    function pointInsideProjectedLassoMask(worldX, _worldY, _worldZ, _mask) {{
      return Number(worldX) > 0;
    }}
    function activeVolumeLassoSelectionMask() {{
      return currentLassoSelectionMask;
    }}
    function volumeStateKeyForLayer(layer) {{
      return String((layer && (layer.state_key || layer.key)) || "");
    }}
    function volumeStateNameForLayer(layer) {{
      return String((layer && (layer.state_name || layer.name)) || "");
    }}
    function volumeBaseNameForLayer(layer) {{
      return volumeStateNameForLayer(layer);
    }}
    function frameVolumeLayers(_frame) {{
      return volumeLayers;
    }}
    function frameVolumeLayerForStateKey(stateKey, _frame) {{
      return volumeLayers.find((layer) => volumeStateKeyForLayer(layer) === String(stateKey)) || null;
    }}
    function volumeLayerForKey(layerKey) {{
      return frameVolumeLayerForStateKey(layerKey, null)
        || volumeLayers.find((layer) => String(layer.key) === String(layerKey))
        || null;
    }}
    function volumeVisibleForFrame(_layer, state, _frame) {{
      return state && state.visible !== false;
    }}
    function volumeScalarArrayFor(_layer) {{
      scalarReadCount += 1;
      throw new Error("full volume data should not be read when an AR proxy exists");
    }}
    function base64ToUint8Array(encoded) {{
      return new Uint8Array(Buffer.from(String(encoded || ""), "base64"));
    }}
    function normalizedVolumeWindowFor(_layer, _state) {{
      return {{ low: 0.0, high: 1.0 }};
    }}
    function volumeColormapOptionFor(layer, _colormapName) {{
      return layer.colormap_options[0];
    }}
    function volumeColorBytesForOption(_option) {{
      return new Uint8Array([
        0, 0, 0, 255,
        64, 64, 180, 255,
        180, 120, 64, 255,
        255, 255, 255, 255,
      ]);
    }}
    function focusViewer() {{}}

    {THREEJS_AR_RUNTIME_JS}

    const snapshot = collectOvizArSnapshot("3d");
    renderArSnapshotButtonState();
    const selectedButtonState = {{
      disabled: mobileArButtonEl.disabled,
      hasSelection: mobileArButtonEl.dataset.hasSelection,
      title: mobileArButtonEl.title,
    }};
    const lon0 = ovizArSkyDirectionForLonLatDeg(0, 0, 1);
    const lon90 = ovizArSkyDirectionForLonLatDeg(90, 0, 1);
    const lat90 = ovizArSkyDirectionForLonLatDeg(0, 90, 1);
    ovizTestCoordsys = "icrs";
    const icrsPoint = ovizArSkyDirectionForPoint({{ ra: 90, dec: 0, l: 0, b: 90 }}, 1);
    const activeSkyUrl = ovizArSkyTextureUrl(2048, 1024);
    const arVector = ovizArVectorFromPoint(
      {{ x: 1, y: 2, z: 3 }},
      {{ center: {{ x: 0, y: 0, z: 0 }}, scale: 1 }}
    );
    selectedClusterKeys = new Set();
    currentSelections = [];
    currentSelection = null;
    const emptySnapshot = collectOvizArSnapshot("3d");
    renderArSnapshotButtonState();
    currentLassoSelectionMask = {{ mask: true }};
    const maskOnlySnapshot = collectOvizArSnapshot("3d");
    const volumeResult = await ovizArCollectVolumeSamples(maskOnlySnapshot, {{ maxSamples: 3, scanTarget: 100 }});
    currentLassoSelectionMask = null;
    const sourceBuffer = new Uint8Array([0, 1, 2, 3, 4, 5]).buffer;
    const slicedViewBuffer = await ovizArArrayBufferFromExporterResult(new Uint8Array(sourceBuffer, 2, 3));
    const validUsdZBuffer = ovizArValidateUsdZArrayBuffer(new ArrayBuffer(OVIZ_AR_MIN_USDZ_BYTES + 32));
    let emptyUsdZError = "";
    try {{
      ovizArValidateUsdZArrayBuffer(new ArrayBuffer(0));
    }} catch (err) {{
      emptyUsdZError = err && err.message ? String(err.message) : String(err);
    }}
    const zipBytes = ovizUsdZWriteStoredZip([{{
      name: "model.usda",
      data: new TextEncoder().encode("#usda 1.0\\n"),
    }}]);
    const zipView = new DataView(zipBytes.buffer, zipBytes.byteOffset, zipBytes.byteLength);
    const zipNameLength = zipView.getUint16(26, true);
    const zipExtraLength = zipView.getUint16(28, true);
    const zipDataOffset = 30 + zipNameLength + zipExtraLength;

    process.stdout.write(JSON.stringify({{
      presentTimeMyr: snapshot.presentTimeMyr,
      presentFrameIndex: snapshot.presentFrameIndex,
      pointCount: snapshot.points.length,
      firstX: snapshot.points[0].x,
      firstY: snapshot.points[0].y,
      firstZ: snapshot.points[0].z,
      trailPointCount: snapshot.trails[0].points.length,
      lon0,
      lon90,
      lat90,
      icrsPoint,
      activeSkyUrl,
      arVector,
      emptyPointCount: emptySnapshot.points.length,
      emptySelectionMode: emptySnapshot.selectionMode,
      maskOnlyPointCount: maskOnlySnapshot.points.length,
      maskOnlySelectionMode: maskOnlySnapshot.selectionMode,
      volumeSampleCount: volumeResult.samples.length,
      firstVolumeSample: volumeResult.samples[0],
      volumeLayerCount: volumeResult.layers.length,
      scalarReadCount,
      emptyCanExport: ovizArCanExportSelection(),
      selectedButtonState,
      emptyButtonState: {{
        disabled: mobileArButtonEl.disabled,
        hasSelection: mobileArButtonEl.dataset.hasSelection,
        title: mobileArButtonEl.title,
      }},
      zipSignature: zipView.getUint32(0, true),
      zipDataOffset,
      slicedViewBytes: Array.from(new Uint8Array(slicedViewBuffer)),
      validUsdZBytes: validUsdZBuffer.byteLength,
      emptyUsdZError,
    }}));
    }})().catch((err) => {{
      console.error(err && err.stack ? err.stack : err);
      process.exit(1);
    }});
    """
    result = subprocess.run(
        ["node"],
        input=textwrap.dedent(script),
        text=True,
        capture_output=True,
        check=True,
    )
    payload = json.loads(result.stdout)

    assert payload["presentTimeMyr"] == 0
    assert payload["presentFrameIndex"] == 1
    assert payload["pointCount"] == 1
    assert payload["firstX"] == 5
    assert payload["firstY"] == 6
    assert payload["firstZ"] == 7
    assert payload["trailPointCount"] == 2
    assert payload["emptyPointCount"] == 1
    assert payload["emptySelectionMode"] == "present-day-scene"
    assert payload["maskOnlyPointCount"] == 1
    assert payload["maskOnlySelectionMode"] == "volume-lasso"
    assert payload["volumeSampleCount"] == 3
    assert payload["volumeLayerCount"] == 1
    assert payload["scalarReadCount"] == 0
    assert payload["firstVolumeSample"]["x"] == pytest.approx(7.5)
    assert payload["firstVolumeSample"]["y"] == pytest.approx(2.5)
    assert payload["firstVolumeSample"]["z"] == pytest.approx(5.0)
    assert payload["arVector"] == {"x": 1, "y": 3, "z": -2}
    assert payload["emptyCanExport"] is True
    assert payload["selectedButtonState"]["disabled"] is False
    assert payload["selectedButtonState"]["hasSelection"] == "true"
    assert payload["emptyButtonState"]["disabled"] is False
    assert payload["emptyButtonState"]["hasSelection"] == "false"
    assert "present-day scene" in payload["emptyButtonState"]["title"]
    assert payload["zipSignature"] == 0x04034B50
    assert payload["zipDataOffset"] % 64 == 0
    assert payload["slicedViewBytes"] == [2, 3, 4]
    assert payload["validUsdZBytes"] > 1024
    assert "empty or incomplete" in payload["emptyUsdZError"]

    assert payload["lon0"]["x"] == pytest.approx(-1)
    assert payload["lon0"]["y"] == pytest.approx(0, abs=1e-12)
    assert payload["lon0"]["z"] == pytest.approx(0, abs=1e-12)
    assert payload["lon90"]["x"] == pytest.approx(0, abs=1e-12)
    assert payload["lon90"]["z"] == pytest.approx(1)
    assert payload["lat90"]["y"] == pytest.approx(1)
    assert payload["icrsPoint"]["z"] == pytest.approx(1)
    assert "hips=P%2FMellinger%2Fcolor" in payload["activeSkyUrl"]
    assert "width=2048" in payload["activeSkyUrl"]
    assert "min_cut" not in payload["activeSkyUrl"]


def test_ar_runtime_uses_self_contained_usdz_exporter():
    assert "class OvizUSDZExporter" in THREEJS_AR_RUNTIME_JS
    assert "defaultPrim = \"${defaultPrim}\"" in THREEJS_AR_RUNTIME_JS
    assert 'prepend apiSchemas = ["MaterialBindingAPI"]' in THREEJS_AR_RUNTIME_JS
    assert "sceneAr.updateMatrixWorld(true)" in THREEJS_AR_RUNTIME_JS
    assert 'string inputs:varname = "st"' in THREEJS_AR_RUNTIME_JS
    assert "fflate" not in THREEJS_AR_RUNTIME_JS
    assert "examples/js/exporters/USDZExporter.js" not in THREEJS_AR_RUNTIME_JS
