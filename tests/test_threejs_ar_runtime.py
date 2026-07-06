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
    let minimalModeEnabled = false;
    let ovizTestCoordsys = "galactic";
    const mobileArButtonEl = null;
    const root = {{ appendChild: () => {{}} }};

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
    function focusViewer() {{}}

    {THREEJS_AR_RUNTIME_JS}

    const snapshot = collectOvizArSnapshot("3d");
    const lon0 = ovizArSkyDirectionForLonLatDeg(0, 0, 1);
    const lon90 = ovizArSkyDirectionForLonLatDeg(90, 0, 1);
    const lat90 = ovizArSkyDirectionForLonLatDeg(0, 90, 1);
    ovizTestCoordsys = "icrs";
    const icrsPoint = ovizArSkyDirectionForPoint({{ ra: 90, dec: 0, l: 0, b: 90 }}, 1);
    selectedClusterKeys = new Set();
    currentSelections = [];
    currentSelection = null;
    const emptySnapshot = collectOvizArSnapshot("3d");

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
      emptyPointCount: emptySnapshot.points.length,
      emptyCanExport: ovizArCanExportSelection(),
    }}));
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
    assert payload["emptyPointCount"] == 0
    assert payload["emptyCanExport"] is False

    assert payload["lon0"]["x"] == pytest.approx(-1)
    assert payload["lon0"]["y"] == pytest.approx(0, abs=1e-12)
    assert payload["lon0"]["z"] == pytest.approx(0, abs=1e-12)
    assert payload["lon90"]["x"] == pytest.approx(0, abs=1e-12)
    assert payload["lon90"]["z"] == pytest.approx(1)
    assert payload["lat90"]["y"] == pytest.approx(1)
    assert payload["icrsPoint"]["z"] == pytest.approx(1)
