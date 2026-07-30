from __future__ import annotations

import json
import shutil
import subprocess
import unittest

from oviz.threejs_figure import ThreeJSFigure


def _runtime_html() -> str:
    return ThreeJSFigure(
        {
            "width": 640,
            "height": 480,
            "frames": [],
            "initial_state": {},
        }
    ).to_html(compress_scene_spec=False)


class ThreeJSRetainedSceneTests(unittest.TestCase):
    @unittest.skipIf(shutil.which("node") is None, "node is not available")
    def test_retained_live_motion_cannot_mutate_source_frame_points(self):
        html = _runtime_html()
        interpolation_helpers = (
            "function interpolateNumber(fromValue, toValue, alpha, fallbackValue = 0.0)"
            + html.split(
                "function interpolateNumber(fromValue, toValue, alpha, fallbackValue = 0.0)", 1
            )[1].split("function cloneTraceLabel", 1)[0]
        )
        retained_clone = (
            "function ovizRetainedCloneLivePoint(pair, alpha, displayedTimeMyr)"
            + html.split(
                "function ovizRetainedCloneLivePoint(pair, alpha, displayedTimeMyr)", 1
            )[1].split("function ovizPointBirthVisibility", 1)[0]
        )
        script = f"""
        {interpolation_helpers}
        {retained_clone}
        const source = {{ motion: {{ key: "cluster-a", time_myr: 0 }} }};
        const pair = {{
          from: {{ metadata: {{ point: source }} }},
          to: null,
          livePoint: cloneTracePoint(source),
        }};
        const rendered = ovizRetainedCloneLivePoint(pair, 0.5, -12.5);
        process.stdout.write(JSON.stringify({{
          sourceTime: source.motion.time_myr,
          renderedTime: rendered.motion.time_myr,
          motionIsOwned: rendered.motion !== source.motion,
        }}));
        """
        result = subprocess.run(
            ["node"], input=script, text=True, capture_output=True, check=True
        )
        payload = json.loads(result.stdout)
        self.assertEqual(payload["sourceTime"], 0)
        self.assertEqual(payload["renderedTime"], -12.5)
        self.assertTrue(payload["motionIsOwned"])

    @unittest.skipIf(shutil.which("node") is None, "node is not available")
    def test_exact_frame_birth_visibility_uses_motion_time_when_override_is_missing(self):
        html = _runtime_html()
        helper_source = (
            "function animatedPointState(point, trace = null, timeOverride = null)"
            + html.split(
                "function animatedPointState(point, trace = null, timeOverride = null)", 1
            )[1].split("function pointHoverText", 1)[0]
        )
        script = f"""
        let fadeOpacityByBirthTimeEnabled = false;
        let fadeInTimeMyr = 8;
        let fadeInAndOutEnabled = false;
        function pointSizeForTrace() {{ return 10; }}
        function pointOpacityForTrace() {{ return 0.75; }}
        function birthOpacityFadePointSize() {{ return 10; }}
        function fadeVisibilityFactor(time) {{ return time === -20 ? 0.25 : (time === -10 ? 0.5 : 1); }}
        {helper_source}
        const point = {{ motion: {{ time_myr: -20 }} }};
        const implicit = animatedPointState(point, {{}});
        const explicit = animatedPointState(point, {{}}, -10);
        process.stdout.write(JSON.stringify({{ implicit, explicit }}));
        """
        result = subprocess.run(
            ["node"],
            input=script,
            text=True,
            capture_output=True,
            check=True,
        )
        payload = json.loads(result.stdout)
        self.assertEqual(payload["implicit"], {"size": 2.5, "opacity": 0.75})
        self.assertEqual(payload["explicit"], {"size": 5, "opacity": 0.75})

    @unittest.skipIf(shutil.which("node") is None, "node is not available")
    def test_birth_time_presence_reaches_zero_before_exact_frame(self):
        html = _runtime_html()
        helper_source = (
            "function ovizPointBirthVisibility(pointState, point, trace)"
            + html.split("function ovizPointBirthVisibility(pointState, point, trace)", 1)[1].split(
                "function ovizRetainedPointVisual", 1
            )[0]
        )
        script = f"""
        let fadeOpacityByBirthTimeEnabled = false;
        function clamp01(value) {{ return Math.min(Math.max(Number(value) || 0, 0), 1); }}
        function pointSizeForTrace(point) {{ return Number(point.baseSize); }}
        {helper_source}
        const point = {{ baseSize: 10 }};
        const sizeFade = [10, 5, 0].map((size) =>
          ovizPointBirthVisibility({{ size }}, point, {{}})
        );
        fadeOpacityByBirthTimeEnabled = true;
        const opacityFade = ovizPointBirthVisibility({{ size: 0 }}, point, {{}});
        process.stdout.write(JSON.stringify({{ sizeFade, opacityFade }}));
        """
        result = subprocess.run(
            ["node"],
            input=script,
            text=True,
            capture_output=True,
            check=True,
        )
        payload = json.loads(result.stdout)
        self.assertEqual(payload["sizeFade"], [1, 0.5, 0])
        self.assertEqual(payload["opacityFade"], 1)

    @unittest.skipIf(shutil.which("node") is None, "node is not available")
    def test_fractional_points_preserve_omitted_size_and_opacity_defaults(self):
        html = _runtime_html()
        interpolation_source = (
            "function interpolateNumber(fromValue, toValue, alpha, fallbackValue = 0.0)"
            + html.split(
                "function interpolateNumber(fromValue, toValue, alpha, fallbackValue = 0.0)", 1
            )[1].split("function interpolateTraceLabel", 1)[0]
        )
        script = f"""
        function clampRange(value, minimum, maximum) {{
          return Math.min(Math.max(Number(value) || 0, minimum), maximum);
        }}
        function clusterFilterSelectionKeyForPoint() {{ return ""; }}
        {interpolation_source}
        const pointA = {{
          x: 0, y: 2, z: 4,
          motion: {{ key: "cluster", age_now_myr: 10, time_myr: -1 }},
        }};
        const pointB = {{
          x: 2, y: 4, z: 6,
          motion: {{ key: "cluster", age_now_myr: 10, time_myr: 0 }},
        }};
        const blended = interpolateTracePoint(pointA, pointB, 0.5, -0.5);
        const explicit = interpolateTracePoint(
          {{ ...pointA, size: 4, opacity: 0.2 }},
          {{ ...pointB, size: 8, opacity: 0.6 }},
          0.5,
          -0.5
        );
        process.stdout.write(JSON.stringify({{
          blended,
          explicit,
          resolvedSize: blended.size ?? 7,
          resolvedOpacity: blended.opacity ?? 0.6,
        }}));
        """
        result = subprocess.run(
            ["node"],
            input=script,
            text=True,
            capture_output=True,
            check=True,
        )
        payload = json.loads(result.stdout)
        self.assertNotIn("size", payload["blended"])
        self.assertNotIn("opacity", payload["blended"])
        self.assertEqual(payload["resolvedSize"], 7)
        self.assertEqual(payload["resolvedOpacity"], 0.6)
        self.assertEqual(payload["explicit"]["size"], 6)
        self.assertAlmostEqual(payload["explicit"]["opacity"], 0.4)
        self.assertEqual(payload["blended"]["motion"]["time_myr"], -0.5)

    def test_transition_updates_do_not_rebuild_the_plot_group(self):
        html = _runtime_html()
        update_body = html.split(
            "function ovizUpdateRetainedTransitionScene(", 1
        )[1].split("function renderInterpolatedFrameValue(", 1)[0]
        prepare_body = html.split(
            "function ovizPrepareRetainedTransitionScene(", 1
        )[1].split("function ovizUpdateRetainedTransitionScene(", 1)[0]

        self.assertNotIn("clearGroup(plotGroup)", update_body)
        self.assertNotIn("clearGroup(runtime.overlayRoot)", update_body)
        self.assertNotIn("ovizRefreshRetainedSelectionOverlay", html)
        self.assertIn("ovizApplyRetainedSelectionOverlay(runtime, alpha)", update_body)
        self.assertIn("actionHeldAppearanceRollback", html)
        self.assertIn("sourceAppearanceComposite", html)
        self.assertIn("retainedOverlayLineMaterials", html)
        self.assertIn("actualOpacityByPoint", html)
        self.assertIn("rendered_selection_extra", html)
        self.assertIn("expectedVisiblePointCount", html)
        self.assertIn("clearGroup(plotGroup)", prepare_body)
        self.assertIn("intervalKey", update_body)
        self.assertIn("retainedSceneMetrics", update_body)
        self.assertIn("* birthVisibility", html)
        self.assertIn("* birthVisibility;", html)
        self.assertIn("ovizPointBirthVisibility(pointState, point, trace)", html)
        self.assertIn("presenceOpacity * birthVisibility", html)
        self.assertIn("material.opacity > 0.001", html)
        self.assertIn("expectedPointOpacity > 0.001", html)

    def test_galactic_references_preserve_retained_endpoint_crossfade(self):
        html = _runtime_html()
        endpoint_body = html.split(
            "function ovizApplyRetainedEndpointWeight(", 1
        )[1].split("function ovizPrepareRetainedSelectionOverlay(", 1)[0]
        timeline_body = html.split(
            "function updateGalacticReferenceMotionOpacity()", 1
        )[1].split("function updateMilkyWayTimelineOpacity()", 1)[0]

        self.assertIn(
            "material.userData.ovizRetainedEndpointOpacityScale = endpointWeight",
            endpoint_body,
        )
        self.assertIn("ovizRetainedEndpointOpacityScale", timeline_body)
        self.assertIn("? clamp01(retainedEndpointOpacity)", timeline_body)

    @unittest.skipIf(shutil.which("node") is None, "node is not available")
    def test_galactic_reference_refresh_does_not_reveal_both_retained_endpoints(self):
        html = _runtime_html()
        helper_source = (
            "function updateGalacticReferenceMotionOpacity()"
            + html.split(
                "function updateGalacticReferenceMotionOpacity()", 1
            )[1].split("function updateMilkyWayTimelineOpacity()", 1)[0]
        )
        script = f"""
        function clamp01(value) {{
          return Math.min(Math.max(Number(value) || 0, 0), 1);
        }}
        function galacticReferenceMotionVisible() {{ return true; }}
        function galacticPresentDayGridOpacity() {{ return 1; }}
        const root = {{ dataset: {{}} }};
        const materials = [0.7, 0.3].map((endpointWeight) => ({{
          opacity: 0,
          transparent: false,
          userData: {{
            ovizTimelineBaseOpacity: 0.5,
            ovizRetainedEndpointOpacityScale: endpointWeight,
          }},
        }}));
        const galacticReferenceOpacityGroups = materials.map((material) => ({{
          visible: false,
          userData: {{ ovizGalacticPresentDayReference: false }},
          traverse(callback) {{ callback({{ material }}); }},
        }}));
        {helper_source}
        updateGalacticReferenceMotionOpacity();
        process.stdout.write(JSON.stringify(materials.map((material) => material.opacity)));
        """
        result = subprocess.run(
            ["node"],
            input=script,
            text=True,
            capture_output=True,
            check=True,
        )
        self.assertEqual(json.loads(result.stdout), [0.35, 0.15])

    @unittest.skipIf(shutil.which("node") is None, "node is not available")
    def test_galactic_quadrant_grid_fades_away_from_present_day(self):
        html = _runtime_html()
        helper_source = (
            "function presentDayOnlyTimelineOpacity()"
            + html.split(
                "function presentDayOnlyTimelineOpacity()", 1
            )[1].split("function milkyWayTimelineOpacity()", 1)[0]
        )
        script = f"""
        function clampRange(value, minimum, maximum) {{
          return Math.min(Math.max(Number(value) || 0, minimum), maximum);
        }}
        function smoothstep01(value) {{
          const x = clampRange(value, 0, 1);
          return x * x * (3 - (2 * x));
        }}
        function frameTimeForValue(value) {{ return value; }}
        let displayedFrameValue = 0;
        {helper_source}
        const values = [0, 1, -0.5, -1, -10].map((time) => {{
          displayedFrameValue = time;
          return galacticPresentDayGridOpacity();
        }});
        process.stdout.write(JSON.stringify(values));
        """
        result = subprocess.run(
            ["node"],
            input=script,
            text=True,
            capture_output=True,
            check=True,
        )
        self.assertEqual(json.loads(result.stdout), [1, 0, 0.5, 0, 0])

    @unittest.skipIf(shutil.which("node") is None, "node is not available")
    def test_milky_way_uses_the_same_present_day_fade_as_longitude_markers(self):
        html = _runtime_html()
        helper_source = (
            "function presentDayOnlyTimelineOpacity()"
            + html.split(
                "function presentDayOnlyTimelineOpacity()", 1
            )[1].split("function setSkyDomeTimelineOpacityScale", 1)[0]
        )
        script = f"""
        function clampRange(value, minimum, maximum) {{
          return Math.min(Math.max(Number(value) || 0, minimum), maximum);
        }}
        function smoothstep01(value) {{
          const x = clampRange(value, 0, 1);
          return x * x * (3 - (2 * x));
        }}
        function frameTimeForValue(value) {{ return value; }}
        let displayedFrameValue = 0;
        {helper_source}
        const values = [0, -0.25, -0.5, -0.75, -1].map((time) => {{
          displayedFrameValue = time;
          return [galacticPresentDayGridOpacity(), milkyWayTimelineOpacity()];
        }});
        process.stdout.write(JSON.stringify(values));
        """
        result = subprocess.run(
            ["node"],
            input=script,
            text=True,
            capture_output=True,
            check=True,
        )
        values = json.loads(result.stdout)
        self.assertEqual([pair[0] for pair in values], [pair[1] for pair in values])
        self.assertEqual(values[0], [1, 1])
        self.assertEqual(values[-1], [0, 0])

    @unittest.skipIf(shutil.which("node") is None, "node is not available")
    def test_galactic_circle_geometry_interpolates_in_place_without_crossfade(self):
        html = _runtime_html()
        helper_source = (
            "function ovizApplyRetainedGalacticCircleLinePair("
            + html.split(
                "function ovizApplyRetainedGalacticCircleLinePair(", 1
            )[1].split("function ovizRetainedCloneLivePoint(", 1)[0]
        )
        script = f"""
        function clamp01(value) {{
          return Math.min(Math.max(Number(value) || 0, 0), 1);
        }}
        function interpolateNumber(from, to, alpha, fallback) {{
          const a = Number(from);
          const b = Number(to);
          return Number.isFinite(a) && Number.isFinite(b)
            ? a + ((b - a) * alpha)
            : fallback;
        }}
        function traceStyleStateForKey() {{ return null; }}
        function traceVisible() {{ return true; }}
        function traceVisibilityOpacityMultiplier() {{ return 1; }}
        function galacticReferenceMotionVisible() {{ return true; }}
        const buffer = {{ array: new Float32Array(6), needsUpdate: false }};
        const fromMaterial = {{
          opacity: 0,
          transparent: false,
          color: {{ set(value) {{ this.value = value; }} }},
          userData: {{}},
        }};
        const toMaterial = {{ opacity: 0.2, userData: {{}} }};
        const pair = {{
          from: {{
            object: {{ visible: false }},
            trace: {{
              key: "ring",
              opacity: 0.34,
              oviz_presence_opacity: 1,
              line: {{ color: "#94a3b8" }},
            }},
            materials: [fromMaterial],
          }},
          to: {{
            object: {{ visible: true }},
            trace: {{
              key: "ring",
              opacity: 0.34,
              oviz_presence_opacity: 1,
              line: {{ color: "#94a3b8" }},
            }},
            materials: [toMaterial],
          }},
          fromSegments: [[0, 0, 0, 2, 2, 2]],
          toSegments: [[4, 8, 12, 6, 10, 14]],
          positionBuffer: buffer,
        }};
        {helper_source}
        ovizApplyRetainedGalacticCircleLinePair(pair, 0.25);
        process.stdout.write(JSON.stringify({{
          positions: Array.from(buffer.array),
          needsUpdate: buffer.needsUpdate,
          fromVisible: pair.from.object.visible,
          toVisible: pair.to.object.visible,
          fromOpacity: fromMaterial.opacity,
          toOpacity: toMaterial.opacity,
        }}));
        """
        result = subprocess.run(
            ["node"],
            input=script,
            text=True,
            capture_output=True,
            check=True,
        )
        self.assertEqual(
            json.loads(result.stdout),
            {
                "positions": [1, 2, 3, 3, 4, 5],
                "needsUpdate": True,
                "fromVisible": True,
                "toVisible": False,
                "fromOpacity": 0.34,
                "toOpacity": 0,
            },
        )

    @unittest.skipIf(shutil.which("node") is None, "node is not available")
    def test_galactic_labels_interpolate_as_one_sprite_without_doubling(self):
        html = _runtime_html()
        helper_source = (
            "function ovizApplyRetainedGalacticLabelPair("
            + html.split(
                "function ovizApplyRetainedGalacticLabelPair(", 1
            )[1].split("function ovizRetainedCloneLivePoint(", 1)[0]
        )
        script = f"""
        function clamp01(value) {{
          return Math.min(Math.max(Number(value) || 0, 0), 1);
        }}
        function interpolateNumber(from, to, alpha, fallback) {{
          const a = Number(from);
          const b = Number(to);
          return Number.isFinite(a) && Number.isFinite(b)
            ? a + ((b - a) * alpha)
            : fallback;
        }}
        function traceStyleStateForKey() {{ return null; }}
        function traceVisible() {{ return true; }}
        function traceVisibilityOpacityMultiplier() {{ return 1; }}
        function isGalacticPresentDayReferenceTrace() {{ return false; }}
        function galacticPresentDayGridOpacity() {{ return 1; }}
        function galacticReferenceMotionVisible() {{ return true; }}
        const theme = {{ axis_color: "#ffffff" }};
        const fromParent = {{ visible: false }};
        const toParent = {{ visible: true }};
        const fromObject = {{
          visible: false,
          parent: fromParent,
          position: {{ set(x, y, z) {{ this.x = x; this.y = y; this.z = z; }} }},
        }};
        const toObject = {{ visible: true, parent: toParent }};
        const fromMaterial = {{
          opacity: 0,
          transparent: false,
          color: {{ set(value) {{ this.value = value; }} }},
          userData: {{}},
        }};
        const toMaterial = {{ opacity: 0.5, userData: {{}} }};
        const trace = {{
          key: "gc-label",
          oviz_presence_opacity: 1,
        }};
        const pair = {{
          from: {{
            object: fromObject,
            materials: [fromMaterial],
            labelMetadata: {{ trace }},
          }},
          to: {{
            object: toObject,
            materials: [toMaterial],
            labelMetadata: {{ trace }},
          }},
          fromLabel: {{
            text: "GALACTIC CENTER",
            x: 0,
            y: 10,
            z: 20,
            color: "#94a3b8",
            oviz_presence_opacity: 1,
          }},
          toLabel: {{
            text: "GALACTIC CENTER",
            x: 4,
            y: 18,
            z: 28,
            color: "#94a3b8",
            oviz_presence_opacity: 1,
          }},
        }};
        {helper_source}
        ovizApplyRetainedGalacticLabelPair(pair, 0.25);
        process.stdout.write(JSON.stringify({{
          position: [fromObject.position.x, fromObject.position.y, fromObject.position.z],
          fromVisible: fromObject.visible,
          toVisible: toObject.visible,
          fromParentVisible: fromParent.visible,
          toParentVisible: toParent.visible,
          fromOpacity: fromMaterial.opacity,
          toOpacity: toMaterial.opacity,
        }}));
        """
        result = subprocess.run(
            ["node"],
            input=script,
            text=True,
            capture_output=True,
            check=True,
        )
        self.assertEqual(
            json.loads(result.stdout),
            {
                "position": [1, 12, 22],
                "fromVisible": True,
                "toVisible": False,
                "fromParentVisible": True,
                "toParentVisible": False,
                "fromOpacity": 1,
                "toOpacity": 0,
            },
        )

    @unittest.skipIf(shutil.which("node") is None, "node is not available")
    def test_bottom_controls_stack_and_clear_the_timeline_on_collision(self):
        html = _runtime_html()
        helper_source = (
            "function ovizRectsApproach("
            + html.split("function ovizRectsApproach(", 1)[1].split(
                "function resize()", 1
            )[0]
        )
        script = f"""
        const root = {{
          clientWidth: 900,
          clientHeight: 700,
          dataset: {{}},
          getBoundingClientRect() {{
            return {{ left: 0, right: this.clientWidth, top: 0, bottom: this.clientHeight }};
          }},
        }};
        const footerEl = {{
          offsetParent: {{}},
          getBoundingClientRect() {{
            return root.clientWidth < 1000
              ? {{ left: 80, right: 820, top: 640, bottom: 688 }}
              : {{ left: 420, right: 1180, top: 840, bottom: 888 }};
          }},
        }};
        const bottomSwitchesEl = {{
          offsetParent: {{}},
          dataset: {{}},
          style: {{}},
          getBoundingClientRect() {{
            if (this.dataset.layout === "stacked") {{
              return this.style.bottom
                ? {{ left: 760, right: 886, top: 500, bottom: 616 }}
                : {{ left: 760, right: 886, top: 530, bottom: 632 }};
            }}
            return root.clientWidth < 1000
              ? {{ left: 600, right: 860, top: 642, bottom: 676 }}
              : {{ left: 1220, right: 1528, top: 848, bottom: 882 }};
          }},
        }};
        const timelineEnabled = true;
        {helper_source}
        updateBottomControlLayout();
        const compact = {{
          layout: bottomSwitchesEl.dataset.layout,
          clear: root.dataset.bottomControlsTimelineClear,
          bottom: bottomSwitchesEl.style.bottom,
        }};
        root.clientWidth = 1600;
        root.clientHeight = 900;
        updateBottomControlLayout();
        const wide = {{
          layout: bottomSwitchesEl.dataset.layout,
          clear: root.dataset.bottomControlsTimelineClear,
          bottom: bottomSwitchesEl.style.bottom,
        }};
        process.stdout.write(JSON.stringify({{ compact, wide }}));
        """
        result = subprocess.run(
            ["node"],
            input=script,
            text=True,
            capture_output=True,
            check=True,
        )
        self.assertEqual(
            json.loads(result.stdout),
            {
                "compact": {
                    "layout": "stacked",
                    "clear": "true",
                    "bottom": "72px",
                },
                "wide": {
                    "layout": "inline",
                    "clear": "true",
                    "bottom": "",
                },
            },
        )

    @unittest.skipIf(shutil.which("node") is None, "node is not available")
    def test_member_star_reveal_is_forced_off_in_3d(self):
        html = _runtime_html()
        helper_source = (
            "function setSkyMemberRevealProgress("
            + html.split("function setSkyMemberRevealProgress(", 1)[1].split(
                "function animateSkyMemberReveal(", 1
            )[0]
        )
        script = f"""
        function clampRange(value, minimum, maximum) {{
          return Math.min(Math.max(Number(value) || 0, minimum), maximum);
        }}
        const cameraViewMode = "free";
        let skyMemberRevealProgress = 1.0;
        let skyMemberBatchesEnabled = true;
        const starMaterial = {{ opacity: 1.0 }};
        const bulkMaterial = {{ opacity: 0.0 }};
        const skyMemberBatchOpacityEntries = [{{
          material: starMaterial,
          baseOpacity: 0.8,
          overflowActive: false,
        }}];
        const skyMemberBulkOpacityEntries = [{{
          material: bulkMaterial,
          baseOpacity: 0.6,
        }}];
        const root = {{ dataset: {{}} }};
        {helper_source}
        setSkyMemberRevealProgress(1.0);
        process.stdout.write(JSON.stringify({{
          reveal: skyMemberRevealProgress,
          batchesEnabled: skyMemberBatchesEnabled,
          starOpacity: starMaterial.opacity,
          bulkOpacity: bulkMaterial.opacity,
          visibleIn3d: root.dataset.skyMemberStarsVisibleIn3d,
        }}));
        """
        result = subprocess.run(
            ["node"],
            input=script,
            text=True,
            capture_output=True,
            check=True,
        )
        self.assertEqual(
            json.loads(result.stdout),
            {
                "reveal": 0,
                "batchesEnabled": False,
                "starOpacity": 0,
                "bulkOpacity": 0.6,
                "visibleIn3d": "false",
            },
        )

    def test_trace_and_point_interpolation_use_stable_keys(self):
        html = _runtime_html()
        interpolation_body = html.split(
            "function interpolatedFrameSpecForValue(", 1
        )[1].split("function makeVectorObject(", 1)[0]

        self.assertIn("upperTraceByKey", interpolation_body)
        self.assertIn("cloneTraceWithPresence", interpolation_body)
        self.assertNotIn("index < upperTraces.length", interpolation_body)
        self.assertIn("ovizTransitionStablePointKey", html)
        self.assertIn("motion:", html)
        self.assertIn("selection:", html)

    def test_retained_endpoints_keep_hidden_items_resident(self):
        html = _runtime_html()
        render_body = html.split("function renderFrameScene(", 1)[1].split(
            "let ovizRetainedTransitionScene", 1
        )[0]
        marker_body = html.split("function addMarkerTrace(parent, trace)", 1)[1].split(
            "function addTextTrace(parent, trace)", 1
        )[0]

        self.assertIn("forceResident", render_body)
        self.assertIn("addMarkerTrace(traceParent, trace, { forceResident })", render_body)
        self.assertIn("effectiveOpacity <= 0.001 && !forceResident", marker_body)
        self.assertIn("ovizRetainedPoint", marker_body)

    def test_retained_transition_filters_nonresident_traces_without_touching_frame_renderer(self):
        html = _runtime_html()
        prepare_body = html.split(
            "function ovizPrepareRetainedTransitionScene(", 1
        )[1].split("function ovizRetainedDebugSnapshot", 1)[0]

        self.assertIn("ovizRetainedTraceKeySet(sourceFromFrame, sourceToFrame)", prepare_body)
        self.assertIn("ovizRetainedFrameWithResidentTraces", prepare_body)
        self.assertIn("residentTraceCount", prepare_body)
        self.assertIn("livePoint: cloneTracePoint(fromEntry.metadata.point || {}) || {}", html)
        self.assertIn("const point = livePoint || metadata.point || {}", html)

    def test_adjacent_retained_intervals_reuse_the_shared_endpoint(self):
        html = _runtime_html()
        prepare_body = html.split(
            "function ovizPrepareRetainedTransitionScene(", 1
        )[1].split("function ovizRetainedDebugSnapshot", 1)[0]

        self.assertIn("previousRuntime.toEndpoint.frameIndex === frameState.lowerIndex", prepare_body)
        self.assertIn("previousRuntime.fromEndpoint.frameIndex === frameState.upperIndex", prepare_body)
        self.assertIn("ovizDisposeRetainedRoot", prepare_body)
        self.assertIn("ovizRetainedEndpointReuseSerial += 1", prepare_body)
        self.assertIn("ovizBuildRetainedEndpoint", prepare_body)
        self.assertNotIn("clearGroup(plotGroup);\n        plotGroup.position", prepare_body)
        self.assertIn('transitionOwnerToken: "timeline-playback"', html)
        self.assertIn('transitionOwnerToken: "timeline-scrub"', html)
        self.assertIn(
            "ovizPruneRetainedArray(milkyWayModelGroups, (entry) => entry, rootGroup)",
            html,
        )

    def test_retained_hot_loop_shares_point_evaluations_and_avoids_traversal(self):
        html = _runtime_html()
        update_body = html.split(
            "function ovizUpdateRetainedTransitionScene(", 1
        )[1].split("function renderInterpolatedFrameValue", 1)[0]
        endpoint_body = html.split(
            "function ovizApplyRetainedEndpointWeight(", 1
        )[1].split("function ovizPrepareRetainedSelectionOverlay", 1)[0]

        self.assertIn("runtime.livePointByIdentity.clear()", update_body)
        self.assertIn("runtime.commonVisualByEndpointIdentity.clear()", update_body)
        self.assertIn("ovizRetainedCommonPointVisual", update_body)
        self.assertIn("pointEvaluations", update_body)
        self.assertNotIn("updateCameraResponsiveImagePlanes()", update_body)
        self.assertIn("endpoint.nonPointEntries.forEach", endpoint_body)
        self.assertNotIn("endpoint.root.traverse", endpoint_body)
        self.assertIn("now - runtime.lastDiagnosticsAt < 100.0", update_body)

    def test_retained_point_visual_cache_is_separate_for_each_endpoint(self):
        html = _runtime_html()
        update_body = html.split(
            "function ovizUpdateRetainedTransitionScene(", 1
        )[1].split("function renderInterpolatedFrameValue", 1)[0]

        self.assertIn("const fromVisualIdentity = logicalIdentity", update_body)
        self.assertIn("const toVisualIdentity = logicalIdentity", update_body)
        self.assertIn(
            "ovizRetainedCommonPointVisual(pair.from, livePoint, displayedTimeMyr)",
            update_body,
        )
        self.assertIn(
            "ovizRetainedCommonPointVisual(pair.to, livePoint, displayedTimeMyr)",
            update_body,
        )
        self.assertNotIn("let common = logicalIdentity", update_body)

    def test_cached_textures_are_not_disposed_with_retained_materials(self):
        html = _runtime_html()
        clear_body = html.split("function clearGroup(group)", 1)[1].split(
            "function markerTextureFor(", 1
        )[0]
        text_body = html.split("function makeTextSprite(", 1)[1].split(
            "function updateScreenStableTextSprite(", 1
        )[0]

        self.assertIn("ovizSharedTexture", clear_body)
        self.assertIn("ovizSharedTexture", text_body)

    def test_manual_labels_and_selection_boxes_remain_resident(self):
        html = _runtime_html()
        prepare_body = html.split(
            "function ovizPrepareRetainedSelectionOverlay(", 1
        )[1].split("function ovizApplyRetainedSelectionOverlay(", 1)[0]
        apply_body = html.split(
            "function ovizApplyRetainedSelectionOverlay(", 1
        )[1].split("function ovizPrepareRetainedTransitionScene(", 1)[0]

        self.assertIn("createManualLabelEndpoint", prepare_body)
        self.assertIn("createSelectionBoxEndpoint", prepare_body)
        self.assertIn("fromLower", prepare_body)
        self.assertIn("toUpper", prepare_body)
        self.assertIn("overlay.manualFrom", apply_body)
        self.assertIn("overlay.manualTo", apply_body)
        self.assertIn("boxFromWeight * lowerWeight", apply_body)
        self.assertIn("boxToWeight * upperWeight", apply_body)
        self.assertIn("boxCompositeWeight * lowerWeight", apply_body)


if __name__ == "__main__":
    unittest.main()
