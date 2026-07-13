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
        self.assertIn(
            "addMarkerTrace(traceParent, trace, { forceResident, retainedPointComponents })",
            render_body,
        )
        self.assertIn("effectiveOpacity <= 0.001 && !forceResident", marker_body)
        self.assertIn("ovizRetainedPoint", marker_body)

    def test_retained_scene_filters_hidden_traces_and_unused_point_components(self):
        html = _runtime_html()
        prepare_body = html.split(
            "function ovizPrepareRetainedTransitionScene(", 1
        )[1].split("function ovizRetainedDebugSnapshot", 1)[0]
        marker_body = html.split("function addMarkerTrace(parent, trace)", 1)[1].split(
            "function addTextTrace(parent, trace)", 1
        )[0]

        self.assertIn("ovizRetainedTraceKeySet(fromFrame, toFrame)", prepare_body)
        self.assertIn("residentTraceKeys", prepare_body)
        self.assertIn("ovizRetainedPointComponentPlan()", prepare_body)
        self.assertIn("retainGlowComponents", marker_body)
        self.assertIn("retainMarkerComponent", marker_body)
        self.assertIn("residentTraceCount", html)
        self.assertIn("runtime.livePointByIdentity.clear()", html)
        self.assertIn("runtime.commonVisualByIdentity.clear()", html)
        self.assertIn("ovizRetainedCommonPointVisual", html)
        self.assertIn("pointEvaluations", html)
        self.assertIn("ageKdeLastRenderSignature", html)
        self.assertIn("ovizRetainedEndpointReuseSerial", html)
        self.assertIn("ovizRetainedSetsEqual", prepare_body)
        self.assertIn("previousRuntime.toEndpoint.frameIndex === frameState.lowerIndex", prepare_body)
        self.assertIn("ovizDisposeRetainedRoot", prepare_body)
        self.assertIn("ovizResetSceneRegistriesFromRoots", prepare_body)

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
