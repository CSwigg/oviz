from __future__ import annotations

import json
import re
import unittest

from oviz.threejs_figure import ThreeJSFigure
from oviz.threejs_actions import normalize_threejs_actions
from oviz.threejs_states import (
    DEFAULT_STATE_TRANSITION,
    deduplicate_state_assets,
    normalize_states_spec,
    normalize_transition,
)
from oviz.threejs_scene import _round_scene_floats


class ThreeJSStatesSchemaTests(unittest.TestCase):
    def test_action_can_delegate_exactly_one_step_to_a_viewer_state(self):
        actions = normalize_threejs_actions(
            [{"key": "sky", "label": "Sky", "steps": [{"type": "state", "state": "state-sky"}]}],
            group_names=["All"],
            trace_key_by_name={},
            playback_interval_ms=160,
        )

        self.assertEqual(actions[0]["steps"], [{
            "type": "state",
            "start": "after_previous",
            "delay_ms": 0,
            "state": "state-sky",
        }])
        with self.assertRaisesRegex(ValueError, "exactly one state step"):
            normalize_threejs_actions(
                [{
                    "key": "bad",
                    "label": "Bad",
                    "steps": [
                        {"type": "state", "state": "state-sky"},
                        {"type": "time", "direction": "backward"},
                    ],
                }],
                group_names=["All"],
                trace_key_by_name={},
                playback_interval_ms=160,
            )
        with self.assertRaisesRegex(ValueError, "cannot overlap two 'camera' steps"):
            normalize_threejs_actions(
                [{
                    "key": "overlap",
                    "label": "Overlap",
                    "steps": [
                        {"type": "camera", "target": {"kind": "point", "x": 0, "y": 0, "z": 0}},
                        {
                            "type": "camera",
                            "start": "with_previous",
                            "target": {"kind": "point", "x": 1, "y": 1, "z": 1},
                        },
                    ],
                }],
                group_names=["All"],
                trace_key_by_name={},
                playback_interval_ms=160,
            )

    def test_compact_rounding_preserves_control_range_precision(self):
        rounded = _round_scene_floats({"x": 1.234, "vmax": 0.07, "cut_max": 0.07123}, 1)

        self.assertEqual(rounded["x"], 1.2)
        self.assertEqual(rounded["vmax"], 0.07)
        self.assertEqual(rounded["cut_max"], 0.0712)

    def test_missing_states_gets_stable_project_schema(self):
        scene = {"width": 640, "height": 480, "frames": [], "initial_state": {}}
        figure = ThreeJSFigure(scene)

        first = figure.to_dict()["states"]
        second = figure.to_dict()["states"]

        self.assertEqual(first["project_id"], second["project_id"])
        self.assertEqual(first["schema_version"], 1)
        self.assertEqual(first["default_mode"], "edit")
        self.assertEqual(first["default_transition"], DEFAULT_STATE_TRANSITION)
        self.assertEqual(first["items"], [])

        same_scene = ThreeJSFigure({"width": 640, "height": 480, "frames": [], "initial_state": {}})
        self.assertEqual(first["project_id"], same_scene.to_dict()["states"]["project_id"])

    def test_normalization_preserves_order_and_repairs_duplicate_ids(self):
        states = normalize_states_spec({
            "project_id": "project-fixed",
            "default_mode": "present",
            "items": [
                {"id": "same", "name": "First", "snapshot": {"current_frame_index": 1}},
                {"id": "same", "name": "Second", "snapshot": {"current_frame_index": 2}},
            ],
        })

        self.assertEqual([item["name"] for item in states["items"]], ["First", "Second"])
        self.assertEqual(states["items"][0]["id"], "same")
        self.assertNotEqual(states["items"][1]["id"], "same")
        self.assertEqual(states["project_id"], "project-fixed")
        self.assertEqual(states["default_mode"], "present")

    def test_target_transition_override_is_normalized_independently(self):
        states = normalize_states_spec({
            "default_transition": {"duration_ms": 1400, "easing": "linear"},
            "items": [
                {"name": "Default", "snapshot": {}},
                {"name": "Fast", "transition": {"duration_ms": 250, "easing": "easeOutCubic"}, "snapshot": {}},
            ],
        })

        self.assertIsNone(states["items"][0]["transition"])
        self.assertEqual(states["items"][1]["transition"], {"duration_ms": 250, "easing": "easeOutCubic"})
        self.assertEqual(normalize_transition({"duration_ms": -10})["duration_ms"], 0)

    def test_large_state_assets_are_content_addressed_and_deduplicated(self):
        data_url = "data:image/png;base64," + ("A" * 5000)
        compact, assets = deduplicate_state_assets({"one": data_url, "nested": [data_url]})

        first_ref = compact["one"]["__oviz_asset_ref__"]
        second_ref = compact["nested"][0]["__oviz_asset_ref__"]
        self.assertEqual(first_ref, second_ref)
        self.assertEqual(len(assets), 1)
        self.assertEqual(assets[first_ref], data_url)


class ThreeJSStatesRuntimeTests(unittest.TestCase):
    def test_runtime_and_non_dom_api_are_embedded(self):
        figure = ThreeJSFigure({"width": 640, "height": 480, "frames": [], "initial_state": {}})
        html = figure.to_html(compress_scene_spec=False)

        self.assertIn("const OVIZ_STATES_VERSION = 1", html)
        self.assertIn("window.Oviz.get =", html)
        self.assertIn('ovizStateEvent("states-ready"', html)
        self.assertIn('ovizStateEvent("transition-progress"', html)
        self.assertIn('data.source !== "oviz-command"', html)
        self.assertNotIn("function ovizSwapTransitionScene(", html)
        self.assertNotIn("const canvasOpacity =", html)
        self.assertNotIn("renderer.domElement.style.opacity = String(canvasOpacity)", html)
        self.assertIn("restoreSkyLayerStateFromSnapshot(initialState, {", html)
        self.assertIn("postToAladin: options.postSkyLayersToAladin !== false", html)
        self.assertIn("ovizStateCameraSignature", html)
        self.assertIn("lockEarthViewCameraToTarget()", html)
        self.assertIn("viewFromEarth();", html)
        self.assertIn("exitEarthView();", html)
        self.assertIn("transition.nativeViewTransition", html)
        self.assertIn("destinationCameraState", html)
        self.assertIn("preApplyFovError", html)
        self.assertIn("now - transition.lastProgressEventAt >= 100", html)
        self.assertIn("root.dataset.stateTransitionMetrics", html)
        self.assertIn("root.dataset.stateFidelity", html)
        self.assertIn("ovizStateFidelityDifferences", html)
        self.assertIn("ovizApplyCapturedCameraState", html)
        self.assertIn("ovizCancelStateTransitionWithoutSnap", html)
        self.assertIn("ovizStateTransition.targetIndex", html)
        self.assertIn("startStateAction", html)
        self.assertIn("lassoSelectionRestoreSerial", html)
        self.assertIn("restoreSerial !== lassoSelectionRestoreSerial", html)
        self.assertIn("previousStep.finished", html)
        self.assertIn("completeCamera: false", html)
        self.assertIn("ovizSkyLayerTransitionWaiters", html)
        self.assertIn('data.type === "oviz-aladin-sky-layer-transition-complete"', html)
        self.assertIn("stateOwnsCameraAndTime", html)
        self.assertIn("actionHeldTraceOpacityByKey", html)
        self.assertIn("captureActionTraceOpacityMap", html)
        self.assertIn("preserveLegendTransitionFrame: true", html)
        cancel_action_body = html.split(
            "function ovizCancelActionWithoutSnap()", 1
        )[1].split("function ovizCancelStateTransitionWithoutSnap", 1)[0]
        self.assertIn("actionHeldTraceOpacityByKey = null", cancel_action_body)
        self.assertNotIn("__STATE_RUNTIME_JS__", html)

    def test_state_transition_runtime_exposes_ordered_phase_instrumentation(self):
        html = ThreeJSFigure({
            "width": 640,
            "height": 480,
            "frames": [],
            "initial_state": {},
        }).to_html(compress_scene_spec=False)

        self.assertIn("phaseProgress", html)
        self.assertIn("effectiveDurationMs", html)
        self.assertIn("stateTransitionPhase", html)
        self.assertIn("stateTransitionPhaseProgress", html)
        self.assertIn("stateTransitionEffectiveDurationMs", html)
        self.assertRegex(html, r"(?i)phase[^\n]{0,80}(?:minimum|min)[^\n]{0,80}800")
        self.assertRegex(html, r'["\']camera["\']')
        self.assertRegex(html, r'["\']appearance["\']')
        self.assertRegex(html, r'["\']time["\']')

    def test_legacy_time_actions_advance_fractional_frames_each_animation_frame(self):
        html = ThreeJSFigure({
            "width": 640,
            "height": 480,
            "frames": [],
            "initial_state": {},
        }).to_html(compress_scene_spec=False)
        update_body = html.split("function updateTimeAction(now)", 1)[1].split(
            "function updateViewerActions(now)", 1
        )[0]

        self.assertIn("displayedFrameValue", update_body)
        self.assertIn("renderInterpolatedFrameValue", update_body)
        self.assertIn("preserveCamera: true", update_body)
        self.assertNotIn("currentFrameIndex + timeActionTrack.direction", update_body)
        self.assertNotIn("for (let idx = 0; idx < steps; idx += 1)", update_body)

    def test_appearance_phase_crossfades_trace_and_lasso_membership(self):
        html = ThreeJSFigure({
            "width": 640,
            "height": 480,
            "frames": [],
            "initial_state": {},
        }).to_html(compress_scene_spec=False)

        self.assertIn("function ovizTraceCandidatesForSnapshots(", html)
        self.assertIn("ovizStateSelectionTransition", html)
        self.assertIn("function ovizSelectionMembershipOpacity(", html)
        self.assertIn("applyLassoSelectionTransitionUniforms", html)
        self.assertIn("selectionTransitionProgress", html)
        self.assertIn("selectionSourceSecondaryMaskTexture", html)
        self.assertIn("float sourceSelectionWeight = mix(", html)
        self.assertIn("selectionWeight = mix(", html)
        self.assertIn("function ovizTransitionOpacityBucket(", html)
        opacity_bucket_body = html.split(
            "function ovizTransitionOpacityBucket(", 1
        )[1].split("function markerMaterialFor(", 1)[0]
        self.assertIn("256", opacity_bucket_body)
        self.assertIn("for (let index = firstIndex; index <= lastIndex", html)
        self.assertIn("1.0 - selectionProgress", html)
        marker_material_body = html.split("function markerMaterialFor(", 1)[1].split(
            "function smoothstep(", 1
        )[0]
        self.assertIn("ovizTransitionOpacityBucket(", marker_material_body)
        glow_material_body = html.split("function starGlowMaterialFor(", 1)[1].split(
            "function starCoreMaterialFor(", 1
        )[0]
        core_material_body = html.split("function starCoreMaterialFor(", 1)[1].split(
            "function glowScaleForPoint(", 1
        )[0]
        self.assertIn("ovizTransitionOpacityBucket(", glow_material_body)
        self.assertIn("ovizTransitionOpacityBucket(", core_material_body)

    def test_exact_target_render_drops_transient_visibility_overrides_first(self):
        html = ThreeJSFigure({
            "width": 640,
            "height": 480,
            "frames": [],
            "initial_state": {},
        }).to_html(compress_scene_spec=False)
        finish_body = html.split(
            "async function ovizFinishStateTransition(transition)", 1
        )[1].split("function ovizFailStateTransition", 1)[0]

        exact_apply = finish_body.index(
            "ovizApplyStateImmediately(transition.targetSnapshot"
        )
        self.assertLess(
            finish_body.index("ovizStateTransitionTraceOpacity = null"),
            exact_apply,
        )
        self.assertLess(
            finish_body.index("ovizStateSelectionTransition = null"),
            exact_apply,
        )

    def test_mask_only_lasso_states_keep_destination_points_visible(self):
        html = ThreeJSFigure({
            "width": 640,
            "height": 480,
            "frames": [],
            "initial_state": {},
        }).to_html(compress_scene_spec=False)
        load_target_mask = html.split(
            "transition.targetRuntimeLassoMask = mask", 1
        )[1].split("return mask", 1)[0]
        marker_body = html.split("function addMarkerTrace(parent, trace)", 1)[1].split(
            "function addTextTrace(parent, trace)", 1
        )[0]

        self.assertIn("transition.selectionTransition.toMask = mask", load_target_mask)
        self.assertNotIn("destination.lasso_volume_selection_enabled", load_target_mask)
        self.assertIn(
            'if (typeof ovizSelectionMembershipOpacity === "function")',
            marker_body,
        )

    def test_embedded_scene_state_schema_round_trips_to_payload(self):
        scene = {
            "width": 640,
            "height": 480,
            "frames": [],
            "initial_state": {"current_frame_index": 0},
            "states": {
                "project_id": "project-roundtrip",
                "default_mode": "present",
                "items": [{
                    "id": "state-one",
                    "name": "One",
                    "transition": {"duration_ms": 333, "easing": "linear"},
                    "snapshot": {"current_frame_value": 0.25},
                }],
            },
        }
        html = ThreeJSFigure(scene).to_html(compress_scene_spec=False)
        match = re.search(
            r"/\*__SCENE_SPEC_START__\*/const sceneSpec = (.*?);/\*__SCENE_SPEC_END__\*/",
            html,
            re.DOTALL,
        )
        self.assertIsNotNone(match)
        payload = json.loads(match.group(1))
        self.assertEqual(payload["states"]["project_id"], "project-roundtrip")
        self.assertEqual(payload["states"]["items"][0]["id"], "state-one")
        self.assertEqual(payload["states"]["items"][0]["snapshot"]["current_frame_value"], 0.25)


if __name__ == "__main__":
    unittest.main()
