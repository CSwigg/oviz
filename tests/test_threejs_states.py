from __future__ import annotations

import json
import re
import shutil
import subprocess
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
    def test_sky_member_display_mode_round_trips_through_states(self):
        html = ThreeJSFigure({
            "width": 640,
            "height": 480,
            "frames": [],
            "initial_state": {
                "global_controls": {
                    "camera_view_mode": "earth",
                    "sky_member_display_mode": "clusters",
                },
            },
        }).to_html(compress_scene_spec=False)

        self.assertIn(
            "sky_member_display_mode: skyMemberDisplayMode",
            html,
        )
        self.assertIn(
            "savedGlobalControls.sky_member_display_mode",
            html,
        )
        self.assertIn(
            "skyMemberDisplayTargetProgress(fromGlobal.sky_member_display_mode)",
            html,
        )
        self.assertIn(
            "skyMemberDisplayTargetProgress(toGlobal.sky_member_display_mode)",
            html,
        )
        self.assertIn(
            "(destination.global_controls || {}).sky_member_display_mode",
            html,
        )

    @unittest.skipIf(shutil.which("node") is None, "node is not available")
    def test_state_only_presentation_navigation_cycles_authored_states(self):
        html = ThreeJSFigure({
            "width": 640,
            "height": 480,
            "frames": [],
            "initial_state": {},
        }).to_html(compress_scene_spec=False)
        cycle_source = (
            "function ovizSyncPresentationStateNavigationDiagnostics()"
            + html.split("function ovizSyncPresentationStateNavigationDiagnostics()", 1)[1].split(
                "function ovizAssertEditable", 1
            )[0]
        )
        script = f"""
        let ovizStatesProject = {{ items: [{{ id: "state-1" }}, {{ id: "state-2" }}] }};
        let ovizActiveStateId = null;
        let ovizStateTransition = null;
        let ovizStateControllerReady = true;
        let ovizPresentationStateNavigationQueue = Promise.resolve({{ idle: true }});
        let ovizPresentationStateNavigationGeneration = 0;
        let ovizPresentationStateNavigationPending = 0;
        const root = {{ dataset: {{}} }};
        const visited = [];
        function ovizGoToState(target) {{
          if (target === "original") ovizActiveStateId = null;
          else ovizActiveStateId = ovizStatesProject.items[Number(target) - 1].id;
          visited.push(target);
          return Promise.resolve(target);
        }}
        {cycle_source}
        (async () => {{
          await ovizStatesPresentationNext();
          await ovizStatesPresentationNext();
          await ovizStatesPresentationNext();
          await ovizStatesPresentationPrevious();
          await ovizStatesPresentationPrevious();
          await ovizStatesPresentationPrevious();
          process.stdout.write(JSON.stringify(visited));
        }})().catch((error) => {{ throw error; }});
        """
        result = subprocess.run(
            ["node"],
            input=script,
            text=True,
            capture_output=True,
            check=True,
        )

        self.assertEqual(json.loads(result.stdout), [1, 2, 1, 2, 1, 2])

    @unittest.skipIf(shutil.which("node") is None, "node is not available")
    def test_state_only_presentation_navigation_waits_for_ready_and_serializes_rapid_taps(self):
        html = ThreeJSFigure({
            "width": 640,
            "height": 480,
            "frames": [],
            "initial_state": {},
        }).to_html(compress_scene_spec=False)
        cycle_source = (
            "function ovizSyncPresentationStateNavigationDiagnostics()"
            + html.split("function ovizSyncPresentationStateNavigationDiagnostics()", 1)[1].split(
                "function ovizAssertEditable", 1
            )[0]
        )
        script = f"""
        let ovizStatesProject = {{ items: [
          {{ id: "state-1" }}, {{ id: "state-2" }}, {{ id: "state-3" }}
        ] }};
        let ovizActiveStateId = null;
        let ovizStateTransition = null;
        let ovizStateControllerReady = false;
        let ovizPresentationStateNavigationQueue = Promise.resolve({{ idle: true }});
        let ovizPresentationStateNavigationGeneration = 0;
        let ovizPresentationStateNavigationPending = 0;
        const root = new EventTarget();
        root.dataset = {{}};
        const window = {{
          setTimeout: globalThis.setTimeout,
          clearTimeout: globalThis.clearTimeout,
        }};
        const visited = [];
        let inFlight = 0;
        let maxInFlight = 0;
        function ovizGoToState(target) {{
          inFlight += 1;
          maxInFlight = Math.max(maxInFlight, inFlight);
          visited.push(target);
          const promise = new Promise((resolve) => globalThis.setTimeout(() => {{
            ovizActiveStateId = ovizStatesProject.items[Number(target) - 1].id;
            ovizStateTransition = null;
            inFlight -= 1;
            resolve(target);
          }}, 5));
          ovizStateTransition = {{ targetIndex: Number(target) - 1, promise }};
          return promise;
        }}
        {cycle_source}
        (async () => {{
          const first = ovizStatesPresentationNext();
          const second = ovizStatesPresentationNext();
          const third = ovizStatesPresentationNext();
          globalThis.setTimeout(() => {{
            ovizStateControllerReady = true;
            root.dispatchEvent(new Event("states-ready"));
          }}, 5);
          await Promise.all([first, second, third]);
          process.stdout.write(JSON.stringify({{
            visited,
            maxInFlight,
            active: ovizActiveStateId,
            pending: ovizPresentationStateNavigationPending,
          }}));
        }})().catch((error) => {{ throw error; }});
        """
        result = subprocess.run(
            ["node"],
            input=script,
            text=True,
            capture_output=True,
            check=True,
        )
        payload = json.loads(result.stdout)

        self.assertEqual(payload["visited"], [1, 2, 3])
        self.assertEqual(payload["maxInFlight"], 1)
        self.assertEqual(payload["active"], "state-3")
        self.assertEqual(payload["pending"], 0)

    @unittest.skipIf(shutil.which("node") is None, "node is not available")
    def test_exact_state_restore_preserves_captured_legend_panel_geometry(self):
        html = ThreeJSFigure({
            "width": 640,
            "height": 480,
            "frames": [],
            "initial_state": {},
        }).to_html(compress_scene_spec=False)
        legend_source = (
            "function setLegendPanelOpen(isOpen, options = {})"
            + html.split("function setLegendPanelOpen(isOpen, options = {})", 1)[1].split(
                "function visibleLegendItemsForCurrentGroup", 1
            )[0]
        )
        script = f"""
        let minimalModeEnabled = false;
        let legendPanelOpen = false;
        let legendPanelRectState = {{ left: 14, top: 22, width: 220, height: 300 }};
        const legendPanelEl = {{ dataset: {{}}, style: {{}} }};
        const legendPanelBodyEl = {{ style: {{}} }};
        const legendPanelToggleEl = {{ setAttribute() {{}}, textContent: "" }};
        const mobileLegendButtonEl = {{ dataset: {{}}, setAttribute() {{}}, title: "" }};
        const legendResizeEls = [];
        let receivedOptions = null;
        function defaultLegendPanelRect() {{ return legendPanelRectState; }}
        function applyLegendPanelRect(_rect, options) {{ receivedOptions = options; }}
        function closeLegendPopover() {{}}
        {legend_source}
        setLegendPanelOpen(true, {{ allowAutoCap: false }});
        process.stdout.write(JSON.stringify(receivedOptions));
        """
        result = subprocess.run(
            ["node"],
            input=script,
            text=True,
            capture_output=True,
            check=True,
        )
        self.assertEqual(json.loads(result.stdout), {"allowAutoCap": False})
        self.assertIn("allowAutoCap: !hasSavedLegendPanelRect", html)
        exact_apply_source = html.split("function ovizApplyStateImmediately", 1)[1].split(
            "function ovizPointFrom", 1
        )[0]
        self.assertIn("applyLegendPanelRect(hydrated.legend_panel_state", exact_apply_source)
        self.assertIn("allowAutoCap: false", exact_apply_source)
        self.assertLess(
            exact_apply_source.index("resize();"),
            exact_apply_source.index("applyLegendPanelRect(hydrated.legend_panel_state"),
        )

    @unittest.skipIf(shutil.which("node") is None, "node is not available")
    def test_free_view_null_sky_camera_fields_restore_as_null(self):
        html = ThreeJSFigure({
            "width": 640,
            "height": 480,
            "frames": [],
            "initial_state": {},
        }).to_html(compress_scene_spec=False)
        restore_source = (
            "function restoreEarthViewDerivedState(savedGlobalControls)"
            + html.split("function restoreEarthViewDerivedState(savedGlobalControls)", 1)[1].split(
                "function applyViewerStateSyncInternal", 1
            )[0]
        )
        script = f"""
        let earthViewFocusDistance = 8122;
        let earthViewReturnCameraState = {{ stale: true }};
        function cameraReturnStateFromPlainObject(value) {{
          return value && typeof value === "object" ? value : null;
        }}
        {restore_source}
        restoreEarthViewDerivedState({{
          earth_view_focus_distance_pc: null,
          earth_view_return_camera_state: null,
        }});
        process.stdout.write(JSON.stringify({{
          focusDistance: earthViewFocusDistance,
          returnCamera: earthViewReturnCameraState,
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
        self.assertIsNone(payload["focusDistance"])
        self.assertIsNone(payload["returnCamera"])

    @unittest.skipIf(shutil.which("node") is None, "node is not available")
    def test_state_navigation_clears_cancelled_action_trace_hold(self):
        html = ThreeJSFigure({
            "width": 640,
            "height": 480,
            "frames": [],
            "initial_state": {},
        }).to_html(compress_scene_spec=False)
        cancel_source = (
            "function ovizCancelActionWithoutSnap(reason"
            + html.split("function ovizCancelActionWithoutSnap(reason", 1)[1].split(
                "function ovizCancelStateTransitionWithoutSnap", 1
            )[0]
        )
        script = f"""
        let activeActionRun = {{ id: "active" }};
        let actionHeldTraceOpacityByKey = {{ young: 0.0 }};
        function cancelActionRun() {{
          actionHeldTraceOpacityByKey = {{ young: 0.0 }};
          return true;
        }}
        {cancel_source}
        const result = ovizCancelActionWithoutSnap("state-navigation");
        process.stdout.write(JSON.stringify({{
          result,
          heldTraceOpacity: actionHeldTraceOpacityByKey,
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
        self.assertTrue(payload["result"])
        self.assertIsNone(payload["heldTraceOpacity"])

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
        self.assertFalse(first["present_only"])
        self.assertFalse(first["preserve_camera_on_navigation"])
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
        self.assertEqual(
            [item["camera_behavior"] for item in states["items"]],
            ["follow", "follow"],
        )

    def test_camera_preservation_policy_round_trips(self):
        states = normalize_states_spec({
            "preserve_camera_on_navigation": True,
            "items": [
                {"name": "Legacy default", "snapshot": {}},
                {"name": "Explicit follow", "camera_behavior": "follow", "snapshot": {}},
                {"name": "Explicit keep", "cameraBehavior": "keep", "snapshot": {}},
            ],
        })

        self.assertTrue(states["preserve_camera_on_navigation"])
        self.assertEqual(
            [item["camera_behavior"] for item in states["items"]],
            ["keep", "follow", "keep"],
        )

    def test_present_only_export_mode_round_trips(self):
        states = normalize_states_spec({"default_mode": "present", "present_only": True})

        self.assertEqual(states["default_mode"], "present")
        self.assertTrue(states["present_only"])

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
    @unittest.skipIf(shutil.which("node") is None, "node is not available")
    def test_keep_camera_policy_replaces_only_camera_domains(self):
        html = ThreeJSFigure({
            "width": 640,
            "height": 480,
            "frames": [],
            "initial_state": {},
        }).to_html(compress_scene_spec=False)
        helper_source = (
            "function ovizPreserveDestinationCamera(destination, source)"
            + html.split("function ovizPreserveDestinationCamera(destination, source)", 1)[1].split(
                "function ovizBeginStateTransition", 1
            )[0]
        )
        script = f"""
        function ovizStatesClone(value, fallback = null) {{
          try {{ return JSON.parse(JSON.stringify(value)); }} catch (_err) {{ return fallback; }}
        }}
        {helper_source}
        const destination = {{
          camera: {{ position: {{ x: 100, y: 200, z: 300 }} }},
          global_controls: {{
            camera_view_mode: "earth",
            camera_fov: 20,
            camera_auto_orbit_enabled: true,
            point_opacity_scale: 0.35,
          }},
          current_frame_value: 7.5,
        }};
        const source = {{
          camera: {{ position: {{ x: 1, y: 2, z: 3 }}, target: {{ x: 4, y: 5, z: 6 }} }},
          global_controls: {{
            camera_view_mode: "free",
            camera_fov: 61,
            camera_auto_orbit_enabled: false,
            camera_auto_orbit_speed: 0.2,
            earth_view_focus_distance_pc: 8122,
          }},
        }};
        const result = ovizPreserveDestinationCamera(destination, source);
        process.stdout.write(JSON.stringify(result));
        """
        result = subprocess.run(
            ["node"], input=script, text=True, capture_output=True, check=True
        )
        payload = json.loads(result.stdout)

        self.assertEqual(payload["camera"], {
            "position": {"x": 1, "y": 2, "z": 3},
            "target": {"x": 4, "y": 5, "z": 6},
        })
        self.assertEqual(payload["global_controls"]["camera_view_mode"], "free")
        self.assertEqual(payload["global_controls"]["camera_fov"], 61)
        self.assertFalse(payload["global_controls"]["camera_auto_orbit_enabled"])
        self.assertEqual(payload["global_controls"]["point_opacity_scale"], 0.35)
        self.assertEqual(payload["current_frame_value"], 7.5)

    def test_runtime_and_non_dom_api_are_embedded(self):
        figure = ThreeJSFigure({"width": 640, "height": 480, "frames": [], "initial_state": {}})
        html = figure.to_html(compress_scene_spec=False)

        self.assertIn("const OVIZ_STATES_VERSION = 1", html)
        self.assertIn("window.Oviz.get =", html)
        self.assertIn("quickAdd: ovizQuickAddState", html)
        self.assertIn("setPreserveCamera: ovizSetPreserveCameraOnNavigation", html)
        self.assertIn("setCameraBehavior: ovizSetStateCameraBehavior", html)
        self.assertIn("exportPresentOnlyHtml:", html)
        self.assertIn("function ovizPreserveDestinationCamera(destination, source)", html)
        self.assertIn("ovizPreserveDestinationCamera(destination, from);", html)
        self.assertIn('"camera_view_mode"', html)
        self.assertIn('"earth_view_return_camera_state"', html)
        self.assertIn("Cam: Keep", html)
        self.assertIn("Cam: Follow", html)
        self.assertIn('target.camera_behavior || target.cameraBehavior || ""', html)
        self.assertIn("row.dataset.cameraGroup = String(cameraGroupIndex)", html)
        self.assertIn("const allowLiveCameraOrbit = Boolean(", html)
        self.assertIn("if (!transition.preserveCamera) {", html)
        self.assertIn("transition.targetSnapshot,\n            ovizCurrentTransitionSnapshot()", html)
        self.assertIn("&& !activeStateTransition.finishing", html)
        self.assertIn("stateAllowsLiveCameraOrbit && cameraAutoOrbitEnabled", html)
        self.assertIn("updateGalacticSimpleDefaultOrbit(deltaSeconds);", html)
        self.assertIn('.oviz-three-topbar:has(.oviz-states-shell[data-open=true]){z-index:90 !important}', html)
        self.assertIn("states-preload-complete", html)
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
        self.assertIn('differences.push("transient_action_trace_opacity")', html)
        self.assertIn("ovizStateFidelityDifferences", html)
        self.assertIn("ovizApplyCapturedCameraState", html)
        self.assertIn("savedGlobalControls.earth_view_focus_distance_pc === null", html)
        self.assertIn(
            'Object.prototype.hasOwnProperty.call(savedGlobalControls, "earth_view_return_camera_state")',
            html,
        )
        self.assertIn("ovizCancelStateTransitionWithoutSnap", html)
        self.assertIn("ovizStateTransition.targetIndex", html)
        self.assertIn("startStateAction", html)
        self.assertIn("lassoSelectionRestoreSerial", html)
        self.assertIn("restoreSerial !== lassoSelectionRestoreSerial", html)
        self.assertIn("const OvizTransitionRuntime = (() =>", html)
        self.assertIn("OvizTransitionRuntime.createFrameCoordinator", html)
        self.assertIn("updateOvizUnifiedTransitionSession(now)", html)
        self.assertIn("previous.startedAt === null ? null", html)
        self.assertIn("previous.finishedAt === null ? null", html)
        self.assertIn("completeCamera: false", html)
        self.assertIn("ovizSkyLayerTransitionWaiters", html)
        self.assertIn('data.type === "oviz-aladin-sky-layer-transition-complete"', html)
        self.assertIn("transitionOwnsCameraAndTime", html)
        self.assertIn('typeof activeActionRun !== "undefined" && activeActionRun', html)
        self.assertIn("actionHeldTraceOpacityByKey", html)
        self.assertIn("captureActionTraceOpacityMap", html)
        self.assertIn("preserveLegendTransitionFrame: true", html)
        self.assertIn("actionHeldAppearanceRollback", html)
        self.assertIn("sourceAppearanceComposite", html)
        self.assertIn("preserveRenderedSelection: true", html)
        self.assertIn("restorePresentation: false", html)
        self.assertIn('destinationViewMode === "earth"', html)
        self.assertIn("if (transition !== ovizStateTransition)", html)
        self.assertIn("finishActionHeldAppearanceRollback", html)
        cancel_action_body = html.split(
            "function ovizCancelActionWithoutSnap(reason", 1
        )[1].split("function ovizCancelStateTransitionWithoutSnap", 1)[0]
        self.assertIn("actionHeldTraceOpacityByKey = null", cancel_action_body)
        active_action_cancel_body = cancel_action_body.split(
            'if (typeof cancelActionRun === "function" && activeActionRun)', 1
        )[1].split("activeActionRun = null", 1)[0]
        self.assertIn("const cancelled = cancelActionRun", active_action_cancel_body)
        self.assertIn("actionHeldTraceOpacityByKey = null", active_action_cancel_body)
        self.assertLess(
            active_action_cancel_body.index("actionHeldTraceOpacityByKey = null"),
            active_action_cancel_body.index("return cancelled"),
        )
        self.assertNotIn("return cancelActionRun", active_action_cancel_body)
        self.assertNotIn("__STATE_RUNTIME_JS__", html)
        fail_state_body = html.split(
            "function ovizFailStateTransition(transition, err)", 1
        )[1].split("function ovizGoToState", 1)[0]
        stale_guard_index = fail_state_body.index("if (transition !== ovizStateTransition)")
        sky_cancel_index = fail_state_body.index("cancelSkyViewTransitionAnimations")
        self.assertLess(stale_guard_index, sky_cancel_index)
        trace_visible_body = html.split("function traceVisible(trace)", 1)[1].split(
            "function isGalacticReferenceTrace", 1
        )[0]
        self.assertIn("Object.prototype.hasOwnProperty.call(legendState, trace.key)", trace_visible_body)
        self.assertIn(": mode === true", trace_visible_body)

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
        self.assertIn('name: "camera+time"', html)
        self.assertIn('domains: ["camera", "time"]', html)

    def test_states_drawer_mount_removes_serialized_or_stale_copy(self):
        html = ThreeJSFigure({
            "width": 640,
            "height": 480,
            "frames": [],
            "initial_state": {},
        }).to_html(compress_scene_spec=False)
        build_body = html.split(
            "function ovizBuildStatesDrawer()", 1
        )[1].split("function ovizRenderStatesDrawer()", 1)[0]

        remove_index = build_body.index('querySelectorAll(".oviz-states-shell")')
        create_index = build_body.index('document.createElement("div")')
        prepend_index = build_body.index("widgetMenuEl.prepend(ovizStatesShellEl)")
        self.assertLess(remove_index, create_index)
        self.assertLess(create_index, prepend_index)
        self.assertIn(
            'ovizStatesMode === "edit" ? "Present" : "Edit"',
            html,
        )

    def test_states_export_button_prompts_for_filename(self):
        html = ThreeJSFigure({
            "width": 640,
            "height": 480,
            "frames": [],
            "initial_state": {},
        }).to_html(compress_scene_spec=False)
        prompt_body = html.split(
            "function ovizPromptExportStatesHtml()", 1
        )[1].split("async function ovizExportStatesHtml", 1)[0]

        self.assertIn("window.prompt(", prompt_body)
        self.assertIn('"Name the exported HTML file"', prompt_body)
        self.assertIn("ovizDefaultStatesExportFilename()", prompt_body)
        self.assertIn("if (filename === null)", prompt_body)
        self.assertIn("ovizNormalizeStatesExportFilename(filename)", prompt_body)
        self.assertIn(
            'ovizMakeButton("Export HTML", ovizPromptExportStatesHtml)',
            html,
        )
        self.assertIn(
            'ovizMakeButton("Export Present Only", ovizPromptExportStatesPresentOnlyHtml)',
            html,
        )
        self.assertIn('"Name the present-only HTML file"', html)
        self.assertIn("compact.present_only = presentOnly", html)
        self.assertIn("delete exportSceneSpec.deck", html)

    def test_states_export_filename_is_safe_and_gets_html_extension(self):
        html = ThreeJSFigure({
            "width": 640,
            "height": 480,
            "frames": [],
            "initial_state": {},
        }).to_html(compress_scene_spec=False)
        normalize_body = html.split(
            "function ovizNormalizeStatesExportFilename(value, options = {})", 1
        )[1].split("function ovizPromptExportStatesHtml", 1)[0]

        self.assertIn('.replace(/[\\\\/:*?"<>|]+/g, "_")', normalize_body)
        self.assertIn(r"/\.html?$/i.test(cleaned)", normalize_body)
        self.assertIn('cleaned + ".html"', normalize_body)

    @unittest.skipIf(shutil.which("node") is None, "node is not available")
    def test_present_only_export_locks_states_and_excludes_slides(self):
        html = ThreeJSFigure({
            "width": 640,
            "height": 480,
            "frames": [],
            "initial_state": {},
        }).to_html(compress_scene_spec=False)
        export_source = (
            "async function ovizExportStatesHtml(options = {})"
            + html.split("async function ovizExportStatesHtml(options = {})", 1)[1].split(
                "function ovizOpenDraftDb", 1
            )[0]
        )
        script = f"""
        const sceneSpec = {{
          title: "Fixture",
          width: 640,
          height: 480,
          initial_state: {{ old: true }},
          deck: {{ available: true, slides: [{{ id: "slide-1" }}] }},
        }};
        const ovizOriginalSceneInitialState = {{ original: true }};
        const root = {{ clientWidth: 800, clientHeight: 600 }};
        const ovizStatesProject = {{ revision: 3 }};
        function ovizStatesPublicProject() {{
          return {{ revision: 3, items: [{{ id: "state-1" }}], assets: {{}} }};
        }}
        async function ovizCompactProjectForStorage(value) {{
          return JSON.parse(JSON.stringify(value));
        }}
        function safeJsonClone(value, fallback) {{
          try {{ return JSON.parse(JSON.stringify(value)); }} catch (_err) {{ return fallback; }}
        }}
        function ovizStatesClone(value, fallback) {{ return safeJsonClone(value, fallback); }}
        function ovizDeckExportSpec() {{
          return {{ available: true, slides: [{{ id: "exported-slide" }}] }};
        }}
        async function buildExportHtml(value) {{ return JSON.stringify(value); }}
        {export_source}
        (async () => {{
          const presentOnly = JSON.parse(await ovizExportStatesHtml({{ download: false, presentOnly: true }}));
          const editable = JSON.parse(await ovizExportStatesHtml({{ download: false }}));
          process.stdout.write(JSON.stringify({{ presentOnly, editable }}));
        }})().catch((error) => {{ throw error; }});
        """
        result = subprocess.run(
            ["node"],
            input=script,
            text=True,
            capture_output=True,
            check=True,
        )
        payload = json.loads(result.stdout)

        self.assertTrue(payload["presentOnly"]["states"]["present_only"])
        self.assertEqual(payload["presentOnly"]["states"]["default_mode"], "present")
        self.assertNotIn("deck", payload["presentOnly"])
        self.assertEqual(payload["presentOnly"]["initial_state"], {"original": True})
        self.assertFalse(payload["editable"]["states"]["present_only"])
        self.assertEqual(
            payload["editable"]["deck"]["slides"],
            [{"id": "exported-slide"}],
        )

    @unittest.skipIf(shutil.which("node") is None, "node is not available")
    def test_present_only_initialization_applies_first_authored_state(self):
        html = ThreeJSFigure({
            "width": 640,
            "height": 480,
            "frames": [],
            "initial_state": {},
        }).to_html(compress_scene_spec=False)
        helper_source = (
            "async function ovizApplyPresentOnlyInitialState()"
            + html.split("async function ovizApplyPresentOnlyInitialState()", 1)[1].split(
                "async function initializeOvizStates()", 1
            )[0]
        )
        script = f"""
        const ovizStatesProject = {{
          present_only: true,
          items: [{{ id: "first", snapshot: {{ marker: "first-snapshot" }} }}],
        }};
        let ovizActiveStateId = null;
        let ovizStateDirty = true;
        const root = {{ dataset: {{}} }};
        const applied = [];
        function ovizStateTargetFor() {{
          return {{ id: "first", snapshot: ovizStatesProject.items[0].snapshot }};
        }}
        function ovizHydrateAssets(value) {{ return value; }}
        async function ovizApplyStateImmediately(snapshot, options) {{
          applied.push({{ snapshot, options }});
          return snapshot;
        }}
        {helper_source}
        (async () => {{
          const result = await ovizApplyPresentOnlyInitialState();
          process.stdout.write(JSON.stringify({{
            result,
            active: ovizActiveStateId,
            dirty: ovizStateDirty,
            dataset: root.dataset,
            applied,
          }}));
        }})().catch((error) => {{ throw error; }});
        """
        result = subprocess.run(
            ["node"],
            input=script,
            text=True,
            capture_output=True,
            check=True,
        )
        payload = json.loads(result.stdout)

        self.assertTrue(payload["result"])
        self.assertEqual(payload["active"], "first")
        self.assertFalse(payload["dirty"])
        self.assertEqual(payload["dataset"]["activeStatePosition"], "1")
        self.assertEqual(payload["dataset"]["presentOnlyInitialState"], "first")
        self.assertEqual(payload["applied"][0]["snapshot"]["marker"], "first-snapshot")
        self.assertTrue(payload["applied"][0]["options"]["forceSkyBackground"])

    def test_3d_camera_and_time_share_one_phase_by_default(self):
        html = ThreeJSFigure({
            "width": 640,
            "height": 480,
            "frames": [],
            "initial_state": {},
        }).to_html(compress_scene_spec=False)
        phase_builder = html.split(
            "function ovizBuildTransitionPhases(", 1
        )[1].split("function ovizTransitionPhaseState", 1)[0]

        self.assertIn("const concurrent3dCameraTime = Boolean(", phase_builder)
        self.assertIn('fromViewMode !== "earth"', phase_builder)
        self.assertIn('toViewMode !== "earth"', phase_builder)
        self.assertIn('{ name: "camera+time", domains: ["camera", "time"] }', phase_builder)
        self.assertIn("const minimumDurationMs = phaseSpecs.length * phaseMinimumDurationMs", phase_builder)

    def test_shared_3d_camera_time_phase_drives_both_domains_with_same_progress(self):
        html = ThreeJSFigure({
            "width": 640,
            "height": 480,
            "frames": [],
            "initial_state": {},
        }).to_html(compress_scene_spec=False)
        progress_body = html.split(
            "function ovizTransitionPhaseProgress(", 1
        )[1].split("function ovizTraceCandidatesForFrameValue", 1)[0]

        self.assertIn("Array.isArray(item.domains)", progress_body)
        self.assertIn("item.domains.includes(phaseName)", progress_body)
        self.assertIn(
            'const cameraRaw = ovizTransitionPhaseProgress(transition, "camera", elapsedMs)',
            html,
        )
        self.assertIn(
            'const timeRaw = ovizTransitionPhaseProgress(transition, "time", elapsedMs)',
            html,
        )

    def test_sky_to_3d_time_change_exits_sky_before_timeline_motion(self):
        html = ThreeJSFigure({
            "width": 640,
            "height": 480,
            "frames": [],
            "initial_state": {},
        }).to_html(compress_scene_spec=False)
        phase_builder = html.split(
            "function ovizBuildTransitionPhases(", 1
        )[1].split("function ovizRenderStateTimelineFrameLikeSlider", 1)[0]

        self.assertIn('fromViewMode === "earth" && toViewMode !== "earth"', phase_builder)
        self.assertIn('["camera", "time", "appearance"]', phase_builder)
        self.assertLess(
            phase_builder.index('["camera", "time", "appearance"]'),
            phase_builder.index('name: "camera+time"'),
        )
        self.assertIn("function ovizStatePhasePlanAfterViewButtonHandoff", phase_builder)
        self.assertIn('.filter((domain) => domain !== "camera")', phase_builder)

    @unittest.skipIf(shutil.which("node") is None, "node is not available")
    def test_view_button_handoff_removes_duplicate_state_camera_phase(self):
        html = ThreeJSFigure({
            "width": 640,
            "height": 480,
            "frames": [],
            "initial_state": {},
        }).to_html(compress_scene_spec=False)
        helper_source = (
            "function ovizStatePhasePlanAfterViewButtonHandoff("
            + html.split(
                "function ovizStatePhasePlanAfterViewButtonHandoff(", 1
            )[1].split("function ovizRenderStateTimelineFrameLikeSlider", 1)[0]
        )
        script = f"""
        {helper_source}
        const result = ovizStatePhasePlanAfterViewButtonHandoff({{
          changed: {{ camera: true, time: true, appearance: true }},
          requestedDurationMs: 2400,
          effectiveDurationMs: 2400,
          phases: [
            {{ name: "camera", domains: ["camera"], startMs: 0, endMs: 800 }},
            {{ name: "time", domains: ["time"], startMs: 800, endMs: 1600 }},
            {{ name: "appearance", domains: ["appearance"], startMs: 1600, endMs: 2400 }},
          ],
        }});
        process.stdout.write(JSON.stringify(result));
        """
        result = subprocess.run(
            ["node"],
            input=script,
            text=True,
            capture_output=True,
            check=True,
        )
        payload = json.loads(result.stdout)
        self.assertFalse(payload["changed"]["camera"])
        self.assertEqual(
            [(phase["name"], phase["startMs"], phase["endMs"]) for phase in payload["phases"]],
            [("time", 0, 800), ("appearance", 800, 1600)],
        )
        self.assertEqual(payload["effectiveDurationMs"], 1600)

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
        self.assertIn(
            "if (useSelectionPolygon || selectionTransitionActive)",
            html,
        )
        self.assertIn("float sourceSelectionWeight = dot(", html)
        self.assertIn("selectionWeight = mix(", html)
        self.assertIn("ovizHeldSelectionTransition", html)
        self.assertIn("preserveRenderedSelection: true", html)
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
        self.assertIn("ovizRenderedSceneFidelityDifferences", finish_body)

    def test_exact_state_completion_defers_saved_motion_until_after_fidelity(self):
        html = ThreeJSFigure({
            "width": 640,
            "height": 480,
            "frames": [],
            "initial_state": {},
        }).to_html(compress_scene_spec=False)
        apply_body = html.split(
            "function ovizApplyStateImmediately(snapshot, options = {})", 1
        )[1].split("function ovizPointFrom", 1)[0]
        finish_body = html.split(
            "async function ovizFinishStateTransition(transition)", 1
        )[1].split("function ovizFailStateTransition", 1)[0]

        self.assertLess(
            apply_body.index("setCameraAutoOrbitEnabled(false);"),
            apply_body.index("applyViewerStateSyncInternal(hydrated, options);"),
        )
        self.assertIn("options.deferMotionActivation === true", apply_body)
        self.assertEqual(finish_body.count("deferMotionActivation: true"), 2)
        self.assertLess(
            finish_body.index("ovizStateTransition = null;"),
            finish_body.index("syncCameraAutoOrbitMode();"),
        )
        self.assertLess(
            finish_body.index("State fidelity check failed"),
            finish_body.index("syncCameraAutoOrbitMode();"),
        )

    def test_exact_camera_restore_drains_orbit_damping_before_saved_pose(self):
        html = ThreeJSFigure({
            "width": 640,
            "height": 480,
            "frames": [],
            "initial_state": {},
        }).to_html(compress_scene_spec=False)
        helper = html.split(
            "function ovizApplyCapturedCameraState(snapshot)", 1
        )[1].split("function ovizStateFidelityDifferences", 1)[0]

        drain_index = helper.index("controls.enableDamping = false;")
        update_index = helper.index("controls.update();", drain_index)
        restore_index = helper.index("camera.position.set", update_index)
        matrix_index = helper.index("camera.updateMatrixWorld(true);", restore_index)

        self.assertLess(drain_index, update_index)
        self.assertLess(update_index, restore_index)
        self.assertLess(restore_index, matrix_index)
        self.assertNotIn("controls.update();", helper[matrix_index:])

    def test_final_retained_target_is_painted_before_exact_restoration(self):
        html = ThreeJSFigure({
            "width": 640,
            "height": 480,
            "frames": [],
            "initial_state": {},
        }).to_html(compress_scene_spec=False)
        update_body = html.split(
            "function updateOvizStateTransition(now)", 1
        )[1].split("async function ovizFinishStateTransition", 1)[0]

        self.assertIn("targetFrameLatched: false", html)
        self.assertIn("if (!transition.targetFrameLatched)", update_body)
        self.assertLess(
            update_body.index("transition.targetFrameLatched = true"),
            update_body.index("ovizFinishStateTransition(transition)"),
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

    def test_transition_only_work_is_dirty_driven_and_diagnostics_are_throttled(self):
        html = ThreeJSFigure({
            "width": 640,
            "height": 480,
            "frames": [],
            "initial_state": {},
        }).to_html(compress_scene_spec=False)
        update_body = html.split("function updateOvizStateTransition(now)", 1)[1].split(
            "async function ovizFinishStateTransition", 1
        )[0]
        diagnostics_body = html.split("function ovizWriteTransitionDiagnostics", 1)[1].split(
            "function ovizStatesClone", 1
        )[0]

        self.assertIn("lastAppliedAppearanceProgress", update_body)
        self.assertIn(
            "Math.abs(transition.lastAppliedAppearanceProgress - appearanceProgress)",
            update_body,
        )
        self.assertIn("const timeSceneDirty = Boolean(", update_body)
        self.assertIn("transition.phasePlan.changed.time", update_body)
        self.assertIn("transition.lastRenderedTimeProgress - timeRaw", update_body)
        self.assertIn("const appearanceSceneDirty = Boolean(", update_body)
        self.assertIn("transition.phasePlan.changed.appearance", update_body)
        self.assertIn(
            "transition.lastRenderedAppearanceProgress - appearanceRaw",
            update_body,
        )
        self.assertIn("ovizRenderStateTimelineFrameLikeSlider", update_body)
        self.assertIn("if (timeSceneDirty && !appearanceSceneDirty)", update_body)
        self.assertIn("transition.skippedSceneUpdateCount += 1", update_body)
        self.assertNotIn("if (appearanceRaw > 0.0 || timeRaw > 0.0)", update_body)
        self.assertEqual(update_body.count("updateTimelineUi("), 0)
        self.assertEqual(update_body.count("updateTimelineMotionOpacity();"), 0)
        self.assertIn("< 100.0", diagnostics_body)

    def test_state_time_phase_uses_slider_style_stable_frame_rendering(self):
        html = ThreeJSFigure({
            "width": 640,
            "height": 480,
            "frames": [],
            "initial_state": {},
        }).to_html(compress_scene_spec=False)
        timeline_body = html.split(
            "function ovizRenderStateTimelineFrameLikeSlider(", 1
        )[1].split("function ovizTransitionPhaseState", 1)[0]

        self.assertIn("const stableFrameIndex = clampFrameIndex(clampedValue)", timeline_body)
        self.assertIn("transition.lastRenderedTimelineFrameIndex !== stableFrameIndex", timeline_body)
        self.assertIn("renderFrameScene(frame, stableTimeMyr", timeline_body)
        self.assertIn("updateWidgets: false", timeline_body)
        self.assertIn("updateTimelineUi(clampedValue, frameTimeForValue(clampedValue))", timeline_body)
        self.assertNotIn("renderInterpolatedFrameValue", timeline_body)

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
