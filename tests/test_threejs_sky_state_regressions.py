from __future__ import annotations

import re
import unittest

from oviz.threejs_figure import ThreeJSFigure


def _function_region(source: str, name: str, next_name: str) -> str:
    start_marker = f"function {name}("
    start = source.find(start_marker)
    if start < 0:
        raise AssertionError(f"Missing JavaScript function {name}")
    end_marker = f"function {next_name}("
    end = source.find(end_marker, start + len(start_marker))
    if end < 0:
        raise AssertionError(f"Missing JavaScript function {next_name} after {name}")
    return source[start:end]


class ThreeJSSkyStateRegressionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.html = ThreeJSFigure({
            "width": 640,
            "height": 480,
            "frames": [],
            "initial_state": {},
        }).to_html(compress_scene_spec=False)

    def assertHtmlContains(self, text: str) -> None:
        self.assertTrue(text in self.html, f"Missing rendered JavaScript contract: {text}")

    def test_sky_to_sky_uses_spherical_path_not_cartesian_camera_lerp(self):
        begin_transition = _function_region(
            self.html,
            "ovizBeginStateTransition",
            "ovizApplyTransitionNumericControls",
        )
        update_transition = _function_region(
            self.html,
            "updateOvizStateTransition",
            "ovizFinishStateTransition",
        )

        self.assertIn("viewTransitionKind", begin_transition)
        self.assertIn('"earth-to-earth"', begin_transition)
        self.assertHtmlContains('type: "oviz-sky-background-transition"')

        earth_branch = update_transition.find(
            'if (transition.viewTransitionKind === "earth-to-earth")'
        )
        spherical_apply = update_transition.find(
            "ovizApplyEarthCameraTrack(transition.earthCameraTrack, progress);",
            earth_branch,
        )
        native_mode_apply = update_transition.find(
            "ovizApplyNativeViewCameraTrack(transition.nativeCameraTrack, progress);",
            spherical_apply,
        )
        generic_guard = update_transition.find(
            "else if (!transition.nativeViewTransition && !targetEarthViewLocked)",
            native_mode_apply,
        )
        cartesian_apply = update_transition.find("camera.position.set(", generic_guard)
        self.assertTrue(
            0 <= earth_branch < spherical_apply < native_mode_apply < generic_guard < cartesian_apply,
            "Sky and mode-changing tracks must be mutually exclusive from generic Cartesian interpolation.",
        )
        self.assertNotIn("enterEarthViewFromCurrentCamera(", begin_transition)
        self.assertNotIn("exitEarthViewToCameraState(", begin_transition)
        self.assertIn("const synchronizedStart = performance.now();", begin_transition)

    def test_states_forward_the_rendered_camera_instead_of_starting_a_second_sky_animation(self):
        begin_transition = _function_region(
            self.html,
            "ovizBeginStateTransition",
            "ovizApplyTransitionNumericControls",
        )
        update_background = _function_region(
            self.html,
            "updateSkyDomeBackgroundFrame",
            "normalizeSkyAperturePreset",
        )

        self.assertNotIn("startSkyDomeBackgroundProgrammaticTransition({", begin_transition)
        self.assertIn("ovizStateTransition.skyBackgroundPromise = null;", begin_transition)
        self.assertIn("stateCameraTransitionActive", update_background)
        self.assertIn("skyDomeBackgroundUserCameraActive || stateCameraTransitionActive", update_background)
        self.assertIn(") ? 16.0 : 50.0", update_background)
        self.assertIn('type: "oviz-sky-background-view"', update_background)
        self.assertIn(
            "signature === skyDomeBackgroundViewSignature",
            update_background,
        )
        self.assertNotIn("< 500.0", update_background)

    def test_live_aladin_background_keeps_the_same_north_up_roll_as_threejs(self):
        self.assertHtmlContains("lockNorthUp: skyDomeBackgroundOnly")
        self.assertHtmlContains("northPoleOrientation: 0")
        self.assertHtmlContains("inertia: !skyDomeBackgroundOnly")

    def test_live_aladin_background_uses_same_origin_direct_pose_updates(self):
        direct_apply = _function_region(
            self.html,
            "applySkyBackgroundViewDirect",
            "setHoveredClusterKey",
        )
        update_background = _function_region(
            self.html,
            "updateSkyDomeBackgroundFrame",
            "normalizeSkyAperturePreset",
        )
        self.assertIn("applySkyBackgroundViewNow(data)", direct_apply)
        self.assertIn("window.OvizSkyBackgroundBridge", direct_apply)
        self.assertIn("skyDomeFrameEl.contentWindow.OvizSkyBackgroundBridge", update_background)
        self.assertIn("bridge.applyView(viewPayload)", update_background)
        self.assertIn("if (!appliedDirectly)", update_background)
        self.assertIn("postMessage(viewPayload", update_background)

    def test_pan_does_not_force_redundant_aladin_fov_redraws(self):
        apply_view = _function_region(
            self.html,
            "applySkyBackgroundViewNow",
            "postSkyBackgroundViewAppliedAfterPaint",
        )
        self.assertIn("lastAppliedSkyBackgroundFovDeg", apply_view)
        self.assertIn("Math.abs(clampedFovDeg - Number(lastAppliedSkyBackgroundFovDeg)) > 1e-5", apply_view)
        self.assertIn("lastAppliedSkyBackgroundFovDeg = clampedFovDeg", apply_view)

    def test_state_finish_does_not_force_duplicate_aladin_redraw(self):
        finish_transition = _function_region(
            self.html,
            "ovizFinishStateTransition",
            "ovizFailStateTransition",
        )
        self.assertIn(
            "ovizApplyStateImmediately(transition.targetSnapshot, {",
            finish_transition,
        )
        self.assertIn("forceSkyBackground: false", finish_transition)
        self.assertIn("postSkyLayersToAladin: false", finish_transition)

    def test_identical_sky_layers_do_not_crossfade_or_reload(self):
        start_layers = _function_region(
            self.html,
            "ovizStartSkyLayerTransition",
            "ovizCreateEarthCameraTrack",
        )
        self.assertIn(
            "JSON.stringify(fromLayers) === JSON.stringify(toLayers)",
            start_layers,
        )
        self.assertIn("return false;", start_layers)
        self.assertIn(
            "postToAladin: options.postSkyLayersToAladin !== false",
            self.html,
        )

    def test_sky_layers_are_resident_and_crossfaded_in_the_iframe(self):
        self.assertHtmlContains('type: "oviz-sky-layer-transition"')
        self.assertHtmlContains('type: "oviz-sky-layer-transition-cancel"')
        self.assertHtmlContains("function startSkyLayerSemanticTransition(data)")
        self.assertHtmlContains("function applySkyLayerTransitionFrame(")
        self.assertHtmlContains("applySkyImageLayerEffectiveOpacity(")
        self.assertHtmlContains("const stackLayers = residentStack")
        self.assertHtmlContains("skyLayerSemanticTransitionSerial")

    def test_live_aladin_iframe_signature_is_independent_of_survey_and_layers(self):
        restore_layers = _function_region(
            self.html,
            "restoreSkyLayerStateFromSnapshot",
            "postSkyLayerStateToAladin",
        )
        update_capture = _function_region(
            self.html,
            "updateSkyDomeCaptureFrame",
            "updateSkyPanel",
        )
        apply_layers = _function_region(
            self.html,
            "applySkyLayerState",
            "scheduleSkyBackgroundView",
        )

        self.assertIn("postSkyLayerStateToAladin();", restore_layers)
        self.assertNotIn("srcdoc", restore_layers)
        self.assertNotIn("srcdoc", apply_layers)

        live_signature = re.search(
            r"const\s+signature\s*=\s*liveAladinBackground\s*"
            r"\?\s*(?P<identity>\"[^\"]+\")\s*:\s*JSON\.stringify\(\{",
            update_capture,
            re.DOTALL,
        )
        self.assertIsNotNone(live_signature)
        self.assertNotRegex(live_signature.group("identity"), r"\b(?:survey|layer)\b")
        self.assertEqual(update_capture.count("skyDomeFrameEl.srcdoc ="), 1)

    def test_aladin_frames_boot_directly_into_the_saved_base_layer(self):
        initial_survey = _function_region(
            self.html,
            "initialAladinFrameSurvey",
            "buildEmptySkySrcdoc",
        )
        update_capture = _function_region(
            self.html,
            "updateSkyDomeCaptureFrame",
            "updateSkyPanel",
        )

        self.assertIn("layers[layers.length - 1]", initial_survey)
        self.assertIn("{ survey: initialAladinFrameSurvey() }", update_capture)

    def test_empty_sky_layer_list_is_an_explicit_restorable_state(self):
        restore_layers = _function_region(
            self.html,
            "restoreSkyLayerStateFromSnapshot",
            "postSkyLayerStateToAladin",
        )

        self.assertIn("skyLayerState = restoredLayers;", restore_layers)
        self.assertNotIn("if (!restoredLayers.length)", restore_layers)
        self.assertRegex(
            restore_layers,
            r"(?s)activeSkyLayerKey\s*=.*?\(skyLayerState\[0\]\s*\?\s*"
            r"skyLayerState\[0\]\.key\s*:\s*\"\"\)",
        )
        self.assertIn("return true;", restore_layers)

    def test_semantic_transition_serial_cancels_and_ignores_stale_callbacks(self):
        self.assertHtmlContains("skyBackgroundSemanticTransitionSerial")
        self.assertHtmlContains("skyBackgroundSemanticTransitionFrame")
        self.assertHtmlContains("activeSkyBackgroundSemanticTransition")
        self.assertHtmlContains('type: "oviz-sky-background-transition"')
        self.assertHtmlContains('type: "oviz-sky-background-transition-cancel"')
        self.assertHtmlContains('data.type === "oviz-sky-background-transition-cancel"')
        self.assertHtmlContains(
            'type: "oviz-aladin-sky-background-transition-complete"',
        )
        self.assertTrue(re.search(
            r"(?:const|let)\s+\w*[Ss]erial\w*\s*=\s*"
            r"\+\+skyBackgroundSemanticTransitionSerial",
            self.html,
        ), "Semantic transition start must advance its cancellation serial.")
        self.assertTrue(re.search(
            r"\w*[Ss]erial\w*\s*!==\s*skyBackgroundSemanticTransitionSerial",
            self.html,
        ), "Stale semantic-transition callbacks must compare their captured serial.")
        self.assertTrue(
            re.search(r"transitionId\s*:\s*[^,}\n]+", self.html),
            "Semantic transition results must preserve the correlated transitionId.",
        )
        self.assertTrue(
            re.search(r"cancelled\s*:\s*true", self.html),
            "Superseded semantic transitions must emit an explicit cancellation result.",
        )

        accept_completion = _function_region(
            self.html,
            "markSkyDomeBackgroundTransitionComplete",
            "finishSkyDomeBackgroundProgrammaticTransition",
        )
        self.assertRegex(
            accept_completion,
            r"String\(transitionId\s*\|\|\s*\"\"\)\s*!==\s*active\.id",
        )
        self.assertRegex(
            accept_completion,
            r"Math\.round\(Number\(seq\)\s*\|\|\s*0\)\s*!==\s*active\.seq",
        )
        self.assertIn("return false;", accept_completion)

    def test_aladin_fallback_script_cannot_terminate_its_own_srcdoc_script(self):
        self.assertHtmlContains(
            "document.write('<script src=\"https://aladin.cds.unistra.fr/AladinLite/"
        )
        self.assertHtmlContains("charset=\"utf-8\"><' + '/script>');")
        self.assertNotIn(
            "charset=\"utf-8\"></script>');",
            self.html,
        )

    def test_parent_child_layer_transitions_share_one_epoch(self):
        begin_transition = _function_region(
            self.html,
            "ovizBeginStateTransition",
            "ovizApplyTransitionNumericControls",
        )
        child_camera = _function_region(
            self.html,
            "startSkyBackgroundSemanticTransition",
            "skyLayerNameFor",
        )
        child_layers = _function_region(
            self.html,
            "startSkyLayerSemanticTransition",
            "scheduleSkyBackgroundView",
        )

        self.assertIn("startedAtEpochMs: ovizTransitionEpochMs(now)", begin_transition)
        self.assertIn("startedAtEpochMs: transition.startedAtEpochMs", self.html)
        self.assertIn("skyBackgroundTransitionStartedAt(data, transitionNow)", child_camera)
        self.assertIn("skyBackgroundTransitionStartedAt(data, transitionNow)", child_layers)

    def test_view_offset_is_interpolated_before_exact_state_restoration(self):
        create_track = _function_region(
            self.html,
            "ovizCreateNativeViewCameraTrack",
            "ovizApplyNativeViewCameraTrack",
        )
        apply_track = _function_region(
            self.html,
            "ovizApplyNativeViewCameraTrack",
            "ovizPrepareEarthViewModeForStateTransition",
        )

        self.assertIn("startViewOffset", create_track)
        self.assertIn("endViewOffset", create_track)
        self.assertIn("ovizApplyInterpolatedViewOffset", apply_track)

    def test_parent_forces_exact_child_target_if_semantic_completion_times_out(self):
        finish_background = _function_region(
            self.html,
            "finishSkyDomeBackgroundProgrammaticTransition",
            "markSkyDomeBackgroundViewApplied",
        )

        self.assertIn("!active.completed && active.sent", finish_background)
        self.assertIn('type: "oviz-sky-background-transition-cancel"', finish_background)
        self.assertIn('type: "oviz-sky-layer-transition-cancel"', finish_background)
        self.assertIn("updateSkyDomeBackgroundFrame(", finish_background)

    def test_live_sky_background_never_uses_a_planar_css_prediction(self):
        predictive_update = _function_region(
            self.html,
            "updateSkyDomeBackgroundPredictiveTransform",
            "recordSkyDomeBackgroundSentView",
        )
        ordinary_update = _function_region(
            self.html,
            "scheduleSkyBackgroundView",
            "setHoveredClusterKey",
        )

        self.assertNotIn("translate3d", predictive_update)
        self.assertNotIn("scale(", predictive_update)
        self.assertIn("clearSkyDomeBackgroundPredictiveTransform();", predictive_update)
        self.assertIn("pendingSkyBackgroundView = data;", ordinary_update)
        self.assertIn("applySkyBackgroundViewNow(nextView);", ordinary_update)
        self.assertIn("requestAnimationFrame", ordinary_update)


if __name__ == "__main__":
    unittest.main()
