import base64
import gzip
import json
import os
import re
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from astropy.io import fits

os.environ.setdefault("MPLCONFIGDIR", "/tmp")
os.environ.setdefault("XDG_CACHE_HOME", "/tmp")
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from oviz import Layer, LayerCollection, Scene3D, galactic_lite_profile, website_background_profile  # noqa: E402
from oviz.threejs_figure import ThreeJSFigure  # noqa: E402
from oviz.viz import Animate3D  # noqa: E402


class _FakeCluster:
    def __init__(self, show_tracks=False):
        self.data_name = "Cluster A"
        self.show_tracks = show_tracks
        self.color = "cyan"
        self.opacity = 0.8
        self.marker_style = "circle"
        self.integrated = True
        self.size_by_n_stars = False
        self.min_size = 2
        self.max_size = 6
        self.colormap = None
        self.cmin = None
        self.cmax = None
        self.df = pd.DataFrame(
            {
                "x": [0.0],
                "y": [0.0],
                "z": [0.0],
                "U": [0.0],
                "V": [0.0],
                "W": [0.0],
                "n_stars": [42.0],
            }
        )
        self.df_int = pd.DataFrame(
            {
                "x": [0.0, 10.0],
                "y": [0.0, -5.0],
                "z": [0.0, 2.0],
                "x_rot": [0.0, 10.0],
                "y_rot": [0.0, -5.0],
                "z_rot": [0.0, 2.0],
                "x_helio": [0.0, 10.0],
                "y_helio": [0.0, -5.0],
                "z_helio": [0.0, 2.0],
                "age_myr": [12.0, 12.0],
                "name": ["member_1", "member_1"],
                "time": [0.0, -1.0],
                "size": [6.0, 4.0],
                "n_stars": [42.0, 42.0],
            }
        )


class _FakeCollection:
    def __init__(self, show_tracks=False):
        self.cluster = _FakeCluster(show_tracks=show_tracks)
        self.time = None

    def get_all_clusters(self):
        return [self.cluster]

    def get_cluster(self, _name):
        return self.cluster

    def integrate_all_orbits(self, time, **_kwargs):
        self.time = np.array(time, dtype=float)

    def set_all_cluster_sizes(self, *_args, **_kwargs):
        return None


class _FakeFamilyCollection(_FakeCollection):
    def __init__(self, show_tracks=False):
        super().__init__(show_tracks=show_tracks)
        self.family = _FakeCluster(show_tracks=show_tracks)
        self.family.data_name = "Family A"
        self.family.color = "violet"

    def get_all_clusters(self):
        return [self.cluster, self.family]


class ThreeJSRendererTests(unittest.TestCase):
    def _embedded_scene_spec_marker(self, html):
        start_marker = "/*__SCENE_SPEC_START__*/"
        end_marker = "/*__SCENE_SPEC_END__*/"
        start = html.find(start_marker)
        self.assertGreaterEqual(start, 0)
        end = html.find(end_marker, start)
        self.assertGreater(end, start)
        return html[start : end + len(end_marker)]

    def _decode_embedded_compressed_scene_spec(self, html):
        pattern = (
            r"/\*__SCENE_SPEC_START__\*/const sceneSpec = await "
            r"inflateOvizGzipBase64SceneSpec\((.*?)\);/\*__SCENE_SPEC_END__\*/"
        )
        match = re.search(pattern, html, re.S)
        self.assertIsNotNone(match)
        encoded_expr = match.group(1).strip()
        payload_match = re.fullmatch(r"readOvizSceneSpecPayload\((.*?)\)", encoded_expr, re.S)
        if payload_match:
            payload_id = json.loads(payload_match.group(1))
            payload_pattern = (
                rf'<script type="application/octet-stream" id="{re.escape(payload_id)}">\n'
                r"(.*?)\n</script>"
            )
            payload_script_match = re.search(payload_pattern, html, re.S)
            self.assertIsNotNone(payload_script_match)
            encoded = re.sub(r"\s+", "", payload_script_match.group(1))
        elif encoded_expr.endswith('.join("")'):
            encoded = "".join(json.loads(encoded_expr[: -len('.join("")')]))
        else:
            encoded = json.loads(encoded_expr)
        return json.loads(gzip.decompress(base64.b64decode(encoded)).decode("utf-8"))

    def test_make_plot_keeps_plotly_renderer_default(self):
        viz = Animate3D(_FakeCollection(), figure_theme="dark")
        fig = viz.make_plot(time=np.array([0.0, -1.0]), renderer="plotly", show=False)

        self.assertIsInstance(fig, go.Figure)
        step_labels = [step["label"] for step in viz.fig_dict["layout"]["sliders"][0]["steps"]]
        self.assertEqual(step_labels, ["-1.0", "0.0"])

    def test_make_plot_can_return_threejs_html_wrapper(self):
        viz = Animate3D(_FakeCollection(show_tracks=True), figure_theme="dark")
        fig = viz.make_plot(time=np.array([0.0, -1.0]), renderer="threejs", show=False)

        html = fig.to_html()
        repr_html = fig._repr_html_()

        self.assertEqual(fig.__class__.__name__, "ThreeJSFigure")
        self.assertIn("three.min.js", html)
        self.assertIn("oviz-three-legend-popover", html)
        self.assertIn(">Legend<", html)
        self.assertIn(">Group<", html)
        self.assertNotIn("border: 1px solid rgba(170, 176, 185, 0.42)", html)
        self.assertIn("oviz-three-group-dropdown", html)
        self.assertIn("oviz-three-group-option", html)
        self.assertIn("Cluster A", html)
        self.assertIn("Shift+drag or use Lasso", html)
        self.assertIn("Enable click select", html)
        self.assertIn("Lasso volumetric data", html)
        self.assertIn("let lassoVolumeSelectionEnabled = true;", html)
        self.assertTrue(fig.scene_spec["initial_state"]["lasso_volume_selection_enabled"])
        self.assertIn("oviz-three-controls-toggle", html)
        self.assertIn("oviz-three-sky-controls-toggle", html)
        self.assertIn("oviz-three-zen-mode", html)
        self.assertIn("oviz-three-reset-view", html)
        self.assertIn("data-zen=\"false\"", html)
        self.assertIn('toggleButton.addEventListener("dblclick"', html)
        self.assertIn("activeLegendEditorKey", html)
        self.assertIn("focusSelectionKey", html)
        self.assertIn("Camera FOV", html)
        self.assertIn('class="oviz-three-camera-fov" type="range" min="0.05"', html)
        self.assertIn('class="oviz-three-camera-fov" type="range" min="0.05" max="120"', html)
        self.assertIn("new THREE.PerspectiveCamera(60", html)
        self.assertIn("Focus group", html)
        self.assertIn("Star glow", html)
        self.assertIn("let globalPointGlowStrength = 0.60;", html)
        self.assertIn("const pointScaleBaseline = Math.max(Number(sceneSpec.point_size_baseline_scale) || (4.0 / 3.0), 0.01);", html)
        self.assertIn("const pointScale = (Math.max(sceneSpec.max_span || 1, 1) / 2600.0) * pointScaleBaseline;", html)
        self.assertIn("function starCoreMaterialFor", html)
        self.assertIn("function starCoreScaleForPoint", html)
        self.assertNotIn("function starBloomMaterialFor", html)
        self.assertNotIn("function starBloomScaleForPoint", html)
        self.assertIn("Glow keeps its old apparent radius while the regular marker stays on the smaller baseline.", html)
        self.assertIn("6.45 + 1.26 * strength", html)
        self.assertIn("function writeStarPsfTexture", html)
        self.assertIn("const airyCore = 0.82", html)
        self.assertIn("const saturatedCore = 1.00", html)
        self.assertIn("blending: THREE.AdditiveBlending", html)
        self.assertNotIn("function diffractionSpikeAlpha", html)
        self.assertIn("THREE.LinearMipmapLinearFilter", html)
        self.assertNotIn("screenPixelsToWorldScale", html)
        self.assertIn('registerCameraResponsivePointSprite(glowSprite, "glow"', html)
        self.assertIn('registerCameraResponsivePointSprite(coreSprite, "core"', html)
        self.assertNotIn('registerCameraResponsivePointSprite(bloomSprite, "bloom"', html)
        self.assertIn('registerCameraResponsivePointSprite(sprite, "marker"', html)
        self.assertIn("Size points by n_stars", html)
        self.assertIn('"n_stars": 42.0', html)
        self.assertIn('"has_n_stars"', html)
        self.assertIn("Fade time (Myr)", html)
        self.assertIn("Fade in and out", html)
        self.assertIn("Cluster Filter", html)
        self.assertIn("Dendrogram", html)
        self.assertIn("Parameter", html)
        self.assertIn("Birth to older track", html)
        self.assertIn("Birth to birth", html)
        self.assertIn("Keyboard help", html)
        self.assertIn(">Graphite<", html)
        self.assertIn(">Aurora<", html)
        self.assertIn("Keyboard controls are active as soon as the viewer loads.", html)
        self.assertIn("Shift + W / A / S / D", html)
        self.assertIn("Toggle between 3D View and Sky View.", html)
        self.assertIn(">3D View<", html)
        self.assertIn("oviz-three-earth-view-toggle", html)
        self.assertIn("Reset camera", html)
        self.assertIn("oviz-three-scale-bar", html)
        self.assertIn("oviz-three-save-state", html)
        self.assertIn('canvas.addEventListener("dblclick", onCanvasDoubleClick);', html)
        self.assertIn('window.addEventListener("keydown", onKeyDown);', html)
        self.assertIn('window.addEventListener("keyup", onKeyUp);', html)
        self.assertIn("updateKeyboardMotion(deltaSeconds);", html)
        self.assertIn("data:text/html", repr_html)
        self.assertIn("alphaTest: 0.15", html)
        self.assertIn("width: 100vw", html)
        self.assertIn("height: 100vh", html)
        self.assertEqual(viz.fig_dict["renderer"], "threejs")
        self.assertEqual(viz.fig_dict["camera_up"], {"x": 0.0, "y": 0.0, "z": 1.0})
        self.assertEqual(viz.fig_dict["animation"]["fade_in_time_default"], 5.0)
        self.assertFalse(viz.fig_dict["animation"]["fade_in_and_out_default"])
        self.assertTrue(viz.fig_dict["cluster_filter"]["enabled"])
        self.assertEqual(viz.fig_dict["cluster_filter"]["default_parameter_key"], "age_now_myr")
        self.assertTrue(viz.fig_dict["dendrogram"]["enabled"])
        cluster_legend = next(item for item in viz.fig_dict["legend"]["items"] if item["name"] == "Cluster A")
        self.assertEqual(cluster_legend["color"], "#00ffff")

        with tempfile.TemporaryDirectory() as tmp_dir:
            out_file = Path(tmp_dir) / "threejs_scene.html"
            fig.write_html(out_file)
            self.assertTrue(out_file.exists())
            html_text = out_file.read_text()
            self.assertIn("oviz-three-legend-panel", html_text)
            self.assertNotIn("Three.js modules are loaded from a CDN", html_text)

    def test_threejs_html_can_force_compressed_scene_spec_payload(self):
        scene_spec = {
            "renderer": "threejs",
            "width": 900,
            "height": 700,
            "initial_state": {"compact_payload_enabled": True},
            "frames": [
                {
                    "name": "0",
                    "time": 0.0,
                    "traces": [
                        {
                            "name": "Dense Trace",
                            "points": [
                                {
                                    "x": 1.0,
                                    "y": 2.0,
                                    "z": 3.0,
                                    "payload": "abcdefghij" * 200,
                                }
                            ],
                        }
                    ],
                }
            ],
        }

        html = ThreeJSFigure(scene_spec).to_html(compress_scene_spec=True)

        self.assertIn("/*__SCENE_SPEC_START__*/const sceneSpec = await inflateOvizGzipBase64SceneSpec(", html)
        self.assertIn('type="application/octet-stream"', html)
        self.assertIn("readOvizSceneSpecPayload(", html)
        self.assertIn('new DecompressionStream("gzip")', html)
        self.assertIn("support DecompressionStream('gzip')", html)
        self.assertNotIn('support DecompressionStream("gzip")', html)
        self.assertIn("compress_scene_spec=False", html)
        self.assertIn('"compression_method":"gzip+base64"', html)
        self.assertIn('"raw_scene_spec_size_bytes":', html)
        self.assertIn('"compressed_size_bytes":', html)
        self.assertIn('"embedded_base64_size_bytes":', html)
        self.assertNotIn('"payload":"abcdefghij', html)
        self.assertEqual(self._decode_embedded_compressed_scene_spec(html), scene_spec)

    def test_threejs_html_auto_compresses_scene_spec_above_threshold(self):
        scene_spec = {
            "renderer": "threejs",
            "width": 900,
            "height": 700,
            "initial_state": {},
            "frames": [{"name": "0", "time": 0.0, "traces": []}],
        }
        fig = ThreeJSFigure(scene_spec)

        raw_html = fig.to_html(scene_spec_compression_threshold_bytes=10**9)
        auto_compressed_html = fig.to_html(scene_spec_compression_threshold_bytes=1)
        forced_raw_html = fig.to_html(
            compress_scene_spec=False,
            scene_spec_compression_threshold_bytes=1,
        )

        self.assertIn("/*__SCENE_SPEC_START__*/const sceneSpec = {", self._embedded_scene_spec_marker(raw_html))
        self.assertNotIn("await inflateOvizGzipBase64SceneSpec(", self._embedded_scene_spec_marker(raw_html))
        self.assertIn(
            "/*__SCENE_SPEC_START__*/const sceneSpec = await inflateOvizGzipBase64SceneSpec(",
            self._embedded_scene_spec_marker(auto_compressed_html),
        )
        self.assertEqual(self._decode_embedded_compressed_scene_spec(auto_compressed_html), scene_spec)
        self.assertIn("/*__SCENE_SPEC_START__*/const sceneSpec = {", self._embedded_scene_spec_marker(forced_raw_html))
        self.assertNotIn("await inflateOvizGzipBase64SceneSpec(", self._embedded_scene_spec_marker(forced_raw_html))
        self.assertIn('"compression_method":"none"', forced_raw_html)

    def test_threejs_save_state_recompresses_large_scene_spec_exports(self):
        scene_spec = {
            "renderer": "threejs",
            "width": 900,
            "height": 700,
            "initial_state": {"compact_payload_enabled": True},
            "frames": [
                {
                    "name": "0",
                    "time": 0.0,
                    "traces": [
                        {
                            "type": "points",
                            "name": "large",
                            "points": [
                                {
                                    "x": 1.0,
                                    "y": 2.0,
                                    "z": 3.0,
                                    "payload": "abcdefghij" * 200,
                                }
                            ],
                        }
                    ],
                }
            ],
        }

        html = ThreeJSFigure(scene_spec).to_html(compress_scene_spec=True)

        self.assertIn("async function buildExportHtml(exportSceneSpec)", html)
        self.assertIn("const currentHtml = removeExistingSceneSpecPayloadHtml(", html)
        self.assertIn("await gzipBase64EncodeText(sceneJsonText)", html)
        self.assertIn('new CompressionStream("gzip")', html)
        self.assertIn("sceneSpecMarkerWrappedCompressedPayload(payloadId)", html)
        self.assertIn("insertSceneSpecPayloadHtml(exportHtml, payloadHtml)", html)
        self.assertIn("const htmlText = await buildExportHtml(exportSceneSpec);", html)
        self.assertIn("compression_mode: sceneSpecPayloadCompressionMode()", html)
        self.assertIn("compression_method: \"gzip+base64\"", html)
        self.assertIn("<\\/script>`", html)

    def test_threejs_renderer_rejects_actions_for_plotly(self):
        viz = Animate3D(_FakeCollection(show_tracks=True), figure_theme="dark")

        with self.assertRaisesRegex(ValueError, "actions are currently only supported"):
            viz.make_plot(
                time=np.array([0.0, -1.0]),
                renderer="plotly",
                show=False,
                actions=[
                    {
                        "key": "paper-1",
                        "label": "Paper 1",
                        "steps": [{"type": "legend_group", "group": "All"}],
                    }
                ],
            )

    def test_threejs_renderer_serializes_actions(self):
        viz = Animate3D(_FakeCollection(show_tracks=True), figure_theme="dark")
        fig = viz.make_plot(
            time=np.array([0.0, -1.0]),
            renderer="threejs",
            show=False,
            threejs_initial_state={"lite_mode_enabled": True},
            actions=[
                {
                    "key": "paper-1",
                    "label": "Paper 1",
                    "description": "Focus on the main cluster trace",
                    "steps": [
                        {
                            "type": "legend_group",
                            "group": "All",
                            "duration_ms": 320,
                        },
                        {
                            "type": "camera",
                            "target": {"kind": "trace", "name": "Cluster A"},
                            "anchor_target": {"kind": "point", "x": 0.0, "y": 0.0, "z": 0.0},
                            "duration_ms": 900,
                            "distance_scale": 0.6,
                            "camera": {"view_offset": {"x": 0.22, "y": 0.0}},
                            "orbit": {"enabled": True, "speed_multiplier": 1.4, "direction": -1, "persist": True},
                        },
                        {
                            "type": "time",
                            "start": "with_previous",
                            "direction": "backward",
                            "interval_ms": 180,
                            "stop_after_frames": 1,
                        },
                    ],
                }
            ],
        )

        html = fig.to_html()
        actions_payload = viz.fig_dict["actions"]

        self.assertTrue(actions_payload["enabled"])
        self.assertEqual(len(actions_payload["items"]), 1)
        self.assertEqual(actions_payload["items"][0]["label"], "Paper 1")
        self.assertEqual(actions_payload["items"][0]["steps"][0]["group"], "All")
        self.assertEqual(actions_payload["items"][0]["steps"][1]["target"]["key"], "trace-0")
        self.assertEqual(actions_payload["items"][0]["steps"][1]["anchor_target"]["kind"], "point")
        self.assertEqual(actions_payload["items"][0]["steps"][1]["distance_scale"], 0.6)
        self.assertEqual(actions_payload["items"][0]["steps"][1]["orbit"]["direction"], -1.0)
        self.assertEqual(actions_payload["items"][0]["steps"][1]["camera"]["view_offset"]["x"], 0.22)
        self.assertEqual(actions_payload["items"][0]["steps"][2]["direction"], "backward")
        self.assertIn("oviz-three-action-bar", html)
        self.assertIn("oviz-three-action-button", html)
        self.assertIn("startViewerAction", html)
        self.assertIn('controls.addEventListener("start"', html)
        self.assertIn("handleManualCameraInteractionStart", html)
        self.assertIn("clearSelectedAction: options.clearSelectedAction === true", html)
        self.assertIn("function applyActionCameraViewOffset(viewOffset)", html)
        self.assertIn("function restoreTimeIntervalMsForAction(action)", html)
        self.assertIn("intervalMs: restoreTimeIntervalMs", html)
        self.assertIn("Paper 1", html)

    def test_threejs_renderer_keeps_actions_in_minimal_mode(self):
        viz = Animate3D(_FakeCollection(show_tracks=True), figure_theme="dark")
        fig = viz.make_plot(
            time=np.array([0.0, -1.0]),
            renderer="threejs",
            show=False,
            threejs_initial_state={"lite_mode_enabled": True},
            actions=[
                {
                    "key": "paper-1",
                    "label": "Paper 1",
                    "steps": [{"type": "legend_group", "group": "All"}],
                }
            ],
        )

        html = fig.to_html()
        self.assertEqual(viz.fig_dict["export_profile"], "lite")
        self.assertTrue(viz.fig_dict["actions"]["enabled"])
        self.assertIn("oviz-three-action-bar", html)
        self.assertIn("Paper 1", html)
        self.assertIn("initialState.lite_mode_enabled", html)

    def test_threejs_renderer_keeps_grouped_family_traces_in_galactic_lite_mode(self):
        viz = Animate3D(
            _FakeFamilyCollection(show_tracks=True),
            figure_theme="dark",
            trace_grouping_dict={
                "Clusters": ["Cluster A"],
                "Families": ["Family A"],
            },
        )
        viz.make_plot(
            time=np.array([0.0, -1.0]),
            renderer="threejs",
            show=False,
            galactic_mode=True,
            threejs_initial_state={
                "lite_mode_enabled": True,
                "galactic_lite_mode_enabled": True,
            },
        )

        trace_names = [trace["name"] for trace in viz.fig_dict["frames"][0]["traces"]]
        legend_names = [item["name"] for item in viz.fig_dict["legend"]["items"]]

        self.assertIn("Family A", trace_names)
        self.assertIn("Family A", legend_names)

    def test_threejs_renderer_renders_clickable_group_dropdown(self):
        viz = Animate3D(
            _FakeFamilyCollection(show_tracks=True),
            figure_theme="dark",
            trace_grouping_dict={
                "Clusters": ["Cluster A"],
                "Families": ["Family A"],
            },
        )
        fig = viz.make_plot(
            time=np.array([0.0, -1.0]),
            renderer="threejs",
            show=False,
        )

        html = fig.to_html()

        self.assertEqual(viz.fig_dict["group_order"], ["Clusters", "Families", "All"])
        self.assertIn('class="oviz-three-group-dropdown"', html)
        self.assertIn('class="oviz-three-group-trigger"', html)
        self.assertIn('class="oviz-three-group-menu"', html)
        self.assertIn('class="oviz-three-group-menu-list"', html)
        self.assertIn('className = "oviz-three-group-option"', html)
        self.assertIn('option.addEventListener("click"', html)
        self.assertIn('option.addEventListener("keydown"', html)
        self.assertIn("--group-option-delay", html)
        self.assertIn("setLegendGroupDropdownOpen", html)
        self.assertIn('groupDropdownTriggerEl.addEventListener("click"', html)
        self.assertIn("max-height 340ms", html)
        self.assertIn("@keyframes oviz-three-group-option-drop", html)
        self.assertIn("animation: oviz-three-group-option-drop 340ms", html)
        self.assertIn("justify-content: flex-start", html)
        self.assertIn("order: -1", html)
        self.assertIn("font: 760 15px/1.16", html)
        self.assertNotIn("text-decoration: underline", html)
        self.assertNotIn("border-bottom: 1px solid rgba(255, 255, 255, 0.58)", html)
        self.assertIn("transform: translateY(-10px)", html)
        self.assertIn('controls.dataset.visible = "false"', html)
        self.assertIn('controls.dataset.visible = "true"', html)
        self.assertIn('collapseLegendEntryControls', html)
        self.assertIn('max-height 300ms cubic-bezier', html)
        self.assertIn("normalizedLegendGroupOrder", html)
        self.assertIn('groupName.toLowerCase() === "all"', html)
        self.assertIn("setLegendGroup(groupName", html)
        self.assertIn('target.closest(".oviz-three-group-dropdown")', html)
        self.assertNotIn("oviz-three-group-carousel", html)

    def test_threejs_renderer_rejects_invalid_actions(self):
        viz = Animate3D(_FakeCollection(show_tracks=True), figure_theme="dark")

        with self.assertRaisesRegex(ValueError, "lite threejs exports"):
            viz.make_plot(
                time=np.array([0.0, -1.0]),
                renderer="threejs",
                show=False,
                actions=[
                    {
                        "key": "paper-1",
                        "label": "Paper 1",
                        "steps": [{"type": "legend_group", "group": "All"}],
                    }
                ],
            )

    def test_scene3d_supports_profile_helpers(self):
        viz = Scene3D(_FakeCollection(show_tracks=True), figure_theme="dark")
        fig = viz.make_plot(
            time=np.array([0.0, -1.0]),
            renderer="threejs",
            show=False,
            preset="galactic_lite",
            threejs_initial_state={
                "galaxy_image": True,
                "galaxy_image_path": "/tmp/galaxy.png",
            },
        )

        self.assertEqual(fig.__class__.__name__, "ThreeJSFigure")
        self.assertEqual(viz.fig_dict["export_profile"], "lite")
        self.assertTrue(viz.fig_dict["galactic_simple"]["enabled"])
        self.assertEqual(viz.fig_dict["initial_state"]["lite_mode_enabled"], True)
        self.assertEqual(viz.threejs_initial_state["galaxy_image_path"], "/tmp/galaxy.png")

    def test_profile_helpers_normalize_legacy_aliases(self):
        profile = galactic_lite_profile(
            galaxy_image=True,
            galaxy_image_path="/tmp/galaxy.png",
            minimal_mode_enabled=False,
            galactic_plane_opacity=0.4,
        )

        self.assertTrue(profile["lite_mode_enabled"])
        self.assertTrue(profile["galactic_lite_mode_enabled"])
        self.assertEqual(profile["galaxy_image_path"], "/tmp/galaxy.png")
        self.assertEqual(profile["galaxy_image_opacity"], 0.4)

    def test_website_background_profile_has_expected_defaults(self):
        profile = website_background_profile(galaxy_image=True, galaxy_image_path="/tmp/galaxy.png")

        self.assertTrue(profile["lite_mode_enabled"])
        self.assertTrue(profile["galactic_lite_mode_enabled"])
        self.assertTrue(profile["track_orbit_target_to_sun"])
        self.assertFalse(profile["click_selection_enabled"])
        self.assertEqual(profile["camera"]["view_offset"]["x"], 0.22)

    def test_standard_scene_supports_standalone_galaxy_image_overlay(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            image_path = Path(tmp_dir) / "galaxy.jpg"
            image_path.write_bytes(b"\xff\xd8\xff\xd9")
            viz = Animate3D(_FakeCollection(show_tracks=True), figure_theme="dark")
            viz.make_plot(
                time=np.array([0.0, -1.0]),
                renderer="threejs",
                show=False,
                threejs_initial_state={
                    "galaxy_image": True,
                    "galaxy_image_path": str(image_path),
                    "galaxy_image_size_pc": 40000.0,
                    "galaxy_image_opacity": 0.35,
                },
            )

        self.assertFalse(viz.fig_dict["galactic_simple"]["enabled"])
        self.assertEqual(len(viz.fig_dict["image_planes"]), 1)
        image_plane = viz.fig_dict["image_planes"][0]
        self.assertEqual(image_plane["key"], "galaxy-image-overlay")
        self.assertTrue(image_plane["image_data_url"].startswith("data:image/jpeg;base64,"))
        self.assertEqual(image_plane["width_pc"], 40000.0)

        zero_frame = next(frame for frame in viz.fig_dict["frames"] if frame["time"] == 0.0)
        past_frame = next(frame for frame in viz.fig_dict["frames"] if frame["time"] == -1.0)
        zero_image = next(
            item for item in zero_frame["decorations"] if item.get("key") == "galaxy-image-overlay"
        )
        past_image = next(
            item for item in past_frame["decorations"] if item.get("key") == "galaxy-image-overlay"
        )
        self.assertEqual(zero_image["kind"], "image_plane")
        self.assertAlmostEqual(zero_image["opacity"], 0.35)
        self.assertEqual(past_image["opacity"], 0.0)
        self.assertIn("Cluster A", [trace["name"] for trace in zero_frame["traces"]])

    def test_static_lite_scene_disables_timeline_footer(self):
        static_sun_df = pd.DataFrame(
            {
                "x": [0.0],
                "y": [0.0],
                "z": [27.0],
                "name": ["Sun"],
                "age_myr": [8000.0],
            }
        )
        layers = LayerCollection(
            [
                Layer.from_dataframe(
                    static_sun_df,
                    layer_name="Sun",
                    min_size=5.0,
                    max_size=5.0,
                    color="yellow",
                    opacity=1.0,
                    marker_style="circle",
                    assume_stationary=True,
                )
            ]
        )
        viz = Scene3D(
            layers,
            figure_theme="dark",
            xyz_widths=(1000, 1000, 400),
            trace_grouping_dict={"Sun": ["Sun"]},
        )
        fig = viz.make_plot(
            time=np.array([0.0, -1.0]),
            renderer="threejs",
            show=False,
            galactic_mode=True,
            preset="galactic_lite",
        )

        self.assertFalse(viz.fig_dict["timeline"]["enabled"])
        html = fig.to_html()
        self.assertIn('const timelineSpec = sceneSpec.timeline || { enabled: frameSpecs.length > 1 };', html)
        self.assertIn('footerEl.style.display = "none";', html)

    def test_threejs_milky_way_model_is_t0_only(self):
        viz = Animate3D(_FakeCollection(show_tracks=False), figure_theme="dark")
        viz.make_plot(
            time=np.array([0.0, -1.0]),
            renderer="threejs",
            show=False,
            show_milky_way_model=True,
            galactic_mode=True,
        )

        frames = viz.fig_dict["frames"]
        zero_frame = next(frame for frame in frames if frame["time"] == 0.0)
        past_frame = next(frame for frame in frames if frame["time"] == -1.0)

        self.assertEqual(len(zero_frame["decorations"]), 1)
        self.assertEqual(zero_frame["decorations"][0]["kind"], "milky_way_model")
        self.assertEqual(past_frame["decorations"], [])

    def test_threejs_galactic_circles_include_gc_marker_and_drop_radius_labels(self):
        viz = Animate3D(_FakeCollection(show_tracks=False), figure_theme="dark")
        viz.make_plot(
            time=np.array([0.0, -1.0]),
            renderer="threejs",
            show=False,
            galactic_mode=True,
            show_gc_line=True,
        )

        zero_frame = next(frame for frame in viz.fig_dict["frames"] if frame["time"] == 0.0)
        trace_names = [trace.get("name") for trace in zero_frame["traces"]]
        label_texts = [
            label.get("text")
            for trace in zero_frame["traces"]
            for label in trace.get("labels", [])
        ]
        gc_labels = [
            label
            for trace in zero_frame["traces"]
            if trace.get("name") == "GC"
            for label in trace.get("labels", [])
        ]

        self.assertIn("R = 4 kpc", trace_names)
        self.assertIn("R = 8.12 kpc", trace_names)
        self.assertIn("R = 12 kpc", trace_names)
        self.assertIn("GC", trace_names)
        self.assertIn("GC Ring", trace_names)
        self.assertIn("G.C.", label_texts)
        self.assertTrue(gc_labels)
        self.assertGreater(float(gc_labels[0].get("z", 0.0)), 0.0)
        self.assertNotIn("R = 4 kpc", label_texts)
        self.assertNotIn("R = 8.12 kpc", label_texts)
        self.assertNotIn("R = 12 kpc", label_texts)

    def test_threejs_text_trace_can_be_a_single_legend_item(self):
        viz = Animate3D(_FakeCollection(show_tracks=False), figure_theme="dark")
        text_trace = go.Scatter3d(
            x=[10.0, -20.0],
            y=[5.0, 15.0],
            z=[2.0, 8.0],
            mode="text",
            text=["Orion", "Sco-Cen"],
            textfont={"color": "#d7dce2", "size": 30, "family": "Helvetica"},
            name="Nearby Region Labels",
            showlegend=True,
            meta={"screen_stable_text": True, "screen_px": 28.0},
            hoverinfo="skip",
        )

        viz.make_plot(
            time=np.array([0.0, -1.0]),
            renderer="threejs",
            show=False,
            static_traces=[text_trace],
            static_traces_times=[[0.0]],
        )

        legend_item = next(
            item for item in viz.fig_dict["legend"]["items"] if item["name"] == "Nearby Region Labels"
        )
        zero_frame = next(frame for frame in viz.fig_dict["frames"] if frame["time"] == 0.0)
        label_trace = next(trace for trace in zero_frame["traces"] if trace["name"] == "Nearby Region Labels")

        self.assertTrue(legend_item["has_labels"])
        self.assertFalse(legend_item["has_points"])
        self.assertEqual([label["text"] for label in label_trace["labels"]], ["Orion", "Sco-Cen"])
        self.assertTrue(all(label["screen_stable"] for label in label_trace["labels"]))

    def test_threejs_renderer_scale_bar_is_fixed_and_not_saved(self):
        viz = Animate3D(_FakeCollection(), figure_theme="dark")
        fig = viz.make_plot(time=np.array([0.0, -1.0]), renderer="threejs", show=False)
        html = fig.to_html()

        self.assertIn('class="oviz-three-scale-bar" data-dragging="false"', html)
        self.assertIn("pointer-events: none;", html)
        self.assertIn("function onScaleBarPointerStart(event) {\n        return false;\n      }", html)
        self.assertIn("function captureScaleBarState()", html)
        self.assertIn("function captureScaleBarState() {\n        return null;\n      }", html)
        self.assertIn("scale_bar_state: captureScaleBarState()", html)

    def test_threejs_renderer_can_serialize_sky_panel_state(self):
        viz = Animate3D(_FakeCollection(), figure_theme="dark")
        viz.data_collection.cluster.df_int["x_helio"] = [10.0, 10.0]
        viz.data_collection.cluster.df_int["y_helio"] = [-5.0, -5.0]
        viz.data_collection.cluster.df_int["z_helio"] = [2.0, 2.0]

        with tempfile.TemporaryDirectory() as tmp_dir:
            members_file = Path(tmp_dir) / "members.csv"
            members_file.write_text("name,l,b\nCluster A,120.0,-20.0\nCluster A,120.1,-20.1\n", encoding="utf-8")

            fig = viz.make_plot(
                time=np.array([0.0, -1.0]),
                renderer="threejs",
                show=False,
                enable_sky_panel=True,
                sky_radius_deg=1.25,
                sky_survey="P/DSS2/color",
                cluster_members_file=str(members_file),
            )

        html = fig.to_html()
        scene_spec = viz.fig_dict
        zero_frame = next(frame for frame in scene_spec["frames"] if frame["time"] == 0.0)
        cluster_trace = next(trace for trace in zero_frame["traces"] if trace["name"] == "Cluster A")
        selection = cluster_trace["points"][0]["selection"]

        self.assertTrue(scene_spec["sky_panel"]["enabled"])
        self.assertAlmostEqual(scene_spec["sky_panel"]["radius_deg"], 1.25)
        self.assertEqual(scene_spec["sky_panel"]["survey"], "P/DSS2/color")
        self.assertIn("Cluster A", scene_spec["sky_panel"]["members_by_cluster"])
        self.assertEqual(len(scene_spec["sky_panel"]["members_by_cluster"]["Cluster A"]), 2)
        self.assertEqual(selection["cluster_name"], "member_1")
        self.assertEqual(selection["trace_name"], "Cluster A")
        self.assertTrue(np.isfinite(selection["ra_deg"]))
        self.assertTrue(np.isfinite(selection["dec_deg"]))
        self.assertIn("AladinLite/api/v3/latest/aladin.js", html)
        self.assertIn("oviz-three-sky-panel", html)
        self.assertIn('aladinOptions.projection = skyDomeBackgroundOnly ? "TAN" : "MOL";', html)
        self.assertIn("expandLayersControl: false", html)
        self.assertIn(".aladin-stack-box", html)
        self.assertIn("View: Mollweide all-sky", html)
        self.assertIn('type: "oviz-sky-hover-cluster"', html)
        self.assertIn('type: "oviz-parent-hover-cluster"', html)
        self.assertIn('aladin.on("objectHovered"', html)

    def test_threejs_renderer_can_enable_local_sky_dome_from_image_path(self):
        viz = Animate3D(_FakeCollection(), figure_theme="dark")
        with tempfile.TemporaryDirectory() as tmp_dir:
            image_path = Path(tmp_dir) / "local_sky.png"
            image_path.write_bytes(b"\x89PNG\r\n\x1a\n")

            viz.make_plot(
                time=np.array([0.0, -1.0]),
                renderer="threejs",
                show=False,
                enable_sky_panel=False,
                threejs_initial_state={
                    "sky_dome": {
                        "enabled": True,
                        "image_path": str(image_path),
                        "projection": "MOL",
                        "projection_metadata": {
                            "coordinate_frame": "galactic",
                            "lon_center_deg": 0.0,
                            "lat_center_deg": 0.0,
                        },
                        "radius_pc": 32000.0,
                        "opacity": 0.42,
                        "full_opacity_scale_bar_pc": 180.0,
                        "fade_out_scale_bar_pc": 900.0,
                        "capture_width_px": 2048,
                    },
                },
            )

        scene_spec = viz.fig_dict
        sky_dome = scene_spec["sky_dome"]
        self.assertFalse(scene_spec["sky_panel"]["enabled"])
        self.assertTrue(sky_dome["enabled"])
        self.assertEqual(sky_dome["source"], "local_image")
        self.assertTrue(sky_dome["image_data_url"].startswith("data:image/png;base64,"))
        self.assertEqual(sky_dome["projection"], "MOL")
        self.assertEqual(
            sky_dome["projection_metadata"],
            {
                "coordinate_frame": "galactic",
                "lon_center_deg": 0.0,
                "lat_center_deg": 0.0,
            },
        )
        self.assertEqual(sky_dome["radius_pc"], 32000.0)
        self.assertEqual(sky_dome["opacity"], 0.42)
        self.assertEqual(sky_dome["full_opacity_scale_bar_pc"], 180.0)
        self.assertEqual(sky_dome["fade_out_scale_bar_pc"], 900.0)
        self.assertNotIn("image_path", sky_dome)
        self.assertNotIn("capture_width_px", sky_dome)

    def test_threejs_renderer_can_enable_local_sky_dome_from_data_url(self):
        viz = Animate3D(_FakeCollection(), figure_theme="dark")
        viz.make_plot(
            time=np.array([0.0, -1.0]),
            renderer="threejs",
            show=False,
            enable_sky_panel=False,
            threejs_initial_state={
                "sky_dome_image_data_url": "data:image/jpeg;base64,abcd",
                "sky_dome_projection": "CAR",
                "sky_dome_projection_metadata": {"coordinate_frame": "icrs"},
                "sky_dome_radius_pc": 25000.0,
                "sky_dome_opacity": 0.65,
            },
        )

        sky_dome = viz.fig_dict["sky_dome"]
        self.assertTrue(sky_dome["enabled"])
        self.assertEqual(sky_dome["source"], "local_image")
        self.assertEqual(sky_dome["image_data_url"], "data:image/jpeg;base64,abcd")
        self.assertEqual(sky_dome["projection"], "CAR")
        self.assertEqual(sky_dome["projection_metadata"], {"coordinate_frame": "icrs"})
        self.assertEqual(sky_dome["radius_pc"], 25000.0)
        self.assertEqual(sky_dome["opacity"], 0.65)

    def test_threejs_renderer_exposes_embedded_sky_dome_source_selector(self):
        image_data_url = "data:image/png;base64,iVBORw0KGgo="
        scene_spec = {
            "renderer": "threejs",
            "width": 900,
            "height": 700,
            "initial_state": {},
            "theme": {},
            "center": {"x": 0.0, "y": 0.0, "z": 0.0},
            "ranges": {"x": [-1.0, 1.0], "y": [-1.0, 1.0], "z": [-1.0, 1.0]},
            "max_span": 2.0,
            "show_axes": False,
            "layout": {"scene": {}},
            "frames": [{"name": "0", "time": 0.0, "traces": []}],
            "sky_panel": {"enabled": False},
            "sky_dome": {
                "enabled": True,
                "default_source_key": "main",
                "radius_pc": 32000.0,
                "opacity": 0.5,
                "sources": [
                    {
                        "key": "main",
                        "label": "Main sky",
                        "projection": "MOL",
                        "data_url": image_data_url,
                        "main": True,
                    },
                    {
                        "key": "alternate",
                        "label": "Alternate sky",
                        "projection": "CAR",
                        "image_data_url": image_data_url,
                    },
                ],
            },
        }

        html = ThreeJSFigure(scene_spec).to_html(compress_scene_spec=False)

        self.assertIn("oviz-three-sky-source-select", html)
        self.assertIn("function normalizeSkyDomeSourceOptions", html)
        self.assertIn("function setSkyDomeSourceByKey", html)
        self.assertIn("function applyInitialSkyDomeSource", html)
        self.assertIn('startsWith("data:image/")', html)
        self.assertIn("skyDomeHasLocalSources()", html)
        self.assertIn("sky_dome_source_key", html)
        self.assertIn('"default_source_key": "main"', html)

    def test_threejs_renderer_can_enable_aladin_sky_dome_without_local_image(self):
        viz = Animate3D(_FakeCollection(), figure_theme="dark")
        fig = viz.make_plot(
            time=np.array([0.0, -1.0]),
            renderer="threejs",
            show=False,
            enable_sky_panel=True,
            sky_survey="P/DSS2/color",
            threejs_initial_state={
                "sky_dome_enabled": True,
                "sky_dome_projection": "MOL",
                "sky_dome_capture_width_px": 2048,
                "sky_dome_capture_height_px": 1024,
                "sky_dome_capture_format": "image/png",
                "sky_dome_capture_quality": 0.8,
            },
        )

        sky_dome = viz.fig_dict["sky_dome"]
        self.assertTrue(sky_dome["enabled"])
        self.assertEqual(sky_dome["source"], "aladin")
        self.assertEqual(sky_dome["projection"], "MOL")
        self.assertEqual(sky_dome["capture_width_px"], 2048)
        self.assertEqual(sky_dome["capture_height_px"], 1024)
        self.assertEqual(sky_dome["capture_format"], "image/png")
        self.assertEqual(sky_dome["capture_quality"], 0.8)
        self.assertNotIn("image_data_url", sky_dome)
        self.assertNotIn("capture_face_px", sky_dome)

        html = fig.to_html()
        self.assertIn('class="oviz-three-sky-dome-frame"', html)
        self.assertIn('type: "oviz-aladin-sky-snapshot"', html)
        self.assertIn("getViewDataURL({", html)
        self.assertIn("const skyDomeCaptureQuality = ${skyDomeCaptureQuality};", html)
        self.assertIn("quality: skyDomeCaptureQuality", html)
        self.assertIn("setSkyDomeTextureFromDataUrl", html)
        self.assertIn('aladin.setProjection(skyDomeBackgroundOnly ? "TAN" : "MOL")', html)
        self.assertIn("skyDomeUsesAladinBackground()", html)
        self.assertIn("captureOnlyMode\n          && skyDomeUsesAladinBackground()", html)
        self.assertIn('buildAladinSrcdoc([], [], "overview", null, true)', html)
        self.assertIn("galacticAxisMatrix", html)
        self.assertIn("skyDomeAxisTransformMatrix", html)
        self.assertIn("root.dataset.skyAxisTransform", html)
        self.assertIn('mode === "live_aladin"', html)
        self.assertNotIn("oviz-aladin-sky-cubemap", html)

    def test_threejs_renderer_can_use_live_aladin_sky_dome_background(self):
        viz = Animate3D(_FakeCollection(), figure_theme="dark")
        fig = viz.make_plot(
            time=np.array([0.0, -1.0]),
            renderer="threejs",
            show=False,
            enable_sky_panel=True,
            sky_survey="P/DSS2/color",
            threejs_initial_state={
                "sky_dome_enabled": True,
                "sky_dome_background_mode": "live_aladin",
            },
        )

        sky_dome = viz.fig_dict["sky_dome"]
        self.assertEqual(sky_dome["source"], "aladin")
        self.assertEqual(sky_dome["background_mode"], "live_aladin")

        html = fig.to_html()
        self.assertIn('mode === "live_aladin"', html)
        self.assertIn("captureOnlyMode\n          && skyDomeUsesAladinBackground()", html)
        self.assertIn('type: "oviz-sky-background-view"', html)
        self.assertIn("controls.enableDamping = !liveAladinSkyBackground", html)
        self.assertIn('controls.addEventListener("change"', html)
        self.assertIn("{ force: skyDomeBackgroundUserCameraActive }", html)
        self.assertIn('skyDomeSurvey = skyDomeLocalSourceName(data.survey || (skySpec && skySpec.survey) || "aladin")', html)
        self.assertIn('root.dataset.skyDomeMotion = skyDomeBackgroundUserCameraActive ? "camera-moving" : "settled"', html)
        self.assertIn("function scheduleSkyBackgroundView(data)", html)
        self.assertIn("pendingSkyBackgroundView = data", html)
        self.assertIn("function updateSkyDomeBackgroundPredictiveTransform", html)
        self.assertIn("markSkyDomeBackgroundViewApplied(data.seq)", html)
        self.assertIn('type: "oviz-aladin-sky-background-view-applied"', html)
        self.assertIn('if (skyDomeBackgroundFrameReady) {\n              setSkyDomeSnapshotStatus("loaded", "");', html)
        self.assertIn("const cameraFovDeg = Math.min(Math.max(Number(camera.fov) || 60.0, 0.05), 120.0);", html)
        self.assertIn("const horizontalFovDeg = 2.0 * Math.atan(", html)
        self.assertIn("const fovDeg = Math.min(Math.max(horizontalFovDeg, 0.05), 179.0);", html)
        self.assertNotIn("cameraFovDeg: fovDeg", html)
        self.assertIn("cameraFovDeg,", html)
        self.assertIn("const alignedFov = Math.min(Math.max(Number(skyDomeBackgroundAlignedView.fovDeg) || 1.0, 0.05), 179.0);", html)
        self.assertIn("const alignedTan = Math.max(Math.tan(THREE.MathUtils.degToRad(alignedFov * 0.5)), 1e-6);", html)
        self.assertIn("aladinInstance.setFoV(Math.min(Math.max(fovDeg, 0.05), 179.0));", html)
        self.assertIn("const scale = Math.min(Math.max(alignedTan / currentTan, 0.55), 1.9);", html)
        self.assertIn("function zoomEarthViewByWheelDelta(deltaY)", html)
        self.assertIn("function exitEarthView()", html)
        self.assertIn("function toggleEarthView()", html)
        self.assertIn("function resetToSunReferenceFrameForSkyView()", html)
        self.assertIn("resetToSunReferenceFrameForSkyView();", html)
        self.assertIn("oviz-three-earth-view-toggle", html)
        self.assertIn("earthViewReturnCameraState", html)
        self.assertIn('earthViewToggleButtonEl.textContent = isEarthView ? "Sky View" : "3D View"', html)
        self.assertIn('if (cameraViewMode === "earth")', html)
        self.assertNotIn("function exitEarthViewToFree", html)
        self.assertNotIn("function maybeEnterEarthViewAfterScrollZoom", html)
        self.assertNotIn("function enterEarthViewFromScroll", html)
        self.assertIn('if (cameraViewMode === "earth") {\n            fadeFactor = 0.0;', html)
        self.assertIn('if (cameraViewMode !== "earth") {\n          return 0.0;', html)
        self.assertIn("function animateSkyDomeViewOpacity", html)
        self.assertIn("setSkyDomeViewOpacityScale(0.0", html)
        self.assertIn("preserveDirection: true", html)
        self.assertIn("preserveFov: false", html)
        self.assertIn("function leveledSkyViewDirectionForEarthView", html)
        self.assertIn("function updateControlSensitivityForView", html)
        self.assertIn("function startSkyViewCameraDrag", html)
        self.assertIn("function rotateSkyViewCameraByPixels", html)
        self.assertIn("const xOffset = (2.0 * dx / width) * Math.tan(horizontalFovRad * 0.5);", html)
        self.assertIn(".addScaledVector(right, xOffset)", html)
        self.assertIn(".addScaledVector(up, yOffset)", html)
        self.assertIn("controls.enableRotate = !isEarthView", html)
        self.assertIn("controls.enabled = !isEarthView", html)
        self.assertIn("controls.rotateSpeed = 0.0", html)
        self.assertIn("function lockEarthViewCameraToTarget()", html)
        self.assertIn("lockEarthViewCameraToTarget();", html)
        self.assertIn('canvas.addEventListener("pointerdown", onCanvasPointerDown, { capture: true });', html)
        self.assertIn("startLockedDirection", html)
        self.assertIn("startDirection,", html)
        self.assertIn("exitTargetDirection", html)
        self.assertIn("exitStartDirection", html)
        self.assertIn("endDistance: exitTargetDistance", html)
        self.assertIn("(!galacticReferenceVisible || cameraViewMode === \"earth\")", html)
        self.assertIn("lockDirection: preserveDirection", html)
        self.assertIn("oviz-three-sky-layer-preset-select", html)
        self.assertIn("oviz-three-sky-layer-list", html)
        self.assertIn("oviz-three-sky-layer-summary", html)
        self.assertIn("oviz-three-sky-layer-body", html)
        self.assertIn("oviz-three-sky-layer-opacity", html)
        self.assertIn("oviz-three-sky-layer-stretch", html)
        self.assertIn("oviz-three-sky-layer-colormap", html)
        self.assertIn("oviz-three-sky-layer-cut-min", html)
        self.assertIn("oviz-three-sky-layer-cut-max", html)
        self.assertIn("oviz-three-sky-hips-search", html)
        self.assertIn("function moveActiveSkyLayer", html)
        self.assertIn("function visibleSkyLayers", html)
        self.assertIn("Search or paste HiPS ID", html)
        self.assertIn("function fetchSkyLayerHipsRegistry()", html)
        self.assertIn("MocServer/query", html)
        self.assertIn("MAXREC\", \"5000\"", html)
        self.assertIn("skyLayerStretchOptions", html)
        self.assertIn("skyLayerColormapOptions", html)
        self.assertIn("P/Mellinger/color", html)
        self.assertIn("Mellinger Color", html)
        self.assertIn("P/PLANCK/R2/HFI/color", html)
        self.assertIn("Planck Dust Emission Color", html)
        self.assertIn('value: "log", label: "Log10"', html)
        self.assertIn('value: "asinh", label: "Asinh"', html)
        self.assertIn('value: "grayscale", label: "Gray"', html)
        self.assertNotIn("oviz-three-sky-layer-active-select", html)
        self.assertNotIn("oviz-three-legend-sky-section", html)
        self.assertIn('type: "oviz-sky-layer-state"', html)
        self.assertIn("fallbackSkyLayerPresetOptions", html)
        self.assertIn("function postAladinHipsFavorites()", html)
        self.assertIn("aladinInstance.hipsFavorites", html)
        self.assertIn('type: "oviz-aladin-hips-favorites"', html)
        self.assertIn("setAladinDefaultSkyLayerPresets(data.favorites)", html)
        self.assertIn("function setOverlaySkyImageLayer", html)
        self.assertIn("aladinInstance.newImageSurvey(survey)", html)
        self.assertIn("aladinInstance.setOverlayImageLayer(overlaySurvey, layerName)", html)
        self.assertIn("aladinInstance.setImageSurvey(survey)", html)
        self.assertIn("applySkyLayerState({ forceTiles: false, renderLegend: false, syncControls: false });", html)
        self.assertIn("function skyDomeFrameVisualFilterCss", html)
        self.assertNotIn("Show sky dome", html)
        self.assertNotIn("Keep sky visible", html)
        self.assertIn("aladinInstance.stopAnimation()", html)
        self.assertIn("aladinInstance.gotoPosition(lDeg, bDeg)", html)
        self.assertIn("const minUpdateIntervalMs = skyDomeBackgroundUserCameraActive ? 16.0 : 50.0", html)
        self.assertIn("updateSkyDomeBackgroundFrame(\n            (typeof performance", html)

    def test_threejs_renderer_can_use_native_hips_sky_dome_background(self):
        viz = Animate3D(_FakeCollection(), figure_theme="dark")
        fig = viz.make_plot(
            time=np.array([0.0, -1.0]),
            renderer="threejs",
            show=False,
            enable_sky_panel=True,
            sky_survey="P/DSS2/color",
            threejs_initial_state={
                "sky_dome_enabled": True,
                "sky_dome_background_mode": "native_hips",
                "sky_dome_source": "hips",
                "sky_dome_hips_base_url": "https://alasky.cds.unistra.fr/DSS/DSSColor",
                "sky_dome_hips_tile_order": 9,
                "sky_dome_hips_allsky_tile_subdivisions": 18,
                "sky_dome_hips_max_active_tiles": 96,
                "sky_dome_hips_max_concurrent_tile_loads": 6,
                "sky_dome_hips_startup_preload_tiles": 24,
                "sky_dome_hips_startup_wait_ms": 1200,
                "sky_dome_hips_brightness": 2.5,
                "sky_dome_hips_contrast": 1.2,
                "sky_dome_hips_gamma": 1.4,
                "sky_dome_force_visible": True,
            },
        )

        sky_dome = viz.fig_dict["sky_dome"]
        self.assertEqual(sky_dome["source"], "hips")
        self.assertEqual(sky_dome["background_mode"], "native_hips")
        self.assertEqual(sky_dome["hips_base_url"], "https://alasky.cds.unistra.fr/DSS/DSSColor")
        self.assertEqual(sky_dome["hips_tile_order"], 9)
        self.assertEqual(sky_dome["hips_allsky_tile_subdivisions"], 18)
        self.assertEqual(sky_dome["hips_max_active_tiles"], 96)
        self.assertEqual(sky_dome["hips_max_concurrent_tile_loads"], 6)
        self.assertEqual(sky_dome["hips_startup_preload_tiles"], 24)
        self.assertEqual(sky_dome["hips_startup_wait_ms"], 1200.0)
        self.assertEqual(sky_dome["hips_brightness"], 2.5)
        self.assertEqual(sky_dome["hips_contrast"], 1.2)
        self.assertEqual(sky_dome["hips_gamma"], 1.4)
        self.assertTrue(sky_dome["force_visible"])
        self.assertNotIn("capture_width_px", sky_dome)
        self.assertNotIn("image_data_url", sky_dome)

        html = fig.to_html()
        self.assertIn("function skyDomeUsesNativeHips()", html)
        self.assertIn("function initializeNativeHipsSkyDome()", html)
        self.assertIn("function healpixNestedVectorFromXyf", html)
        self.assertIn("function healpixNestedVectorFromFaceXY", html)
        self.assertIn("function healpixNestedPixelFromVector", html)
        self.assertIn("function nativeHipsTilesForCurrentView", html)
        self.assertIn("function preloadNativeHipsStartupTiles", html)
        self.assertIn("preloadNativeHipsStartupTiles();", html)
        self.assertIn("function skyDomeWorldCenterForCurrentView", html)
        self.assertIn("return camera.position.clone();", html)
        self.assertIn("function skyDomeHipsAllskyTileSubdivisions", html)
        self.assertIn("function skyDomeHipsTileOrderForFov", html)
        self.assertIn("Norder${safeOrder}/Dir${dirIndex}/Npix${safePix}.${skyDomeHipsTileFormat()}", html)
        self.assertIn("Allsky.${skyDomeHipsTileFormat()}", html)
        self.assertIn("function pumpNativeHipsTileLoads()", html)
        self.assertIn("function skyDomeHipsBrightness()", html)
        self.assertIn("function applySkyDomeFrameVisualSettings()", html)
        self.assertIn("side: THREE.DoubleSide", html)
        self.assertIn("oviz-three-sky-layer-preset-select", html)
        self.assertIn("oviz-three-sky-layer-list", html)
        self.assertIn("oviz-three-sky-layer-summary", html)
        self.assertIn("oviz-three-sky-layer-stretch", html)
        self.assertIn("oviz-three-sky-layer-colormap", html)
        self.assertNotIn("oviz-three-sky-layer-active-select", html)
        self.assertIn("function syncSkyLayerControls()", html)
        self.assertNotIn("Show sky dome", html)
        self.assertNotIn("Keep sky visible", html)
        self.assertIn("function syncSkyDomeControls()", html)
        self.assertIn("function refreshSkyDomeControlStatus()", html)
        self.assertIn('root.dataset.skyDomeMode = "native-hips"', html)
        self.assertIn('skyDomeFrameEl.removeAttribute("srcdoc")', html)

    def test_threejs_renderer_can_use_hips2fits_sky_dome_background(self):
        viz = Animate3D(_FakeCollection(), figure_theme="dark")
        fig = viz.make_plot(
            time=np.array([0.0, -1.0]),
            renderer="threejs",
            show=False,
            enable_sky_panel=True,
            sky_survey="P/DSS2/color",
            threejs_initial_state={
                "sky_dome_enabled": True,
                "sky_dome_background_mode": "hips2fits",
                "sky_dome_source": "hips2fits",
                "sky_dome_projection": "CAR",
                "sky_dome_hips_survey": "P/DSS2/color",
                "sky_dome_hips_frame": "galactic",
                "sky_dome_hips2fits_width": 8192,
                "sky_dome_hips2fits_height": 4096,
                "sky_dome_hips2fits_preview_width": 4096,
                "sky_dome_hips2fits_preview_height": 2048,
                "sky_dome_hips2fits_medium_width": 4096,
                "sky_dome_hips2fits_medium_height": 2048,
                "sky_dome_hips2fits_projection": "CAR",
                "sky_dome_hips2fits_coordsys": "galactic",
                "sky_dome_hips2fits_format": "jpg",
                "sky_dome_hips2fits_center_frame": "galactic",
                "sky_dome_hips2fits_l_deg": 0.0,
                "sky_dome_hips2fits_b_deg": 0.0,
                "sky_dome_hips2fits_fov_deg": 360.0,
                "sky_dome_flip_y": True,
            },
        )

        sky_dome = viz.fig_dict["sky_dome"]
        self.assertEqual(sky_dome["source"], "hips2fits")
        self.assertEqual(sky_dome["background_mode"], "hips2fits")
        self.assertEqual(sky_dome["hips_survey"], "P/DSS2/color")
        self.assertEqual(sky_dome["hips_frame"], "galactic")
        self.assertEqual(sky_dome["hips2fits_width"], 8192)
        self.assertEqual(sky_dome["hips2fits_height"], 4096)
        self.assertEqual(sky_dome["hips2fits_preview_width"], 4096)
        self.assertEqual(sky_dome["hips2fits_preview_height"], 2048)
        self.assertEqual(sky_dome["hips2fits_medium_width"], 4096)
        self.assertEqual(sky_dome["hips2fits_medium_height"], 2048)
        self.assertEqual(sky_dome["hips2fits_projection"], "CAR")
        self.assertEqual(sky_dome["hips2fits_coordsys"], "galactic")
        self.assertEqual(sky_dome["hips2fits_format"], "jpg")
        self.assertEqual(sky_dome["hips2fits_center_frame"], "galactic")
        self.assertEqual(sky_dome["hips2fits_l_deg"], 0.0)
        self.assertEqual(sky_dome["hips2fits_b_deg"], 0.0)
        self.assertEqual(sky_dome["hips2fits_fov_deg"], 360.0)
        self.assertTrue(sky_dome["flip_y"])
        self.assertNotIn("image_data_url", sky_dome)

        html = fig.to_html()
        self.assertIn("function skyDomeUsesHips2Fits()", html)
        self.assertIn("function skyDomeHips2FitsUrl(width = null, height = null)", html)
        self.assertIn("function initializeHips2FitsSkyDome()", html)
        self.assertIn("function skyDomeTextureCoordinateFrameName()", html)
        self.assertIn("uniform mat3 galacticToIcrsMatrix;", html)
        self.assertIn("uniform int textureFrameMode;", html)
        self.assertIn("function skyDomeHips2FitsCenterIcrsDeg()", html)
        self.assertIn("icrsDegFromGalacticDeg(", html)
        self.assertIn("function skyDomeHips2FitsPreviewWidth()", html)
        self.assertIn("function skyDomeHips2FitsMediumWidth()", html)
        self.assertIn("loaded-preview", html)
        self.assertIn("loading-medium", html)
        self.assertIn("loading-hires", html)
        self.assertIn("hips-image-services/hips2fits", html)
        self.assertIn('root.dataset.skyDomeMode = String(modeName || "remote-image")', html)
        self.assertIn('skyDomeCaptureFrameSignature = skyDomeUsesHips2Fits() ? "hips2fits" : "native-hips";', html)

    def test_threejs_renderer_keeps_selection_metadata_without_sky_panel(self):
        viz = Animate3D(_FakeCollection(), figure_theme="dark")
        viz.data_collection.cluster.df_int["x_helio"] = [10.0, 10.0]
        viz.data_collection.cluster.df_int["y_helio"] = [-5.0, -5.0]
        viz.data_collection.cluster.df_int["z_helio"] = [2.0, 2.0]

        viz.make_plot(
            time=np.array([0.0, -1.0]),
            renderer="threejs",
            show=False,
            enable_sky_panel=False,
        )

        zero_frame = next(frame for frame in viz.fig_dict["frames"] if frame["time"] == 0.0)
        cluster_trace = next(trace for trace in zero_frame["traces"] if trace["name"] == "Cluster A")
        selection = cluster_trace["points"][0]["selection"]

        self.assertFalse(viz.fig_dict["sky_panel"]["enabled"])
        self.assertEqual(selection["trace_name"], "Cluster A")
        self.assertEqual(selection["cluster_name"], "member_1")
        self.assertTrue(np.isfinite(selection["ra_deg"]))
        self.assertTrue(np.isfinite(selection["dec_deg"]))

    def test_threejs_renderer_can_compact_repeated_frame_payloads(self):
        viz = Animate3D(_FakeCollection(), figure_theme="dark")

        fig = viz.make_plot(
            time=np.array([0.0, -1.0]),
            renderer="threejs",
            show=False,
            threejs_initial_state={"compact_payload_enabled": True},
        )

        frame_by_time = {float(frame["time"]): frame for frame in viz.fig_dict["frames"]}
        zero_point = frame_by_time[0.0]["traces"][0]["points"][0]
        past_point = frame_by_time[-1.0]["traces"][0]["points"][0]

        self.assertIn("selection", zero_point)
        self.assertIn("hovertext", zero_point)
        self.assertNotIn("selection", past_point)
        self.assertNotIn("hovertext", past_point)
        if "motion" in past_point:
            self.assertLessEqual(
                set(past_point["motion"]),
                {"key", "age_now_myr", "time_myr"},
            )
        self.assertIn('"compact_payload_enabled":true', fig.to_html())

    def test_threejs_renderer_builds_trace_catalog_without_members_file(self):
        viz = Animate3D(_FakeCollection(), figure_theme="dark")
        viz.data_collection.cluster.df_int["x_helio"] = [10.0, 20.0]
        viz.data_collection.cluster.df_int["y_helio"] = [-5.0, -10.0]
        viz.data_collection.cluster.df_int["z_helio"] = [2.0, 4.0]

        viz.make_plot(
            time=np.array([0.0, -1.0]),
            renderer="threejs",
            show=False,
            enable_sky_panel=True,
        )

        sky_panel = viz.fig_dict["sky_panel"]
        self.assertTrue(sky_panel["enabled"])
        self.assertIn("Cluster A", sky_panel["members_by_cluster"])
        self.assertEqual(len(sky_panel["members_by_cluster"]["Cluster A"]), 1)
        point = sky_panel["members_by_cluster"]["Cluster A"][0]
        self.assertTrue(np.isfinite(point["ra"]))
        self.assertTrue(np.isfinite(point["dec"]))

    def test_threejs_renderer_allows_initial_state_overrides(self):
        viz = Animate3D(_FakeCollection(), figure_theme="dark")

        fig = viz.make_plot(
            time=np.array([0.0, -1.0]),
            renderer="threejs",
            show=False,
            enable_sky_panel=True,
            threejs_initial_state={
                "lasso_volume_selection_enabled": True,
                "camera": {"view_offset": {"x": 0.28, "y": 0.0}},
                "global_controls": {
                    "camera_auto_orbit_speed": 0.45,
                    "camera_auto_orbit_direction": -1,
                },
            },
        )

        html = fig.to_html()
        self.assertTrue(viz.fig_dict["initial_state"]["lasso_volume_selection_enabled"])
        self.assertIn('"lasso_volume_selection_enabled": true', html)
        self.assertIn('"view_offset": {"x": 0.28, "y": 0.0}', html)
        self.assertIn('"camera_auto_orbit_speed": 0.45', html)
        self.assertIn('"camera_auto_orbit_direction": -1', html)
        self.assertIn("Volume lasso", html)
        self.assertNotIn("Dust lasso", html)

    def test_threejs_renderer_can_serialize_age_kde_widget(self):
        viz = Animate3D(_FakeCollection(), figure_theme="dark")
        viz.data_collection.cluster.df_int["age_myr"] = [0.5, 0.5]

        fig = viz.make_plot(
            time=np.array([0.0, -1.0]),
            renderer="threejs",
            show=False,
            show_age_kde_inset=True,
            age_kde_bandwidth_myr=0.4,
        )

        html = fig.to_html()
        age_kde = viz.fig_dict["age_kde"]

        self.assertTrue(age_kde["enabled"])
        self.assertEqual(age_kde["title"], "Relative SFH")
        self.assertEqual(age_kde["bandwidth_myr"], 0.4)
        self.assertEqual(len(age_kde["traces"]), 1)
        self.assertEqual(age_kde["traces"][0]["trace_name"], "Cluster A")
        self.assertEqual(age_kde["traces"][0]["trace_key"], "trace-0")
        self.assertIn("oviz-three-widget-select", html)
        self.assertIn("oviz-three-age-panel", html)
        self.assertIn("Age KDE", html)
        self.assertIn("Relative SFH", html)
        self.assertIn("oviz-three-age-filter-range-min", html)
        self.assertIn("Age filter for clusters contributing to the current KDE view", html)
        self.assertIn("oviz-three-widget-window-controls", html)
        self.assertIn("oviz-three-window-button-max", html)
        self.assertIn("oviz-three-window-button-min", html)
        self.assertNotIn("oviz-three-age-readout", html)
        self.assertNotIn(">Full</button>", html)
        self.assertNotIn(">Hide</button>", html)

    def test_threejs_renderer_can_serialize_dendrogram_widget(self):
        viz = Animate3D(_FakeCollection(), figure_theme="dark")
        viz.data_collection.cluster.df_int["age_myr"] = [8.0, 2.0]
        viz.data_collection.cluster.df_int["name"] = ["member_1", "member_2"]

        fig = viz.make_plot(
            time=np.array([0.0, -1.0]),
            renderer="threejs",
            show=False,
        )

        html = fig.to_html()
        dendrogram = viz.fig_dict["dendrogram"]

        self.assertTrue(dendrogram["enabled"])
        self.assertEqual(dendrogram["title"], "Birth Tree")
        self.assertEqual(dendrogram["default_trace_key"], "trace-0")
        self.assertGreaterEqual(dendrogram["traces"][0]["max_age_myr"], 8.0)
        self.assertGreaterEqual(len(dendrogram["entries"]), 2)
        self.assertIn("birth_time_myr", dendrogram["entries"][0])
        self.assertIn("x_birth", dendrogram["entries"][0])
        self.assertIn("time_samples", dendrogram["entries"][0])
        self.assertIn("oviz-three-dendrogram-panel", html)
        self.assertIn("Threshold (pc)", html)
        self.assertIn("buildDendrogramModel", html)
        self.assertIn("interpolateDendrogramPosition", html)
        self.assertIn("setDendrogramPinnedSelectionKeys", html)
        self.assertIn("traceMaxAge + 5.0", html)
        self.assertIn("Click to pin that branch highlight.", html)
        self.assertNotIn("oviz-three-dendrogram-readout", html)

    def test_threejs_renderer_can_serialize_volume_layer(self):
        viz = Animate3D(_FakeCollection(), figure_theme="dark")

        with tempfile.TemporaryDirectory() as tmp_dir:
            cube_path = Path(tmp_dir) / "dust_cube.fits"
            cube = np.linspace(0.0, 1.0, 3 * 4 * 5, dtype=np.float32).reshape(3, 4, 5)
            cube[0, 0, 0] = np.nan

            image_hdu = fits.ImageHDU(cube, name="MEAN")
            image_hdu.header["CTYPE1"] = "X"
            image_hdu.header["CTYPE2"] = "Y"
            image_hdu.header["CTYPE3"] = "Z"
            image_hdu.header["CUNIT1"] = "pc"
            image_hdu.header["CUNIT2"] = "pc"
            image_hdu.header["CUNIT3"] = "pc"
            image_hdu.header["CRVAL1"] = -4.0
            image_hdu.header["CRVAL2"] = -3.0
            image_hdu.header["CRVAL3"] = -2.0
            image_hdu.header["CRPIX1"] = 1.0
            image_hdu.header["CRPIX2"] = 1.0
            image_hdu.header["CRPIX3"] = 1.0
            image_hdu.header["CDELT1"] = 2.0
            image_hdu.header["CDELT2"] = 2.0
            image_hdu.header["CDELT3"] = 2.0
            fits.HDUList([fits.PrimaryHDU(), image_hdu]).writeto(cube_path)

            fig = viz.make_plot(
                time=np.array([0.0, -1.0]),
                renderer="threejs",
                show=False,
                enable_sky_panel=True,
                volumes=[{
                    "path": str(cube_path),
                    "name": "Dust Cube",
                    "hdu": "MEAN",
                    "max_resolution": 8,
                    "opacity": 0.22,
                    "samples": 128,
                    "alpha_coef": 42,
                    "stretch": "log10",
                    "colormap": "gist_heat_r",
                }],
            )

        html = fig.to_html()
        scene_spec = viz.fig_dict
        layer = scene_spec["volumes"]["layers"][0]
        zero_frame = next(frame for frame in scene_spec["frames"] if frame["time"] == 0.0)

        self.assertTrue(scene_spec["volumes"]["enabled"])
        self.assertEqual(layer["name"], "Dust Cube")
        self.assertEqual(layer["shape"], {"x": 5, "y": 4, "z": 3})
        self.assertEqual(layer["source_shape"], {"x": 5, "y": 4, "z": 3})
        self.assertEqual(layer["downsample_method"], "scipy_zoom")
        self.assertEqual(layer["downsample_step"], {"x": 1.0, "y": 1.0, "z": 1.0})
        self.assertEqual(layer["bounds"]["x"], [-5.0, 5.0])
        self.assertEqual(layer["bounds"]["y"], [-4.0, 4.0])
        self.assertAlmostEqual(layer["bounds"]["z"][0], -3.0)
        self.assertAlmostEqual(layer["bounds"]["z"][1], 3.0)
        self.assertEqual(layer["default_controls"]["steps"], 128)
        self.assertEqual(layer["default_controls"]["samples"], 128)
        self.assertEqual(layer["default_controls"]["alpha_coef"], 42.0)
        self.assertEqual(layer["default_controls"]["stretch"], "log10")
        self.assertAlmostEqual(layer["default_controls"]["opacity"], 0.22)
        self.assertEqual(layer["default_controls"]["colormap"], "gist_heat_r")
        self.assertTrue(layer["legend_color"])
        self.assertTrue(layer["data_b64"])
        self.assertTrue(layer["sky_overlay_data_b64"])
        self.assertEqual(layer["sky_overlay_data_encoding"], "png_atlas_uint8")
        self.assertEqual(layer["sky_overlay_atlas_tiles"]["x"], 2)
        self.assertEqual(layer["sky_overlay_atlas_tiles"]["y"], 2)
        self.assertGreaterEqual(len(layer["colormap_options"]), 3)
        self.assertIn("gist_heat", [option["name"] for option in layer["colormap_options"]])
        self.assertIn("gist_heat_r", [option["name"] for option in layer["colormap_options"]])
        self.assertEqual(zero_frame["decorations"][0]["kind"], "volume_layer")
        self.assertIn("Dust Cube", [item["name"] for item in scene_spec["legend"]["items"]])
        self.assertIn("oviz-three-volume", html)
        self.assertIn("Show volume", html)
        self.assertIn("Stretch", html)
        self.assertIn("DataTexture3D", html)
        self.assertIn("Alpha coef", html)
        self.assertIn("stretch_mode", html)
        self.assertIn("precision highp sampler3D;", html)
        self.assertIn("useSelectionPolygon", html)
        self.assertIn("selectionMaskTexture", html)
        self.assertIn("selectionViewProjectionMatrix", html)
        self.assertIn("selectionDimOutside", html)
        self.assertIn("buildVolumeSkyImageOverlaySpec", html)
        self.assertIn('CTYPE1: "GLON-CAR"', html)
        self.assertIn('CTYPE2: "GLAT-CAR"', html)
        self.assertIn("deriveVolumeSkyAxisTransform", html)
        self.assertIn("applyVolumeSkyAxisTransform", html)
        self.assertIn("const volumeSkyAxisTransform = deriveVolumeSkyAxisTransform();", html)
        self.assertIn("CDELT2: latSpan / height", html)
        self.assertIn("setOverlayImageLayer(imageLayer", html)
        self.assertIn("Rendered at t=0 only as a WebGL2 ray-marched volume.", html)
        self.assertIn("function volumeWindowStepForLayer", html)
        self.assertIn("syncVolumeWindowInput(vminInput, state.vmin", html)
        self.assertIn('vminInput.addEventListener("input"', html)
        self.assertIn('vmaxInput.addEventListener("input"', html)
        self.assertNotIn('vminInput.step = "any"', html)
        self.assertNotIn('vmaxInput.step = "any"', html)

    def test_threejs_renderer_volume_defaults_to_100_samples(self):
        viz = Animate3D(_FakeCollection(), figure_theme="dark")

        with tempfile.TemporaryDirectory() as tmp_dir:
            cube_path = Path(tmp_dir) / "default_samples_cube.fits"
            cube = np.linspace(0.0, 1.0, 2 * 3 * 4, dtype=np.float32).reshape(2, 3, 4)
            fits.HDUList([
                fits.PrimaryHDU(),
                fits.ImageHDU(cube, name="MEAN"),
            ]).writeto(cube_path)

            fig = viz.make_plot(
                time=np.array([0.0, -1.0]),
                renderer="threejs",
                show=False,
                volumes=[{
                    "path": str(cube_path),
                    "name": "Default Samples Cube",
                    "max_resolution": 8,
                }],
            )

        layer = viz.fig_dict["volumes"]["layers"][0]
        html = fig.to_html()
        self.assertEqual(layer["default_controls"]["steps"], 100)
        self.assertEqual(layer["default_controls"]["samples"], 100)
        self.assertIn('samplesInput.step = "1"', html)

    def test_threejs_renderer_inline_volume_honors_max_resolution_cap(self):
        viz = Animate3D(_FakeCollection(), figure_theme="dark")
        cube = np.linspace(0.0, 1.0, 12 * 12 * 12, dtype=np.float32).reshape(12, 12, 12)

        viz.make_plot(
            time=np.array([0.0, -1.0]),
            renderer="threejs",
            show=False,
            volumes=[{
                "data": cube,
                "name": "Capped Inline Cube",
                "max_resolution": 10,
                "max_resolution_cap": 9,
            }],
        )

        layer = viz.fig_dict["volumes"]["layers"][0]
        self.assertEqual(layer["shape"], {"x": 9, "y": 9, "z": 9})
        self.assertEqual(layer["source_shape"], {"x": 12, "y": 12, "z": 12})

    def test_threejs_renderer_volume_auto_detects_hdu_and_centered_bounds(self):
        viz = Animate3D(_FakeCollection(), figure_theme="dark")

        with tempfile.TemporaryDirectory() as tmp_dir:
            cube_path = Path(tmp_dir) / "auto_volume.fits"
            std_cube = np.ones((3, 4, 5), dtype=np.float32)
            mean_cube = np.arange(3 * 4 * 5, dtype=np.float32).reshape(3, 4, 5)
            fits.HDUList([
                fits.PrimaryHDU(),
                fits.ImageHDU(std_cube, name="STD"),
                fits.ImageHDU(mean_cube, name="MEAN"),
            ]).writeto(cube_path)

            viz.make_plot(
                time=np.array([0.0, -1.0]),
                renderer="threejs",
                show=False,
                volumes=[{
                    "path": str(cube_path),
                    "name": "Auto Volume",
                    "max_resolution": 8,
                    "opacity": 0.18,
                    "samples": 96,
                    "alpha_coef": 40,
                    "colormap": "inferno",
                }],
            )

        layer = viz.fig_dict["volumes"]["layers"][0]
        self.assertEqual(layer["hdu"], "MEAN")
        self.assertEqual(layer["bounds"]["x"], [-2.5, 2.5])
        self.assertEqual(layer["bounds"]["y"], [-2.0, 2.0])
        self.assertAlmostEqual(layer["bounds"]["z"][0], -1.5)
        self.assertAlmostEqual(layer["bounds"]["z"][1], 1.5)


if __name__ == "__main__":
    unittest.main()
