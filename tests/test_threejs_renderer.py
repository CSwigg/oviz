import os
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

from oviz.viz import Animate3D  # noqa: E402


class _FakeCluster:
    def __init__(self, show_tracks=False):
        self.data_name = "Cluster A"
        self.show_tracks = show_tracks
        self.color = "cyan"
        self.opacity = 0.8
        self.marker_style = "circle"
        self.integrated = True
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


class ThreeJSRendererTests(unittest.TestCase):
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
        self.assertIn("Cluster A", html)
        self.assertIn("Shift+drag or use Lasso", html)
        self.assertIn("oviz-three-lasso-button", html)
        self.assertIn("Enable click select", html)
        self.assertIn("Lasso volumetric data", html)
        self.assertIn("let lassoVolumeSelectionEnabled = false;", html)
        self.assertIn("oviz-three-tools-toggle", html)
        self.assertIn("oviz-three-controls-toggle", html)
        self.assertIn("oviz-three-zen-mode", html)
        self.assertIn("oviz-three-reset-view", html)
        self.assertIn("data-zen=\"false\"", html)
        self.assertIn('toggleButton.addEventListener("dblclick"', html)
        self.assertIn("activeLegendEditorKey", html)
        self.assertIn("focusSelectionKey", html)
        self.assertIn("Camera FOV", html)
        self.assertIn("Focus group", html)
        self.assertIn("Fade time (Myr)", html)
        self.assertIn("Fade in and out", html)
        self.assertIn("Cluster Filter", html)
        self.assertIn("Dendrogram", html)
        self.assertIn("Parameter", html)
        self.assertIn("Keyboard help", html)
        self.assertIn(">Graphite<", html)
        self.assertIn(">Aurora<", html)
        self.assertIn("Keyboard controls are active as soon as the viewer loads.", html)
        self.assertIn("Shift + W / A / S / D", html)
        self.assertIn("View from Earth", html)
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
            self.assertIn("Three.js modules are loaded from a CDN", out_file.read_text())

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

    def test_threejs_galactic_circles_keep_circles_but_drop_labels(self):
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

        self.assertIn("R = 4 kpc", trace_names)
        self.assertIn("R = 8.12 kpc", trace_names)
        self.assertIn("R = 12 kpc", trace_names)
        self.assertNotIn("GC", trace_names)
        self.assertNotIn("GC", label_texts)
        self.assertNotIn("R = 4 kpc", label_texts)
        self.assertNotIn("R = 8.12 kpc", label_texts)
        self.assertNotIn("R = 12 kpc", label_texts)

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
        self.assertIn('aladinOptions.projection = "MOL";', html)
        self.assertIn("expandLayersControl: false", html)
        self.assertIn(".aladin-stack-box", html)
        self.assertIn("View: Mollweide all-sky", html)
        self.assertIn('type: "oviz-sky-hover-cluster"', html)
        self.assertIn('type: "oviz-parent-hover-cluster"', html)
        self.assertIn('aladin.on("objectHovered"', html)

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
                    "colormap": "inferno",
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
        self.assertEqual(layer["default_controls"]["colormap"], "inferno")
        self.assertTrue(layer["legend_color"])
        self.assertTrue(layer["data_b64"])
        self.assertTrue(layer["sky_overlay_data_b64"])
        self.assertEqual(layer["sky_overlay_data_encoding"], "png_atlas_uint8")
        self.assertEqual(layer["sky_overlay_atlas_tiles"]["x"], 2)
        self.assertEqual(layer["sky_overlay_atlas_tiles"]["y"], 2)
        self.assertGreaterEqual(len(layer["colormap_options"]), 3)
        self.assertEqual(zero_frame["decorations"][0]["kind"], "volume_layer")
        self.assertIn("Dust Cube", [item["name"] for item in scene_spec["legend"]["items"]])
        self.assertIn("oviz-three-volume", html)
        self.assertIn("Show at t=0", html)
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
        self.assertIn("Dust sky pixels:", html)
        self.assertIn("setOverlayImageLayer(imageLayer", html)
        self.assertIn("Rendered at t=0 only as a WebGL2 ray-marched volume.", html)

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
