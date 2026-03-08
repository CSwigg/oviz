import os
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go

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
        self.assertIn("Click to toggle traces on/off", html)
        self.assertIn("Cluster A", html)
        self.assertIn("data:text/html", repr_html)
        self.assertIn("alphaTest: 0.15", html)
        self.assertIn("width: 100vw", html)
        self.assertIn("height: 100vh", html)
        self.assertEqual(viz.fig_dict["renderer"], "threejs")
        self.assertEqual(viz.fig_dict["camera_up"], {"x": 0.0, "y": 0.0, "z": 1.0})

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
        self.assertIn("Click a cluster member at t=0 to open Aladin Lite.", html)

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

    def test_threejs_renderer_rejects_age_kde_inset(self):
        viz = Animate3D(_FakeCollection(), figure_theme="dark")

        with self.assertRaises(NotImplementedError):
            viz.make_plot(
                time=np.array([0.0, -1.0]),
                renderer="threejs",
                show=False,
                show_age_kde_inset=True,
            )


if __name__ == "__main__":
    unittest.main()
