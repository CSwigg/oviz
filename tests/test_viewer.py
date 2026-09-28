"""Tests for the Oviz WebGL2 viewer: compiler, bundle, HTML writer, upgrade."""

from __future__ import annotations

import base64
import gzip
import json
import re
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

import numpy as np

from oviz import Animate3D, OvizFigure
from oviz.viewer import build
from oviz.viewer.bundle import BundleBuilder, Bundle, byte_shuffle, byte_unshuffle
from oviz.viewer.compile import compile_scene_spec
from oviz.viewer.figure import render_bundle_html
from oviz.viewer.upgrade import read_legacy_scene_spec, upgrade_html

REPO = Path(__file__).resolve().parents[1]
NODE = shutil.which("node")

try:  # Reuse the renderer tests' synthetic collection (no catalogue files).
    from tests.test_threejs_renderer import _FakeCollection
except ImportError:  # pragma: no cover - direct invocation from tests/
    from test_threejs_renderer import _FakeCollection


def _spec_with_rings():
    """A small hand-written legacy spec exercising every trace kind."""

    frames = []
    ring = [[np.cos(a) * 100, np.sin(a) * 100, 0.0] for a in np.linspace(0, 2 * np.pi, 9)]
    for f, t in enumerate([-2.0, -1.0, 0.0]):
        shift = np.array([f * 5.0, -f * 3.0, 0.0])
        segs = []
        for a, b in zip(ring[:-1], ring[1:]):
            segs.append([*(np.array(a) + shift).round(2), *(np.array(b) + shift).round(2)])
        traces = [
            {
                "key": "trace-0", "name": "Clusters", "showlegend": True, "default_color": "#ff8800",
                "default_opacity": 0.8, "default_point_size": 6.0, "has_n_stars": True,
                "points": [
                    {"x": 1.0 + f, "y": 2.0, "z": 3.0, "size": 6.0 if f else 0.0, "n_stars": 40.0,
                     "motion": {"key": "A", "age_now_myr": 1.5, "time_myr": t},
                     **({"selection": {"cluster_name": "A", "l_deg": 10.0, "b_deg": -1.0, "dist_pc": 3.7, "name_all": "A,Alpha"}} if t == 0 else {})},
                    {"x": -4.0, "y": 5.0 - f, "z": 6.0, "size": 7.0, "n_stars": 160.0, "color": "#0088ff",
                     "motion": {"key": "B", "age_now_myr": 30.0, "time_myr": t}},
                ],
            },
            {"key": "trace-1", "name": "Ring", "showlegend": False, "segments": segs,
             "line": {"color": "#94a3b8", "width": 1.0, "dash": "dash"}, "opacity": 0.4},
            {"key": "trace-2", "name": "Ring label", "showlegend": False,
             "labels": [{"text": "R = 100 pc", "x": 100.0 + shift[0], "y": shift[1], "z": 0.0, "color": "#e2e8f0", "size": 24, "screen_px": 24}]},
        ]
        if t == 0:
            traces.append({"key": "trace-3", "name": "Track", "showlegend": True, "default_color": "#aaaaaa",
                           "points": [{"x": 0.0, "y": 0.0, "z": float(i), "size": 2.0} for i in range(3)]})
        frames.append({"name": str(t), "time": t, "traces": traces, "decorations": []})
    return {
        "title": "Test </script> figure",
        "frames": frames,
        "initial_frame_index": 2,
        "center": {"x": 0, "y": 0, "z": 0},
        "max_span": 2600.0,
        "legend": {"items": [
            {"key": "trace-0", "name": "Clusters", "kind": "trace", "color": "#ff8800"},
            {"key": "trace-3", "name": "Track", "kind": "trace", "color": "#aaaaaa"},
        ]},
        "group_order": ["All"],
        "default_group": "All",
        "group_visibility": {"All": {"trace-0": True, "trace-1": True, "trace-2": True, "trace-3": "legendonly"}},
        "initial_state": {},
        "states": {"items": [{"id": "s1", "name": "Start", "snapshot": {"current_frame_index": 0}}]},
    }


class BundleTests(unittest.TestCase):
    def test_shuffle_roundtrip(self):
        data = np.arange(1000, dtype=np.float32).tobytes()
        self.assertEqual(byte_unshuffle(byte_shuffle(data, 4), 4), data)

    def test_arrays_roundtrip_and_deduplicate(self):
        b = BundleBuilder()
        values = np.linspace(-5, 5, 300, dtype=np.float32).reshape(100, 3)
        first = b.array(values, "f32", hint="pos")
        second = b.array(values.copy(), "f32", hint="again")
        self.assertEqual(first, second)
        bundle = Bundle(manifest={}, blobs=b.blobs)
        np.testing.assert_array_equal(bundle.decode(first), values)
        j = b.json({"name": ["α", "β"]})
        self.assertEqual(bundle.decode(j), {"name": ["α", "β"]})


class CompileTests(unittest.TestCase):
    def setUp(self):
        self.spec = _spec_with_rings()
        self.bundle = compile_scene_spec(self.spec)
        self.m = self.bundle.manifest
        self.traces = {t["key"]: t for t in self.m["traces"]}

    def test_time_and_camera(self):
        self.assertEqual(self.m["time"]["values"], [-2.0, -1.0, 0.0])
        self.assertEqual(self.m["time"]["initialIndex"], 2)
        self.assertEqual(self.m["time"]["zeroIndex"], 2)
        self.assertAlmostEqual(self.m["world"]["pointScale"], 4.0 / 3.0)

    def test_points_are_exact_and_identities_stable(self):
        pts = self.traces["trace-0"]["points"]
        self.assertEqual(pts["count"], 2)
        pos = self.bundle.decode(pts["position"]["blob"])
        self.assertEqual(pos.shape, (3, 2, 3))
        for f, frame in enumerate(self.spec["frames"]):
            for i, p in enumerate(frame["traces"][0]["points"]):
                np.testing.assert_allclose(pos[f, i], [p["x"], p["y"], p["z"]], rtol=0, atol=1e-5)
        sizes = self.bundle.decode(pts["size"]["blob"])
        self.assertEqual(sizes[0, 0], 0.0)  # pre-birth size carried over
        rgb = self.bundle.decode(pts["rgb"]["blob"])
        self.assertEqual(rgb[1].tolist(), [0, 136, 255])
        stars = self.bundle.decode(pts["starsFactor"]["blob"])
        np.testing.assert_allclose(stars, np.clip(np.sqrt(np.array([40, 160]) / 100.0), 0.35, 2.75), rtol=1e-6)
        meta = self.bundle.decode(pts["meta"]["blob"])
        self.assertEqual(meta["name"], ["A", "B"])
        self.assertEqual(meta["aliases"][0], "A,Alpha")
        self.assertEqual(meta["l"][0], 10.0)

    def test_rigid_lines_and_labels_use_offsets(self):
        lines = self.traces["trace-1"]["lines"]
        self.assertEqual(lines["count"], 8)
        self.assertEqual(lines["position"]["frames"], 1)
        offsets = self.bundle.decode(lines["position"]["offset"])
        np.testing.assert_allclose(offsets[2] - offsets[0], [10.0, -6.0, 0.0], atol=0.02)
        self.assertEqual(lines["dash"], "dash")
        labels = self.traces["trace-2"]["labels"]
        self.assertEqual(labels["text"], ["R = 100 pc"])
        label_pos = self.bundle.decode(labels["position"]["blob"])
        np.testing.assert_allclose(label_pos[:, 0, 0], [100.0, 105.0, 110.0])

    def test_single_frame_traces_collapse(self):
        track = self.traces["trace-3"]
        self.assertEqual(track["presence"], [2])
        self.assertEqual(track["points"]["position"]["frames"], 1)
        self.assertTrue(track["showInLegend"])
        self.assertFalse(self.traces["trace-1"]["showInLegend"])

    def test_states_pass_through(self):
        self.assertEqual(self.m["states"]["items"][0]["name"], "Start")


class CompileEdgeCaseTests(unittest.TestCase):
    def test_descending_legacy_frames_remap_state_indices(self):
        spec = _spec_with_rings()
        spec["frames"] = list(reversed(spec["frames"]))  # legacy order 0, -1, -2
        spec["initial_frame_index"] = 0
        spec["initial_state"] = {"current_frame_index": 0}
        spec["states"] = {"items": [
            {"id": "a", "name": "t=-2", "snapshot": {"current_frame_index": 2, "current_frame_value": 2.0}},
            {"id": "b", "name": "t=-0.5", "snapshot": {"current_frame_value": 0.5}},
        ]}
        m = compile_scene_spec(spec).manifest
        self.assertEqual(m["time"]["values"], [-2.0, -1.0, 0.0])
        self.assertEqual(m["time"]["initialIndex"], 2)
        self.assertEqual(m["initialState"]["current_frame_index"], 2)
        snaps = [it["snapshot"] for it in m["states"]["items"]]
        self.assertEqual(snaps[0]["current_frame_index"], 0)  # t = -2 is first after sorting
        self.assertAlmostEqual(snaps[0]["current_frame_value"], 0.0)
        self.assertAlmostEqual(snaps[1]["current_frame_value"], 1.5)  # halfway between t=0 and t=-1

    def test_points_never_collapse_to_rigid_offsets(self):
        spec = _spec_with_rings()
        for f, frame in enumerate(spec["frames"]):
            for p in frame["traces"][0]["points"]:
                p["x"] += 10.0 * f
        bundle = compile_scene_spec(spec)
        pts = next(t for t in bundle.manifest["traces"] if t["key"] == "trace-0")["points"]
        self.assertNotIn("offset", pts["position"])
        self.assertEqual(pts["position"]["frames"], 3)

    def test_html_comment_markers_in_names_stay_valid_json(self):
        spec = _spec_with_rings()
        spec["title"] = "Clusters <!-- draft --> & </script>"
        html = render_bundle_html(compile_scene_spec(spec))
        text = re.search(r'<script type="application/json" id="oviz-manifest">(.*?)</script>', html, re.S).group(1)
        self.assertNotIn("<", text)
        self.assertEqual(json.loads(text)["title"], "Clusters <!-- draft --> & </script>")


class HtmlTests(unittest.TestCase):
    def test_html_is_self_contained_and_script_safe(self):
        bundle = compile_scene_spec(_spec_with_rings())
        html = render_bundle_html(bundle, title="__CSS__ & <b>")
        self.assertIn("<title>__CSS__ &amp; &lt;b&gt;</title>", html)
        head = html.split('<script id="oviz-runtime">')[0]
        self.assertEqual(head.count('<style id="oviz-style">'), 1)
        manifest_text = re.search(r'<script type="application/json" id="oviz-manifest">(.*?)</script>', html, re.S).group(1)
        self.assertNotIn("</script", manifest_text)
        manifest = json.loads(manifest_text)
        self.assertEqual(manifest["format"], "oviz-bundle/1")
        blob_ids = re.findall(r'data-oviz-blob="([^"]+)"', head)
        self.assertEqual(sorted(blob_ids), sorted(b["id"] for b in manifest["blobs"]))
        self.assertNotIn("three.min.js", html)
        self.assertNotIn("cdn.jsdelivr", html)  # no runtime CDN dependency for 3D

    def test_browser_export_skeleton_matches_python_template(self):
        template = build.template_html()
        boot = re.search(r'(<div class="ov-boot".*?</div></div>)', template, re.S).group(1)
        export_js = (build.SRC_ROOT / "app" / "export.js").read_text()
        self.assertIn(boot, export_js)
        for token in ('id="oviz-manifest"', 'id="oviz-runtime"', 'id="oviz-style"', 'id="oviz-root"'):
            self.assertIn(token, template)
            self.assertIn(token, export_js)

    def test_figure_recompiles_when_states_are_attached_later(self):
        spec = _spec_with_rings()
        spec["states"] = {}
        fig = OvizFigure(spec)
        self.assertEqual(fig.bundle.manifest["states"], {})
        fig.scene_spec["states"] = {"items": [{"id": "late", "name": "Late", "snapshot": {}}]}
        self.assertEqual(fig.bundle.manifest["states"]["items"][0]["id"], "late")


class RuntimeBuildTests(unittest.TestCase):
    def test_bundle_has_unique_exports_and_valid_syntax(self):
        js = build.bundle_js()
        self.assertTrue(js.startswith("(() => {"))
        self.assertNotRegex(js, r"(?m)^\s*import\s+\{")
        self.assertNotRegex(js, r"(?m)^export\s")
        if not NODE:
            self.skipTest("node is not installed")
        with tempfile.NamedTemporaryFile("w", suffix=".js", delete=False) as fh:
            fh.write(js)
        try:
            r = subprocess.run([NODE, "--check", fh.name], capture_output=True, text=True)
            self.assertEqual(r.returncode, 0, r.stderr[:2000])
        finally:
            Path(fh.name).unlink(missing_ok=True)

    def test_bundler_rejects_unsupported_module_forms(self):
        with self.assertRaises(build.ViewerBuildError):
            build._strip_module("x.js", "export default 1;\n")
        with self.assertRaises(build.ViewerBuildError):
            build._strip_module("x.js", "import * as z from './z.js';\n")

    @unittest.skipUnless(NODE, "node is not installed")
    def test_node_unit_tests(self):
        files = sorted(str(p) for p in (REPO / "tests" / "viewer_js").glob("*.test.mjs"))
        r = subprocess.run([NODE, "--test", *files], capture_output=True, text=True, cwd=REPO)
        self.assertEqual(r.returncode, 0, (r.stdout + r.stderr)[-4000:])


class IntegrationTests(unittest.TestCase):
    def test_make_plot_writes_the_new_viewer(self):
        viz = Animate3D(_FakeCollection(), figure_theme="dark")
        fig = viz.make_plot(time=np.array([0.0, -1.0, -2.0]), show=False)
        self.assertIsInstance(fig, OvizFigure)
        html = fig.to_html()
        self.assertIn('id="oviz-manifest"', html)
        manifest = fig.bundle.manifest
        self.assertEqual(manifest["time"]["values"], [-2.0, -1.0, 0.0])
        names = {t["name"] for t in manifest["traces"]}
        self.assertTrue(names)

    def test_upgrade_reads_compressed_and_inline_legacy_figures(self):
        viz = Animate3D(_FakeCollection(), figure_theme="dark")
        classic = viz.make_plot(time=np.array([0.0, -1.0]), show=False, viewer="classic")
        with tempfile.TemporaryDirectory() as tmp:
            for mode in (True, False):
                src = Path(tmp) / f"classic_{mode}.html"
                classic.write_html(src, compress_scene_spec=mode)
                spec = read_legacy_scene_spec(src)
                self.assertEqual(len(spec["frames"]), 2)
                out = Path(tmp) / f"new_{mode}.html"
                report = upgrade_html(src, out)
                self.assertTrue(out.exists())
                self.assertGreater(report["output_bytes"], 0)
            # Single-element payloads (older figures) are accepted too.
            payload = base64.b64encode(gzip.compress(json.dumps(spec).encode())).decode()
            old = Path(tmp) / "old.html"
            old.write_text(
                '<script type="application/octet-stream" id="p1">' + payload + "</script>"
                '<script>const sceneSpec = await inflateOvizGzipBase64SceneSpec(readOvizSceneSpecPayload("p1"));</script>'
            )
            self.assertEqual(len(read_legacy_scene_spec(old)["frames"]), 2)


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
