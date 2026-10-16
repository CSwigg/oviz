"""Tests for animated flow layers: the streamline tracer and its viewer plumbing."""

from __future__ import annotations

import json
import re
import unittest

import numpy as np

from oviz import Animate3D
from oviz.viewer import build
from oviz.viewer.compile import compile_scene_spec
from oviz.viewer.flow import KMS_TO_PC_PER_MYR, colormap_lut, flow_layer, trace_streamlines

try:
    from tests.test_threejs_renderer import _FakeCollection
    from tests.test_viewer import _spec_with_rings
except ImportError:  # pragma: no cover - direct invocation from tests/
    from test_threejs_renderer import _FakeCollection
    from test_viewer import _spec_with_rings


def _rotation_field(omega=0.02, half=500.0, cell=25.0, zhalf=100.0):
    """Rigid rotation about +z: v = omega * (-y, x, 0) km/s, speed = omega * r."""

    x = np.arange(-half + cell / 2, half, cell)
    z = np.arange(-zhalf + cell / 2, zhalf, cell)
    zz, yy, xx = np.meshgrid(z, x, x, indexing="ij")
    vx, vy, vz = -omega * yy, omega * xx, np.zeros_like(xx)
    r = np.hypot(xx, yy)
    support = (r > 60) & (r < half - 2 * cell)
    return (vx, vy, vz), (x, x, z), support


class TraceStreamlinesTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.velocity, cls.axes, cls.support = _rotation_field()
        cls.lines = trace_streamlines(
            cls.velocity, cls.axes, support=cls.support, spacing_pc=40.0,
            max_length_pc=600.0, max_lines=400, seed=3,
        )

    def test_lines_follow_the_field(self):
        s = self.lines
        self.assertGreater(s.count, 20)
        for line, (a, b) in zip(s.lines(), zip(s.offsets[:-1], s.offsets[1:])):
            r = np.hypot(line[:, 0], line[:, 1])
            # RK4 keeps a rigid-rotation path on its circle.
            self.assertLess(np.ptp(r), 1.0, "streamline drifted off its circle")
            np.testing.assert_allclose(line[:, 2], line[0, 2], atol=1e-3)
            # Lines run downstream: counter-clockwise about +z.
            t = np.diff(line[:, :2], axis=0)
            cross = line[:-1, 0] * t[:, 1] - line[:-1, 1] * t[:, 0]
            self.assertTrue(np.all(cross > 0))

    def test_integration_time_matches_arc_length_over_speed(self):
        s = self.lines
        omega = 0.02
        for a, b in zip(s.offsets[:-1], s.offsets[1:]):
            p = s.positions[a:b].astype(float)
            arc = np.concatenate([[0.0], np.cumsum(np.linalg.norm(np.diff(p, axis=0), axis=1))])
            r = np.hypot(p[:, 0], p[:, 1]).mean()
            expected = arc / (omega * r * KMS_TO_PC_PER_MYR)
            np.testing.assert_allclose(s.tau_myr[a:b], expected, rtol=0.02, atol=0.05)
            np.testing.assert_allclose(s.speed_kms[a:b], omega * np.hypot(p[:, 0], p[:, 1]), rtol=1e-3)

    def test_lines_are_evenly_spaced_and_stay_in_support(self):
        s = self.lines
        from scipy.spatial import cKDTree

        line_id = np.repeat(np.arange(s.count), np.diff(s.offsets))
        dist, idx = cKDTree(s.positions).query(s.positions, k=24)
        other = line_id[idx] != line_id[:, None]
        nearest = np.where(other, dist, np.inf).min(axis=1)
        nearest = nearest[np.isfinite(nearest)]
        # Half-spacing occupancy cells, one-cell neighbourhood: at least ~0.5 spacing apart.
        self.assertGreater(np.percentile(nearest, 1), 0.5 * 40.0)
        x, y, z = self.axes
        ix = np.clip(np.rint((s.positions[:, 0] - x[0]) / 25.0).astype(int), 0, len(x) - 1)
        iy = np.clip(np.rint((s.positions[:, 1] - y[0]) / 25.0).astype(int), 0, len(y) - 1)
        iz = np.clip(np.rint((s.positions[:, 2] - z[0]) / 25.0).astype(int), 0, len(z) - 1)
        self.assertTrue(self.support[iz, iy, ix].all())

    def test_edges_fade_and_columns_align(self):
        s = self.lines
        V = len(s.positions)
        for col in (s.tau_myr, s.speed_kms, s.trust, s.edge):
            self.assertEqual(col.shape, (V,))
        starts, ends = s.offsets[:-1], s.offsets[1:] - 1
        np.testing.assert_array_equal(s.edge[starts], 0.0)
        np.testing.assert_array_equal(s.edge[ends], 0.0)
        np.testing.assert_array_equal(s.tau_myr[starts], 0.0)
        self.assertTrue(np.all(s.trust == 1.0))

    def test_deterministic_for_a_seed(self):
        again = trace_streamlines(self.velocity, self.axes, support=self.support, spacing_pc=40.0,
                                  max_length_pc=600.0, max_lines=400, seed=3)
        np.testing.assert_array_equal(again.positions, self.lines.positions)
        np.testing.assert_array_equal(again.offsets, self.lines.offsets)

    def test_rejects_bad_inputs(self):
        (vx, vy, vz), axes, support = _rotation_field()
        with self.assertRaises(ValueError):
            trace_streamlines((vx, vy, vz[:, :, :-1]), axes)
        with self.assertRaises(ValueError):
            trace_streamlines((vx, vy, vz), (axes[0][::-1], axes[1], axes[2]))
        with self.assertRaises(ValueError):
            trace_streamlines((vx, vy, vz), axes, support=np.zeros_like(support))

    def test_trust_is_interpolated(self):
        (vx, vy, vz), axes, support = _rotation_field()
        trust = np.full(vx.shape, 0.25)
        s = trace_streamlines((vx, vy, vz), axes, support=support, trust=trust, spacing_pc=60.0, max_lines=30, seed=1)
        np.testing.assert_allclose(s.trust, 0.25, atol=1e-6)


class FlowLayerTests(unittest.TestCase):
    def setUp(self):
        velocity, axes, support = _rotation_field()
        self.lines = trace_streamlines(velocity, axes, support=support, spacing_pc=60.0, max_lines=40, seed=2)

    def test_colormap_luts(self):
        lut = colormap_lut("flow-ice")
        self.assertEqual(lut.shape, (256, 4))
        self.assertEqual(lut.dtype, np.uint8)
        self.assertTrue(np.all(lut[:, 3] == 255))
        custom = colormap_lut([(0.0, "#000000"), (1.0, "#ffffff")], width=3)
        np.testing.assert_array_equal(custom[:, 0], [0, 128, 255])

    def test_compiles_into_a_flow_trace(self):
        spec = _spec_with_rings()
        spec["flows"] = [
            flow_layer(self.lines, key="gas", name="Gas flow", speed_range=(1.0, 9.0)),
            flow_layer(self.lines, key="gas-fixed", name="Fixed", color="#ff00aa", visible=False),
        ]
        bundle = compile_scene_spec(spec)
        m = bundle.manifest
        traces = {t["key"]: t for t in m["traces"]}
        F = traces["gas"]["flow"]
        V = len(self.lines.positions)
        self.assertEqual(F["vertices"], V)
        self.assertEqual(F["lines"], self.lines.count)
        self.assertEqual(F["speedRange"], [1.0, 9.0])
        self.assertFalse(F["fixedColor"])
        self.assertIn(F["colormap"], m["colormaps"])
        pos = bundle.decode(F["position"]["blob"])
        attr = bundle.decode(F["attr"]["blob"])
        self.assertEqual(pos.shape, (V, 4))
        self.assertEqual(attr.shape, (V, 4))
        np.testing.assert_allclose(pos[:, :3], self.lines.positions, atol=1e-4)
        np.testing.assert_array_equal(pos[:, 3], np.repeat(np.arange(self.lines.count), np.diff(self.lines.offsets)))
        np.testing.assert_allclose(attr[:, 0], self.lines.tau_myr, rtol=1e-6)
        np.testing.assert_allclose(attr[:, 1], self.lines.speed_kms, rtol=1e-6)
        self.assertTrue(traces["gas"]["showInLegend"])
        self.assertNotIn("points", traces["gas"])
        # Flows join every trace group; a hidden flow starts hidden.
        self.assertIs(m["groups"]["visibility"]["All"]["gas"], True)
        self.assertIs(m["initialState"]["legend_state"]["gas-fixed"], False)
        self.assertNotIn("gas", m["initialState"]["legend_state"])
        self.assertTrue(traces["gas-fixed"]["flow"]["fixedColor"])
        self.assertEqual(traces["gas-fixed"]["color"], "#ff00aa")
        # 3D context: out of Sky view; shown at every time unless asked.
        self.assertTrue(traces["gas"]["skyHidden"])
        self.assertNotIn("presentDayOnly", traces["gas"])

    def test_values_colour_the_lines_and_present_day_only(self):
        spec = _spec_with_rings()
        signed = np.where(self.lines.positions[:, 0] > 0, 1.0, -1.0) * self.lines.speed_kms
        spec["flows"] = [flow_layer(self.lines, key="vlos", values=signed, value_range=(-5.0, 5.0),
                                    value_label="v_los", present_day_only=True, hide_in_sky=False)]
        bundle = compile_scene_spec(spec)
        trace = next(t for t in bundle.manifest["traces"] if t["key"] == "vlos")
        attr = bundle.decode(trace["flow"]["attr"]["blob"])
        np.testing.assert_allclose(attr[:, 1], signed, rtol=1e-6)
        self.assertEqual(trace["flow"]["speedRange"], [-5.0, 5.0])
        self.assertEqual(trace["flow"]["colorLabel"], "v_los")
        self.assertTrue(trace["presentDayOnly"])
        self.assertFalse(trace["skyHidden"])
        with self.assertRaises(ValueError):
            flow_layer(self.lines, values=signed[:-1])

    def test_bad_flow_columns_fail_loudly(self):
        spec = _spec_with_rings()
        flow = flow_layer(self.lines, key="gas")
        flow["tau_myr"] = flow["tau_myr"][:-1]
        spec["flows"] = [flow]
        with self.assertRaises(ValueError):
            compile_scene_spec(spec)

    def test_make_plot_passes_flows_to_the_viewer(self):
        viz = Animate3D(_FakeCollection(), figure_theme="dark")
        fig = viz.make_plot(time=np.array([0.0, -1.0]), show=False, flows=[flow_layer(self.lines, key="gas")])
        self.assertIn("gas", {t["key"] for t in fig.bundle.manifest["traces"] if t.get("flow")})
        html = fig.to_html()
        manifest = json.loads(re.search(r'id="oviz-manifest">(.*?)</script>', html, re.S).group(1))
        self.assertTrue(any(t.get("flow") for t in manifest["traces"]))

    def test_runtime_has_the_flow_layer_and_controls(self):
        js = build.bundle_js()
        for needle in ("class FlowLayer", "class FlowPlugin", "_stepFlow", "flowPaused", "drawArraysInstanced", "function flowStyle", "flowEditor(trace)"):
            self.assertIn(needle, js)


if __name__ == "__main__":
    unittest.main()
