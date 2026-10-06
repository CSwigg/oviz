"""Equivalence tests for the make_plot hot-path optimisations."""

from __future__ import annotations

import gc
import unittest

import astropy.units as u
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord

from oviz import Scene3D, Trace, TraceCollection
from oviz.point_sizes import set_cluster_point_sizes, size_easing
from oviz.viz import _galactic_to_icrs_deg, _radius_circle_xyz_pc, _rows_at_time


class SceneBuildSpeedupTests(unittest.TestCase):
    def test_galactic_to_icrs_matches_astropy(self):
        rng = np.random.default_rng(11)
        lon = rng.uniform(0, 360, 400)
        lat = np.degrees(np.arcsin(rng.uniform(-1, 1, 400)))
        ref = SkyCoord(l=lon * u.deg, b=lat * u.deg, frame="galactic").icrs
        ours = np.array([_galactic_to_icrs_deg(float(a), float(b)) for a, b in zip(lon, lat)])
        dra = ((ours[:, 0] - ref.ra.deg + 180.0) % 360.0) - 180.0
        self.assertLess(np.max(np.abs(dra * np.cos(np.radians(ref.dec.deg)))), 1e-9)
        self.assertLess(np.max(np.abs(ours[:, 1] - ref.dec.deg)), 1e-9)
        self.assertTrue(all(np.isnan(v) for v in _galactic_to_icrs_deg(float("nan"), 1.0)))

    def test_radius_circle_cache_matches_direct_transform(self):
        n = 1000
        circle = SkyCoord(
            rho=-8.122 * np.ones(n) * u.kpc,
            phi=np.linspace(-180, 180, n) * u.deg,
            z=[0.0] * n * u.pc,
            frame="galactocentric",
            representation_type="cylindrical",
        )
        x, y, z = _radius_circle_xyz_pc(8.122, "galactic")
        np.testing.assert_array_equal(x, circle.galactic.cartesian.x.value * 1000.0)
        np.testing.assert_array_equal(y, circle.galactic.cartesian.y.value * 1000.0)
        np.testing.assert_array_equal(z, circle.galactic.cartesian.z.value * 1000.0)
        with self.assertRaises(ValueError):
            x[0] = 1.0  # cached arrays are read-only

    def test_rows_at_time_matches_isclose_mask(self):
        times = np.repeat(np.arange(-5.0, 0.5, 0.5), 7)
        rng = np.random.default_rng(3)
        rng.shuffle(times)
        df = pd.DataFrame({"time": times, "v": np.arange(times.size)})
        for t in (-5.0, -2.5, 0.0, 0.0000000001, 1.0):
            expected = df[np.isclose(df["time"].to_numpy(float), t, rtol=0.0, atol=1e-9)]
            pd.testing.assert_frame_equal(_rows_at_time(df, t), expected)
        # A modified time column must not return stale rows.
        df2 = df.copy()
        df2.loc[:, "time"] = df2["time"] - 0.5
        expected = df2[np.isclose(df2["time"].to_numpy(float), -5.5, rtol=0.0, atol=1e-9)]
        pd.testing.assert_frame_equal(_rows_at_time(df2, -5.5), expected)

    def test_point_sizes_are_eased_per_cluster_in_index_order(self):
        times = [0.0, -1.0, -2.0, -3.0]
        df = pd.DataFrame({
            "name": ["b"] * 4 + ["a"] * 4 + [None],
            "time": times * 2 + [0.0],
            "age_myr": [1.5] * 4 + [30.0] * 4 + [2.0],
            "n_stars": [40] * 4 + [900] * 4 + [10],
        }, index=[10, 11, 12, 13, 20, 21, 22, 23, 30])
        sized = set_cluster_point_sizes(df, 1, 5, 2.0, False, True)

        # Rows without a name are dropped, as by groupby; the others keep their order.
        self.assertEqual(list(sized.index), [10, 11, 12, 13, 20, 21, 22, 23])
        self.assertNotIn("size", df.columns)
        for name in ("a", "b"):
            expected = size_easing(df[df["name"] == name].copy(), 1, 5, True, 2.0)["size"]
            np.testing.assert_array_equal(sized.loc[expected.index, "size"].to_numpy(), expected.to_numpy())

    def test_make_plot_restores_the_garbage_collector_and_copies_colormaps(self):
        table = pd.DataFrame({
            "x": [10.0, -40.0], "y": [5.0, 30.0], "z": [1.0, -2.0],
            "U": [-10.0, -12.0], "V": [-15.0, -14.0], "W": [-7.0, -6.0],
            "name": ["a", "b"], "age_myr": [3.0, 20.0],
        })
        enabled = gc.isenabled()
        try:
            for start_enabled in (True, False):
                gc.enable() if start_enabled else gc.disable()
                scene = Scene3D(TraceCollection([Trace(table.copy(), "Clusters", color="#e34a4a")]), figure_theme="dark")
                spec = scene.make_plot(time=np.array([0.0, -1.0]), show_sun=False).to_dict()
                self.assertEqual(gc.isenabled(), start_enabled)
        finally:
            gc.enable() if enabled else gc.disable()

        options = [trace["color_by"]["colormap_options"] for frame in spec["frames"] for trace in frame["traces"] if "color_by" in trace]
        self.assertGreaterEqual(len(options), 2)
        self.assertEqual(options[0], options[1])
        self.assertIsNot(options[0][0], options[1][0])


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
