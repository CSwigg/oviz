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
from html import unescape
from pathlib import Path
from unittest import mock

import numpy as np

from oviz import Animate3D, OvizFigure
from oviz.viewer import build
from oviz.viewer.bundle import BundleBuilder, Bundle, byte_shuffle, byte_unshuffle
from oviz.viewer.compile import _radec_to_gal, compile_scene_spec
from oviz.viewer.figure import _spec_fingerprint, render_bundle_html
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


def _lut_option(name, alpha=255):
    """A colormap option with a 16-entry LUT; volumes often bake alpha in."""

    ramp = np.linspace(0, 255, 16).astype(np.uint8)
    lut = np.stack([ramp, ramp // 2, 255 - ramp, np.broadcast_to(np.asarray(alpha, np.uint8), 16)], 1)
    return {"name": name, "lut_b64": base64.b64encode(lut.tobytes()).decode()}


def _member_rows(ra, dec, n=3, **extra):
    return [{"ra": ra + 0.1 * k, "dec": dec, "is_cluster_member": True, "source_id": f"{ra}-{k}", **extra}
            for k in range(n)]


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

    def test_members_without_distances_sit_at_the_parent_distance(self):
        # A member table with only a name and ra/dec is enough (CLAUDE.md):
        # such stars go on the sky at the parent cluster's present distance.
        spec = _spec_with_rings()
        spec["sky_panel"] = {"enabled": True, "members_by_cluster": {
            "A": _member_rows(83.8, -5.4),
            "B": _member_rows(56.7, 24.1, distance_pc=135.0)[:1] + _member_rows(56.9, 24.1)[:1],
        }}
        bundle = compile_scene_spec(spec)
        members = bundle.manifest["sky"]["members"]
        self.assertEqual(members["clusters"], ["A", "B"])
        self.assertEqual(bundle.decode(members["valid"]["blob"]).tolist(), [1] * 5)
        xyz = bundle.decode(members["position"]["blob"]).astype(np.float64)
        # Object A is at (3, 2, 3) pc at t = 0; B's rows keep their own median.
        np.testing.assert_allclose(np.linalg.norm(xyz, axis=1), [np.sqrt(22.0)] * 3 + [135.0] * 2, rtol=1e-5)
        lon, lat = _radec_to_gal(np.array([83.8]), np.array([-5.4]))
        l, b = np.radians(lon[0]), np.radians(lat[0])
        np.testing.assert_allclose(xyz[0] / np.sqrt(22.0), [np.cos(b) * np.cos(l), np.cos(b) * np.sin(l), np.sin(b)], atol=1e-5)
        pts = next(t for t in bundle.manifest["traces"] if t["key"] == "trace-0")["points"]
        self.assertEqual(bundle.decode(pts["members"]["blob"]).tolist(), [0, 1])

    def test_member_links_match_the_loaders_names_and_alias_separators(self):
        spec = _spec_with_rings()
        sel = spec["frames"][2]["traces"][0]["points"][0]["selection"]
        sel.update(cluster_name="Melotte_22", name_all="Melotte_22;Pleiades")
        # The loader keeps the catalogue's spelling and matches it by a
        # normalized key; the classic runtime splits aliases on , ; and |.
        for key in ("melotte 22", "MELOTTE-22", "Pleiades"):
            spec["sky_panel"] = {"enabled": True, "members_by_cluster": {key: _member_rows(56.75, 24.12, distance_pc=135.0)}}
            bundle = compile_scene_spec(spec)
            pts = next(t for t in bundle.manifest["traces"] if t["key"] == "trace-0")["points"]
            self.assertIn("members", pts, key)
            self.assertEqual(bundle.decode(pts["members"]["blob"]).tolist(), [0, -1], key)

    def test_colormap_name_clashes_follow_the_registered_name(self):
        # Two traces and a volume each define "inferno" with a different
        # table; each reference must resolve to its own table.
        ramp = np.linspace(0, 255, 16).astype(np.uint8)
        tables = {"trace-0": _lut_option("inferno"), "trace-3": _lut_option("inferno", ramp),
                  "dust": _lut_option("inferno", ramp[::-1])}
        spec = _spec_with_rings()
        for item in spec["legend"]["items"]:
            item["color_by"] = {"mode": "age", "colormap": "inferno", "colormap_options": [tables[item["key"]]]}
        spec["volumes"] = {"layers": [{
            "key": "dust", "shape": {"x": 4, "y": 4, "z": 4},
            "data_b64": base64.b64encode(np.arange(64, dtype=np.uint8).tobytes()).decode(),
            "bounds": {"x": [-10, 10], "y": [-10, 10], "z": [-10, 10]},
            "colormap_options": [tables["dust"]],
            "default_controls": {"colormap": "inferno"},
        }]}
        legacy = {
            "trace_style_state": {"trace-0": {"colormap": "inferno"}, "trace-3": {"colormap": "inferno"}},
            "volume_state_by_key": {"dust": {"colormap": "inferno"}},
        }
        spec["initial_state"] = json.loads(json.dumps(legacy))
        spec["states"] = {"items": [{"id": "s1", "name": "Start", "snapshot": json.loads(json.dumps(legacy))}]}
        bundle = compile_scene_spec(spec)
        m = bundle.manifest
        traces = {t["key"]: t for t in m["traces"]}
        self.assertEqual(sorted(m["colormaps"]), ["inferno", "inferno~1", "inferno~2"])
        vol = m["volumes"][0]
        chosen = {
            "trace-0": traces["trace-0"]["colorBy"]["colormap"],
            "trace-3": traces["trace-3"]["colorBy"]["colormap"],
            "dust": vol["defaults"]["colormap"],
        }
        self.assertEqual(traces["trace-3"]["colorBy"]["colormaps"], [chosen["trace-3"]])
        self.assertEqual(vol["colormaps"], [chosen["dust"]])
        for key, name in chosen.items():
            expected = np.frombuffer(base64.b64decode(tables[key]["lut_b64"]), np.uint8).reshape(-1, 4)
            np.testing.assert_array_equal(bundle.decode(m["colormaps"][name]["blob"]), expected, key)
        for snap in (m["initialState"], m["states"]["items"][0]["snapshot"]):
            self.assertEqual(snap["trace_style_state"]["trace-0"]["colormap"], chosen["trace-0"])
            self.assertEqual(snap["trace_style_state"]["trace-3"]["colormap"], chosen["trace-3"])
            self.assertEqual(snap["volume_state_by_key"]["dust"]["colormap"], chosen["dust"])
        self.assertEqual(spec["initial_state"], legacy)  # the scene spec itself is untouched


class LegacyVolumeCompileTests(unittest.TestCase):
    """The ccSN figure's volume kinds: gzip voxels, event KDEs, colour fields."""

    BOUNDS = {"x": [-10, 10], "y": [-10, 10], "z": [-10, 10]}

    def test_gzip_uint8_voxels_decode(self):
        vox = np.arange(64, dtype=np.uint8)
        spec = _spec_with_rings()
        spec["volumes"] = {"layers": [{
            "key": "dust", "shape": {"x": 4, "y": 4, "z": 4}, "bounds": self.BOUNDS,
            "data_encoding": "gzip_uint8",
            "data_b64": base64.b64encode(gzip.compress(vox.tobytes())).decode(),
        }]}
        bundle = compile_scene_spec(spec)
        np.testing.assert_array_equal(bundle.decode(bundle.manifest["volumes"][0]["data"]["blob"]).ravel(), vox)

    def test_event_kde_ships_sorted_events_recipe_and_float_initial_volume(self):
        events = np.array([[-1.0, 0, 0, 0], [-5.0, 1, 1, 1], [0.0, 2, 2, 2]], dtype="<f4")
        initial = np.linspace(0, 1, 64, dtype="<f4")
        spec = _spec_with_rings()
        spec["procedural_kde_sources"] = {"sn": {
            "encoding": "float32_le_interleaved_txyz", "stride": 4,
            "events_b64": base64.b64encode(events.tobytes()).decode(),
            "shape": {"x": 4, "y": 4, "z": 4}, "bounds": self.BOUNDS, "data_max": 8.0,
            "surface_rate_per_scalar_myr_kpc2": 1.6, "rolling_window_myr": 10.0,
            "initial_frame_time_myr": 0.0, "initial_volume_b64": base64.b64encode(initial.tobytes()).decode(),
        }}
        spec["volumes"] = {"layers": [{
            "key": "kde", "shape": {"x": 4, "y": 4, "z": 4}, "bounds": self.BOUNDS, "data_b64": "",
            "data_encoding": "procedural_kde_float32", "data_range": [0, 8],
            "procedural_kde_source_key": "sn", "procedural_kde_sigma_vox": 0.8,
            "procedural_kde_halo_sigma_vox": 2.4, "procedural_kde_core_weight": 0.78,
            "procedural_kde_default_rolling_window_myr": 10.0, "procedural_kde_time_sample_myr": 2.0,
            "procedural_kde_window_options_myr": [1, 5, 10], "procedural_kde_data_max": 8.0,
        }]}
        bundle = compile_scene_spec(spec)
        vol = bundle.manifest["volumes"][0]
        proc = vol["procedural"]
        self.assertEqual(proc["kind"], "event-kde")
        self.assertEqual((proc["coreSigmaVox"], proc["haloSigmaVox"], proc["coreWeight"]), (0.8, 2.4, 0.78))
        self.assertEqual((proc["windowMyr"], proc["sampleMyr"], proc["dataMax"]), (10.0, 2.0, 8.0))
        self.assertEqual(proc["windowOptions"], [1.0, 5.0, 10.0])
        np.testing.assert_array_equal(bundle.decode(proc["events"]["blob"])[:, 0], [-5.0, -1.0, 0.0])
        np.testing.assert_array_equal(bundle.decode(vol["data"]["blob"]).ravel(), initial)
        self.assertEqual(proc["initialTimeMyr"], 0.0)
        self.assertEqual(vol["colorbar"]["scale"], 1.6)
        self.assertTrue(vol["colorbar"]["perWindow"])

    def test_colour_field_note_and_segment_colours_carry_over(self):
        spec = _spec_with_rings()
        vox = base64.b64encode(np.arange(64, dtype=np.uint8).tobytes()).decode()
        spec["volumes"] = {"layers": [{
            "key": "hi", "shape": {"x": 4, "y": 4, "z": 4}, "bounds": self.BOUNDS, "data_b64": vox,
            "color_data_b64": vox, "color_range_kms": [-20, 20], "display_note": "Present day.",
            "fixed_display_encoding": True,
        }]}
        for frame in spec["frames"]:
            ring = frame["traces"][1]
            ring["segment_colors"] = ["#ff0000"] * (len(ring["segments"]) - 1) + ["#0000ff"]
        bundle = compile_scene_spec(spec)
        vol = bundle.manifest["volumes"][0]
        self.assertIn("colorData", vol)
        self.assertEqual(vol["colorbar"]["range"], [-20.0, 20.0])
        self.assertEqual((vol["note"], vol["fixedDisplay"]), ("Present day.", True))
        ring = next(t for t in bundle.manifest["traces"] if t["key"] == "trace-1")
        rgb = bundle.decode(ring["lines"]["segmentColors"]["blob"])
        np.testing.assert_array_equal(rgb[0], [255, 0, 0])
        np.testing.assert_array_equal(rgb[-1], [0, 0, 255])


class ParityCompileTests(unittest.TestCase):
    """Classic features the viewer rebuilds: widgets, Actions, focus group, Sky thumbnails."""

    def test_widget_settings_actions_and_focus_group_carry_over(self):
        spec = _spec_with_rings()
        spec["dendrogram"] = {
            "enabled": True, "title": "Birth Tree", "default_trace_key": "trace-0",
            "default_threshold_mode": "distance_pc", "default_connection_mode": "birth_to_birth",
            "default_threshold_pc": 80.0, "default_threshold_age_myr": 3.0, "entries": [{"heavy": True}],
        }
        spec["age_kde"] = {"enabled": True, "title": "Relative SFH", "bandwidth_myr": 1.5, "x_range": [-2.0, 0.0], "cluster_points": [1, 2]}
        spec["cluster_filter"] = {
            "enabled": True, "default_parameter_key": "age_now_myr", "entries": [{"heavy": True}],
            "parameters": [{"key": "age_now_myr", "min": 1.0, "max": 90.0}, {"key": "n_stars", "min": 5.0, "max": None}],
        }
        step = {"type": "time", "direction": "forward", "stop_after_frames": 1, "start": "after_previous", "delay_ms": 0, "interval_ms": 240}
        spec["actions"] = {"enabled": True, "items": [{"key": "a", "label": "Go", "description": "", "steps": [step]}]}
        spec["initial_state"] = {"global_controls": {"focus_trace_key": "trace-0"}}
        m = compile_scene_spec(spec).manifest
        self.assertEqual(m["widgets"]["birthTree"], {
            "title": "Birth Tree", "defaultTrace": "trace-0", "mode": "distance_pc",
            "links": "birth_to_birth", "thresholdPc": 80.0, "thresholdAgeMyr": 3.0,
        })
        # Only settings travel: the viewer rebuilds the data from the traces.
        self.assertEqual(m["widgets"]["ageKde"], {"title": "Relative SFH", "bandwidthMyr": 1.5, "xRange": [-2.0, 0.0]})
        # Classic filter extents, so full-range classic filters convert to none.
        self.assertEqual(m["widgets"]["clusterFilter"], {"extents": {"age_now_myr": [1.0, 90.0]}})
        self.assertEqual(m["actions"][0]["label"], "Go")
        self.assertEqual(m["actions"][0]["steps"][0]["stop_after_frames"], 1)
        self.assertEqual(m["animation"]["focusTrace"], "trace-0")
        spec["actions"]["enabled"] = False
        spec.pop("dendrogram")
        m2 = compile_scene_spec(spec).manifest
        self.assertEqual(m2["actions"], [])
        self.assertNotIn("birthTree", m2["widgets"])

    def test_longitude_grid_is_marked_present_day_only(self):
        spec = _spec_with_rings()
        for frame in spec["frames"]:
            for trace in frame["traces"]:
                if trace["key"] == "trace-1":
                    trace["name"] = "Galactic Quadrants"
        m = compile_scene_spec(spec).manifest
        flags = {t["name"]: t.get("presentDayOnly", False) for t in m["traces"]}
        self.assertTrue(flags["Galactic Quadrants"])  # drawn from today's Sun
        self.assertFalse(flags["Ring label"])  # radius circles stay at every time

    def test_sky_figures_ship_a_thumbnail_for_every_curated_survey(self):
        spec = _spec_with_rings()
        spec["sky_panel"] = {"enabled": True, "survey": "P/DSS2/color"}
        html = render_bundle_html(compile_scene_spec(spec))
        manifest = json.loads(unescape(re.search(r'id="oviz-manifest">(.*?)</script>', html, re.S).group(1)).replace("\\u003c", "<"))
        thumbs = manifest["sky"]["thumbs"]
        catalog = (REPO / "oviz/viewer/web/src/sky/catalog.js").read_text()
        ids = re.findall(r'\{ id: "([^"]+)", label:', catalog)
        self.assertGreaterEqual(len(ids), 12)
        key = lambda i: re.sub(r"[^a-z0-9.-]+", "_", re.sub(r"^cds/", "", i.lower()))
        for i in ids:
            self.assertIn(key(i), thumbs, i)
            self.assertTrue(thumbs[key(i)].startswith("data:image/jpeg;base64,"))
        # Figures without Sky carry none.
        plain = render_bundle_html(compile_scene_spec(_spec_with_rings()))
        self.assertNotIn("data:image/jpeg;base64", plain)


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

    def test_figure_recompiles_after_in_place_edits_of_large_members(self):
        spec = _spec_with_rings()
        fig = OvizFigure(spec)
        first = fig.to_html()
        self.assertIs(fig.bundle, fig.bundle)  # unchanged content reuses the bundle
        trace = spec["frames"][2]["traces"][0]
        trace["default_color"] = "#00ff00"
        trace["points"][0]["x"] = 999.0
        self.assertNotEqual(fig.to_html(), first)
        compiled = next(t for t in fig.bundle.manifest["traces"] if t["key"] == "trace-0")
        self.assertEqual(compiled["color"], "#00ff00")
        self.assertEqual(float(fig.bundle.decode(compiled["points"]["position"]["blob"])[2, 0, 0]), 999.0)
        # numpy content is fingerprinted by value; unhashable content disables the cache.
        a = _spec_fingerprint({"frames": [{"data": np.arange(4.0)}]})
        self.assertNotEqual(a, _spec_fingerprint({"frames": [{"data": np.arange(4.0) + 1}]}))
        self.assertEqual(a, _spec_fingerprint({"frames": [{"data": np.arange(4.0)}]}))
        self.assertIsNone(_spec_fingerprint({"note": object()}))

    def test_notebook_display_uses_srcdoc(self):
        # Chromium does not load data: URL documents over 2 MiB.
        fig = OvizFigure(_spec_with_rings())
        frame = fig._repr_html_()
        self.assertTrue(frame.startswith('<iframe srcdoc="'))
        self.assertNotIn("data:text/html", frame)
        srcdoc = re.match(r'<iframe srcdoc="([^"]*)"', frame).group(1)
        self.assertEqual(unescape(srcdoc), fig.to_html())
        # Isolated from the notebook's origin, as the data: URL was.
        sandbox = re.search(r'sandbox="([^"]*)"', frame).group(1).split()
        self.assertIn("allow-scripts", sandbox)
        self.assertNotIn("allow-same-origin", sandbox)


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
        with self.assertRaises(build.ViewerBuildError):
            build._strip_module("x.js", "export var z = 1;\n")
        with self.assertRaisesRegex(build.ViewerBuildError, "aliases"):
            build._strip_module("x.js", 'import { a as c } from "./m.js";\n')

    def _bundle_sources(self, files):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for rel, text in files.items():
                (root / rel).parent.mkdir(parents=True, exist_ok=True)
                (root / rel).write_text(text)
            with mock.patch.object(build, "SRC_ROOT", root):
                return build.bundle_js("app/main.js")

    def test_bundler_checks_imported_names_against_exports(self):
        math = "export function clamp(v, lo, hi) { return Math.min(Math.max(v, lo), hi); }\n"
        main = 'import { clamp, lerp } from "../core/math.js";\nexport const y = lerp(clamp(1, 0, 2));\n'
        with self.assertRaisesRegex(build.ViewerBuildError, r"app/main\.js imports `lerp` from core/math\.js"):
            self._bundle_sources({"core/math.js": math, "app/main.js": main})

    def test_bundler_leaves_template_literals_and_comments_alone(self):
        help_text = 'Usage:\nexport const figure = oviz.save();\nimport { x } from "./missing.js";\n'
        main = (
            'import { clamp } from "../core/math.js";\n'
            f"const help = `{help_text}${{clamp(5, 0, 2)}}`;\n"
            '/* export default 1;\nimport * as q from "./q.js";\n*/\n'
            "const quotes = /[`\"']/g;\n"
            "export const shown = help + quotes.source;\n"
            "console.log(JSON.stringify(shown));\n"
        )
        math = "export function clamp(v, lo, hi) { return Math.min(Math.max(v, lo), hi); }\n"
        js = self._bundle_sources({"core/math.js": math, "app/main.js": main})
        self.assertIn(help_text, js)
        self.assertIn("return { shown };", js)
        if not NODE:
            self.skipTest("node is not installed")
        r = subprocess.run([NODE, "-e", js], capture_output=True, text=True)
        self.assertEqual(r.returncode, 0, r.stderr[-2000:])
        self.assertEqual(json.loads(r.stdout), help_text + "2" + "[`\"']")

    @unittest.skipUnless(NODE, "node is not installed")
    def test_shaders_avoid_glsl_reserved_identifiers(self):
        # GLSL ES 3.0 reserves words such as `half`; declaring one fails the
        # shader at runtime only, so check every shader source statically.
        reserved = (
            "half|fixed|input|output|filter|sample|sizeof|cast|namespace|using|packed|superp|"
            "common|partition|active|external|interface|long|short|double|unsigned|goto|inline|"
            "noinline|volatile|public|static|extern|template|this|class|union|enum|typedef|asm"
        )
        decl = re.compile(r"\b(?:float|int|uint|bool|[iu]?vec[234]|mat[234])\s+(" + reserved + r")\b")
        shaders = 0
        for js in sorted((REPO / "oviz" / "viewer" / "web" / "src").rglob("*.js")):
            for source in re.findall(r"`([^`]*void main\(\)[^`]*)`", js.read_text()):
                shaders += 1
                self.assertIsNone(decl.search(source), f"reserved GLSL identifier in {js.name}")
        self.assertGreater(shaders, 5)

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

    def test_viewer_mode_is_written_and_validated(self):
        viz = Animate3D(_FakeCollection(), figure_theme="dark")
        fig = viz.make_plot(time=np.array([0.0, -1.0, -2.0]), show=False, viewer_mode="detailed")
        self.assertEqual(fig.mode, "detailed")
        html = fig.to_html()
        self.assertIn('data-oviz-mode="detailed"', html)
        manifest = json.loads(re.search(r'id="oviz-manifest">(.*?)</script>', html, re.S).group(1))
        self.assertEqual(manifest["viewer"]["mode"], "detailed")
        self.assertIn('data-oviz-mode="focus"', fig.to_html(mode="focus"))
        # Focus is the default.
        self.assertIn('data-oviz-mode="focus"', OvizFigure(bundle=fig.bundle).to_html())
        with self.assertRaises(ValueError):
            OvizFigure(bundle=fig.bundle, mode="studio")
        # Every mode the Python side accepts has a stylesheet and a runtime entry.
        css = build.bundle_css()
        runtime = build.bundle_js()
        from oviz.viewer.figure import VIEWER_MODES
        for mode in VIEWER_MODES:
            self.assertIn(f'id: "{mode}"', runtime)
            self.assertIn(f'[data-oviz-mode="{mode}"]', css)

    def test_guide_traces_get_their_own_switches_not_rows_of_the_key(self):
        guide = {
            "type": "scatter3d", "name": "Length guide", "mode": "lines+text",
            "x": [0.0, 100.0, None, 50.0], "y": [0.0, 0.0, None, 0.0], "z": [0.0, 0.0, None, -20.0],
            "text": ["", "", "", "≈100 pc"], "line": {"color": "#ffffff", "width": 1.0, "dash": "longdash"},
            "meta": {"oviz_role": "guide"},
        }
        data = {"type": "scatter3d", "name": "Cloud", "mode": "markers", "x": [1.0], "y": [2.0], "z": [3.0]}
        viz = Animate3D(_FakeCollection(), figure_theme="dark")
        fig = viz.make_plot(time=np.array([0.0, -1.0]), show=False, static_traces=[guide, data],
                            static_traces_times=[[0.0, -1.0], [0.0, -1.0]])
        items = {i["name"]: i for i in fig.scene_spec["legend"]["items"]}
        self.assertEqual(items["Length guide"]["role"], "guide")
        self.assertNotIn("role", items["Cloud"])  # other traces' scene specs are unchanged
        traces = {t["name"]: t for t in fig.bundle.manifest["traces"]}
        self.assertEqual(traces["Length guide"]["role"], "guide")
        self.assertTrue(traces["Length guide"]["showInLegend"])
        self.assertEqual(traces["Length guide"]["labels"]["text"], ["≈100 pc"])
        self.assertNotIn("role", traces["Cloud"])
        # The layers panel lists guides under Guides; the key and T leave them out.
        runtime = build.bundle_js()
        self.assertIn('t.role !== "guide"', runtime)
        self.assertIn("guideTraces()", runtime)

    def test_surface_samples_are_not_picked_and_3d_context_leaves_sky(self):
        shell = {
            "type": "scatter3d", "name": "Shell", "mode": "markers",
            "x": [100.0, -100.0], "y": [0.0, 0.0], "z": [0.0, 0.0], "hoverinfo": "skip",
            "meta": {"oviz_pickable": False, "oviz_hide_in_sky": True},
        }
        data = {"type": "scatter3d", "name": "Cloud", "mode": "markers", "x": [1.0], "y": [2.0], "z": [3.0]}
        viz = Animate3D(_FakeCollection(), figure_theme="dark")
        fig = viz.make_plot(time=np.array([0.0, -1.0]), show=False, static_traces=[shell, data],
                            static_traces_times=[[0.0, -1.0], [0.0, -1.0]])
        traces = {t["name"]: t for t in fig.bundle.manifest["traces"]}
        self.assertIs(traces["Shell"]["pickable"], False)
        self.assertIs(traces["Shell"]["skyHidden"], True)
        self.assertNotIn("pickable", traces["Cloud"])
        self.assertNotIn("skyHidden", traces["Cloud"])
        runtime = build.bundle_js()
        # Not picked (no hover, click, lasso, search or distribution filter) ...
        self.assertIn("pickable: trace.pickable !== false", runtime)
        self.assertIn("trace.pickable === false", runtime)
        # ... and, like the Sun, left out of Sky view and its key.
        self.assertIn("trace.key === this.world.sunTrace || trace.skyHidden", runtime)
        self.assertIn("!(sky && t.skyHidden)", runtime)

    def test_figures_without_a_timeline_get_a_compact_bar(self):
        css = build.bundle_css()
        # Narrow desktops and phones: the bare bar hugs its buttons.
        self.assertGreaterEqual(css.count(".ov-dock.ov-bar-only {"), 2)
        self.assertIn(".ov-dock.ov-bar-only { justify-self: center; }", css)
        self.assertIn(".ov-dock.ov-bar-only { width: auto; justify-self: center;", css)
        runtime = build.bundle_js()
        # Only figures with clusters invite a click for details.
        self.assertIn('clusters ? " · click a cluster for details" : ""', runtime)
        # Sky-start figures open in Sky (see tests/viewer_js/state.test.mjs).
        self.assertIn("this.applySkyStart();", runtime)

    def test_make_plot_adds_the_sun_by_default(self):
        import pandas as pd
        from oviz import Trace, TraceCollection

        def rows(**cols):
            base = {"name": ["A", "B"], "age_myr": [10.0, 30.0], "x": [100.0, -50.0], "y": [20.0, 80.0],
                    "z": [5.0, -10.0], "U": [-10.0, -12.0], "V": [-15.0, -14.0], "W": [-7.0, -6.0]}
            base.update(cols)
            return pd.DataFrame(base)

        time = np.array([0.0, -5.0, -10.0])
        tc = TraceCollection([Trace(rows(), "Clusters", color="#e34a4a")])
        viz = Animate3D(tc, figure_theme="dark", trace_grouping_dict={"Clusters": ["Clusters"]})
        m = viz.make_plot(time=time, show=False).bundle.manifest
        sun = next(t for t in m["traces"] if t["name"] == "Sun")
        self.assertEqual(m["world"]["sunTrace"], sun["key"])
        self.assertEqual(m["traces"][0]["name"], "Sun")  # first in the key, like the scripts' Suns
        # The Sun is in every layer group, though no group lists it.
        for group in m["groups"]["order"]:
            self.assertIs(m["groups"]["visibility"][group][sun["key"]], True)
        # The user's collection is not changed: same traces, order and indices.
        self.assertEqual([c.data_name for c in tc.clusters], ["Clusters"])
        self.assertIs(viz.data_collection, tc)
        # Never twice; not with show_sun=False; and the classic viewer is unchanged.
        again = viz.make_plot(time=time, show=False).bundle.manifest
        self.assertEqual([t["name"] for t in again["traces"]].count("Sun"), 1)
        off = viz.make_plot(time=time, show=False, show_sun=False).bundle.manifest
        self.assertNotIn("Sun", [t["name"] for t in off["traces"]])
        viz.make_plot(time=time, show=False, viewer="classic")
        self.assertNotIn("Sun", [t["name"] for t in viz.fig["frames"][0]["traces"]])
        self.assertEqual([c.data_name for c in tc.clusters], ["Clusters"])
        # Data that bring their own Sun keep theirs.
        own = Trace(rows(name=["Sun", "Sun"], z=[27.0, 27.0]).iloc[:1], "Sun", color="yellow")
        tc2 = TraceCollection([own, Trace(rows(), "Clusters")])
        m2 = Animate3D(tc2, figure_theme="dark").make_plot(time=time, show=False).bundle.manifest
        self.assertEqual([t["name"] for t in m2["traces"]].count("Sun"), 1)
        # The runtime leaves the Sun out in Sky view (the eye is at the Sun).
        runtime = build.bundle_js()
        self.assertIn("(trace.key === this.world.sunTrace || trace.skyHidden)) visible = false", runtime)

    def test_camera_anchor_defaults_to_the_lsr_and_is_validated(self):
        viz = Animate3D(_FakeCollection(), figure_theme="dark")
        time = np.array([0.0, -1.0, -2.0])
        fig = viz.make_plot(time=time, show=False)
        manifest_of = lambda html: json.loads(re.search(r'id="oviz-manifest">(.*?)</script>', html, re.S).group(1))
        # Unset, the viewer anchors the camera to the LSR itself (the home
        # view orbits the origin of the LSR-centred frame).
        viewer = manifest_of(fig.to_html())["viewer"]
        self.assertNotIn("cameraAnchor", viewer)
        self.assertNotIn("lsrOrigin", viewer)
        self.assertEqual(fig.bundle.manifest["camera"]["target"], [0.0, 0.0, 0.0])
        self.assertEqual(manifest_of(fig.to_html(camera_anchor="sun"))["viewer"]["cameraAnchor"], "sun")
        pinned = viz.make_plot(time=time, show=False, camera_anchor="LSR")
        self.assertEqual(pinned.camera_anchor, "lsr")
        self.assertEqual(manifest_of(pinned.to_html())["viewer"]["cameraAnchor"], "lsr")
        with self.assertRaises(ValueError):
            OvizFigure(bundle=fig.bundle, camera_anchor="galactic-centre")
        # A scene centred on another reference orbit has no LSR at its origin.
        custom = viz.make_plot(time=time, show=False, reference_frame_center=[0, 0, 0, -11.1, -12.24, -7.25])
        self.assertIs(manifest_of(custom.to_html())["viewer"]["lsrOrigin"], False)
        # The runtime knows every anchor the Python side accepts.
        from oviz.viewer.figure import CAMERA_ANCHORS
        runtime = build.bundle_js()
        self.assertIn("cameraAnchor", runtime)
        for anchor in CAMERA_ANCHORS:
            self.assertIn(f'"{anchor}"', runtime)

    def test_ar_model_option_points_view_in_ar_at_a_hosted_usdz(self):
        viz = Animate3D(_FakeCollection(), figure_theme="dark")
        fig = viz.make_plot(time=np.array([0.0, -1.0, -2.0]), show=False)
        manifest_of = lambda html: json.loads(re.search(r'id="oviz-manifest">(.*?)</script>', html, re.S).group(1))
        # Without a hosted model, "View in AR" builds one from the live view.
        self.assertNotIn("arModel", manifest_of(fig.to_html())["viewer"])
        hosted = OvizFigure(bundle=fig.bundle, ar_model="figure.usdz")
        self.assertEqual(manifest_of(hosted.to_html())["viewer"]["arModel"], "figure.usdz")
        self.assertEqual(manifest_of(fig.to_html(ar_model="m.usdz"))["viewer"]["arModel"], "m.usdz")
        # The runtime carries the Quick Look hand-off and the USDZ writer.
        runtime = build.bundle_js()
        self.assertIn('relList.supports("ar")', runtime)
        self.assertIn('rel: "ar"', runtime)
        self.assertIn("#usda 1.0", runtime)

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
