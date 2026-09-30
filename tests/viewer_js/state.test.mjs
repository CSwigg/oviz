// State model tests: legacy conversion, group presets, volume clamping.
import test from "node:test";
import assert from "node:assert/strict";

import { initialViewerState, legacySnapshotToState, clampVolumeState, groupMode } from "../../oviz/viewer/web/src/app/state.js";
import { stateDifferences, importStates, exportStatesBlock } from "../../oviz/viewer/web/src/app/states.js";

const manifest = {
  time: { values: [-2, -1, 0], initialIndex: 2 },
  animation: { fadeTimeMyr: 3, fadeInOut: false, fadeByOpacity: true },
  groups: {
    order: ["All", "Young"],
    default: "All",
    visibility: {
      All: { a: true, b: true, ring: true, dust: "legendonly" },
      Young: { a: true, b: false, ring: true, dust: "legendonly" },
    },
  },
  traces: [
    { key: "a", name: "A", color: "#ff0000", opacity: 0.7, showInLegend: true },
    { key: "b", name: "B", color: "#00ff00", opacity: 1, showInLegend: true },
    { key: "ring", name: "R = 8 kpc", color: "#888", opacity: 0.3, showInLegend: false },
  ],
  volumes: [{ key: "dust", stateKey: "dust", name: "Dust", dataRange: [0, 0.5], colormaps: ["Greys"], defaults: { vmin: 0, vmax: 0.1, opacity: 0.9, steps: 200, alpha_coef: 105, stretch: "asinh", colormap: "Greys" }, visible: true }],
  images: [],
  initialState: {
    current_group: "All",
    legend_state: { dust: true },
    global_controls: { size_points_by_stars_enabled: true },
  },
};

test("initial state honours group presets and legacy legend state", () => {
  const s = initialViewerState(manifest);
  assert.equal(s.group, "All");
  assert.equal(s.traces.a.visible, true);
  assert.equal(s.traces.b.visible, true);
  assert.equal(s.volumes.dust.visible, true); // legend_state overrides legendonly
  assert.equal(s.global.sizeByStars, true);
  assert.equal(s.global.fadeByOpacity, true);
  assert.equal(s.time.frame, 2);
  assert.equal(groupMode(manifest, "Young", "b"), false);
  assert.equal(groupMode(manifest, "Young", "missing"), false);
});

test("legacy snapshots convert: camera, time, group, styles, Sky", () => {
  const base = initialViewerState(manifest);
  const snap = {
    current_group: "Young",
    current_frame_value: 0.5,
    legend_state: { a: false, dust: false },
    trace_style_state: { a: { opacity: 0.4, sizeScale: 2, color: "#ff0000", colorMode: "fixed" } },
    volume_state_by_key: { dust: { vmax: 0.2, alphaCoef: 50 } },
    global_controls: { camera_view_mode: "earth", camera_fov: 30, point_glow_strength: 1.5 },
    camera: { position: { x: 0, y: 0, z: 0 }, target: { x: 0, y: 100, z: 0 } },
    sky_layers: [{ key: "P/DSS2/color", survey: "P/DSS2/color", opacity: 1, visible: true }],
  };
  const s = legacySnapshotToState(snap, manifest, base);
  assert.equal(s.group, "Young");
  assert.equal(s.traces.b.visible, false); // group preset
  assert.equal(s.traces.a.visible, false); // legend_state wins
  assert.equal(s.traces.a.opacity, 0.4);
  assert.equal(s.traces.a.sizeScale, 2);
  assert.equal(s.traces.a.color, undefined); // unchanged colour is not an override
  assert.equal(s.volumes.dust.visible, false);
  assert.equal(s.volumes.dust.vmax, 0.2);
  assert.equal(s.volumes.dust.alphaCoef, 50);
  assert.equal(s.global.glow, 1.5);
  assert.equal(s.view.mode, "sky");
  assert.equal(s.view.pose.distance, 0);
  assert.ok(Math.abs(s.view.pose.yaw - Math.PI / 2) < 1e-9); // looking toward l = 90°
  assert.equal(s.view.pose.fov, 30);
  assert.equal(s.time.frame, 0.5);
  assert.equal(s.sky.layers.length, 1);
});

test("volume clamping follows the legacy rules", () => {
  const spec = manifest.volumes[0];
  const c = clampVolumeState(spec, { vmin: -1, vmax: 9, opacity: 2, steps: 5000, alphaCoef: 0 });
  assert.equal(c.vmin, 0);
  assert.equal(c.vmax, 0.5);
  assert.equal(c.opacity, 1);
  assert.equal(c.steps, 768);
  assert.equal(c.alphaCoef, 1);
  const inverted = clampVolumeState(spec, { vmin: 0.3, vmax: 0.1, opacity: 1, steps: 100, alphaCoef: 10 });
  assert.equal(inverted.vmax, 0.5);
});

test("States import/export round-trip and legacy items", () => {
  const base = initialViewerState(manifest);
  const withLegacy = {
    ...manifest,
    states: {
      present_only: true,
      default_transition: { duration_ms: 900, easing: "linear" },
      items: [{ id: "x", name: "First", camera_behavior: "keep", snapshot: { current_group: "Young" } }],
    },
  };
  const project = importStates(withLegacy, base);
  assert.equal(project.presentOnly, true);
  assert.equal(project.items.length, 1);
  assert.equal(project.items[0].camera, "keep");
  assert.equal(project.items[0].state.group, "Young");
  const block = exportStatesBlock(project, { presentOnly: true, defaultMode: "present" });
  assert.equal(block.schema_version, 2);
  const again = importStates({ ...manifest, states: block }, base);
  assert.equal(again.items[0].name, "First");
  assert.deepEqual(stateDifferences(again.items[0].state, project.items[0].state), []);
});

test("stateDifferences reports numeric and structural changes", () => {
  assert.deepEqual(stateDifferences({ a: 1, b: { c: 2 } }, { a: 1, b: { c: 2 + 1e-9 } }), []);
  assert.deepEqual(stateDifferences({ a: 1 }, { a: 2 }), ["state.a"]);
  assert.deepEqual(stateDifferences({ a: "x" }, { a: "y", b: 1 }), ["state.a", "state.b"]);
});

test("legacy snapshots leave the camera anchor to be inferred on apply", () => {
  const base = initialViewerState(manifest);
  assert.deepEqual(base.view.anchor, { kind: "free" });
  base.view.anchor = { kind: "sun" };
  const s = legacySnapshotToState({ camera: { position: { x: 0, y: -10, z: 5 }, target: { x: 40, y: 0, z: 0 } } }, manifest, base);
  assert.equal(s.view.anchor, undefined);
});
