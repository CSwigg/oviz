// State model tests: legacy conversion, group presets, volume clamping.
import test from "node:test";
import assert from "node:assert/strict";

import { initialViewerState, legacySnapshotToState, clampVolumeState, groupMode, toggleAllPlan, soloPlan, applyLegacySkyGroup, autoOrbitRate, skyStartView } from "../../oviz/viewer/web/src/app/state.js";
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

test("sky-start figures open in Sky with their gaze, and V leads to their 3D view", () => {
  const back = { position: { x: 800, y: -1000, z: 650 }, target: { x: 0, y: 0, z: 0 }, up: { x: 0, y: 0, z: 1 }, fov: 60 };
  const sky = {
    ...manifest,
    sky: { enabled: true },
    initialState: {
      ...manifest.initialState,
      // The classic export: the camera just behind the Sun, looking toward l = 90°.
      camera: { position: { x: 0, y: -0.008, z: 0 }, target: { x: 0, y: 0, z: 0 } },
      global_controls: { camera_view_mode: "earth", camera_fov: 110, earth_view_return_camera_state: back },
    },
  };
  const start = skyStartView(sky, initialViewerState(sky));
  assert.equal(start.pose.distance, 0);
  assert.deepEqual(start.pose.target, [0, 0, 0]);
  assert.ok(Math.abs(start.pose.yaw - Math.PI / 2) < 1e-9);
  assert.equal(start.pose.fov, 110);
  assert.ok(Math.abs(start.returnPose.distance - Math.hypot(800, 1000, 650)) < 1e-6);
  assert.equal(start.returnPose.fov, 60);
  // A 3D opening view, or a figure without Sky, opens in 3D.
  assert.equal(skyStartView(manifest, initialViewerState(manifest)), null);
  const noSky = { ...sky, sky: null };
  assert.equal(skyStartView(noSky, initialViewerState(noSky)), null);
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

test("show/hide all takes volumes along and restores the mix that was on", () => {
  const items = ["trace:a", "trace:b", "volume:dust", "volume:ha"];
  // Something is on: everything goes off and the visible layers are remembered.
  const off = toggleAllPlan(items, new Set(["trace:a", "volume:dust"]), null);
  assert.deepEqual([...off.visible.values()], [false, false, false, false]);
  assert.deepEqual([...off.snapshot], ["trace:a", "volume:dust"]);
  // Showing again brings back exactly that mix: the dust volume, not Hα.
  const on = toggleAllPlan(items, new Set(), off.snapshot);
  assert.deepEqual(Object.fromEntries(on.visible), { "trace:a": true, "trace:b": false, "volume:dust": true, "volume:ha": false });
  assert.equal(on.snapshot, null);
  // Nothing remembered (or nothing of these): everything comes on.
  assert.deepEqual([...toggleAllPlan(items, new Set(), null).visible.values()], [true, true, true, true]);
  assert.deepEqual([...toggleAllPlan(["volume:ha"], new Set(), new Set(["trace:a"])).visible.values()], [true]);
  // A section restores its part of a remembered hide-all.
  const vols = toggleAllPlan(["volume:dust", "volume:ha"], new Set(), off.snapshot);
  assert.deepEqual(Object.fromEntries(vols.visible), { "volume:dust": true, "volume:ha": false });
});

test("solo hides volumes too and a second solo restores the earlier mix", () => {
  const items = ["trace:a", "trace:b", "volume:dust", "volume:ha"];
  const before = new Set(["trace:a", "trace:b", "volume:dust"]);
  const solo = soloPlan(items, before, "trace:b", null);
  assert.deepEqual(Object.fromEntries(solo.visible), { "trace:a": false, "trace:b": true, "volume:dust": false, "volume:ha": false });
  assert.deepEqual([...solo.snapshot], ["trace:a", "trace:b", "volume:dust"]);
  const back = soloPlan(items, new Set(["trace:b"]), "trace:b", solo.snapshot);
  assert.deepEqual(Object.fromEntries(back.visible), { "trace:a": true, "trace:b": true, "volume:dust": true, "volume:ha": false });
  assert.equal(back.snapshot, null);
  // Nothing remembered: soloing the lone visible layer again shows everything.
  assert.deepEqual([...soloPlan(items, new Set(["volume:ha"]), "volume:ha", null).visible.values()], [true, true, true, true]);
});

test("classic volume visibility needs both the state and the legend entry", () => {
  const m = { ...manifest, initialState: { ...manifest.initialState, legend_state: { dust: true }, volume_state_by_key: { dust: { visible: false } } } };
  assert.equal(initialViewerState(m).volumes.dust.visible, false);
  const base = initialViewerState(manifest);
  const snap = { legend_state: { dust: true }, volume_state_by_key: { dust: { visible: false, vmax: 0.2 } } };
  assert.equal(legacySnapshotToState(snap, manifest, base).volumes.dust.visible, false);
  const on = { legend_state: { dust: true }, volume_state_by_key: { dust: { visible: true } } };
  assert.equal(legacySnapshotToState(on, manifest, base).volumes.dust.visible, true);
});

test("classic sky groups draw only the active group's layers", () => {
  const groups = [
    { key: "all", label: "All", surveys: [] },
    { key: "all-sky", label: "All-Sky", surveys: ["P/Mellinger/color", "P/DSS2/color"] },
    { key: "infrared", label: "IR", surveys: ["P/2MASS/color"] },
  ];
  const layers = [
    { key: "P/2MASS/color", survey: "P/2MASS/color", visible: true, opacity: 1 },
    { key: "P/DSS2/color", survey: "P/DSS2/color", visible: true, opacity: 1 },
    { key: "mine", survey: "P/Fermi/color", group_key: "infrared", visible: true, opacity: 1 },
  ];
  const vis = (ls) => ls.filter((l) => l.visible).map((l) => l.key);
  assert.deepEqual(vis(applyLegacySkyGroup(layers, groups, "infrared")), ["P/2MASS/color", "mine"]);
  assert.deepEqual(vis(applyLegacySkyGroup(layers, groups, null)), ["P/DSS2/color"]); // classic default: all-sky
  assert.deepEqual(vis(applyLegacySkyGroup(layers, groups, "all")), ["P/2MASS/color", "P/DSS2/color", "mine"]);
  assert.deepEqual(vis(applyLegacySkyGroup(layers, [], "infrared")), ["P/2MASS/color", "P/DSS2/color", "mine"]);
  const m = { ...manifest, sky: { layers, layerGroups: groups } };
  const st = legacySnapshotToState({ sky_layers: layers, active_sky_layer_group_key: "infrared" }, m, initialViewerState(m));
  assert.deepEqual(vis(st.sky.layers), ["P/2MASS/color", "mine"]);
});

test("auto-orbit turns like classic OrbitControls: 1.2 × 2π/60 rad/s, clockwise from north", () => {
  assert.equal(autoOrbitRate({ autoOrbit: false }), 0);
  const r = autoOrbitRate({ autoOrbit: true, autoOrbitSpeed: 1, autoOrbitDirection: 1 });
  assert.ok(Math.abs(r + (2 * Math.PI / 60) * 1.2) < 1e-12); // yaw decreases: clockwise
  assert.ok(autoOrbitRate({ autoOrbit: true, autoOrbitSpeed: 0.45, autoOrbitDirection: -1 }) > 0);
  const s = initialViewerState({ ...manifest, initialState: { global_controls: { camera_auto_orbit_enabled: true, camera_auto_orbit_speed: 0.45, camera_auto_orbit_direction: -1, scroll_speed: 2.5 } } });
  assert.equal(s.global.autoOrbit, true);
  assert.equal(s.global.autoOrbitSpeed, 0.45);
  assert.equal(s.global.autoOrbitDirection, -1);
  assert.equal(s.global.navSpeed, 2.5);
});

test("action steps group into phases at each after_previous step", async () => {
  const { actionPhases } = await import("../../oviz/viewer/web/src/ui/actions.js");
  const steps = [{ type: "camera" }, { type: "time", start: "with_previous" }, { type: "legend_group", start: "after_previous" }, { type: "camera" }];
  assert.deepEqual(actionPhases(steps).map((p) => p.map((s) => s.type)), [["camera", "time"], ["legend_group"], ["camera"]]);
  assert.deepEqual(actionPhases([]), []);
});

test("classic lasso selections, filters, widgets and apertures convert", async () => {
  const { selectionKey } = await import("../../oviz/viewer/web/src/app/state.js");
  assert.equal(selectionKey("  Alpha_Per  cluster "), "alpha per cluster");
  const m = { ...manifest, widgets: { clusterFilter: { extents: { age_now_myr: [1, 11880], n_stars: [1, 1e5] } } } };
  const base = initialViewerState(m);
  const vp = Array.from({ length: 16 }, (_, i) => (i % 5 === 0 ? 1 : 0));
  const snap = {
    current_selection_mode: "lasso",
    selected_cluster_keys: ["Alpha_Per", "alpha per", "IC 2602"],
    lasso_selection_filter_enabled: true,
    lasso_selection_mask: { data_url: "data:,", view_projection_matrix: vp, polygon_ndc: [{ x: -0.5, y: -0.5 }, { x: 0.5, y: -0.5 }, { x: 0, y: 0.5 }] },
    cluster_filter_state: { parameter_key: "n_stars", ranges_by_key: { n_stars: { min: 10, max: 1000 }, age_now_myr: { min: 1, max: 11880 } } },
    widgets: { dendrogram: { mode: "normal", rect: { left: 40, top: 60, width: 500, height: 320 } }, age_kde: { mode: "hidden", rect: {} } },
    dendrogram_state: { trace_key: "a", connection_mode: "birth_to_birth", threshold_mode: "birth_age_myr", threshold_pc: 80, threshold_age_myr: 4 },
    zen_mode_enabled: true,
    current_selection: null,
    global_controls: { sky_aperture: { open: true, sky_center: { l: 80, b: -2 }, size_deg: 9, opacity: 0.8, blend: { lower_survey: "P/DSS2/color", upper_survey: "P/PLANCK/R2/HFI/color", fraction: 0.7 } } },
  };
  const s = legacySnapshotToState(snap, m, base);
  assert.deepEqual(s.ext.lasso.names, ["alpha per", "ic 2602"]);
  assert.equal(s.ext.lasso.isolate, true);
  assert.equal(s.ext.lasso.filterOff, false);
  assert.deepEqual(s.ext.lasso.mask.polygon, [[-0.5, -0.5], [0.5, -0.5], [0, 0.5]]);
  assert.equal(s.ext.lasso.mask.vp.length, 16);
  assert.equal(s.ext.filter.param, "nStars");
  assert.deepEqual(s.ext.filter.range, [1, 3]); // log10 of 10 and 1000
  assert.equal(s.ext.filter.mode, "hide");
  assert.deepEqual(s.ext.widgets, { open: { "birth-tree": { rect: { x: 40, y: 60, w: 500, h: 320 }, state: { trace: "a", mode: "birth_age_myr", links: "birth_to_birth", thresholdPc: 80, thresholdAgeMyr: 4 } } } });
  assert.equal(s.ext.ui.zen, true);
  assert.equal(s.ext.ui.selection, null);
  assert.deepEqual(s.sky.lens, { l: 80, b: -2, sizeDeg: 9, survey: "P/PLANCK/R2/HFI/color", opacity: 0.8 });

  // A full-range filter, no lasso, closed widgets and no aperture clear them.
  const off = legacySnapshotToState({
    current_selection_mode: "click",
    selected_cluster_keys: ["alpha per"],
    cluster_filter_state: { parameter_key: "age_now_myr", ranges_by_key: { age_now_myr: { min: 0, max: 20000 } } },
    widgets: { dendrogram: { mode: "hidden" } },
    global_controls: { sky_aperture: null, sky_apertures: [] },
  }, m, { ...s, ext: { ...s.ext }, sky: { ...s.sky } });
  assert.equal(off.ext.lasso, null);
  assert.equal(off.ext.filter, null);
  assert.equal(off.ext.widgets, null);
  assert.equal(off.sky.lens, undefined);
  // An age range converts unscaled; with the filter off the selection only stays.
  const age = legacySnapshotToState({
    selected_cluster_keys: ["x"], lasso_selection_filter_enabled: false,
    cluster_filter_state: { parameter_key: "age_now_myr", ranges_by_key: { age_now_myr: { min: 5, max: 50 } } },
  }, m, base);
  assert.deepEqual(age.ext.filter.range, [5, 50]);
  assert.equal(age.ext.lasso.filterOff, true);
  // Snapshots without these fields leave them as they were.
  const bare = legacySnapshotToState({ current_group: "All" }, m, s);
  assert.deepEqual(bare.ext.lasso, s.ext.lasso);
  assert.deepEqual(bare.sky.lens, s.sky.lens);
});

test("volumes are layers like the traces: group presets show and hide them", async () => {
  const { volumeGroupMode, applyGroupDefaults } = await import("../../oviz/viewer/web/src/app/state.js");
  const m = {
    ...manifest,
    groups: { order: ["All", "Dust"], default: "All", visibility: { All: { a: true, b: true, dust: "legendonly" }, Dust: { a: false, dust: true } } },
  };
  const spec = m.volumes[0];
  assert.equal(volumeGroupMode(m, "All", spec), "legendonly");
  assert.equal(volumeGroupMode(m, "Dust", spec), true);
  const st = initialViewerState({ ...m, initialState: { legend_state: { dust: true } } });
  assert.equal(st.volumes.dust.visible, true); // the opening legend state wins
  applyGroupDefaults(st, m, "All");
  assert.equal(st.volumes.dust.visible, false); // listed but hidden, like a legendonly trace
  applyGroupDefaults(st, m, "Dust");
  assert.equal(st.volumes.dust.visible, true);
  assert.equal(st.traces.a.visible, false);
  // A group listing another volume but not this one leaves it out.
  const m2 = { ...m, volumes: [...m.volumes, { ...spec, key: "ha", stateKey: "ha", name: "Ha" }], groups: { ...m.groups, visibility: { ...m.groups.visibility, Ha: { a: true, ha: true } } } };
  assert.equal(volumeGroupMode(m2, "Ha", spec), false);
  // Groups that list no volume at all leave volumes as they are.
  const bare = { ...manifest, groups: { order: ["All"], default: "All", visibility: { All: { a: true, b: false } } } };
  assert.equal(volumeGroupMode(bare, "All", spec), undefined);
  const s2 = initialViewerState(bare);
  const before = s2.volumes.dust.visible;
  applyGroupDefaults(s2, bare, "All");
  assert.equal(s2.volumes.dust.visible, before);
  // Classic States: the legend entry and the volume's own state decide, not the group preset.
  const snap = { current_group: "All", legend_state: { dust: true } };
  assert.equal(legacySnapshotToState(snap, m, initialViewerState(m)).volumes.dust.visible, true);
  assert.equal(legacySnapshotToState({ current_group: "All" }, m, initialViewerState(m)).volumes.dust.visible, false);
});
