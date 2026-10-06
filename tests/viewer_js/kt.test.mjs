// KT maps: volumes coloured by one of several fields, and flow-layer styles.
import test from "node:test";
import assert from "node:assert/strict";

import { initialViewerState, legacySnapshotToState, clampVolumeState, volumeColorField, volumeColormap } from "../../oviz/viewer/web/src/app/state.js";
import { flowStyle } from "../../oviz/viewer/web/src/layers/flow.js";

const kt = {
  key: "kt", stateKey: "kt", name: "KT map", dataRange: [0, 1], visible: true,
  colormaps: ["magma", "rdbu_r", "viridis"],
  defaults: { vmin: 0.1, vmax: 1, colormap: "magma" },
  colorFields: [
    { key: "velocity", label: "vlos", unit: "km/s", range: [-70, 70], colormap: "rdbu_r" },
    { key: "residual", label: "Residual", unit: "km/s", range: [-25, 25], colormap: "rdbu_r" },
  ],
  colorDefault: "velocity",
};
const manifest = {
  time: { values: [-1, 0], initialIndex: 1 },
  traces: [],
  volumes: [kt, { key: "dust", stateKey: "dust", name: "Dust", dataRange: [0, 1], colormaps: ["Greys"], defaults: { colormap: "Greys" }, visible: true }],
  images: [],
  initialState: {},
};

test("a volume with colour fields opens coloured by its default field", () => {
  const s = initialViewerState(manifest);
  assert.equal(s.volumes.kt.colorBy, "velocity");
  assert.deepEqual(s.volumes.kt.fieldColormaps, {});
  assert.equal(volumeColorField(kt, s.volumes.kt).key, "velocity");
  assert.equal(volumeColormap(kt, s.volumes.kt), "rdbu_r");
  // Plain volumes are untouched.
  assert.equal(s.volumes.dust.colorBy, undefined);
  assert.equal(volumeColorField(manifest.volumes[1], s.volumes.dust), null);
  assert.equal(volumeColormap(manifest.volumes[1], s.volumes.dust), "Greys");
  // A figure can open coloured by density.
  const byDensity = initialViewerState({ ...manifest, volumes: [{ ...kt, colorDefault: "value" }] });
  assert.equal(byDensity.volumes.kt.colorBy, "value");
});

test("each colour mode keeps its own colormap", () => {
  const vs = { ...initialViewerState(manifest).volumes.kt, fieldColormaps: { residual: "viridis" } };
  assert.equal(volumeColormap(kt, { ...vs, colorBy: "residual" }), "viridis");
  assert.equal(volumeColormap(kt, { ...vs, colorBy: "velocity" }), "rdbu_r");
  assert.equal(volumeColormap(kt, { ...vs, colorBy: "value" }), "magma");
  assert.equal(volumeColorField(kt, { ...vs, colorBy: "value" }), null);
});

test("clamping repairs unknown fields and colormaps", () => {
  const vs = initialViewerState(manifest).volumes.kt;
  const c = clampVolumeState(kt, { ...vs, colorBy: "metallicity", fieldColormaps: { residual: "nope", velocity: "viridis", gone: "magma" } });
  assert.equal(c.colorBy, "velocity");
  assert.deepEqual(c.fieldColormaps, { velocity: "viridis" });
  assert.equal(clampVolumeState(kt, { ...vs, colorBy: "value" }).colorBy, "value");
});

test("saved snapshots carry the colour mode", () => {
  const base = initialViewerState(manifest);
  const s = legacySnapshotToState({ volume_state_by_key: { kt: { colorBy: "residual", fieldColormaps: { residual: "viridis" } } } }, manifest, base);
  assert.equal(s.volumes.kt.colorBy, "residual");
  assert.deepEqual(s.volumes.kt.fieldColormaps, { residual: "viridis" });
});

test("flow styles: authored defaults under the reader's choices", () => {
  const trace = {
    key: "f", color: "#eef3ff",
    flow: { colormap: "RdBu_r", fixedColor: true, speedRange: [-70, 70], railOpacity: 0.06, animation: { rate: 1.5, period: 6, tail: 1.6 } },
  };
  const d = flowStyle(trace);
  assert.equal(d.byValue, false); // pale lines over a coloured map, like wind particles
  assert.deepEqual(d.range, [-70, 70]);
  assert.equal(d.rate, 1.5);
  assert.equal(d.period, 6);
  assert.equal(d.tail, 1.6);
  assert.equal(d.rail, 0.06);
  assert.equal(d.color, "#eef3ff");
  const r = flowStyle(trace, { colorMode: "by_value", colormap: "viridis", cmin: -20, cmax: "", flowRate: 4, flowPeriod: 0, flowTail: 3, flowRail: 2 });
  assert.equal(r.byValue, true);
  assert.equal(r.colormap, "viridis");
  assert.deepEqual(r.range, [-20, 70]); // an empty field keeps the authored end
  assert.equal(r.rate, 4);
  assert.equal(r.period, 0.05); // never zero: pulses repeat every period
  assert.equal(r.tail, 3);
  assert.equal(r.rail, 1);
  // A flow coloured by speed opens that way; a solid colour wins when chosen.
  const bySpeed = { ...trace, flow: { ...trace.flow, fixedColor: false } };
  assert.equal(flowStyle(bySpeed).byValue, true);
  assert.equal(flowStyle(bySpeed, { colorMode: "fixed", color: "#ff0000" }).color, "#ff0000");
  assert.equal(flowStyle({ ...trace, flow: { speedRange: [0, 1] } }).byValue, false); // no colormap: solid
});
