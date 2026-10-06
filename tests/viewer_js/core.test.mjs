// Unit tests for the Oviz viewer's pure modules (run with `node --test`).
import test from "node:test";
import assert from "node:assert/strict";

import { galToIcrs, icrsToGal, m4invert, m4mul, m4perspective, m4lookAt, niceFloor, formatDistance } from "../../oviz/viewer/web/src/core/math.js";
import { Timeline } from "../../oviz/viewer/web/src/app/timeline.js";
import { cpuFramePosition, frameOffset } from "../../oviz/viewer/web/src/engine/frames.js";
import { unshuffle } from "../../oviz/viewer/web/src/core/loader.js";
import { poseFromEyeTarget, poseEye, lerpPose, makePose } from "../../oviz/viewer/web/src/engine/camera.js";
import { encodeViewHash, decodeViewHash, encodeStatePart, decodeStatePart, encodeViewsPart, decodeViewsPart } from "../../oviz/viewer/web/src/app/viewhash.js";
import { normalizeAnchor, encodeAnchor, decodeAnchor, sameAnchor, inferAnchor, LSR_POINT } from "../../oviz/viewer/web/src/engine/anchor.js";
import { fuzzyScore } from "../../oviz/viewer/web/src/ui/dom.js";
import { lassoBlend } from "../../oviz/viewer/web/src/layers/points.js";
import { lassoFade } from "../../oviz/viewer/web/src/ui/lasso.js";

const close = (a, b, eps = 1e-6) => assert.ok(Math.abs(a - b) <= eps, `${a} ≉ ${b}`);

test("galactic ↔ ICRS matches the IAU definition", () => {
  const [ra, dec] = galToIcrs(0, 0);
  close(ra, 266.40499, 1e-4);
  close(dec, -28.93617, 1e-4);
  const [ngpRa, ngpDec] = galToIcrs(0, 90);
  close(ngpRa, 192.85948, 1e-4);
  close(ngpDec, 27.12825, 1e-4);
  const [l, b] = icrsToGal(ra, dec);
  close(((l + 180) % 360) - 180, 0, 1e-5);
  close(b, 0, 1e-5);
});

test("matrix inverse round-trips a view-projection", () => {
  const proj = m4perspective(0.9, 1.6, 0.1, 1e5);
  const view = m4lookAt([100, -800, 400], [0, 0, 0], [0, 0, 1]);
  const vp = m4mul(proj, view);
  const inv = m4invert(vp);
  const id = m4mul(vp, inv);
  for (let i = 0; i < 16; i++) close(id[i], i % 5 === 0 ? 1 : 0, 1e-3);
});

test("nice numbers and distances", () => {
  assert.equal(niceFloor(733), 500);
  assert.equal(niceFloor(0.034), 0.02);
  assert.equal(formatDistance(1234), "1.23 kpc");
  assert.equal(formatDistance(56.4), "56.4 pc");
});

test("timeline maps time ↔ fractional frame on a non-uniform grid", () => {
  const tl = new Timeline([-60, -30, -10, 0], { initialIndex: 3 });
  assert.equal(tl.time, 0);
  close(tl.timeToFrame(-20), 1.5);
  close(tl.frameToTime(2.5), -5);
  assert.equal(tl.zeroFrame(), 3);
  tl.setTime(-45);
  close(tl.frame, 0.5);
  tl.step(1);
  assert.equal(tl.frame, 2); // step snaps to integer frames
});

test("timeline playback advances continuously and loops", () => {
  const tl = new Timeline([0, 1, 2, 3], { initialIndex: 0, intervalMs: 250 });
  tl.play(1);
  tl.tick(1000); // first tick only arms the clock
  tl.tick(1100); // 0.1 s at 4 frames/s
  close(tl.frame, 0.4, 1e-9);
  tl.frame = 2.9;
  tl.tick(1200); // wraps past the end
  assert.ok(tl.frame < 1 && tl.playing);
  tl.loop = false;
  tl.frame = 2.95;
  tl.tick(1300);
  assert.equal(tl.frame, 3);
  assert.equal(tl.playing, false);
});

test("CPU frame interpolation mirrors the shader (Catmull-Rom, absence, offsets)", () => {
  // 4 frames, 2 objects on straight lines; object 1 missing at frame 3.
  const N = 2, F = 4;
  const pos = new Float32Array(F * N * 3);
  for (let f = 0; f < F; f++) {
    pos.set([f * 10, 0, 0], (f * N) * 3);
    pos.set(f === 3 ? [NaN, NaN, NaN] : [0, f * 5, 0], (f * N + 1) * 3);
  }
  const a = cpuFramePosition(pos, F, N, 0, 1.5);
  close(a[0], 15, 1e-5); // Catmull-Rom is exact on linear motion
  const b = cpuFramePosition(pos, F, N, 1, 2.5);
  close(b[1], 10); // next frame absent → holds the present frame
  assert.equal(cpuFramePosition(new Float32Array([NaN, NaN, NaN]), 1, 1, 0, 0), null);
  const off = frameOffset(new Float32Array([0, 0, 0, 10, 20, 30]), 0.5);
  assert.deepEqual(off.map((x) => Math.round(x)), [5, 10, 15]);
});

test("byte un-shuffle inverts the Python shuffle", () => {
  const values = new Float32Array([1.5, -2.25, 3.125, 1e6]);
  const bytes = new Uint8Array(values.buffer);
  const shuffled = new Uint8Array(bytes.length);
  const n = values.length;
  for (let i = 0; i < n; i++) for (let b = 0; b < 4; b++) shuffled[b * n + i] = bytes[i * 4 + b];
  const back = new Float32Array(unshuffle(shuffled, 4).buffer);
  assert.deepEqual([...back], [...values]);
});

test("camera pose ↔ eye/target round-trip and log-distance lerp", () => {
  const p = poseFromEyeTarget([300, -400, 500], [10, 20, 30], 45);
  const eye = poseEye(p);
  close(eye[0], 300, 1e-6);
  close(eye[1], -400, 1e-6);
  close(eye[2], 500, 1e-6);
  const a = makePose({ distance: 10 });
  const b = makePose({ distance: 1000 });
  close(lerpPose(a, b, 0.5).distance, 100, 1e-6);
});

test("view links encode and decode", () => {
  const viewer = {
    pose: { target: [1.23, 4.56, 7.89], distance: 1500, yaw: 1.2345, pitch: -0.25, fov: 45 },
    timeline: { time: -12.5 },
    state: { view: { mode: "sky" }, group: "Clusters & Families" },
  };
  const hash = encodeViewHash(viewer);
  const v = decodeViewHash(`#${hash}`);
  assert.equal(v.time, -12.5);
  assert.equal(v.mode, "sky");
  assert.equal(v.group, "Clusters & Families");
  assert.deepEqual(v.pose.target, [1.2, 4.6, 7.9]);
  assert.equal(decodeViewHash("#nothing=1"), null);
  assert.equal(v.selection, undefined); // plain links keep no selection
});

test("view links keep layers, lasso selections with their outline, and the selected object", () => {
  const viewer = { pose: { target: [0, 0, 0], distance: 900, yaw: 0.5, pitch: 0.3, fov: 60 }, timeline: { time: 0 }, state: { view: { mode: "3d" }, group: "All" } };
  const dense = Array.from({ length: 700 }, (_, i) => i * 1.6 | 0); // most of a 1163-object trace
  const sparse = [5, 9, 300, 4000, 70000];
  // A hand-drawn loop: 240 points on a circle, which the link simplifies.
  const loop = Array.from({ length: 240 }, (_, i) => [0.3 * Math.cos(i / 240 * 2 * Math.PI) - 0.1, 0.4 * Math.sin(i / 240 * 2 * Math.PI) + 0.05]);
  const vp = [1.2, 0, 0, 0, 0, 1.8, 0.1, 0.1, 0, -0.3, -1.0002, -1, 5, -7, 899.5, 900];
  const extra = { layers: [true, false, true, true, false, false], lasso: { selection: [[1, dense], [7, sparse]], isolate: true, filterOff: false, mask: { polygon: loop, vp } }, object: { trace: 1, index: 42 } };
  const hash = encodeViewHash(viewer, extra);
  assert.match(hash, /&s=h,1:b[\w-]+,7:i[\w-]+/); // a bitmap for the dense set, gaps for the sparse one
  assert.ok(hash.length < 1200, `link is ${hash.length} characters`);
  const v = decodeViewHash(`#${hash}`);
  assert.deepEqual(v.layers, extra.layers);
  assert.equal(v.isolate, true);
  assert.equal(v.filterOff, false);
  assert.deepEqual(v.selection, [[1, [...new Set(dense)]], [7, sparse]]);
  assert.deepEqual(v.object, { trace: 1, index: 42 });
  // The outline comes back within a pixel or so (simplified, 16-bit), the camera exactly (float32).
  assert.ok(v.mask.polygon.length >= 12 && v.mask.polygon.length < 120, `${v.mask.polygon.length} outline points`);
  for (const [x, y] of v.mask.polygon) close(Math.hypot((x + 0.1) / 0.3, (y - 0.05) / 0.4), 1, 0.02);
  vp.forEach((x, i) => close(v.mask.vp[i], Math.fround(x), 1e-9));
  // A lasso around dust alone (no objects), dimming the rest with its filter off.
  const dust = decodeViewHash(`#${encodeViewHash(viewer, { lasso: { selection: [], isolate: false, filterOff: true, mask: { polygon: [[-0.5, -0.5], [0.5, -0.5], [0, 0.5]], vp } } })}`);
  assert.deepEqual(dust.selection, []);
  assert.equal(dust.isolate, false);
  assert.equal(dust.filterOff, true);
  assert.deepEqual(dust.mask.polygon.map(([x, y]) => [+x.toFixed(3), +y.toFixed(3)]), [[-0.5, -0.5], [0.5, -0.5], [0, 0.5]]);
  // Garbled parts are dropped, not fatal: bad sets one by one, bad escapes whole.
  const bad = decodeViewHash("#t=0&c=0,0,0,900,0.5,0.3,60&s=h,1:b***,x:i00,2:iAQ&k=zz&o=a:b&l=12");
  assert.deepEqual(bad.selection, [[2, [1]]]); // "AQ" is one gap of 1
  assert.equal(bad.mask, undefined);
  assert.equal(bad.object, undefined);
  assert.equal(bad.layers, undefined);
  const escapes = decodeViewHash("#t=0&c=0,0,0,900,0.5,0.3,60&s=h,1:b%E0%A4%A&k=%%");
  assert.equal(escapes.selection, undefined);
  assert.equal(escapes.mask, undefined);
  assert.equal(escapes.time, 0);
});

test("view links carry the whole State: every property comes back as sent", async () => {
  const layer = (key, visible, opacity) => ({ key, survey: key, label: key.split("/")[1], opacity, visible });
  const base = {
    version: 2, group: "All",
    view: { mode: "3d", pose: { target: [0, 0, 0], distance: 6053.098765, yaw: 1.5707963, pitch: -1.43818, fov: 60 }, anchor: { kind: "lsr" } },
    time: { frame: 60, speed: 1, playing: true, direction: 1 },
    traces: { a: { visible: true, inGroup: true }, b: { visible: true, inGroup: true, opacity: 0.8 } },
    volumes: { dust: { visible: true, vmin: 0, vmax: 0.1, opacity: 0.9, colormap: "Greys", stretch: "asinh" } },
    images: { mw: { visible: true, opacity: 1 } },
    global: { pointSize: 1, glow: 0.6, grid: true, labels: true, trails: 0, lassoVolumes: true },
    sky: { layers: [layer("P/PLANCK", false, 1), layer("P/Mellinger", true, 0.23), layer("P/DSS2", true, 1), layer("P/2MASS", false, 1)], backgroundVisible: true, members: "stars" },
    ext: { lasso: null, filter: null, notes: [], widgets: null, ui: { selection: null, theme: "dark" } },
  };
  const cur = JSON.parse(JSON.stringify(base));
  Object.assign(cur.view, { mode: "sky", pose: { target: [0, 0, 0], distance: 0, yaw: 1.2345678, pitch: 0.0123, fov: 24.5 }, returnPose: { target: [10, 20, 0], distance: 900, yaw: 0.5, pitch: 0.4, fov: 60 } });
  cur.time = { frame: 47.5, speed: 2 }; // paused: playing and direction are gone
  cur.traces.a.color = "#33ccff";
  cur.traces.b.visible = false;
  Object.assign(cur.volumes.dust, { vmax: 0.07, colormap: "magma" });
  cur.global.glow = 1.25;
  cur.sky.layers[0] = { ...cur.sky.layers[0], visible: true, opacity: 0.4, stretch: "log" };
  cur.sky.layers[1].visible = false;
  cur.sky.layers.push(layer("P/XMM/PN/color", true, 0.5)); // a survey added from the picker
  cur.sky.lens = { l: 80, b: -2, sizeDeg: 9, survey: "P/2MASS", opacity: 0.8 };
  cur.sky.members = "clusters";
  const selected = Array.from({ length: 700 }, (_, i) => i * 1.6 | 0);
  const loop = Array.from({ length: 200 }, (_, i) => [0.3 * Math.cos(i / 200 * 2 * Math.PI), 0.4 * Math.sin(i / 200 * 2 * Math.PI)]);
  cur.ext.lasso = { selection: { "trace-3": selected, "trace-6": [0] }, isolate: true, filterOff: false, mask: { polygon: loop, vp: Array.from({ length: 16 }, (_, i) => (i % 5 ? 0.0125 * i : 1.0001)) } };
  cur.ext.filter = { param: "nStars", range: [1, 3], mode: "hide" };
  cur.ext.notes = [{ id: "n1", text: "Orion", anchor: { pos: [100, -20, 5] }, color: "#ffffff" }];
  cur.ext.ui = { selection: { trace: "trace-3", index: 42 }, theme: "dark", layers: true };

  const part = await encodeStatePart(cur, base);
  assert.match(part, /^1z[\w-]+$/);
  assert.ok(part.length < 1100, `state part is ${part.length} characters`);
  const got = await decodeStatePart(part, base);
  // Equal up to 7 significant digits, except the outline (simplified to about a pixel).
  const lassoGot = got.ext.lasso, lassoCur = cur.ext.lasso;
  for (const [x, y] of lassoGot.mask.polygon) close(Math.hypot(x / 0.3, y / 0.4), 1, 0.02);
  delete lassoGot.mask.polygon;
  const want = JSON.parse(JSON.stringify(cur));
  delete want.ext.lasso.mask.polygon;
  want.ext.lasso.selection["trace-3"] = [...new Set(lassoCur.selection["trace-3"])];
  const near = (a, b, path = "") => {
    if (typeof a === "number" && typeof b === "number") return assert.ok(Math.abs(a - b) <= 1e-6 * Math.max(1, Math.abs(b)), `${path}: ${a} vs ${b}`);
    if (Array.isArray(b)) { assert.ok(Array.isArray(a), `${path} is an array`); assert.equal(a.length, b.length, `${path} length`); return b.forEach((x, i) => near(a[i], x, `${path}[${i}]`)); }
    if (b && typeof b === "object") { assert.deepEqual(Object.keys(a).sort(), Object.keys(b).sort(), `${path} keys`); return Object.keys(b).forEach((k) => near(a[k], b[k], `${path}.${k}`)); }
    assert.equal(a, b, path);
  };
  near(got, want);
  assert.equal(got.time.playing, undefined); // a property the view no longer has is gone
  // The receiver's base is not changed in place.
  assert.equal(base.view.mode, "3d");
  // Nothing changed: no state part at all.
  assert.equal(await encodeStatePart(JSON.parse(JSON.stringify(base)), base), "");
  // A link carries it as x= and keeps the readable camera and time.
  const viewer = { pose: cur.view.pose, timeline: { time: -12.5 }, state: { view: { mode: "sky" }, group: "All" } };
  const hv = decodeViewHash(`#${encodeViewHash(viewer, { state: part })}`);
  assert.equal(hv.state, part);
  assert.equal(hv.mode, "sky");
  // Unreadable parts are refused, not fatal.
  assert.equal(await decodeStatePart("1zAAAA", base), null);
  assert.equal(await decodeStatePart("nope", base), null);
  assert.equal(decodeViewHash("#t=0&c=0,0,0,1,0,0,60&x=1z***").state, undefined);
});

test("the longitude grid fades out within a Myr of the present day", async () => {
  const { presentDayFade } = await import("../../oviz/viewer/web/src/app/timeline.js");
  assert.equal(presentDayFade(0), 1);
  close(presentDayFade(-0.5), 0.5);
  close(presentDayFade(0.5), 0.5);
  assert.equal(presentDayFade(-1), 0);
  assert.equal(presentDayFade(-12), 0);
  assert.ok(presentDayFade(-0.2) > presentDayFade(-0.4));
});

test("fuzzy search prefers prefix and word matches", () => {
  const q = "orion";
  const s1 = fuzzyScore(q, "Orion Nebula Cluster");
  const s2 = fuzzyScore(q, "Proto Orion Family");
  const s3 = fuzzyScore(q, "Collinder 135");
  assert.ok(s1 > s2 && s2 > 0);
  assert.equal(s3, 0);
  assert.ok(fuzzyScore("ngc2516", "NGC 2516") > 0);
});

test("motion trails point back along playback and respect the time axis", async () => {
  const { trailParams } = await import("../../oviz/viewer/web/src/engine/frames.js");
  const times = [-60, -50, -40, -30, -20, -10, 0];
  const fwd = trailParams(times, 1, 6);
  close(fwd.myrPerFrame, 10);
  close(fwd.span, 0.6); // 6 Myr at 10 Myr per frame, toward earlier frames
  close(trailParams(times, -1, 6).span, -0.6);
  // A descending time axis flips the time step, not the frame direction.
  const desc = trailParams([...times].reverse(), 1, 6);
  close(desc.span, 0.6);
  close(desc.myrPerFrame, -10);
  assert.equal(trailParams(times, 1, 0), null);
  assert.equal(trailParams([0], 1, 5), null);
  assert.equal(trailParams([3, 3, 3], 1, 5), null);
});

test("spring easings settle, and livelier springs overshoot", async () => {
  const { springEasing, SPRINGS } = await import("../../oviz/viewer/web/src/ui/motion.js");
  const pop = springEasing(SPRINGS.pop);
  assert.ok(pop.duration > 150 && pop.duration < 2500, `duration ${pop.duration}`);
  assert.equal(springEasing(SPRINGS.pop), pop); // cached
  const stiff = springEasing({ stiffness: 900, damping: 60 });
  assert.ok(stiff.duration < springEasing({ stiffness: 120, damping: 14 }).duration);
});

test("the birth profile counts births inside the timeline only", async () => {
  const { Timeline, birthProfile } = await import("../../oviz/viewer/web/src/app/timeline.js");
  const tl = new Timeline([-100, -75, -50, -25, 0]);
  const profile = birthProfile(tl, [Float64Array.from([10, 12, 90, 200, -5])], 4);
  assert.equal(profile.length, 4);
  // Ages 200 (born before the timeline) and −5 (in the future) are ignored.
  const inside = birthProfile(tl, [Float64Array.from([10, 12, 90])], 4);
  assert.deepEqual(Array.from(profile), Array.from(inside));
  // Two births near the present outweigh the one at −90 Myr.
  assert.ok(profile[3] > profile[0] && profile[0] > 0);
  assert.equal(birthProfile(tl, [], 4).every((x) => x === 0), true);
});

test("the CPU birth fade matches the point shader", async () => {
  const { birthFadeAt } = await import("../../oviz/viewer/web/src/app/timeline.js");
  // A 20 Myr old cluster is born at t = −20.
  assert.equal(birthFadeAt(-30, 20), 0);
  assert.equal(birthFadeAt(-20, 20), 1);
  assert.equal(birthFadeAt(0, 20), 1);
  // With a fade it appears over the 8 Myr before birth.
  assert.equal(birthFadeAt(-28.5, 20, 8), 0);
  close(birthFadeAt(-24, 20, 8), 0.5);
  assert.equal(birthFadeAt(-20, 20, 8), 1);
  // Fade in and out: gone again once the fade after birth has passed.
  close(birthFadeAt(-16, 20, 8, true), 0.5);
  assert.equal(birthFadeAt(-5, 20, 8, true), 0);
  // No age (NaN or the shader's sentinel): always visible.
  assert.equal(birthFadeAt(-100, NaN, 8), 1);
  assert.equal(birthFadeAt(-100, -1e30, 8), 1);
});

test("looping playback wraps however far one tick goes", () => {
  const tl = new Timeline([-1, 0], { initialIndex: 0 });
  tl.loop = true;
  tl.speed = 8;
  tl.play(1);
  let t = 1000;
  tl.tick(t);
  const seen = new Set();
  for (let i = 0; i < 40; i++) { t += 100; tl.tick(t); seen.add(Math.round(tl.frame * 1000)); }
  assert.ok(tl.playing);
  // It keeps cycling instead of sticking on the last frame.
  assert.ok(seen.size >= 2, `stuck at ${[...seen]}`);
  for (const f of seen) assert.ok(f >= 0 && f <= 1000);
});

test("view links survive malformed escapes", () => {
  const v = decodeViewHash("#g=All&t=-1%&c=0,0,0,100,0,0,60");
  assert.equal(v.group, "All");
  assert.equal(v.time, undefined);
  assert.equal(v.pose.distance, 100);
  assert.deepEqual(decodeViewHash("#t=-12.5&v=sky"), { time: -12.5, mode: "sky" });
});

test("CPU positions vanish exactly where the GPU draws nothing", () => {
  // One object, three frames; absent at frame 0.
  const P = Float64Array.from([NaN, NaN, NaN, 1, 2, 3, 4, 5, 6]);
  assert.equal(cpuFramePosition(P, 3, 1, 0, 0), null);
  assert.deepEqual(Array.from(cpuFramePosition(P, 3, 1, 0, 0.5)), [1, 2, 3]);
  // Leaving: absent at the last frame.
  const Q = Float64Array.from([1, 2, 3, 4, 5, 6, NaN, NaN, NaN]);
  assert.deepEqual(Array.from(cpuFramePosition(Q, 3, 1, 0, 1)), [4, 5, 6]);
  assert.equal(cpuFramePosition(Q, 3, 1, 0, 2), null);
});

test("a straight-down eye keeps +x to the right and stays inside the orbit limit", async () => {
  const { poseFromEyeTarget, PITCH_LIMIT } = await import("../../oviz/viewer/web/src/engine/camera.js");
  const p = poseFromEyeTarget([0, 0, 2000], [0, 0, 0], 60);
  close(p.pitch, -PITCH_LIMIT);
  close(p.yaw, Math.PI / 2);
  const up = poseFromEyeTarget([0, 0, -10], [0, 0, 0], 60);
  close(up.pitch, PITCH_LIMIT);
});

test("USDZ packages are aligned, stored ZIPs with the scene first", async () => {
  const { usdzPackage, usdaScene } = await import("../../oviz/viewer/web/src/ar/usdz.js");
  const scene = {
    name: "Test \"figure\"",
    materials: [{ name: "Points_0", diffuse: [0.3, 0.2, 0.1], emissive: [1, 0.5, 0] },
      { name: "Dust_1", texture: "textures/cloud.png", diffuse: [1, 1, 1] }],
    meshes: [{ name: "Points_0", material: "Points_0", points: Float32Array.from([0, 0, 0, 1, 0, 0, 0, 1, 0]), normals: Float32Array.from([0, 0, 1, 0, 0, 1, 0, 0, 1]), indices: Uint32Array.from([0, 1, 2]) }],
  };
  const text = usdaScene(scene);
  assert.match(text, /^#usda 1\.0/);
  assert.match(text, /defaultPrim = "Root"/);
  assert.match(text, /upAxis = "Y"/);
  assert.match(text, /rel material:binding = <\/Root\/Materials\/Points_0>/);
  assert.match(text, /asset inputs:file = @textures\/cloud\.png@/);
  assert.ok(!text.includes('Test "figure"'), "quotes are stripped from the scene name");
  const png = new Uint8Array([137, 80, 78, 71, 1, 2, 3]);
  const zip = usdzPackage(scene, { "textures/cloud.png": png });
  const dv = new DataView(zip.buffer, zip.byteOffset, zip.byteLength);
  // Walk the local headers: stored, every data start on a 64-byte boundary.
  const names = [];
  let p = 0;
  while (dv.getUint32(p, true) === 0x04034b50) {
    assert.equal(dv.getUint16(p + 8, true), 0, "stored, not deflated");
    const size = dv.getUint32(p + 18, true);
    const nameLen = dv.getUint16(p + 26, true), extra = dv.getUint16(p + 28, true);
    names.push(new TextDecoder().decode(zip.subarray(p + 30, p + 30 + nameLen)));
    const data = p + 30 + nameLen + extra;
    assert.equal(data % 64, 0, `${names.at(-1)} data is 64-byte aligned`);
    p = data + size;
  }
  assert.deepEqual(names, ["scene.usda", "textures/cloud.png"]);
  assert.equal(dv.getUint32(p, true), 0x02014b50, "central directory follows");
  const eocd = zip.length - 22;
  assert.equal(dv.getUint32(eocd, true), 0x06054b50);
  assert.equal(dv.getUint16(eocd + 10, true), 2);
  assert.equal(dv.getUint32(eocd + 16, true), p, "central directory offset");
});

test("AR spheres are closed, outward-facing icospheres", async () => {
  const { icosphere } = await import("../../oviz/viewer/web/src/ar/model.js");
  for (const level of [0, 2]) {
    const { verts, faces } = icosphere(level);
    assert.equal(faces.length / 3, 20 * 4 ** level);
    for (let f = 0; f < faces.length; f += 3) {
      const [a, b, c] = [faces[f], faces[f + 1], faces[f + 2]].map((i) => [verts[i * 3], verts[i * 3 + 1], verts[i * 3 + 2]]);
      const u = [b[0] - a[0], b[1] - a[1], b[2] - a[2]], w = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
      const n = [u[1] * w[2] - u[2] * w[1], u[2] * w[0] - u[0] * w[2], u[0] * w[1] - u[1] * w[0]];
      const centre = [(a[0] + b[0] + c[0]) / 3, (a[1] + b[1] + c[1]) / 3, (a[2] + b[2] + c[2]) / 3];
      assert.ok(n[0] * centre[0] + n[1] * centre[1] + n[2] * centre[2] > 0, "counter-clockwise from outside");
    }
  }
});

test("AR scenes carry time samples, shared geometry and the timeline", async () => {
  const { usdaScene, usdaGeometry } = await import("../../oviz/viewer/web/src/ar/usdz.js");
  const text = usdaScene({
    name: "T", time: { start: 0, end: 24, perSecond: 12 },
    materials: [{ name: "Points_0", diffuse: [1, 0, 0] }],
    nodes: [
      { name: "O_trace-1_3", geometry: "sphere", material: "Points_0", translateSamples: [[0, [0, 0, 0]], [12, [1, 0, 0]]], scaleSamples: [[0, 0], [12, 0.01]] },
      { name: "Plate", translate: [0, 0.002, 0], children: [{ name: "Label", mesh: { points: Float32Array.from([0, 0, 0, 1, 0, 0, 0, 0, 1]), normals: Float32Array.from([0, 1, 0, 0, 1, 0, 0, 1, 0]), indices: Uint32Array.from([0, 1, 2]) } }] },
    ],
  });
  assert.match(text, /startTimeCode = 0\n    endTimeCode = 24\n    timeCodesPerSecond = 12/);
  assert.match(text, /def Xform "O_trace_1_3" \(\n\s+prepend references = @\.\/geometries\/sphere\.usda@<\/Geometry>/);
  assert.match(text, /double3 xformOp:translate\.timeSamples = \{\n\s+0: \(0, 0, 0\),\n\s+12: \(1, 0, 0\),\n\s+\}/);
  assert.match(text, /float3 xformOp:scale\.timeSamples = \{\n\s+0: \(0, 0, 0\),\n\s+12: \(0\.01, 0\.01, 0\.01\),/);
  assert.match(text, /uniform token\[\] xformOpOrder = \["xformOp:translate", "xformOp:scale"\]/);
  assert.match(text, /def Xform "Plate"\n\s+\{\n\s+double3 xformOp:translate = \(0, 0\.002, 0\)[\s\S]*def Mesh "Label"/);
  assert.match(usdaGeometry({ points: Float32Array.from([0, 0, 0]), normals: Float32Array.from([0, 1, 0]), indices: Uint32Array.from([]) }), /defaultPrim = "Geometry"[\s\S]*def Xform "Geometry"[\s\S]*def Mesh "Mesh"/);
});

test("AR time samples and labels follow the timeline", async () => {
  const { timeSamples, windowKeys, timeText } = await import("../../oviz/viewer/web/src/ar/model.js");
  const tl = new Timeline(Array.from({ length: 301 }, (_, i) => -300 + i), { initialIndex: 300 });
  const s = timeSamples(tl, 150);
  // Every third frame (at most 150 samples), plus the present.
  assert.equal(s.length, 101);
  assert.equal(s.at(-1).time, 0, "ends at the present");
  for (let i = 1; i < s.length; i++) assert.ok(s[i].time > s[i - 1].time, "in time order");
  assert.equal(timeText(0), "Today");
  assert.equal(timeText(-12.5), "−12.5 Myr");
  assert.equal(timeText(-120), "−120 Myr");
  assert.deepEqual(windowKeys(0, 10, 30), [[0, 1], [9.99, 1], [10, 0], [30, 0]]);
  assert.deepEqual(windowKeys(10, 30, 30), [[0, 0], [9.99, 0], [10, 1], [29.99, 1]]);
});

test("camera anchors normalise, encode and infer", () => {
  assert.deepEqual(normalizeAnchor("LSR"), { kind: "lsr" });
  assert.deepEqual(normalizeAnchor({ kind: "object", trace: "trace-3", index: 12 }), { kind: "object", trace: "trace-3", index: 12 });
  assert.equal(normalizeAnchor({ kind: "object", trace: "", index: 1 }), null);
  assert.equal(normalizeAnchor({ kind: "object", trace: "t", index: -1 }), null);
  assert.equal(normalizeAnchor({ kind: "galaxy" }), null);
  assert.equal(encodeAnchor({ kind: "object", trace: "trace-3", index: 12 }), "o:trace-3:12");
  assert.deepEqual(decodeAnchor("o:trace-3:12"), { kind: "object", trace: "trace-3", index: 12 });
  assert.deepEqual(decodeAnchor("o:a:b:7"), { kind: "object", trace: "a:b", index: 7 });
  assert.equal(decodeAnchor("o:trace-3:x"), null);
  assert.equal(decodeAnchor("o::3"), null);
  assert.equal(decodeAnchor("elsewhere"), null);
  assert.ok(sameAnchor("sun", { kind: "sun" }));
  assert.ok(!sameAnchor("sun", "lsr"));
  // States and links from before anchors: the LSR if the camera orbits it.
  const home = { target: [0, 0, 0], distance: 6000, yaw: 0, pitch: 1.2, fov: 60 };
  assert.deepEqual(inferAnchor(home, LSR_POINT, 2), { kind: "lsr" });
  assert.deepEqual(inferAnchor({ ...home, target: [120, 0, 0] }, LSR_POINT, 2), { kind: "free" });
  assert.deepEqual(inferAnchor({ ...home, distance: 0 }, LSR_POINT, 2, { kind: "sun" }), { kind: "sun" });
  assert.deepEqual(inferAnchor(home, null, 2), { kind: "free" }, "no LSR anchor when the origin is not the LSR");
});

test("view links carry the camera anchor", () => {
  const viewer = {
    pose: { target: [10, 20, 30], distance: 500, yaw: 0.5, pitch: 0.2, fov: 50 },
    timeline: { time: -3 },
    state: { view: { mode: "3d", anchor: { kind: "object", trace: "trace-2", index: 5 } }, group: "" },
  };
  const v = decodeViewHash(`#${encodeViewHash(viewer)}`);
  assert.deepEqual(v.anchor, { kind: "object", trace: "trace-2", index: 5 });
  viewer.state.view.anchor = { kind: "lsr" };
  assert.ok(encodeViewHash(viewer).includes("a=lsr"));
  assert.equal(decodeViewHash("#t=0&a=bogus").anchor, undefined);
  assert.equal(decodeViewHash("#t=0").anchor, undefined);
});

test("frame stepping stops at the ends even when playback loops (classic)", () => {
  const tl = new Timeline([0, 1, 2], { initialIndex: 2 });
  tl.loop = true;
  tl.step(1);
  assert.equal(tl.frame, 2);
  tl.step(-5);
  assert.equal(tl.frame, 0);
});

test("view links carry the saved views and can open presenting them", async () => {
  const base = {
    version: 2, group: "All", time: { frame: 60, speed: 1 },
    view: { mode: "3d", pose: { target: [0, 0, 0], distance: 6053.1, yaw: 1.5708, pitch: -1.4382, fov: 60 }, anchor: { kind: "lsr" } },
    traces: { a: { visible: true }, b: { visible: true } },
    volumes: { dust: { visible: true, opacity: 1, vmax: 0.07 } },
    global: { pointSize: 1, glow: 0.6 }, sky: { layers: [], members: "stars" }, ext: {},
  };
  const views = [0, 1, 2, 3, 4, 5].map((i) => {
    const state = JSON.parse(JSON.stringify(base));
    state.time.frame = 60 - 10 * i;
    state.view.pose = { ...state.view.pose, distance: 6053.1 / (1 + i), yaw: 1.5708 + 0.3 * i };
    if (i % 2) state.traces.b.visible = false;
    if (i === 3) state.view.mode = "sky";
    return { name: `View ${i + 1}`, caption: i === 2 ? "The Radcliffe Wave forms" : "", camera: i === 4 ? "keep" : "follow", transition: i === 1 ? { duration_ms: 2500, easing: "linear" } : null, thumb: "data:image/jpeg;base64,AAAA", state };
  });
  const part = await encodeViewsPart(views, base, { duration_ms: 1800, easing: "easeInOutCubic" });
  assert.match(part, /^1z[\w-]+$/);
  assert.ok(part.length < 700, `views part is ${part.length} characters`);
  const got = await decodeViewsPart(part, base);
  // States arrive exactly, up to the links' 7 significant digits.
  const near = (a, b, path = "") => {
    if (typeof a === "number" && typeof b === "number") return assert.ok(Math.abs(a - b) <= 1e-6 * Math.max(1, Math.abs(b)), `${path}: ${a} vs ${b}`);
    if (Array.isArray(b)) { assert.equal(a.length, b.length, `${path} length`); return b.forEach((x, i) => near(a[i], x, `${path}[${i}]`)); }
    if (b && typeof b === "object") { assert.deepEqual(Object.keys(a).sort(), Object.keys(b).sort(), `${path} keys`); return Object.keys(b).forEach((k) => near(a[k], b[k], `${path}.${k}`)); }
    assert.equal(a, b, path);
  };
  assert.equal(got.items.length, 6);
  assert.deepEqual(got.defaultTransition, { duration_ms: 1800, easing: "easeInOutCubic" });
  got.items.forEach((it, i) => {
    assert.equal(it.name, views[i].name);
    assert.equal(it.caption, views[i].caption);
    assert.equal(it.camera, views[i].camera);
    assert.deepEqual(it.transition, views[i].transition);
    near(it.state, views[i].state, `view ${i + 1}`);
    assert.equal(it.thumb, undefined); // thumbnails stay behind and are redrawn on visiting
  });
  // A presentation link: the views travel as w= and p= says where to start.
  const viewer = { pose: base.view.pose, timeline: { time: 0 }, state: { view: { mode: "3d" }, group: "All" } };
  const hv = decodeViewHash(`#${encodeViewHash(viewer, { views: part, present: 1 })}`);
  assert.equal(hv.views, part);
  assert.equal(hv.present, 1);
  assert.equal(decodeViewHash(`#${encodeViewHash(viewer, { views: part })}`).present, undefined);
  // No views: no w= part; unreadable parts are refused, not fatal.
  assert.equal(await encodeViewsPart([], base), "");
  assert.equal(await decodeViewsPart("1zNOT-DEFLATE", base), null);
  assert.equal(decodeViewHash("#t=0&c=0,0,0,1,0,0,60&w=bad!").views, undefined);
});

test("a State change crossfades one lasso selection into the next", () => {
  const dim = 0.16;
  // A fade starts from the weight drawn (0–255) toward the level (255
  // shown, 128 dimmed, 0 hidden); from what the level already draws,
  // nothing moves whatever the mix, and at rest only the level counts.
  assert.equal(lassoBlend(255, 255, [0, 0], dim), 1);
  close(lassoBlend(Math.round(255 * dim), 128, [0.3, 0.7], dim), dim, 2e-3);
  assert.equal(lassoBlend(0, 0, [0.5, 0.5], dim), 0);
  assert.equal(lassoBlend(0, 128, [1, 1], dim), dim);
  assert.equal(lassoBlend(255, 0, [1, 1], dim), 0);
  // Fades: out over the first three quarters, in over the last three,
  // eased, never both gone mid-way, and exact at both ends.
  assert.deepEqual(lassoFade(0), [0, 0]);
  assert.deepEqual(lassoFade(1), [1, 1]);
  const [out, inn] = lassoFade(0.5);
  close(out, 0.7407, 1e-3);
  close(inn, 0.2593, 1e-3);
  let prevOut = -1, prevIn = -1;
  for (let i = 0; i <= 20; i++) {
    const [o, n] = lassoFade(i / 20);
    assert.ok(o >= prevOut && n >= prevIn, "fades only move forward");
    assert.ok(o >= n - 1e-9, "the old selection leaves before the new one arrives");
    prevOut = o; prevIn = n;
  }
  // An object leaving the selection follows the fade-out, one joining the
  // fade-in; one dimmed in the next view stops at the dim level.
  const mid = lassoFade(0.5);
  close(lassoBlend(255, 0, mid, dim), 1 - mid[0]);
  close(lassoBlend(0, 255, mid, dim), mid[1]);
  close(lassoBlend(0, 128, lassoFade(1), dim), dim);
  assert.ok(lassoBlend(255, 0, mid, dim) + lassoBlend(0, 255, mid, dim) > 0.5, "no blink mid-way");
});

test("a lasso fade cut short carries on from what it drew", () => {
  const dim = 0.16;
  // Shown → hidden, hidden → shown and shown → dimmed, interrupted at 40 %
  // of the flight by a change to another selection.
  const cases = [[255, 0, 255], [0, 255, 0], [255, 128, 0], [0, 255, 128]];
  const cut = lassoFade(0.4);
  for (const [from, to, next] of cases) {
    const drawn = lassoBlend(from, to, cut, dim);
    const start = Math.round(255 * drawn); // what PointBatch.lassoDrawn hands the next fade
    // The next fade begins where the last one was (within a byte) …
    close(lassoBlend(start, next, lassoFade(0), dim), drawn, 1 / 255);
    // … moves only toward its own level …
    const target = next > 191 ? 1 : next > 63 ? dim : 0;
    let prev = lassoBlend(start, next, lassoFade(0), dim);
    for (let i = 1; i <= 20; i++) {
      const w = lassoBlend(start, next, lassoFade(i / 20), dim);
      assert.ok(Math.abs(target - w) <= Math.abs(target - prev) + 1e-9, `${from}→${to}→${next} heads for its level`);
      prev = w;
    }
    // … and lands on it exactly.
    assert.equal(lassoBlend(start, next, lassoFade(1), dim), target);
  }
});
