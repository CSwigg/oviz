// Unit tests for the Oviz viewer's pure modules (run with `node --test`).
import test from "node:test";
import assert from "node:assert/strict";

import { galToIcrs, icrsToGal, m4invert, m4mul, m4perspective, m4lookAt, niceFloor, formatDistance } from "../../oviz/viewer/web/src/core/math.js";
import { Timeline } from "../../oviz/viewer/web/src/app/timeline.js";
import { cpuFramePosition, frameOffset } from "../../oviz/viewer/web/src/engine/frames.js";
import { unshuffle } from "../../oviz/viewer/web/src/core/loader.js";
import { poseFromEyeTarget, poseEye, lerpPose, makePose } from "../../oviz/viewer/web/src/engine/camera.js";
import { encodeViewHash, decodeViewHash } from "../../oviz/viewer/web/src/app/viewhash.js";
import { fuzzyScore } from "../../oviz/viewer/web/src/ui/dom.js";

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
