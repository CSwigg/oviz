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
