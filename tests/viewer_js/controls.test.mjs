// Keyboard camera motion must keep the classic viewer's directions.
import test from "node:test";
import assert from "node:assert/strict";

import { keyMotion } from "../../oviz/viewer/web/src/engine/controls.js";
import { makePose, poseEye, forwardFromAngles } from "../../oviz/viewer/web/src/engine/camera.js";

const hold = (pose, key, opts = {}) => keyMotion(pose, new Set([key]), { dt: 0.25, maxSpan: 4000, ...opts });
const pose3d = () => makePose({ target: [0, 0, 0], distance: 1000, yaw: -Math.PI / 2, pitch: 0.4, fov: 50 });
const azimuth = (p) => { const e = poseEye(p); return Math.atan2(e[1] - p.target[1], e[0] - p.target[0]); };

test("W tilts toward overhead and S away, like the classic viewer", () => {
  const up = hold(pose3d(), "w");
  const down = hold(pose3d(), "s");
  assert.ok(poseEye(up)[2] > poseEye(pose3d())[2]);
  assert.ok(poseEye(down)[2] < poseEye(pose3d())[2]);
});

test("A orbits the camera azimuth up and D down", () => {
  const start = azimuth(pose3d());
  assert.ok(azimuth(hold(pose3d(), "a")) > start);
  assert.ok(azimuth(hold(pose3d(), "d")) < start);
});

test("Q zooms out, E zooms in, R/F move along Galactic z", () => {
  assert.ok(hold(pose3d(), "q").distance > 1000);
  assert.ok(hold(pose3d(), "e").distance < 1000);
  assert.ok(hold(pose3d(), "r").target[2] > 0);
  assert.ok(hold(pose3d(), "f").target[2] < 0);
});

test("Shift+W flies forward and Shift+D to the right without orbiting", () => {
  const p = pose3d();
  const f = forwardFromAngles(p.yaw, p.pitch);
  const w = hold(pose3d(), "w", { fast: true });
  assert.ok(w.target.reduce((s, x, i) => s + x * f[i], 0) > 0);
  assert.equal(w.yaw, p.yaw);
  assert.equal(w.pitch, p.pitch);
  const d = hold(pose3d(), "d", { fast: true });
  const right = [f[1], -f[0], 0];
  assert.ok(d.target.reduce((s, x, i) => s + x * right[i], 0) > 0);
});

test("Sky view: keys turn the gaze and Q/E change the field of view", () => {
  const sky = () => makePose({ target: [0, 0, 0], distance: 0, yaw: 0, pitch: 0, fov: 60 });
  assert.ok(hold(sky(), "a", { sky: true }).yaw > 0);
  assert.ok(hold(sky(), "q", { sky: true }).fov > 60);
  assert.ok(hold(sky(), "e", { sky: true }).fov < 60);
  const r = hold(sky(), "r", { sky: true });
  assert.deepEqual(r.target, [0, 0, 0]);
  assert.equal(r.distance, 0);
});
