// Padding (a host page's sidebar over the figure) must centre the view in
// the uncovered part without moving the pose, and keep Sky registered.
import test from "node:test";
import assert from "node:assert/strict";

import { Camera, makePose, forwardFromAngles } from "../../oviz/viewer/web/src/engine/camera.js";
import { Controls } from "../../oviz/viewer/web/src/engine/controls.js";
import { readPadding } from "../../oviz/viewer/web/src/app/viewer.js";

const DEG = Math.PI / 180;

function camera(padding, pose = {}) {
  const cam = new Camera();
  cam.setViewport(1440, 900);
  cam.pose = makePose({ target: [120, -40, 15], distance: 900, yaw: -1.1, pitch: 0.5, fov: 45, ...pose });
  cam.setPadding(padding);
  cam.update(1e4);
  return cam;
}

test("the orbit target lands in the middle of the uncovered view", () => {
  const cam = camera({ right: 640 });
  assert.deepEqual(cam.center, [400, 450]);
  const s = cam.project(cam.pose.target);
  assert.ok(Math.abs(s[0] - 400) < 1e-3 && Math.abs(s[1] - 450) < 1e-3);
  const both = camera({ top: 100, right: 300, bottom: 200, left: 100 });
  const t = both.project(both.pose.target);
  assert.ok(Math.abs(t[0] - 620) < 1e-3 && Math.abs(t[1] - 400) < 1e-3);
});

test("padding moves the image, not the pose or its scale", () => {
  const plain = camera({});
  const padded = camera({ right: 640 });
  assert.deepEqual(padded.pose, plain.pose);
  const p = [300, 200, -50];
  const a = plain.project(p), b = padded.project(p);
  assert.ok(Math.abs(b[0] - (a[0] - 320)) < 1e-3 && Math.abs(b[1] - a[1]) < 1e-3);
  // Rays through the shifted pixel run where they did before.
  const ra = plain.ray(a[0], a[1]), rb = padded.ray(b[0], b[1]);
  for (let i = 0; i < 3; i++) assert.ok(Math.abs(ra.dir[i] - rb.dir[i]) < 1e-5);
});

test("a fifth of each axis always stays uncovered", () => {
  const cam = camera({ left: 900, right: 900 });
  const used = cam.usedPadding();
  assert.ok(Math.abs(used.left + used.right - 1440 * 0.8) < 1e-3);
  assert.deepEqual(cam.center, [720, 450]);
});

test("the Sky frame keeps a TAN survey registered with the padded view", () => {
  for (const padding of [{}, { right: 640 }, { bottom: 380 }, { left: 120, top: 60 }]) {
    const cam = camera(padding, { target: [0, 0, 0], distance: 0, yaw: 0.7, pitch: 0.2, fov: 50 });
    const f = cam.skyFrame();
    // The frame is centred on the optical centre and covers the window.
    assert.ok(Math.abs(f.left + f.width / 2 - cam.center[0]) < 1e-9 && Math.abs(f.top + f.height / 2 - cam.center[1]) < 1e-9);
    assert.ok(f.left <= 1e-9 && f.top <= 1e-9 && f.left + f.width >= 1440 - 1e-9 && f.top + f.height >= 900 - 1e-9);
    // A TAN view of the frame's width and field puts a direction where the camera does.
    const scale = f.width / 2 / Math.tan((f.hfov * DEG) / 2);
    for (const [dy, dp] of [[0.2, 0.1], [-0.3, -0.15], [0.05, 0.25]]) {
      const d = forwardFromAngles(cam.pose.yaw + dy, cam.pose.pitch + dp);
      const s = cam.project([d[0] * 1e5, d[1] * 1e5, d[2] * 1e5]);
      const fw = cam.forward, rt = cam.right, up = cam.camUp;
      const z = d[0] * fw[0] + d[1] * fw[1] + d[2] * fw[2];
      const x = f.left + f.width / 2 + scale * (d[0] * rt[0] + d[1] * rt[1] + d[2] * rt[2]) / z;
      const y = f.top + f.height / 2 - scale * (d[0] * up[0] + d[1] * up[1] + d[2] * up[2]) / z;
      assert.ok(Math.abs(s[0] - x) < 1e-3 && Math.abs(s[1] - y) < 1e-3, `${JSON.stringify(padding)}: ${s} vs ${x},${y}`);
    }
  }
  const plain = camera({}, { distance: 0 });
  assert.ok(Math.abs(plain.skyFrame().hfov - plain.hfov) < 1e-9);
});

test("Sky zoom keeps the direction under the cursor, padded or not", () => {
  for (const padding of [{}, { right: 640 }]) {
    const cam = camera(padding, { target: [0, 0, 0], distance: 0, yaw: 0.4, pitch: 0.1, fov: 40 });
    const controls = Object.create(Controls.prototype);
    Object.assign(controls, { camera: cam, mode: "sky", minFov: 1, maxFov: 120, _changed: () => {} });
    Object.defineProperty(controls, "pose", { get: () => cam.pose });
    // Zooming about the optical centre leaves the view direction alone.
    const [cx, cy] = cam.center;
    const before = [cam.pose.yaw, cam.pose.pitch];
    controls.zoomAt(0.8, cx, cy);
    assert.ok(Math.abs(cam.pose.yaw - before[0]) < 1e-12 && Math.abs(cam.pose.pitch - before[1]) < 1e-12);
  }
});

test("padding reads from objects and from ?padding=top,right,bottom,left", () => {
  assert.deepEqual(readPadding("0,640,0,0"), { top: 0, right: 640, bottom: 0, left: 0 });
  assert.deepEqual(readPadding("10, 20"), { top: 10, right: 20 });
  assert.deepEqual(readPadding({ right: "300", left: -5, top: "x" }), { right: 300, left: 0 });
  assert.deepEqual(readPadding(""), {});
  assert.deepEqual(readPadding(null), {});
});
