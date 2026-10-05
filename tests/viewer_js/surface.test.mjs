// Double-click re-centring geometry: boxes, volume density, planes, lines.
import test from "node:test";
import assert from "node:assert/strict";

import { rayBox, volumeMedianDepth, rayPlaneRect, screenSegment } from "../../oviz/viewer/web/src/engine/surface.js";

const close = (a, b, eps = 1e-6) => assert.ok(Math.abs(a - b) < eps, `${a} ≉ ${b}`);

test("rays enter and leave boxes, rotated ones too", () => {
  const hit = rayBox([-10, 0, 0], [1, 0, 0], [-1, -1, -1], [1, 1, 1]);
  close(hit.near, 9);
  close(hit.far, 11);
  assert.equal(rayBox([-10, 5, 0], [1, 0, 0], [-1, -1, -1], [1, 1, 1]), null);
  assert.equal(rayBox([10, 0, 0], [1, 0, 0], [-1, -1, -1], [1, 1, 1]), null); // behind the eye
  // A long thin box turned 90° about z now spans y, so a ray along x crosses it thinly.
  const thin = rayBox([-10, 0, 0], [1, 0, 0], [-5, -0.5, -1], [5, 0.5, 1], [0, 0, 0], Math.PI / 2);
  close(thin.far - thin.near, 1, 1e-6);
  // From inside, the entry is the eye itself.
  close(rayBox([0, 0, 0], [0, 0, 1], [-1, -1, -1], [1, 1, 1]).near, 0);
});

test("a volume pick lands in the middle of the dust along the ray", () => {
  // 10 voxels along x, dust only in the last three.
  const dims = [10, 1, 1];
  const vox = new Uint8Array(10);
  vox[7] = vox[8] = vox[9] = 255;
  const hit = rayBox([-5, 0.5, 0.5], [1, 0, 0], [0, 0, 0], [10, 1, 1]);
  const t = volumeMedianDepth(hit, vox, dims, [0, 0, 0], [10, 1, 1], 0, 1, 100);
  close(t, 5 + 8.5, 0.11); // the middle of voxels 7–9
  // An empty ray is not a hit.
  assert.equal(volumeMedianDepth(hit, new Uint8Array(10), dims, [0, 0, 0], [10, 1, 1]), null);
  // The display window counts: dust below the low end is invisible.
  vox.fill(100);
  assert.equal(volumeMedianDepth(hit, vox, dims, [0, 0, 0], [10, 1, 1], 0.5, 1), null);
});

test("image planes and screen segments", () => {
  close(rayPlaneRect([0, 0, 10], [0, 0, -1], [0, 0, 0], 4, 4), 10);
  assert.equal(rayPlaneRect([3, 0, 10], [0, 0, -1], [0, 0, 0], 4, 4), null);
  assert.equal(rayPlaneRect([0, 0, 10], [1, 0, 0], [0, 0, 0], 4, 4), null);
  const s = screenSegment(5, 3, [0, 0], [10, 0]);
  close(s.dist, 3);
  close(s.u, 0.5);
  close(screenSegment(-4, 3, [0, 0], [10, 0]).dist, 5);
});
