// Picking: a click lands on the point itself, not on its glow.
import test from "node:test";
import assert from "node:assert/strict";

import { POINT_PICK_CORE, PICK_MIN_PX, starProfile } from "../../oviz/viewer/web/src/layers/points.js";
import { MEMBER_PICK_CORE } from "../../oviz/viewer/web/src/sky/members.js";

/** Outermost radius (sprite units) where the profile is above `level`. */
function edge(profile, level) {
  let r = 0;
  for (let x = 0; x < 1; x += 0.0005) if (profile(x) > level) r = x;
  return r;
}

test("a cluster point is picked on its bright core, not its glow", () => {
  const core = (r) => starProfile("core", r);
  // The pick disc reaches the core's edge, where it has faded to a few percent.
  assert.ok(core(POINT_PICK_CORE) < 0.05 * core(0));
  assert.ok(core(POINT_PICK_CORE / 2) > 0.15 * core(0));
  for (const glow of [0.2, 0.6, 1.0]) {
    // Sprite geometry from the point shader: the quad is the glow's size, the core is smaller.
    const coreRatio = (3.15 + 0.7 * glow) / (2.65 + 0.18 * glow);
    const pickRadius = POINT_PICK_CORE / coreRatio; // in quad radii
    const haloAlpha = (r) => starProfile("halo", r) * Math.min(0.34 + 0.18 * glow, 0.78);
    // The glow stays visible (above 2% alpha) out to over three pick radii:
    // a click out there misses.
    assert.ok(edge(haloAlpha, 0.02) > 3 * pickRadius, `glow ${glow}`);
  }
  assert.ok(PICK_MIN_PX <= 4, "small points stay precise; pick()'s search radius forgives near misses");
});

test("a Sky member star is picked on its core, not its glow", () => {
  const core = (r) => starProfile("member-core", r);
  assert.ok(core(MEMBER_PICK_CORE) < 0.1 * core(0));
  const haloAlpha = (r) => starProfile("member-halo", r) * 0.57;
  assert.ok(edge(haloAlpha, 0.02) > 1.5 * MEMBER_PICK_CORE);
  assert.ok(MEMBER_PICK_CORE < 0.25, "the old pick disc was the whole sprite (1.0)");
});
