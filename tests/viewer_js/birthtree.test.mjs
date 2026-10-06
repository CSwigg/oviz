// Birth tree (classic dendrogram): links, thresholds and plot order.
import test from "node:test";
import assert from "node:assert/strict";

import { buildBirthTree } from "../../oviz/viewer/web/src/ui/birthtree.js";

// Clusters on straight tracks along x: pos(t) = x0 + vx * (t − birthTime).
const cluster = (key, age, x0, vx = 0) => ({
  key, name: key, age, birthTime: -age, birth: [x0, 0, 0],
  at: (t) => [x0 + vx * (t + age), 0, 0],
});

test("each cluster links to the nearest older one within the distance threshold", () => {
  const A = cluster("A", 30, 0, 2);   // born at x = 0, drifting +2 pc/Myr
  const B = cluster("B", 20, 25);     // born 25 pc away; A is at x = 20 then
  const C = cluster("C", 10, 300);    // too far from everyone
  const D = cluster("D", 5, 22);      // nearest older at its birth: B (22 vs 25)
  const t = buildBirthTree([C, B, D, A], { threshold: 50 });
  const parent = Object.fromEntries(t.nodes.map((n) => [n.key, n.parent]));
  assert.deepEqual(parent, { A: null, B: "A", C: null, D: "B" });
  assert.equal(t.branches.length, 2);
  const bAB = t.branches.find((b) => b.key === "B->A");
  assert.ok(Math.abs(bAB.distance - 5) < 1e-9); // A's track is at x = 20 when B is born
  assert.deepEqual(bAB.keys, ["B", "D"]);        // B and its descendant D
  assert.equal(t.roots.map((n) => n.key).join(), "A,C");
  // Leaves get consecutive orders, parents sit over their children.
  const byKey = Object.fromEntries(t.nodes.map((n) => [n.key, n]));
  assert.equal(byKey.D.order, 0);
  assert.equal(byKey.B.order, byKey.D.order);
  assert.equal(t.leafCount, 2);
  assert.equal(byKey.A.axis, 30);                // distance mode: height = birth age
});

test("birth-to-birth links and the birth-age threshold mode", () => {
  const A = cluster("A", 30, 0, 2);
  const B = cluster("B", 20, 25);
  const tree = buildBirthTree([A, B], { threshold: 10, links: "birth_to_birth" });
  assert.equal(tree.byKey.get("B").parent, null); // 25 pc from A's birth site: over 10 pc
  const byAge = buildBirthTree([A, B, cluster("E", 28, 400)], { threshold: 5, mode: "birth_age_myr" });
  // B is 10 Myr younger than A and 8 Myr younger than E: over 5 Myr, so a root.
  assert.equal(byAge.byKey.get("B").parent, null);
  // E is 2 Myr younger than A, within 5 Myr: linked to A.
  assert.equal(byAge.byKey.get("E").parent, "A");
  assert.ok(byAge.byKey.get("E").axis > 0);       // cumulative birth distance
  assert.equal(byAge.mode, "birth_age_myr");
});

test("relative SFH: a unit-area Gaussian KDE scaled to its own peak", async () => {
  const { gaussianKde, relativeCurve } = await import("../../oviz/viewer/web/src/ui/sfh.js");
  const grid = Array.from({ length: 601 }, (_, i) => -30 + i * 0.05);
  const d = gaussianKde([-10, -10, -20], grid, 1);
  let area = 0;
  for (let i = 1; i < grid.length; i++) area += 0.5 * (d[i] + d[i - 1]) * 0.05;
  assert.ok(Math.abs(area - 1) < 1e-3);
  const r = relativeCurve(d);
  assert.equal(Math.max(...r), 1);
  const at = (t) => r[Math.round((t + 30) / 0.05)];
  assert.ok(Math.abs(at(-10) - 1) < 1e-6);
  assert.ok(Math.abs(at(-20) - 0.5) < 1e-3); // one cluster there against two at −10
  assert.deepEqual(Array.from(gaussianKde([], [0, 1], 2)), [0, 0]);
});
