// Annotations: the document model, curve and shell geometry, grouping and morphs.
import test from "node:test";
import assert from "node:assert/strict";

import {
  normalizeAnnotations, notesToAnnotations, curveSamples, sphereRotation, sphereAxes, rotationQuat,
  insideSphere, polylineDistance, assignGroup, annotGroupName, lerpAnnotations, annotationsTimeDependent,
  resolveAnnotPoint, emptyAnnotations,
} from "../../oviz/viewer/web/src/ui/annotate-model.js";

const close = (a, b, eps = 1e-9) => assert.ok(Math.abs(a - b) <= eps, `${a} vs ${b}`);

test("documents from files, links and Python are cleaned", () => {
  const doc = normalizeAnnotations({
    groups: [{ id: "g1", name: "Bubbles", auto: false }, { id: "empty" }],
    items: [
      { id: "a", group: "g1", kind: "sphere", center: [1, 2, 3], radius: 5, style: "shell", opacity: 7 },
      { id: "a", group: "g1", kind: "text", at: { pos: [0, 0, 0] }, text: "Same id" },
      { kind: "curve", points: [[0, 0, 0]] },              // one point: dropped
      { kind: "curve", points: [[0, 0, 0], [1, 1, 1]], group: "nowhere", dash: "zigzag" },
      { kind: "nonsense", at: [0, 0, 0] },
      { kind: "text", at: { trace: "trace-2", index: 4, offset: [0, 0, 0] } },
    ],
  });
  assert.equal(doc.items.length, 4);
  assert.deepEqual(doc.items[0].radii, [5, 5, 5]);
  assert.equal(doc.items[0].opacity, 1);                 // clamped
  assert.notEqual(doc.items[1].id, "a");                 // ids stay unique
  assert.equal(doc.items[2].dash, "solid");
  assert.deepEqual(doc.items[3].at, { trace: "trace-2", index: 4 });
  // Every item has a group; empty groups go.
  const ids = new Set(doc.groups.map((g) => g.id));
  assert.ok(doc.items.every((i) => ids.has(i.group)));
  assert.ok(!ids.has("empty"));
  assert.deepEqual(normalizeAnnotations(null), emptyAnnotations());
});

test("notes from earlier figures become labels in a Notes entry", () => {
  const doc = notesToAnnotations([
    { id: "n1", text: "Orion", anchor: { pos: [100, -20, 5] }, color: "#ff0000", size: 18 },
    { id: "n2", text: "On a cluster", anchor: { trace: "t", index: 3 } },
    { text: "broken" },
  ]);
  assert.equal(doc.groups.length, 1);
  assert.equal(doc.groups[0].name, "Notes");
  assert.deepEqual(doc.items.map((i) => [i.id, i.text, i.size]), [["n1", "Orion", 18], ["n2", "On a cluster", 13]]);
  assert.equal(doc.items[1].dy, -14);
});

test("curves pass through their points; smoothness 0 is straight", () => {
  const pts = [[0, 0, 0], [10, 5, 0], [20, 0, 3]];
  const smooth = curveSamples(pts, { smooth: 1, perSegment: 16 });
  assert.deepEqual(smooth.points[0], pts[0]);
  assert.deepEqual(smooth.points[16], pts[1]);
  assert.deepEqual(smooth.points.at(-1), pts[2]);
  // Centripetal Catmull-Rom stays between the points (no overshoot here).
  assert.ok(smooth.points.every((p) => p[1] <= 5 + 1e-9 && p[1] >= -1e-9));
  const straight = curveSamples(pts, { smooth: 0 });
  assert.deepEqual(straight.points, pts);
  close(straight.arc.at(-1), Math.hypot(10, 5) + Math.hypot(10, 5, 3));
  const loop = curveSamples(pts, { smooth: 1, closed: true, perSegment: 8 });
  assert.deepEqual(loop.points.at(-1), pts[0]);
  assert.equal(loop.points.length, 3 * 8 + 1);
});

test("ellipsoid rotations, axes and quaternions agree", () => {
  const item = { radii: [3, 2, 1], rot: [0, 0, 90] };
  const axes = sphereAxes(item);
  // The x semi-axis turns onto +y.
  close(axes[0][0], 0, 1e-12); close(axes[0][1], 3, 1e-12);
  close(axes[1][0], -2, 1e-12);
  const m = sphereRotation([20, -35, 110]);
  const q = rotationQuat(m);
  close(Math.hypot(...q), 1, 1e-12);
  // Rotate e_x by q and compare with the matrix's first column.
  const v = [1, 0, 0];
  const t = [2 * (q[1] * v[2] - q[2] * v[1]), 2 * (q[2] * v[0] - q[0] * v[2]), 2 * (q[0] * v[1] - q[1] * v[0])];
  const r = [v[0] + q[3] * t[0] + (q[1] * t[2] - q[2] * t[1]), v[1] + q[3] * t[1] + (q[2] * t[0] - q[0] * t[2]), v[2] + q[3] * t[2] + (q[0] * t[1] - q[1] * t[0])];
  close(r[0], m[0], 1e-9); close(r[1], m[3], 1e-9); close(r[2], m[6], 1e-9);
  assert.ok(insideSphere(item, [0, 0, 0], [0, 2.9, 0]));
  assert.ok(!insideSphere(item, [0, 0, 0], [2.9, 0, 0]));
});

test("screen distance to a polyline", () => {
  close(polylineDistance(5, 3, [[0, 0], [10, 0]]), 3);
  close(polylineDistance(-4, 3, [[0, 0], [10, 0]]), 5);
  assert.equal(polylineDistance(0, 0, [[1, 1]]), Infinity);
});

test("what is drawn groups itself the way people sort it", () => {
  const doc = emptyAnnotations();
  const add = (item) => { item.id = `i${doc.items.length}`; item.group = assignGroup(doc, item); doc.items.push(item); return item; };
  const shell = add({ kind: "sphere", center: { pos: [0, 0, 0] }, radii: [100, 100, 100], rot: [0, 0, 0], style: "shell", color: "#ff9e7a", opacity: 1 });
  assert.equal(annotGroupName(doc, doc.groups[0]), "Shell");
  // A label inside the shell joins it and names it.
  const label = add({ kind: "text", at: { pos: [20, 10, 0] }, text: "Local Bubble\nsecond line", color: "#fff", opacity: 1 });
  assert.equal(label.group, shell.group);
  assert.equal(annotGroupName(doc, doc.groups[0]), "Local Bubble");
  // Arrows of one colour make one entry; another colour starts another.
  const a1 = add({ kind: "curve", points: [{ pos: [500, 0, 0] }, { pos: [600, 0, 0] }], arrow: "end", color: "#ffd27a", opacity: 1 });
  const a2 = add({ kind: "curve", points: [{ pos: [500, 50, 0] }, { pos: [600, 50, 0] }], arrow: "end", color: "#FFD27A", opacity: 1 });
  const a3 = add({ kind: "curve", points: [{ pos: [500, 90, 0] }, { pos: [600, 90, 0] }], arrow: "end", color: "#4d8dff", opacity: 1 });
  assert.equal(a1.group, a2.group);
  assert.notEqual(a3.group, a1.group);
  assert.equal(annotGroupName(doc, doc.groups.find((g) => g.id === a1.group)), "Arrows");
  // A label beside a curve on screen joins it (given a screen measure).
  const near = assignGroup(doc, { kind: "text", at: { pos: [550, 90, 0] }, text: "Flow" }, { screenDistance: (c, p) => (c === a3 ? 4 : 100) });
  assert.equal(near, a3.group);
  // A bubble around an existing label joins the label's entry.
  const tag = add({ kind: "text", at: { pos: [-900, 0, 0] }, text: "Orion", color: "#fff", opacity: 1 });
  const bub = add({ kind: "sphere", center: { pos: [-880, 0, 0] }, radii: [60, 60, 60], rot: [0, 0, 0], style: "bubble", color: "#7cc4ff", opacity: 1 });
  assert.equal(bub.group, tag.group);
  assert.equal(annotGroupName(doc, doc.groups.find((g) => g.id === tag.group)), "Orion");
});

test("view transitions morph matching annotations and fade the rest", () => {
  const a = normalizeAnnotations({ groups: [{ id: "g" }], items: [
    { id: "s", group: "g", kind: "sphere", center: [0, 0, 0], radius: 10, color: "#000000" },
    { id: "gone", group: "g", kind: "text", at: [5, 5, 5], text: "Bye" },
  ] });
  const b = normalizeAnnotations({ groups: [{ id: "g" }], items: [
    { id: "s", group: "g", kind: "sphere", center: [100, 0, 0], radius: 40, color: "#ffffff" },
    { id: "new", group: "g", kind: "curve", points: [[0, 0, 0], [1, 0, 0]] },
  ] });
  const mid = lerpAnnotations(a, b, 0.5);
  const s = mid.items.find((i) => i.id === "s");
  assert.deepEqual(s.center.pos, [50, 0, 0]);
  close(s.radii[0], 20, 1e-9); // geometric mean: radii grow evenly
  assert.equal(s.color, "#808080");
  close(mid.items.find((i) => i.id === "gone").fade, 0.5);
  close(mid.items.find((i) => i.id === "new").fade, 0.5);
  // A hidden group fades its items.
  const hidden = normalizeAnnotations({ ...b, groups: [{ id: "g", visible: false }] });
  close(lerpAnnotations(b, hidden, 0.25).items.find((i) => i.id === "s").fade, 0.75);
});

test("only attached or present-day annotations depend on time", () => {
  const fixed = normalizeAnnotations({ items: [{ kind: "text", at: [0, 0, 0] }] });
  assert.equal(annotationsTimeDependent(fixed), false);
  const riding = normalizeAnnotations({ items: [{ kind: "curve", points: [[0, 0, 0], { trace: "t", index: 1 }] }] });
  assert.equal(annotationsTimeDependent(riding), true);
  const present = normalizeAnnotations({ items: [{ kind: "text", at: [0, 0, 0], present: true }] });
  assert.equal(annotationsTimeDependent(present), true);
  const pos = resolveAnnotPoint({ trace: "t", index: 1, offset: [1, 2, 3] }, () => [10, 10, 10]);
  assert.deepEqual(pos, [11, 12, 13]);
  assert.equal(resolveAnnotPoint({ trace: "t", index: 1 }, () => null), null);
});

test("boxes: half sizes, corners and what lies inside", async () => {
  const { boxCorners, insideBox, BOX_EDGES } = await import("../../oviz/viewer/web/src/ui/annotate-model.js");
  const doc = normalizeAnnotations({ items: [{ kind: "box", center: [10, 0, 0], size: [20, 10, 4], select: "isolate", style: "edges" }] });
  const box = doc.items[0];
  assert.deepEqual(box.radii, [10, 5, 2]);
  assert.equal(box.style, "box");           // unknown styles fall back
  assert.equal(box.select, "isolate");
  const c = boxCorners(box, [10, 0, 0]);
  assert.deepEqual(c[0], [0, -5, -2]);
  assert.deepEqual(c[7], [20, 5, 2]);
  assert.equal(BOX_EDGES.length, 12);
  assert.ok(BOX_EDGES.every(([a, b]) => [1, 2, 4].includes(a ^ b)));
  assert.ok(insideBox(box, [10, 0, 0], [19, 4, 1]));
  assert.ok(!insideBox(box, [10, 0, 0], [21, 0, 0]));
  // Turned 90° about z, the long side lies along y.
  const turned = { ...box, rot: [0, 0, 90] };
  assert.ok(insideBox(turned, [10, 0, 0], [10, 9, 0]));
  assert.ok(!insideBox(turned, [10, 0, 0], [19, 0, 0]));
});

test("selections: inside a shell, a box, or within reach of a curve", async () => {
  const { selectionTest, curveReach } = await import("../../oviz/viewer/web/src/ui/annotate-model.js");
  const shell = normalizeAnnotations({ items: [{ kind: "sphere", center: [0, 0, 0], radius: 10, select: "highlight" }] }).items[0];
  const inShell = selectionTest(shell, { center: [0, 0, 0] });
  assert.ok(inShell([0, 9, 0]) && !inShell([0, 11, 0]));
  assert.equal(selectionTest({ ...shell, select: undefined }, { center: [0, 0, 0] }), null);
  const line = normalizeAnnotations({ items: [{ kind: "curve", points: [[0, 0, 0], [100, 0, 0]], select: "highlight", reach: 5 }] }).items[0];
  const samples = curveSamples([[0, 0, 0], [100, 0, 0]], { smooth: 0 }).points;
  const near = selectionTest(line, { samples });
  assert.ok(near([50, 4, 0]) && near([50, 0, -4.9]));
  assert.ok(!near([50, 6, 0]) && !near([-6, 0, 0]));
  // Without a reach a curve selects within a twentieth of its length.
  assert.equal(curveReach({ ...line, reach: undefined }, samples), 5);
});

test("rotations turn and read back as the same Euler angles", async () => {
  const { rotationToEuler, axisRotation, turnEuler, mul3 } = await import("../../oviz/viewer/web/src/ui/annotate-model.js");
  for (const rot of [[0, 0, 0], [30, -20, 75], [-120, 40, 10], [5, 89, -30]]) {
    const back = rotationToEuler(sphereRotation(rot));
    const m1 = sphereRotation(rot), m2 = sphereRotation(back);
    m1.forEach((x, i) => close(x, m2[i], 1e-6));
  }
  // A quarter turn about z, applied twice, is a half turn.
  const q = axisRotation([0, 0, 1], Math.PI / 2);
  const half = rotationToEuler(mul3(q, q));
  sphereRotation(half).forEach((x, i) => close(x, sphereRotation([0, 0, 180])[i], 1e-9));
  const turned = turnEuler([0, 0, 0], axisRotation([1, 0, 0], Math.PI / 6));
  close(turned[0], 30, 1e-6); close(turned[1], 0, 1e-6); close(turned[2], 0, 1e-6);
});
