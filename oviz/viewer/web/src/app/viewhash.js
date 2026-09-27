// Shareable view links: #t=<time>&v=<mode>&c=<target,distance,yaw,pitch,fov>&g=<group>

/** Compact, human-readable hash: #t=-12.5&v=3d&c=x,y,z,d,yaw,pitch,fov */
export function encodeViewHash(viewer) {
  const p = viewer.pose;
  const r = (x, d = 2) => Number(x.toFixed(d));
  const parts = [
    `t=${r(viewer.timeline.time, 3)}`,
    `v=${viewer.state.view.mode}`,
    `c=${[...p.target.map((x) => r(x, 1)), r(p.distance, 1), r(p.yaw, 4), r(p.pitch, 4), r(p.fov, 2)].join(",")}`,
  ];
  if (viewer.state.group) parts.push(`g=${encodeURIComponent(viewer.state.group)}`);
  return parts.join("&");
}

export function decodeViewHash(hash) {
  const out = {};
  for (const part of String(hash || "").replace(/^#/, "").split("&")) {
    const [k, val] = part.split("=");
    if (!k || val == null) continue;
    out[k] = decodeURIComponent(val);
  }
  if (!out.c && !out.t) return null;
  const view = {};
  if (out.t != null && Number.isFinite(Number(out.t))) view.time = Number(out.t);
  if (out.v) view.mode = out.v === "sky" ? "sky" : "3d";
  if (out.g) view.group = out.g;
  if (out.c) {
    const n = out.c.split(",").map(Number);
    if (n.length === 7 && n.every(Number.isFinite)) {
      view.pose = { target: n.slice(0, 3), distance: n[3], yaw: n[4], pitch: n[5], fov: n[6] };
    }
  }
  return view;
}
