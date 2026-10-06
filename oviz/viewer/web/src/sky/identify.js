// Sky identify and name lookup (classic Oviz): right-click (or long-press)
// the sky to list the catalogued objects there from SIMBAD, clusters and
// nebulae first; and find a named object on the sky from the search box
// through the CDS Sesame resolver.

import { h, icon, iconButton, clear } from "../ui/dom.js";
import { dirToLb, galToIcrs, icrsToGal, clamp, DEG, angleDelta } from "../core/math.js";
import { clonePose } from "../engine/camera.js";

const TAP_URL = "https://simbad.cds.unistra.fr/simbad/sim-tap/sync";
const SESAME_URL = "https://cds.unistra.fr/cgi-bin/nph-sesame/-oI/SNV?";

// Rank by scientific prominence, not raw proximity: clusters and
// associations, then nebulae and ISM structures, then galaxies, then stars;
// within a tier by how often SIMBAD's literature cites them.
const TIERS = [
  { tier: 5, types: ["Cl*", "OpC", "GlC", "As*", "St*", "MGr", "Cl?", "OB?"] },
  { tier: 4, types: ["HII", "SNR", "PN", "RNe", "GNe", "DNe", "MoC", "Cld", "glb", "sh", "SFR", "HH", "bub", "flt", "ISM", "Neb"] },
  { tier: 3, types: ["G", "AGN", "SyG", "Sy1", "Sy2", "LIN", "SBG", "bCG", "GiC", "GiG", "GiP", "IG", "PaG", "GrG", "ClG", "QSO", "Bla", "BLL", "rG"] },
];

export function simbadTier(otype) {
  const code = String(otype ?? "").trim();
  for (const g of TIERS) if (g.types.includes(code)) return g.tier;
  return 1;
}

/** SIMBAD cone rows [main_id, otype, ra, dec, oid, sp_type, refs, sep] → ranked entries. */
export function rankSimbadRows(rows, limit = 8) {
  return (Array.isArray(rows) ? rows : [])
    .map((r) => ({
      name: String(r[0] ?? "").replace(/\s+/g, " ").trim(),
      otype: String(r[1] ?? "").trim(),
      ra: Number(r[2]), dec: Number(r[3]), oid: Number(r[4]),
      spType: String(r[5] ?? "").trim(),
      refs: Number(r[6]) || 0,
      sep: Number(r[7]),
      tier: simbadTier(r[1]),
    }))
    .filter((e) => e.name)
    .sort((a, b) => b.tier - a.tier || b.refs - a.refs || a.sep - b.sep)
    .slice(0, limit);
}

/** Parse Sesame's plain-text answer: `{ ra, dec, type, name }`, or null. */
export function parseSesame(text) {
  const s = String(text || "");
  const m = /^%J\s+([+-]?\d+(?:\.\d+)?)\s+([+-]?\d+(?:\.\d+)?)/m.exec(s);
  if (!m) return null;
  const type = /^%C\.0\s+(\S+)/m.exec(s)?.[1] || "";
  const ident = /^%I\.0\s+(.+)$/m.exec(s)?.[1]?.replace(/\s+/g, " ").trim() || "";
  const name = /^#\s*(\S.*?)\s*\t/m.exec(s)?.[1] || ident;
  return { ra: Number(m[1]), dec: Number(m[2]), type, name };
}

/** "12.3′", "1.4°", "25″" (classic separation readout). */
export function formatSeparation(deg) {
  const arcmin = Number(deg) * 60;
  if (!Number.isFinite(arcmin)) return "";
  if (arcmin >= 60) return `${(arcmin / 60).toFixed(1)}°`;
  if (arcmin >= 1) return `${arcmin.toFixed(1)}′`;
  return `${(arcmin * 60).toFixed(0)}″`;
}

async function tap(adql, signal) {
  const url = `${TAP_URL}?request=doQuery&lang=adql&format=json&query=${encodeURIComponent(adql)}`;
  const res = await fetch(url, { signal });
  if (!res.ok) throw new Error(`SIMBAD answered ${res.status}`);
  return res.json();
}

export async function resolveName(name, signal) {
  const res = await fetch(SESAME_URL + encodeURIComponent(String(name).trim()), { signal });
  if (!res.ok) throw new Error(`Sesame answered ${res.status}`);
  return parseSesame(await res.text());
}

export class SkyIdentify {
  constructor(plugin) {
    this.plugin = plugin;
    this.ui = plugin.ui;
    this.viewer = plugin.viewer;
    this.el = null;
    this.abort = null;
    this.bind();
  }

  bind() {
    const v = this.viewer;
    const canvas = v.canvas;
    let down = null;
    let press = 0;
    canvas.addEventListener("contextmenu", (e) => { if (v.state.view.mode === "sky") e.preventDefault(); });
    canvas.addEventListener("pointerdown", (e) => {
      if (v.state.view.mode !== "sky") return;
      down = { x: e.clientX, y: e.clientY, button: e.button, type: e.pointerType };
      clearTimeout(press);
      // A long press on a touch screen is the right-click.
      if (e.pointerType === "touch") press = setTimeout(() => { if (down && v.controls.pointers.size <= 1) { const d = down; down = null; this.open(d.x, d.y); } }, 620);
    });
    canvas.addEventListener("pointermove", (e) => {
      if (down && Math.hypot(e.clientX - down.x, e.clientY - down.y) > 8) { down = null; clearTimeout(press); }
    });
    canvas.addEventListener("pointerup", (e) => {
      clearTimeout(press);
      const d = down;
      down = null;
      if (!d || d.button !== 2 || v.state.view.mode !== "sky") return;
      if (Math.hypot(e.clientX - d.x, e.clientY - d.y) <= 5) this.open(e.clientX, e.clientY);
    });
    v.on("viewmode", ({ mode }) => { if (mode !== "sky") this.close(); });
    window.addEventListener("keydown", (e) => { if (e.key === "Escape" && this.el) { this.close(); e.stopPropagation(); } }, true);
  }

  close() {
    this.abort?.abort();
    this.abort = null;
    this.el?.remove();
    this.el = null;
    if (this._outside) { this.ui.root.removeEventListener("pointerdown", this._outside, true); this._outside = null; }
  }

  /** List the SIMBAD objects around the sky direction under (clientX, clientY). */
  async open(clientX, clientY) {
    const v = this.viewer;
    const rect = v.canvas.getBoundingClientRect();
    const ray = v.renderer.camera.ray(clientX - rect.left, clientY - rect.top);
    const [l, b] = dirToLb(ray.dir);
    const [ra, dec] = galToIcrs(l, b);
    const radius = clamp(v.pose.fov * 0.02, 0.05, 1.2);
    this.close();
    const list = h("div", { class: "ov-identify-list", role: "list" });
    const status = h("div", { class: "ov-identify-status" }, "Searching SIMBAD…");
    const el = h("div", { class: "ov-pop ov-glass ov-identify", role: "dialog", "aria-label": "Sky objects" },
      h("div", { class: "ov-identify-head" },
        h("strong", null, "Sky objects"),
        h("span", { class: "ov-faint" }, `l ${l.toFixed(2)}°  b ${b >= 0 ? "+" : "−"}${Math.abs(b).toFixed(2)}°`),
        iconButton("close", "Close", () => this.close())),
      list, status);
    this.ui.ui.append(el);
    const rr = this.ui.root.getBoundingClientRect();
    const place = () => {
      const w = el.offsetWidth || 280, ht = el.offsetHeight || 140;
      el.style.left = `${clamp(clientX - rr.left + 12, 8, rr.width - w - 8)}px`;
      el.style.top = `${clamp(clientY - rr.top + 12, 8, rr.height - ht - 8)}px`;
    };
    place();
    this.el = el;
    this._outside = (e) => { if (!el.contains(e.target)) this.close(); };
    setTimeout(() => this.ui.root.addEventListener("pointerdown", this._outside, true), 0);
    const ctrl = new AbortController();
    this.abort = ctrl;
    const timer = setTimeout(() => ctrl.abort(), 9000);
    try {
      const adql = `SELECT TOP 40 b.main_id, b.otype, b.ra, b.dec, b.oid, b.sp_type, b.nbref AS refs, `
        + `DISTANCE(POINT('ICRS', b.ra, b.dec), POINT('ICRS', ${ra}, ${dec})) AS sep FROM basic b `
        + `WHERE CONTAINS(POINT('ICRS', b.ra, b.dec), CIRCLE('ICRS', ${ra}, ${dec}, ${radius})) = 1 AND b.nbref >= 1 ORDER BY refs DESC`;
      const payload = await tap(adql, ctrl.signal);
      if (this.el !== el) return;
      const shown = rankSimbadRows(payload?.data, 8);
      if (!shown.length) {
        status.textContent = `No catalogued objects within ${formatSeparation(radius)}.`;
        return;
      }
      status.textContent = `SIMBAD via CDS · ${formatSeparation(radius)} radius · click a row for details`;
      for (const entry of shown) list.append(this.row(entry));
      place();
    } catch (err) {
      if (this.el === el && !ctrl.signal.aborted) status.textContent = "SIMBAD lookup failed. Check the network connection.";
      else if (this.el === el) status.textContent = "SIMBAD did not answer in time.";
    } finally {
      clearTimeout(timer);
      if (this.abort === ctrl) this.abort = null;
    }
  }

  row(entry) {
    const item = h("div", { class: "ov-identify-row", role: "listitem", "data-open": "false", "data-tier": String(entry.tier) });
    const summary = h("button", { class: "ov-identify-summary", type: "button" },
      // SIMBAD prefixes common names with "NAME "; the link keeps the full id.
      h("span", { class: "ov-identify-name", title: entry.name }, entry.name.replace(/^NAME\s+/, "")),
      h("span", { class: "ov-identify-meta" }, [entry.otype, formatSeparation(entry.sep)].filter(Boolean).join(" · ")));
    summary.addEventListener("click", () => this.expand(item, entry));
    item.append(summary);
    return item;
  }

  async expand(item, entry) {
    const list = item.parentElement;
    const open = item.dataset.open === "true";
    for (const other of list?.querySelectorAll('.ov-identify-row[data-open="true"]') || []) {
      other.dataset.open = "false";
      other.querySelector(".ov-identify-detail")?.remove();
    }
    if (open) return;
    item.dataset.open = "true";
    const detail = h("div", { class: "ov-identify-detail" }, "Loading details…");
    item.append(detail);
    try {
      const adql = `SELECT b.otype, b.sp_type, b.plx_value, b.rvz_radvel, b.pmra, b.pmdec, b.nbref, f.V, f.B, f.G `
        + `FROM basic b LEFT JOIN allfluxes f ON f.oidref = b.oid WHERE b.oid = ${Number(entry.oid)}`;
      const payload = await tap(adql);
      if (!item.isConnected || item.dataset.open !== "true") return;
      const r = payload?.data?.[0] || [];
      const facts = [["Type", String(r[0] || entry.otype || "?")]];
      if (entry.spType || r[1]) facts.push(["Spectral type", String(r[1] || entry.spType)]);
      const mag = Number.isFinite(Number(r[7])) && r[7] != null ? ["V mag", Number(r[7]).toFixed(1)] : Number.isFinite(Number(r[9])) && r[9] != null ? ["G mag", Number(r[9]).toFixed(1)] : null;
      if (mag) facts.push(mag);
      const plx = Number(r[2]);
      if (Number.isFinite(plx) && plx > 0.01) {
        const pc = 1000 / plx;
        facts.push(["Distance", pc >= 1000 ? `${(pc / 1000).toFixed(2)} kpc` : `${pc.toFixed(0)} pc`]);
      }
      if (Number.isFinite(Number(r[4])) && Number.isFinite(Number(r[5])) && r[4] != null) facts.push(["Proper motion", `${Math.hypot(Number(r[4]), Number(r[5])).toFixed(1)} mas/yr`]);
      if (r[3] != null && Number.isFinite(Number(r[3]))) facts.push(["Radial velocity", `${Number(r[3]).toFixed(1)} km/s`]);
      const [gl, gb] = icrsToGal(entry.ra, entry.dec);
      facts.push(["Galactic", `${gl.toFixed(3)}°, ${gb.toFixed(3)}°`]);
      if (Number(r[6]) > 0) facts.push(["References", String(r[6])]);
      clear(detail);
      detail.append(h("div", { class: "ov-identify-facts" }, facts.map(([k, val]) => [h("span", { class: "ov-faint" }, k), h("span", null, val)])));
      detail.append(h("div", { class: "ov-identify-links" },
        h("a", { href: `https://simbad.cds.unistra.fr/simbad/sim-id?Ident=${encodeURIComponent(entry.name)}`, target: "_blank", rel: "noopener noreferrer" }, "Open in SIMBAD ↗"),
        h("button", { class: "ov-btn ov-btn--sm", type: "button", onclick: () => this.plugin.lookAt(gl, gb, { fov: Math.min(this.viewer.pose.fov, 12) }) }, icon("target"), "Centre")));
    } catch (_) {
      if (item.isConnected && item.dataset.open === "true") detail.textContent = "Detail lookup failed.";
    }
  }
}

/** Turn the Sky camera toward galactic (l, b), optionally narrowing the field. */
export function lookAtPose(pose, lDeg, bDeg, { fov } = {}) {
  const p = clonePose(pose);
  p.yaw = pose.yaw + angleDelta(pose.yaw, lDeg * DEG);
  p.pitch = clamp(bDeg * DEG, -Math.PI / 2 + 1e-3, Math.PI / 2 - 1e-3);
  if (Number.isFinite(fov)) p.fov = fov;
  return p;
}
