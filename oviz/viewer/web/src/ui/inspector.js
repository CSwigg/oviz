// Selection inspector + hover card.

import { h, icon, iconButton, clear, fmt, fmtInt, copyText } from "./dom.js";
import { formatDistance } from "../core/math.js";

export class Inspector {
  constructor(ui) {
    this.ui = ui;
    const v = (this.viewer = ui.viewer);
    this.el = h("section", { class: "ov-panel ov-inspector ov-glass ov-chrome", "aria-label": "Selection", "data-open": "false" });
    this.title = h("span", { class: "ov-panel-title" }, "Selection");
    this.head = h("div", { class: "ov-panel-head" }, this.title,
      iconButton("close", "Clear selection", () => ui.select(null), { shortcut: "Esc" }));
    this.body = h("div", { class: "ov-panel-body", style: { padding: "0" } });
    this.el.append(this.head, this.body);
    this.current = null;
    this.dyn = null;
    this.track = null;
    // Time only changes a few values: patch those instead of rebuilding.
    v.on("time", () => this.current && this.updateDynamic());
    this.bindTilt();
  }

  show(hit) {
    this.current = hit;
    this.el.dataset.open = hit ? "true" : "false";
    if (hit) this.refresh();
    else { this.dyn = null; this.track = null; }
  }

  /** The card leans a few degrees toward the pointer, with a moving sheen. */
  bindTilt() {
    const el = this.el;
    const fine = window.matchMedia?.("(hover: hover) and (pointer: fine)");
    const reduce = window.matchMedia?.("(prefers-reduced-motion: reduce)");
    const set = (x, y) => {
      el.style.setProperty("--ov-tilt-y", `${((x - 0.5) * 5).toFixed(2)}deg`);
      el.style.setProperty("--ov-tilt-x", `${((0.5 - y) * 4).toFixed(2)}deg`);
      el.style.setProperty("--ov-sheen-x", `${(x * 100).toFixed(1)}%`);
      el.style.setProperty("--ov-sheen-y", `${(y * 100).toFixed(1)}%`);
    };
    el.addEventListener("pointermove", (e) => {
      if (!fine?.matches || reduce?.matches || this.ui.narrow || e.pointerType !== "mouse") return;
      const r = el.getBoundingClientRect();
      set(clamp01((e.clientX - r.left) / r.width), clamp01((e.clientY - r.top) / r.height));
      el.dataset.tilting = "true";
    });
    el.addEventListener("pointerleave", () => {
      el.dataset.tilting = "false";
      set(0.5, 0.5);
    });
  }

  async refresh() {
    const hit = this.current;
    if (!hit) return;
    const v = this.viewer;
    if (hit.kind === "member") return this.renderMember(hit);
    await v.objectMeta(hit.trace);
    if (hit !== this.current) return;
    const d = v.describeObject(hit.trace, hit.index);
    if (!d) return;
    clear(this.body);
    const kicker = h("div", { class: "ov-insp-kicker" }, h("span", { class: "ov-swatch", style: { color: d.color } }), d.traceName);
    const hero = h("div", { class: "ov-insp-hero" }, kicker, h("div", { class: "ov-insp-name" }, d.name));
    if (d.aliases.length) hero.append(h("div", { class: "ov-insp-aliases" }, `Also: ${d.aliases.slice(0, 6).join(", ")}`));
    // Key facts up front; coordinates fold away until asked for.
    const facts = h("div", { class: "ov-facts" });
    const coords = h("div", { class: "ov-facts ov-facts--coords" });
    this.dyn = { hit, facts: new Map() };
    for (const f of this.factList(d)) {
      const label = h("div", { class: "ov-fact-label" }, f.label);
      const value = h("div", { class: "ov-fact-value" });
      (f.wide ? coords : facts).append(h("div", { class: `ov-fact${f.wide ? " ov-fact--wide" : ""}` }, label, value));
      this.dyn.facts.set(f.key, { label, value });
      setFact(this.dyn.facts.get(f.key), f);
    }
    const more = coords.childElementCount
      ? h("details", { class: "ov-insp-more", open: this.coordsOpen ? true : null },
        h("summary", null, h("span", null, "Coordinates"), icon("chevronDown")), coords)
      : null;
    more?.addEventListener("toggle", () => { this.coordsOpen = more.open; });
    this.track = DistanceTrack.create(this.ui, hit, d.color);
    const following = this.ui.following && this.ui.following.trace === hit.trace && this.ui.following.index === hit.index;
    const trailOn = this.ui.hasTrail(hit);
    const actions = h("div", { class: "ov-insp-actions" },
      h("button", { class: "ov-btn", type: "button", onclick: () => v.flyToObject(hit.trace, hit.index) }, icon("target"), "Fly to"),
      h("button", { class: "ov-btn", type: "button", "aria-pressed": String(!!following), onclick: () => { this.ui.toggleFollow(hit); this.refresh(); } }, icon("follow"), following ? "Following" : "Follow"),
      h("button", { class: "ov-btn", type: "button", "aria-pressed": String(trailOn), disabled: v.timeline.count <= 1 ? true : null, onclick: () => { this.ui.toggleTrail(hit); this.refresh(); } }, icon("trail"), "Orbit trail"),
      h("button", { class: "ov-btn", type: "button", onclick: () => this.copy(v.describeObject(hit.trace, hit.index) || d) }, icon("copy"), "Copy"),
      h("button", { class: "ov-btn", type: "button", onclick: () => this.ui.notes?.addForObject(hit) }, icon("edit"), "Add note"),
    );
    this.body.append(hero);
    if (this.track) this.body.append(this.track.el);
    this.body.append(facts);
    if (more) this.body.append(more);
    this.body.append(actions);
    if (d.ra != null && d.dec != null) {
      const coord = `${d.ra.toFixed(5)} ${d.dec >= 0 ? "+" : ""}${d.dec.toFixed(5)}`;
      const simbad = `https://simbad.cds.unistra.fr/simbad/sim-coo?Coord=${encodeURIComponent(coord)}&CooFrame=ICRS&Radius=10&Radius.unit=arcmin`;
      const aladin = `https://aladin.cds.unistra.fr/AladinLite/?target=${encodeURIComponent(coord)}&fov=1.5&survey=P%2FDSS2%2Fcolor`;
      this.body.append(h("div", { class: "ov-insp-links" },
        h("span", { class: "ov-fact-label" }, "Look up"),
        h("a", { class: "ov-btn ov-btn--sm ov-btn--ghost", href: simbad, target: "_blank", rel: "noopener noreferrer" }, "SIMBAD ↗"),
        h("a", { class: "ov-btn ov-btn--sm ov-btn--ghost", href: aladin, target: "_blank", rel: "noopener noreferrer" }, "Aladin ↗")));
    }
    if (d.hover) {
      const legacy = h("div", { class: "ov-insp-legacy" });
      legacy.innerHTML = sanitize(d.hover);
      this.body.append(legacy);
    }
    this.track?.draw();
  }

  /** Facts shown for an object; `key` identifies the ones that follow time. */
  factList(d) {
    const t = this.viewer.timeline.time;
    const now = Math.abs(t) < 1e-6;
    const out = [];
    if (d.ageNow != null) {
      const digits = d.ageNow >= 100 ? 0 : 1;
      out.push({ key: "age", label: "Age today", value: fmt(d.ageNow, digits), unit: "Myr" });
      const ageAt = d.ageAt;
      out.push({
        key: "ageAt",
        label: now ? "Born" : "Age at t",
        value: now ? fmt(-d.ageNow, digits) : ageAt < 0 ? "unborn" : fmt(ageAt, 1),
        unit: now || ageAt >= 0 ? "Myr" : "",
      });
    }
    if (d.nStars != null) out.push({ key: "members", label: "Members", value: fmtInt(d.nStars), unit: "stars" });
    if (d.dist != null) out.push({ key: "dist", label: "Distance", value: fmt(d.dist, d.dist >= 1000 ? 0 : 1), unit: "pc" });
    if (d.l != null && d.b != null) out.push({ key: "lb", label: "Galactic l, b", value: `${fmt(d.l, 2)}°, ${fmt(d.b, 2)}°`, wide: true });
    if (d.ra != null && d.dec != null) out.push({ key: "radec", label: "ICRS α, δ", value: `${fmt(d.ra, 3)}°, ${fmt(d.dec, 3)}°`, wide: true });
    if (d.position) {
      const [x, y, z] = d.position;
      out.push({ key: "pos", label: now ? "Position" : `Position at ${this.ui.formatTime(t)}`, value: `${fmt(x, 0)}, ${fmt(y, 0)}, ${fmt(z, 0)}`, unit: "pc", wide: true });
    }
    return out;
  }

  updateDynamic() {
    const hit = this.current;
    if (!hit || hit.kind === "member" || !this.dyn || this.dyn.hit !== hit) return;
    const d = this.viewer.describeObject(hit.trace, hit.index);
    if (!d) return;
    for (const f of this.factList(d)) {
      const slot = this.dyn.facts.get(f.key);
      if (slot) setFact(slot, f);
    }
    this.track?.draw();
  }

  renderMember(hit) {
    clear(this.body);
    this.dyn = null;
    this.track = null;
    const m = hit.member;
    const hero = h("div", { class: "ov-insp-hero" },
      h("div", { class: "ov-insp-kicker" }, h("span", { class: "ov-swatch", style: { color: m.color } }), `Member of ${m.cluster}`),
      h("div", { class: "ov-insp-name" }, `Star ${hit.index + 1}`));
    const facts = h("div", { class: "ov-facts" });
    const add = (label, value, unit, wide) => facts.append(h("div", { class: `ov-fact${wide ? " ov-fact--wide" : ""}` },
      h("div", { class: "ov-fact-label" }, label), h("div", { class: "ov-fact-value" }, value, unit ? h("small", null, unit) : null)));
    add("Distance", fmt(m.dist, 1), "pc");
    add("Tangential v", fmt(m.speed, 2), "km/s");
    add("Galactic l, b", `${fmt(m.l, 3)}°, ${fmt(m.b, 3)}°`, "", true);
    this.body.append(hero, facts);
  }

  copy(d) {
    const lines = [d.name];
    if (d.l != null) lines.push(`l, b = ${d.l.toFixed(4)}, ${d.b.toFixed(4)} deg`);
    if (d.ra != null) lines.push(`RA, Dec = ${d.ra.toFixed(5)}, ${d.dec.toFixed(5)} deg`);
    if (d.dist != null) lines.push(`distance = ${d.dist.toFixed(1)} pc`);
    if (d.ageNow != null) lines.push(`age = ${d.ageNow.toFixed(2)} Myr`);
    copyText(lines.join("\n")).then(() => this.ui.toast("Copied to clipboard", { icon: icon("copy") }));
  }
}

function setFact(slot, f) {
  if (slot.label.textContent !== f.label) slot.label.textContent = f.label;
  const text = `${f.value}${f.unit ? ` ${f.unit}` : ""}`;
  if (slot.value.dataset.text === text) return;
  slot.value.dataset.text = text;
  slot.value.replaceChildren(f.value, f.unit ? h("small", null, f.unit) : "");
}

function clamp01(x) {
  return x < 0 ? 0 : x > 1 ? 1 : x;
}

/**
 * Distance from the Sun across the whole timeline, with a cursor at the
 * current time. Dragging along the chart scrubs time.
 */
class DistanceTrack {
  static create(ui, hit, color) {
    const v = ui.viewer;
    const tl = v.timeline;
    const trace = v.traceByKey.get(hit.trace);
    const F = trace?.points?.position?.frames || 1;
    if (tl.count <= 1 || F <= 1 || /^sun$/i.test(String(trace.name).trim())) return null;
    const dist = new Float64Array(F).fill(NaN);
    const eye = v.skyEye();
    let n = 0;
    for (let f = 0; f < F; f++) {
      const p = v.objectPosition(hit.trace, hit.index, f);
      if (!p) continue;
      dist[f] = Math.hypot(p[0] - eye[0], p[1] - eye[1], p[2] - eye[2]);
      n++;
    }
    if (n < 2) return null;
    return new DistanceTrack(ui, hit, dist, color);
  }

  constructor(ui, hit, dist, color) {
    this.ui = ui;
    this.hit = hit;
    this.dist = dist;
    this.color = color;
    const tl = ui.viewer.timeline;
    let lo = Infinity, hi = -Infinity;
    for (const d of dist) if (d === d) { lo = Math.min(lo, d); hi = Math.max(hi, d); }
    const pad = Math.max((hi - lo) * 0.12, 0.5);
    this.lo = lo - pad;
    this.hi = hi + pad;
    this.readout = h("span", { class: "ov-insp-track-read" });
    this.canvas = h("canvas", { class: "ov-spark", role: "slider", "aria-label": "Distance from the Sun over time; drag to change time" });
    this.el = h("div", { class: "ov-insp-track" },
      h("div", { class: "ov-insp-track-head" }, h("span", { class: "ov-fact-label" }, "Distance from the Sun"), this.readout),
      this.canvas,
      h("div", { class: "ov-insp-track-axis" }, h("span", null, ui.formatTime(tl.times[0])), h("span", null, ui.formatTime(tl.times[tl.count - 1]))));
    this.bindScrub();
  }

  bindScrub() {
    const c = this.canvas;
    const tl = this.ui.viewer.timeline;
    let dragging = false;
    const seek = (e) => {
      const r = c.getBoundingClientRect();
      tl.setFrame(clamp01((e.clientX - r.left) / r.width) * (tl.count - 1));
    };
    c.addEventListener("pointerdown", (e) => {
      if (e.button !== 0) return;
      dragging = true;
      tl.pause();
      c.setPointerCapture(e.pointerId);
      seek(e);
    });
    c.addEventListener("pointermove", (e) => { if (dragging) seek(e); });
    const end = () => { dragging = false; };
    c.addEventListener("pointerup", end);
    c.addEventListener("pointercancel", end);
  }

  draw() {
    const c = this.canvas;
    const w = c.clientWidth, hgt = c.clientHeight;
    if (!w || !hgt) return;
    const dpr = Math.min(window.devicePixelRatio || 1, 2);
    if (c.width !== Math.round(w * dpr) || c.height !== Math.round(hgt * dpr)) {
      c.width = Math.round(w * dpr);
      c.height = Math.round(hgt * dpr);
    }
    const ctx = c.getContext("2d");
    const de = document.documentElement.dataset;
    const theme = `${de.ovizTheme || "dark"}/${de.ovizSkin || ""}`;
    if (this._theme !== theme) {
      const css = getComputedStyle(this.ui.root);
      const get = (name, fb) => css.getPropertyValue(name).trim() || fb;
      this._colors = {
        accent: get("--ov-accent", "#f4c46a"),
        line: get("--ov-line-strong", "rgba(255,255,255,0.14)"),
        text: get("--ov-text", "#fff"),
        surface: get("--ov-surface-solid", "#0e1118"),
      };
      this._theme = theme;
    }
    const { accent, line, text, surface } = this._colors;
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.clearRect(0, 0, w, hgt);
    const tl = this.ui.viewer.timeline;
    const F = this.dist.length;
    const X = (f) => (f / (F - 1)) * (w - 8) + 4;
    const Y = (d) => 5 + (1 - (d - this.lo) / (this.hi - this.lo)) * (hgt - 10);
    // Present day.
    const z = tl.zeroFrame?.() ?? -1;
    if (z > 0 && z < tl.count - 1) {
      ctx.strokeStyle = line;
      ctx.setLineDash([2, 3]);
      ctx.beginPath();
      ctx.moveTo(X(z * (F - 1) / (tl.count - 1)), 2);
      ctx.lineTo(X(z * (F - 1) / (tl.count - 1)), hgt - 2);
      ctx.stroke();
      ctx.setLineDash([]);
    }
    const path = () => {
      ctx.beginPath();
      let pen = false;
      for (let f = 0; f < F; f++) {
        const d = this.dist[f];
        if (d !== d) { pen = false; continue; }
        if (pen) ctx.lineTo(X(f), Y(d));
        else { ctx.moveTo(X(f), Y(d)); pen = true; }
      }
    };
    // Soft area under the curve, then the line.
    path();
    ctx.lineTo(X(F - 1), hgt);
    ctx.lineTo(X(0), hgt);
    ctx.closePath();
    ctx.globalAlpha = 0.1;
    ctx.fillStyle = accent;
    ctx.fill();
    ctx.globalAlpha = 1;
    path();
    ctx.strokeStyle = accent;
    ctx.lineWidth = 1.6;
    ctx.lineJoin = "round";
    ctx.stroke();
    // Cursor at the current time.
    const fr = (tl.frame / Math.max(tl.count - 1, 1)) * (F - 1);
    const d = sampleSeries(this.dist, fr);
    const cx = X(fr);
    ctx.strokeStyle = line;
    ctx.lineWidth = 1;
    ctx.beginPath();
    ctx.moveTo(cx, 0);
    ctx.lineTo(cx, hgt);
    ctx.stroke();
    if (d === d) {
      ctx.fillStyle = surface;
      ctx.strokeStyle = text;
      ctx.lineWidth = 1.6;
      ctx.beginPath();
      ctx.arc(cx, Y(d), 3.6, 0, Math.PI * 2);
      ctx.fill();
      ctx.stroke();
    }
    const t = tl.time;
    const when = Math.abs(t) < 1e-6 ? "today" : `at ${this.ui.formatTime(t)}`;
    this.readout.replaceChildren(d === d ? `${fmt(d, d >= 1000 ? 0 : 1)} pc` : "—", h("small", null, when));
  }
}

function sampleSeries(arr, f) {
  const i = Math.floor(f);
  const j = Math.min(i + 1, arr.length - 1);
  const a = arr[i], b = arr[j];
  if (a !== a) return b;
  if (b !== b) return a;
  return a + (b - a) * (f - i);
}

/** Allow only simple formatting tags from legacy hover HTML. */
function sanitize(html) {
  const tpl = document.createElement("template");
  tpl.innerHTML = String(html);
  const allowed = new Set(["B", "I", "EM", "STRONG", "BR", "SPAN", "SUB", "SUP", "SMALL"]);
  const walk = (node) => {
    for (const child of [...node.childNodes]) {
      if (child.nodeType === 1) {
        if (!allowed.has(child.tagName)) {
          child.replaceWith(document.createTextNode(child.textContent));
          continue;
        }
        for (const attr of [...child.attributes]) child.removeAttribute(attr.name);
        walk(child);
      }
    }
  };
  walk(tpl.content);
  return tpl.innerHTML;
}

export class HoverCard {
  constructor(ui) {
    this.ui = ui;
    this.el = h("div", { class: "ov-hovercard ov-glass", role: "status" });
    this.swatch = h("span", { class: "ov-swatch" });
    this.name = h("span");
    this.sub = h("small");
    this.el.append(this.swatch, this.name, this.sub);
  }

  show(hit, x, y) {
    const v = this.ui.viewer;
    if (!hit) {
      this.el.dataset.show = "false";
      return;
    }
    let name, sub, color;
    if (hit.kind === "member") {
      name = `${hit.member.cluster} member`;
      sub = `${fmt(hit.member.dist, 0)} pc`;
      color = hit.member.color;
    } else {
      const meta = v.objectMetaSync(hit.trace);
      if (!meta) v.objectMeta(hit.trace).then(() => this.ui.refreshHover());
      const d = v.describeObject(hit.trace, hit.index);
      if (!d) return;
      name = d.name;
      const parts = [];
      if (d.ageAt != null) parts.push(d.ageAt < 0 ? "unborn" : `${fmt(d.ageAt, d.ageAt >= 100 ? 0 : 1)} Myr`);
      if (d.dist != null && Math.abs(v.timeline.time) < 1e-6) parts.push(formatDistance(d.dist));
      sub = parts.join(" · ");
      color = d.color;
    }
    this.swatch.style.color = color;
    this.name.textContent = name;
    this.sub.textContent = sub;
    this.el.dataset.show = "true";
    const W = this.ui.root.clientWidth;
    const w = this.el.offsetWidth;
    const px = x + 16 + w > W - 8 ? x - w - 14 : x + 16;
    this.el.style.transform = `translate3d(${px}px, ${y + 14}px, 0)`;
  }
}
