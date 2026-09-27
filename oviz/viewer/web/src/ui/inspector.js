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
    v.on("time", () => this.current && this.refresh());
  }

  show(hit) {
    this.current = hit;
    this.el.dataset.open = hit ? "true" : "false";
    if (hit) this.refresh();
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
    const t = v.timeline.time;
    clear(this.body);
    const kicker = h("div", { class: "ov-insp-kicker" }, h("span", { class: "ov-swatch", style: { color: d.color } }), d.traceName);
    const hero = h("div", { class: "ov-insp-hero" }, kicker, h("div", { class: "ov-insp-name" }, d.name));
    if (d.aliases.length) hero.append(h("div", { class: "ov-insp-aliases" }, `Also: ${d.aliases.slice(0, 6).join(", ")}`));
    const facts = h("div", { class: "ov-facts" });
    const fact = (label, value, unit, wide) => facts.append(h("div", { class: `ov-fact${wide ? " ov-fact--wide" : ""}` },
      h("div", { class: "ov-fact-label" }, label),
      h("div", { class: "ov-fact-value" }, value, unit ? h("small", null, unit) : null)));
    if (d.ageNow != null) {
      fact("Age today", fmt(d.ageNow, d.ageNow >= 100 ? 0 : 1), "Myr");
      const ageAt = d.ageAt;
      fact(Math.abs(t) < 1e-6 ? "Born" : "Age at t", Math.abs(t) < 1e-6 ? fmt(-d.ageNow, d.ageNow >= 100 ? 0 : 1) : ageAt < 0 ? "unborn" : fmt(ageAt, 1), Math.abs(t) < 1e-6 || ageAt >= 0 ? "Myr" : "");
    }
    if (d.nStars != null) fact("Members", fmtInt(d.nStars), "stars");
    if (d.dist != null) fact("Distance", fmt(d.dist, d.dist >= 1000 ? 0 : 1), "pc");
    if (d.l != null && d.b != null) fact("Galactic l, b", `${fmt(d.l, 2)}°, ${fmt(d.b, 2)}°`, "", true);
    if (d.ra != null && d.dec != null) fact("ICRS α, δ", `${fmt(d.ra, 3)}°, ${fmt(d.dec, 3)}°`, "", true);
    if (d.position) {
      const [x, y, z] = d.position;
      fact(Math.abs(t) < 1e-6 ? "Position" : `Position at ${this.ui.formatTime(t)}`, `${fmt(x, 0)}, ${fmt(y, 0)}, ${fmt(z, 0)}`, "pc", true);
    }
    const following = this.ui.following && this.ui.following.trace === hit.trace && this.ui.following.index === hit.index;
    const trailOn = this.ui.hasTrail(hit);
    const actions = h("div", { class: "ov-insp-actions" },
      h("button", { class: "ov-btn", type: "button", onclick: () => v.flyToObject(hit.trace, hit.index) }, icon("target"), "Fly to"),
      h("button", { class: "ov-btn", type: "button", "aria-pressed": String(!!following), onclick: () => { this.ui.toggleFollow(hit); this.refresh(); } }, icon("follow"), following ? "Following" : "Follow"),
      h("button", { class: "ov-btn", type: "button", "aria-pressed": String(trailOn), disabled: v.timeline.count <= 1 ? true : null, onclick: () => { this.ui.toggleTrail(hit); this.refresh(); } }, icon("trail"), "Orbit trail"),
      h("button", { class: "ov-btn", type: "button", onclick: () => this.copy(d) }, icon("copy"), "Copy"),
    );
    this.body.append(hero, facts, actions);
    if (d.hover) {
      const legacy = h("div", { class: "ov-insp-legacy" });
      legacy.innerHTML = sanitize(d.hover);
      this.body.append(legacy);
    }
  }

  renderMember(hit) {
    clear(this.body);
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
