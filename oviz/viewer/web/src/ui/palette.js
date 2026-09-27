// Command palette (⌘K / "/"): fuzzy search over objects, layers, States and actions.

import { h, icon, kbd, fuzzyScore, clear } from "./dom.js";

export class Palette {
  constructor(ui) {
    this.ui = ui;
    this.viewer = ui.viewer;
    this.index = null;
    this.open = false;
    this.results = [];
    this.active = 0;
  }

  async buildObjectIndex() {
    if (this.index) return this.index;
    const v = this.viewer;
    const items = [];
    for (const trace of v.traces) {
      if (!trace.points || !trace.showInLegend) continue;
      const meta = await v.objectMeta(trace.key);
      const names = meta.name || [];
      for (let i = 0; i < names.length; i++) {
        const raw = names[i];
        if (!raw) continue;
        const name = raw.replace(/_/g, " ");
        const aliases = (meta.aliases?.[i] || "").split(",").map((s) => s.trim().replace(/_/g, " ")).filter((s) => s && s !== name);
        items.push({ trace: trace.key, index: i, name, aliases, traceName: trace.name, color: trace.color });
      }
    }
    this.index = items;
    return items;
  }

  toggle() {
    if (this.open) this.close();
    else this.show();
  }

  show(initial = "") {
    if (this.open) return;
    this.open = true;
    const ui = this.ui;
    this.scrim = h("div", { class: "ov-scrim", onclick: () => this.close() });
    this.input = h("input", {
      class: "ov-palette-input", type: "text", spellcheck: "false", autocomplete: "off",
      placeholder: "Search clusters, layers, views and actions…", "aria-label": "Command palette",
      role: "combobox", "aria-expanded": "true", "aria-controls": "ov-palette-list",
    });
    this.list = h("div", { class: "ov-palette-list", id: "ov-palette-list", role: "listbox" });
    this.el = h("div", { class: "ov-palette ov-glass", role: "dialog", "aria-label": "Command palette" },
      h("div", { class: "ov-palette-input-row" }, icon("search"), this.input, kbd("Esc")),
      this.list,
      h("div", { class: "ov-palette-foot" },
        h("span", null, kbd("↑"), kbd("↓"), "navigate"),
        h("span", null, kbd("↵"), "open"),
        h("span", null, kbd("⇧↵"), "select without flying")),
    );
    ui.root.append(this.scrim, this.el);
    this.input.value = initial;
    this.input.addEventListener("input", () => this.update());
    this.input.addEventListener("keydown", (e) => this.onKey(e));
    this.input.focus();
    this.update();
    this.buildObjectIndex().then(() => this.open && this.update());
  }

  close() {
    if (!this.open) return;
    this.open = false;
    this.scrim?.remove();
    this.el?.remove();
    this.ui.focusCanvas();
  }

  onKey(e) {
    e.stopPropagation();
    if (e.key === "Escape") { this.close(); e.preventDefault(); }
    else if (e.key === "ArrowDown") { this.setActive(this.active + 1); e.preventDefault(); }
    else if (e.key === "ArrowUp") { this.setActive(this.active - 1); e.preventDefault(); }
    else if (e.key === "Enter") {
      const item = this.results[this.active];
      if (item) { this.close(); item.run({ shift: e.shiftKey }); }
      e.preventDefault();
    }
  }

  setActive(i) {
    const n = this.results.length;
    if (!n) return;
    this.active = (i + n) % n;
    for (const [j, el] of this.itemEls.entries()) el.setAttribute("aria-selected", String(j === this.active));
    this.itemEls[this.active]?.scrollIntoView({ block: "nearest" });
  }

  commands() {
    return this.ui.commands().map((c) => ({ ...c, group: "Actions", score: 0 }));
  }

  update() {
    const q = this.input.value.trim();
    const results = [];
    const push = (item, text) => {
      const s = fuzzyScore(q, text);
      if (s > 0) results.push({ ...item, score: s });
    };
    for (const c of this.commands()) push(c, `${c.title} ${c.keywords || ""}`);
    for (const s of this.ui.stateItems()) push({ group: "Views", ...s }, s.title);
    for (const l of this.ui.layerItems()) push({ group: "Layers", ...l }, `${l.title} layer`);
    if (q && this.index) {
      const objs = [];
      for (const o of this.index) {
        let s = fuzzyScore(q, o.name);
        let via = "";
        if (s <= 0 || s < 30) {
          for (const a of o.aliases) {
            const sa = fuzzyScore(q, a);
            if (sa > s) { s = sa; via = a; }
          }
        }
        if (s > 0) objs.push({ o, s, via });
      }
      objs.sort((a, b) => b.s - a.s);
      for (const { o, s, via } of objs.slice(0, 40)) {
        results.push({
          group: "Objects",
          title: o.name,
          sub: via ? `${o.traceName} · matched “${via}”` : o.traceName,
          swatch: o.color,
          score: s + 5,
          hint: "Fly to",
          run: ({ shift } = {}) => this.ui.selectAndFly({ trace: o.trace, index: o.index }, !shift),
        });
      }
    }
    if (!q) {
      // Curated default view: a few actions, then views and layers.
      results.forEach((r) => { r.score = r.group === "Views" ? 3 : r.group === "Actions" ? (r.pinned ? 4 : 1) : 2; });
    }
    const order = { Objects: 0, Views: 1, Layers: 2, Actions: 3 };
    results.sort((a, b) => (q ? b.score - a.score : (order[a.group] - order[b.group]) || b.score - a.score));
    const top = q ? results.slice(0, 60) : results.filter((r) => r.group !== "Layers" || results.length < 30).slice(0, 40);
    // Group while keeping score order within groups.
    const grouped = [];
    const seen = new Set();
    for (const r of top) {
      if (seen.has(r.group)) continue;
      seen.add(r.group);
      grouped.push(...top.filter((x) => x.group === r.group));
    }
    this.results = grouped;
    this.renderList(q);
  }

  renderList(q) {
    clear(this.list);
    this.itemEls = [];
    if (!this.results.length) {
      this.list.append(h("div", { class: "ov-palette-empty" }, q ? `Nothing matches “${q}”.` : "Loading…"));
      return;
    }
    let group = "";
    this.results.forEach((r, i) => {
      if (r.group !== group) {
        group = r.group;
        this.list.append(h("div", { class: "ov-palette-group" }, group));
      }
      const iconEl = r.swatch
        ? h("span", { class: "ov-pi-icon" }, h("span", { class: "ov-swatch", style: { color: r.swatch } }))
        : h("span", { class: "ov-pi-icon" }, icon(r.icon || "sparkle"));
      const title = h("div", { class: "ov-pi-title" });
      highlight(title, r.title, q);
      const el = h("div", { class: "ov-palette-item", role: "option", "aria-selected": String(i === 0) },
        iconEl,
        h("div", { class: "ov-pi-text" }, title, r.sub ? h("div", { class: "ov-pi-sub" }, r.sub) : null),
        h("div", { class: "ov-pi-hint" }, r.shortcut ? r.shortcut.split(" ").map((k) => kbd(k)) : r.hint || ""),
      );
      el.addEventListener("pointermove", () => { if (this.active !== i) this.setActive(i); });
      el.addEventListener("click", (e) => { this.close(); r.run({ shift: e.shiftKey }); });
      this.itemEls.push(el);
      this.list.append(el);
    });
    this.active = 0;
  }
}

function highlight(el, text, q) {
  if (!q) { el.textContent = text; return; }
  const i = text.toLowerCase().indexOf(q.toLowerCase());
  if (i < 0) { el.textContent = text; return; }
  el.append(text.slice(0, i), h("mark", null, text.slice(i, i + q.length)), text.slice(i + q.length));
}
