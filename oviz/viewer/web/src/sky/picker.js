// Sky background picker: choose the survey behind Sky view the way a maps
// app picks its map type. Thumbnails of the whole sky, in wavelength order;
// a wavelength slider that crossfades between neighbouring surveys; what is
// in view, with a blend strength for each; and a search over every survey
// the CDS publishes.

import { h, icon, iconButton, clear } from "../ui/dom.js";
import { toggle } from "../ui/controls.js";
import { SKY_SURVEYS, surveyInfo, thumbnailUrl, formatWavelength, normSurveyId, loadRegistry, searchRegistry } from "./catalog.js";
import { layersInView, wavelengthStops, wavelengthPosition } from "./stack.js";

export class SkyPicker {
  constructor(plugin) {
    this.plugin = plugin;
    this.ui = plugin.ui;
    this.viewer = plugin.viewer;
    this.el = null;
    this.query = "";
    this.records = null; // CDS registry, once loaded
  }

  get open() {
    return !!this.el?.isConnected;
  }

  get thumbs() {
    return this.plugin.spec?.thumbs || null;
  }

  toggle(anchor) {
    if (this.open) { this.ui.closeMenu(); return; }
    this.show(anchor);
  }

  show(anchor) {
    const body = h("div", { class: "ov-skypick" });
    this.body = body;
    this.render();
    const pop = this.ui.menu(anchor || this.plugin.bgBtn, [{ body }], { align: "right", above: true, cls: "ov-pop--sky" });
    if (!pop) return;
    pop.setAttribute("role", "dialog");
    pop.setAttribute("aria-label", "Sky background");
    this.el = pop;
    // A menu takes focus on its first control; the tiles read better when
    // the current background has it.
    (body.querySelector('.ov-skytile-main[aria-pressed="true"]') || body.querySelector(".ov-skytile-main"))?.focus({ preventScroll: true });
    // Thumbnails may appear after the popover was placed.
    requestAnimationFrame(() => this.keepOnScreen());
  }

  keepOnScreen() {
    const pop = this.el;
    if (!pop) return;
    const rr = this.ui.root.getBoundingClientRect();
    const r = pop.getBoundingClientRect();
    if (r.bottom > rr.bottom - 8) pop.style.top = `${Math.max(8, r.top - rr.top - (r.bottom - rr.bottom + 8))}px`;
  }

  // ---------------------------------------------------------------- render

  render() {
    const body = this.body;
    if (!body) return;
    clear(body);
    const v = this.viewer;
    const st = v.state.sky;
    // Header: title and the on/off switch (the B key).
    this.switch = toggle({ label: "Show", checked: st.backgroundVisible, onChange: (x) => this.plugin.setBackground(x) });
    this.switch.classList.add("ov-skypick-switch");
    body.append(h("div", { class: "ov-skypick-head" },
      h("div", { class: "ov-skypick-title" }, "Sky background"),
      this.switch));
    this.hint = h("div", { class: "ov-skypick-hint", hidden: true });
    body.append(this.hint);
    // Search replaces the grid while it has a query.
    this.searchInput = h("input", {
      class: "ov-skypick-input", type: "search", spellcheck: "false", autocomplete: "off",
      placeholder: "Search every CDS survey: Planck, Hα, WISE…", "aria-label": "Search sky surveys",
      value: this.query,
    });
    this.searchInput.addEventListener("keydown", (e) => {
      e.stopPropagation();
      if (e.key === "Escape" && this.query) { e.preventDefault(); this.setQuery(""); }
      if (e.key === "Enter") this.enter();
    });
    this.searchInput.addEventListener("input", () => this.setQuery(this.searchInput.value, { keepInput: true }));
    this.searchInput.addEventListener("focus", () => this.ensureRegistry());
    this.content = h("div", { class: "ov-skypick-content" });
    body.append(this.content, h("div", { class: "ov-skypick-search" }, icon("search"), this.searchInput));
    body.append(h("a", { class: "ov-skypick-credit", href: "https://aladin.cds.unistra.fr/AladinLite/", target: "_blank", rel: "noopener noreferrer" }, "Imagery: Aladin Lite · CDS, Strasbourg"));
    this.renderContent();
    this.sync();
  }

  renderContent() {
    const c = this.content;
    clear(c);
    if (this.query.trim()) {
      c.append(this.resultsView());
      return;
    }
    c.append(this.gridView(), this.wavelengthView(), this.lensRow(), this.inViewList());
  }

  /** Every survey the picker offers: the figure's, then the curated ones, in wavelength order. */
  tiles() {
    const layers = this.viewer.state.sky.layers || [];
    const seen = new Set();
    const out = [];
    const add = (id, label) => {
      const info = surveyInfo(id, this.recordFor(id));
      const key = normSurveyId(info.id);
      if (seen.has(key)) return;
      seen.add(key);
      out.push({ ...info, id: info.known ? info.id : String(id), label: label && !info.known ? label : info.label });
    };
    for (const l of layers) add(l.survey || l.key, l.label);
    for (const s of SKY_SURVEYS) add(s.id);
    const order = (x) => (Number.isFinite(x.lambda) ? x.lambda : Infinity);
    return out.sort((a, b) => order(a) - order(b));
  }

  recordFor(id) {
    if (!this.records) return null;
    const n = normSurveyId(id);
    return this.records.find((r) => normSurveyId(r.id) === n) || null;
  }

  gridView() {
    this.tileEls = new Map();
    const grid = h("div", { class: "ov-skypick-grid", role: "group", "aria-label": "Backgrounds" });
    for (const t of this.tiles()) {
      const img = h("img", { alt: "", src: thumbnailUrl(t.id, this.thumbs), loading: "lazy", decoding: "async", draggable: "false" });
      img.addEventListener("error", () => { img.remove(); }, { once: true });
      const sub = [t.band, t.sky ? `${Math.round(t.sky * 100)}% of sky` : ""].filter(Boolean).join(" · ");
      const main = h("button", {
        class: "ov-skytile-main", type: "button",
        "data-tip": `${t.label}${t.band ? ` · ${t.band}` : ""}${Number.isFinite(t.lambda) ? ` · ${formatWavelength(t.lambda)}` : ""}${t.note ? `\n${t.note}` : ""}`,
        onclick: () => this.plugin.showOnly(t.id, { label: t.label }),
      },
      h("span", { class: "ov-skytile-map" }, img),
      h("span", { class: "ov-skytile-name" }, t.label),
      h("span", { class: "ov-skytile-band" }, sub || " "),
      h("span", { class: "ov-skytile-badge", "aria-hidden": "true" }));
      const add = iconButton("plus", `Blend ${t.label} over the background`, () => this.plugin.blendOnTop(t.id, { label: t.label }), { cls: "ov-skytile-add" });
      const tile = h("div", { class: "ov-skytile" }, main, add);
      this.tileEls.set(normSurveyId(t.id), { tile, main, add, t });
      grid.append(tile);
    }
    return grid;
  }

  wavelengthView() {
    const stops = (this.stops = wavelengthStops(this.viewer.state.sky.layers || []));
    const wrap = h("div", { class: "ov-skypick-lambda" });
    if (stops.length < 2) { wrap.hidden = true; return wrap; }
    this.lambdaRead = h("span", { class: "ov-field-value" });
    this.lambdaInput = h("input", { class: "ov-range ov-lambda-range", type: "range", min: "0", max: String((stops.length - 1) * 100), step: "1", "aria-label": "Wavelength" });
    this.lambdaInput.addEventListener("input", () => {
      this.plugin.setWavelength(Number(this.lambdaInput.value) / 100, stops);
    });
    this.lambdaInput.addEventListener("keydown", (e) => e.stopPropagation());
    const ticks = h("div", { class: "ov-lambda-ticks", "aria-hidden": "true" },
      stops.map((s, i) => h("span", { class: "ov-lambda-tick", style: { left: `${(i / (stops.length - 1)) * 100}%` }, title: `${s.label} · ${formatWavelength(s.lambda)}` })));
    wrap.append(
      h("div", { class: "ov-field-head" }, h("span", { class: "ov-field-label" }, "Wavelength"), this.lambdaRead),
      h("div", { class: "ov-lambda-track" }, h("div", { class: "ov-lambda-spectrum" }), ticks, this.lambdaInput),
      h("div", { class: "ov-lambda-ends" }, h("span", null, stops[0].band || "Short"), h("span", null, stops[stops.length - 1].band || "Long")),
    );
    return wrap;
  }

  /** The sky lens (classic aperture): another band inside a circle on the sky. */
  lensRow() {
    const lens = this.plugin.lens;
    if (!lens) return h("div", { hidden: true });
    const on = !!this.viewer.state.sky.lens;
    return h("div", { class: "ov-skypick-lens" },
      h("span", { class: "ov-skypick-lens-text" }, h("b", null, "Sky lens"), " · another wavelength in a circle you drag over the sky"),
      h("button", { class: "ov-btn ov-btn--sm", type: "button", "aria-pressed": String(on), onclick: () => { lens.toggle(); this.renderContent(); this.sync(); } }, on ? "Close" : "Open"));
  }

  inViewList() {
    this.viewHost = h("div", { class: "ov-skypick-view" });
    return this.viewHost;
  }

  renderInView() {
    const host = this.viewHost;
    if (!host) return;
    clear(host);
    const st = this.viewer.state.sky;
    const view = layersInView(st.layers || [], st.backgroundVisible).reverse(); // top first
    if (!view.length) return;
    host.append(h("div", { class: "ov-skypick-label" }, "In view"));
    view.forEach((l, i) => {
      const isBase = i === view.length - 1;
      const pct = h("span", { class: "ov-skyrow-pct" }, `${Math.round((l.opacity ?? 1) * 100)}%`);
      const range = h("input", { class: "ov-range", type: "range", min: "0", max: "100", step: "1", value: String(Math.round((l.opacity ?? 1) * 100)), "aria-label": `${l.label || l.key} strength` });
      const paint = () => range.style.setProperty("--p", `${range.value}%`);
      paint();
      range.addEventListener("input", () => {
        paint();
        pct.textContent = `${range.value}%`;
        this.plugin.setLayerOpacity(l.key, Number(range.value) / 100, { quiet: true });
      });
      range.addEventListener("change", () => this.plugin.changed());
      range.addEventListener("keydown", (e) => e.stopPropagation());
      host.append(h("div", { class: "ov-skyrow" },
        h("img", { class: "ov-skyrow-map", alt: "", src: thumbnailUrl(l.survey || l.key, this.thumbs), loading: "lazy" }),
        h("div", { class: "ov-skyrow-text" },
          h("span", { class: "ov-skyrow-name", title: l.survey || l.key }, l.label || l.key),
          h("span", { class: "ov-skyrow-role" }, isBase ? "Background" : "Blended over it")),
        range, pct,
        iconButton("close", `Remove ${l.label || l.key} from view`, () => this.plugin.hideLayer(l.key), { cls: "ov-skyrow-x" })));
    });
  }

  // ---------------------------------------------------------------- search

  setQuery(q, { keepInput = false } = {}) {
    this.query = String(q || "");
    if (!keepInput && this.searchInput) this.searchInput.value = this.query;
    this.renderContent();
    this.sync();
    if (this.query.trim()) this.ensureRegistry();
  }

  ensureRegistry() {
    if (this.records || this._loading) return;
    this._loading = loadRegistry()
      .then((records) => { this.records = records; })
      .catch((err) => { this.registryError = err.message || "The CDS survey list could not be loaded"; })
      .finally(() => {
        this._loading = null;
        if (this.open && this.query.trim()) { this.renderContent(); this.sync(); }
      });
  }

  resultsView() {
    const wrap = h("div", { class: "ov-skypick-results", role: "list" });
    const q = this.query.trim();
    this.resultEls = new Map();
    if (!this.records) {
      wrap.append(h("div", { class: "ov-skypick-empty" }, this.registryError
        ? `${this.registryError}. Press Enter to add “${q}” as a HiPS ID or URL.`
        : "Loading the CDS survey list…"));
      return wrap;
    }
    const found = searchRegistry(this.records, q, 40);
    if (!found.length) {
      wrap.append(h("div", { class: "ov-skypick-empty" }, looksLikeHips(q) ? `No survey matches. Press Enter to add “${q}”.` : `No survey matches “${q}”.`));
      return wrap;
    }
    wrap.append(h("div", { class: "ov-skypick-label" }, `${found.length === 40 ? "Top 40" : found.length} of ${this.records.length.toLocaleString()} surveys`));
    for (const r of found) {
      const info = surveyInfo(r.id, r);
      const meta = [info.band, Number.isFinite(r.sky) && r.sky < 0.995 ? `${Math.max(1, Math.round(r.sky * 100))}% of sky` : "all sky", r.id].filter(Boolean).join(" · ");
      const main = h("button", { class: "ov-skyresult-main", type: "button", title: r.title || r.id, onclick: () => this.plugin.showOnly(r.id, { label: info.label }) },
        h("img", { class: "ov-skyresult-map", alt: "", src: thumbnailUrl(r.id, this.thumbs), loading: "lazy", decoding: "async" }),
        h("span", { class: "ov-skyresult-text" },
          h("span", { class: "ov-skyresult-title" }, r.title || r.id),
          h("span", { class: "ov-skyresult-meta" }, meta)));
      const add = iconButton("plus", `Blend ${info.label} over the background`, () => this.plugin.blendOnTop(r.id, { label: info.label }), { cls: "ov-skyresult-add" });
      const row = h("div", { class: "ov-skyresult", role: "listitem" }, main, add);
      this.resultEls.set(normSurveyId(r.id), { row, main });
      wrap.append(row);
    }
    return wrap;
  }

  /** Enter in the search field: the top result, or a typed HiPS ID or URL. */
  enter() {
    const q = this.query.trim();
    if (!q) return;
    const top = this.records ? searchRegistry(this.records, q, 1)[0] : null;
    if (top && (!looksLikeHips(q) || normSurveyId(top.id) === normSurveyId(q))) {
      this.plugin.showOnly(top.id, { label: surveyInfo(top.id, top).label });
    } else if (looksLikeHips(q)) {
      this.plugin.showOnly(q);
    } else {
      this.ui.toast("Type a survey name, a HiPS ID such as P/2MASS/color, or an https:// URL");
    }
  }

  // ---------------------------------------------------------------- sync

  /** Reflect the current stack without rebuilding (called on every change). */
  sync() {
    if (!this.body) return;
    const v = this.viewer;
    const st = v.state.sky;
    const layers = st.layers || [];
    this.switch?.set(st.backgroundVisible);
    this.body.dataset.off = String(!st.backgroundVisible);
    const view = layersInView(layers, st.backgroundVisible);
    const role = new Map(view.map((l, i) => [normSurveyId(surveyInfo(l.survey || l.key).id), { i, l }]));
    const failed = this.plugin.sky?.failed;
    const failedIds = new Set([...(failed?.keys() || [])].map((k) => {
      const l = layers.find((x) => x.key === k);
      return normSurveyId(surveyInfo(l ? l.survey || l.key : k).id);
    }));
    for (const [id, r] of this.tileEls || []) {
      const hit = role.get(id);
      const state = failedIds.has(id) ? "failed" : !hit ? "off" : hit.i === 0 ? "base" : "blend";
      r.tile.dataset.state = state;
      r.main.setAttribute("aria-pressed", String(state === "base"));
      const badge = r.main.querySelector(".ov-skytile-badge");
      badge.textContent = state === "blend" ? `${Math.round((hit.l.opacity ?? 1) * 100)}%` : state === "failed" ? "!" : "";
      r.add.hidden = state === "base";
    }
    for (const [id, r] of this.resultEls || []) {
      const hit = role.get(id);
      r.row.dataset.state = !hit ? "off" : hit.i === 0 ? "base" : "blend";
    }
    if (this.lambdaInput && this.stops) {
      const pos = wavelengthPosition(layers, this.stops, st.backgroundVisible);
      this.lambdaInput.parentElement.parentElement.dataset.free = String(pos == null);
      if (pos != null && document.activeElement !== this.lambdaInput) this.lambdaInput.value = String(Math.round(pos * 100));
      const p = pos == null ? null : pos;
      const s = p == null ? null : this.stops[Math.round(p)];
      this.lambdaRead.textContent = s ? `${s.label} · ${formatWavelength(s.lambda)}` : "Drag to sweep the spectrum";
      const max = Number(this.lambdaInput.max) || 1;
      this.lambdaInput.style.setProperty("--p", `${(Number(this.lambdaInput.value) / max) * 100}%`);
    }
    this.renderInView();
    // Today's sky: Aladin imagery is the present-day sky, shown at t = 0.
    const t = v.timeline.time;
    const away = st.backgroundVisible && Math.abs(t) >= 1;
    this.hint.hidden = !away;
    if (away) this.hint.textContent = "The imagery is today's sky: it shows at t = 0 (press 0).";
  }
}

function looksLikeHips(q) {
  return /^(https?:\/\/\S+|[A-Za-z0-9_.-]+\/[A-Za-z0-9_./+-]+)$/.test(String(q || "").trim());
}
