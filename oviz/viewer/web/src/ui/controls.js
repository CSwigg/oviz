// Reusable form controls.

import { h } from "./dom.js";

function setProgress(input) {
  const min = Number(input.min), max = Number(input.max), v = Number(input.value);
  const p = max > min ? ((v - min) / (max - min)) * 100 : 0;
  input.style.setProperty("--p", `${p}%`);
}

/**
 * Slider with label + live value. `scale: "log"` maps the slider linearly in
 * log space (for opacity/alpha coefficients spanning decades).
 */
export function slider({ label, min, max, step = 0.01, value, format = (v) => v.toFixed(2), onInput, onChange, scale = "linear" }) {
  const log = scale === "log";
  const toSlider = (v) => (log ? Math.log(Math.max(v, min)) : v);
  const fromSlider = (s) => (log ? Math.exp(s) : s);
  const input = h("input", {
    class: "ov-range", type: "range",
    min: log ? Math.log(min) : min, max: log ? Math.log(max) : max,
    step: log ? (Math.log(max) - Math.log(min)) / 400 : step,
    value: toSlider(value),
    "aria-label": label,
  });
  const out = h("span", { class: "ov-field-value" }, format(value));
  setProgress(input);
  input.addEventListener("input", () => {
    const v = fromSlider(Number(input.value));
    out.textContent = format(v);
    setProgress(input);
    onInput?.(v);
  });
  input.addEventListener("change", () => onChange?.(fromSlider(Number(input.value))));
  input.addEventListener("dblclick", () => {
    // Double-click resets to the value the control was created with.
    input.value = toSlider(value);
    input.dispatchEvent(new Event("input"));
    input.dispatchEvent(new Event("change"));
  });
  const field = h("div", { class: "ov-field" },
    h("div", { class: "ov-field-head" }, h("span", { class: "ov-field-label" }, label), out),
    input,
  );
  field.set = (v) => {
    input.value = toSlider(v);
    out.textContent = format(v);
    setProgress(input);
  };
  return field;
}

export function toggle({ label, hint, checked, onChange }) {
  const input = h("input", { type: "checkbox", role: "switch" });
  input.checked = !!checked;
  input.addEventListener("change", () => onChange?.(input.checked));
  const el = h("label", { class: "ov-switch" },
    h("span", { class: "ov-switch-label" }, label, hint ? h("small", null, hint) : null),
    input,
    h("span", { class: "ov-switch-track", "aria-hidden": "true" }),
  );
  el.set = (v) => { input.checked = !!v; };
  return el;
}

export function select({ label, options, value, onChange }) {
  const sel = h("select", { class: "ov-select", "aria-label": label });
  for (const o of options) {
    const opt = h("option", { value: o.value }, o.label);
    if (o.value === value) opt.selected = true;
    sel.append(opt);
  }
  sel.addEventListener("change", () => onChange?.(sel.value));
  const field = label
    ? h("div", { class: "ov-field" }, h("div", { class: "ov-field-head" }, h("span", { class: "ov-field-label" }, label)), sel)
    : sel;
  field.set = (v) => { sel.value = v; };
  return field;
}

export function miniSeg({ options, value, onChange, label, hint }) {
  const wrap = h("div", { class: "ov-mini-seg", role: "group", "aria-label": label || "" });
  const buttons = options.map((o) => {
    const b = h("button", { type: "button", "aria-pressed": String(o.value === value) }, o.label);
    b.addEventListener("click", () => {
      buttons.forEach((x) => x.setAttribute("aria-pressed", String(x === b)));
      onChange?.(o.value);
    });
    wrap.append(b);
    return b;
  });
  wrap.set = (v) => buttons.forEach((b, i) => b.setAttribute("aria-pressed", String(options[i].value === v)));
  return label
    ? Object.assign(h("div", { class: "ov-field" }, h("div", { class: "ov-field-head" }, h("span", { class: "ov-field-label" }, label)), wrap,
      hint ? h("small", { class: "ov-field-hint" }, hint) : null), { set: wrap.set })
    : wrap;
}

export function numberInput({ label, value, step = "any", onChange }) {
  const input = h("input", { class: "ov-input", type: "number", step, value: String(value), "aria-label": label });
  input.addEventListener("change", () => {
    const v = Number(input.value);
    if (Number.isFinite(v)) onChange?.(v);
  });
  input.addEventListener("keydown", (e) => e.stopPropagation());
  const field = h("div", { class: "ov-field" }, h("div", { class: "ov-field-head" }, h("span", { class: "ov-field-label" }, label)), input);
  field.set = (v) => { input.value = String(v); };
  return field;
}

const PALETTE = ["#f4c46a", "#ff7a59", "#ff4d6d", "#e56bff", "#9b7bff", "#4d8dff", "#2fb7ff", "#35d6c4", "#6fdc6b", "#c8e25a", "#ffffff", "#9aa4b2"];

export function colorPicker({ value, onChange, label = "Colour" }) {
  const wrap = h("div", { class: "ov-colors", role: "group", "aria-label": label });
  const chips = PALETTE.map((c) => {
    const b = h("button", { type: "button", class: "ov-chip", style: { background: c }, "aria-label": c, "aria-pressed": String(c === value) });
    b.addEventListener("click", () => { mark(c); onChange?.(c); });
    wrap.append(b);
    return b;
  });
  const custom = h("input", { type: "color", class: "ov-color-input", value: toHex6(value), "aria-label": "Custom colour" });
  custom.addEventListener("input", () => { mark(custom.value); onChange?.(custom.value); });
  wrap.append(custom);
  function mark(c) {
    chips.forEach((b, i) => b.setAttribute("aria-pressed", String(PALETTE[i].toLowerCase() === String(c).toLowerCase())));
  }
  wrap.set = (c) => { mark(c); custom.value = toHex6(c); };
  return wrap;
}

function toHex6(c) {
  const s = String(c || "#ffffff");
  return /^#[0-9a-f]{6}$/i.test(s) ? s : "#ffffff";
}

// ---- Tooltips ----------------------------------------------------------

export function installTooltips(root) {
  const tip = h("div", { class: "ov-tip", role: "tooltip" });
  root.append(tip);
  let current = null;
  let timer = 0;
  const show = (el) => {
    const text = el.getAttribute("data-tip");
    // No tip over a button whose menu is already open.
    if (!text || el.getAttribute("aria-expanded") === "true") return;
    const [main, keys] = text.split("  ");
    tip.textContent = main;
    if (keys) {
      for (const k of keys.split(" ")) tip.append(h("kbd", { class: "ov-kbd" }, k));
    }
    const r = el.getBoundingClientRect();
    const rr = root.getBoundingClientRect();
    tip.dataset.show = "true";
    const tw = tip.offsetWidth, th = tip.offsetHeight;
    let x = r.left + r.width / 2 - tw / 2;
    let y = r.bottom + 8;
    if (y + th > rr.bottom - 8) y = r.top - th - 8;
    x = Math.max(8, Math.min(x, rr.right - tw - 8));
    tip.style.left = `${x}px`;
    tip.style.top = `${y}px`;
  };
  const hide = () => {
    clearTimeout(timer);
    tip.dataset.show = "false";
    current = null;
  };
  root.addEventListener("pointerover", (e) => {
    const el = e.target.closest?.("[data-tip]");
    if (el === current) return;
    hide();
    if (!el || e.pointerType === "touch") return;
    current = el;
    timer = setTimeout(() => show(el), 280);
  });
  root.addEventListener("pointerdown", hide, true);
  root.addEventListener("pointerleave", hide);
  return { hide };
}

// ---- Toasts --------------------------------------------------------------

export function createToasts(root) {
  const host = h("div", { class: "ov-toasts", "aria-live": "polite" });
  root.append(host);
  return {
    show(text, { icon: iconEl = null, action, ms = 2600 } = {}) {
      const t = h("div", { class: "ov-toast ov-glass" }, iconEl, h("span", null, text));
      if (action) {
        const b = h("button", { type: "button" }, action.label);
        b.addEventListener("click", () => { action.run(); dismiss(); });
        t.append(b);
      }
      host.append(t);
      const dismiss = () => {
        t.dataset.leaving = "true";
        setTimeout(() => t.remove(), 240);
      };
      setTimeout(dismiss, ms);
      return dismiss;
    },
  };
}
