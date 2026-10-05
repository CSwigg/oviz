// Tiny DOM helpers: hyperscript, icons, formatting, focus utilities.

export function h(tag, props, ...children) {
  const el = document.createElement(tag);
  if (props) {
    for (const [k, v] of Object.entries(props)) {
      if (v == null || v === false) continue;
      if (k === "class") el.className = v;
      else if (k === "style" && typeof v === "object") Object.assign(el.style, v);
      else if (k === "dataset") Object.assign(el.dataset, v);
      else if (k.startsWith("on") && typeof v === "function") el.addEventListener(k.slice(2).toLowerCase(), v);
      else if (k === "html") el.innerHTML = v;
      else if (v === true) el.setAttribute(k, "");
      else el.setAttribute(k, String(v));
    }
  }
  append(el, children);
  return el;
}

function append(el, children) {
  for (const c of children) {
    if (c == null || c === false) continue;
    if (Array.isArray(c)) append(el, c);
    else el.append(c instanceof Node ? c : document.createTextNode(String(c)));
  }
}

export function clear(el) {
  while (el.firstChild) el.firstChild.remove();
  return el;
}

// 20×20 stroke icons (1.6 px), designed as one family.
const ICONS = {
  play: '<path d="M7 4.8v10.4a.6.6 0 0 0 .92.5l8.1-5.2a.6.6 0 0 0 0-1L7.92 4.3A.6.6 0 0 0 7 4.8Z" fill="currentColor" stroke="none"/>',
  pause: '<rect x="5.5" y="4.5" width="3" height="11" rx="1" fill="currentColor" stroke="none"/><rect x="11.5" y="4.5" width="3" height="11" rx="1" fill="currentColor" stroke="none"/>',
  stepBack: '<path d="M5.5 5v10"/><path d="M14.5 5.4 8 10l6.5 4.6V5.4Z" fill="currentColor" stroke="none"/>',
  stepFwd: '<path d="M14.5 5v10"/><path d="M5.5 5.4 12 10l-6.5 4.6V5.4Z" fill="currentColor" stroke="none"/>',
  loop: '<path d="M4 9.5V8a3 3 0 0 1 3-3h8.5"/><path d="m13 2.8 2.5 2.2L13 7.2"/><path d="M16 10.5V12a3 3 0 0 1-3 3H4.5"/><path d="m7 17.2-2.5-2.2L7 12.8"/>',
  search: '<circle cx="9" cy="9" r="5.2"/><path d="m13 13 3.5 3.5"/>',
  layers: '<path d="m10 3 7 3.8-7 3.8-7-3.8L10 3Z"/><path d="m3 10.2 7 3.8 7-3.8"/><path d="m3 13.6 7 3.8 7-3.8"/>',
  bookmark: '<path d="M6 3.5h8a1 1 0 0 1 1 1V17l-5-3.2L5 17V4.5a1 1 0 0 1 1-1Z"/>',
  camera: '<path d="M3.5 7a1.5 1.5 0 0 1 1.5-1.5h1.6l1.2-1.8h4.4l1.2 1.8H15A1.5 1.5 0 0 1 16.5 7v7.5A1.5 1.5 0 0 1 15 16H5a1.5 1.5 0 0 1-1.5-1.5V7Z"/><circle cx="10" cy="10.6" r="2.8"/>',
  download: '<path d="M10 3.5v9"/><path d="m6.5 9.2 3.5 3.5 3.5-3.5"/><path d="M4 15.5h12"/>',
  share: '<path d="M8.5 11.5 11.5 8.5"/><path d="M9.2 6.3 10.6 5a3 3 0 0 1 4.3 4.3l-1.4 1.4"/><path d="M10.8 13.7 9.4 15a3 3 0 0 1-4.3-4.3l1.4-1.4"/>',
  sliders: '<path d="M4 6h7"/><path d="M15 6h1"/><circle cx="13" cy="6" r="1.8"/><path d="M4 14h1"/><path d="M9 14h7"/><circle cx="7" cy="14" r="1.8"/>',
  help: '<circle cx="10" cy="10" r="6.8"/><path d="M8.2 8.1a1.9 1.9 0 1 1 2.6 1.8c-.5.2-.8.6-.8 1.1v.6"/><circle cx="10" cy="13.9" r=".5" fill="currentColor"/>',
  expand: '<path d="M4 8V4h4"/><path d="M16 8V4h-4"/><path d="M4 12v4h4"/><path d="M16 12v4h-4"/>',
  collapse: '<path d="M8 4v4H4"/><path d="M12 4v4h4"/><path d="M8 16v-4H4"/><path d="M12 16v-4h4"/>',
  eye: '<path d="M2.8 10s2.6-5 7.2-5 7.2 5 7.2 5-2.6 5-7.2 5-7.2-5-7.2-5Z"/><circle cx="10" cy="10" r="2.2"/>',
  eyeOff: '<path d="M4.5 4.5l11 11"/><path d="M8.2 5.3A7 7 0 0 1 10 5c4.6 0 7.2 5 7.2 5a12 12 0 0 1-2 2.6"/><path d="M12.6 14.4A6.8 6.8 0 0 1 10 15c-4.6 0-7.2-5-7.2-5a12.3 12.3 0 0 1 2.5-3"/>',
  chevron: '<path d="m7.5 5 5 5-5 5"/>',
  chevronDown: '<path d="m5 7.5 5 5 5-5"/>',
  close: '<path d="m5.5 5.5 9 9"/><path d="m14.5 5.5-9 9"/>',
  plus: '<path d="M10 4.5v11"/><path d="M4.5 10h11"/>',
  target: '<circle cx="10" cy="10" r="5.8"/><circle cx="10" cy="10" r="1.6" fill="currentColor"/><path d="M10 1.8v2.4M10 15.8v2.4M1.8 10h2.4M15.8 10h2.4"/>',
  home: '<path d="M3.8 9.2 10 4l6.2 5.2"/><path d="M5.5 8v7.5h9V8"/>',
  globe: '<circle cx="10" cy="10" r="6.8"/><path d="M3.2 10h13.6"/><path d="M10 3.2c2 2 2.8 4.3 2.8 6.8S12 14.8 10 16.8C8 14.8 7.2 12.5 7.2 10S8 5.2 10 3.2Z"/>',
  cube: '<path d="m10 2.8 6.2 3.5v7.4L10 17.2l-6.2-3.5V6.3L10 2.8Z"/><path d="m3.9 6.3 6.1 3.5 6.1-3.5"/><path d="M10 9.8v7.3"/>',
  ar: '<path d="M3 7V4.5A1.5 1.5 0 0 1 4.5 3H7M13 3h2.5A1.5 1.5 0 0 1 17 4.5V7M17 13v2.5a1.5 1.5 0 0 1-1.5 1.5H13M7 17H4.5A1.5 1.5 0 0 1 3 15.5V13"/><path d="m10 5.8 3.7 2.1v4.2L10 14.2l-3.7-2.1V7.9L10 5.8Z"/><path d="M6.4 7.9 10 10l3.6-2.1M10 10v4.2"/>',
  sparkle: '<path d="M10 3v3.2M10 13.8V17M3 10h3.2M13.8 10H17"/><path d="m5.1 5.1 2 2m5.8 5.8 2 2m0-9.8-2 2m-5.8 5.8-2 2"/>',
  orbit: '<ellipse cx="10" cy="10" rx="7" ry="3.2" transform="rotate(-25 10 10)"/><circle cx="10" cy="10" r="1.8" fill="currentColor"/>',
  grid: '<circle cx="10" cy="10" r="6.8"/><circle cx="10" cy="10" r="3.4"/><path d="M10 3.2v13.6M3.2 10h13.6"/>',
  present: '<rect x="3" y="4" width="14" height="9.5" rx="1.2"/><path d="M10 13.5V17M7 17h6"/><path d="m8.6 7 3 1.8-3 1.8V7Z" fill="currentColor" stroke="none"/>',
  record: '<circle cx="10" cy="10" r="6.5"/><circle cx="10" cy="10" r="3" fill="currentColor" stroke="none"/>',
  stop: '<rect x="5.5" y="5.5" width="9" height="9" rx="1.6" fill="currentColor" stroke="none"/>',
  ruler: '<path d="m3.5 13.3 9.8-9.8 3.2 3.2-9.8 9.8Z"/><path d="m6.2 10.6 1.4 1.4M8.4 8.4l1 1M10.6 6.2l1.4 1.4"/>',
  follow: '<circle cx="10" cy="10" r="2.2"/><path d="M10 3.5a6.5 6.5 0 0 1 6.5 6.5"/><path d="M10 16.5A6.5 6.5 0 0 1 3.5 10"/><path d="m15 7.6 1.5 2.4 2-1.9M5 12.4 3.5 10l-2 1.9"/>',
  trail: '<path d="M3.5 15.5c2.2-5.6 5.2-8.4 9-9"/><circle cx="14.8" cy="6.2" r="2" fill="currentColor"/><circle cx="6.2" cy="11.7" r=".7" fill="currentColor"/><circle cx="9" cy="8.6" r=".9" fill="currentColor"/>',
  copy: '<rect x="6.5" y="6.5" width="9" height="9" rx="1.5"/><path d="M13.5 4.5v-.2A.8.8 0 0 0 12.7 3.5H4.3a.8.8 0 0 0-.8.8v8.4a.8.8 0 0 0 .8.8h.2"/>',
  edit: '<path d="M4 16h3.2l8.3-8.3a1.6 1.6 0 0 0-2.2-2.2L5 13.8V16Z"/>',
  trash: '<path d="M4.5 6h11"/><path d="M8 6V4.5h4V6"/><path d="M6 6l.7 9.5h6.6L14 6"/>',
  drag: '<circle cx="8" cy="6" r=".9" fill="currentColor"/><circle cx="12" cy="6" r=".9" fill="currentColor"/><circle cx="8" cy="10" r=".9" fill="currentColor"/><circle cx="12" cy="10" r=".9" fill="currentColor"/><circle cx="8" cy="14" r=".9" fill="currentColor"/><circle cx="12" cy="14" r=".9" fill="currentColor"/>',
  info: '<circle cx="10" cy="10" r="6.8"/><path d="M10 9v4.6"/><circle cx="10" cy="6.4" r=".6" fill="currentColor"/>',
  keyboard: '<rect x="2.8" y="5.2" width="14.4" height="9.6" rx="1.5"/><path d="M5.5 8h.01M8.2 8h.01M11 8h.01M13.8 8h.01M5.5 10.6h.01M13.8 10.6h.01M7.6 12.2h4.8"/>',
  moon: '<path d="M15.8 12.3A6.5 6.5 0 0 1 7.7 4.2a6.5 6.5 0 1 0 8.1 8.1Z"/>',
  sun: '<circle cx="10" cy="10" r="3.2"/><path d="M10 2.8v1.6M10 15.6v1.6M2.8 10h1.6M15.6 10h1.6M4.9 4.9l1.1 1.1M14 14l1.1 1.1M15.1 4.9 14 6M6 14l-1.1 1.1"/>',
  stars: '<path d="m8 3.5 1.1 3.1 3.1 1.1-3.1 1.1L8 11.9 6.9 8.8 3.8 7.7l3.1-1.1L8 3.5Z"/><path d="m14.2 10.8.6 1.6 1.6.6-1.6.6-.6 1.6-.6-1.6-1.6-.6 1.6-.6.6-1.6Z"/>',
  wand: '<path d="m4 16 8.5-8.5"/><path d="m11 6 1.5 1.5"/><path d="M14.5 2.5v2M13.5 3.5h2M16.5 7v1.4M15.8 7.7h1.4"/>',
  more: '<circle cx="5" cy="10" r="1.4" fill="currentColor" stroke="none"/><circle cx="10" cy="10" r="1.4" fill="currentColor" stroke="none"/><circle cx="15" cy="10" r="1.4" fill="currentColor" stroke="none"/>',
  minus: '<path d="M4.5 10h11"/>',
  menu: '<path d="M4 6h12"/><path d="M4 10h12"/><path d="M4 14h12"/>',
  widgets: '<rect x="3" y="3" width="6" height="6" rx="1.5"/><rect x="11" y="3" width="6" height="6" rx="1.5"/><rect x="3" y="11" width="6" height="6" rx="1.5"/><rect x="11" y="11" width="6" height="6" rx="1.5"/>',
  sidebar: '<rect x="3.5" y="4.5" width="13" height="11" rx="2"/><path d="M8 4.5v11"/>',
  chart: '<path d="M3.5 16.5h13"/><path d="m3.8 13.6 3.4-4.1 3 2.6 3.4-5.6 2.9 3.6"/>',
};

export function icon(name, cls = "") {
  const span = document.createElement("span");
  span.className = `ov-icon ${cls}`;
  span.setAttribute("aria-hidden", "true");
  span.innerHTML = `<svg viewBox="0 0 20 20" width="20" height="20" fill="none" stroke="currentColor" stroke-width="1.6" stroke-linecap="round" stroke-linejoin="round">${ICONS[name] || ""}</svg>`;
  return span;
}

export function iconButton(name, label, onClick, { cls = "", shortcut = "", pressed } = {}) {
  const b = h("button", {
    class: `ov-ibtn ${cls}`,
    type: "button",
    "aria-label": label,
    "data-tip": shortcut ? `${label}  ${shortcut}` : label,
    onclick: onClick,
  }, icon(name));
  if (pressed !== undefined) b.setAttribute("aria-pressed", String(!!pressed));
  return b;
}

export function fmt(n, digits = 1) {
  if (n == null || !Number.isFinite(n)) return "—";
  const s = Math.abs(n) >= 1e4 ? Math.round(n).toLocaleString() : n.toFixed(digits);
  return s.replace("-", "−");
}

export function fmtInt(n) {
  if (n == null || !Number.isFinite(n)) return "—";
  return Math.round(n).toLocaleString();
}

export function kbd(text) {
  return h("kbd", { class: "ov-kbd" }, text);
}

const isMac = typeof navigator !== "undefined" && /Mac|iPhone|iPad/.test(navigator.platform || navigator.userAgent);
export const MOD = isMac ? "⌘" : "Ctrl";

/** Is the event target a text-editing element (so shortcuts should yield)? */
export function isEditable(target) {
  if (!target || !(target instanceof Element)) return false;
  if (target.isContentEditable) return true;
  const tag = target.tagName;
  if (tag === "TEXTAREA" || tag === "SELECT") return true;
  if (tag === "INPUT") {
    const type = (target.getAttribute("type") || "text").toLowerCase();
    return !["range", "checkbox", "radio", "button", "color"].includes(type);
  }
  return false;
}

/**
 * Does the focused control handle this key itself? Sliders own arrows,
 * Home/End and paging; buttons, switches, options and links own Space and
 * Enter; selects own navigation keys. Global shortcuts must yield to them.
 */
export function controlConsumesKey(e) {
  const t = e.target;
  if (!t || !(t instanceof Element)) return false;
  const k = e.key;
  const role = t.getAttribute("role") || "";
  const type = t.tagName === "INPUT" ? (t.getAttribute("type") || "text").toLowerCase() : "";
  if (type === "range" || role === "slider") {
    return ["ArrowLeft", "ArrowRight", "ArrowUp", "ArrowDown", "Home", "End", "PageUp", "PageDown"].includes(k);
  }
  if (t.tagName === "SELECT") return ["ArrowUp", "ArrowDown", " ", "Enter"].includes(k);
  const activatable = t.tagName === "BUTTON" || t.tagName === "A" || t.tagName === "SUMMARY"
    || type === "checkbox" || type === "radio" || type === "color"
    || ["button", "switch", "option", "menuitem", "tab", "listitem"].includes(role);
  if (activatable) return k === " " || k === "Enter";
  return false;
}

/**
 * Fuzzy score: substring matches rank highest (earlier and word-start
 * better); otherwise a subsequence match must be compact — its span may be
 * at most ~2.5× the query length — so scattered letters never match.
 */
export function fuzzyScore(query, text) {
  if (!query) return 1;
  const q = query.toLowerCase(), t = text.toLowerCase();
  const idx = t.indexOf(q);
  if (idx >= 0) return 100 - idx * 0.5 + (idx === 0 || /[\s_\-(]/.test(t[idx - 1]) ? 40 : 0) - t.length * 0.02;
  if (q.length < 2) return 0;
  let best = 0;
  // Try each possible start of the first character for the tightest span.
  for (let start = t.indexOf(q[0]); start >= 0; start = t.indexOf(q[0], start + 1)) {
    let ti = start + 1, score = 3, streak = 0, ok = true;
    for (let k = 1; k < q.length; k++) {
      const j = t.indexOf(q[k], ti);
      if (j < 0) { ok = false; break; }
      streak = j === ti ? streak + 1 : 0;
      score += 1 + streak * 2 + (/[\s_\-(]/.test(t[j - 1]) ? 3 : 0);
      ti = j + 1;
    }
    if (!ok) break;
    const span = ti - start;
    if (span > Math.max(q.length * 2.5, q.length + 3)) continue;
    const s = score - span * 0.3 + (start === 0 || /[\s_\-(]/.test(t[start - 1]) ? 4 : 0);
    if (s > best) best = s;
  }
  return best > 0 ? best - t.length * 0.02 : 0;
}

export function copyText(text) {
  if (navigator.clipboard?.writeText) return navigator.clipboard.writeText(text);
  const ta = h("textarea", { style: { position: "fixed", opacity: "0" } }, text);
  document.body.append(ta);
  ta.select();
  document.execCommand("copy");
  ta.remove();
  return Promise.resolve();
}

/**
 * Copy text that is still being made (a promise) without losing the click:
 * Safari only lets the page write to the clipboard during the gesture, so
 * the clipboard item is handed over now and filled in when the text is ready.
 */
export function copyTextLater(textPromise) {
  if (typeof ClipboardItem === "function" && navigator.clipboard?.write) {
    try {
      const item = new ClipboardItem({ "text/plain": textPromise.then((t) => new Blob([t], { type: "text/plain" })) });
      return navigator.clipboard.write([item]).catch(() => textPromise.then(copyText));
    } catch (_) { /* older ClipboardItem: copy when ready */ }
  }
  return textPromise.then(copyText);
}

export function downloadBlob(blob, filename) {
  const url = URL.createObjectURL(blob);
  const a = h("a", { href: url, download: filename });
  document.body.append(a);
  a.click();
  a.remove();
  setTimeout(() => URL.revokeObjectURL(url), 4000);
}

export function localStorageGet(k) {
  try { return localStorage.getItem(k); } catch (_) { return null; }
}

export function localStorageSet(k, v) {
  try { localStorage.setItem(k, v); } catch (_) { /* private mode or quota */ }
}
