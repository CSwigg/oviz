// Fluid motion: physically based springs, shared-element reveals (a menu
// grows out of the button that opened it), cascading items, pointer
// parallax, a gaze highlight, and drag-to-dismiss with flick physics.
// Everything collapses to instant changes under prefers-reduced-motion.

const reduceQuery = typeof window !== "undefined" ? window.matchMedia?.("(prefers-reduced-motion: reduce)") : null;
export const reducedMotion = () => !!reduceQuery?.matches;

let linearSupported = null;
function supportsLinear() {
  if (linearSupported == null) {
    try { linearSupported = CSS.supports("animation-timing-function", "linear(0, 1)"); } catch (_) { linearSupported = false; }
  }
  return linearSupported;
}

const cache = new Map();

/**
 * Simulate a unit spring (0 → 1) and sample it into a CSS `linear()`
 * easing, so springs run on the compositor via the Web Animations API.
 */
export function springEasing({ stiffness = 320, damping = 26, mass = 1 } = {}) {
  const key = `${stiffness}|${damping}|${mass}`;
  if (cache.has(key)) return cache.get(key);
  const dt = 1 / 240;
  let x = 0, v = 0, t = 0;
  const pts = [0];
  while (t < 2.5) {
    const a = (stiffness * (1 - x) - damping * v) / mass;
    v += a * dt;
    x += v * dt;
    t += dt;
    pts.push(x);
    if (t > 0.08 && Math.abs(1 - x) < 0.0015 && Math.abs(v) < 0.02) break;
  }
  const n = Math.min(pts.length - 1, 60);
  const out = [];
  for (let i = 0; i <= n; i++) out.push(+pts[Math.round((i / n) * (pts.length - 1))].toFixed(4));
  out[0] = 0;
  out[out.length - 1] = 1;
  const res = {
    easing: supportsLinear() ? `linear(${out.join(", ")})` : "cubic-bezier(0.2, 0.9, 0.25, 1.08)",
    duration: Math.round(t * 1000),
  };
  cache.set(key, res);
  return res;
}

export const SPRINGS = {
  pop: { stiffness: 380, damping: 26 },     // menus, palette
  soft: { stiffness: 220, damping: 24 },    // sheets, cards
  snappy: { stiffness: 520, damping: 34 },  // closing
  wobble: { stiffness: 300, damping: 16 },  // playful settles
};

function rectIn(el) {
  const r = el.getBoundingClientRect();
  return { left: r.left, top: r.top, width: r.width, height: r.height, right: r.right, bottom: r.bottom };
}

/**
 * Keyframes that start as a pill sitting exactly on `from` (a DOMRect of
 * the trigger) and end as `el` in place: the element appears to grow out of
 * the button, with no text distortion (clip + translate, never scale).
 */
function growFrames(el, from) {
  // Compose with the element's own transform (centred dialogs use one).
  el.style.animation = "none";
  const cs = getComputedStyle(el);
  const base = cs.transform === "none" ? "" : cs.transform;
  const at = (t) => (t || base ? `${t} ${base}`.trim() : "none");
  const r = rectIn(el);
  const radius = cs.borderRadius || "16px";
  if (!from || !from.width) {
    const w = Math.min(120, r.width), hh = Math.min(40, r.height);
    const cx = (r.width - w) / 2, cy = (r.height - hh) / 2;
    return [
      { clipPath: `inset(${cy}px ${cx}px ${cy}px ${cx}px round ${hh / 2}px)`, transform: at("translateY(10px) scale(0.98)"), opacity: 0 },
      { clipPath: `inset(0px 0px 0px 0px round ${radius})`, transform: at(""), opacity: 1 },
    ];
  }
  const below = from.top < r.top + r.height / 2;
  const w = Math.min(from.width, r.width);
  const hh = Math.min(from.height, r.height);
  const left = Math.max(0, Math.min(from.left - r.left, r.width - w));
  const top = below ? 0 : r.height - hh;
  const dx = from.left - (r.left + left);
  const dy = from.top - (r.top + top);
  const pill = Math.min(w, hh) / 2;
  return [
    {
      clipPath: `inset(${top}px ${r.width - left - w}px ${r.height - top - hh}px ${left}px round ${pill}px)`,
      transform: at(`translate(${dx}px, ${dy}px)`),
      opacity: 0.6,
    },
    { clipPath: `inset(0px 0px 0px 0px round ${radius})`, transform: at(""), opacity: 1 },
  ];
}

/** Grow `el` out of the trigger rect `from`, then cascade its items in. */
export function reveal(el, from, { spring = SPRINGS.pop, items = null } = {}) {
  if (reducedMotion() || !el.animate) return;
  const s = springEasing(spring);
  el.animate(growFrames(el, from), { duration: s.duration, easing: s.easing });
  if (items) cascade(el.querySelectorAll(items));
}

/** Shrink `el` back into `from`; resolves when it is safe to remove. */
export function conceal(el, from) {
  if (reducedMotion() || !el.animate || !el.isConnected) return Promise.resolve();
  const s = springEasing(SPRINGS.snappy);
  const frames = growFrames(el, from).reverse();
  frames[1].opacity = 0;
  const anim = el.animate(frames, { duration: Math.min(s.duration, 320), easing: "cubic-bezier(0.4, 0, 0.6, 1)", fill: "forwards" });
  el.style.pointerEvents = "none";
  return anim.finished.catch(() => {});
}

/** Items fall into place one after another. */
export function cascade(nodes, { step = 22, max = 14, from = "translateY(-6px) scale(0.98)" } = {}) {
  if (reducedMotion()) return;
  const s = springEasing(SPRINGS.pop);
  let i = 0;
  for (const n of nodes) {
    if (i >= max) break;
    n.animate?.([{ opacity: 0, transform: from }, { opacity: 1, transform: "none" }],
      { duration: s.duration, easing: s.easing, delay: 40 + i * step, fill: "backwards" });
    i++;
  }
}

/**
 * Pointer parallax: smoothed pointer position in [-1, 1] published as
 * --ov-px / --ov-py on the root, for skins that float chrome in depth.
 */
export function installParallax(root) {
  let tx = 0, ty = 0, x = 0, y = 0, raf = 0;
  const step = () => {
    x += (tx - x) * 0.12;
    y += (ty - y) * 0.12;
    root.style.setProperty("--ov-px", x.toFixed(4));
    root.style.setProperty("--ov-py", y.toFixed(4));
    raf = Math.abs(tx - x) + Math.abs(ty - y) > 0.001 ? requestAnimationFrame(step) : 0;
  };
  root.addEventListener("pointermove", (e) => {
    if (e.pointerType !== "mouse" || reducedMotion()) return;
    const r = root.getBoundingClientRect();
    tx = ((e.clientX - r.left) / r.width) * 2 - 1;
    ty = ((e.clientY - r.top) / r.height) * 2 - 1;
    if (!raf) raf = requestAnimationFrame(step);
  }, { passive: true });
}

/** Gaze highlight: interactive elements light up where the pointer is. */
const GAZE = ".ov-ibtn, .ov-btn, .ov-row, .ov-pop-item, .ov-palette-item, .ov-card, .ov-skin-card, .ov-seg button, .ov-search-trigger";
export function installGaze(root) {
  let last = null;
  root.addEventListener("pointermove", (e) => {
    if (e.pointerType !== "mouse") return;
    const el = e.target.closest?.(GAZE);
    if (last && last !== el) last.removeAttribute("data-gaze");
    last = el;
    if (!el) return;
    const r = el.getBoundingClientRect();
    el.style.setProperty("--gx", `${(e.clientX - r.left).toFixed(0)}px`);
    el.style.setProperty("--gy", `${(e.clientY - r.top).toFixed(0)}px`);
    el.dataset.gaze = "true";
  }, { passive: true });
  root.addEventListener("pointerleave", () => { last?.removeAttribute("data-gaze"); last = null; });
}

/**
 * Drag a card along `axis` ("x" or "y") with the pointer; release past a
 * distance or with a flick to throw it away (onDismiss), otherwise it
 * springs home. Only starts from `handle` (and never from controls).
 */
export function flickToDismiss(el, handle, { axis = "x", sign = 1, onDismiss, enabled = () => true } = {}) {
  let start = null;
  const pos = (e) => (axis === "x" ? e.clientX : e.clientY);
  handle.addEventListener("pointerdown", (e) => {
    if (e.button !== 0 || !enabled() || e.target.closest("button, input, select, a, canvas, .ov-range")) return;
    start = { p: pos(e), t: performance.now(), last: pos(e), lt: performance.now(), v: 0 };
    handle.setPointerCapture(e.pointerId);
    el.dataset.dragging = "true";
  });
  handle.addEventListener("pointermove", (e) => {
    if (!start) return;
    const now = performance.now();
    const p = pos(e);
    start.v = (p - start.last) / Math.max(1, now - start.lt);
    start.last = p;
    start.lt = now;
    let d = (p - start.p) * sign;
    if (d < 0) d = -Math.sqrt(-d) * 3; // rubber band the wrong way
    el.style.translate = axis === "x" ? `${d * sign}px 0` : `0 ${d * sign}px`;
    el.style.rotate = axis === "x" ? `${(d * sign) / 60}deg` : "";
    el.style.opacity = String(Math.max(0.35, 1 - Math.max(0, d) / 600));
  });
  const end = () => {
    if (!start) return;
    const d = (start.last - start.p) * sign;
    const v = start.v * sign;
    start = null;
    el.dataset.dragging = "false";
    const s = springEasing(SPRINGS.soft);
    if (d > 110 || v > 0.9) {
      const to = axis === "x" ? `${sign * (innerWidth * 0.6)}px 0` : `0 ${sign * innerHeight * 0.6}px`;
      const anim = el.animate?.([{ translate: el.style.translate || "0 0", opacity: el.style.opacity || 1 }, { translate: to, opacity: 0 }],
        { duration: 260, easing: "cubic-bezier(0.3, 0.6, 0.4, 1)" });
      const done = () => { el.style.translate = ""; el.style.rotate = ""; el.style.opacity = ""; onDismiss?.(); };
      if (anim) anim.finished.then(done, done); else done();
    } else {
      el.animate?.([{ translate: el.style.translate || "0 0", rotate: el.style.rotate || "0deg" }, { translate: "0 0", rotate: "0deg" }],
        { duration: s.duration, easing: s.easing });
      el.style.translate = "";
      el.style.rotate = "";
      el.style.opacity = "";
    }
  };
  handle.addEventListener("pointerup", end);
  handle.addEventListener("pointercancel", end);
}
