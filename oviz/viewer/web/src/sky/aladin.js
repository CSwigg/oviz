// Aladin Lite v3 integration: lazy loading, registered camera sync and the
// HiPS layer stack.
//
// Registration: the Sky camera sits at the Sun with galactic +z up. A
// pinhole camera is a gnomonic (TAN) projection, so the Aladin TAN view
// registers exactly when its centre (l, b), horizontal field of view and
// "galactic north up" orientation match the WebGL camera. We push the pose
// to Aladin synchronously from the render loop, every frame the camera moves.

const ALADIN_URLS = [
  "https://aladin.cds.unistra.fr/AladinLite/api/v3/latest/aladin.js",
  "https://aladin.u-strasbg.fr/AladinLite/api/v3/latest/aladin.js",
];

// IRSA serves HiPS without CORS headers; route through the CDS proxy.
const CORS_PROXY = "https://alaskybis.cds.unistra.fr/cgi/JSONProxy?url=";
const CORS_BLOCKED = ["https://irsa.ipac.caltech.edu/data/hips/", "http://irsa.ipac.caltech.edu/data/hips/"];
const SURVEY_ALIASES = new Map([
  ["P/GLIMPSE360", "https://irsa.ipac.caltech.edu/data/hips/Spitzer/GLIMPSE360"],
  ["IPAC/P/GLIMPSE360", "https://irsa.ipac.caltech.edu/data/hips/Spitzer/GLIMPSE360"],
  ["GLIMPSE360", "https://irsa.ipac.caltech.edu/data/hips/Spitzer/GLIMPSE360"],
]);

let loading = null;
let proxyInstalled = false;

function installCorsProxy() {
  if (proxyInstalled) return;
  proxyInstalled = true;
  const blocked = (url) => CORS_BLOCKED.some((p) => String(url || "").startsWith(p));
  const origFetch = window.fetch.bind(window);
  window.fetch = (input, init) => {
    const url = typeof input === "string" ? input : input?.url || "";
    if (blocked(url)) return origFetch(CORS_PROXY + encodeURIComponent(url), typeof input === "string" ? init : undefined);
    return origFetch(input, init);
  };
  const origOpen = XMLHttpRequest.prototype.open;
  XMLHttpRequest.prototype.open = function open(method, url, ...rest) {
    return origOpen.call(this, method, blocked(url) ? CORS_PROXY + encodeURIComponent(String(url)) : url, ...rest);
  };
}

export function loadAladin() {
  if (window.A?.aladin) return Promise.resolve(window.A);
  if (loading) return loading;
  installCorsProxy();
  loading = new Promise((resolve, reject) => {
    const tryUrl = (i) => {
      if (i >= ALADIN_URLS.length) {
        loading = null;
        reject(new Error("Aladin Lite could not be loaded (offline?)."));
        return;
      }
      const s = document.createElement("script");
      s.src = ALADIN_URLS[i];
      s.charset = "utf-8";
      s.async = true;
      s.onload = () => (window.A ? window.A.init.then(() => resolve(window.A), reject) : tryUrl(i + 1));
      s.onerror = () => { s.remove(); tryUrl(i + 1); };
      document.head.append(s);
    };
    tryUrl(0);
  });
  return loading;
}

export function resolveSurvey(id) {
  const key = String(id || "").trim();
  return SURVEY_ALIASES.get(key) || key;
}

export class AladinSky {
  constructor(host) {
    this.host = host;
    this.aladin = null;
    this.ready = false;
    this.base = null;
    this.overlays = new Map(); // key → {name, hips, opacity}
    this.order = [];
    this.lastSig = "";
  }

  async init({ survey = "P/DSS2/color", fov = 60 } = {}) {
    const A = await loadAladin();
    if (this.aladin) return this.aladin;
    const div = document.createElement("div");
    div.className = "ov-aladin";
    div.style.cssText = "position:absolute;inset:0;";
    this.host.append(div);
    this.aladin = A.aladin(div, {
      survey: resolveSurvey(survey),
      fov,
      target: "0 0",
      cooFrame: "galactic",
      projection: "TAN",
      showReticle: false,
      showZoomControl: false,
      showFullscreenControl: false,
      showLayersControl: false,
      showGotoControl: false,
      showShareControl: false,
      showSimbadPointerControl: false,
      showCooGridControl: false,
      showSettingsControl: false,
      showProjectionControl: false,
      showStatusBar: false,
      showCooLocation: false,
      showFrame: false,
      showFov: false,
      showContextMenu: false,
      showCooGrid: false,
      lockNorthUp: true,
      inertia: false,
      realFullscreen: false,
    });
    try { this.aladin.setProjection?.("TAN"); } catch (_) { /* older API */ }
    // Hook Aladin's animation loop: `redraw` re-requests `this.redrawClbk`
    // every frame, so wrapping that property lets Oviz render in the same
    // callback, immediately before Aladin draws.
    const view = this.aladin.view;
    const loop = view && view.redrawClbk;
    if (typeof loop === "function") {
      view.redrawClbk = (t) => {
        try { this.onBeforeDraw?.(t); } catch (err) { console.error(err); }
        return loop(t);
      };
      this.synchronized = true;
    }
    this.base = resolveSurvey(survey);
    this.ready = true;
    this.host.dataset.ready = "true";
    return this.aladin;
  }

  /** Point Aladin at galactic (l, b) with a horizontal field of view. */
  setView(lDeg, bDeg, hfovDeg) {
    const a = this.aladin;
    if (!a) return;
    const sig = `${lDeg.toFixed(6)}|${bDeg.toFixed(6)}|${hfovDeg.toFixed(5)}`;
    if (sig === this.lastSig) return;
    try {
      a.gotoPosition(((lDeg % 360) + 360) % 360, bDeg);
      a.setFoV(Math.min(Math.max(hfovDeg, 1e-4), 179));
      a.view?.requestRedraw?.();
      // Only now is the pose Aladin's: a failed push is retried next frame.
      this.lastSig = sig;
    } catch (err) {
      // Aladin throws while its WebGL context is still warming up.
    }
  }

  _makeHiPS(id, opts) {
    const A = window.A;
    const a = this.aladin;
    const url = resolveSurvey(id);
    if (typeof a.newImageSurvey === "function") return a.newImageSurvey(url, { name: id, ...opts });
    if (typeof A.imageHiPS === "function") return A.imageHiPS(url, { name: id, ...opts });
    if (typeof A.HiPS === "function") return A.HiPS(url, { name: id, ...opts });
    return url;
  }

  /**
   * Apply a layer stack (top first). The bottom visible layer is the base;
   * the rest are overlays with their own opacity. Attach/detach is
   * incremental so visible layers never blink during State transitions.
   */
  applyStack(layers, { backgroundVisible = true } = {}) {
    const a = this.aladin;
    if (!a) return;
    const visible = layers.filter((l) => l.visible !== false && (l.opacity ?? 1) > 0.001 && backgroundVisible);
    const bottom = visible[visible.length - 1];
    const baseId = bottom ? resolveSurvey(bottom.survey || bottom.key) : null;
    if (baseId && baseId !== this.base) {
      try {
        (a.setImageSurvey || a.setBaseImageLayer).call(a, baseId);
        this.base = baseId;
      } catch (err) {
        console.warn("Oviz sky: base survey failed", baseId, err);
      }
    }
    try {
      const baseLayer = a.getBaseImageLayer?.();
      baseLayer?.setOpacity?.(bottom ? Math.min(1, bottom.opacity ?? 1) : 0);
    } catch (_) { /* ignore */ }
    const overlays = visible.slice(0, -1).reverse(); // bottom → top
    const wantKeys = overlays.map((l) => l.key);
    // Aladin always adds an overlay on top of the stack. Keep the attached
    // layers that already sit in the wanted order, and detach everything
    // above the first difference so it is re-added in order below.
    let keep = 0;
    while (keep < this.order.length && keep < wantKeys.length && this.order[keep] === wantKeys[keep]) keep++;
    for (const key of this.order.slice(keep)) {
      const o = this.overlays.get(key);
      if (o) this._detach(key, o);
    }
    for (const [key, o] of this.overlays) {
      if (!wantKeys.includes(key)) this._detach(key, o);
    }
    const attached = this.order.slice(0, keep);
    for (const l of overlays) {
      let o = this.overlays.get(l.key);
      const opacity = Math.min(1, Math.max(0, l.opacity ?? 1));
      if (!o) {
        const name = `ov-${l.key}`;
        const hips = this._makeHiPS(l.survey || l.key, {});
        try { hips?.setOpacity?.(opacity); } catch (_) { /* not ready */ }
        try {
          if (typeof a.addImageLayer === "function") a.addImageLayer(hips, name);
          else a.setOverlayImageLayer(hips, name);
        } catch (err) {
          console.warn("Oviz sky: overlay failed", l.key, err);
          continue;
        }
        o = { name, hips, opacity: -1 };
        this.overlays.set(l.key, o);
        attached.push(l.key);
      }
      if (Math.abs(o.opacity - opacity) > 1e-4) {
        o.opacity = opacity;
        this._setOpacity(o, opacity);
      }
      this._applyStyle(o, l);
    }
    // The real stack (an overlay that failed to attach is not in it).
    this.order = attached;
  }

  _setOpacity(o, opacity) {
    const a = this.aladin;
    let layer = null;
    try { layer = a.getOverlayImageLayer?.(o.name); } catch (_) { layer = null; }
    const target = layer || o.hips;
    try { target?.setOpacity?.(opacity); } catch (_) { /* retried next apply */ }
    if (!layer) setTimeout(() => {
      try { a.getOverlayImageLayer?.(o.name)?.setOpacity?.(opacity); } catch (_) { /* ignore */ }
    }, 250);
  }

  _applyStyle(o, l) {
    const a = this.aladin;
    let layer = null;
    try { layer = a.getOverlayImageLayer?.(o.name); } catch (_) { return; }
    if (!layer) return;
    const sig = `${l.stretch || ""}|${l.colormap || ""}|${l.cut_min ?? ""}|${l.cut_max ?? ""}`;
    if (o.styleSig === sig) return;
    o.styleSig = sig;
    try {
      if (l.colormap) layer.setColormap?.(l.colormap, { stretch: l.stretch || "linear" });
      else if (l.stretch) layer.setColormap?.("native", { stretch: l.stretch });
      if (Number.isFinite(l.cut_min) && Number.isFinite(l.cut_max) && l.cut_max > l.cut_min) layer.setCuts?.(l.cut_min, l.cut_max);
    } catch (_) { /* optional features */ }
  }

  _detach(key, o) {
    try { this.aladin.removeImageLayer?.(o.name); } catch (_) { /* ignore */ }
    this.overlays.delete(key);
  }

  /** Draw the current Aladin view into a 2D canvas (for screenshots). */
  async drawInto(ctx, w, h) {
    const a = this.aladin;
    if (!a) return;
    try {
      const url = await a.getViewDataURL?.({ format: "image/png", width: w, height: h });
      if (!url) return;
      const img = new Image();
      img.src = url;
      await img.decode();
      ctx.drawImage(img, 0, 0, w, h);
    } catch (err) {
      console.warn("Oviz sky: background capture failed", err);
    }
  }
}
