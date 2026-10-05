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

// Hidden overlays kept attached so showing them again is instant.
const MAX_HIDDEN = 6;

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
  if (loading) return loading;
  // The script can be on the page before its WebAssembly core is ready:
  // always wait for A.init, or Aladin starts half-initialised.
  if (window.A?.aladin) {
    loading = Promise.resolve(window.A.init).then(() => window.A);
    loading.catch(() => { loading = null; });
    return loading;
  }
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
    this.overlays = new Map(); // key → {name, hips, opacity, used}
    this.order = []; // attached overlays, bottom → top
    this.failed = new Map(); // layer key → load error
    this.tries = new Map(); // layer key → failed attempts so far
    this.lastSig = "";
  }

  /** Start Aladin once, however many callers ask at the same time. */
  init(opts = {}) {
    if (!this._init) {
      this._init = this._start(opts);
      this._init.catch(() => { this._init = null; });
    }
    return this._init;
  }

  async _start({ survey = "P/DSS2/color", fov = 60 } = {}) {
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
    if (typeof A.imageHiPS === "function") return A.imageHiPS(url, { name: id, ...opts });
    if (typeof a.newImageSurvey === "function") return a.newImageSurvey(url, { name: id, ...opts });
    if (typeof A.HiPS === "function") return A.HiPS(url, { name: id, ...opts });
    return url;
  }

  /**
   * Apply a layer stack (top first). Aladin's base layer keeps the survey it
   * opened with and stands in for that survey whenever it is the bottom
   * visible layer; every other survey is an overlay. Overlays stay attached
   * (at zero opacity) once used, so switching, blending and State
   * crossfades change opacities instead of re-downloading tiles. Only an
   * order change re-attaches the layers above it.
   */
  applyStack(layers, { backgroundVisible = true } = {}) {
    const a = this.aladin;
    if (!a) return;
    const shown = (l) => backgroundVisible && l.visible !== false && (l.opacity ?? 1) > 0.001;
    const visible = layers.filter(shown).reverse(); // bottom → top
    const bottom = visible[0];
    const onBase = bottom && resolveSurvey(bottom.survey || bottom.key) === this.base ? bottom : null;
    this._applyBase(onBase);
    const want = visible.filter((l) => l !== onBase).map((l) => l.key);
    const byKey = new Map(layers.map((l) => [l.key, l]));
    // Aladin adds overlays on top: keep the attached ones already in the
    // wanted order and re-attach the rest above them.
    let last = -1;
    let from = want.length;
    for (let i = 0; i < want.length; i++) {
      const pos = this.order.indexOf(want[i]);
      if (pos < 0 || pos < last) { from = i; break; }
      last = pos;
    }
    for (const key of want.slice(from)) if (this.overlays.has(key)) this._detach(key);
    for (const key of want.slice(from)) this._attach(byKey.get(key));
    const wanted = new Set(want);
    const now = performance.now();
    for (const key of this.order) {
      const o = this.overlays.get(key);
      const l = byKey.get(key);
      const on = wanted.has(key) && l;
      const opacity = on ? Math.min(1, Math.max(0, l.opacity ?? 1)) : 0;
      if (Math.abs(o.opacity - opacity) > 1e-4) {
        o.opacity = opacity;
        this._setOpacity(o, opacity);
      }
      if (on) {
        o.used = now;
        this._applyStyle(o, l, () => a.getOverlayImageLayer?.(o.name));
      }
    }
    this._trim(wanted);
  }

  /** Opacity and style of the base layer (0 when another survey is at the bottom). */
  _applyBase(l) {
    const a = this.aladin;
    const opacity = l ? Math.min(1, Math.max(0, l.opacity ?? 1)) : 0;
    const layer = (() => { try { return a.getBaseImageLayer?.(); } catch (_) { return null; } })();
    if (!layer) {
      // The base loads asynchronously at start: try again shortly.
      clearTimeout(this._baseRetry);
      this._baseRetry = setTimeout(() => this._applyBase(l), 250);
      return;
    }
    if (Math.abs((this.baseOpacity ?? -1) - opacity) > 1e-4) {
      try { layer.setOpacity?.(opacity); this.baseOpacity = opacity; } catch (_) { /* not ready */ }
    }
    if (l) this._applyStyle((this._baseStyle ||= {}), l, () => layer);
  }

  _attach(l) {
    if (!l) return;
    const a = this.aladin;
    const key = l.key;
    const name = `ov-${key}`;
    const settle = (ok, err) => {
      if (this.overlays.get(key)?.hips !== hips) return; // replaced meanwhile
      if (ok) {
        this.tries.delete(key);
        if (this.failed.delete(key)) this.onLayerStatus?.(key, "ready");
        return;
      }
      if (this.failed.has(key)) return;
      // A slow or busy server: attach again (twice) before calling the
      // survey unavailable.
      const tries = (this.tries.get(key) || 0) + 1;
      this.tries.set(key, tries);
      if (tries <= 2) {
        setTimeout(() => {
          if (this.overlays.get(key)?.hips !== hips) return;
          this._detach(key);
          this.onNeedApply?.();
        }, 1200 * tries);
        return;
      }
      this.failed.set(key, err || new Error("unavailable"));
      this.onLayerStatus?.(key, "error", err);
    };
    // Aladin's metadata query decides (it moves on to a mirror when a server
    // fails, and reports that failure first); the callbacks are a fallback.
    const verdict = (err) => {
      const q = hips?.query;
      if (q && typeof q.then === "function") q.then(() => settle(true), (e) => settle(false, e || err));
      else settle(false, err);
    };
    const hips = this._makeHiPS(l.survey || key, {
      successCallback: () => settle(true),
      errorCallback: (err) => verdict(err),
    });
    const q = hips?.query;
    if (q && typeof q.then === "function") q.then(() => settle(true), () => {});
    try { hips?.setOpacity?.(0); } catch (_) { /* not ready */ }
    try {
      if (typeof a.addImageLayer === "function") a.addImageLayer(hips, name);
      else a.setOverlayImageLayer(hips, name);
    } catch (err) {
      console.warn("Oviz sky: overlay failed", key, err);
      this.failed.set(key, err);
      this.onLayerStatus?.(key, "error", err);
      return;
    }
    this.overlays.set(key, { name, hips, opacity: 0, used: performance.now() });
    this.order.push(key);
  }

  /**
   * Attach a survey hidden (opacity 0) ahead of need, e.g. the backgrounds
   * saved views switch to, so their tiles are in when a transition shows them.
   */
  preload(l) {
    if (!this.aladin || !l?.key || this.overlays.has(l.key)) return;
    if (resolveSurvey(l.survey || l.key) === this.base) return;
    this._attach({ key: l.key, survey: l.survey || l.key });
    const o = this.overlays.get(l.key);
    if (o) o.used = 0; // first to go if the cache is full
  }

  /** Let go of a survey taken out of the figure (hidden ones stay cached). */
  forget(key) {
    if (this.overlays.has(key)) this._detach(key);
  }

  /** Keep at most MAX_HIDDEN hidden overlays attached (the most recently shown). */
  _trim(wanted) {
    const hidden = this.order.filter((k) => !wanted.has(k));
    if (hidden.length <= MAX_HIDDEN) return;
    hidden.sort((x, y) => this.overlays.get(x).used - this.overlays.get(y).used);
    for (const key of hidden.slice(0, hidden.length - MAX_HIDDEN)) this._detach(key);
  }

  _setOpacity(o, opacity) {
    const a = this.aladin;
    let layer = null;
    try { layer = a.getOverlayImageLayer?.(o.name); } catch (_) { layer = null; }
    const target = layer || o.hips;
    try { target?.setOpacity?.(opacity); } catch (_) { /* retried next apply */ }
    if (!layer) setTimeout(() => {
      try { a.getOverlayImageLayer?.(o.name)?.setOpacity?.(o.opacity); } catch (_) { /* ignore */ }
    }, 250);
  }

  /** Stretch, colormap and cuts of a layer (base or overlay). */
  _applyStyle(o, l, getLayer) {
    const sig = `${l.stretch || ""}|${l.colormap || ""}|${l.cut_min ?? ""}|${l.cut_max ?? ""}`;
    if (o.styleSig === sig) return;
    let layer = null;
    try { layer = getLayer(); } catch (_) { layer = null; }
    if (!layer) return;
    o.styleSig = sig;
    try {
      if (l.colormap || l.stretch || o.styled) {
        layer.setColormap?.(l.colormap || "native", { stretch: l.stretch || "linear" });
        o.styled = true;
      }
      // Like the classic viewer, cuts apply only when both ends are set.
      if (Number.isFinite(l.cut_min) && Number.isFinite(l.cut_max) && l.cut_max > l.cut_min) layer.setCuts?.(l.cut_min, l.cut_max);
    } catch (_) { /* optional features */ }
  }

  _detach(key) {
    const o = this.overlays.get(key);
    if (o) {
      try { this.aladin.removeImageLayer?.(o.name); } catch (_) { /* ignore */ }
    }
    this.overlays.delete(key);
    this.order = this.order.filter((k) => k !== key);
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
