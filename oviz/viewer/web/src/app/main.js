// Boot sequence: manifest → critical blobs → first frame → stream the rest.

import { readManifest, BundleStore } from "../core/loader.js";
import { Viewer } from "./viewer.js";
import { mountUI } from "../ui/app-ui.js";
import { installApi } from "./api.js";

async function boot() {
  const root = document.getElementById("oviz-root");
  const bootEl = root.querySelector(".ov-boot");
  const bootText = root.querySelector(".ov-boot-text");
  const bootBar = root.querySelector(".ov-boot-bar i");
  const t0 = performance.now();
  const timing = {};
  const setBoot = (text, frac) => {
    if (bootText && text) bootText.textContent = text;
    if (bootBar && frac != null) bootBar.style.transform = `scaleX(${Math.max(0.02, Math.min(frac, 1))})`;
  };
  try {
    const manifest = readManifest();
    const store = new BundleStore(manifest);
    const criticalBytes = (manifest.blobs || []).filter((b) => (b.priority ?? 1) <= 0).reduce((s, b) => s + b.bytes, 0) || 1;
    let loadedCritical = 0;
    const off = store.onProgress((_, __, d) => {
      if ((d.priority ?? 1) <= 0) {
        loadedCritical += d.bytes;
        setBoot(null, loadedCritical / criticalBytes);
      }
    });
    setBoot("Decoding scene…", 0.05);
    await store.loadPriority(0);
    off();
    timing.critical = performance.now() - t0;
    const stage = document.createElement("div");
    stage.className = "ov-stage";
    root.prepend(stage);
    const viewer = new Viewer(stage, manifest, store);
    await viewer.attachTraces();
    viewer.start();
    timing.firstFrame = performance.now() - t0;
    const ui = mountUI(root, viewer);
    installApi(root, viewer, ui);
    root.dataset.boot = "ready";
    bootEl?.remove();
    // Stream heavier assets without blocking interaction.
    viewer.attachImages();
    ui.onBackgroundLoad?.("start");
    viewer.attachVolumes(() => ui.onBackgroundLoad?.("volume")).then(() => {
      timing.volumes = performance.now() - t0;
      ui.onBackgroundLoad?.("volumes-done");
    });
    ui.afterBoot?.();
    root.dataset.bootMs = timing.firstFrame.toFixed(0);
    viewer.timing = timing;
  } catch (err) {
    console.error(err);
    root.dataset.boot = "error";
    setBoot(`This figure could not start: ${err.message || err}`, 1);
    bootEl?.classList.add("ov-boot--error");
  }
}

if (document.readyState === "loading") {
  document.addEventListener("DOMContentLoaded", boot, { once: true });
} else {
  boot();
}
