// Self-export: rebuild this figure as a standalone HTML file.
//
// The page never regex-splices its own markup. Every piece it needs is kept
// verbatim in the DOM: the stylesheet, the manifest JSON, the base64 blob
// scripts (untouched, so data is never re-encoded) and the runtime script.
// Export swaps in a new manifest and reassembles the same skeleton the
// Python writer produces.

// The runtime itself is inlined in a script element, so its source must never
// contain a literal closing tag or comment opener: build them from parts.
const OPEN = "<" + "script";
const CLOSE = "<" + "/script>";

function jsonForScript(obj) {
  // \u003c keeps data from ending the element early and stays valid JSON.
  return JSON.stringify(obj).replace(/</g, "\\u003c");
}

function escapeHtml(s) {
  return String(s).replace(/[&<>"']/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[c]));
}

export function buildExportHtml(manifest, { title, theme = "dark", mode } = {}) {
  const doc = document;
  const modeId = mode || manifest.viewer?.mode || "focus";
  const style = doc.getElementById("oviz-style")?.textContent || "";
  const runtime = doc.getElementById("oviz-runtime")?.textContent || "";
  const blobs = [...doc.querySelectorAll("script[data-oviz-blob]")].map((el) => el.outerHTML).join("\n");
  const version = manifest.viewer?.version || "";
  const pageTitle = title ?? (doc.title || manifest.title || "Oviz figure");
  return [
    "<!doctype html>",
    `<html lang="en" data-oviz-theme="${escapeHtml(theme)}" data-oviz-mode="${escapeHtml(modeId)}">`,
    "<head>",
    '<meta charset="utf-8">',
    '<meta name="viewport" content="width=device-width, initial-scale=1, viewport-fit=cover">',
    `<meta name="generator" content="oviz ${escapeHtml(version)}">`,
    '<meta name="color-scheme" content="dark light">',
    `<title>${escapeHtml(pageTitle)}</title>`,
    `<style id="oviz-style">\n${style.replace(/^\n/, "")}</style>`,
    "</head>",
    "<body>",
    '<div id="oviz-root" class="ov-root" data-boot="loading">',
    '<div class="ov-boot" role="status" aria-live="polite"><div class="ov-boot-mark"><span></span><span></span><span></span></div><div class="ov-boot-text">Loading figure…</div><div class="ov-boot-bar"><i></i></div></div>',
    "</div>",
    `${OPEN} type="application/json" id="oviz-manifest">${jsonForScript(manifest)}${CLOSE}`,
    blobs,
    `${OPEN} id="oviz-runtime">\n${runtime.replace(/^\n/, "")}${CLOSE}`,
    "</body>",
    "</html>",
    "",
  ].join("\n");
}

let fileHandle = null;

/**
 * Save text as an HTML file: the File System Access API when available
 * (subsequent saves overwrite the same file), otherwise a download.
 */
export async function saveHtml(text, suggestedName, { reuseHandle = false } = {}) {
  if (document.fullscreenElement) {
    try { await document.exitFullscreen(); } catch (_) { /* ignore */ }
  }
  const blob = new Blob([text], { type: "text/html" });
  if (window.showSaveFilePicker && window.isSecureContext !== false) {
    try {
      let handle = reuseHandle ? fileHandle : null;
      if (!handle) {
        handle = await window.showSaveFilePicker({
          suggestedName,
          types: [{ description: "Oviz figure", accept: { "text/html": [".html"] } }],
        });
      }
      const w = await handle.createWritable();
      await w.write(blob);
      await w.close();
      if (reuseHandle) fileHandle = handle;
      return { ok: true, name: handle.name };
    } catch (err) {
      if (err && err.name === "AbortError") return { ok: false, cancelled: true };
      // Fall through to a download (e.g. cross-origin iframes).
    }
  }
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url;
  a.download = suggestedName;
  document.body.append(a);
  a.click();
  a.remove();
  setTimeout(() => URL.revokeObjectURL(url), 5000);
  return { ok: true, name: suggestedName, downloaded: true };
}
