from __future__ import annotations

import base64
import json
import tempfile
import uuid
import webbrowser
from pathlib import Path
from typing import Any


_THREEJS_HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
  <head>
    <meta charset="utf-8" />
    <meta name="viewport" content="width=device-width, initial-scale=1" />
    <title>oviz three.js figure</title>
    <style>
      html, body {
        margin: 0;
        padding: 0;
        background: var(--oviz-paper-bg, #000);
      }
      #__ROOT_ID__ {
        --oviz-paper-bg: #000000;
        --oviz-scene-bg: #000000;
        --oviz-text: #d0d0d0;
        --oviz-axis: #808080;
        --oviz-panel-bg: rgba(0, 0, 0, 0.50);
        --oviz-panel-border: rgba(128, 128, 128, 0.55);
        --oviz-panel-solid: #121212;
        --oviz-footprint: #6ec5ff;
        position: relative;
        width: 100vw;
        height: 100vh;
        margin: 0;
        overflow: hidden;
        background: var(--oviz-paper-bg);
        color: var(--oviz-text);
        font-family: Helvetica, Arial, sans-serif;
      }
      #__ROOT_ID__ .oviz-three-canvas {
        display: block;
        width: 100%;
        height: 100%;
      }
      #__ROOT_ID__ .oviz-three-title {
        position: absolute;
        top: 12px;
        left: 50%;
        transform: translateX(-50%);
        z-index: 5;
        font-size: 20px;
        font-weight: 500;
        color: var(--oviz-text);
        pointer-events: none;
        white-space: nowrap;
      }
      #__ROOT_ID__ .oviz-three-toolbar {
        position: absolute;
        top: 12px;
        left: 12px;
        z-index: 6;
        display: flex;
        flex-direction: column;
        gap: 10px;
        max-width: min(320px, 32vw);
      }
      #__ROOT_ID__ .oviz-three-group-select {
        min-width: 180px;
        padding: 6px 10px;
        border-radius: 6px;
        border: 1px solid var(--oviz-panel-border);
        background: var(--oviz-panel-bg);
        color: var(--oviz-text);
        font-size: 14px;
        font-family: inherit;
      }
      #__ROOT_ID__ .oviz-three-legend {
        display: flex;
        flex-direction: column;
        gap: 6px;
        padding: 10px 12px;
        border-radius: 8px;
        border: 1px solid var(--oviz-panel-border);
        background: var(--oviz-panel-bg);
        max-height: 58vh;
        overflow: auto;
      }
      #__ROOT_ID__ .oviz-three-legend-title {
        font-size: 14px;
        font-weight: 600;
        color: var(--oviz-text);
      }
      #__ROOT_ID__ .oviz-three-legend-item {
        display: flex;
        align-items: center;
        gap: 8px;
        font-size: 13px;
        color: var(--oviz-text);
      }
      #__ROOT_ID__ .oviz-three-legend-item input {
        margin: 0;
      }
      #__ROOT_ID__ .oviz-three-footer {
        position: absolute;
        left: 50%;
        bottom: 14px;
        transform: translateX(-50%);
        z-index: 6;
        display: flex;
        align-items: center;
        gap: 10px;
        width: min(88vw, 1100px);
        padding: 8px 12px;
        border-radius: 999px;
        border: 1px solid var(--oviz-panel-border);
        background: var(--oviz-panel-bg);
        backdrop-filter: blur(6px);
      }
      #__ROOT_ID__ .oviz-three-footer button {
        border: 1px solid var(--oviz-panel-border);
        border-radius: 999px;
        background: transparent;
        color: var(--oviz-text);
        cursor: pointer;
        font-size: 15px;
        width: 34px;
        height: 34px;
        line-height: 1;
      }
      #__ROOT_ID__ .oviz-three-time-label {
        min-width: 136px;
        font-size: 15px;
        color: var(--oviz-text);
        white-space: nowrap;
      }
      #__ROOT_ID__ .oviz-three-slider {
        flex: 1 1 auto;
        accent-color: var(--oviz-axis);
      }
      #__ROOT_ID__ .oviz-three-tooltip {
        position: absolute;
        z-index: 7;
        display: none;
        min-width: 120px;
        max-width: 320px;
        padding: 8px 10px;
        border-radius: 8px;
        border: 1px solid var(--oviz-panel-border);
        background: rgba(0, 0, 0, 0.82);
        color: #ffffff;
        font-size: 12px;
        line-height: 1.4;
        pointer-events: none;
        white-space: normal;
      }
      #__ROOT_ID__ .oviz-three-note {
        position: absolute;
        right: 12px;
        bottom: 64px;
        z-index: 6;
        max-width: min(340px, 34vw);
        padding: 8px 10px;
        border-radius: 8px;
        border: 1px solid var(--oviz-panel-border);
        background: var(--oviz-panel-bg);
        color: var(--oviz-text);
        font-size: 12px;
        display: none;
      }
      #__ROOT_ID__ .oviz-three-sky-panel {
        position: absolute;
        top: 72px;
        right: 12px;
        width: min(38vw, 560px);
        height: min(56vh, 560px);
        z-index: 8;
        display: none;
        overflow: hidden;
        border: 1px solid var(--oviz-panel-border);
        border-radius: 10px;
        background: var(--oviz-panel-solid);
        box-shadow: 0 18px 46px rgba(0, 0, 0, 0.30);
        touch-action: none;
      }
      #__ROOT_ID__ .oviz-three-sky-panel[data-mode="normal"] {
        display: block;
      }
      #__ROOT_ID__ .oviz-three-sky-panel[data-mode="fullscreen"] {
        display: block;
        top: 0;
        right: 0;
        bottom: 0;
        left: 0;
        width: auto;
        height: auto;
        border-radius: 0;
        border-width: 0;
      }
      #__ROOT_ID__ .oviz-three-sky-drag {
        position: absolute;
        top: 0;
        left: 0;
        right: 0;
        height: 28px;
        z-index: 2;
        display: flex;
        align-items: center;
        padding: 0 10px;
        font-size: 11px;
        letter-spacing: 0.04em;
        text-transform: uppercase;
        color: var(--oviz-text);
        background: linear-gradient(to bottom, rgba(0, 0, 0, 0.18), rgba(0, 0, 0, 0.0));
        cursor: grab;
        user-select: none;
      }
      #__ROOT_ID__ .oviz-three-sky-panel[data-mode="fullscreen"] .oviz-three-sky-drag {
        cursor: default;
      }
      #__ROOT_ID__ .oviz-three-sky-frame {
        position: absolute;
        top: 0;
        left: 0;
        width: 100%;
        height: calc(100% - 72px);
        border: 0;
        display: block;
        background: var(--oviz-panel-solid);
      }
      #__ROOT_ID__ .oviz-three-sky-readout {
        position: absolute;
        right: 0;
        bottom: 0;
        left: 0;
        height: 72px;
        margin: 0;
        padding: 8px 10px;
        overflow: auto;
        border-top: 1px solid var(--oviz-panel-border);
        background: rgba(0, 0, 0, 0.18);
        color: var(--oviz-text);
        font: 11px/1.45 Menlo, Monaco, Consolas, monospace;
        white-space: pre-wrap;
      }
      #__ROOT_ID__ .oviz-three-sky-panel button,
      #__ROOT_ID__ .oviz-three-sky-actions button {
        border: 1px solid var(--oviz-panel-border);
        border-radius: 4px;
        background: var(--oviz-panel-bg);
        color: var(--oviz-text);
        cursor: pointer;
        font: 10px Menlo, Monaco, Consolas, monospace;
        padding: 2px 6px;
      }
      #__ROOT_ID__ .oviz-three-sky-full,
      #__ROOT_ID__ .oviz-three-sky-hide {
        position: absolute;
        top: 6px;
        z-index: 3;
      }
      #__ROOT_ID__ .oviz-three-sky-full {
        right: 52px;
      }
      #__ROOT_ID__ .oviz-three-sky-hide {
        right: 6px;
      }
      #__ROOT_ID__ .oviz-three-sky-actions {
        position: absolute;
        right: 12px;
        top: 72px;
        z-index: 8;
        display: none;
        gap: 8px;
      }
      #__ROOT_ID__ .oviz-three-sky-resize {
        position: absolute;
        width: 14px;
        height: 14px;
        z-index: 3;
        background: rgba(120, 120, 120, 0.10);
        display: block;
        touch-action: none;
      }
      #__ROOT_ID__ .oviz-three-sky-panel[data-mode="fullscreen"] .oviz-three-sky-resize {
        display: none;
      }
      #__ROOT_ID__ .oviz-three-sky-resize[data-dir="nw"] {
        top: 0;
        left: 0;
        cursor: nwse-resize;
      }
      #__ROOT_ID__ .oviz-three-sky-resize[data-dir="ne"] {
        top: 0;
        right: 0;
        cursor: nesw-resize;
      }
      #__ROOT_ID__ .oviz-three-sky-resize[data-dir="sw"] {
        bottom: 0;
        left: 0;
        cursor: nesw-resize;
      }
      #__ROOT_ID__ .oviz-three-sky-resize[data-dir="se"] {
        right: 0;
        bottom: 0;
        cursor: nwse-resize;
      }
    </style>
  </head>
  <body>
    <div id="__ROOT_ID__">
      <div class="oviz-three-title"></div>
      <div class="oviz-three-toolbar">
        <select class="oviz-three-group-select"></select>
        <div class="oviz-three-legend"></div>
      </div>
      <canvas class="oviz-three-canvas"></canvas>
      <div class="oviz-three-tooltip"></div>
      <div class="oviz-three-footer">
        <button class="oviz-three-play" type="button" title="Play">▶</button>
        <button class="oviz-three-pause" type="button" title="Pause">⏸</button>
        <span class="oviz-three-time-label"></span>
        <input class="oviz-three-slider" type="range" min="0" max="0" step="1" value="0" />
      </div>
      <div class="oviz-three-note"></div>
      <div class="oviz-three-sky-panel" data-mode="hidden">
        <div class="oviz-three-sky-drag">Sky</div>
        <iframe class="oviz-three-sky-frame" loading="eager" referrerpolicy="no-referrer"></iframe>
        <pre class="oviz-three-sky-readout"></pre>
        <button class="oviz-three-sky-full" type="button" title="Toggle fullscreen sky panel">Full</button>
        <button class="oviz-three-sky-hide" type="button" title="Hide sky panel">Hide</button>
        <div class="oviz-three-sky-resize" data-dir="nw"></div>
        <div class="oviz-three-sky-resize" data-dir="ne"></div>
        <div class="oviz-three-sky-resize" data-dir="sw"></div>
        <div class="oviz-three-sky-resize" data-dir="se"></div>
      </div>
      <div class="oviz-three-sky-actions">
        <button class="oviz-three-sky-show" type="button" title="Show sky panel">Sky</button>
        <button class="oviz-three-sky-show-full" type="button" title="Show fullscreen sky panel">Sky Full</button>
      </div>
    </div>

    <script>
      (async function () {
      const dependencyGroups = [
        [
          "https://cdn.jsdelivr.net/npm/three@0.128.0/build/three.min.js",
          "https://unpkg.com/three@0.128.0/build/three.min.js"
        ],
        [
          "https://cdn.jsdelivr.net/npm/three@0.128.0/examples/js/controls/OrbitControls.js",
          "https://unpkg.com/three@0.128.0/examples/js/controls/OrbitControls.js"
        ],
        [
          "https://cdn.jsdelivr.net/npm/three@0.128.0/examples/js/lines/LineSegmentsGeometry.js",
          "https://unpkg.com/three@0.128.0/examples/js/lines/LineSegmentsGeometry.js"
        ],
        [
          "https://cdn.jsdelivr.net/npm/three@0.128.0/examples/js/lines/LineMaterial.js",
          "https://unpkg.com/three@0.128.0/examples/js/lines/LineMaterial.js"
        ],
        [
          "https://cdn.jsdelivr.net/npm/three@0.128.0/examples/js/lines/LineSegments2.js",
          "https://unpkg.com/three@0.128.0/examples/js/lines/LineSegments2.js"
        ]
      ];

      function renderError(root, message) {
        if (!root) {
          return;
        }
        root.innerHTML = "";
        const box = document.createElement("div");
        box.style.cssText = [
          "height: 100%",
          "display: flex",
          "align-items: center",
          "justify-content: center",
          "padding: 24px",
          "background: #111",
          "color: #eee",
          "font: 14px Helvetica, Arial, sans-serif",
          "text-align: center",
          "white-space: pre-wrap"
        ].join(";");
        box.textContent = message;
        root.appendChild(box);
      }

      function loadScript(url) {
        return new Promise((resolve, reject) => {
          const script = document.createElement("script");
          script.src = url;
          script.async = true;
          script.onload = () => resolve(url);
          script.onerror = () => reject(new Error(url));
          document.head.appendChild(script);
        });
      }

      async function loadDependencyGroup(urls) {
        let lastError = null;
        for (const url of urls) {
          try {
            await loadScript(url);
            return url;
          } catch (err) {
            lastError = err;
          }
        }
        throw lastError || new Error(urls.join(", "));
      }

      try {
        for (const group of dependencyGroups) {
          await loadDependencyGroup(group);
        }
      } catch (err) {
        const root = document.getElementById("__ROOT_ID__");
        renderError(
          root,
          "oviz three.js renderer could not load its browser dependencies.\\n\\n"
          + "Tried jsDelivr and unpkg CDNs.\\n"
          + "If this environment blocks network requests, the Three.js figure cannot render yet.\\n\\n"
          + "Last failed URL: " + (err && err.message ? err.message : String(err))
        );
        return;
      }

      const THREE = window.THREE;
      const OrbitControls = THREE && THREE.OrbitControls;
      const LineMaterial = THREE && THREE.LineMaterial;
      const LineSegments2 = THREE && THREE.LineSegments2;
      const LineSegmentsGeometry = THREE && THREE.LineSegmentsGeometry;
      if (!THREE || !OrbitControls || !LineMaterial || !LineSegments2 || !LineSegmentsGeometry) {
        const root = document.getElementById("__ROOT_ID__");
        renderError(
          root,
          "oviz three.js renderer loaded partial browser dependencies but the expected Three.js helpers were not available."
        );
        return;
      }

      const sceneSpec = __SCENE_JSON__;
      const root = document.getElementById("__ROOT_ID__");
      const canvas = root.querySelector(".oviz-three-canvas");
      const titleEl = root.querySelector(".oviz-three-title");
      const groupSelectEl = root.querySelector(".oviz-three-group-select");
      const legendEl = root.querySelector(".oviz-three-legend");
      const sliderEl = root.querySelector(".oviz-three-slider");
      const timeLabelEl = root.querySelector(".oviz-three-time-label");
      const playButtonEl = root.querySelector(".oviz-three-play");
      const pauseButtonEl = root.querySelector(".oviz-three-pause");
      const tooltipEl = root.querySelector(".oviz-three-tooltip");
      const noteEl = root.querySelector(".oviz-three-note");
      const skyPanelEl = root.querySelector(".oviz-three-sky-panel");
      const skyFrameEl = root.querySelector(".oviz-three-sky-frame");
      const skyReadoutEl = root.querySelector(".oviz-three-sky-readout");
      const skyFullButtonEl = root.querySelector(".oviz-three-sky-full");
      const skyHideButtonEl = root.querySelector(".oviz-three-sky-hide");
      const skyShowButtonEl = root.querySelector(".oviz-three-sky-show");
      const skyShowFullButtonEl = root.querySelector(".oviz-three-sky-show-full");
      const skyActionsEl = root.querySelector(".oviz-three-sky-actions");
      const skyDragHandleEl = root.querySelector(".oviz-three-sky-drag");
      const skyResizeEls = Array.from(root.querySelectorAll(".oviz-three-sky-resize"));

      const theme = sceneSpec.theme;
      root.style.setProperty("--oviz-paper-bg", theme.paper_bgcolor || "#000000");
      root.style.setProperty("--oviz-scene-bg", theme.scene_bgcolor || theme.paper_bgcolor || "#000000");
      root.style.setProperty("--oviz-text", theme.text_color || "#d0d0d0");
      root.style.setProperty("--oviz-axis", theme.axis_color || "#808080");
      root.style.setProperty("--oviz-panel-bg", theme.panel_bg || "rgba(0, 0, 0, 0.45)");
      root.style.setProperty("--oviz-panel-border", theme.panel_border || "rgba(128, 128, 128, 0.50)");
      root.style.setProperty("--oviz-panel-solid", theme.panel_solid || theme.paper_bgcolor || "#121212");
      root.style.setProperty("--oviz-footprint", theme.footprint || "#6ec5ff");
      titleEl.textContent = sceneSpec.title || "";
      titleEl.style.display = sceneSpec.title ? "block" : "none";

      if (sceneSpec.note) {
        noteEl.textContent = sceneSpec.note;
        noteEl.style.display = "block";
      }

      const layoutScene = (sceneSpec.layout || {}).scene || {};
      const axisSpec = sceneSpec.axes || {};
      const frameSpecs = sceneSpec.frames || [];
      const legendItems = (sceneSpec.legend || {}).items || [];
      const groupVisibility = sceneSpec.group_visibility || {};
      const defaultGroup = sceneSpec.default_group || "All";
      const skySpec = sceneSpec.sky_panel || { enabled: false };
      const pointScale = Math.max(sceneSpec.max_span || 1, 1) / 4000.0;
      const hoverTargets = [];
      const axisLineMaterials = [];
      const frameLineMaterials = [];
      const markerTextureCache = new Map();
      const markerMaterialCache = new Map();
      const textTextureCache = new Map();
      const galaxyTextureCache = new Map();
      let legendState = {};
      let currentGroup = defaultGroup;
      let currentFrameIndex = sceneSpec.initial_frame_index || 0;
      let currentSelection = null;
      let playbackTimer = null;
      let skyPanelMode = skySpec.enabled ? "normal" : "hidden";
      let skyPointerState = null;

      const renderer = new THREE.WebGLRenderer({
        canvas,
        antialias: true,
        alpha: false,
        powerPreference: "high-performance",
      });
      renderer.setPixelRatio(window.devicePixelRatio || 1);
      if ("outputColorSpace" in renderer && THREE.SRGBColorSpace) {
        renderer.outputColorSpace = THREE.SRGBColorSpace;
      } else if ("outputEncoding" in renderer && THREE.sRGBEncoding) {
        renderer.outputEncoding = THREE.sRGBEncoding;
      }

      const scene = new THREE.Scene();
      scene.background = new THREE.Color(theme.scene_bgcolor || theme.paper_bgcolor || "#000000");

      const sceneUp = sceneSpec.camera_up || { x: 0.0, y: 0.0, z: 1.0 };
      const camera = new THREE.PerspectiveCamera(38, 1, 0.1, Math.max((sceneSpec.max_span || 1) * 20.0, 10000.0));
      camera.up.set(sceneUp.x ?? 0.0, sceneUp.y ?? 0.0, sceneUp.z ?? 1.0);
      const controls = new OrbitControls(camera, renderer.domElement);
      controls.enableDamping = true;
      controls.dampingFactor = 0.08;
      controls.rotateSpeed = 0.7;
      controls.panSpeed = 0.7;
      controls.minPolarAngle = 0.02;
      controls.maxPolarAngle = Math.PI - 0.02;
      controls.target.set(sceneSpec.center.x, sceneSpec.center.y, sceneSpec.center.z);

      const eye = (layoutScene.camera || {}).eye || { x: 0.0, y: -0.8, z: 1.5 };
      const eyeScale = Math.max(sceneSpec.max_span || 1, 1);
      camera.position.set(
        sceneSpec.center.x + (eye.x ?? 0.0) * eyeScale,
        sceneSpec.center.y + (eye.y ?? -0.8) * eyeScale,
        sceneSpec.center.z + (eye.z ?? 1.5) * eyeScale
      );
      camera.lookAt(controls.target);

      const plotGroup = new THREE.Group();
      scene.add(plotGroup);

      const axisGroup = new THREE.Group();
      scene.add(axisGroup);

      const raycaster = new THREE.Raycaster();
      const pointer = new THREE.Vector2(2, 2);

      function formatTick(value) {
        if (!Number.isFinite(value)) {
          return "";
        }
        const rounded = Math.round(value);
        if (Math.abs(value - rounded) < 1e-6) {
          return String(rounded);
        }
        return value.toFixed(1);
      }

      function approximatelyZero(value) {
        return Number.isFinite(value) && Math.abs(Number(value)) <= 1e-9;
      }

      function currentFrame() {
        return frameSpecs[Math.max(0, Math.min(currentFrameIndex, frameSpecs.length - 1))] || null;
      }

      function clampSkyPosition(left, top, width, height) {
        const margin = 6;
        return {
          left: Math.min(Math.max(margin, left), Math.max(margin, window.innerWidth - width - margin)),
          top: Math.min(Math.max(margin, top), Math.max(margin, window.innerHeight - height - margin)),
        };
      }

      function resizeSkyRect(state, clientX, clientY) {
        const margin = 6;
        const minWidth = Math.min(220, Math.max(80, window.innerWidth - 2 * margin));
        const minHeight = Math.min(220, Math.max(80, window.innerHeight - 2 * margin));
        let left = state.startLeft;
        let right = state.startRight;
        let top = state.startTop;
        let bottom = state.startBottom;
        const dx = clientX - state.startX;
        const dy = clientY - state.startY;
        const dir = state.dir || "se";

        if (dir.includes("w")) {
          left = Math.max(margin, Math.min(state.startLeft + dx, state.startRight - minWidth));
        } else if (dir.includes("e")) {
          right = Math.min(window.innerWidth - margin, Math.max(state.startRight + dx, state.startLeft + minWidth));
        }
        if (dir.includes("n")) {
          top = Math.max(margin, Math.min(state.startTop + dy, state.startBottom - minHeight));
        } else if (dir.includes("s")) {
          bottom = Math.min(window.innerHeight - margin, Math.max(state.startBottom + dy, state.startTop + minHeight));
        }

        return {
          left,
          top,
          width: Math.max(minWidth, right - left),
          height: Math.max(minHeight, bottom - top),
        };
      }

      function skyReadoutText(selection) {
        if (!selection) {
          return "Click a cluster member at t=0 to load sky imagery.";
        }
        const clusterLabel = selection.trace_name || selection.cluster_name || "";
        const clusterName = clusterLabel ? `Cluster: ${clusterLabel}\n` : "";
        return (
          clusterName
          + `Selected direction (t=0): l=${Number(selection.l_deg).toFixed(4)} deg, b=${Number(selection.b_deg).toFixed(4)} deg\n`
          + `Beam radius: ${Number(skySpec.radius_deg || 1.0).toFixed(2)} deg\n`
          + `Survey: ${String(skySpec.survey || "P/DSS2/color")}`
        );
      }

      function buildEmptySkySrcdoc() {
        const bg = String(theme.panel_solid || theme.paper_bgcolor || "#121212");
        const txt = String(theme.text_color || "#d0d0d0");
        return (
          "<!doctype html><html><head><meta charset='utf-8'></head>"
          + "<body style=\\"margin:0;padding:0;background:" + bg + ";color:" + txt + ";"
          + "display:flex;align-items:center;justify-content:center;height:100vh;"
          + "font-family:Helvetica,Arial,sans-serif;font-size:14px;text-align:center;\\">"
          + "Click a cluster member at t=0 to open Aladin Lite."
          + "</body></html>"
        );
      }

      function normalizeMemberKey(value) {
        return String(value || "")
          .trim()
          .toLowerCase()
          .replace(/[_\s]+/g, " ");
      }

      function resolveMemberPoints(selection) {
        const byCluster = skySpec.members_by_cluster || {};
        const traceName = selection && selection.trace_name ? String(selection.trace_name) : "";
        const clusterName = selection && selection.cluster_name ? String(selection.cluster_name) : "";
        const candidates = [clusterName, traceName]
          .filter(Boolean)
          .flatMap((name) => [name, name.replace(/_/g, " "), name.replace(/\s+/g, "_")]);

        for (const key of candidates) {
          if (Array.isArray(byCluster[key]) && byCluster[key].length) {
            return { key, points: byCluster[key] };
          }
        }

        const normalizedCandidates = new Set(candidates.map((name) => normalizeMemberKey(name)));
        for (const key of Object.keys(byCluster)) {
          if (normalizedCandidates.has(normalizeMemberKey(key)) && Array.isArray(byCluster[key]) && byCluster[key].length) {
            return { key, points: byCluster[key] };
          }
        }

        return { key: traceName || clusterName || "Selection", points: null };
      }

      function buildAladinCatalogPayload(selection) {
        if (!selection) {
          return [];
        }
        const resolvedMembers = resolveMemberPoints(selection);
        const lookupName = resolvedMembers.key;
        const clusterColor = selection.cluster_color ? String(selection.cluster_color) : "#ffffff";
        const members = Array.isArray(resolvedMembers.points) ? resolvedMembers.points : null;
        let points = [];
        if (members && members.length) {
          points = members.map((pt) => ({
            l: Number(pt.l),
            b: Number(pt.b),
            ra: Number(pt.ra),
            dec: Number(pt.dec),
            label: pt.label || lookupName,
          }));
        } else if (Number.isFinite(Number(selection.ra_deg)) && Number.isFinite(Number(selection.dec_deg))) {
          points = [{
            l: Number(selection.l_deg),
            b: Number(selection.b_deg),
            ra: Number(selection.ra_deg),
            dec: Number(selection.dec_deg),
            label: lookupName,
          }];
        }
        return [{
          name: lookupName,
          color: clusterColor,
          opacity: 1.0,
          sourceSize: points.length > 1 ? 4 : 7,
          points,
        }];
      }

      function buildAladinSrcdoc(selection, catalogPayload) {
        if (!selection || !Number.isFinite(Number(selection.ra_deg)) || !Number.isFinite(Number(selection.dec_deg))) {
          return buildEmptySkySrcdoc();
        }
        const bg = String(theme.panel_solid || theme.paper_bgcolor || "#121212");
        const txt = String(theme.text_color || "#d0d0d0");
        const beamColor = JSON.stringify(String(theme.footprint || "#6ec5ff"));
        const survey = JSON.stringify(String(skySpec.survey || "P/DSS2/color"));
        const cooFrame = JSON.stringify(String(skySpec.frame || "galactic") === "galactic" ? "galactic" : "equatorial");
        const payloadJson = JSON.stringify(catalogPayload || []);
        const target = JSON.stringify(`${Number(selection.ra_deg).toFixed(6)} ${Number(selection.dec_deg).toFixed(6)}`);
        const radiusDeg = Number(skySpec.radius_deg || 1.0);
        const fovDeg = Math.min(Math.max(radiusDeg * 2.4, 1.2), 180.0);
        const ra = Number(selection.ra_deg);
        const dec = Number(selection.dec_deg);
        return `<!doctype html>
<html>
  <head>
    <meta charset="utf-8" />
    <meta name="viewport" content="width=device-width, initial-scale=1" />
    <style>
      html, body { margin: 0; padding: 0; width: 100%; height: 100%; background: ${bg}; color: ${txt}; overflow: hidden; }
      #oviz-wrap { position: relative; width: 100%; height: 100%; }
      #aladin-lite-div { width: 100%; height: 100%; }
      #oviz-status {
        position: absolute;
        inset: 0;
        display: flex;
        align-items: center;
        justify-content: center;
        font: 14px Helvetica, Arial, sans-serif;
        padding: 18px;
        text-align: center;
      }
    </style>
  </head>
  <body>
    <div id="oviz-wrap">
      <div id="aladin-lite-div"></div>
      <div id="oviz-status">Loading Aladin Lite...</div>
    </div>
    <script src="https://aladin.u-strasbg.fr/AladinLite/api/v3/latest/aladin.js" charset="utf-8"><\/script>
    <script>
      (function () {
        const payload = ${payloadJson};
        const statusEl = document.getElementById("oviz-status");
        function fail(message) {
          if (statusEl) {
            statusEl.textContent = message;
          }
        }
        if (typeof window.A === "undefined") {
          fail("Aladin Lite failed to load.");
          return;
        }
        A.init.then(() => {
          const aladin = A.aladin("#aladin-lite-div", {
            survey: ${survey},
            fov: ${fovDeg},
            target: ${target},
            cooFrame: ${cooFrame},
            showReticle: true,
            showLayersControl: true,
            showGotoControl: true,
            showFrame: true
          });
          if (aladin && typeof aladin.setImageSurvey === "function") {
            aladin.setImageSurvey(${survey});
          }
          if (aladin && typeof aladin.gotoRaDec === "function") {
            aladin.gotoRaDec(${ra.toFixed(8)}, ${dec.toFixed(8)});
          }
          const beam = A.graphicOverlay({ color: ${beamColor}, lineWidth: 2, opacity: 0.95 });
          aladin.addOverlay(beam);
          beam.add(A.circle(${ra.toFixed(8)}, ${dec.toFixed(8)}, ${radiusDeg.toFixed(8)}));
          payload.forEach((catDef) => {
            const cat = A.catalog({
              name: catDef.name,
              color: catDef.color,
              sourceSize: catDef.sourceSize || 7,
              shape: "circle",
              opacity: catDef.opacity
            });
            aladin.addCatalog(cat);
            const sources = [];
            (catDef.points || []).forEach((pt) => {
              sources.push(A.source(Number(pt.ra), Number(pt.dec), { popupTitle: pt.label || catDef.name }));
            });
            if (sources.length) {
              cat.addSources(sources);
            }
            if ((catDef.points || []).length > 1) {
              const membersOverlay = A.graphicOverlay({ color: catDef.color, lineWidth: 1, opacity: 0.90 });
              aladin.addOverlay(membersOverlay);
              const markerRadius = Math.max(Math.min(${fovDeg} / 220.0, 0.04), 0.002);
              (catDef.points || []).forEach((pt) => {
                membersOverlay.add(A.circle(Number(pt.ra), Number(pt.dec), markerRadius));
              });
            }
          });
          if (statusEl) {
            statusEl.style.display = "none";
          }
        }).catch((err) => {
          fail("Aladin initialization error: " + (err && err.message ? err.message : String(err)));
        });
      })();
    <\/script>
  </body>
</html>`;
      }

      function updateSkyPanel() {
        if (!skySpec.enabled) {
          return;
        }
        if (!currentSelection) {
          skyFrameEl.srcdoc = buildEmptySkySrcdoc();
          skyReadoutEl.textContent = skyReadoutText(null);
          return;
        }
        const payload = buildAladinCatalogPayload(currentSelection);
        skyFrameEl.srcdoc = buildAladinSrcdoc(currentSelection, payload);
        skyReadoutEl.textContent = skyReadoutText(currentSelection);
      }

      function applySkyPanelMode() {
        if (!skySpec.enabled) {
          skyPanelEl.dataset.mode = "hidden";
          skyActionsEl.style.display = "none";
          return;
        }
        skyPanelEl.dataset.mode = skyPanelMode;
        skyActionsEl.style.display = skyPanelMode === "hidden" ? "flex" : "none";
        skyFullButtonEl.textContent = skyPanelMode === "fullscreen" ? "Window" : "Full";
      }

      function clearGroup(group) {
        while (group.children.length) {
          const child = group.children[group.children.length - 1];
          if (child.children && child.children.length) {
            clearGroup(child);
          }
          if (child.material) {
            if (Array.isArray(child.material)) {
              child.material.forEach((material) => {
                if (!material || material.userData?.cached) {
                  return;
                }
                if (material.map) {
                  material.map.dispose();
                }
                material.dispose();
              });
            } else {
              if (!child.material.userData?.cached && child.material.map) {
                child.material.map.dispose();
              }
              if (!child.material.userData?.cached) {
                child.material.dispose();
              }
            }
          }
          if (child.geometry) {
            child.geometry.dispose();
          }
          group.remove(child);
        }
      }

      function markerTextureFor(symbol) {
        const key = symbol || "circle";
        if (markerTextureCache.has(key)) {
          return markerTextureCache.get(key);
        }

        const canvasEl = document.createElement("canvas");
        canvasEl.width = 128;
        canvasEl.height = 128;
        const ctx = canvasEl.getContext("2d");
        ctx.clearRect(0, 0, canvasEl.width, canvasEl.height);
        ctx.fillStyle = "#ffffff";
        ctx.strokeStyle = "#ffffff";
        ctx.lineWidth = 10;
        ctx.lineJoin = "round";
        ctx.lineCap = "round";
        ctx.translate(64, 64);

        switch (key) {
          case "square":
            ctx.fillRect(-34, -34, 68, 68);
            break;
          case "diamond":
            ctx.beginPath();
            ctx.moveTo(0, -40);
            ctx.lineTo(40, 0);
            ctx.lineTo(0, 40);
            ctx.lineTo(-40, 0);
            ctx.closePath();
            ctx.fill();
            break;
          case "cross":
            ctx.beginPath();
            ctx.moveTo(-14, -40);
            ctx.lineTo(14, -40);
            ctx.lineTo(14, -14);
            ctx.lineTo(40, -14);
            ctx.lineTo(40, 14);
            ctx.lineTo(14, 14);
            ctx.lineTo(14, 40);
            ctx.lineTo(-14, 40);
            ctx.lineTo(-14, 14);
            ctx.lineTo(-40, 14);
            ctx.lineTo(-40, -14);
            ctx.lineTo(-14, -14);
            ctx.closePath();
            ctx.fill();
            break;
          case "x":
            ctx.rotate(Math.PI / 4.0);
            ctx.fillRect(-14, -40, 28, 80);
            ctx.fillRect(-40, -14, 80, 28);
            break;
          case "triangle-up":
            ctx.beginPath();
            ctx.moveTo(0, -42);
            ctx.lineTo(42, 34);
            ctx.lineTo(-42, 34);
            ctx.closePath();
            ctx.fill();
            break;
          default:
            ctx.beginPath();
            ctx.arc(0, 0, 38, 0, Math.PI * 2.0);
            ctx.fill();
            break;
        }

        const texture = new THREE.CanvasTexture(canvasEl);
        texture.colorSpace = THREE.SRGBColorSpace;
        texture.needsUpdate = true;
        markerTextureCache.set(key, texture);
        return texture;
      }

      function markerMaterialFor(symbol, color, opacity) {
        const cacheKey = [symbol ?? "circle", color ?? "#ffffff", opacity ?? 1.0].join("|");
        if (markerMaterialCache.has(cacheKey)) {
          return markerMaterialCache.get(cacheKey);
        }
        const material = new THREE.SpriteMaterial({
          map: markerTextureFor(symbol),
          color: color ?? "#ffffff",
          transparent: true,
          opacity: opacity ?? 1.0,
          alphaTest: 0.15,
          depthWrite: false,
          depthTest: true,
        });
        material.userData = { cached: true };
        markerMaterialCache.set(cacheKey, material);
        return material;
      }

      function textTextureFor(text, color, size, family) {
        const cacheKey = [text, color ?? "#ffffff", size ?? 12, family ?? "Helvetica"].join("|");
        if (textTextureCache.has(cacheKey)) {
          return textTextureCache.get(cacheKey);
        }
        const fontSize = Math.max(10, size ?? 12);
        const resolutionScale = Math.max(2.0, Math.min(window.devicePixelRatio || 1.0, 4.0));
        const canvasEl = document.createElement("canvas");
        const ctx = canvasEl.getContext("2d");
        ctx.font = `${fontSize}px ${family ?? "Helvetica"}`;
        const metrics = ctx.measureText(text);
        const paddingX = 18;
        const paddingY = 12;
        const logicalWidth = Math.max(4, Math.ceil(metrics.width + 2 * paddingX));
        const logicalHeight = Math.max(4, Math.ceil(fontSize + 2 * paddingY));
        canvasEl.width = Math.max(4, Math.ceil(logicalWidth * resolutionScale));
        canvasEl.height = Math.max(4, Math.ceil(logicalHeight * resolutionScale));
        const drawCtx = canvasEl.getContext("2d");
        drawCtx.clearRect(0, 0, canvasEl.width, canvasEl.height);
        drawCtx.setTransform(resolutionScale, 0, 0, resolutionScale, 0, 0);
        drawCtx.font = `${fontSize}px ${family ?? "Helvetica"}`;
        drawCtx.fillStyle = color ?? "#ffffff";
        drawCtx.textBaseline = "middle";
        drawCtx.fillText(text, paddingX, logicalHeight / 2);
        const texture = new THREE.CanvasTexture(canvasEl);
        texture.colorSpace = THREE.SRGBColorSpace;
        texture.generateMipmaps = false;
        texture.minFilter = THREE.LinearFilter;
        texture.magFilter = THREE.LinearFilter;
        texture.needsUpdate = true;
        texture.userData = { width: logicalWidth, height: logicalHeight };
        textTextureCache.set(cacheKey, texture);
        return texture;
      }

      function makeTextSprite(text, options = {}) {
        const texture = textTextureFor(text, options.color, options.size, options.family);
        const material = new THREE.SpriteMaterial({
          map: texture,
          transparent: true,
          alphaTest: 0.05,
          depthWrite: false,
          depthTest: true,
        });
        const sprite = new THREE.Sprite(material);
        const aspect = texture.userData.width / Math.max(texture.userData.height, 1);
        const baseScale = (sceneSpec.max_span || 1) / 1400.0;
        const sizeScale = Math.max((options.size ?? 12) / 12.0, 0.8);
        sprite.scale.set(baseScale * aspect * 120.0 * sizeScale, baseScale * 120.0 * sizeScale, 1.0);
        return sprite;
      }

      function galaxyTextureFor(kind) {
        if (galaxyTextureCache.has(kind)) {
          return galaxyTextureCache.get(kind);
        }

        const size = 1024;
        const canvasEl = document.createElement("canvas");
        canvasEl.width = size;
        canvasEl.height = size;
        const ctx = canvasEl.getContext("2d");
        ctx.clearRect(0, 0, size, size);
        ctx.translate(size / 2, size / 2);

        function drawGlow(radius, stops) {
          const gradient = ctx.createRadialGradient(0, 0, 0, 0, 0, radius);
          stops.forEach(([offset, color]) => gradient.addColorStop(offset, color));
          ctx.fillStyle = gradient;
          ctx.beginPath();
          ctx.arc(0, 0, radius, 0, Math.PI * 2.0);
          ctx.fill();
        }

        function drawSpiralBand(options) {
          const armCount = options.armCount ?? 4;
          const turns = options.turns ?? 1.15;
          const startRadius = options.startRadius ?? size * 0.10;
          const endRadius = options.endRadius ?? size * 0.46;
          const width = options.width ?? size * 0.055;
          const points = options.points ?? 2600;
          const colorList = options.colors ?? ["rgba(255,255,255,0.12)"];
          const blend = options.blend ?? "screen";
          ctx.save();
          ctx.globalCompositeOperation = blend;
          for (let arm = 0; arm < armCount; arm += 1) {
            const armOffset = arm * (Math.PI * 2.0 / armCount) + (options.phase ?? 0.0);
            for (let i = 0; i < points; i += 1) {
              const t = i / points;
              const radius = startRadius + (endRadius - startRadius) * Math.pow(t, options.radialPower ?? 1.0);
              const theta = armOffset + turns * Math.PI * 2.0 * t + (Math.random() - 0.5) * (options.angularJitter ?? 0.32);
              const spread = width * (0.12 + 0.88 * t) * (0.5 + Math.random());
              const x = Math.cos(theta) * radius + (Math.random() - 0.5) * spread;
              const y = Math.sin(theta) * radius + (Math.random() - 0.5) * spread;
              const dot = (options.minDot ?? 0.35) + Math.random() * (options.maxDot ?? 1.65);
              ctx.fillStyle = colorList[(arm + i) % colorList.length];
              ctx.beginPath();
              ctx.arc(x, y, dot, 0, Math.PI * 2.0);
              ctx.fill();
            }
          }
          ctx.restore();
        }

        if (kind === "stars") {
          drawGlow(size * 0.47, [
            [0.0, "rgba(255, 247, 225, 0.97)"],
            [0.09, "rgba(255, 228, 176, 0.96)"],
            [0.20, "rgba(244, 210, 151, 0.82)"],
            [0.36, "rgba(216, 226, 255, 0.58)"],
            [0.62, "rgba(128, 163, 255, 0.20)"],
            [1.0, "rgba(0, 0, 0, 0.0)"],
          ]);

          const barGradient = ctx.createLinearGradient(-size * 0.24, 0, size * 0.24, 0);
          barGradient.addColorStop(0.0, "rgba(0,0,0,0.0)");
          barGradient.addColorStop(0.20, "rgba(255, 205, 146, 0.16)");
          barGradient.addColorStop(0.50, "rgba(255, 222, 184, 0.40)");
          barGradient.addColorStop(0.80, "rgba(255, 205, 146, 0.16)");
          barGradient.addColorStop(1.0, "rgba(0,0,0,0.0)");
          ctx.save();
          ctx.rotate(0.42);
          ctx.fillStyle = barGradient;
          ctx.beginPath();
          ctx.ellipse(0, 0, size * 0.21, size * 0.070, 0, 0, Math.PI * 2.0);
          ctx.fill();
          ctx.restore();

          drawSpiralBand({
            armCount: 4,
            turns: 0.88,
            startRadius: size * 0.10,
            endRadius: size * 0.46,
            width: size * 0.070,
            points: 2600,
            angularJitter: 0.34,
            minDot: 0.65,
            maxDot: 2.6,
            colors: [
              "rgba(233, 241, 255, 0.16)",
              "rgba(167, 199, 255, 0.20)",
              "rgba(255, 235, 194, 0.18)",
              "rgba(255, 255, 255, 0.12)",
            ],
            blend: "screen",
          });

          drawSpiralBand({
            armCount: 4,
            turns: 0.90,
            startRadius: size * 0.11,
            endRadius: size * 0.46,
            width: size * 0.045,
            points: 1800,
            angularJitter: 0.18,
            minDot: 0.45,
            maxDot: 1.3,
            colors: [
              "rgba(255,255,255,0.18)",
              "rgba(206,224,255,0.16)",
            ],
            blend: "lighter",
          });

          ctx.globalCompositeOperation = "screen";
          for (let i = 0; i < 2400; i += 1) {
            const r = Math.sqrt(Math.random()) * size * 0.47;
            const theta = Math.random() * Math.PI * 2.0;
            const x = Math.cos(theta) * r;
            const y = Math.sin(theta) * r;
            const starSize = 0.30 + Math.random() * 1.45;
            const alpha = (0.025 + Math.random() * 0.12) * (1.0 - 0.55 * (r / (size * 0.47)));
            ctx.fillStyle = i % 5 === 0
              ? `rgba(255, 227, 181, ${alpha})`
              : `rgba(255,255,255,${alpha})`;
            ctx.beginPath();
            ctx.arc(x, y, starSize, 0, Math.PI * 2.0);
            ctx.fill();
          }
          ctx.globalCompositeOperation = "source-over";
        } else if (kind === "dust") {
          const ring = ctx.createRadialGradient(0, 0, size * 0.05, 0, 0, size * 0.48);
          ring.addColorStop(0.0, "rgba(0, 0, 0, 0.0)");
          ring.addColorStop(0.10, "rgba(21, 12, 7, 0.18)");
          ring.addColorStop(0.24, "rgba(35, 20, 12, 0.48)");
          ring.addColorStop(0.52, "rgba(58, 34, 21, 0.34)");
          ring.addColorStop(0.76, "rgba(20, 14, 10, 0.14)");
          ring.addColorStop(1.0, "rgba(0, 0, 0, 0.0)");
          ctx.fillStyle = ring;
          ctx.beginPath();
          ctx.arc(0, 0, size * 0.48, 0, Math.PI * 2.0);
          ctx.fill();

          drawSpiralBand({
            armCount: 4,
            turns: 0.92,
            startRadius: size * 0.15,
            endRadius: size * 0.45,
            width: size * 0.060,
            points: 1600,
            angularJitter: 0.24,
            minDot: 2.4,
            maxDot: 6.8,
            colors: [
              "rgba(31, 18, 12, 0.14)",
              "rgba(58, 34, 21, 0.10)",
            ],
            blend: "multiply",
          });

          ctx.globalCompositeOperation = "destination-out";
          for (let i = 0; i < 820; i += 1) {
            const t = i / 820.0;
            const radius = (0.12 + 0.80 * t) * size * 0.44;
            const theta = 0.34 + (t * Math.PI * 2.2) + (Math.random() - 0.5) * 0.42;
            const x = Math.cos(theta) * radius;
            const y = Math.sin(theta) * radius;
            const blot = 4.0 + Math.random() * 9.0;
            ctx.fillStyle = "rgba(0,0,0,0.06)";
            ctx.beginPath();
            ctx.arc(x, y, blot, 0, Math.PI * 2.0);
            ctx.fill();
          }
          ctx.globalCompositeOperation = "source-over";
        } else if (kind === "halo") {
          drawGlow(size * 0.49, [
            [0.0, "rgba(160, 181, 255, 0.22)"],
            [0.30, "rgba(121, 147, 228, 0.14)"],
            [0.62, "rgba(78, 98, 176, 0.06)"],
            [1.0, "rgba(0, 0, 0, 0.0)"],
          ]);
        }

        const texture = new THREE.CanvasTexture(canvasEl);
        if (THREE.SRGBColorSpace) {
          texture.colorSpace = THREE.SRGBColorSpace;
        } else if (THREE.sRGBEncoding) {
          texture.encoding = THREE.sRGBEncoding;
        }
        texture.needsUpdate = true;
        galaxyTextureCache.set(kind, texture);
        return texture;
      }

      function seededRandom(seed) {
        let state = seed >>> 0;
        return function () {
          state = (1664525 * state + 1013904223) >>> 0;
          return state / 4294967296;
        };
      }

      function gaussianFromRandom(randFn) {
        let u = 0.0;
        let v = 0.0;
        while (u <= 1e-7) {
          u = randFn();
        }
        while (v <= 1e-7) {
          v = randFn();
        }
        return Math.sqrt(-2.0 * Math.log(u)) * Math.cos(2.0 * Math.PI * v);
      }

      function colorAttrArray(count, baseColorHex, jitter, randFn) {
        const arr = new Float32Array(count * 3);
        const base = new THREE.Color(baseColorHex);
        const hsl = { h: 0, s: 0, l: 0 };
        base.getHSL(hsl);
        for (let i = 0; i < count; i += 1) {
          const col = new THREE.Color();
          const dh = (randFn() - 0.5) * (jitter.h ?? 0.0);
          const ds = (randFn() - 0.5) * (jitter.s ?? 0.0);
          const dl = (randFn() - 0.5) * (jitter.l ?? 0.0);
          col.setHSL(
            THREE.MathUtils.clamp(hsl.h + dh, 0.0, 1.0),
            THREE.MathUtils.clamp(hsl.s + ds, 0.0, 1.0),
            THREE.MathUtils.clamp(hsl.l + dl, 0.0, 1.0)
          );
          arr[3 * i] = col.r;
          arr[3 * i + 1] = col.g;
          arr[3 * i + 2] = col.b;
        }
        return arr;
      }

      function makeCircularSpriteTexture(colorStops) {
        const canvasEl = document.createElement("canvas");
        canvasEl.width = 128;
        canvasEl.height = 128;
        const ctx = canvasEl.getContext("2d");
        const gradient = ctx.createRadialGradient(64, 64, 0, 64, 64, 64);
        colorStops.forEach(([offset, color]) => gradient.addColorStop(offset, color));
        ctx.fillStyle = gradient;
        ctx.beginPath();
        ctx.arc(64, 64, 64, 0, Math.PI * 2.0);
        ctx.fill();
        const texture = new THREE.CanvasTexture(canvasEl);
        if (THREE.SRGBColorSpace) {
          texture.colorSpace = THREE.SRGBColorSpace;
        } else if (THREE.sRGBEncoding) {
          texture.encoding = THREE.sRGBEncoding;
        }
        texture.needsUpdate = true;
        return texture;
      }

      function createGalaxyPointCloud(params) {
        const rand = seededRandom(params.seed ?? 1);
        const count = params.count ?? 4000;
        const positions = new Float32Array(count * 3);
        const sizes = new Float32Array(count);
        const colors = params.vertexColors
          ? colorAttrArray(count, params.baseColor ?? "#ffffff", params.colorJitter ?? {}, rand)
          : null;
        const armCount = params.armCount ?? 4;
        const radialScale = params.radius ?? 8000.0;
        const innerRadius = params.innerRadius ?? 0.0;
        const turns = params.turns ?? 0.9;
        const armTightness = params.armTightness ?? 1.0;
        const armWidth = params.armWidth ?? (radialScale * 0.08);
        const verticalScale = params.verticalScale ?? 100.0;
        const phase = params.phase ?? 0.0;
        const flatten = params.flatten ?? 1.0;
        const distribution = params.distribution ?? "spiral";

        for (let i = 0; i < count; i += 1) {
          const radialRand = rand();
          const r = innerRadius + (radialScale - innerRadius) * Math.pow(radialRand, params.radialPower ?? 1.6);
          let x = 0.0;
          let y = 0.0;
          if (distribution === "spiral") {
            const armIndex = Math.floor(rand() * armCount);
            const baseTheta = phase + armIndex * (Math.PI * 2.0 / armCount);
            const theta = baseTheta + turns * Math.PI * 2.0 * Math.pow(r / radialScale, armTightness);
            const tangentJitter = gaussianFromRandom(rand) * (params.angularScatter ?? 0.05);
            const radialJitter = gaussianFromRandom(rand) * armWidth * (0.25 + 0.75 * (r / radialScale));
            const rr = Math.max(innerRadius * 0.5, r + radialJitter);
            x = Math.cos(theta + tangentJitter) * rr;
            y = Math.sin(theta + tangentJitter) * rr;
          } else {
            const theta = rand() * Math.PI * 2.0;
            const rr = innerRadius + (radialScale - innerRadius) * Math.sqrt(rand());
            x = Math.cos(theta) * rr;
            y = Math.sin(theta) * rr;
          }
          const z = gaussianFromRandom(rand) * verticalScale;
          positions[3 * i] = x;
          positions[3 * i + 1] = y * flatten;
          positions[3 * i + 2] = z;
          sizes[i] = (params.minSize ?? 14.0) + rand() * (params.maxSize ?? 30.0);
        }

        const geometry = new THREE.BufferGeometry();
        geometry.setAttribute("position", new THREE.BufferAttribute(positions, 3));
        geometry.setAttribute("size", new THREE.BufferAttribute(sizes, 1));
        if (colors) {
          geometry.setAttribute("color", new THREE.BufferAttribute(colors, 3));
        }

        const spriteTexture = makeCircularSpriteTexture(
          params.textureStops ?? [
            [0.0, "rgba(255,255,255,1.0)"],
            [0.35, "rgba(255,255,255,0.96)"],
            [0.75, "rgba(255,255,255,0.20)"],
            [1.0, "rgba(255,255,255,0.0)"],
          ]
        );
        const material = new THREE.ShaderMaterial({
          uniforms: {
            pointTexture: { value: spriteTexture },
            opacity: { value: params.opacity ?? 1.0 },
            baseColor: { value: new THREE.Color(params.baseColor ?? "#ffffff") },
          },
          vertexShader: `
            attribute float size;
            varying vec3 vColor;
            void main() {
              vColor = color;
              vec4 mvPosition = modelViewMatrix * vec4(position, 1.0);
              gl_PointSize = size * (280.0 / max(-mvPosition.z, 1.0));
              gl_Position = projectionMatrix * mvPosition;
            }
          `,
          fragmentShader: `
            uniform sampler2D pointTexture;
            uniform float opacity;
            uniform vec3 baseColor;
            varying vec3 vColor;
            void main() {
              vec4 tex = texture2D(pointTexture, gl_PointCoord);
              vec3 color = ${params.vertexColors ? "vColor" : "baseColor"};
              gl_FragColor = vec4(color, tex.a * opacity);
              if (gl_FragColor.a < 0.02) discard;
            }
          `,
          transparent: true,
          depthWrite: false,
          depthTest: true,
          blending: params.blending === "additive" ? THREE.AdditiveBlending : THREE.NormalBlending,
          vertexColors: Boolean(params.vertexColors),
        });

        const points = new THREE.Points(geometry, material);
        points.renderOrder = params.renderOrder ?? -8;
        if (params.rotationZ) {
          points.rotation.z = params.rotationZ;
        }
        return points;
      }

      function createMilkyWayModel(decoration) {
        const center = decoration.center || sceneSpec.center || { x: 0, y: 0, z: 0 };
        const discRadius = decoration.disc_radius_pc ?? 10000.0;
        const discThickness = decoration.disc_thickness_pc ?? 200.0;
        const bulgeRadius = decoration.bulge_radius_pc ?? 1800.0;
        const bulgeHeight = decoration.bulge_height_pc ?? 1200.0;
        const haloRadius = decoration.halo_radius_pc ?? (discRadius * 1.15);
        const dustInner = decoration.dust_inner_radius_pc ?? (discRadius * 0.14);
        const dustOuter = decoration.dust_outer_radius_pc ?? (discRadius * 0.92);
        const modelOrientationSign = -1.0;

        const group = new THREE.Group();
        group.position.set(center.x, center.y, center.z);

        const ambient = new THREE.AmbientLight(0x8593a8, 0.55);
        const hemi = new THREE.HemisphereLight(0xbfd0ff, 0x09070c, 0.62);
        const dir = new THREE.DirectionalLight(0xd7e3ff, 0.82);
        dir.position.set(discRadius * 0.35, modelOrientationSign * -discRadius * 0.65, bulgeHeight * 2.2);
        const coreLight = new THREE.PointLight(0xffddb0, 1.7, discRadius * 2.8, 2.0);
        coreLight.position.set(0.0, 0.0, bulgeHeight * 0.28);
        group.add(ambient, hemi, dir, coreLight);

        const stellarCloud = createGalaxyPointCloud({
          seed: 1087,
          count: 5400,
          radius: discRadius,
          innerRadius: 0.08 * discRadius,
          turns: 0.92,
          armCount: 4,
          armWidth: 0.075 * discRadius,
          verticalScale: 0.40 * discThickness,
          minSize: 10.0,
          maxSize: 21.0,
          baseColor: "#dbe7ff",
          vertexColors: true,
          colorJitter: { h: 0.04, s: 0.10, l: 0.18 },
          opacity: 0.56,
          rotationZ: 0.06,
          blending: "additive",
          renderOrder: -7,
        });
        stellarCloud.scale.y = modelOrientationSign;
        group.add(stellarCloud);

        const warmArmCloud = createGalaxyPointCloud({
          seed: 2081,
          count: 2400,
          radius: 0.96 * discRadius,
          innerRadius: 0.14 * discRadius,
          turns: 0.88,
          armCount: 4,
          armWidth: 0.05 * discRadius,
          verticalScale: 0.24 * discThickness,
          minSize: 12.0,
          maxSize: 24.0,
          baseColor: "#ffe1b8",
          vertexColors: true,
          colorJitter: { h: 0.03, s: 0.12, l: 0.16 },
          opacity: 0.40,
          rotationZ: 0.13,
          blending: "additive",
          renderOrder: -6,
        });
        warmArmCloud.scale.y = modelOrientationSign;
        group.add(warmArmCloud);

        const thickDiscCloud = createGalaxyPointCloud({
          seed: 3011,
          count: 2600,
          radius: 0.96 * discRadius,
          innerRadius: 0.12 * discRadius,
          distribution: "disc",
          verticalScale: 1.55 * discThickness,
          minSize: 7.0,
          maxSize: 16.0,
          baseColor: "#c2d1ff",
          vertexColors: true,
          colorJitter: { h: 0.03, s: 0.05, l: 0.11 },
          opacity: 0.14,
          renderOrder: -9,
        });
        group.add(thickDiscCloud);

        const planeSize = haloRadius * 2.0;
        const planeGeometry = new THREE.PlaneGeometry(planeSize, planeSize, 1, 1);
        const stellarTexture = galaxyTextureFor("stars");
        const haloTexture = galaxyTextureFor("halo");
        const dustTexture = galaxyTextureFor("dust");

        const stellarLayers = [
          { z: -0.85 * discThickness, opacity: 0.11, emissive: 0.28 },
          { z: -0.42 * discThickness, opacity: 0.18, emissive: 0.40 },
          { z: 0.0, opacity: 0.34, emissive: 0.55 },
          { z: 0.42 * discThickness, opacity: 0.18, emissive: 0.40 },
          { z: 0.85 * discThickness, opacity: 0.11, emissive: 0.28 },
        ];

        stellarLayers.forEach((layer) => {
          const material = new THREE.MeshStandardMaterial({
            map: stellarTexture,
            transparent: true,
            opacity: layer.opacity,
            alphaTest: 0.02,
            depthWrite: false,
            side: THREE.DoubleSide,
            color: 0xe7ecff,
            emissive: 0x8aa6ff,
            emissiveMap: stellarTexture,
            emissiveIntensity: layer.emissive,
            roughness: 0.92,
            metalness: 0.0,
          });
          const mesh = new THREE.Mesh(planeGeometry, material);
          mesh.position.z = layer.z;
          mesh.scale.y = modelOrientationSign;
          mesh.renderOrder = -12;
          group.add(mesh);
        });

        const dustGeometry = new THREE.RingGeometry(dustInner, dustOuter, 160);
        const dustMaterial = new THREE.MeshStandardMaterial({
          map: dustTexture,
          transparent: true,
          opacity: 0.30,
          alphaTest: 0.03,
          depthWrite: false,
          side: THREE.DoubleSide,
          color: 0x2d1f18,
          emissive: 0x110907,
          emissiveMap: dustTexture,
          emissiveIntensity: 0.10,
          roughness: 1.0,
          metalness: 0.0,
        });
        const dustTop = new THREE.Mesh(dustGeometry, dustMaterial);
        dustTop.position.z = 0.18 * discThickness;
        dustTop.scale.y = modelOrientationSign;
        dustTop.renderOrder = -10;
        group.add(dustTop);

        const dustBottom = new THREE.Mesh(dustGeometry.clone(), dustMaterial.clone());
        dustBottom.position.z = -0.18 * discThickness;
        dustBottom.scale.y = modelOrientationSign;
        dustBottom.renderOrder = -10;
        group.add(dustBottom);

        const haloMaterial = new THREE.MeshBasicMaterial({
          map: haloTexture,
          transparent: true,
          opacity: 0.22,
          alphaTest: 0.01,
          depthWrite: false,
          side: THREE.DoubleSide,
          color: 0x9ab3ff,
        });
        const haloMesh = new THREE.Mesh(planeGeometry.clone(), haloMaterial);
        haloMesh.position.z = 0.0;
        haloMesh.scale.y = modelOrientationSign;
        haloMesh.renderOrder = -14;
        group.add(haloMesh);

        const outerGlow = new THREE.Mesh(
          planeGeometry.clone(),
          new THREE.MeshBasicMaterial({
            map: haloTexture,
            transparent: true,
            opacity: 0.16,
            alphaTest: 0.01,
            depthWrite: false,
            side: THREE.DoubleSide,
            color: 0xc3d4ff,
          })
        );
        outerGlow.scale.set(1.18, 1.18, 1.0);
        outerGlow.position.z = 0.0;
        outerGlow.scale.y *= modelOrientationSign;
        outerGlow.renderOrder = -15;
        group.add(outerGlow);

        const bar = new THREE.Mesh(
          new THREE.BoxGeometry(0.58 * bulgeRadius, 2.25 * bulgeRadius, 0.38 * bulgeHeight),
          new THREE.MeshPhysicalMaterial({
            color: 0xf1d7a7,
            emissive: 0x7a5330,
            emissiveIntensity: 0.24,
            roughness: 0.72,
            metalness: 0.0,
            clearcoat: 0.06,
            transparent: true,
            opacity: 0.74,
          })
        );
        bar.rotation.z = 0.46;
        bar.rotation.z *= modelOrientationSign;
        bar.position.z = 0.0;
        group.add(bar);

        const bulgeMaterial = new THREE.MeshPhysicalMaterial({
          color: 0xf3dfb7,
          emissive: 0x8c6845,
          emissiveIntensity: 0.44,
          roughness: 0.68,
          metalness: 0.0,
          clearcoat: 0.10,
          transparent: true,
          opacity: 0.90,
        });
        const bulge = new THREE.Mesh(
          new THREE.SphereGeometry(bulgeRadius, 56, 40),
          bulgeMaterial
        );
        bulge.scale.set(1.22, 1.22, Math.max(bulgeHeight / Math.max(bulgeRadius, 1.0), 0.65));
        bulge.position.z = 0.0;
        group.add(bulge);

        const bulgeCloud = createGalaxyPointCloud({
          seed: 4109,
          count: 1800,
          radius: 0.65 * bulgeRadius,
          innerRadius: 0.0,
          distribution: "disc",
          verticalScale: 0.55 * bulgeHeight,
          minSize: 10.0,
          maxSize: 20.0,
          baseColor: "#ffe4bd",
          vertexColors: true,
          colorJitter: { h: 0.02, s: 0.08, l: 0.12 },
          opacity: 0.42,
          renderOrder: -5,
        });
        bulgeCloud.scale.set(1.15, 1.15, 0.82);
        group.add(bulgeCloud);

        const core = new THREE.Mesh(
          new THREE.SphereGeometry(0.28 * bulgeRadius, 36, 24),
          new THREE.MeshPhysicalMaterial({
            color: 0xffe1ae,
            emissive: 0xffc06c,
            emissiveIntensity: 1.10,
            roughness: 0.35,
            metalness: 0.0,
            transparent: true,
            opacity: 0.92,
          })
        );
        core.position.z = 0.0;
        group.add(core);

        const innerGlow = new THREE.Mesh(
          new THREE.SphereGeometry(0.60 * bulgeRadius, 28, 22),
          new THREE.MeshBasicMaterial({
            color: 0xffc88a,
            transparent: true,
            opacity: 0.12,
            depthWrite: false,
          })
        );
        innerGlow.scale.set(1.0, 1.0, 0.68);
        innerGlow.position.z = 0.0;
        group.add(innerGlow);

        return group;
      }

      function makeLineObject(trace, materialBucket) {
        if (!trace.segments || !trace.segments.length) {
          return null;
        }
        const positions = [];
        trace.segments.forEach((segment) => {
          positions.push(segment[0], segment[1], segment[2], segment[3], segment[4], segment[5]);
        });
        const geometry = new LineSegmentsGeometry();
        geometry.setPositions(positions);
        const material = new LineMaterial({
          color: (trace.line || {}).color ?? "#ffffff",
          linewidth: Math.max((trace.line || {}).width ?? 1.0, 1.0),
          dashed: ((trace.line || {}).dash || "solid") !== "solid",
          dashSize: 8.0,
          gapSize: 5.0,
          transparent: (trace.opacity ?? 1.0) < 1.0,
          opacity: trace.opacity ?? 1.0,
          worldUnits: false,
        });
        material.resolution.set(root.clientWidth, root.clientHeight);
        const line = new LineSegments2(geometry, material);
        line.computeLineDistances();
        materialBucket.push(material);
        return line;
      }

      function addMarkerTrace(parent, trace) {
        if (!trace.points || !trace.points.length) {
          return;
        }
        const group = new THREE.Group();
        trace.points.forEach((point) => {
          if (!Number.isFinite(point.size) || point.size <= 0) {
            return;
          }
          const sprite = new THREE.Sprite(markerMaterialFor(point.symbol, point.color, point.opacity));
          const scale = Math.max(point.size * pointScale, pointScale * 0.5);
          sprite.position.set(point.x, point.y, point.z);
          sprite.scale.set(scale, scale, scale);
          sprite.userData = {
            hovertext: point.hovertext || trace.name || "",
            selection: point.selection || null,
          };
          group.add(sprite);
          hoverTargets.push(sprite);
        });
        parent.add(group);
      }

      function addTextTrace(parent, trace) {
        if (!trace.labels || !trace.labels.length) {
          return;
        }
        const group = new THREE.Group();
        trace.labels.forEach((label) => {
          if (!label.text) {
            return;
          }
          const sprite = makeTextSprite(label.text, {
            color: label.color ?? theme.axis_color,
            size: label.size ?? 12,
            family: label.family ?? "Helvetica",
          });
          sprite.position.set(label.x, label.y, label.z);
          group.add(sprite);
        });
        parent.add(group);
      }

      function addDecoration(parent, decoration) {
        if (!decoration || !decoration.kind) {
          return;
        }
        if (decoration.kind === "milky_way_model") {
          parent.add(createMilkyWayModel(decoration));
        }
      }

      function buildSelectionFootprint(selection, materialBucket) {
        if (!selection) {
          return null;
        }
        const axis = [
          Number(selection.x0),
          Number(selection.y0),
          Number(selection.z0),
        ];
        if (!axis.every(Number.isFinite)) {
          return null;
        }
        const axisLen = Math.hypot(axis[0], axis[1], axis[2]);
        if (!(axisLen > 0.0)) {
          return null;
        }
        const halfAngleDeg = Number(skySpec.radius_deg || 1.0);
        const halfAngleRad = halfAngleDeg * Math.PI / 180.0;
        if (!(halfAngleRad > 0.0)) {
          return null;
        }

        const wHat = new THREE.Vector3(axis[0], axis[1], axis[2]).normalize();
        let ref = new THREE.Vector3(0, 0, 1);
        if (Math.abs(wHat.dot(ref)) > 0.95) {
          ref = new THREE.Vector3(0, 1, 0);
        }
        const uHat = new THREE.Vector3().crossVectors(wHat, ref).normalize();
        const vHat = new THREE.Vector3().crossVectors(wHat, uHat).normalize();
        const nTheta = 72;
        const nSteps = 24;
        const tanHalf = Math.tan(halfAngleRad);
        const positions = [];
        const index = [];
        const rimSegments = [];
        const rimPoints = [];

        for (let stepI = 0; stepI < nSteps; stepI += 1) {
          const s = axisLen * (stepI / (nSteps - 1));
          const ringRadius = s * tanHalf;
          for (let thetaI = 0; thetaI < nTheta; thetaI += 1) {
            const ang = 2.0 * Math.PI * thetaI / nTheta;
            const offset = uHat.clone().multiplyScalar(Math.cos(ang) * ringRadius)
              .add(vHat.clone().multiplyScalar(Math.sin(ang) * ringRadius));
            const point = wHat.clone().multiplyScalar(s).add(offset);
            positions.push(point.x, point.y, point.z);
            if (stepI === nSteps - 1) {
              rimPoints.push(point);
            }
          }
        }

        for (let stepI = 0; stepI < nSteps - 1; stepI += 1) {
          for (let thetaI = 0; thetaI < nTheta; thetaI += 1) {
            const nextTheta = (thetaI + 1) % nTheta;
            const base0 = stepI * nTheta + thetaI;
            const base1 = stepI * nTheta + nextTheta;
            const top0 = (stepI + 1) * nTheta + thetaI;
            const top1 = (stepI + 1) * nTheta + nextTheta;
            index.push(base0, top0, top1);
            index.push(base0, top1, base1);
          }
        }

        if (rimPoints.length >= 2) {
          for (let i = 0; i < rimPoints.length; i += 1) {
            const a = rimPoints[i];
            const b = rimPoints[(i + 1) % rimPoints.length];
            rimSegments.push([a.x, a.y, a.z, b.x, b.y, b.z]);
          }
        }

        const group = new THREE.Group();
        const geometry = new THREE.BufferGeometry();
        geometry.setAttribute("position", new THREE.Float32BufferAttribute(positions, 3));
        geometry.setIndex(index);
        geometry.computeVertexNormals();
        const mesh = new THREE.Mesh(
          geometry,
          new THREE.MeshBasicMaterial({
            color: theme.footprint || "#6ec5ff",
            transparent: true,
            opacity: 0.18,
            side: THREE.DoubleSide,
            depthWrite: false,
          })
        );
        group.add(mesh);

        const rim = makeLineObject({
          segments: rimSegments,
          line: {
            color: theme.footprint || "#6ec5ff",
            width: 5.0,
            dash: "solid",
          },
          opacity: 1.0,
        }, materialBucket);
        if (rim) {
          group.add(rim);
        }
        return group;
      }

      function buildAxes() {
        clearGroup(axisGroup);
        axisLineMaterials.length = 0;
        if (!sceneSpec.show_axes) {
          return;
        }

        const xAxis = axisSpec.x || {};
        const yAxis = axisSpec.y || {};
        const zAxis = axisSpec.z || {};
        const xRange = sceneSpec.ranges.x;
        const yRange = sceneSpec.ranges.y;
        const zRange = sceneSpec.ranges.z;
        const tickLength = (sceneSpec.max_span || 1) * 0.02;
        const axisColor = xAxis.linecolor ?? theme.axis_color;
        const axisWidth = xAxis.linewidth ?? 2.0;
        const x0 = xRange[0];
        const y0 = yRange[0];
        const z0 = zRange[0];

        const axisSegments = [
          [xRange[0], y0, z0, xRange[1], y0, z0],
          [x0, yRange[0], z0, x0, yRange[1], z0],
          [x0, y0, zRange[0], x0, y0, zRange[1]],
        ];
        const axisLines = makeLineObject({
          segments: axisSegments,
          line: { color: axisColor, width: axisWidth, dash: "solid" },
          opacity: 1.0,
        }, axisLineMaterials);
        if (axisLines) {
          axisGroup.add(axisLines);
        }

        [
          ["x", xAxis, xRange],
          ["y", yAxis, yRange],
          ["z", zAxis, zRange],
        ].forEach(([axisName, axis, range]) => {
          const nTicks = Math.max(Number(axis.nticks ?? 5), 2);
          const step = (range[1] - range[0]) / (nTicks - 1);
          const tickSegments = [];
          for (let i = 0; i < nTicks; i += 1) {
            const value = range[0] + step * i;
            if (axisName === "x") {
              tickSegments.push([value, y0, z0, value, y0 - tickLength, z0]);
            } else if (axisName === "y") {
              tickSegments.push([x0, value, z0, x0 - tickLength, value, z0]);
            } else {
              tickSegments.push([x0, y0, value, x0 - tickLength, y0, value]);
            }

            const label = makeTextSprite(formatTick(value), {
              color: (axis.tickfont || {}).color ?? theme.axis_color,
              size: (axis.tickfont || {}).size || 14,
              family: (axis.tickfont || {}).family ?? "Helvetica",
            });
            if (axisName === "x") {
              label.position.set(value, y0 - tickLength * 2.4, z0);
            } else if (axisName === "y") {
              label.position.set(x0 - tickLength * 2.6, value, z0);
            } else {
              label.position.set(x0 - tickLength * 2.6, y0, value);
            }
            axisGroup.add(label);
          }

          const tickLines = makeLineObject({
            segments: tickSegments,
            line: { color: axisColor, width: Math.max(axisWidth * 0.75, 1.0), dash: "solid" },
            opacity: 1.0,
          }, axisLineMaterials);
          if (tickLines) {
            axisGroup.add(tickLines);
          }

          const titleText = typeof axis.title === "string" ? axis.title : (axis.title || {}).text;
          if (titleText) {
            const titleSprite = makeTextSprite(titleText, {
              color: (axis.title_font || axis.titlefont || {}).color ?? theme.axis_color,
              size: (axis.title_font || axis.titlefont || {}).size || 18,
              family: (axis.title_font || axis.titlefont || {}).family ?? "Helvetica",
            });
            if (axisName === "x") {
              titleSprite.position.set(range[1], y0 - tickLength * 4.5, z0);
            } else if (axisName === "y") {
              titleSprite.position.set(x0 - tickLength * 4.8, range[1], z0);
            } else {
              titleSprite.position.set(x0 - tickLength * 4.8, y0, range[1]);
            }
            axisGroup.add(titleSprite);
          }
        });
      }

      function groupDefaults(groupName) {
        return groupVisibility[groupName] || {};
      }

      function resetLegendState(groupName) {
        legendState = {};
        const defaults = groupDefaults(groupName);
        legendItems.forEach((item) => {
          const mode = defaults[item.key];
          if (mode === true) {
            legendState[item.key] = true;
          } else if (mode === "legendonly") {
            legendState[item.key] = false;
          }
        });
      }

      function traceVisible(trace) {
        const defaults = groupDefaults(currentGroup);
        const mode = defaults[trace.key];
        if (mode === false || mode === undefined) {
          return false;
        }
        if (trace.showlegend) {
          return Boolean(legendState[trace.key]);
        }
        return mode === true;
      }

      function renderLegend() {
        legendEl.innerHTML = "";
        const title = document.createElement("div");
        title.className = "oviz-three-legend-title";
        title.textContent = "Click to toggle traces on/off";
        legendEl.appendChild(title);

        const defaults = groupDefaults(currentGroup);
        legendItems.forEach((item) => {
          const mode = defaults[item.key];
          if (mode === false || mode === undefined) {
            return;
          }

          const label = document.createElement("label");
          label.className = "oviz-three-legend-item";
          const checkbox = document.createElement("input");
          checkbox.type = "checkbox";
          checkbox.checked = Boolean(legendState[item.key]);
          checkbox.addEventListener("change", () => {
            legendState[item.key] = checkbox.checked;
            renderFrame(currentFrameIndex);
          });
          const text = document.createElement("span");
          text.textContent = item.name;
          label.appendChild(checkbox);
          label.appendChild(text);
          legendEl.appendChild(label);
        });
      }

      function renderFrame(index) {
        currentFrameIndex = Math.max(0, Math.min(index, frameSpecs.length - 1));
        sliderEl.value = String(currentFrameIndex);
        const frame = frameSpecs[currentFrameIndex];
        timeLabelEl.textContent = `Time (Myr): ${frame.name}`;
        tooltipEl.style.display = "none";
        hoverTargets.length = 0;
        clearGroup(plotGroup);
        frameLineMaterials.length = 0;

        frame.traces.forEach((trace) => {
          if (!traceVisible(trace)) {
            return;
          }

          if (trace.segments && trace.segments.length) {
            const line = makeLineObject(trace, frameLineMaterials);
            if (line) {
              plotGroup.add(line);
            }
          }
          if (trace.points && trace.points.length) {
            addMarkerTrace(plotGroup, trace);
          }
          if (trace.labels && trace.labels.length) {
            addTextTrace(plotGroup, trace);
          }
        });

        if (currentSelection && approximatelyZero(Number(frame.time))) {
          const footprint = buildSelectionFootprint(currentSelection, frameLineMaterials);
          if (footprint) {
            plotGroup.add(footprint);
          }
        }

        (frame.decorations || []).forEach((decoration) => {
          addDecoration(plotGroup, decoration);
        });

        resize();
      }

      function resize() {
        const width = root.clientWidth;
        const height = root.clientHeight;
        renderer.setSize(width, height, false);
        camera.aspect = width / Math.max(height, 1);
        camera.updateProjectionMatrix();
        [...axisLineMaterials, ...frameLineMaterials].forEach((material) => {
          material.resolution.set(width, height);
        });
      }

      function play() {
        if (playbackTimer || frameSpecs.length <= 1) {
          return;
        }
        playbackTimer = window.setInterval(() => {
          const nextIndex = (currentFrameIndex + 1) % frameSpecs.length;
          renderFrame(nextIndex);
        }, sceneSpec.playback_interval_ms || 500);
      }

      function pause() {
        if (playbackTimer) {
          window.clearInterval(playbackTimer);
          playbackTimer = null;
        }
      }

      function storeSkyPanelRect() {
        const rect = skyPanelEl.getBoundingClientRect();
        skyPanelEl.dataset.normalLeft = String(rect.left);
        skyPanelEl.dataset.normalTop = String(rect.top);
        skyPanelEl.dataset.normalWidth = String(rect.width);
        skyPanelEl.dataset.normalHeight = String(rect.height);
      }

      function restoreSkyPanelRect() {
        const left = Number(skyPanelEl.dataset.normalLeft);
        const top = Number(skyPanelEl.dataset.normalTop);
        const width = Number(skyPanelEl.dataset.normalWidth);
        const height = Number(skyPanelEl.dataset.normalHeight);
        if ([left, top, width, height].every(Number.isFinite)) {
          skyPanelEl.style.left = `${left}px`;
          skyPanelEl.style.top = `${top}px`;
          skyPanelEl.style.right = "auto";
          skyPanelEl.style.bottom = "auto";
          skyPanelEl.style.width = `${width}px`;
          skyPanelEl.style.height = `${height}px`;
          return;
        }
        skyPanelEl.style.left = "";
        skyPanelEl.style.top = "";
        skyPanelEl.style.right = "12px";
        skyPanelEl.style.bottom = "";
        skyPanelEl.style.width = "";
        skyPanelEl.style.height = "";
      }

      function setSkyPanelMode(mode) {
        if (!skySpec.enabled) {
          return;
        }
        const nextMode = ["normal", "fullscreen", "hidden"].includes(mode) ? mode : "normal";
        if (nextMode === "fullscreen" && skyPanelMode === "normal") {
          storeSkyPanelRect();
        }
        skyPanelMode = nextMode;
        applySkyPanelMode();
        if (skyPanelMode === "normal") {
          restoreSkyPanelRect();
        } else if (skyPanelMode === "fullscreen") {
          skyPanelEl.style.left = "0px";
          skyPanelEl.style.top = "0px";
          skyPanelEl.style.right = "0px";
          skyPanelEl.style.bottom = "0px";
          skyPanelEl.style.width = "auto";
          skyPanelEl.style.height = "auto";
        }
        resize();
      }

      function pickSprite(event) {
        const rect = canvas.getBoundingClientRect();
        pointer.x = ((event.clientX - rect.left) / rect.width) * 2.0 - 1.0;
        pointer.y = -((event.clientY - rect.top) / rect.height) * 2.0 + 1.0;
        raycaster.setFromCamera(pointer, camera);
        const hits = raycaster.intersectObjects(hoverTargets, false);
        return hits.length ? hits[0].object : null;
      }

      function onCanvasClick(event) {
        if (!skySpec.enabled || skyPointerState) {
          return;
        }
        const hit = pickSprite(event);
        const selection = hit && hit.userData ? hit.userData.selection : null;
        if (!selection) {
          return;
        }
        const frame = currentFrame();
        if (!frame || !approximatelyZero(Number(selection.click_time_myr)) || !approximatelyZero(Number(frame.time))) {
          return;
        }
        currentSelection = selection;
        updateSkyPanel();
        renderFrame(currentFrameIndex);
      }

      function onSkyPointerStart(event) {
        if (skyPanelMode !== "normal") {
          return;
        }
        const resizeHandle = event.target.closest(".oviz-three-sky-resize");
        const dragHandle = event.target.closest(".oviz-three-sky-drag");
        if (!resizeHandle && !dragHandle) {
          return;
        }
        const rect = skyPanelEl.getBoundingClientRect();
        skyPointerState = {
          mode: resizeHandle ? "resize" : "drag",
          dir: resizeHandle ? (resizeHandle.dataset.dir || "se").toLowerCase() : null,
          startX: event.clientX,
          startY: event.clientY,
          startLeft: rect.left,
          startTop: rect.top,
          startWidth: rect.width,
          startHeight: rect.height,
          startRight: rect.right,
          startBottom: rect.bottom,
          handle: resizeHandle || dragHandle,
        };
        skyPanelEl.style.left = `${rect.left}px`;
        skyPanelEl.style.top = `${rect.top}px`;
        skyPanelEl.style.right = "auto";
        skyPanelEl.style.bottom = "auto";
        skyPanelEl.style.width = `${rect.width}px`;
        skyPanelEl.style.height = `${rect.height}px`;
        controls.enabled = false;
        if (skyPointerState.handle && typeof skyPointerState.handle.setPointerCapture === "function" && event.pointerId !== undefined) {
          try {
            skyPointerState.handle.setPointerCapture(event.pointerId);
          } catch (_err) {
          }
        }
        document.body.style.userSelect = "none";
        event.preventDefault();
        event.stopPropagation();
      }

      function onSkyPointerMove(event) {
        if (!skyPointerState) {
          return;
        }
        if (skyPointerState.mode === "drag") {
          const left = skyPointerState.startLeft + (event.clientX - skyPointerState.startX);
          const top = skyPointerState.startTop + (event.clientY - skyPointerState.startY);
          const next = clampSkyPosition(left, top, skyPointerState.startWidth, skyPointerState.startHeight);
          skyPanelEl.style.left = `${next.left}px`;
          skyPanelEl.style.top = `${next.top}px`;
        } else {
          const next = resizeSkyRect(skyPointerState, event.clientX, event.clientY);
          skyPanelEl.style.left = `${next.left}px`;
          skyPanelEl.style.top = `${next.top}px`;
          skyPanelEl.style.width = `${next.width}px`;
          skyPanelEl.style.height = `${next.height}px`;
        }
        resize();
        event.preventDefault();
        event.stopPropagation();
      }

      function onSkyPointerEnd(event) {
        if (!skyPointerState) {
          return;
        }
        if (skyPointerState.handle && typeof skyPointerState.handle.releasePointerCapture === "function" && event.pointerId !== undefined) {
          try {
            skyPointerState.handle.releasePointerCapture(event.pointerId);
          } catch (_err) {
          }
        }
        controls.enabled = true;
        document.body.style.userSelect = "";
        storeSkyPanelRect();
        skyPointerState = null;
      }

      function initSkyPanel() {
        if (!skySpec.enabled) {
          applySkyPanelMode();
          return;
        }
        skyReadoutEl.textContent = skyReadoutText(null);
        skyFrameEl.srcdoc = buildEmptySkySrcdoc();
        skyShowButtonEl.addEventListener("click", () => setSkyPanelMode("normal"));
        skyShowFullButtonEl.addEventListener("click", () => setSkyPanelMode("fullscreen"));
        skyHideButtonEl.addEventListener("click", () => setSkyPanelMode("hidden"));
        skyFullButtonEl.addEventListener("click", () => {
          setSkyPanelMode(skyPanelMode === "fullscreen" ? "normal" : "fullscreen");
        });
        skyDragHandleEl.addEventListener("pointerdown", onSkyPointerStart);
        skyResizeEls.forEach((handle) => handle.addEventListener("pointerdown", onSkyPointerStart));
        applySkyPanelMode();
      }

      function onPointerMove(event) {
        const rect = canvas.getBoundingClientRect();
        pointer.x = ((event.clientX - rect.left) / rect.width) * 2.0 - 1.0;
        pointer.y = -((event.clientY - rect.top) / rect.height) * 2.0 + 1.0;
        raycaster.setFromCamera(pointer, camera);
        const hits = raycaster.intersectObjects(hoverTargets, false);
        if (!hits.length || !hits[0].object.userData.hovertext) {
          tooltipEl.style.display = "none";
          return;
        }
        tooltipEl.style.display = "block";
        tooltipEl.innerHTML = hits[0].object.userData.hovertext;
        tooltipEl.style.left = `${event.clientX - rect.left + 14}px`;
        tooltipEl.style.top = `${event.clientY - rect.top + 14}px`;
      }

      function onPointerLeave() {
        tooltipEl.style.display = "none";
      }

      function initControls() {
        groupSelectEl.innerHTML = "";
        const groups = sceneSpec.group_order || [defaultGroup];
        groups.forEach((groupName) => {
          const option = document.createElement("option");
          option.value = groupName;
          option.textContent = groupName;
          groupSelectEl.appendChild(option);
        });
        groupSelectEl.value = currentGroup;
        groupSelectEl.style.display = groups.length > 1 ? "block" : "none";
        groupSelectEl.addEventListener("change", () => {
          currentGroup = groupSelectEl.value;
          resetLegendState(currentGroup);
          renderLegend();
          renderFrame(currentFrameIndex);
        });

        sliderEl.max = String(Math.max(frameSpecs.length - 1, 0));
        sliderEl.addEventListener("input", () => {
          pause();
          renderFrame(Number(sliderEl.value));
        });
        playButtonEl.addEventListener("click", play);
        pauseButtonEl.addEventListener("click", pause);
      }

      function animate() {
        window.requestAnimationFrame(animate);
        controls.update();
        renderer.render(scene, camera);
      }

      buildAxes();
      initControls();
      initSkyPanel();
      resetLegendState(currentGroup);
      renderLegend();
      renderFrame(currentFrameIndex);
      resize();
      animate();

      canvas.addEventListener("pointermove", onPointerMove);
      canvas.addEventListener("pointerleave", onPointerLeave);
      canvas.addEventListener("click", onCanvasClick);
      window.addEventListener("pointermove", onSkyPointerMove);
      window.addEventListener("pointerup", onSkyPointerEnd);
      window.addEventListener("pointercancel", onSkyPointerEnd);
      window.addEventListener("resize", resize);
      if (typeof ResizeObserver !== "undefined") {
        const observer = new ResizeObserver(() => resize());
        observer.observe(root);
      }
      })();
    </script>
  </body>
</html>
"""


class ThreeJSFigure:
    """Minimal HTML-backed figure wrapper for standalone three.js exports."""

    def __init__(self, scene_spec: dict[str, Any]):
        self.scene_spec = scene_spec
        self._root_id = f"oviz-three-{uuid.uuid4().hex}"

    def to_dict(self) -> dict[str, Any]:
        return self.scene_spec

    def to_html(self) -> str:
        width = int(self.scene_spec.get("width") or 900)
        height = int(self.scene_spec.get("height") or 700)
        html = _THREEJS_HTML_TEMPLATE.replace("__ROOT_ID__", self._root_id)
        html = html.replace("__WIDTH_PX__", str(width))
        html = html.replace("__HEIGHT_PX__", str(height))
        html = html.replace("__SCENE_JSON__", json.dumps(self.scene_spec))
        return html

    def _data_url(self) -> str:
        encoded = base64.b64encode(self.to_html().encode("utf-8")).decode("ascii")
        return f"data:text/html;charset=utf-8;base64,{encoded}"

    def _iframe_html(self) -> str:
        return (
            f'<iframe src="{self._data_url()}" '
            'style="width: 100vw; height: 100vh; border: 0; display: block; max-width: none; margin: 0 0 0 calc(50% - 50vw);" '
            f'loading="eager" referrerpolicy="no-referrer"></iframe>'
        )

    def _repr_html_(self) -> str:
        return self._iframe_html()

    def write_html(self, file: str | Path, **_: Any) -> None:
        Path(file).write_text(self.to_html(), encoding="utf-8")

    def show(self) -> str | None:
        try:
            from IPython.display import HTML, display

            display(HTML(self._iframe_html()))
            return None
        except Exception:
            pass

        with tempfile.NamedTemporaryFile("w", suffix=".html", delete=False, encoding="utf-8") as tmp:
            tmp.write(self.to_html())
            out_path = tmp.name

        webbrowser.open(Path(out_path).resolve().as_uri())
        return out_path
