from __future__ import annotations

import tempfile
import uuid
import webbrowser
from pathlib import Path
from typing import Any

from .threejs_embed import (
    DEFAULT_SCENE_SPEC_COMPRESSION_THRESHOLD_BYTES,
    render_threejs_html,
    threejs_data_url,
    threejs_iframe_html,
)
from .threejs_shell import THREEJS_SHELL_HTML
from .threejs_runtime_legend import THREEJS_LEGEND_RUNTIME_JS
from .threejs_runtime_interactions import THREEJS_INTERACTION_RUNTIME_JS
from .threejs_runtime_actions import THREEJS_ACTION_RUNTIME_JS
from .threejs_runtime_scene import THREEJS_SCENE_RUNTIME_JS
from .threejs_runtime_sky import THREEJS_SKY_RUNTIME_JS
from .threejs_runtime_viewer import THREEJS_VIEWER_RUNTIME_JS
from .threejs_runtime_widget_content import THREEJS_WIDGET_CONTENT_RUNTIME_JS
from .threejs_runtime_widgets import THREEJS_WIDGET_RUNTIME_JS


_THREEJS_TOPBAR_HTML = """
      <div class="oviz-three-topbar">
        <div class="oviz-three-topbar-brand">
          <span class="oviz-three-topbar-dot" aria-hidden="true"></span>
          <span>oviz</span>
        </div>
        <div class="oviz-three-title"></div>
        <div class="oviz-three-widget-menu">
          <div class="oviz-three-controls-shell" data-open="false">
            <button class="oviz-three-controls-toggle" type="button" title="Show or hide the global scene controls" aria-expanded="false" aria-haspopup="true">Controls ▸</button>
            <div class="oviz-three-controls-drawer" aria-hidden="true" inert>
                <div class="oviz-three-controls">
                  <div class="oviz-three-controls-title">Controls</div>
                <div class="oviz-three-selection">
                  <div class="oviz-three-selection-row">
                    <button class="oviz-three-selection-clear" type="button" title="Clear current cluster selection">Clear selection</button>
                  </div>
                  <label class="oviz-three-selection-toggle">
                    <input class="oviz-three-click-select-toggle" type="checkbox" />
                    <span>Enable click select</span>
                  </label>
                  <label class="oviz-three-selection-toggle">
                    <input class="oviz-three-volume-lasso-toggle" type="checkbox" />
                    <span>Lasso volumetric data</span>
                  </label>
                </div>
                  <label class="oviz-three-controls-field">
                    <span>Theme</span>
                  <select class="oviz-three-theme-select">
                    <option value="default">Default</option>
                    <option value="dark">Dark</option>
                    <option value="light">Light</option>
                    <option value="midnight">Midnight</option>
                    <option value="glacier">Glacier</option>
                    <option value="ember">Ember</option>
                    <option value="graphite">Graphite</option>
                    <option value="aurora">Aurora</option>
                    <option value="paper">Paper</option>
                  </select>
                </label>
                <div class="oviz-three-controls-row">
                  <label class="oviz-three-controls-field">
                    <span class="oviz-three-scroll-speed-label">Scroll speed</span>
                    <input class="oviz-three-scroll-speed" type="range" min="0.2" max="4" step="0.05" />
                  </label>
                  <label class="oviz-three-controls-field">
                    <span class="oviz-three-camera-fov-label">Camera FOV</span>
                    <input class="oviz-three-camera-fov" type="range" min="0.05" max="120" step="0.05" />
                  </label>
                </div>
                <div class="oviz-three-controls-row">
                  <label class="oviz-three-controls-field">
                    <span class="oviz-three-global-point-size-label">Point size</span>
                    <input class="oviz-three-global-point-size" type="range" min="0.25" max="4" step="0.05" />
                  </label>
                  <label class="oviz-three-controls-field">
                    <span class="oviz-three-global-point-opacity-label">Point opacity</span>
                    <input class="oviz-three-global-point-opacity" type="range" min="0" max="2" step="0.02" />
                  </label>
                </div>
                <label class="oviz-three-controls-field">
                  <span class="oviz-three-global-point-glow-label">Star glow</span>
                  <input class="oviz-three-global-point-glow" type="range" min="0" max="4" step="0.02" />
                </label>
                <div class="oviz-three-controls-row">
                  <label class="oviz-three-controls-field">
                    <span>Focus group</span>
                    <select class="oviz-three-focus-group-select"></select>
                  </label>
                  <label class="oviz-three-controls-field">
                    <span>Fade time (Myr)</span>
                    <input class="oviz-three-fade-time" type="number" min="0" step="0.5" />
                  </label>
                </div>
                <label class="oviz-three-controls-toggle-row">
                  <input class="oviz-three-fade-in-out-toggle" type="checkbox" />
                  <span>Fade in and out</span>
                </label>
                <label class="oviz-three-controls-toggle-row">
                  <input class="oviz-three-size-by-stars-toggle" type="checkbox" />
                  <span>Size points by n_stars</span>
                </label>
                <label class="oviz-three-controls-toggle-row">
                  <input class="oviz-three-axes-visible-toggle" type="checkbox" checked />
                  <span>Show axes</span>
                </label>
                <label class="oviz-three-controls-toggle-row">
                  <input class="oviz-three-galactic-reference-toggle" type="checkbox" checked />
                  <span>Show GC/radius overlays</span>
                </label>
                <label class="oviz-three-controls-toggle-row">
                  <input class="oviz-three-region-labels-toggle" type="checkbox" checked />
                  <span>Show region labels</span>
                </label>
                <div class="oviz-three-controls-actions">
                  <button class="oviz-three-key-help-button" type="button" title="Show keyboard controls">Keyboard help</button>
                  <button class="oviz-three-view-from-earth" type="button" title="Move the camera to the sky view from the observer position">3D View</button>
                  <button class="oviz-three-auto-orbit" type="button" title="Rotate around the current camera target" aria-pressed="false">Orbit camera</button>
                  <button class="oviz-three-reset-camera" type="button" title="Reset the camera to the initial view">Reset camera</button>
                  <button class="oviz-three-reset-controls" type="button" title="Reset the global control sliders">Reset controls</button>
                </div>
                <div class="oviz-three-controls-hint">Point size, glow, and opacity act as global multipliers on top of each trace's existing settings.</div>
              </div>
            </div>
          </div>
          <div class="oviz-three-sky-controls-shell" data-open="false">
            <button class="oviz-three-sky-controls-toggle" type="button" title="Show or hide sky controls" aria-expanded="false" aria-haspopup="true">Sky ▸</button>
            <div class="oviz-three-sky-controls-drawer" aria-hidden="true" inert>
              <div class="oviz-three-controls oviz-three-sky-controls">
                <div class="oviz-three-controls-title">Sky</div>
                <label class="oviz-three-controls-field oviz-three-sky-source-field" hidden>
                  <span>Sky image</span>
                  <select class="oviz-three-sky-source-select" aria-label="Sky image"></select>
                </label>
                <div class="oviz-three-sky-dome-controls" hidden>
                  <div class="oviz-three-sky-add-grid">
                  <label class="oviz-three-controls-field">
                      <span>Add</span>
                      <select class="oviz-three-sky-layer-preset-select" aria-label="Add sky image"></select>
                  </label>
                  <label class="oviz-three-controls-field">
                      <span>HiPS ID / URL</span>
                      <input class="oviz-three-sky-layer-custom-input" type="text" placeholder="Search or paste HiPS ID" list="__ROOT_ID__-sky-hips-search" autocomplete="off" />
                      <datalist class="oviz-three-sky-hips-search" id="__ROOT_ID__-sky-hips-search"></datalist>
                  </label>
                  </div>
                  <div class="oviz-three-controls-actions">
                    <button class="oviz-three-sky-layer-add" type="button">Add</button>
                  </div>
                  <div class="oviz-three-sky-menu-heading">Layers</div>
                  <div class="oviz-three-sky-layer-list" aria-label="Sky layer stack"></div>
                  <div class="oviz-three-sky-dome-status"></div>
                </div>
              </div>
            </div>
          </div>
          <button class="oviz-three-zen-mode" type="button" title="Hide interface panels and keep only the time slider visible">Zen</button>
          <button class="oviz-three-reset-view" type="button" title="Reset the camera and clear current lasso and click selections">Reset</button>
          <button class="oviz-three-save-state" type="button" title="Export an HTML copy of the figure with the current state">Save State</button>
          <select class="oviz-three-widget-select">
            <option value="">Widgets</option>
          </select>
        </div>
      </div>
""".strip()

_THREEJS_MINIMAL_TOPBAR_HTML = """
      <div class="oviz-three-topbar">
        <div class="oviz-three-title"></div>
      </div>
""".strip()


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
      #__ROOT_ID__:focus {
        outline: none;
      }
      #__ROOT_ID__ .oviz-three-canvas {
        display: block;
        position: absolute;
        inset: 0;
        z-index: 1;
        width: 100%;
        height: 100%;
        outline: none;
        background: transparent;
      }
      #__ROOT_ID__ .oviz-three-save-state:disabled {
        opacity: 0.55;
        cursor: default;
      }
      #__ROOT_ID__[data-zen="true"] .oviz-three-key-help,
      #__ROOT_ID__[data-zen="true"] .oviz-three-widget-panel {
        display: none !important;
      }
      #__ROOT_ID__[data-zen="true"] .oviz-three-topbar {
        grid-template-columns: auto;
        grid-template-areas: "actions";
        justify-content: end;
        padding-inline: 12px;
        background: transparent;
        border-bottom-color: transparent;
        box-shadow: none;
      }
      #__ROOT_ID__[data-zen="true"] .oviz-three-topbar-brand,
      #__ROOT_ID__[data-zen="true"] .oviz-three-title,
      #__ROOT_ID__[data-zen="true"] .oviz-three-tools-shell,
      #__ROOT_ID__[data-zen="true"] .oviz-three-controls-shell,
      #__ROOT_ID__[data-zen="true"] .oviz-three-sky-controls-shell,
      #__ROOT_ID__[data-zen="true"] .oviz-three-widget-select,
      #__ROOT_ID__[data-zen="true"] .oviz-three-reset-view,
      #__ROOT_ID__[data-zen="true"] .oviz-three-save-state {
        display: none !important;
      }
      #__ROOT_ID__[data-zen="true"] .oviz-three-widget-menu {
        justify-self: end;
      }
      #__ROOT_ID__[data-minimal="true"] .oviz-three-topbar,
      #__ROOT_ID__[data-minimal="true"] .oviz-three-key-help,
      #__ROOT_ID__[data-minimal="true"] .oviz-three-note,
      #__ROOT_ID__[data-minimal="true"] .oviz-three-legend-popover,
      #__ROOT_ID__[data-minimal="true"] .oviz-three-widget-panel,
      #__ROOT_ID__[data-minimal="true"] .oviz-three-legend-resize {
        display: none !important;
      }
      #__ROOT_ID__[data-minimal="true"] .oviz-three-lasso-overlay {
        display: none;
      }
      #__ROOT_ID__[data-minimal="true"] .oviz-three-legend-panel {
        top: 12px;
        left: 12px;
        width: auto;
        min-width: 0;
        min-height: 0;
        border: 0 !important;
        border-radius: 0;
        background: transparent !important;
        box-shadow: none !important;
        backdrop-filter: none !important;
        -webkit-backdrop-filter: none !important;
      }
      #__ROOT_ID__[data-minimal="true"] .oviz-three-legend-panel-head {
        display: none !important;
      }
      #__ROOT_ID__[data-minimal="true"] .oviz-three-legend-panel-body {
        gap: 0;
        padding: 0;
        background: transparent !important;
        border: 0 !important;
        box-shadow: none !important;
      }
      #__ROOT_ID__[data-minimal="true"] .oviz-three-legend-group-field {
        padding: 0 1px 8px;
      }
      #__ROOT_ID__[data-minimal="true"] .oviz-three-legend-section {
        gap: 0;
        background: transparent !important;
        border: 0 !important;
        box-shadow: none !important;
      }
      #__ROOT_ID__[data-minimal="true"] .oviz-three-legend-section-toggle {
        display: none !important;
      }
      #__ROOT_ID__[data-minimal="true"] .oviz-three-legend-section-chevron {
        display: none;
      }
      #__ROOT_ID__[data-minimal="true"] .oviz-three-legend-section-title {
        font-size: 9px;
        letter-spacing: 0.08em;
      }
      #__ROOT_ID__[data-minimal="true"] .oviz-three-legend-trace-section {
        display: block;
      }
      #__ROOT_ID__[data-minimal="true"] .oviz-three-legend-volume-section {
        display: none !important;
      }
      #__ROOT_ID__[data-minimal="true"] .oviz-three-legend,
      #__ROOT_ID__[data-minimal="true"] .oviz-three-legend-volume-list {
        gap: 0;
        padding: 0;
        border: 0;
        border-radius: 0;
        background: transparent;
        max-height: none;
        overflow: visible;
      }
      #__ROOT_ID__[data-minimal="true"] .oviz-three-legend-entry,
      #__ROOT_ID__[data-minimal="true"] .oviz-three-legend-volume-list .oviz-three-legend-entry {
        display: block;
        padding: 0;
        border: 0 !important;
        border-radius: 0;
        background: transparent !important;
        box-shadow: none !important;
      }
      #__ROOT_ID__[data-minimal="true"] .oviz-three-legend-item {
        display: block;
        width: auto;
        gap: 0;
        padding: 1px 0;
        font: 650 13px/1.2 -apple-system, BlinkMacSystemFont, "SF Pro Text", "Helvetica Neue", sans-serif;
        cursor: pointer;
        pointer-events: auto;
        border: 0 !important;
        border-radius: 0 !important;
        background: transparent !important;
        box-shadow: none !important;
        appearance: none;
        -webkit-appearance: none;
        outline: none;
      }
      #__ROOT_ID__[data-minimal="true"] .oviz-three-legend-swatch {
        display: none;
      }
      #__ROOT_ID__[data-minimal="true"] .oviz-three-legend-row,
      #__ROOT_ID__[data-minimal="true"] .oviz-three-legend-meta {
        background: transparent !important;
        border: 0 !important;
        box-shadow: none !important;
      }
      #__ROOT_ID__[data-minimal="true"] .oviz-three-legend-meta {
        display: inline-block !important;
        min-width: auto !important;
        max-width: none !important;
      }
      #__ROOT_ID__[data-minimal="true"] .oviz-three-legend-name {
        display: inline-block !important;
        min-width: auto !important;
        max-width: none !important;
        overflow: visible !important;
        text-overflow: clip !important;
        white-space: normal !important;
        opacity: 1 !important;
        text-shadow: 0 0 8px rgba(0, 0, 0, 0.92), 0 1px 2px rgba(0, 0, 0, 0.88);
      }
      #__ROOT_ID__[data-minimal="true"] .oviz-three-legend-item[data-active="false"] .oviz-three-legend-name {
        opacity: 0.58 !important;
      }
      #__ROOT_ID__[data-minimal="true"] .oviz-three-legend-item:hover,
      #__ROOT_ID__[data-minimal="true"] .oviz-three-legend-entry[data-editor-open="true"] {
        background: transparent !important;
      }
      #__ROOT_ID__[data-minimal="true"] .oviz-three-legend-edit,
      #__ROOT_ID__[data-minimal="true"] .oviz-three-legend-panel-toggle {
        display: none !important;
      }
      #__ROOT_ID__[data-minimal="true"][data-galactic-simple="true"] .oviz-three-footer {
        left: 29%;
        bottom: 22px;
        transform: translateX(-50%);
        width: auto;
        max-width: min(calc(100vw - 32px), 580px);
        gap: 14px;
      }
      #__ROOT_ID__[data-minimal="true"][data-galactic-simple="true"] .oviz-three-footer button {
        width: 46px;
        height: 46px;
        font-size: 18px;
      }
      #__ROOT_ID__[data-minimal="true"][data-galactic-simple="true"] .oviz-three-time-label {
        min-width: 156px;
        font: 600 16px/1.2 -apple-system, BlinkMacSystemFont, "SF Pro Text", "Helvetica Neue", sans-serif;
      }
      #__ROOT_ID__[data-minimal="true"][data-galactic-simple="true"] .oviz-three-slider-shell {
        flex: 0 0 clamp(220px, 26vw, 300px);
        width: clamp(220px, 26vw, 300px);
        max-width: clamp(220px, 26vw, 300px);
        height: 44px;
      }
      #__ROOT_ID__[data-minimal="true"][data-galactic-simple="true"] .oviz-three-slider-track-wrap {
        height: 30px;
      }
      #__ROOT_ID__[data-minimal="true"][data-galactic-simple="true"] .oviz-three-slider::-webkit-slider-runnable-track {
        height: 5px;
      }
      #__ROOT_ID__[data-minimal="true"][data-galactic-simple="true"] .oviz-three-slider::-webkit-slider-thumb {
        width: 16px;
        height: 16px;
        margin-top: -6px;
      }
      #__ROOT_ID__[data-minimal="true"][data-galactic-simple="true"] .oviz-three-slider::-moz-range-track {
        height: 5px;
      }
      #__ROOT_ID__[data-minimal="true"][data-galactic-simple="true"] .oviz-three-slider::-moz-range-thumb {
        width: 16px;
        height: 16px;
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
        flex: 1 1 auto;
        padding: 0;
        border: 0;
        background: transparent;
        text-align: left;
        cursor: pointer;
        font: 13px Helvetica, Arial, sans-serif;
        line-height: 1.35;
      }
      #__ROOT_ID__ .oviz-three-legend-entry {
        display: flex;
        flex-direction: column;
        gap: 6px;
        padding-bottom: 6px;
        border-bottom: 1px solid rgba(255, 255, 255, 0.08);
      }
      #__ROOT_ID__ .oviz-three-legend-row {
        display: flex;
        align-items: center;
        gap: 8px;
      }
      #__ROOT_ID__ .oviz-three-legend-disclosure {
        flex: 0 0 auto;
        width: 24px;
        height: 24px;
        padding: 0;
        border-radius: 6px;
        border: 1px solid var(--oviz-panel-border);
        background: transparent;
        color: var(--oviz-text);
        cursor: pointer;
        font: 12px Helvetica, Arial, sans-serif;
        line-height: 1;
      }
      #__ROOT_ID__ .oviz-three-legend-disclosure[data-open="true"] {
        background: rgba(255, 255, 255, 0.08);
      }
      #__ROOT_ID__ .oviz-three-legend-controls {
        display: none;
        flex-direction: column;
        gap: 8px;
        padding: 8px 10px;
        margin-top: 2px;
        border-radius: 8px;
        background: rgba(255, 255, 255, 0.04);
        border: 1px solid rgba(255, 255, 255, 0.06);
      }
      #__ROOT_ID__ .oviz-three-legend-entry[data-open="true"] .oviz-three-legend-controls {
        display: flex;
      }
      #__ROOT_ID__ .oviz-three-legend-control-row {
        display: grid;
        grid-template-columns: 1fr 1fr;
        gap: 8px;
      }
      #__ROOT_ID__ .oviz-three-legend-field {
        display: flex;
        flex-direction: column;
        gap: 3px;
        min-width: 0;
        color: rgba(255, 255, 255, 0.8);
        font-size: 11px;
      }
      #__ROOT_ID__ .oviz-three-legend-field input[type="range"],
      #__ROOT_ID__ .oviz-three-legend-field input[type="number"],
      #__ROOT_ID__ .oviz-three-legend-field select,
      #__ROOT_ID__ .oviz-three-legend-field input[type="color"] {
        width: 100%;
        min-width: 0;
        box-sizing: border-box;
      }
      #__ROOT_ID__ .oviz-three-legend-field input[type="number"],
      #__ROOT_ID__ .oviz-three-legend-field select {
        padding: 2px 0 3px;
        border-radius: 0;
        border: 0;
        border-bottom: 1px solid rgba(255, 255, 255, 0.16);
        background: transparent;
        color: var(--oviz-text);
        font: 12px Helvetica, Arial, sans-serif;
      }
      #__ROOT_ID__ .oviz-three-legend-field input[type="range"] {
        accent-color: var(--oviz-axis);
      }
      #__ROOT_ID__ .oviz-three-legend-field input[type="color"] {
        height: 26px;
        border: 0;
        border-radius: 0;
        padding: 0;
        background: transparent;
      }
      #__ROOT_ID__ .oviz-three-legend-toggle {
        display: flex;
        align-items: center;
        gap: 8px;
        color: rgba(255, 255, 255, 0.82);
        font-size: 11px;
      }
      #__ROOT_ID__ .oviz-three-legend-toggle input {
        margin: 0;
        accent-color: var(--oviz-axis);
      }
      #__ROOT_ID__ .oviz-three-legend-summary {
        color: rgba(255, 255, 255, 0.72);
        font-size: 10px;
        line-height: 1.4;
        white-space: pre-wrap;
        opacity: 1;
      }
      #__ROOT_ID__ .oviz-three-legend-item[data-active="false"] {
        opacity: 0.38;
        text-decoration: line-through;
      }
      #__ROOT_ID__ .oviz-three-legend-item[data-active="true"] {
        opacity: 1.0;
      }
      #__ROOT_ID__ .oviz-three-legend-item:hover {
        opacity: 1.0;
      }
      #__ROOT_ID__ .oviz-three-selection-toggle {
        display: flex;
        align-items: center;
        gap: 8px;
        font-size: 13px;
        color: var(--oviz-text);
      }
      #__ROOT_ID__ .oviz-three-selection-toggle input {
        margin: 0;
        accent-color: var(--oviz-axis);
      }
      #__ROOT_ID__ .oviz-three-selection {
        display: flex;
        flex-direction: column;
        gap: 8px;
        padding: 10px 12px;
        border-radius: 8px;
        border: 1px solid var(--oviz-panel-border);
        background: var(--oviz-panel-bg);
      }
      #__ROOT_ID__ .oviz-three-selection-row {
        display: flex;
        gap: 8px;
      }
      #__ROOT_ID__ .oviz-three-selection button {
        border: 1px solid var(--oviz-panel-border);
        border-radius: 6px;
        background: transparent;
        color: var(--oviz-text);
        cursor: pointer;
        font: 12px Helvetica, Arial, sans-serif;
        padding: 6px 10px;
      }
      #__ROOT_ID__ .oviz-three-selection button[data-active="true"] {
        background: rgba(110, 140, 255, 0.20);
        border-color: rgba(140, 170, 255, 0.7);
      }
      #__ROOT_ID__ .oviz-three-selection button:disabled {
        opacity: 0.45;
        cursor: default;
      }
      #__ROOT_ID__ .oviz-three-selection-readout {
        color: var(--oviz-text);
        font-size: 12px;
        line-height: 1.35;
        white-space: pre-wrap;
      }
      #__ROOT_ID__ .oviz-three-controls {
        display: flex;
        flex-direction: column;
        gap: 10px;
        padding: 10px 12px;
        border-radius: 8px;
        border: 1px solid var(--oviz-panel-border);
        background: var(--oviz-panel-bg);
      }
      #__ROOT_ID__ .oviz-three-controls-title {
        font-size: 14px;
        font-weight: 600;
        color: var(--oviz-text);
      }
      #__ROOT_ID__ .oviz-three-controls-row {
        display: grid;
        grid-template-columns: 1fr 1fr;
        gap: 8px;
      }
      #__ROOT_ID__ .oviz-three-controls-field {
        display: flex;
        flex-direction: column;
        gap: 4px;
        min-width: 0;
        font-size: 12px;
        color: var(--oviz-text);
      }
      #__ROOT_ID__ .oviz-three-controls-field input[type="range"],
      #__ROOT_ID__ .oviz-three-controls-field select,
      #__ROOT_ID__ .oviz-three-controls-field input[type="number"],
      #__ROOT_ID__ .oviz-three-controls-field input[type="text"] {
        width: 100%;
        min-width: 0;
        box-sizing: border-box;
      }
      #__ROOT_ID__ .oviz-three-controls-field select,
      #__ROOT_ID__ .oviz-three-controls-field input[type="number"],
      #__ROOT_ID__ .oviz-three-controls-field input[type="text"] {
        padding: 6px 8px;
        border-radius: 6px;
        border: 1px solid var(--oviz-panel-border);
        background: rgba(0, 0, 0, 0.18);
        color: var(--oviz-text);
        font: 12px Helvetica, Arial, sans-serif;
      }
      #__ROOT_ID__ .oviz-three-controls-field input[type="range"] {
        accent-color: var(--oviz-axis);
      }
      #__ROOT_ID__ .oviz-three-sky-dome-controls {
        display: flex;
        flex-direction: column;
        gap: 9px;
        padding-top: 8px;
        border-top: 1px solid rgba(255, 255, 255, 0.08);
      }
      #__ROOT_ID__ .oviz-three-sky-menu-heading {
        color: rgba(255, 255, 255, 0.72);
        font-size: 11px;
        font-weight: 700;
        letter-spacing: 0;
      }
      #__ROOT_ID__ .oviz-three-sky-add-grid {
        display: grid;
        grid-template-columns: minmax(0, 1fr);
        gap: 7px;
      }
      #__ROOT_ID__ .oviz-three-sky-stretch {
        display: flex;
        flex-direction: column;
        gap: 7px;
        padding-top: 8px;
        border-top: 1px solid rgba(255, 255, 255, 0.08);
      }
      #__ROOT_ID__ .oviz-three-sky-layer-list {
        display: flex;
        flex-direction: column;
        gap: 7px;
        max-height: min(58vh, 520px);
        overflow: auto;
      }
      #__ROOT_ID__ .oviz-three-sky-layer-row {
        min-width: 0;
        border: 1px solid rgba(255, 255, 255, 0.10);
        border-radius: 7px;
        background: rgba(255, 255, 255, 0.035);
      }
      #__ROOT_ID__ .oviz-three-sky-layer-row[data-active="true"] {
        color: rgba(255, 255, 255, 0.96);
      }
      #__ROOT_ID__ .oviz-three-sky-layer-summary {
        display: grid;
        grid-template-columns: minmax(0, 1fr) auto;
        align-items: center;
        gap: 8px;
        padding: 7px 8px;
        cursor: pointer;
        list-style: none;
      }
      #__ROOT_ID__ .oviz-three-sky-layer-summary::-webkit-details-marker {
        display: none;
      }
      #__ROOT_ID__ .oviz-three-sky-layer-name {
        min-width: 0;
        overflow: hidden;
        text-overflow: ellipsis;
        white-space: nowrap;
        font-weight: 700;
      }
      #__ROOT_ID__ .oviz-three-sky-layer-order {
        color: var(--oviz-muted-text);
        font-size: 11px;
      }
      #__ROOT_ID__ .oviz-three-sky-layer-body {
        display: flex;
        flex-direction: column;
        gap: 7px;
        padding: 0 8px 8px;
      }
      #__ROOT_ID__ .oviz-three-sky-layer-body label {
        display: flex;
        flex-direction: column;
        gap: 4px;
        color: rgba(255, 255, 255, 0.72);
        font-size: 11px;
      }
      #__ROOT_ID__ .oviz-three-sky-layer-body input,
      #__ROOT_ID__ .oviz-three-sky-layer-body select {
        width: 100%;
        min-width: 0;
        box-sizing: border-box;
      }
      #__ROOT_ID__ .oviz-three-sky-layer-grid {
        display: grid;
        grid-template-columns: minmax(0, 1fr) minmax(0, 1fr);
        gap: 7px;
      }
      #__ROOT_ID__ .oviz-three-sky-layer-actions {
        display: flex;
        gap: 6px;
        align-items: center;
      }
      #__ROOT_ID__ .oviz-three-sky-layer-row button {
        min-width: 26px;
        height: 24px;
        padding: 0 7px;
        border: 1px solid var(--oviz-panel-border);
        border-radius: 4px;
        background: transparent;
        color: var(--oviz-text);
        cursor: pointer;
        font: 13px Helvetica, Arial, sans-serif;
        line-height: 1;
      }
      #__ROOT_ID__ .oviz-three-sky-layer-row button:disabled {
        opacity: 0.4;
        cursor: default;
      }
      #__ROOT_ID__ .oviz-three-sky-dome-controls[hidden],
      #__ROOT_ID__ .oviz-three-sky-dome-hips-controls[hidden] {
        display: none !important;
      }
      #__ROOT_ID__ .oviz-three-sky-dome-status {
        color: var(--oviz-muted-text);
        font-size: 11px;
        line-height: 1.35;
        min-height: 15px;
      }
      #__ROOT_ID__ .oviz-three-controls-toggle-row {
        display: flex;
        align-items: center;
        gap: 8px;
        color: var(--oviz-text);
        font-size: 12px;
      }
      #__ROOT_ID__ .oviz-three-controls-toggle-row input {
        margin: 0;
        accent-color: var(--oviz-axis);
      }
      #__ROOT_ID__ .oviz-three-controls-actions {
        display: flex;
        gap: 8px;
        flex-wrap: wrap;
      }
      #__ROOT_ID__ .oviz-three-controls-actions button {
        border: 1px solid var(--oviz-panel-border);
        border-radius: 6px;
        background: transparent;
        color: var(--oviz-text);
        cursor: pointer;
        font: 12px Helvetica, Arial, sans-serif;
        padding: 6px 10px;
      }
      #__ROOT_ID__ .oviz-three-controls-actions button[data-active="true"] {
        background: rgba(110, 140, 255, 0.20);
        border-color: rgba(140, 170, 255, 0.7);
      }
      #__ROOT_ID__ .oviz-three-controls-hint {
        color: var(--oviz-text);
        font-size: 11px;
        line-height: 1.4;
        white-space: pre-wrap;
        opacity: 0.88;
      }
      #__ROOT_ID__ .oviz-three-manual-labels {
        display: flex;
        flex-direction: column;
        gap: 8px;
        padding-top: 8px;
        border-top: 1px solid rgba(255, 255, 255, 0.08);
      }
      #__ROOT_ID__ .oviz-three-manual-label-actions {
        align-items: center;
      }
      #__ROOT_ID__ .oviz-three-key-help {
        position: absolute;
        top: 64px;
        right: 12px;
        z-index: 9;
        display: none;
        width: min(420px, 42vw);
        max-height: min(72vh, 760px);
        overflow: auto;
        padding: 14px 16px 16px;
        border-radius: 12px;
        border: 1px solid var(--oviz-panel-border);
        background: var(--oviz-panel-solid);
        color: var(--oviz-text);
        box-shadow: 0 18px 46px rgba(0, 0, 0, 0.30);
      }
      #__ROOT_ID__ .oviz-three-key-help[data-open="true"] {
        display: block;
      }
      #__ROOT_ID__ .oviz-three-key-help-head {
        display: flex;
        align-items: center;
        justify-content: space-between;
        gap: 12px;
        margin-bottom: 10px;
      }
      #__ROOT_ID__ .oviz-three-key-help-title {
        font-size: 16px;
        font-weight: 600;
      }
      #__ROOT_ID__ .oviz-three-key-help-close {
        border: 1px solid var(--oviz-panel-border);
        border-radius: 6px;
        background: transparent;
        color: var(--oviz-text);
        cursor: pointer;
        font: 12px Helvetica, Arial, sans-serif;
        padding: 5px 8px;
      }
      #__ROOT_ID__ .oviz-three-key-help-text {
        margin: 0 0 12px;
        font-size: 12px;
        line-height: 1.45;
        opacity: 0.9;
      }
      #__ROOT_ID__ .oviz-three-key-help-grid {
        display: grid;
        grid-template-columns: minmax(110px, 150px) 1fr;
        gap: 8px 12px;
        font-size: 12px;
        line-height: 1.45;
      }
      #__ROOT_ID__ .oviz-three-key-help-keys {
        font-family: Menlo, Monaco, Consolas, monospace;
        color: var(--oviz-axis);
        white-space: nowrap;
      }
      #__ROOT_ID__ .oviz-three-volume {
        display: none;
        flex-direction: column;
        gap: 8px;
        padding: 10px 2px 2px;
      }
      #__ROOT_ID__ .oviz-three-volume[data-enabled="true"] {
        display: flex;
      }
      #__ROOT_ID__ .oviz-three-volume-title {
        font-size: 14px;
        font-weight: 600;
        color: var(--oviz-text);
      }
      #__ROOT_ID__ .oviz-three-volume-row {
        display: grid;
        grid-template-columns: 1fr 1fr;
        gap: 8px;
      }
      #__ROOT_ID__ .oviz-three-volume-field {
        display: flex;
        flex-direction: column;
        gap: 4px;
        font-size: 12px;
        color: var(--oviz-text);
      }
      #__ROOT_ID__ .oviz-three-volume-toggle {
        display: flex;
        align-items: center;
        gap: 8px;
        font-size: 12px;
        color: var(--oviz-text);
      }
      #__ROOT_ID__ .oviz-three-volume-toggle input {
        margin: 0;
        accent-color: var(--oviz-axis);
      }
      #__ROOT_ID__ .oviz-three-volume select,
      #__ROOT_ID__ .oviz-three-volume input[type="number"],
      #__ROOT_ID__ .oviz-three-volume input[type="range"] {
        width: 100%;
        min-width: 0;
        box-sizing: border-box;
      }
      #__ROOT_ID__ .oviz-three-volume select,
      #__ROOT_ID__ .oviz-three-volume input[type="number"] {
        padding: 6px 8px;
        border-radius: 6px;
        border: 1px solid var(--oviz-panel-border);
        background: rgba(0, 0, 0, 0.18);
        color: var(--oviz-text);
        font: 12px Helvetica, Arial, sans-serif;
      }
      #__ROOT_ID__ .oviz-three-volume input[type="range"] {
        accent-color: var(--oviz-axis);
      }
      #__ROOT_ID__ .oviz-three-volume-summary {
        color: var(--oviz-text);
        font-size: 11px;
        line-height: 1.45;
        white-space: pre-wrap;
        opacity: 0.88;
      }
      #__ROOT_ID__ .oviz-three-slider {
        position: relative;
        z-index: 2;
        width: 100%;
        margin: 0;
        accent-color: var(--oviz-axis);
        background: transparent;
        -webkit-appearance: none;
        appearance: none;
      }
      #__ROOT_ID__ .oviz-three-slider::-webkit-slider-runnable-track {
        height: 4px;
        border-radius: 999px;
        background: rgba(255, 255, 255, 0.18);
      }
      #__ROOT_ID__ .oviz-three-slider::-webkit-slider-thumb {
        -webkit-appearance: none;
        appearance: none;
        width: 14px;
        height: 14px;
        margin-top: -5px;
        border: 1px solid rgba(255, 255, 255, 0.18);
        border-radius: 50%;
        background: var(--oviz-axis);
        box-shadow: 0 0 0 3px rgba(0, 0, 0, 0.18);
      }
      #__ROOT_ID__ .oviz-three-slider::-moz-range-track {
        height: 4px;
        border-radius: 999px;
        background: rgba(255, 255, 255, 0.18);
      }
      #__ROOT_ID__ .oviz-three-slider::-moz-range-thumb {
        width: 14px;
        height: 14px;
        border: 1px solid rgba(255, 255, 255, 0.18);
        border-radius: 50%;
        background: var(--oviz-axis);
        box-shadow: 0 0 0 3px rgba(0, 0, 0, 0.18);
      }
      #__ROOT_ID__ .oviz-three-slider-ticks,
      #__ROOT_ID__ .oviz-three-slider-labels {
        position: absolute;
        left: 7px;
        right: 7px;
        pointer-events: none;
        display: none;
      }
      #__ROOT_ID__ .oviz-three-slider-ticks {
        inset-block: 0;
      }
      #__ROOT_ID__ .oviz-three-slider-labels {
        position: absolute;
        left: 0;
        right: 0;
        top: 34px;
        height: 16px;
      }
      #__ROOT_ID__ .oviz-three-slider-tick {
        position: absolute;
        top: 50%;
        transform: translate(-50%, -50%);
      }
      #__ROOT_ID__ .oviz-three-slider-tick-label {
        position: absolute;
        transform: translateX(-50%);
      }
      #__ROOT_ID__ .oviz-three-slider-tick {
        width: 2px;
        border-radius: 999px;
        background: rgba(255, 255, 255, 0.34);
      }
      #__ROOT_ID__ .oviz-three-slider-tick-minor {
        height: 6px;
        opacity: 0.78;
      }
      #__ROOT_ID__ .oviz-three-slider-tick-major {
        height: 11px;
        background: rgba(255, 255, 255, 0.56);
      }
      #__ROOT_ID__ .oviz-three-slider-tick[data-active="true"] {
        background: var(--oviz-axis);
        opacity: 0.95;
      }
      #__ROOT_ID__ .oviz-three-slider-tick-label {
        top: 0;
        color: var(--oviz-muted-text);
        font: 500 10px/1 -apple-system, BlinkMacSystemFont, "SF Pro Text", "Helvetica Neue", sans-serif;
        white-space: nowrap;
        letter-spacing: 0.01em;
        font-variant-numeric: tabular-nums;
      }
      #__ROOT_ID__ .oviz-three-slider-tick-label[data-active="true"] {
        color: var(--oviz-text);
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
      #__ROOT_ID__ .oviz-three-scale-label {
        color: var(--oviz-text);
        font: 12px Helvetica, Arial, sans-serif;
        white-space: nowrap;
      }
      #__ROOT_ID__ .oviz-three-scale-line {
        position: relative;
        width: 120px;
        height: 12px;
      }
      #__ROOT_ID__ .oviz-three-scale-line::before {
        content: "";
        position: absolute;
        left: 0;
        right: 0;
        top: 5px;
        border-top: 2px solid var(--oviz-axis);
      }
      #__ROOT_ID__ .oviz-three-scale-line::after {
        content: "";
        position: absolute;
        left: 0;
        top: 1px;
        width: 100%;
        height: 10px;
        border-left: 2px solid var(--oviz-axis);
        border-right: 2px solid var(--oviz-axis);
        box-sizing: border-box;
      }
      #__ROOT_ID__ .oviz-three-lasso-overlay {
        position: absolute;
        inset: 0;
        z-index: 7;
        display: none;
        pointer-events: none;
      }
      #__ROOT_ID__ .oviz-three-lasso-overlay[data-active="true"] {
        display: block;
      }
      #__ROOT_ID__ .oviz-three-lasso-overlay svg {
        width: 100%;
        height: 100%;
        display: block;
      }
      #__ROOT_ID__ .oviz-three-lasso-overlay polyline {
        fill: rgba(110, 197, 255, 0.10);
        stroke: var(--oviz-footprint);
        stroke-width: 2;
        stroke-linejoin: round;
        stroke-linecap: round;
        vector-effect: non-scaling-stroke;
      }
      #__ROOT_ID__ .oviz-three-widget-panel {
        position: absolute;
        z-index: 8;
        display: none;
        overflow: hidden;
        border: 1px solid var(--oviz-panel-border);
        border-radius: 10px;
        background: var(--oviz-panel-solid);
        box-shadow: 0 18px 46px rgba(0, 0, 0, 0.30);
        touch-action: none;
      }
      #__ROOT_ID__ .oviz-three-widget-panel[data-mode="normal"] {
        display: block;
      }
      #__ROOT_ID__ .oviz-three-widget-panel[data-mode="fullscreen"] {
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
      #__ROOT_ID__ .oviz-three-sky-panel {
        top: 72px;
        right: 12px;
        width: min(38vw, 560px);
        height: min(56vh, 560px);
      }
      #__ROOT_ID__ .oviz-three-age-panel {
        top: 96px;
        right: 44px;
        width: min(36vw, 540px);
        height: min(52vh, 460px);
      }
      #__ROOT_ID__ .oviz-three-filter-panel {
        top: 108px;
        right: 72px;
        width: min(34vw, 480px);
        height: min(40vh, 340px);
      }
      #__ROOT_ID__ .oviz-three-box-panel {
        right: 18px;
        bottom: 86px;
        width: min(40vw, 580px);
        height: min(46vh, 420px);
      }
      #__ROOT_ID__ .oviz-three-widget-title {
        flex: 1 1 auto;
        min-width: 0;
        overflow: hidden;
        text-overflow: ellipsis;
        white-space: nowrap;
      }
      #__ROOT_ID__ .oviz-three-widget-window-controls {
        display: inline-flex;
        align-items: center;
        gap: 6px;
        flex: 0 0 auto;
      }
      #__ROOT_ID__ .oviz-three-sky-frame {
        position: absolute;
        top: 34px;
        left: 0;
        width: 100%;
        height: calc(100% - 34px);
        border: 0;
        display: block;
        background: var(--oviz-panel-solid);
      }
      #__ROOT_ID__ .oviz-three-sky-dome-frame {
        position: absolute;
        top: 0;
        left: 0;
        width: 100%;
        height: 100%;
        border: 0;
        opacity: 0;
        pointer-events: none;
        z-index: 0;
        transform-origin: 50% 50%;
        will-change: opacity, transform;
        transition: none;
      }
      #__ROOT_ID__ .oviz-three-age-body {
        position: absolute;
        top: 34px;
        right: 0;
        bottom: 0;
        left: 0;
        display: flex;
        flex-direction: column;
        gap: 8px;
        padding: 10px;
        box-sizing: border-box;
      }
      #__ROOT_ID__ .oviz-three-age-kde-shell {
        position: relative;
        flex: 1 1 auto;
        min-height: 220px;
        border: 1px solid var(--oviz-panel-border);
        border-radius: 8px;
        overflow: hidden;
        background: rgba(0, 0, 0, 0.16);
      }
      #__ROOT_ID__ .oviz-three-age-canvas {
        width: 100%;
        height: 100%;
        display: block;
        background: var(--oviz-panel-solid);
      }
      #__ROOT_ID__ .oviz-three-age-filter {
        display: flex;
        flex-direction: column;
        gap: 6px;
        flex: 0 0 auto;
        position: relative;
        z-index: 1;
      }
      #__ROOT_ID__ .oviz-three-age-filter-hint {
        color: var(--oviz-axis);
        font: 10px Menlo, Monaco, Consolas, monospace;
        white-space: nowrap;
        overflow: hidden;
        text-overflow: ellipsis;
      }
      #__ROOT_ID__ .oviz-three-filter-body {
        position: absolute;
        top: 34px;
        right: 0;
        bottom: 0;
        left: 0;
        display: flex;
        flex-direction: column;
        gap: 8px;
        padding: 10px;
        box-sizing: border-box;
      }
      #__ROOT_ID__ .oviz-three-filter-row {
        display: flex;
        gap: 8px;
        align-items: center;
      }
      #__ROOT_ID__ .oviz-three-filter-field {
        display: flex;
        flex-direction: column;
        gap: 4px;
        flex: 1 1 auto;
        color: var(--oviz-text);
        font: 11px Helvetica, Arial, sans-serif;
      }
      #__ROOT_ID__ .oviz-three-filter-field select {
        width: 100%;
      }
      #__ROOT_ID__ .oviz-three-filter-hist {
        position: relative;
        flex: 1 1 auto;
        min-height: 120px;
        border: 1px solid var(--oviz-panel-border);
        border-radius: 8px;
        overflow: hidden;
        background: rgba(0, 0, 0, 0.16);
      }
      #__ROOT_ID__ .oviz-three-filter-canvas {
        width: 100%;
        height: 100%;
        display: block;
        background: transparent;
      }
      #__ROOT_ID__ .oviz-three-filter-slider-shell,
      #__ROOT_ID__ .oviz-three-age-filter-slider-shell {
        position: relative;
        height: 30px;
        margin-top: 2px;
      }
      #__ROOT_ID__ .oviz-three-filter-slider-shell::before,
      #__ROOT_ID__ .oviz-three-age-filter-slider-shell::before {
        content: "";
        position: absolute;
        left: 10px;
        right: 10px;
        top: 14px;
        height: 2px;
        background: rgba(255, 255, 255, 0.18);
      }
      #__ROOT_ID__ .oviz-three-filter-range,
      #__ROOT_ID__ .oviz-three-age-filter-range {
        position: absolute;
        inset: 0;
        width: 100%;
        margin: 0;
        background: transparent;
        pointer-events: none;
        -webkit-appearance: none;
        appearance: none;
      }
      #__ROOT_ID__ .oviz-three-filter-range::-webkit-slider-runnable-track,
      #__ROOT_ID__ .oviz-three-age-filter-range::-webkit-slider-runnable-track {
        height: 30px;
        background: transparent;
      }
      #__ROOT_ID__ .oviz-three-filter-range::-webkit-slider-thumb,
      #__ROOT_ID__ .oviz-three-age-filter-range::-webkit-slider-thumb {
        -webkit-appearance: none;
        appearance: none;
        width: 12px;
        height: 12px;
        margin-top: 9px;
        border: 1px solid var(--oviz-panel-border);
        border-radius: 50%;
        background: var(--oviz-text);
        cursor: pointer;
        pointer-events: auto;
      }
      #__ROOT_ID__ .oviz-three-filter-range::-moz-range-track,
      #__ROOT_ID__ .oviz-three-age-filter-range::-moz-range-track {
        height: 30px;
        background: transparent;
      }
      #__ROOT_ID__ .oviz-three-filter-range::-moz-range-thumb,
      #__ROOT_ID__ .oviz-three-age-filter-range::-moz-range-thumb {
        width: 12px;
        height: 12px;
        border: 1px solid var(--oviz-panel-border);
        border-radius: 50%;
        background: var(--oviz-text);
        cursor: pointer;
        pointer-events: auto;
      }
      #__ROOT_ID__ .oviz-three-filter-labels,
      #__ROOT_ID__ .oviz-three-age-filter-labels {
        display: flex;
        justify-content: space-between;
        color: var(--oviz-axis);
        font: 10px Menlo, Monaco, Consolas, monospace;
      }
      #__ROOT_ID__ .oviz-three-dendrogram-body {
        position: absolute;
        top: 34px;
        right: 0;
        bottom: 0;
        left: 0;
        display: flex;
        flex-direction: column;
        gap: 8px;
        padding: 10px;
        box-sizing: border-box;
      }
      #__ROOT_ID__ .oviz-three-dendrogram-row {
        display: flex;
        gap: 8px;
        align-items: flex-end;
        flex-wrap: wrap;
      }
      #__ROOT_ID__ .oviz-three-dendrogram-field {
        display: flex;
        flex-direction: column;
        gap: 4px;
        flex: 1 1 auto;
        min-width: 120px;
        color: var(--oviz-text);
        font: 11px Helvetica, Arial, sans-serif;
      }
      #__ROOT_ID__ .oviz-three-dendrogram-field--trace {
        flex: 1.6 1 180px;
      }
      #__ROOT_ID__ .oviz-three-dendrogram-field select,
      #__ROOT_ID__ .oviz-three-dendrogram-field input {
        width: 100%;
      }
      #__ROOT_ID__ .oviz-three-dendrogram-shell {
        position: relative;
        flex: 1 1 auto;
        min-height: 180px;
        border: 1px solid var(--oviz-panel-border);
        border-radius: 8px;
        overflow: hidden;
        background: rgba(0, 0, 0, 0.16);
      }
      #__ROOT_ID__ .oviz-three-dendrogram-canvas {
        width: 100%;
        height: 100%;
        display: block;
        background: transparent;
      }
      #__ROOT_ID__ .oviz-three-dendrogram-hint {
        color: var(--oviz-axis);
        font: 10px Menlo, Monaco, Consolas, monospace;
        white-space: nowrap;
        overflow: hidden;
        text-overflow: ellipsis;
      }
      #__ROOT_ID__ .oviz-three-box-body {
        position: absolute;
        top: 34px;
        right: 0;
        bottom: 0;
        left: 0;
        display: flex;
        flex-direction: column;
        gap: 8px;
        padding: 10px;
        box-sizing: border-box;
      }
      #__ROOT_ID__ .oviz-three-box-toolbar {
        display: flex;
        align-items: flex-start;
        justify-content: space-between;
        gap: 8px;
        flex-wrap: wrap;
      }
      #__ROOT_ID__ .oviz-three-box-summary {
        flex: 1 1 auto;
        color: var(--oviz-text);
        font: 11px Menlo, Monaco, Consolas, monospace;
        line-height: 1.35;
        white-space: pre-wrap;
      }
      #__ROOT_ID__ .oviz-three-box-controls {
        display: flex;
        align-items: flex-end;
        gap: 8px;
        flex-wrap: wrap;
        justify-content: flex-end;
      }
      #__ROOT_ID__ .oviz-three-box-field {
        display: flex;
        flex-direction: column;
        gap: 4px;
        color: var(--oviz-text);
        font: 10px Helvetica, Arial, sans-serif;
      }
      #__ROOT_ID__ .oviz-three-box-field input {
        width: 82px;
      }
      #__ROOT_ID__ .oviz-three-box-reset {
        flex: 0 0 auto;
      }
      #__ROOT_ID__ .oviz-three-box-shell {
        position: relative;
        flex: 1 1 auto;
        min-height: 220px;
        border: 1px solid var(--oviz-panel-border);
        border-radius: 8px;
        overflow: hidden;
        background: rgba(0, 0, 0, 0.16);
      }
      #__ROOT_ID__ .oviz-three-box-canvas {
        width: 100%;
        height: 100%;
        display: block;
        background: transparent;
      }
      #__ROOT_ID__ .oviz-three-box-hint {
        color: var(--oviz-axis);
        font: 10px Menlo, Monaco, Consolas, monospace;
        white-space: nowrap;
        overflow: hidden;
        text-overflow: ellipsis;
      }
      #__ROOT_ID__ .oviz-three-widget-panel button {
        border: 1px solid var(--oviz-panel-border);
        border-radius: 4px;
        background: var(--oviz-panel-bg);
        color: var(--oviz-text);
        cursor: pointer;
        font: 10px Menlo, Monaco, Consolas, monospace;
        padding: 2px 6px;
      }
      #__ROOT_ID__ .oviz-three-window-button {
        position: relative;
        z-index: 3;
        width: 20px;
        height: 20px;
        min-width: 20px;
        padding: 0;
        display: inline-flex;
        align-items: center;
        justify-content: center;
        border-radius: 6px;
        background: rgba(255, 255, 255, 0.03);
        color: var(--oviz-text);
        cursor: pointer;
        font-size: 0;
        line-height: 0;
      }
      #__ROOT_ID__ .oviz-three-window-button:hover {
        background: rgba(255, 255, 255, 0.10);
      }
      #__ROOT_ID__ .oviz-three-window-button-min::before {
        content: "";
        display: block;
        width: 10px;
        height: 1.5px;
        background: currentColor;
        transform: translateY(3px);
      }
      #__ROOT_ID__ .oviz-three-window-button-max::before {
        content: "";
        display: block;
        width: 9px;
        height: 9px;
        border: 1.5px solid currentColor;
        box-sizing: border-box;
      }
      #__ROOT_ID__ .oviz-three-sky-actions {
        display: none !important;
      }
      #__ROOT_ID__ .oviz-three-widget-resize {
        position: absolute;
        width: 14px;
        height: 14px;
        z-index: 3;
        background: rgba(120, 120, 120, 0.10);
        display: block;
        touch-action: none;
      }
      #__ROOT_ID__ .oviz-three-widget-panel[data-mode="fullscreen"] .oviz-three-widget-resize {
        display: none;
      }
      #__ROOT_ID__ .oviz-three-widget-resize[data-dir="nw"] {
        top: 0;
        left: 0;
        cursor: nwse-resize;
      }
      #__ROOT_ID__ .oviz-three-widget-resize[data-dir="ne"] {
        top: 0;
        right: 0;
        cursor: nesw-resize;
      }
      #__ROOT_ID__ .oviz-three-widget-resize[data-dir="sw"] {
        bottom: 0;
        left: 0;
        cursor: nesw-resize;
      }
      #__ROOT_ID__ .oviz-three-widget-resize[data-dir="se"] {
        right: 0;
        bottom: 0;
        cursor: nwse-resize;
      }
      #__ROOT_ID__ {
        --oviz-panel-bg: rgba(18, 22, 28, 0.50);
        --oviz-panel-border: rgba(255, 255, 255, 0.11);
        --oviz-panel-solid: rgba(19, 23, 30, 0.76);
        --oviz-axis: #c4c9d1;
        --oviz-text: #eef2f7;
        --oviz-muted-text: rgba(238, 242, 247, 0.58);
        --oviz-accent: #d8dde6;
        --oviz-accent-strong: #eef1f5;
        --oviz-shadow-lg: 0 24px 70px rgba(0, 0, 0, 0.34);
        --oviz-shadow-md: 0 18px 42px rgba(0, 0, 0, 0.24);
      }
      #__ROOT_ID__::before {
        content: "";
        position: absolute;
        inset: 0;
        background:
          radial-gradient(circle at 16% 12%, rgba(255, 255, 255, 0.07), transparent 24%),
          radial-gradient(circle at 84% 10%, rgba(255, 255, 255, 0.04), transparent 18%),
          linear-gradient(180deg, rgba(255, 255, 255, 0.025), transparent 20%);
        pointer-events: none;
        z-index: 0;
      }
      #__ROOT_ID__ .oviz-three-topbar {
        position: absolute;
        top: 14px;
        left: 14px;
        right: 14px;
        z-index: 7;
        display: grid;
        grid-template-columns: 1fr auto 1fr;
        align-items: start;
        gap: 0;
        padding: 0;
        border: 0;
        border-radius: 0;
        background: transparent;
        box-shadow: none;
        backdrop-filter: none;
        -webkit-backdrop-filter: none;
      }
      #__ROOT_ID__ .oviz-three-topbar-brand {
        display: none;
      }
      #__ROOT_ID__ .oviz-three-topbar-dot {
        width: 10px;
        height: 10px;
        border-radius: 50%;
        background: radial-gradient(circle at 35% 35%, #ffffff, #d3d8df 62%, #9ea7b3);
        box-shadow: 0 0 0 1px rgba(255, 255, 255, 0.18), 0 0 14px rgba(255, 255, 255, 0.12);
      }
      #__ROOT_ID__::before {
        background:
          radial-gradient(circle at 16% 12%, rgba(255, 255, 255, 0.014), transparent 18%),
          linear-gradient(180deg, rgba(255, 255, 255, 0.006), transparent 14%);
      }
      #__ROOT_ID__ .oviz-three-title {
        position: static;
        transform: none;
        min-width: 0;
        text-align: center;
        justify-self: center;
        align-self: start;
        padding: 5px 10px;
        border: 1px solid rgba(255, 255, 255, 0.05);
        border-radius: 4px;
        background: rgba(20, 22, 26, 0.90);
        box-shadow: 0 3px 10px rgba(0, 0, 0, 0.10);
        backdrop-filter: blur(4px) saturate(105%);
        -webkit-backdrop-filter: blur(4px) saturate(105%);
        font: 600 13px/1.15 -apple-system, BlinkMacSystemFont, "SF Pro Text", "Helvetica Neue", sans-serif;
        letter-spacing: -0.015em;
        color: var(--oviz-text);
        pointer-events: none;
      }
      #__ROOT_ID__ .oviz-three-action-bar {
        position: absolute;
        top: 16px;
        right: 18px;
        z-index: 7;
        display: none;
        flex-direction: column;
        align-items: flex-end;
        justify-content: flex-start;
        gap: 8px;
        width: auto;
        max-width: min(calc(100vw - 32px), 320px);
        padding: 0;
        border: 0;
        border-radius: 0;
        background: transparent;
        box-shadow: none;
        backdrop-filter: none;
        -webkit-backdrop-filter: none;
      }
      #__ROOT_ID__ .oviz-three-action-button {
        height: auto;
        padding: 0;
        border: 0;
        border-radius: 0;
        background: transparent;
        color: rgba(255, 255, 255, 0.42);
        cursor: pointer;
        box-shadow: none;
        font: 500 16px/1.2 -apple-system, BlinkMacSystemFont, "SF Pro Text", "Helvetica Neue", sans-serif;
        letter-spacing: 0.01em;
        transition: color 160ms ease, opacity 160ms ease;
      }
      #__ROOT_ID__ .oviz-three-action-button:hover {
        color: rgba(255, 255, 255, 0.96);
      }
      #__ROOT_ID__ .oviz-three-action-button[data-active="true"] {
        color: rgba(255, 255, 255, 0.96);
      }
      #__ROOT_ID__ .oviz-three-widget-menu {
        position: static;
        grid-column: 3;
        display: flex;
        align-items: center;
        justify-self: end;
        justify-content: flex-end;
        flex-wrap: wrap;
        gap: 4px;
        width: fit-content;
        max-width: min(calc(100vw - 28px), 760px);
        padding: 4px;
        border: 1px solid rgba(255, 255, 255, 0.05);
        border-radius: 4px;
        background: rgba(20, 22, 26, 0.90);
        box-shadow: 0 4px 14px rgba(0, 0, 0, 0.10);
        backdrop-filter: blur(4px) saturate(104%);
        -webkit-backdrop-filter: blur(4px) saturate(104%);
      }
      #__ROOT_ID__ .oviz-three-widget-menu button,
      #__ROOT_ID__ .oviz-three-widget-menu select,
      #__ROOT_ID__ .oviz-three-tools-toggle,
      #__ROOT_ID__ .oviz-three-controls-toggle {
        height: 30px;
        border-radius: 4px;
        border: 1px solid rgba(255, 255, 255, 0.06);
        background: rgba(29, 32, 37, 0.98);
        color: var(--oviz-text);
        cursor: pointer;
        box-shadow: none;
        font: 500 11px/1 -apple-system, BlinkMacSystemFont, "SF Pro Text", "Helvetica Neue", sans-serif;
      }
      #__ROOT_ID__ .oviz-three-widget-menu button:hover,
      #__ROOT_ID__ .oviz-three-widget-menu select:hover,
      #__ROOT_ID__ .oviz-three-tools-toggle:hover,
      #__ROOT_ID__ .oviz-three-controls-toggle:hover {
        background: rgba(38, 42, 48, 0.98);
        border-color: rgba(255, 255, 255, 0.10);
      }
      #__ROOT_ID__ .oviz-three-widget-menu .oviz-three-tools-shell,
      #__ROOT_ID__ .oviz-three-widget-menu .oviz-three-controls-shell {
        position: relative;
        min-height: auto;
        border: 0;
        border-radius: 0;
        background: transparent;
        box-shadow: none;
        backdrop-filter: none;
        -webkit-backdrop-filter: none;
        overflow: visible;
      }
      #__ROOT_ID__ .oviz-three-widget-menu .oviz-three-tools-toggle,
      #__ROOT_ID__ .oviz-three-widget-menu .oviz-three-controls-toggle {
        width: auto;
        padding: 0 12px;
        justify-content: center;
        text-align: center;
      }
      #__ROOT_ID__[data-topbar-density="compact"] .oviz-three-widget-menu {
        gap: 3px;
        padding: 3px;
        max-width: min(calc(100vw - 24px), 520px);
      }
      #__ROOT_ID__[data-topbar-density="compact"] .oviz-three-action-bar {
        gap: 7px;
        top: 14px;
        right: 16px;
        max-width: min(calc(100vw - 24px), 280px);
      }
      #__ROOT_ID__[data-topbar-density="compact"] .oviz-three-action-button {
        font-size: 14px;
      }
      #__ROOT_ID__[data-topbar-density="compact"] .oviz-three-widget-menu button,
      #__ROOT_ID__[data-topbar-density="compact"] .oviz-three-widget-menu select,
      #__ROOT_ID__[data-topbar-density="compact"] .oviz-three-tools-toggle,
      #__ROOT_ID__[data-topbar-density="compact"] .oviz-three-controls-toggle {
        height: 28px;
        font-size: 10px;
      }
      #__ROOT_ID__[data-topbar-density="compact"] .oviz-three-widget-menu .oviz-three-tools-toggle,
      #__ROOT_ID__[data-topbar-density="compact"] .oviz-three-widget-menu .oviz-three-controls-toggle {
        padding: 0 10px;
      }
      #__ROOT_ID__[data-topbar-density="compact"] .oviz-three-widget-select {
        min-width: 124px;
        max-width: 148px;
      }
      #__ROOT_ID__[data-topbar-density="stacked"] .oviz-three-topbar {
        grid-template-columns: 1fr;
        grid-template-areas:
          "title"
          "actions";
        align-items: start;
        gap: 8px;
      }
      #__ROOT_ID__[data-topbar-density="stacked"] .oviz-three-title {
        grid-area: title;
        justify-self: center;
        max-width: calc(100vw - 24px);
      }
      #__ROOT_ID__[data-topbar-density="stacked"] .oviz-three-widget-menu {
        grid-area: actions;
        justify-self: end;
        justify-content: flex-end;
        width: min(calc(100vw - 24px), 420px);
        max-width: calc(100vw - 24px);
        gap: 3px;
        padding: 3px;
      }
      #__ROOT_ID__[data-topbar-density="stacked"] .oviz-three-widget-menu button,
      #__ROOT_ID__[data-topbar-density="stacked"] .oviz-three-widget-menu select,
      #__ROOT_ID__[data-topbar-density="stacked"] .oviz-three-tools-toggle,
      #__ROOT_ID__[data-topbar-density="stacked"] .oviz-three-controls-toggle {
        height: 28px;
        font-size: 10px;
      }
      #__ROOT_ID__[data-topbar-density="stacked"] .oviz-three-widget-menu .oviz-three-tools-toggle,
      #__ROOT_ID__[data-topbar-density="stacked"] .oviz-three-widget-menu .oviz-three-controls-toggle {
        padding: 0 9px;
      }
      #__ROOT_ID__[data-topbar-density="stacked"] .oviz-three-widget-select {
        min-width: 112px;
        max-width: 132px;
      }
      #__ROOT_ID__[data-topbar-density="stacked"] .oviz-three-action-bar {
        top: 12px;
        right: 14px;
        max-width: min(calc(100vw - 24px), 240px);
        gap: 6px;
      }
      #__ROOT_ID__[data-topbar-density="stacked"] .oviz-three-action-button {
        font-size: 13px;
      }
      #__ROOT_ID__ .oviz-three-widget-menu .oviz-three-tools-drawer,
      #__ROOT_ID__ .oviz-three-widget-menu .oviz-three-controls-drawer {
        position: absolute;
        top: calc(100% + 10px);
        right: 0;
        left: auto;
        width: min(340px, calc(100vw - 32px));
        max-height: min(74vh, 760px);
        overflow: auto;
        padding: 12px;
        border: 1px solid rgba(255, 255, 255, 0.08);
        border-radius: 20px;
        background: linear-gradient(180deg, rgba(23, 26, 31, 0.54), rgba(12, 15, 19, 0.34));
        box-shadow: var(--oviz-shadow-md);
        backdrop-filter: blur(22px) saturate(124%);
        -webkit-backdrop-filter: blur(22px) saturate(124%);
        opacity: 0;
        pointer-events: none;
        transform: translateY(-4px);
        transition: opacity 0.16s ease, transform 0.16s ease;
      }
      #__ROOT_ID__ .oviz-three-widget-menu .oviz-three-tools-shell[data-open="true"] .oviz-three-tools-drawer,
      #__ROOT_ID__ .oviz-three-widget-menu .oviz-three-controls-shell[data-open="true"] .oviz-three-controls-drawer {
        opacity: 1;
        pointer-events: auto;
        transform: translateY(0);
      }
      #__ROOT_ID__ .oviz-three-widget-menu .oviz-three-selection,
      #__ROOT_ID__ .oviz-three-widget-menu .oviz-three-controls {
        padding: 0;
        border: 0;
        border-radius: 0;
        background: transparent;
      }
      #__ROOT_ID__ .oviz-three-legend-panel,
      #__ROOT_ID__ .oviz-three-legend-popover {
        border: 1px solid rgba(255, 255, 255, 0.08);
        border-radius: 20px;
        background: linear-gradient(180deg, rgba(23, 26, 31, 0.34), rgba(12, 15, 19, 0.16));
        box-shadow: 0 18px 42px rgba(0, 0, 0, 0.16);
        backdrop-filter: blur(22px) saturate(120%);
        -webkit-backdrop-filter: blur(22px) saturate(120%);
      }
      #__ROOT_ID__ .oviz-three-legend-panel,
      #__ROOT_ID__ .oviz-three-legend-popover,
      #__ROOT_ID__ .oviz-three-widget-menu .oviz-three-tools-drawer,
      #__ROOT_ID__ .oviz-three-widget-menu .oviz-three-controls-drawer,
      #__ROOT_ID__ .oviz-three-widget-panel {
        border: 1px solid rgba(255, 255, 255, 0.06);
        border-radius: 10px;
        background: rgba(18, 20, 24, 0.76);
        box-shadow: 0 8px 22px rgba(0, 0, 0, 0.16);
        backdrop-filter: blur(8px) saturate(108%);
        -webkit-backdrop-filter: blur(8px) saturate(108%);
      }
      #__ROOT_ID__ .oviz-three-legend-panel {
        position: absolute;
        top: 14px;
        left: 14px;
        z-index: 8;
        width: auto;
        min-width: 0;
        min-height: 0;
        display: flex;
        flex-direction: column;
        overflow: visible;
        border: 0 !important;
        border-radius: 0;
        background: transparent !important;
        box-shadow: none !important;
        backdrop-filter: none !important;
        -webkit-backdrop-filter: none !important;
      }
      #__ROOT_ID__ .oviz-three-legend-panel[data-open="false"] {
        min-height: 0;
      }
      #__ROOT_ID__ .oviz-three-legend-panel-head {
        display: none !important;
      }
      #__ROOT_ID__ .oviz-three-legend-panel[data-dragging="true"] .oviz-three-legend-panel-head {
        cursor: grabbing;
      }
      #__ROOT_ID__ .oviz-three-legend-panel-title {
        color: var(--oviz-text);
        font: 600 13px/1.2 -apple-system, BlinkMacSystemFont, "SF Pro Text", "Helvetica Neue", sans-serif;
        letter-spacing: -0.01em;
      }
      #__ROOT_ID__ .oviz-three-legend-panel-title {
        font-size: 12px;
      }
      #__ROOT_ID__ .oviz-three-legend-panel-toggle {
        width: 24px;
        height: 24px;
        border: 0;
        border-radius: 999px;
        background: rgba(255, 255, 255, 0.05);
        color: var(--oviz-muted-text);
        font: 600 12px/1 -apple-system, BlinkMacSystemFont, "SF Pro Text", "Helvetica Neue", sans-serif;
        cursor: pointer;
        position: relative;
        z-index: 1;
        touch-action: manipulation;
      }
      #__ROOT_ID__ .oviz-three-legend-panel-body {
        display: flex;
        flex-direction: column;
        flex: 1 1 auto;
        min-height: 0;
        padding: 0;
        overflow: auto;
        transition: max-height 0.18s ease, padding 0.18s ease, opacity 0.18s ease;
        background: transparent !important;
        border: 0 !important;
        box-shadow: none !important;
      }
      #__ROOT_ID__ .oviz-three-legend-panel button,
      #__ROOT_ID__ .oviz-three-legend-panel button:hover,
      #__ROOT_ID__ .oviz-three-legend-panel button:focus,
      #__ROOT_ID__ .oviz-three-legend-panel button:active {
        background: transparent !important;
        border: 0 !important;
        box-shadow: none !important;
        outline: none;
        -webkit-appearance: none;
        appearance: none;
      }
      #__ROOT_ID__ .oviz-three-legend-panel-body {
        padding: 0;
        scroll-padding-bottom: 8px;
      }
      #__ROOT_ID__ .oviz-three-legend-panel[data-open="false"] .oviz-three-legend-panel-body {
        max-height: 0;
        padding-top: 0;
        padding-bottom: 0;
        opacity: 0;
        overflow: hidden;
      }
      #__ROOT_ID__ .oviz-three-legend {
        display: flex;
        flex-direction: column;
        gap: 2px;
        padding: 0;
        border: 0;
        border-radius: 0;
        background: transparent;
        box-shadow: none;
        max-height: 58vh;
        overflow: auto;
      }
      #__ROOT_ID__ .oviz-three-legend-section {
        display: flex;
        flex-direction: column;
        gap: 0;
      }
      #__ROOT_ID__ .oviz-three-legend-section[data-open="false"] {
        gap: 0;
      }
      #__ROOT_ID__ .oviz-three-legend-section[data-empty="true"] {
        display: none;
      }
      #__ROOT_ID__ .oviz-three-legend-section-toggle {
        display: none !important;
      }
      #__ROOT_ID__ .oviz-three-legend-section-title {
        color: var(--oviz-muted-text);
        font: 600 10px/1.2 -apple-system, BlinkMacSystemFont, "SF Pro Text", "Helvetica Neue", sans-serif;
        letter-spacing: 0.04em;
        text-transform: uppercase;
      }
      #__ROOT_ID__ .oviz-three-legend-section-chevron {
        flex: 0 0 auto;
        color: var(--oviz-muted-text);
        font: 600 11px/1 -apple-system, BlinkMacSystemFont, "SF Pro Text", "Helvetica Neue", sans-serif;
        transition: transform 0.16s ease;
      }
      #__ROOT_ID__ .oviz-three-legend-section[data-open="false"] .oviz-three-legend {
        display: none;
      }
      #__ROOT_ID__ .oviz-three-legend-section[data-open="false"] .oviz-three-legend-section-chevron {
        transform: rotate(-90deg);
      }
      #__ROOT_ID__ .oviz-three-legend-volume-list .oviz-three-legend-entry {
        background: transparent;
      }
      #__ROOT_ID__ .oviz-three-legend-volume-section {
        padding-top: 0;
        border-top: 0;
      }
      #__ROOT_ID__ .oviz-three-legend-section:last-child {
        padding-bottom: 0;
      }
      #__ROOT_ID__ .oviz-three-legend-group-field {
        display: flex;
        flex-direction: column;
        gap: 4px;
        padding: 0 1px 8px;
        color: var(--oviz-muted-text);
        font: 500 10px/1.2 -apple-system, BlinkMacSystemFont, "SF Pro Text", "Helvetica Neue", sans-serif;
      }
      #__ROOT_ID__ .oviz-three-legend-group-field > span {
        display: none;
      }
      #__ROOT_ID__ .oviz-three-group-select {
        display: none;
      }
      @keyframes oviz-three-group-option-drop {
        from {
          opacity: 0;
          transform: translateY(-10px);
        }
        to {
          opacity: 1;
          transform: translateY(0);
        }
      }
      #__ROOT_ID__ .oviz-three-group-dropdown {
        width: min(236px, calc(100vw - 28px));
      }
      #__ROOT_ID__ .oviz-three-group-trigger {
        display: inline-flex;
        align-items: center;
        justify-content: flex-start;
        gap: 5px;
        width: auto;
        min-width: 0;
        max-width: min(236px, calc(100vw - 28px));
        padding: 0 0 5px;
        border: 0;
        border-bottom: 0;
        background: transparent;
        color: rgba(255, 255, 255, 0.96);
        cursor: pointer;
        font: 760 15px/1.16 -apple-system, BlinkMacSystemFont, "SF Pro Text", "Helvetica Neue", sans-serif;
        letter-spacing: 0;
        text-align: left;
      }
      #__ROOT_ID__ .oviz-three-group-current {
        display: block;
        flex: 0 1 auto;
        min-width: 0;
        max-width: min(212px, calc(100vw - 52px));
        overflow: hidden;
        text-overflow: ellipsis;
        white-space: nowrap;
      }
      #__ROOT_ID__ .oviz-three-group-chevron {
        flex: 0 0 auto;
        order: -1;
        color: rgba(255, 255, 255, 0.56);
        font: 760 14px/1 -apple-system, BlinkMacSystemFont, "SF Pro Text", "Helvetica Neue", sans-serif;
        transform: translateY(-1px) rotate(0deg);
        transition: transform 220ms ease, color 160ms ease;
      }
      #__ROOT_ID__ .oviz-three-group-dropdown[data-open="true"] .oviz-three-group-chevron {
        color: rgba(255, 255, 255, 0.9);
        transform: translateY(-1px) rotate(180deg);
      }
      #__ROOT_ID__ .oviz-three-group-menu {
        max-height: 0;
        overflow: hidden;
        opacity: 0;
        transform: translateY(-8px);
        transition: max-height 340ms cubic-bezier(0.22, 1, 0.36, 1), opacity 220ms ease, transform 340ms cubic-bezier(0.22, 1, 0.36, 1), margin-top 340ms ease;
        margin-top: 0;
      }
      #__ROOT_ID__ .oviz-three-group-dropdown[data-open="true"] .oviz-three-group-menu {
        max-height: 260px;
        opacity: 1;
        transform: translateY(0);
        margin-top: 6px;
      }
      #__ROOT_ID__ .oviz-three-group-menu-list {
        display: flex;
        flex-direction: column;
        gap: 1px;
        padding: 0 0 4px;
      }
      #__ROOT_ID__ .oviz-three-group-option {
        width: min(236px, calc(100vw - 28px));
        padding: 3px 0 4px;
        border: 0;
        border-bottom: 1px solid transparent;
        background: transparent;
        color: rgba(235, 238, 245, 0.5);
        cursor: pointer;
        font: 650 11px/1.2 -apple-system, BlinkMacSystemFont, "SF Pro Text", "Helvetica Neue", sans-serif;
        letter-spacing: 0;
        opacity: 0;
        text-align: left;
        transform: translateY(-10px);
        transition: color 160ms ease, border-color 160ms ease;
        will-change: opacity, transform;
      }
      #__ROOT_ID__ .oviz-three-group-dropdown[data-open="true"] .oviz-three-group-option {
        opacity: 1;
        transform: translateY(0);
        animation: oviz-three-group-option-drop 340ms cubic-bezier(0.16, 1, 0.3, 1) both;
        animation-delay: var(--group-option-delay, 0ms);
      }
      #__ROOT_ID__ .oviz-three-group-option:hover,
      #__ROOT_ID__ .oviz-three-group-option:focus-visible {
        color: rgba(255, 255, 255, 0.9);
        outline: none;
      }
      #__ROOT_ID__ .oviz-three-group-option[data-active="true"] {
        color: rgba(255, 255, 255, 0.96);
        border-bottom-color: transparent;
      }
      #__ROOT_ID__ .oviz-three-legend-title {
        padding: 0 2px 6px;
        color: var(--oviz-muted-text);
        font: 500 10px/1.3 -apple-system, BlinkMacSystemFont, "SF Pro Text", "Helvetica Neue", sans-serif;
      }
      #__ROOT_ID__ .oviz-three-legend-entry {
        display: block;
        padding: 0;
        border: 0;
        border-radius: 0;
        background: transparent;
        transition: opacity 0.14s ease;
      }
      #__ROOT_ID__ .oviz-three-legend-row {
        display: inline-flex;
        align-items: center;
        gap: 0;
        min-width: 0;
        background: transparent !important;
        border: 0 !important;
        box-shadow: none !important;
      }
      #__ROOT_ID__ .oviz-three-legend-entry[data-active="false"] {
        opacity: 1;
      }
      #__ROOT_ID__ .oviz-three-legend-item {
        display: inline-flex;
        align-items: center;
        gap: 0;
        min-width: 0;
        width: auto;
        border: 0;
        background: transparent;
        padding: 0;
        color: inherit;
        text-align: left;
        font: 650 13px/1.2 -apple-system, BlinkMacSystemFont, "SF Pro Text", "Helvetica Neue", sans-serif;
      }
      #__ROOT_ID__ .oviz-three-legend-item:hover,
      #__ROOT_ID__ .oviz-three-legend-entry[data-editor-open="true"] {
        background: transparent;
      }
      #__ROOT_ID__ .oviz-three-legend-swatch {
        display: none;
      }
      #__ROOT_ID__ .oviz-three-legend-meta {
        display: block;
        min-width: 0;
      }
      #__ROOT_ID__ .oviz-three-legend-name {
        min-width: 0;
        overflow: visible;
        text-overflow: clip;
        white-space: normal;
        color: inherit;
        opacity: 1;
        text-shadow: 0 0 8px rgba(0, 0, 0, 0.92), 0 1px 2px rgba(0, 0, 0, 0.88);
      }
      #__ROOT_ID__ .oviz-three-legend-edit {
        flex: 0 0 auto;
        align-self: center;
        width: auto;
        height: auto;
        margin-left: 0;
        padding: 0 0 0 1px;
        border: 0;
        border-radius: 0;
        background: transparent;
        color: rgba(255, 255, 255, 0.46);
        cursor: pointer;
        font: 700 13px/1 -apple-system, BlinkMacSystemFont, "SF Pro Text", "Helvetica Neue", sans-serif;
      }
      #__ROOT_ID__ .oviz-three-legend-edit:hover,
      #__ROOT_ID__ .oviz-three-legend-edit[data-open="true"] {
        color: rgba(255, 255, 255, 0.94);
      }
      #__ROOT_ID__ .oviz-three-legend-item[data-active="false"] .oviz-three-legend-name {
        opacity: 0.52;
        text-decoration: none;
      }
      #__ROOT_ID__ .oviz-three-legend-controls {
        display: flex;
        flex-direction: column;
        gap: 8px;
        margin: 0 0 0 0;
        padding: 0 0 0 10px;
        border-radius: 0;
        background: transparent;
        border: 0;
        max-height: 0;
        opacity: 0;
        overflow: hidden;
        pointer-events: none;
        transform: translateY(-6px);
        transition: max-height 300ms cubic-bezier(0.22, 1, 0.36, 1), opacity 190ms ease, transform 300ms cubic-bezier(0.22, 1, 0.36, 1), margin 300ms ease;
      }
      #__ROOT_ID__ .oviz-three-legend-controls[data-visible="true"] {
        max-height: 460px;
        opacity: 1;
        pointer-events: auto;
        transform: translateY(0);
        margin: 1px 0 8px 0;
      }
      #__ROOT_ID__ .oviz-three-legend-resize,
      #__ROOT_ID__ .oviz-three-legend-popover {
        display: none !important;
      }
      #__ROOT_ID__ .oviz-three-legend-resize {
        position: absolute;
        width: 16px;
        height: 16px;
        z-index: 1;
        touch-action: none;
      }
      #__ROOT_ID__ .oviz-three-legend-resize[data-dir="nw"] {
        top: -1px;
        left: -1px;
        cursor: nwse-resize;
      }
      #__ROOT_ID__ .oviz-three-legend-resize[data-dir="ne"] {
        top: -1px;
        right: -1px;
        cursor: nesw-resize;
      }
      #__ROOT_ID__ .oviz-three-legend-resize[data-dir="sw"] {
        left: -1px;
        bottom: -1px;
        cursor: nesw-resize;
      }
      #__ROOT_ID__ .oviz-three-legend-resize[data-dir="se"] {
        right: -1px;
        bottom: -1px;
        cursor: nwse-resize;
      }
      #__ROOT_ID__ .oviz-three-legend-panel[data-open="false"] .oviz-three-legend-resize {
        pointer-events: none;
        opacity: 0;
      }
      #__ROOT_ID__ .oviz-three-widget-menu .oviz-three-tools-drawer,
      #__ROOT_ID__ .oviz-three-widget-menu .oviz-three-controls-drawer {
        top: calc(100% + 8px);
        width: min(320px, calc(100vw - 32px));
        padding: 10px;
      }
      #__ROOT_ID__ .oviz-three-widget-drag {
        position: absolute;
        top: 0;
        left: 0;
        right: 0;
        z-index: 2;
        display: flex;
        align-items: center;
        justify-content: space-between;
        gap: 12px;
        height: 34px;
        padding: 0 10px;
        text-transform: none;
        letter-spacing: 0;
        font: 500 12px/1 -apple-system, BlinkMacSystemFont, "SF Pro Text", "Helvetica Neue", sans-serif;
        color: var(--oviz-text);
        background: rgba(24, 26, 30, 0.96);
        border-bottom: 1px solid var(--oviz-panel-border);
        cursor: grab;
        user-select: none;
      }
      #__ROOT_ID__ .oviz-three-widget-panel[data-mode="fullscreen"] .oviz-three-widget-drag {
        cursor: default;
      }
      #__ROOT_ID__ .oviz-three-footer {
        position: absolute;
        left: 50%;
        bottom: 14px;
        transform: translateX(-50%);
        z-index: 6;
        display: flex;
        align-items: center;
        justify-content: center;
        gap: 10px;
        width: auto;
        max-width: calc(100vw - 28px);
        padding: 5px 10px;
        border-radius: 8px;
        background: rgba(18, 20, 24, 0.94);
        box-shadow: 0 5px 14px rgba(0, 0, 0, 0.12);
        backdrop-filter: blur(4px) saturate(103%);
        -webkit-backdrop-filter: blur(4px) saturate(103%);
      }
      #__ROOT_ID__ .oviz-three-footer button {
        width: 32px;
        height: 32px;
        line-height: 1;
        padding: 0;
        border-radius: 6px;
        border: 1px solid rgba(255, 255, 255, 0.04);
        background: rgba(255, 255, 255, 0.02);
        -webkit-appearance: none;
        appearance: none;
      }
      #__ROOT_ID__ .oviz-three-footer button[data-active="true"] {
        color: var(--oviz-axis);
        border-color: rgba(255, 255, 255, 0.14);
        background: rgba(255, 255, 255, 0.08);
      }
      #__ROOT_ID__ .oviz-three-earth-view-toggle {
        position: absolute;
        right: 18px;
        bottom: 18px;
        z-index: 7;
        min-width: 104px;
        height: 34px;
        padding: 0 14px;
        border-radius: 6px;
        border: 1px solid rgba(255, 255, 255, 0.10);
        background: rgba(18, 20, 24, 0.92);
        color: rgba(255, 255, 255, 0.88);
        box-shadow: 0 6px 18px rgba(0, 0, 0, 0.18);
        cursor: pointer;
        font: 650 12px/1 -apple-system, BlinkMacSystemFont, "SF Pro Text", "Helvetica Neue", sans-serif;
        letter-spacing: 0;
        backdrop-filter: blur(4px) saturate(103%);
        -webkit-backdrop-filter: blur(4px) saturate(103%);
        transition: transform 160ms ease, color 160ms ease, background 180ms ease, border-color 180ms ease;
      }
      #__ROOT_ID__ .oviz-three-earth-view-toggle:hover {
        transform: translateY(-1px);
        color: rgba(255, 255, 255, 0.98);
        border-color: rgba(255, 255, 255, 0.18);
        background: rgba(28, 31, 36, 0.96);
      }
      #__ROOT_ID__ .oviz-three-earth-view-toggle[data-active="true"] {
        color: rgba(255, 255, 255, 0.98);
        border-color: rgba(246, 200, 95, 0.34);
        background: rgba(246, 200, 95, 0.12);
      }
      #__ROOT_ID__ .oviz-three-time-label {
        min-width: 124px;
        height: 30px;
        font-size: 13px;
        display: inline-flex;
        align-items: center;
        justify-content: center;
        white-space: nowrap;
        font-variant-numeric: tabular-nums;
        text-align: center;
      }
      #__ROOT_ID__ .oviz-three-slider-shell {
        position: relative;
        flex: 0 0 clamp(196px, 25vw, 288px);
        width: clamp(196px, 25vw, 288px);
        max-width: clamp(196px, 25vw, 288px);
        display: flex;
        align-items: center;
        height: 38px;
        min-width: 0;
      }
      #__ROOT_ID__ .oviz-three-slider-track-wrap {
        position: relative;
        display: flex;
        align-items: center;
        min-width: 0;
        width: 100%;
        height: 24px;
      }
      #__ROOT_ID__ .oviz-three-slider-labels {
        top: 24px;
        height: 12px;
      }
      #__ROOT_ID__ .oviz-three-scale-bar {
        position: absolute;
        left: 18px;
        bottom: 22px;
        z-index: 6;
        display: flex;
        flex-direction: column;
        align-items: flex-start;
        gap: 6px;
        pointer-events: none;
        cursor: default;
        user-select: none;
        touch-action: none;
        border-radius: 4px;
        background: rgba(18, 20, 24, 0.92);
        backdrop-filter: blur(4px) saturate(103%);
        -webkit-backdrop-filter: blur(4px) saturate(103%);
      }
      #__ROOT_ID__ .oviz-three-scale-bar[data-dragging="true"] {
        cursor: default;
      }
      #__ROOT_ID__[data-minimal="true"] .oviz-three-scale-bar {
        cursor: default;
      }
      #__ROOT_ID__[data-minimal="true"] .oviz-three-scale-bar[data-dragging="true"] {
        cursor: default;
      }
      #__ROOT_ID__ .oviz-three-note {
        position: absolute;
        right: 12px;
        bottom: 64px;
        z-index: 6;
        max-width: min(340px, 34vw);
        padding: 8px 10px;
        display: none;
        color: var(--oviz-text);
        font-size: 12px;
        border-radius: 4px;
        background: rgba(18, 20, 24, 0.92);
        box-shadow: 0 5px 14px rgba(0, 0, 0, 0.12);
        backdrop-filter: blur(4px) saturate(103%);
        -webkit-backdrop-filter: blur(4px) saturate(103%);
      }
      #__ROOT_ID__ .oviz-three-sky-frame {
        position: absolute;
        left: 0;
        width: 100%;
        border: 0;
        display: block;
        background: var(--oviz-panel-solid);
        top: 32px;
        height: calc(100% - 32px);
      }
      #__ROOT_ID__ .oviz-three-widget-select,
      #__ROOT_ID__ .oviz-three-group-select,
      #__ROOT_ID__ .oviz-three-controls-field select,
      #__ROOT_ID__ .oviz-three-legend-field select,
      #__ROOT_ID__ .oviz-three-legend-field input[type="number"],
      #__ROOT_ID__ .oviz-three-controls-field input[type="number"],
      #__ROOT_ID__ .oviz-three-controls-field input[type="text"],
      #__ROOT_ID__ .oviz-three-volume select,
      #__ROOT_ID__ .oviz-three-volume input[type="number"],
      #__ROOT_ID__ .oviz-three-filter-field select {
        background: rgba(8, 10, 14, 0.34);
        border: 1px solid rgba(255, 255, 255, 0.10);
        border-radius: 10px;
        color: var(--oviz-text);
        font: 500 12px/1.2 -apple-system, BlinkMacSystemFont, "SF Pro Text", "Helvetica Neue", sans-serif;
      }
      #__ROOT_ID__ .oviz-three-tools-shell,
      #__ROOT_ID__ .oviz-three-controls-shell {
        border: 1px solid rgba(255, 255, 255, 0.08);
        border-radius: 18px;
        background: linear-gradient(180deg, rgba(23, 28, 34, 0.36), rgba(14, 17, 22, 0.18));
        box-shadow: 0 18px 36px rgba(0, 0, 0, 0.18);
        backdrop-filter: blur(18px) saturate(126%);
        -webkit-backdrop-filter: blur(18px) saturate(126%);
        overflow: hidden;
      }
      #__ROOT_ID__ .oviz-three-scene-card {
        display: flex;
        flex-direction: column;
        gap: 10px;
        padding: 12px;
        min-height: 0;
      }
      #__ROOT_ID__ .oviz-three-card-head {
        display: flex;
        flex-direction: column;
        gap: 4px;
      }
      #__ROOT_ID__ .oviz-three-card-eyebrow {
        color: var(--oviz-muted-text);
        font: 600 11px/1 -apple-system, BlinkMacSystemFont, "SF Pro Text", "Helvetica Neue", sans-serif;
        letter-spacing: 0.08em;
        text-transform: uppercase;
      }
      #__ROOT_ID__ .oviz-three-card-caption {
        color: var(--oviz-text);
        font: 600 15px/1.2 -apple-system, BlinkMacSystemFont, "SF Pro Display", "Helvetica Neue", sans-serif;
        letter-spacing: -0.02em;
      }
      #__ROOT_ID__ .oviz-three-scene-meta {
        display: flex;
        align-items: center;
        justify-content: space-between;
        gap: 10px;
      }
      #__ROOT_ID__ .oviz-three-sidebar-field {
        display: flex;
        flex-direction: column;
        gap: 6px;
        color: var(--oviz-muted-text);
        font: 500 11px/1.2 -apple-system, BlinkMacSystemFont, "SF Pro Text", "Helvetica Neue", sans-serif;
      }
      #__ROOT_ID__ .oviz-three-legend-popover {
        position: absolute;
        z-index: 8;
        width: min(244px, calc(100vw - 92px));
        max-height: min(66vh, 560px);
        padding: 10px 10px 11px;
        border: 1px solid rgba(255, 255, 255, 0.10);
        border-radius: 18px;
        background: linear-gradient(180deg, rgba(23, 28, 34, 0.74), rgba(14, 17, 22, 0.50));
        box-shadow: 0 24px 46px rgba(0, 0, 0, 0.26);
        backdrop-filter: blur(22px) saturate(132%);
        -webkit-backdrop-filter: blur(22px) saturate(132%);
        opacity: 0;
        pointer-events: none;
        transform: translateY(6px);
        transition: opacity 0.16s ease, transform 0.16s ease;
        overflow: auto;
      }
      #__ROOT_ID__ .oviz-three-legend-popover[data-open="true"] {
        opacity: 1;
        pointer-events: auto;
        transform: translateY(0);
      }
      #__ROOT_ID__ .oviz-three-legend-popover-head {
        display: flex;
        align-items: flex-start;
        justify-content: space-between;
        gap: 10px;
        margin-bottom: 8px;
      }
      #__ROOT_ID__ .oviz-three-legend-popover-title {
        display: flex;
        flex-direction: column;
        gap: 2px;
      }
      #__ROOT_ID__ .oviz-three-legend-popover-eyebrow {
        color: var(--oviz-muted-text);
        font: 600 10px/1 -apple-system, BlinkMacSystemFont, "SF Pro Text", "Helvetica Neue", sans-serif;
        letter-spacing: 0.08em;
        text-transform: uppercase;
      }
      #__ROOT_ID__ .oviz-three-legend-popover-name {
        color: var(--oviz-text);
        font: 600 14px/1.25 -apple-system, BlinkMacSystemFont, "SF Pro Display", "Helvetica Neue", sans-serif;
        letter-spacing: -0.02em;
      }
      #__ROOT_ID__ .oviz-three-legend-popover-close {
        width: 26px;
        height: 26px;
        border: 0;
        border-radius: 6px;
        background: rgba(255, 255, 255, 0.05);
        color: var(--oviz-muted-text);
        font: 600 12px/1 -apple-system, BlinkMacSystemFont, "SF Pro Text", "Helvetica Neue", sans-serif;
      }
      #__ROOT_ID__ .oviz-three-legend-popover .oviz-three-legend-controls {
        display: flex;
        margin-top: 2px;
      }
      #__ROOT_ID__ .oviz-three-legend {
        gap: 8px;
        padding: 0;
        border: 0;
        border-radius: 0;
        background: transparent;
        max-height: none;
        min-height: 0;
        overflow: visible;
      }
      #__ROOT_ID__ .oviz-three-legend-controls {
        gap: 8px;
        padding: 7px 8px;
        border-radius: 6px;
        background: rgba(255, 255, 255, 0.022);
      }
      #__ROOT_ID__ .oviz-three-legend-title {
        position: sticky;
        top: 0;
        z-index: 1;
        padding: 5px 7px;
        border: 1px solid rgba(255, 255, 255, 0.06);
        border-radius: 6px;
        background: rgba(255, 255, 255, 0.018);
        color: var(--oviz-muted-text);
        font: 500 10px/1.2 -apple-system, BlinkMacSystemFont, "SF Pro Text", "Helvetica Neue", sans-serif;
        letter-spacing: 0.02em;
      }
      #__ROOT_ID__ .oviz-three-legend-entry {
        gap: 5px;
        padding: 5px 7px;
        border: 1px solid rgba(255, 255, 255, 0.05);
        border-radius: 6px;
        background: rgba(255, 255, 255, 0.022);
      }
      #__ROOT_ID__ .oviz-three-legend-entry:last-child {
        padding-bottom: 5px;
        border-bottom-width: 1px;
        border-bottom-style: solid;
        border-bottom-color: rgba(255, 255, 255, 0.05);
      }
      #__ROOT_ID__ .oviz-three-legend-item {
        display: flex;
        align-items: center;
        gap: 7px;
        font: 600 11px/1.18 -apple-system, BlinkMacSystemFont, "SF Pro Text", "Helvetica Neue", sans-serif;
      }
      #__ROOT_ID__ .oviz-three-legend-swatch {
        width: 8px;
        height: 8px;
        border-radius: 999px;
      }
      #__ROOT_ID__ .oviz-three-legend-meta {
        min-width: 0;
        display: flex;
        flex-direction: column;
        gap: 1px;
      }
      #__ROOT_ID__ .oviz-three-legend-name {
        overflow: hidden;
        text-overflow: ellipsis;
        white-space: nowrap;
      }
      #__ROOT_ID__ .oviz-three-legend-kind {
        color: var(--oviz-muted-text);
        font: 500 10px/1.2 -apple-system, BlinkMacSystemFont, "SF Pro Text", "Helvetica Neue", sans-serif;
      }
      #__ROOT_ID__ .oviz-three-legend-controls {
        gap: 8px;
        padding: 7px 8px;
        border-radius: 6px;
        background: rgba(0, 0, 0, 0.16);
      }
      #__ROOT_ID__ .oviz-three-legend-panel .oviz-three-legend-title,
      #__ROOT_ID__ .oviz-three-legend-panel .oviz-three-legend-entry,
      #__ROOT_ID__ .oviz-three-legend-panel .oviz-three-legend-entry:last-child,
      #__ROOT_ID__ .oviz-three-legend-panel .oviz-three-legend-row,
      #__ROOT_ID__ .oviz-three-legend-panel .oviz-three-legend-item,
      #__ROOT_ID__ .oviz-three-legend-panel .oviz-three-legend-item:hover,
      #__ROOT_ID__ .oviz-three-legend-panel .oviz-three-legend-entry[data-editor-open="true"],
      #__ROOT_ID__ .oviz-three-legend-panel .oviz-three-legend-edit,
      #__ROOT_ID__ .oviz-three-legend-panel .oviz-three-legend-edit:hover,
      #__ROOT_ID__ .oviz-three-legend-panel .oviz-three-legend-controls {
        background: transparent !important;
        border: 0 !important;
        box-shadow: none !important;
      }
      #__ROOT_ID__ .oviz-three-legend-panel .oviz-three-legend-entry,
      #__ROOT_ID__ .oviz-three-legend-panel .oviz-three-legend-entry:last-child {
        display: block !important;
        padding: 0 !important;
      }
      #__ROOT_ID__ .oviz-three-legend-panel .oviz-three-legend-row {
        display: inline-flex !important;
        align-items: center !important;
        gap: 0 !important;
      }
      #__ROOT_ID__ .oviz-three-legend-panel .oviz-three-legend-item {
        display: inline-flex !important;
        width: auto !important;
        gap: 0 !important;
        padding: 0 !important;
      }
      #__ROOT_ID__ .oviz-three-legend-panel .oviz-three-legend-edit {
        margin-left: 0 !important;
        padding-left: 1px !important;
      }
      #__ROOT_ID__ .oviz-three-legend-panel .oviz-three-legend-controls {
        margin: 1px 0 8px 0 !important;
        padding: 0 0 0 10px !important;
      }
      #__ROOT_ID__ .oviz-three-tools-shell,
      #__ROOT_ID__ .oviz-three-controls-shell {
        min-height: 0;
      }
      #__ROOT_ID__ .oviz-three-tools-toggle,
      #__ROOT_ID__ .oviz-three-controls-toggle {
        width: 100%;
        padding: 0 12px;
        justify-content: space-between;
        text-align: left;
      }
      #__ROOT_ID__ .oviz-three-tools-drawer,
      #__ROOT_ID__ .oviz-three-controls-drawer {
        position: static;
        left: auto;
        top: auto;
        width: auto;
        transform: none;
        opacity: 1;
        pointer-events: auto;
        max-height: 0;
        overflow: hidden;
        padding: 0 12px;
        transition: max-height 0.22s ease, padding 0.22s ease;
      }
      #__ROOT_ID__ .oviz-three-tools-shell[data-open="true"] .oviz-three-tools-drawer,
      #__ROOT_ID__ .oviz-three-controls-shell[data-open="true"] .oviz-three-controls-drawer {
        max-height: 1200px;
        padding: 0 12px 12px;
      }
      #__ROOT_ID__ .oviz-three-selection,
      #__ROOT_ID__ .oviz-three-controls {
        border: 0;
        border-radius: 14px;
        background: transparent;
        padding: 10px 2px 2px;
      }
      #__ROOT_ID__ .oviz-three-selection-readout,
      #__ROOT_ID__ .oviz-three-controls-hint,
      #__ROOT_ID__ .oviz-three-legend-summary {
        color: var(--oviz-muted-text);
      }
      #__ROOT_ID__ .oviz-three-key-help {
        top: 78px;
        right: 18px;
        width: min(460px, calc(100vw - 36px));
        border-radius: 22px;
        background: linear-gradient(180deg, rgba(25, 30, 37, 0.82), rgba(14, 17, 22, 0.72));
        box-shadow: var(--oviz-shadow-lg);
        backdrop-filter: blur(24px) saturate(140%);
        -webkit-backdrop-filter: blur(24px) saturate(140%);
      }
      #__ROOT_ID__ .oviz-three-footer {
        left: 50%;
        bottom: 18px;
        width: auto;
        max-width: calc(100vw - 36px);
        padding: 0;
        border-radius: 24px;
        background: transparent;
        box-shadow: none;
        backdrop-filter: none;
        -webkit-backdrop-filter: none;
      }
      #__ROOT_ID__ .oviz-three-footer button {
        width: 38px;
        height: 38px;
        border: 0;
        border-radius: 0;
        background: none;
        color: var(--oviz-text);
        cursor: pointer;
        font-size: 15px;
        line-height: 1;
        margin-top: 0;
        padding: 0;
        box-shadow: none;
        -webkit-appearance: none;
        appearance: none;
      }
      #__ROOT_ID__ .oviz-three-time-label {
        min-width: 136px;
        font: 600 14px/1.2 -apple-system, BlinkMacSystemFont, "SF Pro Text", "Helvetica Neue", sans-serif;
      }
      #__ROOT_ID__ .oviz-three-scale-bar {
        left: 18px;
        bottom: 22px;
        padding: 8px 10px;
        border: 1px solid rgba(255, 255, 255, 0.08);
        border-radius: 10px;
        background: rgba(14, 17, 22, 0.26);
        backdrop-filter: blur(16px) saturate(135%);
        -webkit-backdrop-filter: blur(16px) saturate(135%);
      }
      #__ROOT_ID__ .oviz-three-note {
        right: 16px;
        bottom: 86px;
        border-radius: 16px;
        background: rgba(14, 17, 22, 0.40);
        box-shadow: var(--oviz-shadow-md);
        backdrop-filter: blur(18px) saturate(135%);
        -webkit-backdrop-filter: blur(18px) saturate(135%);
      }
      #__ROOT_ID__ .oviz-three-widget-panel {
        border-radius: 20px;
        border: 1px solid rgba(255, 255, 255, 0.12);
        background: linear-gradient(180deg, rgba(23, 28, 34, 0.82), rgba(14, 17, 22, 0.72));
        box-shadow: var(--oviz-shadow-lg);
        backdrop-filter: blur(24px) saturate(140%);
        -webkit-backdrop-filter: blur(24px) saturate(140%);
      }
      #__ROOT_ID__ .oviz-three-widget-drag {
        height: 34px;
        padding: 0 14px;
        font: 600 11px/1 -apple-system, BlinkMacSystemFont, "SF Pro Text", "Helvetica Neue", sans-serif;
        letter-spacing: 0.08em;
        background: linear-gradient(180deg, rgba(255, 255, 255, 0.06), rgba(255, 255, 255, 0.01));
        border-bottom: 1px solid rgba(255, 255, 255, 0.06);
      }
      #__ROOT_ID__ .oviz-three-widget-panel button {
        border-radius: 10px;
        background: rgba(255, 255, 255, 0.04);
        font: 500 11px/1 -apple-system, BlinkMacSystemFont, "SF Pro Text", "Helvetica Neue", sans-serif;
      }
      #__ROOT_ID__ .oviz-three-legend-panel {
        width: auto !important;
        min-width: 0 !important;
        min-height: 0 !important;
        border: 0 !important;
        border-radius: 0 !important;
        background: transparent !important;
        box-shadow: none !important;
        backdrop-filter: none !important;
        -webkit-backdrop-filter: none !important;
        overflow: visible !important;
      }
      #__ROOT_ID__ .oviz-three-legend-panel-head {
        display: none !important;
      }
      #__ROOT_ID__ .oviz-three-legend-panel-body {
        gap: 0 !important;
        padding: 0 !important;
        background: transparent !important;
        border: 0 !important;
        box-shadow: none !important;
      }
      #__ROOT_ID__ .oviz-three-legend-panel button,
      #__ROOT_ID__ .oviz-three-legend-panel button:hover,
      #__ROOT_ID__ .oviz-three-legend-panel button:focus,
      #__ROOT_ID__ .oviz-three-legend-panel button:active {
        background: transparent !important;
        border: 0 !important;
        box-shadow: none !important;
        outline: none !important;
        -webkit-appearance: none !important;
        appearance: none !important;
      }
      #__ROOT_ID__ .oviz-three-legend-group-field {
        gap: 4px !important;
        padding: 0 1px 8px !important;
      }
      #__ROOT_ID__ .oviz-three-legend-group-field > span {
        display: none !important;
      }
      #__ROOT_ID__ .oviz-three-group-select {
        display: none !important;
      }
      #__ROOT_ID__ .oviz-three-group-dropdown {
        width: min(236px, calc(100vw - 28px)) !important;
      }
      #__ROOT_ID__ .oviz-three-group-trigger {
        display: inline-flex !important;
        align-items: center !important;
        justify-content: flex-start !important;
        gap: 5px !important;
        width: auto !important;
        min-width: 0 !important;
        max-width: min(236px, calc(100vw - 28px)) !important;
        padding: 0 0 5px !important;
        border: 0 !important;
        border-bottom: 0 !important;
        background: transparent !important;
        color: rgba(255, 255, 255, 0.96) !important;
        cursor: pointer !important;
        font: 760 15px/1.16 -apple-system, BlinkMacSystemFont, "SF Pro Text", "Helvetica Neue", sans-serif !important;
        letter-spacing: 0 !important;
        text-align: left !important;
      }
      #__ROOT_ID__ .oviz-three-group-current {
        display: block !important;
        flex: 0 1 auto !important;
        min-width: 0 !important;
        max-width: min(212px, calc(100vw - 52px)) !important;
        overflow: hidden !important;
        text-overflow: ellipsis !important;
        white-space: nowrap !important;
      }
      #__ROOT_ID__ .oviz-three-group-chevron {
        flex: 0 0 auto !important;
        order: -1 !important;
        color: rgba(255, 255, 255, 0.56) !important;
        font: 760 14px/1 -apple-system, BlinkMacSystemFont, "SF Pro Text", "Helvetica Neue", sans-serif !important;
        transform: translateY(-1px) rotate(0deg) !important;
        transition: transform 220ms ease, color 160ms ease !important;
      }
      #__ROOT_ID__ .oviz-three-group-dropdown[data-open="true"] .oviz-three-group-chevron {
        color: rgba(255, 255, 255, 0.9) !important;
        transform: translateY(-1px) rotate(180deg) !important;
      }
      #__ROOT_ID__ .oviz-three-group-menu {
        max-height: 0 !important;
        overflow: hidden !important;
        opacity: 0 !important;
        transform: translateY(-8px) !important;
        transition: max-height 340ms cubic-bezier(0.22, 1, 0.36, 1), opacity 220ms ease, transform 340ms cubic-bezier(0.22, 1, 0.36, 1), margin-top 340ms ease !important;
        margin-top: 0 !important;
      }
      #__ROOT_ID__ .oviz-three-group-dropdown[data-open="true"] .oviz-three-group-menu {
        max-height: 260px !important;
        opacity: 1 !important;
        transform: translateY(0) !important;
        margin-top: 6px !important;
      }
      #__ROOT_ID__ .oviz-three-group-menu-list {
        display: flex !important;
        flex-direction: column !important;
        gap: 1px !important;
        padding: 0 0 4px !important;
      }
      #__ROOT_ID__ .oviz-three-group-option {
        width: min(236px, calc(100vw - 28px)) !important;
        padding: 3px 0 4px !important;
        border: 0 !important;
        border-bottom: 1px solid transparent !important;
        background: transparent !important;
        color: rgba(235, 238, 245, 0.5) !important;
        cursor: pointer !important;
        font: 650 11px/1.2 -apple-system, BlinkMacSystemFont, "SF Pro Text", "Helvetica Neue", sans-serif !important;
        letter-spacing: 0 !important;
        opacity: 0 !important;
        text-align: left !important;
        transform: translateY(-10px) !important;
        transition: color 160ms ease, border-color 160ms ease !important;
        will-change: opacity, transform !important;
      }
      #__ROOT_ID__ .oviz-three-group-dropdown[data-open="true"] .oviz-three-group-option {
        opacity: 1 !important;
        transform: translateY(0) !important;
        animation: oviz-three-group-option-drop 340ms cubic-bezier(0.16, 1, 0.3, 1) both !important;
        animation-delay: var(--group-option-delay, 0ms) !important;
      }
      #__ROOT_ID__ .oviz-three-group-option:hover,
      #__ROOT_ID__ .oviz-three-group-option:focus-visible {
        color: rgba(255, 255, 255, 0.9) !important;
        outline: none !important;
      }
      #__ROOT_ID__ .oviz-three-group-option[data-active="true"] {
        color: rgba(255, 255, 255, 0.96) !important;
        border-bottom-color: transparent !important;
      }
      #__ROOT_ID__ .oviz-three-legend-section {
        gap: 0 !important;
        background: transparent !important;
        border: 0 !important;
        box-shadow: none !important;
      }
      #__ROOT_ID__ .oviz-three-legend-section-toggle {
        display: none !important;
      }
      #__ROOT_ID__ .oviz-three-legend,
      #__ROOT_ID__ .oviz-three-legend-volume-list {
        gap: 0 !important;
        padding: 0 !important;
        border: 0 !important;
        border-radius: 0 !important;
        background: transparent !important;
        max-height: none !important;
        overflow: visible !important;
      }
      #__ROOT_ID__ .oviz-three-legend-entry,
      #__ROOT_ID__ .oviz-three-legend-volume-list .oviz-three-legend-entry,
      #__ROOT_ID__ .oviz-three-legend-entry:last-child {
        display: block !important;
        padding: 0 !important;
        border: 0 !important;
        border-radius: 0 !important;
        background: transparent !important;
        box-shadow: none !important;
      }
      #__ROOT_ID__ .oviz-three-legend-row,
      #__ROOT_ID__ .oviz-three-legend-meta {
        background: transparent !important;
        border: 0 !important;
        box-shadow: none !important;
      }
      #__ROOT_ID__ .oviz-three-legend-row {
        display: inline-grid !important;
        grid-template-columns: 20px auto !important;
        align-items: center !important;
        column-gap: 4px !important;
        min-width: 0 !important;
      }
      #__ROOT_ID__ .oviz-three-legend-item {
        display: inline-flex !important;
        width: auto !important;
        gap: 0 !important;
        padding: 1px 0 !important;
        font: 650 13px/1.2 -apple-system, BlinkMacSystemFont, "SF Pro Text", "Helvetica Neue", sans-serif !important;
        cursor: pointer !important;
        pointer-events: auto !important;
        text-align: left !important;
      }
      #__ROOT_ID__ .oviz-three-legend-swatch {
        display: none !important;
      }
      #__ROOT_ID__ .oviz-three-legend-meta {
        display: inline-block !important;
        min-width: auto !important;
        max-width: none !important;
      }
      #__ROOT_ID__ .oviz-three-legend-name {
        display: inline-block !important;
        min-width: auto !important;
        max-width: none !important;
        overflow: visible !important;
        text-overflow: clip !important;
        white-space: normal !important;
        opacity: 1 !important;
        text-shadow: 0 0 8px rgba(0, 0, 0, 0.92), 0 1px 2px rgba(0, 0, 0, 0.88) !important;
      }
      #__ROOT_ID__ .oviz-three-legend-item[data-active="false"] .oviz-three-legend-name {
        opacity: 0.58 !important;
        text-decoration: none !important;
      }
      #__ROOT_ID__ .oviz-three-legend-item:hover,
      #__ROOT_ID__ .oviz-three-legend-entry[data-editor-open="true"] {
        background: transparent !important;
      }
      #__ROOT_ID__ .oviz-three-legend-edit {
        display: inline-flex !important;
        align-self: center !important;
        justify-content: center !important;
        width: 18px !important;
        height: 18px !important;
        margin-right: 0 !important;
        margin-left: 0 !important;
        padding: 0 !important;
        border: 1px solid currentColor !important;
        border-radius: 999px !important;
        color: rgba(255, 255, 255, 0.46) !important;
        font: 400 16px/1 -apple-system, BlinkMacSystemFont, "SF Pro Text", "Helvetica Neue", sans-serif !important;
      }
      #__ROOT_ID__ .oviz-three-legend-edit:hover,
      #__ROOT_ID__ .oviz-three-legend-edit[data-open="true"] {
        color: currentColor !important;
      }
      #__ROOT_ID__ .oviz-three-legend-controls {
        display: flex !important;
        flex-direction: column !important;
        gap: 8px !important;
        margin: 0 0 0 0 !important;
        padding: 0 0 0 24px !important;
        border: 0 !important;
        border-radius: 0 !important;
        background: transparent !important;
        box-shadow: none !important;
        max-height: 0 !important;
        opacity: 0 !important;
        overflow: hidden !important;
        pointer-events: none !important;
        transform: translateY(-6px) !important;
        transition: max-height 300ms cubic-bezier(0.22, 1, 0.36, 1), opacity 190ms ease, transform 300ms cubic-bezier(0.22, 1, 0.36, 1), margin 300ms ease !important;
      }
      #__ROOT_ID__ .oviz-three-legend-controls[data-visible="true"] {
        max-height: 460px !important;
        opacity: 1 !important;
        pointer-events: auto !important;
        transform: translateY(0) !important;
        margin: 1px 0 8px 0 !important;
      }
      #__ROOT_ID__ .oviz-three-legend-resize,
      #__ROOT_ID__ .oviz-three-legend-popover {
        display: none !important;
      }
      #__ROOT_ID__[data-zen="true"] .oviz-three-legend-panel,
      #__ROOT_ID__[data-zen="true"] .oviz-three-key-help,
      #__ROOT_ID__[data-zen="true"] .oviz-three-widget-panel,
      #__ROOT_ID__[data-zen="true"] .oviz-three-note,
      #__ROOT_ID__[data-zen="true"] .oviz-three-scale-bar {
        display: none !important;
      }
      #__ROOT_ID__ {
        --oviz-instrument-bg: rgba(6, 8, 12, 0.54);
        --oviz-instrument-bg-strong: rgba(8, 10, 15, 0.74);
        --oviz-instrument-border: rgba(238, 242, 247, 0.10);
        --oviz-instrument-border-strong: rgba(246, 200, 95, 0.48);
        --oviz-instrument-hover: rgba(255, 255, 255, 0.055);
        --oviz-instrument-text: rgba(238, 242, 247, 0.88);
        --oviz-instrument-muted: rgba(238, 242, 247, 0.56);
        --oviz-instrument-dim: rgba(238, 242, 247, 0.38);
        --oviz-instrument-accent: #f6c85f;
        --oviz-instrument-ease: cubic-bezier(0.22, 1, 0.36, 1);
      }
      #__ROOT_ID__ .oviz-three-topbar,
      #__ROOT_ID__ .oviz-three-widget-menu,
      #__ROOT_ID__ .oviz-three-widget-menu .oviz-three-tools-shell,
      #__ROOT_ID__ .oviz-three-widget-menu .oviz-three-controls-shell {
        overflow: visible !important;
      }
      #__ROOT_ID__ .oviz-three-widget-menu {
        gap: 12px !important;
        padding: 0 !important;
        border: 0 !important;
        border-radius: 0 !important;
        background: transparent !important;
        box-shadow: none !important;
        backdrop-filter: none !important;
        -webkit-backdrop-filter: none !important;
        max-width: min(calc(100vw - 28px), 920px) !important;
      }
      #__ROOT_ID__ .oviz-three-widget-menu button,
      #__ROOT_ID__ .oviz-three-widget-menu select,
      #__ROOT_ID__ .oviz-three-tools-toggle,
      #__ROOT_ID__ .oviz-three-controls-toggle {
        min-height: 30px !important;
        padding: 0 2px !important;
        border: 0 !important;
        border-radius: 0 !important;
        background: transparent !important;
        color: rgba(238, 242, 247, 0.70) !important;
        box-shadow: none !important;
        font: 720 11.5px/1 -apple-system, BlinkMacSystemFont, "SF Pro Text", "Helvetica Neue", sans-serif !important;
        letter-spacing: 0 !important;
        text-shadow: 0 0 9px rgba(0, 0, 0, 0.84) !important;
        transition: color 160ms ease, box-shadow 180ms ease !important;
      }
      #__ROOT_ID__ .oviz-three-widget-menu button:hover,
      #__ROOT_ID__ .oviz-three-widget-menu select:hover,
      #__ROOT_ID__ .oviz-three-tools-toggle:hover,
      #__ROOT_ID__ .oviz-three-controls-toggle:hover {
        background: transparent !important;
        border-color: transparent !important;
        color: rgba(255, 255, 255, 0.96) !important;
      }
      #__ROOT_ID__ .oviz-three-widget-menu button:focus-visible,
      #__ROOT_ID__ .oviz-three-widget-menu select:focus-visible,
      #__ROOT_ID__ .oviz-three-tools-toggle:focus-visible,
      #__ROOT_ID__ .oviz-three-controls-toggle:focus-visible,
      #__ROOT_ID__ .oviz-three-widget-menu .oviz-three-tools-shell[data-open="true"] > .oviz-three-tools-toggle,
      #__ROOT_ID__ .oviz-three-widget-menu .oviz-three-controls-shell[data-open="true"] > .oviz-three-controls-toggle,
      #__ROOT_ID__ .oviz-three-widget-menu button[data-active="true"],
      #__ROOT_ID__ .oviz-three-widget-menu .oviz-three-auto-orbit[aria-pressed="true"] {
        outline: none !important;
        color: rgba(255, 255, 255, 0.96) !important;
        box-shadow: inset 0 -1px 0 rgba(246, 200, 95, 0.66) !important;
      }
      #__ROOT_ID__ .oviz-three-widget-select {
        min-width: 92px !important;
        padding-right: 17px !important;
        -webkit-appearance: none !important;
        appearance: none !important;
        background:
          linear-gradient(45deg, transparent 50%, rgba(238, 242, 247, 0.58) 50%) right 7px center / 4px 4px no-repeat,
          linear-gradient(135deg, rgba(238, 242, 247, 0.58) 50%, transparent 50%) right 4px center / 4px 4px no-repeat,
          transparent !important;
      }
      #__ROOT_ID__ .oviz-three-widget-menu .oviz-three-tools-drawer,
      #__ROOT_ID__ .oviz-three-widget-menu .oviz-three-controls-drawer {
        top: calc(100% + 8px) !important;
        box-sizing: border-box !important;
        max-height: min(76vh, 720px) !important;
        overflow-y: auto !important;
        overflow-x: visible !important;
        padding: 14px 16px 16px !important;
        border: 1px solid var(--oviz-instrument-border) !important;
        border-radius: 9px !important;
        background: linear-gradient(180deg, rgba(10, 12, 17, 0.74), rgba(6, 8, 12, 0.54)) !important;
        box-shadow: 0 16px 44px rgba(0, 0, 0, 0.28) !important;
        backdrop-filter: blur(14px) saturate(118%) !important;
        -webkit-backdrop-filter: blur(14px) saturate(118%) !important;
        transform: translate3d(0, -7px, 0) scale(0.99) !important;
        transform-origin: top right !important;
        transition: opacity 190ms ease, transform 260ms var(--oviz-instrument-ease) !important;
      }
      #__ROOT_ID__ .oviz-three-widget-menu .oviz-three-tools-drawer {
        width: min(316px, calc(100vw - 32px)) !important;
      }
      #__ROOT_ID__ .oviz-three-widget-menu .oviz-three-controls-drawer {
        width: min(400px, calc(100vw - 32px)) !important;
      }
      #__ROOT_ID__ .oviz-three-widget-menu .oviz-three-tools-shell[data-open="true"] .oviz-three-tools-drawer,
      #__ROOT_ID__ .oviz-three-widget-menu .oviz-three-controls-shell[data-open="true"] .oviz-three-controls-drawer {
        transform: translate3d(0, 0, 0) scale(1) !important;
      }
      #__ROOT_ID__ .oviz-three-widget-menu .oviz-three-selection,
      #__ROOT_ID__ .oviz-three-widget-menu .oviz-three-controls {
        gap: 10px !important;
        padding: 0 !important;
        border: 0 !important;
        border-radius: 0 !important;
        background: transparent !important;
      }
      #__ROOT_ID__ .oviz-three-controls-title,
      #__ROOT_ID__ .oviz-three-controls-hint,
      #__ROOT_ID__ .oviz-three-selection-readout,
      #__ROOT_ID__ .oviz-three-volume-summary,
      #__ROOT_ID__ .oviz-three-legend-summary {
        display: none !important;
      }
      #__ROOT_ID__ .oviz-three-controls-row {
        gap: 10px !important;
      }
      #__ROOT_ID__ .oviz-three-controls-field,
      #__ROOT_ID__ .oviz-three-controls-toggle-row,
      #__ROOT_ID__ .oviz-three-selection-toggle,
      #__ROOT_ID__ .oviz-three-volume-field,
      #__ROOT_ID__ .oviz-three-volume-toggle,
      #__ROOT_ID__ .oviz-three-legend-toggle {
        gap: 5px !important;
        color: var(--oviz-instrument-muted) !important;
        font: 620 11px/1.24 -apple-system, BlinkMacSystemFont, "SF Pro Text", "Helvetica Neue", sans-serif !important;
        letter-spacing: 0 !important;
      }
      #__ROOT_ID__ .oviz-three-controls-field,
      #__ROOT_ID__ .oviz-three-volume-field,
      #__ROOT_ID__ .oviz-three-legend-field,
      #__ROOT_ID__ .oviz-three-controls-row,
      #__ROOT_ID__ .oviz-three-volume-row,
      #__ROOT_ID__ .oviz-three-legend-control-row {
        overflow: visible !important;
      }
      #__ROOT_ID__ .oviz-three-controls-field input[type="number"],
      #__ROOT_ID__ .oviz-three-controls-field input[type="text"],
      #__ROOT_ID__ .oviz-three-controls-field select,
      #__ROOT_ID__ .oviz-three-volume input[type="number"],
      #__ROOT_ID__ .oviz-three-volume select,
      #__ROOT_ID__ .oviz-three-filter-field select,
      #__ROOT_ID__ .oviz-three-legend-field input[type="number"],
      #__ROOT_ID__ .oviz-three-legend-field select {
        min-height: 28px !important;
        padding: 3px 0 4px !important;
        border: 0 !important;
        border-bottom: 1px solid rgba(238, 242, 247, 0.16) !important;
        border-radius: 0 !important;
        background: transparent !important;
        color: rgba(238, 242, 247, 0.84) !important;
        font: 650 11.5px/1.16 -apple-system, BlinkMacSystemFont, "SF Pro Text", "Helvetica Neue", sans-serif !important;
        box-shadow: none !important;
        transition: border-color 170ms ease, color 170ms ease !important;
      }
      #__ROOT_ID__ .oviz-three-controls-field input[type="number"]:focus,
      #__ROOT_ID__ .oviz-three-controls-field input[type="text"]:focus,
      #__ROOT_ID__ .oviz-three-controls-field select:focus,
      #__ROOT_ID__ .oviz-three-volume input[type="number"]:focus,
      #__ROOT_ID__ .oviz-three-volume select:focus,
      #__ROOT_ID__ .oviz-three-filter-field select:focus,
      #__ROOT_ID__ .oviz-three-legend-field input[type="number"]:focus,
      #__ROOT_ID__ .oviz-three-legend-field select:focus {
        outline: none !important;
        border-bottom-color: rgba(246, 200, 95, 0.58) !important;
        color: rgba(255, 255, 255, 0.95) !important;
        box-shadow: none !important;
      }
      #__ROOT_ID__ .oviz-three-controls-toggle-row input[type="checkbox"],
      #__ROOT_ID__ .oviz-three-selection-toggle input[type="checkbox"],
      #__ROOT_ID__ .oviz-three-legend-toggle input[type="checkbox"],
      #__ROOT_ID__ .oviz-three-volume-toggle input[type="checkbox"],
      #__ROOT_ID__ .oviz-three-box-visible {
        position: relative !important;
        width: 12px !important;
        height: 12px !important;
        flex: 0 0 12px !important;
        margin: 0 !important;
        border: 1px solid rgba(238, 242, 247, 0.42) !important;
        border-radius: 3px !important;
        background: rgba(255, 255, 255, 0.025) !important;
        accent-color: var(--oviz-instrument-accent) !important;
        -webkit-appearance: none !important;
        appearance: none !important;
        transition: border-color 160ms ease, background 160ms ease, box-shadow 160ms ease !important;
      }
      #__ROOT_ID__ .oviz-three-controls-toggle-row input[type="checkbox"]:checked,
      #__ROOT_ID__ .oviz-three-selection-toggle input[type="checkbox"]:checked,
      #__ROOT_ID__ .oviz-three-legend-toggle input[type="checkbox"]:checked,
      #__ROOT_ID__ .oviz-three-volume-toggle input[type="checkbox"]:checked,
      #__ROOT_ID__ .oviz-three-box-visible:checked {
        border-color: rgba(246, 200, 95, 0.82) !important;
        background: rgba(246, 200, 95, 0.82) !important;
      }
      #__ROOT_ID__ .oviz-three-controls-toggle-row input[type="checkbox"]:checked::after,
      #__ROOT_ID__ .oviz-three-selection-toggle input[type="checkbox"]:checked::after,
      #__ROOT_ID__ .oviz-three-legend-toggle input[type="checkbox"]:checked::after,
      #__ROOT_ID__ .oviz-three-volume-toggle input[type="checkbox"]:checked::after,
      #__ROOT_ID__ .oviz-three-box-visible:checked::after {
        content: "" !important;
        position: absolute !important;
        left: 3px !important;
        top: 1px !important;
        width: 4px !important;
        height: 7px !important;
        border-right: 2px solid rgba(0, 0, 0, 0.78) !important;
        border-bottom: 2px solid rgba(0, 0, 0, 0.78) !important;
        transform: rotate(40deg) !important;
      }
      #__ROOT_ID__ .oviz-three-controls-toggle-row input[type="checkbox"]:focus-visible,
      #__ROOT_ID__ .oviz-three-selection-toggle input[type="checkbox"]:focus-visible,
      #__ROOT_ID__ .oviz-three-legend-toggle input[type="checkbox"]:focus-visible,
      #__ROOT_ID__ .oviz-three-volume-toggle input[type="checkbox"]:focus-visible,
      #__ROOT_ID__ .oviz-three-box-visible:focus-visible {
        outline: none !important;
        box-shadow: 0 0 0 3px rgba(246, 200, 95, 0.14) !important;
      }
      #__ROOT_ID__ .oviz-three-controls-field input[type="range"],
      #__ROOT_ID__ .oviz-three-volume input[type="range"],
      #__ROOT_ID__ .oviz-three-filter-range,
      #__ROOT_ID__ .oviz-three-age-filter-range,
      #__ROOT_ID__ .oviz-three-legend-field input[type="range"],
      #__ROOT_ID__ .oviz-three-slider {
        height: 30px !important;
        min-height: 30px !important;
        overflow: visible !important;
        padding: 0 !important;
        border: 0 !important;
        background: transparent !important;
        accent-color: var(--oviz-instrument-accent) !important;
        -webkit-appearance: none !important;
        appearance: none !important;
        cursor: pointer !important;
      }
      #__ROOT_ID__ .oviz-three-controls-field input[type="range"]::-webkit-slider-runnable-track,
      #__ROOT_ID__ .oviz-three-volume input[type="range"]::-webkit-slider-runnable-track,
      #__ROOT_ID__ .oviz-three-filter-range::-webkit-slider-runnable-track,
      #__ROOT_ID__ .oviz-three-age-filter-range::-webkit-slider-runnable-track,
      #__ROOT_ID__ .oviz-three-legend-field input[type="range"]::-webkit-slider-runnable-track,
      #__ROOT_ID__ .oviz-three-slider::-webkit-slider-runnable-track,
      #__ROOT_ID__[data-minimal="true"][data-galactic-simple="true"] .oviz-three-slider::-webkit-slider-runnable-track {
        height: 2px !important;
        border: 0 !important;
        border-radius: 999px !important;
        background: rgba(238, 242, 247, 0.22) !important;
      }
      #__ROOT_ID__ .oviz-three-controls-field input[type="range"]::-webkit-slider-thumb,
      #__ROOT_ID__ .oviz-three-volume input[type="range"]::-webkit-slider-thumb,
      #__ROOT_ID__ .oviz-three-filter-range::-webkit-slider-thumb,
      #__ROOT_ID__ .oviz-three-age-filter-range::-webkit-slider-thumb,
      #__ROOT_ID__ .oviz-three-legend-field input[type="range"]::-webkit-slider-thumb,
      #__ROOT_ID__ .oviz-three-slider::-webkit-slider-thumb,
      #__ROOT_ID__[data-minimal="true"][data-galactic-simple="true"] .oviz-three-slider::-webkit-slider-thumb {
        -webkit-appearance: none !important;
        appearance: none !important;
        width: 12px !important;
        height: 12px !important;
        margin-top: -5px !important;
        border: 1px solid rgba(255, 255, 255, 0.66) !important;
        border-radius: 50% !important;
        background: rgba(246, 200, 95, 0.92) !important;
        box-shadow: 0 0 0 4px rgba(246, 200, 95, 0.055) !important;
      }
      #__ROOT_ID__ .oviz-three-controls-field input[type="range"]::-moz-range-track,
      #__ROOT_ID__ .oviz-three-volume input[type="range"]::-moz-range-track,
      #__ROOT_ID__ .oviz-three-filter-range::-moz-range-track,
      #__ROOT_ID__ .oviz-three-age-filter-range::-moz-range-track,
      #__ROOT_ID__ .oviz-three-legend-field input[type="range"]::-moz-range-track,
      #__ROOT_ID__ .oviz-three-slider::-moz-range-track,
      #__ROOT_ID__[data-minimal="true"][data-galactic-simple="true"] .oviz-three-slider::-moz-range-track {
        height: 2px !important;
        border: 0 !important;
        border-radius: 999px !important;
        background: rgba(238, 242, 247, 0.22) !important;
      }
      #__ROOT_ID__ .oviz-three-controls-field input[type="range"]::-moz-range-thumb,
      #__ROOT_ID__ .oviz-three-volume input[type="range"]::-moz-range-thumb,
      #__ROOT_ID__ .oviz-three-filter-range::-moz-range-thumb,
      #__ROOT_ID__ .oviz-three-age-filter-range::-moz-range-thumb,
      #__ROOT_ID__ .oviz-three-legend-field input[type="range"]::-moz-range-thumb,
      #__ROOT_ID__ .oviz-three-slider::-moz-range-thumb,
      #__ROOT_ID__[data-minimal="true"][data-galactic-simple="true"] .oviz-three-slider::-moz-range-thumb {
        width: 12px !important;
        height: 12px !important;
        border: 1px solid rgba(255, 255, 255, 0.66) !important;
        border-radius: 50% !important;
        background: rgba(246, 200, 95, 0.92) !important;
        box-shadow: 0 0 0 4px rgba(246, 200, 95, 0.055) !important;
      }
      #__ROOT_ID__ .oviz-three-filter-slider-shell,
      #__ROOT_ID__ .oviz-three-age-filter-slider-shell {
        height: 38px !important;
        margin-top: 6px !important;
        padding: 2px 0 4px !important;
        overflow: visible !important;
      }
      #__ROOT_ID__ .oviz-three-filter-slider-shell::before,
      #__ROOT_ID__ .oviz-three-age-filter-slider-shell::before {
        left: 10px !important;
        right: 10px !important;
        top: 18px !important;
        height: 2px !important;
      }
      #__ROOT_ID__ .oviz-three-filter-range,
      #__ROOT_ID__ .oviz-three-age-filter-range {
        position: absolute !important;
        inset: 0 !important;
        height: 38px !important;
        min-height: 38px !important;
        margin: 0 !important;
        pointer-events: none !important;
      }
      #__ROOT_ID__ .oviz-three-filter-range::-webkit-slider-runnable-track,
      #__ROOT_ID__ .oviz-three-age-filter-range::-webkit-slider-runnable-track {
        height: 38px !important;
        background: transparent !important;
      }
      #__ROOT_ID__ .oviz-three-filter-range::-webkit-slider-thumb,
      #__ROOT_ID__ .oviz-three-age-filter-range::-webkit-slider-thumb {
        width: 14px !important;
        height: 14px !important;
        margin-top: 12px !important;
        pointer-events: auto !important;
      }
      #__ROOT_ID__ .oviz-three-filter-range::-moz-range-track,
      #__ROOT_ID__ .oviz-three-age-filter-range::-moz-range-track {
        height: 38px !important;
        background: transparent !important;
      }
      #__ROOT_ID__ .oviz-three-filter-range::-moz-range-thumb,
      #__ROOT_ID__ .oviz-three-age-filter-range::-moz-range-thumb {
        width: 14px !important;
        height: 14px !important;
        pointer-events: auto !important;
      }
      #__ROOT_ID__ .oviz-three-selection button,
      #__ROOT_ID__ .oviz-three-controls-actions button,
      #__ROOT_ID__ .oviz-three-action-button {
        min-height: 28px !important;
        padding: 0 2px !important;
        border: 0 !important;
        border-radius: 0 !important;
        background: transparent !important;
        color: rgba(238, 242, 247, 0.64) !important;
        box-shadow: none !important;
        font: 680 11px/1.2 -apple-system, BlinkMacSystemFont, "SF Pro Text", "Helvetica Neue", sans-serif !important;
        transition: color 160ms ease, box-shadow 180ms ease !important;
      }
      #__ROOT_ID__ .oviz-three-selection button:hover,
      #__ROOT_ID__ .oviz-three-controls-actions button:hover,
      #__ROOT_ID__ .oviz-three-action-button:hover,
      #__ROOT_ID__ .oviz-three-selection button[data-active="true"],
      #__ROOT_ID__ .oviz-three-controls-actions button[data-active="true"],
      #__ROOT_ID__ .oviz-three-action-button[data-active="true"] {
        background: transparent !important;
        border-color: transparent !important;
        color: rgba(255, 255, 255, 0.94) !important;
        box-shadow: inset 0 -1px 0 rgba(246, 200, 95, 0.58) !important;
      }
      #__ROOT_ID__ .oviz-three-widget-panel,
      #__ROOT_ID__ .oviz-three-key-help {
        border: 1px solid var(--oviz-instrument-border) !important;
        border-radius: 10px !important;
        background: linear-gradient(180deg, rgba(10, 12, 17, 0.74), rgba(6, 8, 12, 0.58)) !important;
        box-shadow: 0 18px 54px rgba(0, 0, 0, 0.30) !important;
        backdrop-filter: blur(14px) saturate(120%) !important;
        -webkit-backdrop-filter: blur(14px) saturate(120%) !important;
      }
      #__ROOT_ID__ .oviz-three-age-panel,
      #__ROOT_ID__ .oviz-three-filter-panel {
        min-width: min(360px, calc(100vw - 24px)) !important;
        min-height: min(360px, calc(100vh - 24px)) !important;
      }
      #__ROOT_ID__ .oviz-three-age-body,
      #__ROOT_ID__ .oviz-three-filter-body {
        gap: 10px !important;
        padding: 12px 12px 18px !important;
      }
      #__ROOT_ID__ .oviz-three-age-filter,
      #__ROOT_ID__ .oviz-three-filter-slider-shell {
        padding-bottom: 4px !important;
        overflow: visible !important;
      }
      #__ROOT_ID__ .oviz-three-widget-drag {
        height: 32px !important;
        padding: 0 12px !important;
        border-bottom: 1px solid rgba(238, 242, 247, 0.08) !important;
        background: rgba(255, 255, 255, 0.018) !important;
        color: rgba(255, 255, 255, 0.88) !important;
        font: 760 11.5px/1 -apple-system, BlinkMacSystemFont, "SF Pro Text", "Helvetica Neue", sans-serif !important;
        letter-spacing: 0 !important;
      }
      #__ROOT_ID__ .oviz-three-widget-panel button:not(.oviz-three-window-button) {
        min-height: 28px !important;
        border-radius: 6px !important;
        font: 650 11px/1.2 -apple-system, BlinkMacSystemFont, "SF Pro Text", "Helvetica Neue", sans-serif !important;
      }
      #__ROOT_ID__ .oviz-three-window-button {
        width: 22px !important;
        height: 22px !important;
        border-radius: 999px !important;
      }
      #__ROOT_ID__ .oviz-three-footer {
        gap: 10px !important;
        padding: 0 !important;
        background: transparent !important;
      }
      #__ROOT_ID__ .oviz-three-footer button {
        width: 32px !important;
        height: 32px !important;
        border: 1px solid transparent !important;
        border-radius: 999px !important;
        background: transparent !important;
        color: rgba(238, 242, 247, 0.78) !important;
        transition: color 160ms ease, background 180ms ease, border-color 180ms ease !important;
      }
      #__ROOT_ID__ .oviz-three-footer button:hover,
      #__ROOT_ID__ .oviz-three-footer button[data-active="true"] {
        color: rgba(255, 255, 255, 0.96) !important;
        border-color: rgba(246, 200, 95, 0.22) !important;
        background: rgba(246, 200, 95, 0.055) !important;
      }
      #__ROOT_ID__ .oviz-three-time-label {
        color: rgba(255, 255, 255, 0.94) !important;
        font: 760 14px/1.15 -apple-system, BlinkMacSystemFont, "SF Pro Text", "Helvetica Neue", sans-serif !important;
        text-shadow: 0 0 10px rgba(0, 0, 0, 0.82) !important;
      }
      #__ROOT_ID__ .oviz-three-slider-shell {
        height: 42px !important;
      }
      #__ROOT_ID__ .oviz-three-slider-track-wrap {
        height: 30px !important;
        overflow: visible !important;
      }
      #__ROOT_ID__ .oviz-three-scale-bar {
        gap: 5px !important;
        padding: 7px 10px !important;
        border: 1px solid rgba(238, 242, 247, 0.08) !important;
        border-radius: 8px !important;
        background: rgba(6, 8, 12, 0.38) !important;
        box-shadow: 0 10px 30px rgba(0, 0, 0, 0.16) !important;
        backdrop-filter: blur(10px) saturate(116%) !important;
        -webkit-backdrop-filter: blur(10px) saturate(116%) !important;
        cursor: default !important;
        pointer-events: none !important;
      }
      #__ROOT_ID__ .oviz-three-scale-label {
        color: rgba(255, 255, 255, 0.9) !important;
        font: 760 11.5px/1.1 -apple-system, BlinkMacSystemFont, "SF Pro Text", "Helvetica Neue", sans-serif !important;
      }
      #__ROOT_ID__ .oviz-three-scale-line::before {
        border-top-color: rgba(238, 242, 247, 0.66) !important;
        border-top-width: 1px !important;
      }
      #__ROOT_ID__ .oviz-three-scale-line::after {
        border-left-color: rgba(238, 242, 247, 0.66) !important;
        border-right-color: rgba(238, 242, 247, 0.66) !important;
        border-left-width: 1px !important;
        border-right-width: 1px !important;
      }
      #__ROOT_ID__ .oviz-three-legend-name {
        max-inline-size: min(38vw, 360px) !important;
        white-space: normal !important;
        overflow-wrap: anywhere !important;
        font-size: 15px !important;
        line-height: 1.24 !important;
      }
      #__ROOT_ID__ .oviz-three-legend-row {
        grid-template-columns: 18px minmax(0, auto) !important;
        column-gap: 5px !important;
        align-items: center !important;
      }
      #__ROOT_ID__ .oviz-three-legend-edit {
        width: 21px !important;
        height: 21px !important;
        min-width: 21px !important;
        min-height: 21px !important;
        border: 0 !important;
        border-radius: 0 !important;
        color: rgba(255, 255, 255, 0.42) !important;
        background: transparent !important;
        font: 720 15px/1 -apple-system, BlinkMacSystemFont, "SF Pro Text", "Helvetica Neue", sans-serif !important;
        box-shadow: none !important;
      }
      #__ROOT_ID__ .oviz-three-legend-edit:hover,
      #__ROOT_ID__ .oviz-three-legend-edit:focus-visible,
      #__ROOT_ID__ .oviz-three-legend-edit[data-open="true"] {
        outline: none !important;
        color: rgba(255, 255, 255, 0.88) !important;
        box-shadow: inset 0 -1px 0 rgba(246, 200, 95, 0.58) !important;
      }
      #__ROOT_ID__ .oviz-three-legend-controls {
        display: grid !important;
        gap: 8px !important;
        width: min(292px, calc(100vw - 36px)) !important;
        max-width: min(292px, calc(100vw - 36px)) !important;
        max-height: 0 !important;
        margin: 0 !important;
        padding: 0 8px 0 18px !important;
        overflow: hidden !important;
        opacity: 0 !important;
        transform: translateY(-8px) !important;
        transform-origin: top left !important;
        transition:
          max-height 300ms var(--oviz-instrument-ease),
          margin 240ms var(--oviz-instrument-ease),
          padding 240ms var(--oviz-instrument-ease),
          opacity 200ms ease,
          transform 260ms var(--oviz-instrument-ease) !important;
      }
      #__ROOT_ID__ .oviz-three-legend-controls[data-visible="true"] {
        max-height: min(72vh, 680px) !important;
        margin: 5px 0 10px 0 !important;
        padding: 6px 8px 10px 18px !important;
        overflow-y: visible !important;
        overflow-x: visible !important;
        opacity: 1 !important;
        transform: translateY(0) !important;
      }
      #__ROOT_ID__ .oviz-three-legend-control-row {
        display: grid !important;
        grid-template-columns: 1fr !important;
        gap: 8px !important;
        overflow: visible !important;
      }
      #__ROOT_ID__ .oviz-three-legend-field {
        display: grid !important;
        grid-template-columns: minmax(70px, max-content) 154px !important;
        align-items: center !important;
        gap: 10px !important;
        min-height: 32px !important;
        overflow: visible !important;
        color: rgba(238, 242, 247, 0.50) !important;
        font: 620 10.5px/1.18 -apple-system, BlinkMacSystemFont, "SF Pro Text", "Helvetica Neue", sans-serif !important;
      }
      #__ROOT_ID__ .oviz-three-legend-field input[type="range"] {
        width: 140px !important;
        min-width: 140px !important;
        max-width: 140px !important;
        height: 34px !important;
        min-height: 34px !important;
        margin: 0 7px !important;
        overflow: visible !important;
      }
      #__ROOT_ID__ .oviz-three-legend-field input[type="range"]::-webkit-slider-runnable-track {
        height: 2px !important;
      }
      #__ROOT_ID__ .oviz-three-legend-field input[type="range"]::-webkit-slider-thumb {
        width: 12px !important;
        height: 12px !important;
        margin-top: -5px !important;
        box-shadow: 0 0 0 4px rgba(246, 200, 95, 0.06) !important;
      }
      #__ROOT_ID__ .oviz-three-legend-field input[type="range"]::-moz-range-track {
        height: 2px !important;
      }
      #__ROOT_ID__ .oviz-three-legend-field input[type="range"]::-moz-range-thumb {
        width: 12px !important;
        height: 12px !important;
        box-shadow: 0 0 0 4px rgba(246, 200, 95, 0.06) !important;
      }
      #__ROOT_ID__ .oviz-three-legend-field input[type="number"],
      #__ROOT_ID__ .oviz-three-legend-field select {
        width: 154px !important;
        max-width: 154px !important;
        min-height: 30px !important;
      }
      #__ROOT_ID__ .oviz-three-legend-field input[type="color"] {
        width: 44px !important;
        height: 24px !important;
        padding: 0 !important;
        border: 1px solid rgba(238, 242, 247, 0.18) !important;
        border-radius: 2px !important;
        background: transparent !important;
      }
      #__ROOT_ID__ .oviz-three-group-dropdown[data-open="true"] .oviz-three-group-menu {
        max-height: 282px !important;
        overflow-y: auto !important;
        overflow-x: hidden !important;
      }
      #__ROOT_ID__ .oviz-three-group-menu-list {
        max-height: 282px !important;
      }
      #__ROOT_ID__ .oviz-three-sky-controls-shell {
        position: relative;
        min-height: auto;
        border: 0;
        border-radius: 0;
        background: transparent;
        box-shadow: none;
        backdrop-filter: none;
        -webkit-backdrop-filter: none;
        overflow: visible;
      }
      #__ROOT_ID__ .oviz-three-sky-controls-toggle {
        width: auto;
        height: 30px;
        padding: 0 12px;
        justify-content: center;
        text-align: center;
        border-radius: 4px;
        border: 1px solid rgba(255, 255, 255, 0.06);
        background: rgba(29, 32, 37, 0.98);
        color: var(--oviz-text);
        cursor: pointer;
        box-shadow: none;
        font: 500 11px/1 -apple-system, BlinkMacSystemFont, "SF Pro Text", "Helvetica Neue", sans-serif;
      }
      #__ROOT_ID__ .oviz-three-sky-controls-toggle:hover,
      #__ROOT_ID__ .oviz-three-sky-controls-toggle:focus-visible,
      #__ROOT_ID__ .oviz-three-sky-controls-shell[data-open="true"] > .oviz-three-sky-controls-toggle {
        background: rgba(38, 42, 48, 0.98);
        border-color: rgba(255, 255, 255, 0.10);
        outline: none;
      }
      #__ROOT_ID__ .oviz-three-sky-controls-drawer {
        position: absolute;
        top: calc(100% + 10px);
        right: 0;
        left: auto;
        width: min(340px, calc(100vw - 32px));
        max-height: min(70vh, 620px);
        overflow: auto;
        padding: 12px;
        border: 1px solid var(--oviz-instrument-border);
        border-radius: 8px;
        background: var(--oviz-instrument-bg-strong);
        box-shadow: var(--oviz-shadow-md);
        backdrop-filter: blur(14px) saturate(112%);
        -webkit-backdrop-filter: blur(14px) saturate(112%);
        opacity: 0;
        pointer-events: none;
        transform: translateY(-4px);
        transition: opacity 0.16s ease, transform 0.16s ease;
      }
      #__ROOT_ID__ .oviz-three-sky-controls-shell[data-open="true"] .oviz-three-sky-controls-drawer {
        opacity: 1;
        pointer-events: auto;
        transform: translateY(0);
      }
      #__ROOT_ID__ .oviz-three-sky-source-field[hidden] {
        display: none !important;
      }
      #__ROOT_ID__ .oviz-three-sky-controls {
        gap: 10px !important;
      }
      #__ROOT_ID__ .oviz-three-sky-dome-controls {
        gap: 10px !important;
        padding-top: 0 !important;
        border-top: 0 !important;
      }
      #__ROOT_ID__ .oviz-three-sky-add-grid {
        gap: 8px !important;
      }
      #__ROOT_ID__ .oviz-three-sky-layer-row {
        border: 1px solid var(--oviz-instrument-border) !important;
        border-radius: 8px !important;
        background: rgba(255, 255, 255, 0.035) !important;
      }
      #__ROOT_ID__ .oviz-three-sky-layer-row[open],
      #__ROOT_ID__ .oviz-three-sky-layer-row[data-active="true"] {
        border-color: rgba(246, 200, 95, 0.24) !important;
        background: rgba(246, 200, 95, 0.045) !important;
      }
      #__ROOT_ID__ .oviz-three-sky-layer-summary {
        min-height: 34px !important;
        padding: 8px 9px !important;
      }
      #__ROOT_ID__ .oviz-three-sky-layer-name {
        color: rgba(238, 242, 247, 0.92) !important;
        font: 700 12px/1.18 -apple-system, BlinkMacSystemFont, "SF Pro Text", "Helvetica Neue", sans-serif !important;
      }
      #__ROOT_ID__ .oviz-three-sky-layer-order,
      #__ROOT_ID__ .oviz-three-sky-menu-heading {
        color: rgba(238, 242, 247, 0.54) !important;
        font: 650 10.5px/1.2 -apple-system, BlinkMacSystemFont, "SF Pro Text", "Helvetica Neue", sans-serif !important;
      }
      #__ROOT_ID__ .oviz-three-sky-layer-body {
        gap: 8px !important;
        padding: 0 9px 9px !important;
      }
      #__ROOT_ID__ .oviz-three-sky-layer-body label {
        color: rgba(238, 242, 247, 0.58) !important;
        font: 620 10.5px/1.2 -apple-system, BlinkMacSystemFont, "SF Pro Text", "Helvetica Neue", sans-serif !important;
      }
      #__ROOT_ID__ .oviz-three-sky-layer-body input,
      #__ROOT_ID__ .oviz-three-sky-layer-body select {
        border: 1px solid rgba(238, 242, 247, 0.12) !important;
        border-radius: 6px !important;
        background: rgba(0, 0, 0, 0.18) !important;
        color: var(--oviz-text) !important;
      }
      #__ROOT_ID__ .oviz-three-sky-layer-body input[type="range"] {
        border: 0 !important;
        background: transparent !important;
        accent-color: rgba(246, 200, 95, 0.88) !important;
      }
      #__ROOT_ID__ .oviz-three-sky-layer-actions button {
        border-color: rgba(238, 242, 247, 0.12) !important;
        border-radius: 6px !important;
        background: rgba(255, 255, 255, 0.035) !important;
        font: 650 11px/1 -apple-system, BlinkMacSystemFont, "SF Pro Text", "Helvetica Neue", sans-serif !important;
      }
      #__ROOT_ID__ .oviz-three-group-button:focus-visible,
      #__ROOT_ID__ .oviz-three-group-option:focus-visible {
        outline: none !important;
        box-shadow: inset 0 -1px 0 rgba(246, 200, 95, 0.62) !important;
      }
      @media (max-width: 980px) {
        #__ROOT_ID__ .oviz-three-topbar {
          grid-template-columns: 1fr;
          grid-template-areas:
            "title"
            "actions";
          align-items: start;
          gap: 8px;
        }
        #__ROOT_ID__ .oviz-three-title {
          grid-area: title;
          justify-self: center;
          max-width: calc(100vw - 28px);
        }
        #__ROOT_ID__ .oviz-three-widget-menu {
          grid-area: actions;
          justify-self: end;
          justify-content: flex-end;
          width: min(calc(100vw - 28px), 920px);
          max-width: calc(100vw - 28px);
        }
        #__ROOT_ID__ .oviz-three-legend-panel {
          width: min(244px, calc(100vw - 28px));
        }
        #__ROOT_ID__ .oviz-three-footer {
          width: min(calc(100vw - 28px), 920px);
        }
      }
      @media (max-width: 860px) {
        #__ROOT_ID__ .oviz-three-widget-menu {
          row-gap: 4px;
          justify-self: end;
          justify-content: flex-end;
        }
        #__ROOT_ID__ .oviz-three-widget-menu button,
        #__ROOT_ID__ .oviz-three-widget-menu select,
        #__ROOT_ID__ .oviz-three-tools-toggle,
        #__ROOT_ID__ .oviz-three-controls-toggle {
          height: 30px;
          font-size: 10.5px;
        }
        #__ROOT_ID__ .oviz-three-widget-menu .oviz-three-tools-toggle,
        #__ROOT_ID__ .oviz-three-widget-menu .oviz-three-controls-toggle {
          padding: 0 10px;
        }
        #__ROOT_ID__ .oviz-three-widget-select,
        #__ROOT_ID__ .oviz-three-group-select {
          min-width: 138px;
        }
        #__ROOT_ID__ .oviz-three-footer {
          width: min(calc(100vw - 24px), 760px);
        }
      }
      @media (max-width: 720px) {
        #__ROOT_ID__ .oviz-three-title {
          font-size: 12px;
          padding: 5px 9px;
        }
        #__ROOT_ID__ .oviz-three-widget-menu {
          width: min(calc(100vw - 20px), 640px);
          max-width: calc(100vw - 20px);
          padding: 4px;
          gap: 4px;
          justify-self: end;
          justify-content: flex-end;
        }
        #__ROOT_ID__ .oviz-three-widget-menu button,
        #__ROOT_ID__ .oviz-three-widget-menu select,
        #__ROOT_ID__ .oviz-three-tools-toggle,
        #__ROOT_ID__ .oviz-three-controls-toggle {
          height: 28px;
          font-size: 10px;
        }
        #__ROOT_ID__ .oviz-three-widget-select,
        #__ROOT_ID__ .oviz-three-group-select {
          min-width: 120px;
        }
        #__ROOT_ID__ .oviz-three-legend-panel {
          top: 74px;
          left: 12px;
          width: min(220px, calc(100vw - 24px));
        }
        #__ROOT_ID__ .oviz-three-footer {
          bottom: 12px;
          width: min(calc(100vw - 20px), 620px);
        }
        #__ROOT_ID__ .oviz-three-earth-view-toggle {
          right: 12px;
          bottom: 68px;
          min-width: 96px;
          height: 32px;
          padding: 0 12px;
          font-size: 11px;
        }
        #__ROOT_ID__ .oviz-three-scale-bar {
          display: none;
        }
      }
      @media (max-width: 1100px), (max-height: 760px) {
        #__ROOT_ID__[data-minimal="true"][data-galactic-simple="true"] .oviz-three-footer {
          left: 29%;
          bottom: 14px;
          max-width: calc(100vw - 28px);
        }
      }
      @media (max-width: 560px) {
        #__ROOT_ID__ .oviz-three-footer {
          display: grid;
          grid-template-columns: 32px 32px auto minmax(0, 1fr);
          align-items: center;
          gap: 8px;
          width: min(calc(100vw - 16px), 520px);
        }
        #__ROOT_ID__ .oviz-three-footer button {
          width: 32px;
          height: 32px;
        }
        #__ROOT_ID__ .oviz-three-time-label {
          min-width: 96px;
          font-size: 12px;
        }
        #__ROOT_ID__ .oviz-three-slider-shell {
          width: auto;
          max-width: none;
          min-width: 0;
        }
        #__ROOT_ID__ .oviz-three-legend-panel {
          width: min(206px, calc(100vw - 20px));
        }
      }
      @media (max-height: 760px) {
        #__ROOT_ID__ .oviz-three-topbar {
          top: 10px;
          left: 10px;
          right: 10px;
        }
        #__ROOT_ID__ .oviz-three-widget-menu {
          max-width: min(calc(100vw - 20px), 760px);
        }
        #__ROOT_ID__ .oviz-three-widget-menu .oviz-three-tools-drawer,
        #__ROOT_ID__ .oviz-three-widget-menu .oviz-three-controls-drawer {
          max-height: min(62vh, 560px);
        }
        #__ROOT_ID__ .oviz-three-legend-panel {
          top: 10px;
        }
        #__ROOT_ID__ .oviz-three-footer {
          bottom: 10px;
        }
        #__ROOT_ID__ .oviz-three-note {
          display: none;
        }
      }
      @media (max-height: 620px) {
        #__ROOT_ID__ .oviz-three-title {
          font-size: 12px;
          padding: 5px 9px;
        }
        #__ROOT_ID__ .oviz-three-widget-menu button,
        #__ROOT_ID__ .oviz-three-widget-menu select,
        #__ROOT_ID__ .oviz-three-tools-toggle,
        #__ROOT_ID__ .oviz-three-controls-toggle {
          height: 28px;
          font-size: 10px;
        }
        #__ROOT_ID__ .oviz-three-legend-panel {
          width: min(212px, calc(100vw - 20px));
        }
        #__ROOT_ID__ .oviz-three-footer {
          bottom: 8px;
        }
      }
    </style>
  </head>
  <body>
    <div id="__ROOT_ID__" tabindex="0" data-zen="false" __ROOT_MINIMAL_ATTR__ __ROOT_GALACTIC_SIMPLE_ATTR__>
      __TOPBAR_HTML__
      <button class="oviz-three-earth-view-toggle" type="button" title="Move to the sky view from the observer position" aria-pressed="false">3D View</button>
      __SHELL_HTML__
    </div>
    __SCENE_SPEC_PAYLOAD_HTML__

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

      const root = document.getElementById("__ROOT_ID__");
      /*__SCENE_SPEC_METADATA_START__*/const sceneSpecPayloadMetadata = __SCENE_SPEC_PAYLOAD_METADATA__;/*__SCENE_SPEC_METADATA_END__*/
      window.__OVIZ_SCENE_SPEC_PAYLOAD_METADATA__ = sceneSpecPayloadMetadata;

      function attachSceneSpecPayloadMetadata(loadedSceneSpec) {
        if (!loadedSceneSpec || typeof loadedSceneSpec !== "object") {
          return loadedSceneSpec;
        }
        const existingMetadata = loadedSceneSpec.export_metadata && typeof loadedSceneSpec.export_metadata === "object"
          ? loadedSceneSpec.export_metadata
          : {};
        existingMetadata.scene_spec_payload = sceneSpecPayloadMetadata;
        loadedSceneSpec.export_metadata = existingMetadata;
        return loadedSceneSpec;
      }

      function decodeOvizBase64ToUint8Array(encoded) {
        const binary = window.atob(String(encoded || "").replace(/\s+/g, ""));
        const bytes = new Uint8Array(binary.length);
        for (let idx = 0; idx < binary.length; idx += 1) {
          bytes[idx] = binary.charCodeAt(idx);
        }
        return bytes;
      }

      function readOvizSceneSpecPayload(payloadId) {
        const payloadEl = document.getElementById(payloadId);
        if (!payloadEl) {
          renderError(
            root,
            "oviz could not find the embedded compressed scene payload in this HTML file."
          );
          return "";
        }
        return payloadEl.textContent || "";
      }

      async function inflateOvizGzipBase64SceneSpec(encoded) {
        if (typeof DecompressionStream === "undefined") {
          renderError(
            root,
            "This oviz Three.js export uses a compressed scene payload, but this browser does not support DecompressionStream('gzip').\\n\\n"
            + "Open the file in a current browser with gzip DecompressionStream support, or regenerate it with compress_scene_spec=False."
          );
          return null;
        }

        try {
          const compressedBytes = decodeOvizBase64ToUint8Array(encoded);
          const stream = new Blob([compressedBytes])
            .stream()
            .pipeThrough(new DecompressionStream("gzip"));
          const buffer = await new Response(stream).arrayBuffer();
          const jsonText = new TextDecoder().decode(buffer);
          return JSON.parse(jsonText);
        } catch (err) {
          renderError(
            root,
            "oviz could not load the compressed scene payload.\\n\\n"
            + "The embedded gzip/base64 sceneSpec data may be incomplete or unsupported by this browser.\\n\\n"
            + String(err && err.message ? err.message : err)
          );
          return null;
        }
      }

      /*__SCENE_SPEC_START__*/const sceneSpec = __SCENE_SPEC_EXPR__;/*__SCENE_SPEC_END__*/
      if (!sceneSpec) {
        return;
      }
      attachSceneSpecPayloadMetadata(sceneSpec);
      const initialState = sceneSpec.initial_state || {};
      const minimalModeEnabled = Boolean(
        initialState.lite_mode_enabled
        || initialState.minimal_mode_enabled
        || sceneSpec.export_profile === "lite"
        || sceneSpec.export_profile === "minimal"
      );
      const galacticSimpleSpec = sceneSpec.galactic_simple && typeof sceneSpec.galactic_simple === "object"
        ? sceneSpec.galactic_simple
        : {};
      const galacticSimpleModeEnabled = Boolean(galacticSimpleSpec.enabled);
      const galacticSimpleTracksOrbitTargetToSun = Boolean(galacticSimpleSpec.track_orbit_target_to_sun);
      root.dataset.minimal = minimalModeEnabled ? "true" : "false";
      root.dataset.galacticSimple = galacticSimpleModeEnabled ? "true" : "false";
      const initialZoomAnchorSpec = initialState.initial_zoom_anchor && typeof initialState.initial_zoom_anchor === "object"
        ? initialState.initial_zoom_anchor
        : null;
      const canvas = root.querySelector(".oviz-three-canvas");
      const actionBarEl = root.querySelector(".oviz-three-action-bar");
      const titleEl = root.querySelector(".oviz-three-title");
      const widgetMenuEl = root.querySelector(".oviz-three-widget-menu");
      const zenModeButtonEl = root.querySelector(".oviz-three-zen-mode");
      const resetViewButtonEl = root.querySelector(".oviz-three-reset-view");
      const saveStateButtonEl = root.querySelector(".oviz-three-save-state");
      const legendGroupFieldEl = root.querySelector(".oviz-three-legend-group-field");
      const groupDropdownEl = root.querySelector(".oviz-three-group-dropdown");
      const groupDropdownTriggerEl = root.querySelector(".oviz-three-group-trigger");
      const groupDropdownLabelEl = root.querySelector(".oviz-three-group-current");
      const groupDropdownListEl = root.querySelector(".oviz-three-group-menu-list");
      const groupSelectEl = root.querySelector(".oviz-three-group-select");
      const widgetSelectEl = root.querySelector(".oviz-three-widget-select");
      const legendPanelEl = root.querySelector(".oviz-three-legend-panel");
      const legendPanelToggleEl = root.querySelector(".oviz-three-legend-panel-toggle");
      const legendPanelBodyEl = root.querySelector(".oviz-three-legend-panel-body");
      const legendTraceSectionEl = root.querySelector(".oviz-three-legend-trace-section");
      const legendVolumeSectionEl = root.querySelector(".oviz-three-legend-volume-section");
      const legendTraceSectionToggleEl = root.querySelector(".oviz-three-legend-trace-section-toggle");
      const legendVolumeSectionToggleEl = root.querySelector(".oviz-three-legend-volume-section-toggle");
      const legendTraceListEl = root.querySelector(".oviz-three-legend-trace-list");
      const legendVolumeListEl = root.querySelector(".oviz-three-legend-volume-list");
      const legendDragHandleEl = root.querySelector(".oviz-three-legend-panel-drag");
      const legendResizeEls = Array.from(root.querySelectorAll(".oviz-three-legend-resize"));
      const legendPopoverEl = root.querySelector(".oviz-three-legend-popover");
      const keyHelpEl = root.querySelector(".oviz-three-key-help");
      const keyHelpButtonEl = root.querySelector(".oviz-three-key-help-button");
      const keyHelpCloseEl = root.querySelector(".oviz-three-key-help-close");
      const toolsShellEl = root.querySelector(".oviz-three-tools-shell");
      const toolsToggleEl = root.querySelector(".oviz-three-tools-toggle");
      const controlsShellEl = root.querySelector(".oviz-three-controls-shell");
      const controlsToggleEl = root.querySelector(".oviz-three-controls-toggle");
      const skyControlsShellEl = root.querySelector(".oviz-three-sky-controls-shell");
      const skyControlsToggleEl = root.querySelector(".oviz-three-sky-controls-toggle");
      const sliderEl = root.querySelector(".oviz-three-slider");
      const sliderTrackWrapEl = root.querySelector(".oviz-three-slider-track-wrap");
      const sliderMinorTicksEl = root.querySelector(".oviz-three-slider-ticks-minor");
      const sliderMajorTicksEl = root.querySelector(".oviz-three-slider-ticks-major");
      const sliderLabelsEl = root.querySelector(".oviz-three-slider-labels");
      const timeLabelEl = root.querySelector(".oviz-three-time-label");
      const playBackwardButtonEl = root.querySelector(".oviz-three-play-backward");
      const playForwardButtonEl = root.querySelector(".oviz-three-play-forward");
      const footerEl = root.querySelector(".oviz-three-footer");
      const earthViewToggleButtonEl = root.querySelector(".oviz-three-earth-view-toggle");
      const tooltipEl = root.querySelector(".oviz-three-tooltip");
      const scaleBarEl = root.querySelector(".oviz-three-scale-bar");
      const scaleLabelEl = root.querySelector(".oviz-three-scale-label");
      const noteEl = root.querySelector(".oviz-three-note");
      const clearSelectionButtonEl = root.querySelector(".oviz-three-selection-clear");
      const clickSelectToggleEl = root.querySelector(".oviz-three-click-select-toggle");
      const volumeLassoToggleEl = root.querySelector(".oviz-three-volume-lasso-toggle");
      const themeSelectEl = root.querySelector(".oviz-three-theme-select");
      const scrollSpeedEl = root.querySelector(".oviz-three-scroll-speed");
      const scrollSpeedLabelEl = root.querySelector(".oviz-three-scroll-speed-label");
      const cameraFovEl = root.querySelector(".oviz-three-camera-fov");
      const cameraFovLabelEl = root.querySelector(".oviz-three-camera-fov-label");
      const globalPointSizeEl = root.querySelector(".oviz-three-global-point-size");
      const globalPointSizeLabelEl = root.querySelector(".oviz-three-global-point-size-label");
      const globalPointOpacityEl = root.querySelector(".oviz-three-global-point-opacity");
      const globalPointOpacityLabelEl = root.querySelector(".oviz-three-global-point-opacity-label");
      const globalPointGlowEl = root.querySelector(".oviz-three-global-point-glow");
      const globalPointGlowLabelEl = root.querySelector(".oviz-three-global-point-glow-label");
      const skyDomeSourceFieldEl = root.querySelector(".oviz-three-sky-source-field");
      const skyDomeSourceSelectEl = root.querySelector(".oviz-three-sky-source-select");
      const skyDomeControlsEl = root.querySelector(".oviz-three-sky-dome-controls");
      const skyDomeVisibleToggleEl = root.querySelector(".oviz-three-sky-dome-visible-toggle");
      const skyDomeForceVisibleToggleEl = root.querySelector(".oviz-three-sky-dome-force-visible-toggle");
      const skyDomeOpacityEl = root.querySelector(".oviz-three-sky-dome-opacity");
      const skyDomeOpacityLabelEl = root.querySelector(".oviz-three-sky-dome-opacity-label");
      const skyDomeBrightnessEl = root.querySelector(".oviz-three-sky-dome-brightness");
      const skyDomeBrightnessLabelEl = root.querySelector(".oviz-three-sky-dome-brightness-label");
      const skyDomeContrastEl = root.querySelector(".oviz-three-sky-dome-contrast");
      const skyDomeContrastLabelEl = root.querySelector(".oviz-three-sky-dome-contrast-label");
      const skyDomeGammaEl = root.querySelector(".oviz-three-sky-dome-gamma");
      const skyDomeGammaLabelEl = root.querySelector(".oviz-three-sky-dome-gamma-label");
      const skyDomeStatusEl = root.querySelector(".oviz-three-sky-dome-status");
      const skyDomeHipsControlEls = Array.from(root.querySelectorAll(".oviz-three-sky-dome-hips-controls"));
      const skyLayerPresetSelectEl = root.querySelector(".oviz-three-sky-layer-preset-select");
      const skyLayerCustomInputEl = root.querySelector(".oviz-three-sky-layer-custom-input");
      const skyLayerSearchDatalistEl = root.querySelector(".oviz-three-sky-hips-search");
      const skyLayerAddButtonEl = root.querySelector(".oviz-three-sky-layer-add");
      const skyLayerListEl = root.querySelector(".oviz-three-sky-layer-list");
      const sizeByStarsToggleEl = root.querySelector(".oviz-three-size-by-stars-toggle");
      const focusGroupSelectEl = root.querySelector(".oviz-three-focus-group-select");
      const fadeTimeEl = root.querySelector(".oviz-three-fade-time");
      const fadeInOutToggleEl = root.querySelector(".oviz-three-fade-in-out-toggle");
      const axesVisibleToggleEl = root.querySelector(".oviz-three-axes-visible-toggle");
      const galacticReferenceToggleEl = root.querySelector(".oviz-three-galactic-reference-toggle");
      const nearbyRegionLabelsToggleEls = Array.from(root.querySelectorAll(".oviz-three-region-labels-toggle"));
      const manualLabelSelectEl = root.querySelector(".oviz-three-manual-label-select");
      const manualLabelTextEl = root.querySelector(".oviz-three-manual-label-text");
      const manualLabelSizeEl = root.querySelector(".oviz-three-manual-label-size");
      const manualLabelAddButtonEl = root.querySelector(".oviz-three-manual-label-add");
      const manualLabelApplyButtonEl = root.querySelector(".oviz-three-manual-label-apply");
      const manualLabelDeleteButtonEl = root.querySelector(".oviz-three-manual-label-delete");
      const manualLabelReadoutEl = root.querySelector(".oviz-three-manual-label-readout");
      const orbitCameraButtons = Array.from(root.querySelectorAll(".oviz-three-auto-orbit"));
      const viewFromEarthButtonEl = root.querySelector(".oviz-three-view-from-earth");
      const resetCameraButtonEl = root.querySelector(".oviz-three-reset-camera");
      const resetControlsButtonEl = root.querySelector(".oviz-three-reset-controls");
      const volumePanelEl = root.querySelector(".oviz-three-volume");
      const volumeSelectEl = root.querySelector(".oviz-three-volume-select");
      const volumeSmoothingFieldEl = root.querySelector(".oviz-three-volume-smoothing-field");
      const volumeSmoothingEl = root.querySelector(".oviz-three-volume-smoothing");
      const volumeVisibleEl = root.querySelector(".oviz-three-volume-visible");
      const volumeColormapEl = root.querySelector(".oviz-three-volume-colormap");
      const volumeStretchEl = root.querySelector(".oviz-three-volume-stretch");
      const volumeVMinEl = root.querySelector(".oviz-three-volume-vmin");
      const volumeVMaxEl = root.querySelector(".oviz-three-volume-vmax");
      const volumeOpacityEl = root.querySelector(".oviz-three-volume-opacity");
      const volumeAlphaEl = root.querySelector(".oviz-three-volume-alpha");
      const volumeStepsEl = root.querySelector(".oviz-three-volume-steps");
      const volumeOpacityLabelEl = root.querySelector(".oviz-three-volume-opacity-label");
      const volumeAlphaLabelEl = root.querySelector(".oviz-three-volume-alpha-label");
      const volumeStepsLabelEl = root.querySelector(".oviz-three-volume-steps-label");
      const volumeSummaryEl = root.querySelector(".oviz-three-volume-summary");
      const lassoOverlayEl = root.querySelector(".oviz-three-lasso-overlay");
      const lassoPolylineEl = root.querySelector(".oviz-three-lasso-overlay polyline");
      const skyPanelEl = root.querySelector(".oviz-three-sky-panel");
      const skyFrameEl = root.querySelector(".oviz-three-sky-frame");
      const skyDomeFrameEl = root.querySelector(".oviz-three-sky-dome-frame");
      const skyBodyEl = null;
      const skyCanvasEl = root.querySelector(".oviz-three-sky-canvas");
      const skyImageSelectEl = root.querySelector(".oviz-three-sky-image-select");
      const skyReadoutEl = root.querySelector(".oviz-three-sky-readout");
      const skyStatusEl = root.querySelector(".oviz-three-sky-status");
      const skyFullButtonEl = root.querySelector(".oviz-three-sky-full");
      const skyHideButtonEl = root.querySelector(".oviz-three-sky-hide");
      const ageKdePanelEl = root.querySelector(".oviz-three-age-panel");
      const ageKdeCanvasEl = root.querySelector(".oviz-three-age-canvas");
      const ageKdeFilterRangeMinEl = root.querySelector(".oviz-three-age-filter-range-min");
      const ageKdeFilterRangeMaxEl = root.querySelector(".oviz-three-age-filter-range-max");
      const ageKdeFilterRangeReadoutMinEl = root.querySelector(".oviz-three-age-filter-range-readout-min");
      const ageKdeFilterRangeReadoutMaxEl = root.querySelector(".oviz-three-age-filter-range-readout-max");
      const ageKdeFullButtonEl = root.querySelector(".oviz-three-age-full");
      const ageKdeHideButtonEl = root.querySelector(".oviz-three-age-hide");
      const clusterFilterPanelEl = root.querySelector(".oviz-three-filter-panel");
      const clusterFilterCanvasEl = root.querySelector(".oviz-three-filter-canvas");
      const clusterFilterParameterEl = root.querySelector(".oviz-three-filter-parameter");
      const clusterFilterRangeMinEl = root.querySelector(".oviz-three-filter-range-min");
      const clusterFilterRangeMaxEl = root.querySelector(".oviz-three-filter-range-max");
      const clusterFilterRangeReadoutMinEl = root.querySelector(".oviz-three-filter-range-readout-min");
      const clusterFilterRangeReadoutMaxEl = root.querySelector(".oviz-three-filter-range-readout-max");
      const clusterFilterFullButtonEl = root.querySelector(".oviz-three-filter-full");
      const clusterFilterHideButtonEl = root.querySelector(".oviz-three-filter-hide");
      const boxMetricsPanelEl = root.querySelector(".oviz-three-box-panel");
      const boxMetricsCanvasEl = root.querySelector(".oviz-three-box-canvas");
      const boxMetricsSummaryEl = root.querySelector(".oviz-three-box-summary");
      const boxMetricsVisibleEl = root.querySelector(".oviz-three-box-visible");
      const boxMetricsFullButtonEl = root.querySelector(".oviz-three-box-full");
      const boxMetricsHideButtonEl = root.querySelector(".oviz-three-box-hide");
      const boxMetricsResetButtonEl = root.querySelector(".oviz-three-box-reset");
      const dendrogramPanelEl = root.querySelector(".oviz-three-dendrogram-panel");
      const dendrogramCanvasEl = root.querySelector(".oviz-three-dendrogram-canvas");
      const dendrogramTraceEl = root.querySelector(".oviz-three-dendrogram-trace");
      const dendrogramConnectionEl = root.querySelector(".oviz-three-dendrogram-connection");
      const dendrogramModeEl = root.querySelector(".oviz-three-dendrogram-mode");
      const dendrogramThresholdLabelEl = root.querySelector(".oviz-three-dendrogram-threshold-label");
      const dendrogramThresholdEl = root.querySelector(".oviz-three-dendrogram-threshold");
      const dendrogramFullButtonEl = root.querySelector(".oviz-three-dendrogram-full");
      const dendrogramHideButtonEl = root.querySelector(".oviz-three-dendrogram-hide");
      const widgetPanels = Array.from(root.querySelectorAll(".oviz-three-widget-panel"));
      const widgetDragHandles = Array.from(root.querySelectorAll(".oviz-three-widget-drag"));
      const widgetResizeEls = Array.from(root.querySelectorAll(".oviz-three-widget-resize"));

      const baseTheme = safeJsonClone(sceneSpec.theme || {}, {});
      const theme = safeJsonClone(baseTheme, {});
      if (titleEl) {
        titleEl.textContent = sceneSpec.title || "";
      }

      function buildThemePresets(sourceTheme) {
        const defaultTheme = safeJsonClone(sourceTheme || {}, {});
        return {
          default: Object.assign({}, defaultTheme),
          dark: Object.assign({}, defaultTheme, {
            paper_bgcolor: "#000000",
            scene_bgcolor: "#000000",
            text_color: "#d0d0d0",
            axis_color: "#808080",
            panel_bg: "rgba(0, 0, 0, 0.50)",
            panel_border: "rgba(128, 128, 128, 0.55)",
            panel_solid: "#121212",
            footprint: defaultTheme.footprint || "#6ec5ff",
          }),
          light: Object.assign({}, defaultTheme, {
            paper_bgcolor: "#f4f1eb",
            scene_bgcolor: "#fbfaf7",
            text_color: "#1a1a1a",
            axis_color: "#666666",
            panel_bg: "rgba(255, 255, 255, 0.84)",
            panel_border: "rgba(0, 0, 0, 0.18)",
            panel_solid: "#ffffff",
            footprint: "#1d70b8",
          }),
          midnight: Object.assign({}, defaultTheme, {
            paper_bgcolor: "#050913",
            scene_bgcolor: "#050913",
            text_color: "#e7eefc",
            axis_color: "#7f90b8",
            panel_bg: "rgba(8, 13, 32, 0.56)",
            panel_border: "rgba(126, 150, 214, 0.26)",
            panel_solid: "#0a1126",
            footprint: "#8ecbff",
          }),
          glacier: Object.assign({}, defaultTheme, {
            paper_bgcolor: "#ecf5f9",
            scene_bgcolor: "#f8fcfd",
            text_color: "#173342",
            axis_color: "#59727e",
            panel_bg: "rgba(255, 255, 255, 0.80)",
            panel_border: "rgba(27, 66, 82, 0.16)",
            panel_solid: "#ffffff",
            footprint: "#0b91c1",
          }),
          ember: Object.assign({}, defaultTheme, {
            paper_bgcolor: "#140906",
            scene_bgcolor: "#190b07",
            text_color: "#f7e7de",
            axis_color: "#d39a75",
            panel_bg: "rgba(41, 18, 12, 0.54)",
            panel_border: "rgba(233, 165, 118, 0.22)",
            panel_solid: "#1d0d0a",
            footprint: "#ffb273",
          }),
          graphite: Object.assign({}, defaultTheme, {
            paper_bgcolor: "#0d1016",
            scene_bgcolor: "#0d1016",
            text_color: "#eef1f4",
            axis_color: "#8d96a4",
            panel_bg: "rgba(22, 27, 36, 0.54)",
            panel_border: "rgba(194, 204, 220, 0.14)",
            panel_solid: "#161b24",
            footprint: "#9fc5ff",
          }),
          aurora: Object.assign({}, defaultTheme, {
            paper_bgcolor: "#071413",
            scene_bgcolor: "#091917",
            text_color: "#e6faf5",
            axis_color: "#7dbdb2",
            panel_bg: "rgba(10, 29, 27, 0.56)",
            panel_border: "rgba(121, 214, 196, 0.18)",
            panel_solid: "#0d2421",
            footprint: "#9df6cf",
          }),
          paper: Object.assign({}, defaultTheme, {
            paper_bgcolor: "#f3ede2",
            scene_bgcolor: "#fcfaf4",
            text_color: "#2e2821",
            axis_color: "#7a6d60",
            panel_bg: "rgba(255, 252, 247, 0.80)",
            panel_border: "rgba(92, 74, 54, 0.16)",
            panel_solid: "#fffdf8",
            footprint: "#3f7f97",
          }),
        };
      }

      function applyThemeCssVars() {
        root.style.setProperty("--oviz-paper-bg", theme.paper_bgcolor || "#000000");
        root.style.setProperty("--oviz-scene-bg", theme.scene_bgcolor || theme.paper_bgcolor || "#000000");
        root.style.setProperty("--oviz-text", theme.text_color || "#d0d0d0");
        root.style.setProperty("--oviz-axis", theme.axis_color || "#808080");
        root.style.setProperty("--oviz-muted-text", theme.muted_text || "rgba(238, 242, 247, 0.58)");
        root.style.setProperty("--oviz-accent", theme.accent || "#d8dde6");
        root.style.setProperty("--oviz-accent-strong", theme.accent_strong || "#eef1f5");
        root.style.setProperty("--oviz-panel-bg", theme.panel_bg || "rgba(0, 0, 0, 0.45)");
        root.style.setProperty("--oviz-panel-border", theme.panel_border || "rgba(128, 128, 128, 0.50)");
        root.style.setProperty("--oviz-panel-solid", theme.panel_solid || theme.paper_bgcolor || "#121212");
        root.style.setProperty("--oviz-shadow-md", theme.shadow_md || "0 18px 42px rgba(0, 0, 0, 0.24)");
        root.style.setProperty("--oviz-shadow-lg", theme.shadow_lg || "0 24px 70px rgba(0, 0, 0, 0.34)");
        root.style.setProperty("--oviz-footprint", theme.footprint || "#6ec5ff");
      }

      applyThemeCssVars();
      if (titleEl) {
        titleEl.style.display = sceneSpec.title ? "block" : "none";
      }

      if (sceneSpec.note) {
        noteEl.textContent = sceneSpec.note;
        noteEl.style.display = "block";
      }

      const layoutScene = (sceneSpec.layout || {}).scene || {};
      const axisSpec = sceneSpec.axes || {};
      const frameSpecs = sceneSpec.frames || [];
      const timelineSpec = sceneSpec.timeline || { enabled: frameSpecs.length > 1 };
      const timelineEnabled = timelineSpec.enabled !== false;
      const legendItems = (sceneSpec.legend || {}).items || [];
      const groupVisibility = sceneSpec.group_visibility || {};
      const defaultGroup = sceneSpec.default_group || "All";
      const skySpec = sceneSpec.sky_panel || { enabled: false };
      const skyDomeSpec = sceneSpec.sky_dome || { enabled: false };
      const selectionBoxSpec = sceneSpec.selection_box || { enabled: false };
      const ageKdeSpec = sceneSpec.age_kde || { enabled: false };
      const clusterFilterSpec = sceneSpec.cluster_filter || { enabled: false };
      const dendrogramSpec = sceneSpec.dendrogram || { enabled: false };
      const volumeSpec = sceneSpec.volumes || { enabled: false, layers: [] };
      const volumeLayers = (volumeSpec.layers || []);
      const volumeCoRotationRateRadPerMyr = Number(volumeSpec.co_rotation_rate_rad_per_myr);
      const volumeLayersByKey = new Map(volumeLayers.map((layer) => [String(layer.key), layer]));
      const imagePlaneSpecs = Array.isArray(sceneSpec.image_planes) ? sceneSpec.image_planes : [];
      const imagePlaneSpecByKey = new Map(
        imagePlaneSpecs
          .map((spec) => {
            const key = String((spec && spec.key) || "");
            return key ? [key, spec] : null;
          })
          .filter(Boolean)
      );
      const volumeStateKeyForRawLayer = (layer) => String((layer && (layer.state_key || layer.key)) || "");
      const volumeStateKeys = Array.from(
        new Set(
          volumeLayers
            .map((layer) => volumeStateKeyForRawLayer(layer))
            .filter((key) => key)
        )
      );
      if (footerEl && !timelineEnabled) {
        footerEl.style.display = "none";
      }
      const volumeStateKeySet = new Set(volumeStateKeys);
      const animationSpec = sceneSpec.animation || {};
      const clusterFilterParameters = Array.isArray(clusterFilterSpec.parameters) ? clusterFilterSpec.parameters : [];
      const clusterFilterEntries = Array.isArray(clusterFilterSpec.entries) ? clusterFilterSpec.entries : [];
      const clusterFilterEntryByKey = new Map(
        clusterFilterEntries
          .map((entry) => {
            const key = normalizeMemberKey(entry && entry.selection_key ? entry.selection_key : "");
            return key ? [key, entry] : null;
          })
          .filter(Boolean)
      );
      const selectionBoxMetricsSpec = selectionBoxSpec.metrics || { enabled: false };
      const selectionBoxTransformSpec = selectionBoxSpec.transform || {};
      const selectionBoxReferenceOrbitSpec = selectionBoxTransformSpec.reference_orbit || {};
      const selectionBoxDefaultBandSpec = (() => {
        const explicitDefault = Array.isArray(selectionBoxMetricsSpec.twopcf_default_band_pc)
          ? selectionBoxMetricsSpec.twopcf_default_band_pc
          : null;
        const fallbackBand = Array.isArray(selectionBoxMetricsSpec.twopcf_bands_pc)
          && selectionBoxMetricsSpec.twopcf_bands_pc.length
          ? selectionBoxMetricsSpec.twopcf_bands_pc[0]
          : [30.0, 80.0];
        const sourceBand = explicitDefault || fallbackBand;
        const minValue = Number(sourceBand && sourceBand[0]);
        const maxValue = Number(sourceBand && sourceBand[1]);
        if (Number.isFinite(minValue) && Number.isFinite(maxValue) && maxValue > minValue) {
          return [minValue, maxValue];
        }
        return [30.0, 80.0];
      })();
      const galacticReferenceTraceNames = new Set([
        "GC",
        "GC Ring",
        "R = 4 kpc",
        "R = 8.12 kpc",
        "R = 12 kpc",
        "Galactic Quadrants",
        "Galactic l Ticks",
        "Galactic l Labels",
        "Galactic Z Axis",
      ]);
      const nearbyRegionLabelTraceNames = new Set([
        "Nearby Region Labels",
      ]);
      const selectionBoxMetricTimeCenters = buildSelectionBoxMetricTimeCenters();
      const selectionBoxMetricEvents = normalizeSelectionBoxMetricEvents(selectionBoxMetricsSpec.events || {});
      const selectionBoxRandomCatalogSize = Math.max(
        32,
        Math.round(Number(selectionBoxMetricsSpec.random_catalog_size) || 256)
      );
      const selectionBoxRandomUnitCatalog = buildSelectionBoxRandomUnitCatalog(
        selectionBoxRandomCatalogSize,
        734921
      );
      const themePresets = buildThemePresets(baseTheme);
      const pointScaleBaseline = Math.max(Number(sceneSpec.point_size_baseline_scale) || (4.0 / 3.0), 0.01);
      const pointScale = (Math.max(sceneSpec.max_span || 1, 1) / 2600.0) * pointScaleBaseline;
      const hoverTargets = [];
      const cameraResponsivePointEntries = [];
      const cameraResponsiveImagePlaneEntries = [];
      const selectionSpriteEntriesByKey = new Map();
      const axisLineMaterials = [];
      const frameLineMaterials = [];
      const screenStableTextSprites = [];
      const markerTextureCache = new Map();
      const markerMaterialCache = new Map();
      const starCoreTextureCache = new Map();
      const starCoreMaterialCache = new Map();
      const starGlowTextureCache = new Map();
      const starGlowMaterialCache = new Map();
      const textTextureCache = new Map();
      const galaxyTextureCache = new Map();
      const imagePlaneTextureCache = new Map();
      const volumeScalarDataCache = new Map();
      const volumeScalarDataPendingCache = new Map();
      const volumeTextureCache = new Map();
      const volumeColorTextureCache = new Map();
      const volumeColorBytesCache = new Map();
      const volumeRuntimeByKey = new Map();
      const DEFAULT_MANUAL_LABEL_SIZE = 24.0;
      const LEGEND_MAX_VISIBLE_ITEMS = 5;
      const MIN_MANUAL_LABEL_SIZE = 8.0;
      const MAX_MANUAL_LABEL_SIZE = 120.0;
      const MAX_MANUAL_LABEL_TEXT_LENGTH = 80;
      let volumeJitterTexture = null;
      let legendState = {};
      let currentGroup = defaultGroup;
      let currentFrameIndex = sceneSpec.initial_frame_index || 0;
      let currentSelection = null;
      let currentSelections = [];
      let currentSelectionMode = "none";
      let clickSelectionEnabled = false;
      let lassoVolumeSelectionEnabled = true;
      let playbackDirection = 0;
      let displayedFrameValue = currentFrameIndex;
      let frameTransitionState = null;
      let lastPlaybackAdvanceTimestamp = null;
      let lastTransitionRenderTimestamp = null;
      let pendingSliderFrameIndex = null;
      let sliderScrubRenderHandle = null;
      let skyPanelMode = "hidden";
      let boxMetricsPanelMode = "hidden";
      let ageKdePanelMode = "hidden";
      let clusterFilterPanelMode = "hidden";
      let dendrogramPanelMode = "hidden";
      let widgetPointerState = null;
      let legendPointerState = null;
      let scaleBarPointerState = null;
      let widgetZIndexCounter = 8;
      let selectionBoxPointerState = null;
      let lassoState = null;
      let lassoArmed = false;
      let suppressNextCanvasClick = false;
      let selectedClusterKeys = new Set();
      let localHoveredClusterKey = "";
      let skyHoveredClusterKey = "";
      let lastSentSkyHoverClusterKey = null;
      let scaleBarPosition = null;
      let currentScaleBarLengthPc = NaN;
      let skyDomeMesh = null;
      let skyDomeTexture = null;
      let skyDomeTexturePriority = -1;
      let skyDomeHipsState = null;
      let skyDomeHipsWarmupImages = [];
      let skyDomeSurvey = "";
      let skyDomeProjection = String(skyDomeSpec.projection || "MOL");
      let skyDomeSnapshotStatus = "idle";
      const skyDomeDefaultEnabled = Boolean(skyDomeSpec.enabled);
      const skyDomeDefaultForceVisible = Boolean(skyDomeSpec.force_visible || skyDomeSpec.always_visible);
      const skyDomeDefaultOpacity = Math.min(Math.max(Number(skyDomeSpec.opacity ?? 0.55), 0.0), 1.0);
      const skyDomeDefaultHipsBrightness = Math.min(Math.max(Number(skyDomeSpec.hips_brightness) || 2.4, 0.1), 8.0);
      const skyDomeDefaultHipsContrast = Math.min(Math.max(Number(skyDomeSpec.hips_contrast) || 1.25, 0.1), 4.0);
      const skyDomeDefaultHipsGamma = Math.min(Math.max(Number(skyDomeSpec.hips_gamma) || 1.35, 0.2), 4.0);
      let skyDomeForceVisible = skyDomeDefaultForceVisible;
      let skyDomeCaptureFrameSignature = "";
      let skyDomeBackgroundFrameReady = false;
      let skyDomeBackgroundViewSignature = "";
      let skyDomeBackgroundLatestViewSignature = "";
      let skyDomeBackgroundLastSentAt = 0.0;
      let skyDomeBackgroundMotionHoldUntil = 0.0;
      let skyDomeBackgroundUserCameraActive = false;
      let skyDomeBackgroundViewSequence = 0;
      let skyDomeBackgroundAlignedView = null;
      const skyDomeBackgroundSentViews = new Map();
      const skyDomeSourceOptions = normalizeSkyDomeSourceOptions(skyDomeSpec);
      const skyDomeSourceOptionByKey = new Map(skyDomeSourceOptions.map((source) => [source.key, source]));
      let activeSkyDomeSourceKey = "";
      let skyDomeLocalSourceActive = false;
      const fallbackSkyLayerPresetOptions = [
        { key: "P/Mellinger/color", label: "Mellinger Color", survey: "P/Mellinger/color", color: "#f2c178" },
        { key: "P/DSS2/color", label: "DSS2 Color", survey: "P/DSS2/color", color: "#8fbfff" },
        { key: "P/PLANCK/R2/HFI/color", label: "Planck Dust Emission Color", survey: "P/PLANCK/R2/HFI/color", color: "#ffc56e" },
        { key: "P/PanSTARRS/DR1/color-z-zg-g", label: "Pan-STARRS DR1 Color", survey: "P/PanSTARRS/DR1/color-z-zg-g", color: "#f2c178" },
        { key: "P/2MASS/color", label: "2MASS Color", survey: "P/2MASS/color", color: "#ff9a7a" },
        { key: "P/allWISE/color", label: "AllWISE Color", survey: "P/allWISE/color", color: "#d6a4ff" },
        { key: "P/GALEXGR6/AIS/color", label: "GALEX AIS Color", survey: "P/GALEXGR6/AIS/color", color: "#94f4ff" },
        { key: "P/SDSS9/color", label: "SDSS9 Color", survey: "P/SDSS9/color", color: "#a7e0ff" },
        { key: "P/DECaLS/DR5/color", label: "DECaLS DR5 Color", survey: "P/DECaLS/DR5/color", color: "#c7ef9f" },
      ];
      const skyLayerStretchOptions = [
        { value: "linear", label: "Linear" },
        { value: "log", label: "Log10" },
        { value: "asinh", label: "Asinh" },
      ];
      const skyLayerColormapOptions = [
        { value: "native", label: "Native" },
        { value: "grayscale", label: "Gray" },
        { value: "viridis", label: "Viridis" },
        { value: "plasma", label: "Plasma" },
        { value: "magma", label: "Magma" },
        { value: "inferno", label: "Inferno" },
        { value: "cividis", label: "Cividis" },
        { value: "rainbow", label: "Rainbow" },
        { value: "cubehelix", label: "Cubehelix" },
      ];
      let skyLayerPresetOptions = fallbackSkyLayerPresetOptions.slice();
      let skyLayerPresetBySurvey = buildSkyLayerPresetBySurvey(skyLayerPresetOptions);
      let skyLayerHipsRegistry = [];
      let skyLayerHipsRegistryBySurvey = new Map();
      let skyLayerHipsRegistryPromise = null;
      let skyLayerHipsRegistryLoaded = false;
      let skyLayerHipsRegistryFetchAttempted = false;
      let skyLayerHipsRegistryError = "";
      let skyLayerSearchRenderSignature = "";
      let skyLayerSearchUpdateTimer = 0;
      let skyLayerState = [];
      let skyLayerStateInitialized = false;
      let activeSkyLayerKey = "";
      const skyPanelImageCache = new Map();
      let skyPanelRenderSerial = 0;
      let currentZoomAnchorPoint = null;
      let galacticSimpleOrbitTargetTrackingActive = false;
      let legendPanelRectState = null;
      let legendPanelUserSized = false;
      let legendSectionOpenState = { traces: true, volumes: true };
      let activeThemeKey = "default";
      let axesVisible = Boolean(sceneSpec.show_axes);
      let galacticReferenceVisible = true;
      let nearbyRegionLabelsVisible = true;
      let manualLabels = [];
      let activeManualLabelId = "";
      let manualLabelDraftText = "";
      let manualLabelDraftSize = DEFAULT_MANUAL_LABEL_SIZE;
      let manualLabelPointerState = null;
      let manualLabelIdCounter = 0;
      let globalPointSizeScale = 1.0;
      let globalPointOpacityScale = 1.0;
      let globalPointGlowStrength = 0.60;
      let sizePointsByStarsEnabled = false;
      let globalScrollSpeed = 1.0;
      let cameraAutoOrbitBaseSpeed = 1.0;
      let cameraAutoOrbitDirection = 1.0;
      let cameraAutoOrbitSpeedMultiplier = 1.0;
      let cameraAutoOrbitEnabled = Boolean(
        initialState.global_controls && initialState.global_controls.camera_auto_orbit_enabled
      );
      let zenModeEnabled = Boolean(initialState.zen_mode_enabled);
      let legendPanelOpen = initialState.legend_open === undefined ? true : Boolean(initialState.legend_open);
      let fadeInTimeMyr = Number(animationSpec.fade_in_time_default);
      let fadeInAndOutEnabled = Boolean(animationSpec.fade_in_and_out_default);
      let focusTraceKey = String(animationSpec.focus_trace_key_default || "");
      let focusSelectionKey = "";
      let cameraViewMode = "free";
      let earthViewFocusDistance = null;
      let earthViewReturnCameraState = null;
      let cameraTransitionAnimationFrame = 0;
      let skyDomeOpacityAnimationFrame = 0;
      let skyDomeViewOpacityScale = 0.0;
      let skyViewTransitionSerial = 0;
      let skyViewDragState = null;
      let currentLassoSelectionMask = null;
      const pressedKeys = new Set();
      let lastAnimationTimestamp = null;
      const playbackIntervalMs = Math.max(Number(sceneSpec.playback_interval_ms) || 240, 80);
      let clusterFilterParameterKey = String(clusterFilterSpec.default_parameter_key || "");
      const clusterFilterRangeStateByKey = {};
      let dendrogramTraceKey = String(dendrogramSpec.default_trace_key || "");
      let dendrogramConnectionMode = String(dendrogramSpec.default_connection_mode || "birth_to_older_track");
      let dendrogramThresholdMode = String(dendrogramSpec.default_threshold_mode || "distance_pc");
      let dendrogramThresholdPc = Number(dendrogramSpec.default_threshold_pc || 100.0);
      let dendrogramThresholdAgeMyr = Number(dendrogramSpec.default_threshold_age_myr || 5.0);
      let dendrogramHoveredSelectionKeys = new Set();
      let dendrogramPinnedSelectionKeys = new Set();
      let dendrogramHoveredBranchLabel = "";
      let dendrogramHoveredBranchCount = 0;
      let dendrogramPinnedBranchLabel = "";
      let dendrogramPinnedBranchCount = 0;
      let dendrogramHoveredRegionKey = "";
      let dendrogramPinnedRegionKey = "";
      let dendrogramHitRegions = [];
      let activeVolumeKey = volumeStateKeys.length ? String(volumeStateKeys[0]) : null;
      const volumeStateByKey = {};
      const traceStyleStateByKey = {};
      const legendEditButtonByKey = new Map();
      let activeLegendEditorKey = "";
      let selectionBoxState = buildDefaultSelectionBoxState();
      let selectionBoxHitObjects = [];
      let selectionBoxMetricsCache = null;
      let selectionBoxMetricsTimer = null;
      let selectionBoxMetricsVersion = 0;
      let selectionBoxMetricsPending = false;

      const renderer = new THREE.WebGLRenderer({
        canvas,
        antialias: true,
        alpha: true,
        powerPreference: "high-performance",
      });
      renderer.setPixelRatio(window.devicePixelRatio || 1);
      if ("outputColorSpace" in renderer && THREE.SRGBColorSpace) {
        renderer.outputColorSpace = THREE.SRGBColorSpace;
      } else if ("outputEncoding" in renderer && THREE.sRGBEncoding) {
        renderer.outputEncoding = THREE.sRGBEncoding;
      }
      const volumeSupported = volumeStateKeys.length > 0 && Boolean(renderer.capabilities && renderer.capabilities.isWebGL2);

      const scene = new THREE.Scene();
      function skyDomeUsesAladinBackground() {
        const mode = String(
          (skyDomeSpec && (
            skyDomeSpec.background_mode
            || skyDomeSpec.mode
            || skyDomeSpec.render_mode
          ))
          || ""
        ).toLowerCase();
        return Boolean(
          skyDomeSpec
          && skyDomeSpec.enabled
          && String(skyDomeSpec.source || "").toLowerCase() === "aladin"
          && (mode === "live_aladin" || mode === "aladin-background")
        );
      }

      function skyDomeUsesNativeHips() {
        const mode = String(
          (skyDomeSpec && (
            skyDomeSpec.background_mode
            || skyDomeSpec.mode
            || skyDomeSpec.render_mode
          ))
          || ""
        ).toLowerCase();
        const source = String(skyDomeSpec && skyDomeSpec.source || "").toLowerCase();
        return Boolean(
          skyDomeSpec
          && skyDomeSpec.enabled
          && (
            source === "hips"
            || source === "native_hips"
            || source === "native-hips"
            || mode === "native_hips"
            || mode === "native-hips"
            || mode === "hips"
          )
        );
      }

      function skyDomeUsesHips2Fits() {
        const mode = String(
          (skyDomeSpec && (
            skyDomeSpec.background_mode
            || skyDomeSpec.mode
            || skyDomeSpec.render_mode
          ))
          || ""
        ).toLowerCase();
        const source = String(skyDomeSpec && skyDomeSpec.source || "").toLowerCase();
        return Boolean(
          skyDomeSpec
          && skyDomeSpec.enabled
          && (
            source === "hips2fits"
            || source === "hips-2-fits"
            || mode === "hips2fits"
            || mode === "hips-2-fits"
          )
        );
      }

      function applySceneBackground() {
        const bgColor = new THREE.Color(theme.scene_bgcolor || theme.paper_bgcolor || "#000000");
        if (skyDomeUsesAladinBackground()) {
          scene.background = null;
          renderer.setClearColor(bgColor, 0.0);
        } else {
          scene.background = bgColor;
          renderer.setClearColor(bgColor, 1.0);
        }
      }

      const volumeSkyAxisTransform = deriveVolumeSkyAxisTransform();
      if (root && root.dataset) {
        root.dataset.skyAxisTransform = JSON.stringify(volumeSkyAxisTransform);
      }

      applySceneBackground();

      const sceneUp = sceneSpec.camera_up || { x: 0.0, y: 0.0, z: 1.0 };
      const sceneUpVector = new THREE.Vector3(sceneUp.x ?? 0.0, sceneUp.y ?? 0.0, sceneUp.z ?? 1.0).normalize();
      const CAMERA_AUTO_ORBIT_SPEED = 1.2;
      function cameraFarPlanePc() {
        const farCandidates = [Math.max((sceneSpec.max_span || 1) * 20.0, 10000.0)];
        if (skyDomeSpec && skyDomeSpec.enabled !== false) {
          farCandidates.push((Number(skyDomeSpec.radius_pc) || 0.0) * 2.5);
        }
        imagePlaneSpecs.forEach((imagePlaneSpec) => {
          farCandidates.push(Number(imagePlaneSpec && imagePlaneSpec.width_pc) || 0.0);
          farCandidates.push(Number(imagePlaneSpec && imagePlaneSpec.height_pc) || 0.0);
        });
        return Math.max(...farCandidates.filter((value) => Number.isFinite(value) && value > 0.0));
      }
      const camera = new THREE.PerspectiveCamera(60, 1, 0.1, cameraFarPlanePc());
      camera.up.set(sceneUp.x ?? 0.0, sceneUp.y ?? 0.0, sceneUp.z ?? 1.0);
      const controls = new OrbitControls(camera, renderer.domElement);
      const liveAladinSkyBackground = skyDomeUsesAladinBackground();
      controls.enableDamping = !liveAladinSkyBackground;
      controls.dampingFactor = liveAladinSkyBackground ? 0.0 : 0.08;
      controls.rotateSpeed = 0.7;
      controls.panSpeed = 0.7;
      controls.zoomSpeed = globalScrollSpeed;
      controls.autoRotate = cameraAutoOrbitEnabled;
      controls.autoRotateSpeed = CAMERA_AUTO_ORBIT_SPEED;
      controls.minPolarAngle = 0.02;
      controls.maxPolarAngle = Math.PI - 0.02;
      controls.target.set(sceneSpec.center.x, sceneSpec.center.y, sceneSpec.center.z);
      controls.addEventListener("start", () => {
        handleManualCameraInteractionStart();
        setSkyDomeBackgroundCameraActive(true);
      });
      controls.addEventListener("end", () => {
        setSkyDomeBackgroundCameraActive(false);
      });
      controls.addEventListener("change", () => {
        if (!skyDomeUsesAladinBackground()) {
          return;
        }
        updateSkyDomeBackgroundFrame(
          (typeof performance !== "undefined" && performance.now) ? performance.now() : Date.now(),
          { force: skyDomeBackgroundUserCameraActive }
        );
      });

      const eye = (layoutScene.camera || {}).eye || { x: 0.0, y: -0.8, z: 1.5 };
      const eyeScale = Math.max(sceneSpec.max_span || 1, 1);
      camera.position.set(
        sceneSpec.center.x + (eye.x ?? 0.0) * eyeScale,
        sceneSpec.center.y + (eye.y ?? -0.8) * eyeScale,
        sceneSpec.center.z + (eye.z ?? 1.5) * eyeScale
      );
      camera.lookAt(controls.target);
      const initialCameraState = {
        position: camera.position.clone(),
        target: controls.target.clone(),
        up: camera.up.clone(),
        fov: Number(camera.fov),
        viewOffset: { x: 0.0, y: 0.0 },
      };

      const plotGroup = new THREE.Group();
      scene.add(plotGroup);

      const axisGroup = new THREE.Group();
      scene.add(axisGroup);

      volumeLayers.forEach((layer) => {
        const stateKey = volumeStateKeyForRawLayer(layer);
        if (!stateKey || Object.prototype.hasOwnProperty.call(volumeStateByKey, stateKey)) {
          return;
        }
        const defaults = layer.default_controls || {};
        volumeStateByKey[stateKey] = {
          visible: layer.visible !== false,
          vmin: Number(defaults.vmin),
          vmax: Number(defaults.vmax),
          opacity: Number(defaults.opacity),
          steps: Number(defaults.steps),
          alphaCoef: Number(defaults.alpha_coef),
          gradientStep: Number(defaults.gradient_step),
          stretch: normalizeVolumeStretch(defaults.stretch),
          colormap: String(defaults.colormap || (((layer.colormap_options || [])[0] || {}).name || "inferno")),
          showAllTimes: Boolean(defaults.show_all_times),
        };
      });

      legendItems.forEach((item) => {
        const itemKey = String(item.key);
        if (volumeStateKeySet.has(itemKey)) {
          return;
        }
        traceStyleStateByKey[itemKey] = {
          color: String(item.default_color || item.color || theme.axis_color || "#ffffff"),
          opacity: Math.min(Math.max(Number(item.default_opacity ?? 1.0), 0.0), 1.0),
          sizeScale: 1.0,
          hasPoints: Boolean(item.has_points),
          hasSegments: Boolean(item.has_segments),
          hasLabels: Boolean(item.has_labels),
          hasNStars: Boolean(item.has_n_stars),
          sizeByNStarsDefault: Boolean(item.size_by_n_stars_default),
        };
      });

      clusterFilterParameters.forEach((parameter) => {
        const key = String(parameter.key || "");
        if (!key) {
          return;
        }
        clusterFilterRangeStateByKey[key] = {
          min: Number(parameter.min),
          max: Number(parameter.max),
        };
      });
      if (!clusterFilterParameterKey && clusterFilterParameters.length) {
        clusterFilterParameterKey = String(clusterFilterParameters[0].key || "");
      }

      const raycaster = new THREE.Raycaster();
      const pointer = new THREE.Vector2(2, 2);
      const initialZoomAnchorPoint = (
        initialZoomAnchorSpec
        && Number.isFinite(Number(initialZoomAnchorSpec.x))
        && Number.isFinite(Number(initialZoomAnchorSpec.y))
        && Number.isFinite(Number(initialZoomAnchorSpec.z))
      )
        ? new THREE.Vector3(
          Number(initialZoomAnchorSpec.x),
          Number(initialZoomAnchorSpec.y),
          Number(initialZoomAnchorSpec.z)
        )
        : null;
      currentZoomAnchorPoint = initialZoomAnchorPoint ? initialZoomAnchorPoint.clone() : null;
      const initialZoomAnchorTargetSpec = initialState.camera
        && typeof initialState.camera === "object"
        && initialState.camera.target
        && typeof initialState.camera.target === "object"
        ? initialState.camera.target
        : null;
      const initialZoomAnchorCameraSpec = initialState.camera
        && typeof initialState.camera === "object"
        && initialState.camera.position
        && typeof initialState.camera.position === "object"
        ? initialState.camera.position
        : null;
      const initialZoomAnchorTargetPoint = (
        initialZoomAnchorTargetSpec
        && Number.isFinite(Number(initialZoomAnchorTargetSpec.x))
        && Number.isFinite(Number(initialZoomAnchorTargetSpec.y))
        && Number.isFinite(Number(initialZoomAnchorTargetSpec.z))
      )
        ? new THREE.Vector3(
          Number(initialZoomAnchorTargetSpec.x),
          Number(initialZoomAnchorTargetSpec.y),
          Number(initialZoomAnchorTargetSpec.z)
        )
        : null;
      const initialZoomAnchorCameraPoint = (
        initialZoomAnchorCameraSpec
        && Number.isFinite(Number(initialZoomAnchorCameraSpec.x))
        && Number.isFinite(Number(initialZoomAnchorCameraSpec.y))
        && Number.isFinite(Number(initialZoomAnchorCameraSpec.z))
      )
        ? new THREE.Vector3(
          Number(initialZoomAnchorCameraSpec.x),
          Number(initialZoomAnchorCameraSpec.y),
          Number(initialZoomAnchorCameraSpec.z)
        )
        : null;

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

      function clampFrameIndex(value) {
        return Math.max(0, Math.min(Math.round(Number(value) || 0.0), Math.max(frameSpecs.length - 1, 0)));
      }

      function clampFrameValue(value) {
        return Math.max(0.0, Math.min(Number(value) || 0.0, Math.max(frameSpecs.length - 1, 0)));
      }

      function wrapFrameValue(value) {
        const frameCount = Math.max(frameSpecs.length, 0);
        if (frameCount <= 1) {
          return 0.0;
        }
        let wrapped = Number(value) || 0.0;
        wrapped %= frameCount;
        if (wrapped < 0.0) {
          wrapped += frameCount;
        }
        if (wrapped >= frameCount) {
          wrapped = 0.0;
        }
        return wrapped;
      }

      function frameValueState(frameValue) {
        const clampedValue = clampFrameValue(frameValue);
        const lowerIndex = Math.floor(clampedValue);
        const upperIndex = Math.min(lowerIndex + 1, Math.max(frameSpecs.length - 1, 0));
        const alpha = Math.max(0.0, Math.min(clampedValue - lowerIndex, 1.0));
        return {
          value: clampedValue,
          lowerIndex,
          upperIndex,
          alpha,
        };
      }

      function frameTimeForValue(frameValue) {
        const state = frameValueState(frameValue);
        const lowerFrame = frameSpecs[state.lowerIndex] || null;
        const upperFrame = frameSpecs[state.upperIndex] || lowerFrame;
        const lowerTime = Number(lowerFrame && lowerFrame.time);
        const upperTime = Number(upperFrame && upperFrame.time);
        if (!Number.isFinite(lowerTime)) {
          return Number.isFinite(upperTime) ? upperTime : 0.0;
        }
        if (!Number.isFinite(upperTime) || state.upperIndex === state.lowerIndex) {
          return lowerTime;
        }
        return lowerTime + (upperTime - lowerTime) * state.alpha;
      }

      function sceneSpecForRawExport(exportSceneSpec) {
        const payloadSceneSpec = safeJsonClone(exportSceneSpec, {});
        const metadata = payloadSceneSpec.export_metadata && typeof payloadSceneSpec.export_metadata === "object"
          ? payloadSceneSpec.export_metadata
          : null;
        if (metadata && Object.prototype.hasOwnProperty.call(metadata, "scene_spec_payload")) {
          delete metadata.scene_spec_payload;
          if (Object.keys(metadata).length === 0) {
            delete payloadSceneSpec.export_metadata;
          }
        }
        return payloadSceneSpec;
      }

      function rawSceneSpecSizeBytes(jsonText) {
        if (typeof TextEncoder !== "undefined") {
          return new TextEncoder().encode(jsonText).length;
        }
        return String(jsonText || "").length;
      }

      function sceneSpecPayloadMetadataForJsonText(jsonText) {
        return {
          compression_mode: sceneSpecPayloadCompressionMode(),
          compression_method: "none",
          compressed: false,
          raw_scene_spec_size_bytes: rawSceneSpecSizeBytes(jsonText),
          compressed_size_bytes: null,
          embedded_base64_size_bytes: null,
          compression_threshold_bytes: sceneSpecPayloadCompressionThresholdBytes(),
        };
      }

      function sceneSpecPayloadCompressionMode() {
        const mode = sceneSpecPayloadMetadata && sceneSpecPayloadMetadata.compression_mode;
        if (mode === true || mode === false || mode === "auto") {
          return mode;
        }
        if (mode === "true") {
          return true;
        }
        if (mode === "false" || mode === "none") {
          return false;
        }
        return "auto";
      }

      function sceneSpecPayloadCompressionThresholdBytes() {
        const threshold = Number(sceneSpecPayloadMetadata && sceneSpecPayloadMetadata.compression_threshold_bytes);
        if (Number.isFinite(threshold) && threshold >= 0) {
          return Math.floor(threshold);
        }
        return 8 * 1024 * 1024;
      }

      function shouldCompressExportSceneSpec(jsonText) {
        const mode = sceneSpecPayloadCompressionMode();
        if (mode === true) {
          return true;
        }
        if (mode === false) {
          return false;
        }
        return rawSceneSpecSizeBytes(jsonText) > sceneSpecPayloadCompressionThresholdBytes();
      }

      function sceneSpecPayloadMetadataForCompressedJsonText(jsonText, compressedSizeBytes, embeddedBase64SizeBytes) {
        return {
          compression_mode: sceneSpecPayloadCompressionMode(),
          compression_method: "gzip+base64",
          compressed: true,
          raw_scene_spec_size_bytes: rawSceneSpecSizeBytes(jsonText),
          compressed_size_bytes: compressedSizeBytes,
          embedded_base64_size_bytes: embeddedBase64SizeBytes,
          compression_threshold_bytes: sceneSpecPayloadCompressionThresholdBytes(),
        };
      }

      function sceneSpecMarkerWrappedJsonText(jsonText) {
        return `/*__SCENE_SPEC_START__*/const sceneSpec = ${jsonText};/*__SCENE_SPEC_END__*/`;
      }

      function sceneSpecMarkerWrappedCompressedPayload(payloadId) {
        return `/*__SCENE_SPEC_START__*/const sceneSpec = await inflateOvizGzipBase64SceneSpec(readOvizSceneSpecPayload(${JSON.stringify(payloadId)}));/*__SCENE_SPEC_END__*/`;
      }

      function sceneSpecMetadataMarkerWrappedJson(metadata) {
        return `/*__SCENE_SPEC_METADATA_START__*/const sceneSpecPayloadMetadata = ${JSON.stringify(metadata)};/*__SCENE_SPEC_METADATA_END__*/`;
      }

      function sceneSpecPayloadElementId() {
        const idBase = String((root && root.id) || "__ROOT_ID__" || "oviz-three");
        return `${idBase}-scene-spec-payload`;
      }

      function base64EncodeUint8Array(bytes) {
        const chunkSize = 0x6000;
        let encoded = "";
        for (let offset = 0; offset < bytes.length; offset += chunkSize) {
          const chunk = bytes.subarray(offset, offset + chunkSize);
          let binary = "";
          for (let idx = 0; idx < chunk.length; idx += 1) {
            binary += String.fromCharCode(chunk[idx]);
          }
          encoded += window.btoa(binary);
        }
        return encoded;
      }

      async function gzipBase64EncodeText(text) {
        if (typeof CompressionStream === "undefined") {
          return null;
        }
        const stream = new Blob([text])
          .stream()
          .pipeThrough(new CompressionStream("gzip"));
        const buffer = await new Response(stream).arrayBuffer();
        const bytes = new Uint8Array(buffer);
        const encoded = base64EncodeUint8Array(bytes);
        return {
          encoded,
          compressedSizeBytes: bytes.length,
          embeddedBase64SizeBytes: encoded.length,
        };
      }

      function sceneSpecPayloadScriptHtml(payloadId, encodedText) {
        const encodedChunks = String(encodedText || "").match(/.{1,65536}/g) || [""];
        return `<script type="application/octet-stream" id="${payloadId}">\n${encodedChunks.join("\\n")}\n<\\/script>`;
      }

      function removeExistingSceneSpecPayloadHtml(htmlText) {
        return String(htmlText || "").replace(
          /\\s*<script\\b(?=[^>]*\\btype=["']application\\/octet-stream["'])(?=[^>]*\\bid=["'][^"']*-scene-spec-payload["'])[^>]*>[\\s\\S]*?<\\/script>\\s*/gi,
          "\\n"
        );
      }

      function insertSceneSpecPayloadHtml(htmlText, payloadHtml) {
        if (!payloadHtml) {
          return htmlText;
        }
        const markerIndex = htmlText.indexOf("/*__SCENE_SPEC_START__*/");
        const scriptStartIndex = markerIndex >= 0 ? htmlText.lastIndexOf("<script", markerIndex) : -1;
        if (scriptStartIndex < 0) {
          throw new Error("Could not locate the main viewer script for the compressed scene payload.");
        }
        return htmlText.slice(0, scriptStartIndex) + payloadHtml + "\\n" + htmlText.slice(scriptStartIndex);
      }

      function safeJsonClone(value, fallback = null) {
        try {
          return JSON.parse(JSON.stringify(value));
        } catch (_err) {
          return fallback;
        }
      }

      function escapeHtml(value) {
        return String(value ?? "")
          .replace(/&/g, "&amp;")
          .replace(/</g, "&lt;")
          .replace(/>/g, "&gt;")
          .replace(/"/g, "&quot;")
          .replace(/'/g, "&#39;");
      }

      function measuredScaleBarSize() {
        if (!scaleBarEl) {
          return { width: 0.0, height: 0.0 };
        }
        const rect = scaleBarEl.getBoundingClientRect();
        const fallbackWidth = Number(scaleBarEl.offsetWidth) || 160.0;
        const fallbackHeight = Number(scaleBarEl.offsetHeight) || 44.0;
        return {
          width: Math.max(Number(rect.width) || fallbackWidth, 1.0),
          height: Math.max(Number(rect.height) || fallbackHeight, 1.0),
        };
      }

      function clampScaleBarPosition(left, top, width, height) {
        const rootRect = root.getBoundingClientRect();
        const margin = 6.0;
        return {
          left: Math.min(Math.max(margin, Number(left) || 0.0), Math.max(margin, rootRect.width - width - margin)),
          top: Math.min(Math.max(margin, Number(top) || 0.0), Math.max(margin, rootRect.height - height - margin)),
        };
      }

      function defaultScaleBarPosition() {
        if (!scaleBarEl) {
          return { left: 18.0, top: 18.0 };
        }
        const rootRect = root.getBoundingClientRect();
        const size = measuredScaleBarSize();
        const left = 18.0;
        const top = Math.max(6.0, rootRect.height - size.height - 22.0);
        return clampScaleBarPosition(left, top, size.width, size.height);
      }

      function clampLegendPanelRect(left, top, width, height) {
        const margin = 8.0;
        const minWidth = minimalModeEnabled
          ? Math.min(140.0, Math.max(96.0, root.clientWidth - 2 * margin))
          : Math.min(196.0, Math.max(160.0, root.clientWidth - 2 * margin));
        const minHeight = minimalModeEnabled
          ? Math.min(48.0, Math.max(18.0, root.clientHeight - 2 * margin))
          : Math.min(132.0, Math.max(72.0, root.clientHeight - 2 * margin));
        const nextWidth = clampRange(Number(width) || minWidth, minWidth, Math.max(minWidth, root.clientWidth - 2 * margin));
        const nextHeight = clampRange(Number(height) || minHeight, minHeight, Math.max(minHeight, root.clientHeight - 2 * margin));
        return {
          left: clampRange(Number(left) || margin, margin, Math.max(margin, root.clientWidth - nextWidth - margin)),
          top: clampRange(Number(top) || margin, margin, Math.max(margin, root.clientHeight - nextHeight - margin)),
          width: nextWidth,
          height: nextHeight,
        };
      }

      function defaultLegendPanelRect() {
        const actionBarHeight = actionBarEl && actionBarEl.dataset.visible === "true"
          ? Math.max(actionBarEl.getBoundingClientRect().height, 0.0)
          : 0.0;
        const topInset = Math.max(14.0 + actionBarHeight + (actionBarHeight > 0.0 ? 10.0 : 0.0), 12.0);
        if (minimalModeEnabled) {
          const visibleItems = visibleLegendItemsForCurrentGroup();
          const estimatedHeight = Math.max(visibleItems.length * 16 + 8, 24);
          const width = Math.min(Math.max(root.clientWidth * 0.18, 140.0), 240.0);
          const height = Math.min(estimatedHeight, Math.max(root.clientHeight - topInset - 12.0, 24.0));
          return clampLegendPanelRect(12.0, topInset, width, height);
        }
        const width = Math.min(Math.max(root.clientWidth * 0.18, 204.0), 240.0);
        const visibleItems = visibleLegendItemsForCurrentGroup();
        const traceCount = visibleItems.filter((item) => !volumeLayerForKey(item.key)).length;
        const volumeCount = visibleItems.length - traceCount;
        const preferredTraceRows = traceCount ? clampRange(traceCount, 4, 8) : 0;
        const preferredVolumeRows = volumeCount ? clampRange(volumeCount, 1, 3) : 0;
        const estimatedHeight = 108
          + preferredTraceRows * 34
          + preferredVolumeRows * 32
          + ((traceCount && volumeCount) ? 14 : 0);
        const minHeight = volumeCount ? 304.0 : 272.0;
        const height = Math.min(Math.max(estimatedHeight, minHeight), Math.min(root.clientHeight - topInset - 12.0, 520.0));
        return clampLegendPanelRect(14.0, topInset, width, height);
      }

      function legendSectionDefs() {
        return [
          {
            key: "traces",
            sectionEl: legendTraceSectionEl,
            toggleEl: legendTraceSectionToggleEl,
            listEl: legendTraceListEl,
            label: "Traces",
          },
          {
            key: "volumes",
            sectionEl: legendVolumeSectionEl,
            toggleEl: legendVolumeSectionToggleEl,
            listEl: legendVolumeListEl,
            label: "Volumes",
          },
        ];
      }

      function applyLegendSectionState(sectionKey) {
        const section = legendSectionDefs().find((item) => item.key === sectionKey);
        if (!section || !section.sectionEl || !section.toggleEl) {
          return;
        }
        const isOpen = minimalModeEnabled ? true : (legendSectionOpenState[sectionKey] !== false);
        section.sectionEl.dataset.open = isOpen ? "true" : "false";
        section.toggleEl.setAttribute("aria-expanded", isOpen ? "true" : "false");
        section.toggleEl.setAttribute(
          "title",
          `${isOpen ? "Collapse" : "Expand"} ${String(section.label || sectionKey).toLowerCase()}`
        );
        const chevronEl = section.toggleEl.querySelector(".oviz-three-legend-section-chevron");
        if (chevronEl) {
          chevronEl.textContent = isOpen ? "▾" : "▸";
        }
      }

      function setLegendSectionOpen(sectionKey, isOpen) {
        legendSectionOpenState[sectionKey] = Boolean(isOpen);
        applyLegendSectionState(sectionKey);
        if (legendPanelOpen) {
          applyLegendPanelRect(legendPanelRectState || defaultLegendPanelRect());
        }
      }

      function legendPanelPreferredExpandedHeight(maxVisibleEntries = LEGEND_MAX_VISIBLE_ITEMS) {
        if (!legendPanelEl || !legendPanelBodyEl || !legendPanelOpen) {
          return 0.0;
        }
        const headerHeight = legendDragHandleEl
          ? Math.max(legendDragHandleEl.getBoundingClientRect().height || 0.0, 36.0)
          : 36.0;
        const hiddenEntries = [];
        const hiddenSections = [];
        let remainingEntries = Math.max(0, Number(maxVisibleEntries) || 0);

        legendSectionDefs().forEach(({ sectionEl, listEl }) => {
          if (!sectionEl || sectionEl.dataset.empty === "true") {
            return;
          }
          if (sectionEl.dataset.open === "false" || !listEl) {
            return;
          }
          const entries = Array.from(listEl.children).filter(
            (child) => child && child.classList && child.classList.contains("oviz-three-legend-entry")
          );
          const visibleCount = Math.min(entries.length, remainingEntries);
          if (visibleCount <= 0) {
            hiddenSections.push([sectionEl, sectionEl.style.display]);
            sectionEl.style.display = "none";
            return;
          }
          remainingEntries -= visibleCount;
          entries.slice(visibleCount).forEach((entry) => {
            hiddenEntries.push([entry, entry.style.display]);
            entry.style.display = "none";
          });
        });

        const bodyHeight = Math.max(legendPanelBodyEl.scrollHeight || 0.0, 0.0);

        hiddenEntries.forEach(([entry, previousDisplay]) => {
          entry.style.display = previousDisplay;
        });
        hiddenSections.forEach(([sectionEl, previousDisplay]) => {
          sectionEl.style.display = previousDisplay;
        });

        return headerHeight + bodyHeight + 2.0;
      }

      function legendPanelOverflowHeight(expandedHeight) {
        if (!legendPanelEl || !legendPanelBodyEl || !legendPanelOpen) {
          return 0.0;
        }
        const headerHeight = legendDragHandleEl
          ? Math.max(legendDragHandleEl.getBoundingClientRect().height || 0.0, 36.0)
          : 36.0;
        const availableBodyHeight = Math.max(Number(expandedHeight) - headerHeight, 0.0);
        const contentBodyHeight = Math.max(legendPanelBodyEl.scrollHeight || 0.0, 0.0);
        return Math.max(contentBodyHeight - availableBodyHeight, 0.0);
      }

      function applyLegendPanelRect(rectState, options = {}) {
        if (!legendPanelEl) {
          return;
        }
        const allowAutoCap = options.allowAutoCap !== false;
        const requestedHeight = rectState && rectState.height;
        const base = clampLegendPanelRect(
          rectState && rectState.left,
          rectState && rectState.top,
          rectState && rectState.width,
          requestedHeight,
        );
        let next = base;
        if (legendPanelOpen && allowAutoCap && !legendPanelUserSized) {
          const preferredExpandedHeight = legendPanelPreferredExpandedHeight(LEGEND_MAX_VISIBLE_ITEMS);
          let targetHeight = preferredExpandedHeight > 0.0
            ? Math.min(base.height, preferredExpandedHeight)
            : base.height;
          const overflowHeight = legendPanelOverflowHeight(targetHeight);
          if (overflowHeight > 1.0) {
            targetHeight = preferredExpandedHeight > 0.0
              ? Math.min(targetHeight + overflowHeight + 2.0, preferredExpandedHeight)
              : targetHeight + overflowHeight + 2.0;
          }
          if (Math.abs(targetHeight - base.height) > 0.5) {
            next = clampLegendPanelRect(
              base.left,
              base.top,
              base.width,
              targetHeight,
            );
          }
        }
        legendPanelRectState = next;
        legendPanelEl.style.left = `${next.left}px`;
        legendPanelEl.style.top = `${next.top}px`;
        legendPanelEl.style.width = `${next.width}px`;
        if (legendPanelOpen) {
          legendPanelEl.style.height = `${next.height}px`;
          legendPanelEl.dataset.expandedHeight = String(next.height);
        } else {
          const headerHeight = legendDragHandleEl
            ? Math.max(legendDragHandleEl.getBoundingClientRect().height || 0.0, 36.0)
            : 36.0;
          legendPanelEl.style.height = `${headerHeight}px`;
          legendPanelEl.dataset.expandedHeight = String(next.height);
        }
        legendPanelEl.style.right = "auto";
        legendPanelEl.style.bottom = "auto";
      }

      function captureLegendPanelState() {
        if (!legendPanelEl) {
          return null;
        }
        if (!legendPanelRectState) {
          const rootRect = root.getBoundingClientRect();
          const rect = legendPanelEl.getBoundingClientRect();
          const expandedHeight = Number(legendPanelEl.dataset.expandedHeight);
          legendPanelRectState = clampLegendPanelRect(
            rect.left - rootRect.left,
            rect.top - rootRect.top,
            rect.width,
            Number.isFinite(expandedHeight) ? expandedHeight : rect.height,
          );
        }
        return safeJsonClone(legendPanelRectState, null);
      }

      function resizeLegendRect(state, clientX, clientY) {
        const margin = 8.0;
        const minWidth = Math.min(196.0, Math.max(160.0, root.clientWidth - 2 * margin));
        const minHeight = Math.min(112.0, Math.max(64.0, root.clientHeight - 2 * margin));
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
          right = Math.min(root.clientWidth - margin, Math.max(state.startRight + dx, state.startLeft + minWidth));
        }
        if (dir.includes("n")) {
          top = Math.max(margin, Math.min(state.startTop + dy, state.startBottom - minHeight));
        } else if (dir.includes("s")) {
          bottom = Math.min(root.clientHeight - margin, Math.max(state.startBottom + dy, state.startTop + minHeight));
        }

        return clampLegendPanelRect(left, top, right - left, bottom - top);
      }

      function niceFrameStep(rawStep) {
        const numericStep = Math.max(Number(rawStep) || 1, 1);
        const exponent = Math.floor(Math.log10(numericStep));
        const magnitude = Math.pow(10, exponent);
        const normalized = numericStep / magnitude;
        if (normalized <= 1) {
          return 1 * magnitude;
        }
        if (normalized <= 2) {
          return 2 * magnitude;
        }
        if (normalized <= 5) {
          return 5 * magnitude;
        }
        return 10 * magnitude;
      }

      function formatTimelineTickLabel(frame) {
        if (!frame) {
          return "";
        }
        const numericTime = Number(frame.time);
        if (Number.isFinite(numericTime)) {
          return formatTick(numericTime).replace(/^-/, "−");
        }
        return String(frame.name || "");
      }

      function beginLegendPanelInteraction(state, event) {
        legendPointerState = state;
        if (legendPanelEl) {
          legendPanelEl.dataset.dragging = "true";
        }
        document.body.style.userSelect = "none";
        if (event) {
          event.preventDefault();
        }
      }

      function onLegendPointerStart(event) {
        if (!legendPanelEl || !legendDragHandleEl || event.button !== 0) {
          return false;
        }
        const target = event.target;
        if (!target) {
          return false;
        }
        if (target.closest && target.closest(".oviz-three-legend-panel-toggle")) {
          return false;
        }
        const panelRect = legendPanelEl.getBoundingClientRect();
        const rootRect = root.getBoundingClientRect();
        const startLeft = panelRect.left - rootRect.left;
        const startTop = panelRect.top - rootRect.top;
        const startWidth = panelRect.width;
        const startHeight = Number(legendPanelEl.dataset.expandedHeight) || panelRect.height;
        const resizeHandle = target.closest ? target.closest(".oviz-three-legend-resize") : null;
        if (resizeHandle) {
          if (!legendPanelOpen) {
            return false;
          }
          beginLegendPanelInteraction({
            mode: "resize",
            dir: String(resizeHandle.dataset.dir || "se"),
            startX: event.clientX,
            startY: event.clientY,
            startLeft,
            startTop,
            startRight: startLeft + startWidth,
            startBottom: startTop + startHeight,
          }, event);
          return true;
        }
        if (target !== legendDragHandleEl && !(legendDragHandleEl.contains && legendDragHandleEl.contains(target))) {
          return false;
        }
        beginLegendPanelInteraction({
          mode: "move",
          startX: event.clientX,
          startY: event.clientY,
          offsetX: event.clientX - startLeft,
          offsetY: event.clientY - startTop,
          width: startWidth,
          height: startHeight,
        }, event);
        return true;
      }

      function updateLegendPanelInteraction(event) {
        if (!legendPointerState) {
          return false;
        }
        let nextRect = null;
        if (legendPointerState.mode === "move") {
          nextRect = clampLegendPanelRect(
            event.clientX - legendPointerState.offsetX,
            event.clientY - legendPointerState.offsetY,
            legendPointerState.width,
            legendPointerState.height,
          );
        } else if (legendPointerState.mode === "resize") {
          legendPanelUserSized = true;
          nextRect = resizeLegendRect(legendPointerState, event.clientX, event.clientY);
        }
        if (nextRect) {
          applyLegendPanelRect(nextRect, { allowAutoCap: false });
          if (typeof positionLegendPopover === "function" && activeLegendEditorKey && legendEditButtonByKey.has(activeLegendEditorKey)) {
            positionLegendPopover(legendEditButtonByKey.get(activeLegendEditorKey));
          }
        }
        event.preventDefault();
        return true;
      }

      function finishLegendPanelInteraction(event) {
        if (!legendPointerState) {
          return false;
        }
        legendPointerState = null;
        if (legendPanelEl) {
          legendPanelEl.dataset.dragging = "false";
        }
        document.body.style.userSelect = "";
        if (event) {
          event.preventDefault();
        }
        return true;
      }

      function applyScaleBarPosition() {
        if (!scaleBarEl) {
          return;
        }
        const size = measuredScaleBarSize();
        const next = defaultScaleBarPosition();
        scaleBarEl.style.left = `${next.left}px`;
        scaleBarEl.style.top = `${next.top}px`;
        scaleBarEl.style.right = "auto";
        scaleBarEl.style.bottom = "auto";
        scaleBarPosition = null;
      }

      function captureScaleBarState() {
        return null;
      }

      function buildSelectionBoxMetricTimeCenters() {
        if (!selectionBoxMetricsSpec.enabled) {
          return [];
        }
        const lookbackMax = Math.max(Number(selectionBoxMetricsSpec.lookback_max_myr) || 0.0, 0.0);
        const stepMyr = Math.max(Number(selectionBoxMetricsSpec.step_myr) || 0.0, 0.0);
        if (!(lookbackMax > 0.0) || !(stepMyr > 0.0)) {
          return [];
        }
        const centers = [];
        for (let center = 0.0; center <= lookbackMax + 1e-9; center += stepMyr) {
          centers.push(Number(center.toFixed(6)));
        }
        if (!centers.length || Math.abs(centers[centers.length - 1] - lookbackMax) > 1e-6) {
          centers.push(Number(lookbackMax.toFixed(6)));
        }
        return centers;
      }

      function normalizeSelectionBoxMetricEvents(eventSpec) {
        const rawTimes = Array.isArray(eventSpec.time_of_death_myr) ? eventSpec.time_of_death_myr : [];
        const rawX = Array.isArray(eventSpec.x_pc) ? eventSpec.x_pc : [];
        const rawY = Array.isArray(eventSpec.y_pc) ? eventSpec.y_pc : [];
        const rawZ = Array.isArray(eventSpec.z_pc) ? eventSpec.z_pc : [];
        const rawSampleIds = Array.isArray(eventSpec.imf_sample_idx) ? eventSpec.imf_sample_idx : [];
        const limit = Math.min(rawTimes.length, rawX.length, rawY.length, rawZ.length, rawSampleIds.length);
        const lookbackMax = Math.max(Number(selectionBoxMetricsSpec.lookback_max_myr) || 0.0, 0.0);
        const rows = [];
        const sampleIdSet = new Set();
        for (let index = 0; index < limit; index += 1) {
          const time = Number(rawTimes[index]);
          const lookback = -time;
          const x = Number(rawX[index]);
          const y = Number(rawY[index]);
          const z = Number(rawZ[index]);
          const sampleId = Math.round(Number(rawSampleIds[index]));
          if (
            !Number.isFinite(lookback)
            || !Number.isFinite(x)
            || !Number.isFinite(y)
            || !Number.isFinite(z)
            || !Number.isFinite(sampleId)
          ) {
            continue;
          }
          if (lookback < -1e-9 || lookback > lookbackMax + 1e-9) {
            continue;
          }
          rows.push({ lookback, x, y, z, sampleId });
          sampleIdSet.add(sampleId);
        }
        rows.sort((a, b) => a.lookback - b.lookback);
        const sampleValues = Array.from(sampleIdSet).sort((a, b) => a - b);
        const sampleIndexById = new Map(sampleValues.map((value, index) => [value, index]));
        const lookback = new Float32Array(rows.length);
        const x = new Float32Array(rows.length);
        const y = new Float32Array(rows.length);
        const z = new Float32Array(rows.length);
        const sampleIds = new Int32Array(rows.length);
        const sampleIndex = new Int32Array(rows.length);
        rows.forEach((row, index) => {
          lookback[index] = row.lookback;
          x[index] = row.x;
          y[index] = row.y;
          z[index] = row.z;
          sampleIds[index] = row.sampleId;
          sampleIndex[index] = sampleIndexById.get(row.sampleId);
        });
        return {
          count: rows.length,
          lookback,
          x,
          y,
          z,
          sampleIds,
          sampleIndex,
          sampleValues,
        };
      }

      function buildSelectionBoxRandomUnitCatalog(size, seed) {
        const nPoints = Math.max(0, Math.round(Number(size) || 0));
        const x = new Float32Array(nPoints);
        const y = new Float32Array(nPoints);
        const z = new Float32Array(nPoints);
        let state = (Math.round(Number(seed) || 1) >>> 0) || 1;
        function nextUnit() {
          state = (state + 0x6D2B79F5) >>> 0;
          let t = state;
          t = Math.imul(t ^ (t >>> 15), t | 1);
          t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
          return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
        }
        for (let index = 0; index < nPoints; index += 1) {
          x[index] = 2.0 * nextUnit() - 1.0;
          y[index] = 2.0 * nextUnit() - 1.0;
          z[index] = 2.0 * nextUnit() - 1.0;
        }
        return { x, y, z };
      }

      function selectionBoxHalfWidthBounds() {
        const minHalfWidthPc = Math.max(Number(selectionBoxSpec.min_half_width_pc) || 50.0, 10.0);
        const maxHalfWidthPc = Math.max(
          Number(selectionBoxSpec.max_half_width_pc) || Math.max(minHalfWidthPc, sceneSpec.max_span || 1000.0),
          minHalfWidthPc
        );
        return { minHalfWidthPc, maxHalfWidthPc };
      }

      function clampSelectionBoxHalfWidth(value) {
        const bounds = selectionBoxHalfWidthBounds();
        return clampRange(Number(value) || bounds.minHalfWidthPc, bounds.minHalfWidthPc, bounds.maxHalfWidthPc);
      }

      function selectionBoxBandBounds() {
        const { maxHalfWidthPc } = selectionBoxHalfWidthBounds();
        const minBandPc = Math.max(Number(selectionBoxMetricsSpec.twopcf_min_pc) || 1.0, 0.0);
        const maxBandPc = Math.max(
          Number(selectionBoxMetricsSpec.twopcf_max_pc) || (2.0 * Math.sqrt(3.0) * maxHalfWidthPc),
          minBandPc + 1.0
        );
        return { minBandPc, maxBandPc };
      }

      function normalizeSelectionBoxBand(bandLike) {
        const bounds = selectionBoxBandBounds();
        const defaultMin = Number(selectionBoxDefaultBandSpec[0]) || bounds.minBandPc;
        const defaultMax = Number(selectionBoxDefaultBandSpec[1]) || Math.max(defaultMin + 1.0, bounds.maxBandPc);
        let minValue;
        let maxValue;
        if (Array.isArray(bandLike)) {
          minValue = Number(bandLike[0]);
          maxValue = Number(bandLike[1]);
        } else {
          minValue = Number(bandLike && bandLike.min);
          maxValue = Number(bandLike && bandLike.max);
        }
        if (!Number.isFinite(minValue)) {
          minValue = defaultMin;
        }
        if (!Number.isFinite(maxValue)) {
          maxValue = defaultMax;
        }
        if (maxValue < minValue) {
          const swappedMin = maxValue;
          maxValue = minValue;
          minValue = swappedMin;
        }
        minValue = clampRange(minValue, bounds.minBandPc, bounds.maxBandPc - 1.0);
        maxValue = clampRange(maxValue, minValue + 1.0, bounds.maxBandPc);
        return [minValue, maxValue];
      }

      function activeSelectionBoxBandPc() {
        return normalizeSelectionBoxBand(selectionBoxState.twopcfBandPc || selectionBoxDefaultBandSpec);
      }

      function formatSelectionBoxBandInputValue(value) {
        const numericValue = Number(value);
        if (!Number.isFinite(numericValue)) {
          return "";
        }
        return Number.isInteger(numericValue) ? String(numericValue) : numericValue.toFixed(1);
      }

      function syncSelectionBoxBandInputs(force = false) {
        const band = activeSelectionBoxBandPc();
      }

      function syncSelectionBoxVisibilityInput(force = false) {
        if (!boxMetricsVisibleEl) {
          return;
        }
        if (force || document.activeElement !== boxMetricsVisibleEl) {
          boxMetricsVisibleEl.checked = selectionBoxState.visible !== false;
        }
      }

      function updateSelectionBoxBandFromInputs(immediate = false) {
        const currentBand = activeSelectionBoxBandPc();
        const proposedBand = normalizeSelectionBoxBand([
          boxMetricsBandMinEl ? Number(boxMetricsBandMinEl.value) : currentBand[0],
          boxMetricsBandMaxEl ? Number(boxMetricsBandMaxEl.value) : currentBand[1],
        ]);
        selectionBoxState.twopcfBandPc = proposedBand;
        syncSelectionBoxBandInputs();
        renderBoxMetricsWidget();
        scheduleSelectionBoxMetricsRecompute(immediate);
      }

      function buildDefaultSelectionBoxState() {
        const defaultCenter = selectionBoxSpec.default_center_local_pc || {};
        return {
          center: {
            x: Number(defaultCenter.x) || 0.0,
            y: Number(defaultCenter.y) || 0.0,
            z: Number(defaultCenter.z) || 0.0,
          },
          visible: selectionBoxSpec.default_visible !== false,
          halfWidthPc: clampSelectionBoxHalfWidth(
            Number(selectionBoxSpec.default_half_width_pc) || 500.0
          ),
          twopcfBandPc: normalizeSelectionBoxBand(selectionBoxDefaultBandSpec),
        };
      }

      function selectionBoxPanelIsOpen() {
        return widgetModeForKey("box_metrics") !== "hidden";
      }

      function selectionBoxShouldRender() {
        return !minimalModeEnabled
          && Boolean(selectionBoxSpec.enabled)
          && selectionBoxPanelIsOpen()
          && selectionBoxState.visible !== false;
      }

      function selectionBoxReferenceOrbitAtTime(timeMyr) {
        const times = Array.isArray(selectionBoxReferenceOrbitSpec.time_myr)
          ? selectionBoxReferenceOrbitSpec.time_myr
          : [];
        const xs = Array.isArray(selectionBoxReferenceOrbitSpec.x_gc_pc)
          ? selectionBoxReferenceOrbitSpec.x_gc_pc
          : [];
        const ys = Array.isArray(selectionBoxReferenceOrbitSpec.y_gc_pc)
          ? selectionBoxReferenceOrbitSpec.y_gc_pc
          : [];
        const zs = Array.isArray(selectionBoxReferenceOrbitSpec.z_gc_pc)
          ? selectionBoxReferenceOrbitSpec.z_gc_pc
          : [];
        const nSamples = Math.min(times.length, xs.length, ys.length, zs.length);
        if (!nSamples) {
          return { x: 0.0, y: 0.0, z: 0.0 };
        }
        const targetTime = Number(timeMyr) || 0.0;
        const firstTime = Number(times[0]);
        const lastTime = Number(times[nSamples - 1]);
        if (targetTime <= firstTime) {
          return {
            x: Number(xs[0]) || 0.0,
            y: Number(ys[0]) || 0.0,
            z: Number(zs[0]) || 0.0,
          };
        }
        if (targetTime >= lastTime) {
          return {
            x: Number(xs[nSamples - 1]) || 0.0,
            y: Number(ys[nSamples - 1]) || 0.0,
            z: Number(zs[nSamples - 1]) || 0.0,
          };
        }
        let rightIndex = 1;
        while (rightIndex < nSamples && Number(times[rightIndex]) < targetTime) {
          rightIndex += 1;
        }
        const leftIndex = Math.max(0, rightIndex - 1);
        const leftTime = Number(times[leftIndex]);
        const rightTime = Number(times[rightIndex]);
        const fraction = Math.abs(rightTime - leftTime) <= 1e-12
          ? 0.0
          : (targetTime - leftTime) / (rightTime - leftTime);
        function interp(values) {
          const leftValue = Number(values[leftIndex]) || 0.0;
          const rightValue = Number(values[rightIndex]) || leftValue;
          return leftValue + (rightValue - leftValue) * fraction;
        }
        return {
          x: interp(xs),
          y: interp(ys),
          z: interp(zs),
        };
      }

      function selectionBoxRotatingLocalToDisplay(localPoint, timeMyr) {
        const roKpc = Number(selectionBoxTransformSpec.ro_kpc) || 8.0;
        const voKms = Number(selectionBoxTransformSpec.vo_kms) || 220.0;
        const timeCodeFactor = Number(selectionBoxTransformSpec.time_code_factor) || 0.01022;
        const rSunPc = 1000.0 * roKpc;
        const omegaPerMyr = (voKms / Math.max(roKpc, 1e-9)) / 10.0;
        const xRot = Number(localPoint && localPoint.x) || 0.0;
        const yRot = Number(localPoint && localPoint.y) || 0.0;
        const zRot = Number(localPoint && localPoint.z) || 0.0;
        const localX = rSunPc - xRot;
        const localY = yRot;
        const rPc = Math.hypot(localX, localY);
        const phi = Math.atan2(localY, localX);
        const theta = phi + omegaPerMyr * (Number(timeMyr) || 0.0) * timeCodeFactor - 0.5 * Math.PI;
        const xGc = rPc * Math.sin(theta);
        const yGc = rPc * Math.cos(theta);
        const refOrbit = selectionBoxReferenceOrbitAtTime(timeMyr);
        return new THREE.Vector3(
          xGc - refOrbit.x,
          yGc - refOrbit.y,
          zRot - refOrbit.z
        );
      }

      function selectionBoxDisplayToRotatingLocal(displayPoint, timeMyr) {
        const roKpc = Number(selectionBoxTransformSpec.ro_kpc) || 8.0;
        const voKms = Number(selectionBoxTransformSpec.vo_kms) || 220.0;
        const timeCodeFactor = Number(selectionBoxTransformSpec.time_code_factor) || 0.01022;
        const rSunPc = 1000.0 * roKpc;
        const omegaPerMyr = (voKms / Math.max(roKpc, 1e-9)) / 10.0;
        const refOrbit = selectionBoxReferenceOrbitAtTime(timeMyr);
        const xGc = (Number(displayPoint && displayPoint.x) || 0.0) + refOrbit.x;
        const yGc = (Number(displayPoint && displayPoint.y) || 0.0) + refOrbit.y;
        const zGc = (Number(displayPoint && displayPoint.z) || 0.0) + refOrbit.z;
        const rPc = Math.hypot(xGc, yGc);
        const theta = Math.atan2(xGc, yGc);
        const phi = theta - omegaPerMyr * (Number(timeMyr) || 0.0) * timeCodeFactor + 0.5 * Math.PI;
        const localX = rPc * Math.cos(phi);
        const localY = rPc * Math.sin(phi);
        return {
          x: rSunPc - localX,
          y: localY,
          z: zGc,
        };
      }

      function selectionBoxDisplayStateAtTime(timeMyr) {
        if (!selectionBoxSpec.enabled) {
          return null;
        }
        const centerLocal = selectionBoxState.center || { x: 0.0, y: 0.0, z: 0.0 };
        const halfWidthPc = clampSelectionBoxHalfWidth(selectionBoxState.halfWidthPc);
        const signs = [
          [-1, -1, -1], [1, -1, -1], [1, 1, -1], [-1, 1, -1],
          [-1, -1, 1], [1, -1, 1], [1, 1, 1], [-1, 1, 1],
        ];
        const cornersDisplay = signs.map((signsRow) => selectionBoxRotatingLocalToDisplay(
          {
            x: centerLocal.x + signsRow[0] * halfWidthPc,
            y: centerLocal.y + signsRow[1] * halfWidthPc,
            z: centerLocal.z + signsRow[2] * halfWidthPc,
          },
          timeMyr
        ));
        const centerDisplay = selectionBoxRotatingLocalToDisplay(centerLocal, timeMyr);
        const basisXDisplay = selectionBoxRotatingLocalToDisplay(
          { x: centerLocal.x + 1.0, y: centerLocal.y, z: centerLocal.z },
          timeMyr
        ).sub(centerDisplay.clone());
        const angle = Math.atan2(basisXDisplay.y, basisXDisplay.x);
        const edgePairs = [
          [0, 1], [1, 2], [2, 3], [3, 0],
          [4, 5], [5, 6], [6, 7], [7, 4],
          [0, 4], [1, 5], [2, 6], [3, 7],
        ];
        const segments = edgePairs.map(([startIndex, endIndex]) => {
          const start = cornersDisplay[startIndex];
          const end = cornersDisplay[endIndex];
          return [start.x, start.y, start.z, end.x, end.y, end.z];
        });
        return {
          timeMyr: Number(timeMyr) || 0.0,
          centerLocal: { x: centerLocal.x, y: centerLocal.y, z: centerLocal.z },
          halfWidthPc,
          centerDisplay,
          cornersDisplay,
          segments,
          angle,
        };
      }

      function selectionBoxSummaryText(cache, isPending = false) {
        const state = cache || {
          centerLocal: selectionBoxState.center,
          halfWidthPc: selectionBoxState.halfWidthPc,
          selectedEventCount: 0,
          rate: null,
          clustering: null,
        };
        const center = state.centerLocal || { x: 0.0, y: 0.0, z: 0.0 };
        const sideKpc = (2.0 * (Number(state.halfWidthPc) || 0.0)) / 1000.0;
        const volumeKpc3 = Number(state.volumeKpc3) > 0.0 ? Number(state.volumeKpc3) : Math.pow(sideKpc, 3);
        const frame = currentFrame();
        const currentLookback = Math.abs(Number(frame && frame.time) || 0.0);
        let currentRateText = "";
        let currentClusteringText = "";
        if (state.rate && Array.isArray(state.rate.centers) && Array.isArray(state.rate.mean) && state.rate.centers.length) {
          let bestIndex = 0;
          let bestDelta = Infinity;
          state.rate.centers.forEach((centerValue, index) => {
            const delta = Math.abs(Number(centerValue) - currentLookback);
            if (delta < bestDelta) {
              bestDelta = delta;
              bestIndex = index;
            }
          });
          const rateValue = Number(state.rate.mean[bestIndex]);
          if (Number.isFinite(rateValue)) {
            currentRateText = `\\nCurrent-frame mean ccSN density: ${rateValue.toFixed(2)} Myr^-1 kpc^-3`;
          }
        }
        if (
          state.clustering
          && Array.isArray(state.clustering.centers)
          && Array.isArray(state.clustering.median)
          && state.clustering.centers.length
        ) {
          let bestIndex = 0;
          let bestDelta = Infinity;
          state.clustering.centers.forEach((centerValue, index) => {
            const delta = Math.abs(Number(centerValue) - currentLookback);
            if (delta < bestDelta) {
              bestDelta = delta;
              bestIndex = index;
            }
          });
          const clusteringValue = Number(state.clustering.median[bestIndex]);
          if (Number.isFinite(clusteringValue)) {
            currentClusteringText = `\\nCurrent-frame NN enhancement: ${clusteringValue.toFixed(2)}`;
          }
        }
        const pendingText = isPending ? "\\nUpdating metrics..." : "";
        return [
          `Center (pc): (${formatTick(center.x)}, ${formatTick(center.y)}, ${formatTick(center.z)})`,
          `Half-width: ${formatTick(state.halfWidthPc)} pc | Volume: ${volumeKpc3.toFixed(2)} kpc^3 | Events: ${formatCompactNumber(state.selectedEventCount || 0)}`,
        ].join("\\n") + currentRateText + currentClusteringText + pendingText;
      }

      function selectionBoxRayHitFromEvent(event) {
        if (!selectionBoxSpec.enabled || !selectionBoxHitObjects.length) {
          return null;
        }
        pointerRayFromEvent(event);
        const hits = raycaster.intersectObjects(selectionBoxHitObjects, false);
        if (!hits.length) {
          return null;
        }
        const hit = hits[0];
        const handle = hit && hit.object && hit.object.userData
          ? hit.object.userData.selectionBoxHandle
          : null;
        return handle ? { hit, handle } : null;
      }

      function selectionBoxPlanePointFromEvent(event, planeZ) {
        pointerRayFromEvent(event);
        const plane = new THREE.Plane().setFromNormalAndCoplanarPoint(
          new THREE.Vector3(0.0, 0.0, 1.0),
          new THREE.Vector3(0.0, 0.0, Number(planeZ) || 0.0)
        );
        const point = new THREE.Vector3();
        return raycaster.ray.intersectPlane(plane, point) ? point : null;
      }

      function startSelectionBoxInteraction(event) {
        if (!selectionBoxSpec.enabled || widgetPointerState || lassoState || event.button !== 0) {
          return false;
        }
        const hitInfo = selectionBoxRayHitFromEvent(event);
        if (!hitInfo) {
          return false;
        }
        const frame = currentFrame();
        const timeMyr = Number(frame && frame.time) || 0.0;
        const displayState = selectionBoxDisplayStateAtTime(timeMyr);
        if (!displayState) {
          return false;
        }
        const sceneCenter = displayState.centerDisplay.clone().add(plotGroup.position);
        const planePointScene = selectionBoxPlanePointFromEvent(event, sceneCenter.z) || hitInfo.hit.point.clone();
        const planePointDisplay = planePointScene.clone().sub(plotGroup.position);
        const planePointLocal = selectionBoxDisplayToRotatingLocal(planePointDisplay, timeMyr);
        selectionBoxPointerState = {
          pointerId: event.pointerId,
          mode: String(hitInfo.handle.kind || "drag"),
          timeMyr,
          planeZ: sceneCenter.z,
          startClientX: Number(event.clientX) || 0.0,
          startClientY: Number(event.clientY) || 0.0,
          startCenter: {
            x: Number(selectionBoxState.center.x) || 0.0,
            y: Number(selectionBoxState.center.y) || 0.0,
            z: Number(selectionBoxState.center.z) || 0.0,
          },
          startHalfWidthPc: clampSelectionBoxHalfWidth(selectionBoxState.halfWidthPc),
          dragOffsetLocal: {
            x: planePointLocal.x - (Number(selectionBoxState.center.x) || 0.0),
            y: planePointLocal.y - (Number(selectionBoxState.center.y) || 0.0),
          },
        };
        controls.enabled = false;
        document.body.style.userSelect = "none";
        if (typeof canvas.setPointerCapture === "function" && event.pointerId !== undefined) {
          try {
            canvas.setPointerCapture(event.pointerId);
          } catch (_err) {
          }
        }
        suppressNextCanvasClick = true;
        event.preventDefault();
        return true;
      }

      function updateSelectionBoxInteraction(event) {
        if (!selectionBoxPointerState) {
          return false;
        }
        const planePointScene = selectionBoxPlanePointFromEvent(event, selectionBoxPointerState.planeZ);
        if (!planePointScene) {
          return true;
        }
        const planePointDisplay = planePointScene.clone().sub(plotGroup.position);
        const planePointLocal = selectionBoxDisplayToRotatingLocal(planePointDisplay, selectionBoxPointerState.timeMyr);
        if (selectionBoxPointerState.mode === "drag") {
          selectionBoxState.center = {
            x: planePointLocal.x - selectionBoxPointerState.dragOffsetLocal.x,
            y: planePointLocal.y - selectionBoxPointerState.dragOffsetLocal.y,
            z: selectionBoxPointerState.startCenter.z,
          };
        } else {
          const dx = Math.abs(planePointLocal.x - selectionBoxPointerState.startCenter.x);
          const dy = Math.abs(planePointLocal.y - selectionBoxPointerState.startCenter.y);
          const resizePcPerPixel = Math.max(Number(selectionBoxSpec.resize_pc_per_pixel) || 4.0, 0.1);
          const mouseDeltaPc = resizePcPerPixel * Math.max(
            Math.abs(Number(event.clientX) - Number(selectionBoxPointerState.startClientX || event.clientX)),
            Math.abs(Number(event.clientY) - Number(selectionBoxPointerState.startClientY || event.clientY))
          );
          const proposedHalfWidth = Math.max(dx, dy, 0.25 * mouseDeltaPc);
          selectionBoxState.halfWidthPc = clampSelectionBoxHalfWidth(proposedHalfWidth);
        }
        renderFrame(currentFrameIndex);
        scheduleSelectionBoxMetricsRecompute(false);
        event.preventDefault();
        return true;
      }

      function finishSelectionBoxInteraction(event) {
        if (!selectionBoxPointerState) {
          return false;
        }
        if (typeof canvas.releasePointerCapture === "function" && selectionBoxPointerState.pointerId !== undefined) {
          try {
            canvas.releasePointerCapture(selectionBoxPointerState.pointerId);
          } catch (_err) {
          }
        }
        selectionBoxPointerState = null;
        controls.enabled = true;
        document.body.style.userSelect = "";
        scheduleSelectionBoxMetricsRecompute(true);
        if (event) {
          event.preventDefault();
        }
        return true;
      }

      function selectionBoxSubsampleWindow(relX, relY, relZ, startIndex, endIndex, maxPoints, seedOffset = 0) {
        const count = Math.max(0, endIndex - startIndex);
        const take = Math.min(count, Math.max(1, Math.round(Number(maxPoints) || 1)));
        const outX = new Float32Array(take);
        const outY = new Float32Array(take);
        const outZ = new Float32Array(take);
        if (!take) {
          return { x: outX, y: outY, z: outZ };
        }
        if (count <= take) {
          for (let index = 0; index < count; index += 1) {
            outX[index] = relX[startIndex + index];
            outY[index] = relY[startIndex + index];
            outZ[index] = relZ[startIndex + index];
          }
          return { x: outX, y: outY, z: outZ };
        }
        const stride = count / take;
        const jitter = (0.5 + 0.5 * Math.sin(12.9898 * (seedOffset + 1))) * stride;
        let lastIndex = startIndex - 1;
        for (let index = 0; index < take; index += 1) {
          let sourceIndex = startIndex + Math.floor(jitter + index * stride);
          if (sourceIndex <= lastIndex) {
            sourceIndex = lastIndex + 1;
          }
          if (sourceIndex >= endIndex) {
            sourceIndex = endIndex - 1;
          }
          lastIndex = sourceIndex;
          outX[index] = relX[sourceIndex];
          outY[index] = relY[sourceIndex];
          outZ[index] = relZ[sourceIndex];
        }
        return { x: outX, y: outY, z: outZ };
      }

      function selectionBoxScaledRandomCatalog(halfWidthPc) {
        const scale = Math.max(Number(halfWidthPc) || 0.0, 0.0);
        const x = new Float32Array(selectionBoxRandomUnitCatalog.x.length);
        const y = new Float32Array(selectionBoxRandomUnitCatalog.y.length);
        const z = new Float32Array(selectionBoxRandomUnitCatalog.z.length);
        for (let index = 0; index < x.length; index += 1) {
          x[index] = selectionBoxRandomUnitCatalog.x[index] * scale;
          y[index] = selectionBoxRandomUnitCatalog.y[index] * scale;
          z[index] = selectionBoxRandomUnitCatalog.z[index] * scale;
        }
        return { x, y, z };
      }

      function selectionBoxCountPairsInBands(x, y, z, bandSq) {
        const counts = new Float32Array(bandSq.length);
        for (let i = 0; i < x.length - 1; i += 1) {
          const xi = x[i];
          const yi = y[i];
          const zi = z[i];
          for (let j = i + 1; j < x.length; j += 1) {
            const dx = xi - x[j];
            const dy = yi - y[j];
            const dz = zi - z[j];
            const distanceSq = dx * dx + dy * dy + dz * dz;
            for (let bandIndex = 0; bandIndex < bandSq.length; bandIndex += 1) {
              if (distanceSq >= bandSq[bandIndex][0] && distanceSq < bandSq[bandIndex][1]) {
                counts[bandIndex] += 1.0;
              }
            }
          }
        }
        return counts;
      }

      function selectionBoxCountCrossPairsInBands(dataX, dataY, dataZ, randX, randY, randZ, bandSq) {
        const counts = new Float32Array(bandSq.length);
        for (let i = 0; i < dataX.length; i += 1) {
          const xi = dataX[i];
          const yi = dataY[i];
          const zi = dataZ[i];
          for (let j = 0; j < randX.length; j += 1) {
            const dx = xi - randX[j];
            const dy = yi - randY[j];
            const dz = zi - randZ[j];
            const distanceSq = dx * dx + dy * dy + dz * dz;
            for (let bandIndex = 0; bandIndex < bandSq.length; bandIndex += 1) {
              if (distanceSq >= bandSq[bandIndex][0] && distanceSq < bandSq[bandIndex][1]) {
                counts[bandIndex] += 1.0;
              }
            }
          }
        }
        return counts;
      }

      function percentileSorted(values, fraction) {
        if (!values.length) {
          return NaN;
        }
        const clampedFraction = clampRange(Number(fraction) || 0.0, 0.0, 1.0);
        const position = clampedFraction * (values.length - 1);
        const lowIndex = Math.floor(position);
        const highIndex = Math.ceil(position);
        if (lowIndex === highIndex) {
          return values[lowIndex];
        }
        const weight = position - lowIndex;
        return values[lowIndex] + (values[highIndex] - values[lowIndex]) * weight;
      }

      function selectionBoxAxisBoundsWithFloor(values, startIndex, endIndex, minSpanPc) {
        let lo = Infinity;
        let hi = -Infinity;
        for (let index = startIndex; index < endIndex; index += 1) {
          const value = Number(values[index]);
          if (!Number.isFinite(value)) {
            continue;
          }
          if (value < lo) {
            lo = value;
          }
          if (value > hi) {
            hi = value;
          }
        }
        if (!Number.isFinite(lo) || !Number.isFinite(hi)) {
          const half = 0.5 * Math.max(Number(minSpanPc) || 25.0, 1.0);
          return [-half, half];
        }
        const span = hi - lo;
        if (!Number.isFinite(span) || span < minSpanPc) {
          const center = 0.5 * (lo + hi);
          const half = 0.5 * Math.max(Number(minSpanPc) || 25.0, 1.0);
          return [center - half, center + half];
        }
        return [lo, hi];
      }

      function selectionBoxMedianNearestNeighborDistance(x, y, z, startIndex = 0, endIndex = x.length) {
        const nPoints = Math.max(0, endIndex - startIndex);
        if (nPoints < 2) {
          return NaN;
        }
        const minDistances = new Array(nPoints).fill(Infinity);
        for (let i = startIndex; i < endIndex - 1; i += 1) {
          const xi = Number(x[i]);
          const yi = Number(y[i]);
          const zi = Number(z[i]);
          for (let j = i + 1; j < endIndex; j += 1) {
            const dx = xi - Number(x[j]);
            const dy = yi - Number(y[j]);
            const dz = zi - Number(z[j]);
            const distance = Math.hypot(dx, dy, dz);
            const relI = i - startIndex;
            const relJ = j - startIndex;
            if (distance < minDistances[relI]) {
              minDistances[relI] = distance;
            }
            if (distance < minDistances[relJ]) {
              minDistances[relJ] = distance;
            }
          }
        }
        const finite = minDistances.filter((value) => Number.isFinite(value));
        if (!finite.length) {
          return NaN;
        }
        finite.sort((a, b) => a - b);
        return percentileSorted(finite, 0.5);
      }

      function createSelectionBoxMetricRng(seed) {
        let state = (Math.round(Number(seed) || 1) >>> 0) || 1;
        return function nextUnit() {
          state = (state + 0x6D2B79F5) >>> 0;
          let t = state;
          t = Math.imul(t ^ (t >>> 15), t | 1);
          t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
          return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
        };
      }

      function selectionBoxRandomMedianNearestNeighborDistance(
        nEvents,
        xRangePc,
        yRangePc,
        zRangePc,
        nRealizations,
        seedBase
      ) {
        const nPoints = Math.max(0, Math.round(Number(nEvents) || 0));
        const realizations = Math.max(1, Math.round(Number(nRealizations) || 1));
        if (nPoints < 2) {
          return NaN;
        }
        const medians = [];
        for (let realizationIndex = 0; realizationIndex < realizations; realizationIndex += 1) {
          const rng = createSelectionBoxMetricRng(
            (Number(seedBase) || 1) + 7919 * (realizationIndex + 1)
          );
          const randX = new Float32Array(nPoints);
          const randY = new Float32Array(nPoints);
          const randZ = new Float32Array(nPoints);
          for (let pointIndex = 0; pointIndex < nPoints; pointIndex += 1) {
            randX[pointIndex] = Number(xRangePc[0]) + rng() * (Number(xRangePc[1]) - Number(xRangePc[0]));
            randY[pointIndex] = Number(yRangePc[0]) + rng() * (Number(yRangePc[1]) - Number(yRangePc[0]));
            randZ[pointIndex] = Number(zRangePc[0]) + rng() * (Number(zRangePc[1]) - Number(zRangePc[0]));
          }
          const medianDistance = selectionBoxMedianNearestNeighborDistance(randX, randY, randZ);
          if (Number.isFinite(medianDistance)) {
            medians.push(medianDistance);
          }
        }
        if (!medians.length) {
          return NaN;
        }
        medians.sort((a, b) => a - b);
        return percentileSorted(medians, 0.5);
      }

      function computeSelectionBoxMetrics() {
        const halfWidthPc = clampSelectionBoxHalfWidth(selectionBoxState.halfWidthPc);
        const sideKpc = (2.0 * halfWidthPc) / 1000.0;
        const volumeKpc3 = Math.max(Math.pow(sideKpc, 3), 1e-6);
        const centerLocal = {
          x: Number(selectionBoxState.center.x) || 0.0,
          y: Number(selectionBoxState.center.y) || 0.0,
          z: Number(selectionBoxState.center.z) || 0.0,
        };
        const source = selectionBoxMetricEvents;
        const relLookback = [];
        const relX = [];
        const relY = [];
        const relZ = [];
        const relSampleIndex = [];
        for (let index = 0; index < source.count; index += 1) {
          const dx = source.x[index] - centerLocal.x;
          const dy = source.y[index] - centerLocal.y;
          const dz = source.z[index] - centerLocal.z;
          if (
            Math.abs(dx) > halfWidthPc
            || Math.abs(dy) > halfWidthPc
            || Math.abs(dz) > halfWidthPc
          ) {
            continue;
          }
          relLookback.push(source.lookback[index]);
          relX.push(dx);
          relY.push(dy);
          relZ.push(dz);
          relSampleIndex.push(source.sampleIndex[index]);
        }

        const nCenters = selectionBoxMetricTimeCenters.length;
        const sampleCount = Math.max(
          source.sampleValues.length,
          Math.round(Number(selectionBoxMetricsSpec.n_imf_samples) || 0),
          1
        );
        const perSample = Array.from({ length: sampleCount }, () => ({
          lookback: [],
          x: [],
          y: [],
          z: [],
        }));
        for (let index = 0; index < relLookback.length; index += 1) {
          const sampleIndex = Number(relSampleIndex[index]);
          if (!(sampleIndex >= 0 && sampleIndex < sampleCount)) {
            continue;
          }
          perSample[sampleIndex].lookback.push(Number(relLookback[index]));
          perSample[sampleIndex].x.push(Number(relX[index]));
          perSample[sampleIndex].y.push(Number(relY[index]));
          perSample[sampleIndex].z.push(Number(relZ[index]));
        }

        const rateCountsBySample = Array.from({ length: sampleCount }, () => new Float32Array(nCenters));
        const windowMyr = Math.max(Number(selectionBoxMetricsSpec.window_myr) || 0.0, 1.0);
        const lookbackMax = Math.max(Number(selectionBoxMetricsSpec.lookback_max_myr) || 0.0, 0.0);
        const halfWindowMyr = 0.5 * windowMyr;
        const rateMean = new Float32Array(nCenters);
        const rateLo = new Float32Array(nCenters);
        const rateHi = new Float32Array(nCenters);
        const clusteringMedian = new Float32Array(nCenters).fill(NaN);
        const clusteringLo = new Float32Array(nCenters).fill(NaN);
        const clusteringHi = new Float32Array(nCenters).fill(NaN);
        const rateScratch = new Array(sampleCount).fill(0.0);
        const clusteringScratch = [];
        const leftIndices = new Int32Array(sampleCount);
        const rightIndices = new Int32Array(sampleCount);
        const randomRealizations = Math.max(
          1,
          Math.round(Number(selectionBoxMetricsSpec.random_realizations) || 16)
        );
        const minAxisSpanPc = Math.max(Number(selectionBoxMetricsSpec.min_axis_span_pc) || 25.0, 1.0);

        for (let centerIndex = 0; centerIndex < nCenters; centerIndex += 1) {
          const centerTime = Number(selectionBoxMetricTimeCenters[centerIndex]);
          const minTime = Math.max(0.0, centerTime - halfWindowMyr);
          const maxTime = Math.min(lookbackMax, centerTime + halfWindowMyr);
          const effectiveWindowMyr = Math.max(maxTime - minTime, 1e-6);
          let rateSum = 0.0;
          clusteringScratch.length = 0;

          for (let sampleIndex = 0; sampleIndex < sampleCount; sampleIndex += 1) {
            const sample = perSample[sampleIndex];
            const lookback = sample.lookback;
            let leftIndex = leftIndices[sampleIndex];
            let rightIndex = rightIndices[sampleIndex];
            while (leftIndex < lookback.length && Number(lookback[leftIndex]) < minTime) {
              leftIndex += 1;
            }
            if (rightIndex < leftIndex) {
              rightIndex = leftIndex;
            }
            while (rightIndex < lookback.length && Number(lookback[rightIndex]) < maxTime) {
              rightIndex += 1;
            }
            leftIndices[sampleIndex] = leftIndex;
            rightIndices[sampleIndex] = rightIndex;

            const nEvents = Math.max(0, rightIndex - leftIndex);
            const rateValue = nEvents / (effectiveWindowMyr * volumeKpc3);
            rateCountsBySample[sampleIndex][centerIndex] = rateValue;
            rateScratch[sampleIndex] = rateValue;
            rateSum += rateValue;

            if (nEvents < 2) {
              continue;
            }
            const dataMedian = selectionBoxMedianNearestNeighborDistance(
              sample.x,
              sample.y,
              sample.z,
              leftIndex,
              rightIndex
            );
            if (!(Number.isFinite(dataMedian) && dataMedian > 0.0)) {
              continue;
            }
            const xRange = selectionBoxAxisBoundsWithFloor(sample.x, leftIndex, rightIndex, minAxisSpanPc);
            const yRange = selectionBoxAxisBoundsWithFloor(sample.y, leftIndex, rightIndex, minAxisSpanPc);
            const zRange = selectionBoxAxisBoundsWithFloor(sample.z, leftIndex, rightIndex, minAxisSpanPc);
            const randMedian = selectionBoxRandomMedianNearestNeighborDistance(
              nEvents,
              xRange,
              yRange,
              zRange,
              randomRealizations,
              734921 + 104729 * (sampleIndex + 1) + 7919 * (centerIndex + 1)
            );
            if (Number.isFinite(randMedian) && randMedian > 0.0) {
              clusteringScratch.push(randMedian / dataMedian);
            }
          }

          rateScratch.sort((a, b) => a - b);
          rateMean[centerIndex] = rateSum / sampleCount;
          rateLo[centerIndex] = percentileSorted(rateScratch, 0.16);
          rateHi[centerIndex] = percentileSorted(rateScratch, 0.84);

          if (clusteringScratch.length) {
            clusteringScratch.sort((a, b) => a - b);
            clusteringMedian[centerIndex] = percentileSorted(clusteringScratch, 0.5);
            clusteringLo[centerIndex] = percentileSorted(clusteringScratch, 0.16);
            clusteringHi[centerIndex] = percentileSorted(clusteringScratch, 0.84);
          }
        }

        return {
          centerLocal,
          halfWidthPc,
          volumeKpc3,
          selectedEventCount: relLookback.length,
          rate: {
            centers: selectionBoxMetricTimeCenters.slice(),
            mean: Array.from(rateMean),
            lo: Array.from(rateLo),
            hi: Array.from(rateHi),
          },
          clustering: {
            centers: selectionBoxMetricTimeCenters.slice(),
            median: Array.from(clusteringMedian),
            lo: Array.from(clusteringLo),
            hi: Array.from(clusteringHi),
          },
        };
      }

      function scheduleSelectionBoxMetricsRecompute(immediate = false) {
        if (!selectionBoxMetricsSpec.enabled) {
          return;
        }
        selectionBoxMetricsVersion += 1;
        const currentVersion = selectionBoxMetricsVersion;
        selectionBoxMetricsPending = true;
        if (selectionBoxMetricsTimer) {
          window.clearTimeout(selectionBoxMetricsTimer);
        }
        const delayMs = immediate ? 0 : 90;
        selectionBoxMetricsTimer = window.setTimeout(() => {
          selectionBoxMetricsTimer = null;
          const nextMetrics = computeSelectionBoxMetrics();
          if (currentVersion !== selectionBoxMetricsVersion) {
            return;
          }
          selectionBoxMetricsCache = nextMetrics;
          selectionBoxMetricsPending = false;
          renderBoxMetricsWidget();
        }, delayMs);
        renderBoxMetricsWidget();
      }

      function restoreInitialLassoSelectionMask() {
        if (minimalModeEnabled) {
          disposeLassoSelectionMask(currentLassoSelectionMask);
          currentLassoSelectionMask = null;
          return Promise.resolve();
        }
        const savedMask = initialState && typeof initialState === "object"
          ? initialState.lasso_selection_mask
          : null;
        if (
          !savedMask
          || typeof savedMask !== "object"
          || !savedMask.data_url
          || !Array.isArray(savedMask.view_projection_matrix)
          || savedMask.view_projection_matrix.length !== 16
        ) {
          disposeLassoSelectionMask(currentLassoSelectionMask);
          currentLassoSelectionMask = null;
          return Promise.resolve();
        }

        return new Promise((resolve) => {
          const image = new Image();
          image.onload = () => {
            const maskCanvas = document.createElement("canvas");
            const width = Math.max(1, Number(image.naturalWidth || image.width || 0));
            const height = Math.max(1, Number(image.naturalHeight || image.height || 0));
            maskCanvas.width = width;
            maskCanvas.height = height;
            const maskCtx = maskCanvas.getContext("2d");
            if (maskCtx) {
              maskCtx.drawImage(image, 0, 0, width, height);
            }
            const texture = new THREE.CanvasTexture(maskCanvas);
            texture.minFilter = THREE.NearestFilter;
            texture.magFilter = THREE.NearestFilter;
            texture.wrapS = THREE.ClampToEdgeWrapping;
            texture.wrapT = THREE.ClampToEdgeWrapping;
            texture.flipY = false;
            texture.generateMipmaps = false;
            texture.needsUpdate = true;
            const maskImageData = maskCtx ? maskCtx.getImageData(0, 0, width, height) : null;
            const viewProjectionMatrix = new THREE.Matrix4().fromArray(
              savedMask.view_projection_matrix.map((value) => Number(value))
            );
            const polygonNdc = Array.isArray(savedMask.polygon_ndc)
              ? savedMask.polygon_ndc
                .map((point) => ({
                  x: Number(point && point.x),
                  y: Number(point && point.y),
                }))
                .filter((point) => Number.isFinite(point.x) && Number.isFinite(point.y))
              : [];
            disposeLassoSelectionMask(currentLassoSelectionMask);
            currentLassoSelectionMask = {
              maskTexture: texture,
              viewProjectionMatrix,
              polygonNdc,
              maskSize: width,
              maskAlphaData: maskImageData ? maskImageData.data : null,
            };
            resolve();
          };
          image.onerror = () => {
            disposeLassoSelectionMask(currentLassoSelectionMask);
            currentLassoSelectionMask = null;
            resolve();
          };
          image.src = String(savedMask.data_url);
        });
      }

      function applyInitialStateSync() {
        resetLegendState(currentGroup);
        if (!initialState || typeof initialState !== "object") {
          if (minimalModeEnabled) {
            clickSelectionEnabled = false;
            lassoVolumeSelectionEnabled = false;
            lassoArmed = false;
          currentSelection = null;
          currentSelections = [];
          selectedClusterKeys = new Set();
          focusSelectionKey = "";
          legendPanelRectState = null;
          legendPanelUserSized = false;
          ["sky", "box_metrics", "age_kde", "cluster_filter", "dendrogram"].forEach((widgetKey) => {
            setWidgetModeValue(widgetKey, "hidden");
          });
            setToolsDrawerOpen(false);
            setControlsDrawerOpen(false);
            setSkyControlsDrawerOpen(false);
            legendPanelOpen = true;
          }
          return;
        }

        const requestedGroup = String(initialState.current_group || "");
        if (requestedGroup && Object.prototype.hasOwnProperty.call(groupVisibility, requestedGroup)) {
          currentGroup = requestedGroup;
        }
        resetLegendState(currentGroup);

        const requestedFrameIndex = Number(initialState.current_frame_index);
        if (Number.isFinite(requestedFrameIndex)) {
          currentFrameIndex = Math.max(0, Math.min(Math.round(requestedFrameIndex), Math.max(frameSpecs.length - 1, 0)));
        }

        if (typeof initialState.click_selection_enabled === "boolean") {
          clickSelectionEnabled = initialState.click_selection_enabled;
        }
        if (typeof initialState.lasso_volume_selection_enabled === "boolean") {
          lassoVolumeSelectionEnabled = initialState.lasso_volume_selection_enabled;
        }
        if (typeof initialState.lasso_armed === "boolean") {
          lassoArmed = initialState.lasso_armed;
        }
        loadManualLabels(initialState.manual_labels, initialState.active_manual_label_id);

        const savedGlobalControls = initialState.global_controls;
        if (savedGlobalControls && typeof savedGlobalControls === "object") {
          if (typeof savedGlobalControls.theme_key === "string" && savedGlobalControls.theme_key) {
            activeThemeKey = String(savedGlobalControls.theme_key);
          }
          if (Number.isFinite(Number(savedGlobalControls.scroll_speed))) {
            globalScrollSpeed = Number(savedGlobalControls.scroll_speed);
          }
          if (Number.isFinite(Number(savedGlobalControls.point_size_scale))) {
            globalPointSizeScale = Number(savedGlobalControls.point_size_scale);
          }
          if (Number.isFinite(Number(savedGlobalControls.point_opacity_scale))) {
            globalPointOpacityScale = Number(savedGlobalControls.point_opacity_scale);
          }
          if (Number.isFinite(Number(savedGlobalControls.point_glow_strength))) {
            globalPointGlowStrength = Number(savedGlobalControls.point_glow_strength);
          }
          if (typeof savedGlobalControls.size_points_by_stars_enabled === "boolean") {
            sizePointsByStarsEnabled = savedGlobalControls.size_points_by_stars_enabled;
          }
          if (Number.isFinite(Number(savedGlobalControls.fade_in_time_myr))) {
            fadeInTimeMyr = Number(savedGlobalControls.fade_in_time_myr);
          }
          if (typeof savedGlobalControls.fade_in_and_out_enabled === "boolean") {
            fadeInAndOutEnabled = savedGlobalControls.fade_in_and_out_enabled;
          }
          if (typeof savedGlobalControls.focus_trace_key === "string") {
            focusTraceKey = String(savedGlobalControls.focus_trace_key);
          }
          if (typeof savedGlobalControls.focus_selection_key === "string") {
            focusSelectionKey = normalizeMemberKey(savedGlobalControls.focus_selection_key);
          }
          if (Number.isFinite(Number(savedGlobalControls.camera_fov))) {
            camera.fov = Number(savedGlobalControls.camera_fov);
          }
          if (typeof savedGlobalControls.axes_visible === "boolean") {
            axesVisible = savedGlobalControls.axes_visible;
          }
          if (typeof savedGlobalControls.galactic_reference_visible === "boolean") {
            galacticReferenceVisible = savedGlobalControls.galactic_reference_visible;
          }
          if (typeof savedGlobalControls.nearby_region_labels_visible === "boolean") {
            nearbyRegionLabelsVisible = savedGlobalControls.nearby_region_labels_visible;
          }
          if (
            typeof savedGlobalControls.sky_dome_source_key === "string"
            && skyDomeSourceOptionByKey.has(savedGlobalControls.sky_dome_source_key)
          ) {
            activeSkyDomeSourceKey = String(savedGlobalControls.sky_dome_source_key);
          }
          if (typeof savedGlobalControls.sky_dome_enabled === "boolean") {
            skyDomeSpec.enabled = savedGlobalControls.sky_dome_enabled;
          }
          if (typeof savedGlobalControls.sky_dome_force_visible === "boolean") {
            skyDomeForceVisible = savedGlobalControls.sky_dome_force_visible;
          }
          if (Number.isFinite(Number(savedGlobalControls.sky_dome_opacity))) {
            skyDomeSpec.opacity = Math.min(Math.max(Number(savedGlobalControls.sky_dome_opacity), 0.0), 1.0);
          }
          if (Number.isFinite(Number(savedGlobalControls.sky_dome_hips_brightness))) {
            skyDomeSpec.hips_brightness = Math.min(Math.max(Number(savedGlobalControls.sky_dome_hips_brightness), 0.1), 8.0);
          }
          if (Number.isFinite(Number(savedGlobalControls.sky_dome_hips_contrast))) {
            skyDomeSpec.hips_contrast = Math.min(Math.max(Number(savedGlobalControls.sky_dome_hips_contrast), 0.1), 4.0);
          }
          if (Number.isFinite(Number(savedGlobalControls.sky_dome_hips_gamma))) {
            skyDomeSpec.hips_gamma = Math.min(Math.max(Number(savedGlobalControls.sky_dome_hips_gamma), 0.2), 4.0);
          }
          if (typeof savedGlobalControls.camera_auto_orbit_enabled === "boolean") {
            cameraAutoOrbitEnabled = savedGlobalControls.camera_auto_orbit_enabled;
          }
          if (savedGlobalControls.camera_auto_orbit_speed !== undefined) {
            cameraAutoOrbitBaseSpeed = normalizeCameraAutoOrbitSpeed(savedGlobalControls.camera_auto_orbit_speed);
          }
          if (savedGlobalControls.camera_auto_orbit_direction !== undefined) {
            cameraAutoOrbitDirection = normalizeCameraAutoOrbitDirection(savedGlobalControls.camera_auto_orbit_direction);
          }
          if (typeof savedGlobalControls.camera_view_mode === "string" && savedGlobalControls.camera_view_mode) {
            cameraViewMode = String(savedGlobalControls.camera_view_mode);
          }
          if (Number.isFinite(Number(savedGlobalControls.earth_view_focus_distance_pc))) {
            earthViewFocusDistance = Number(savedGlobalControls.earth_view_focus_distance_pc);
          }
        }
        skyDomeViewOpacityScale = cameraViewMode === "earth" ? 1.0 : 0.0;

        const savedScaleBarState = initialState.scale_bar_state;
        if (
          savedScaleBarState
          && typeof savedScaleBarState === "object"
          && Number.isFinite(Number(savedScaleBarState.left))
          && Number.isFinite(Number(savedScaleBarState.top))
        ) {
          scaleBarPosition = {
            left: Number(savedScaleBarState.left),
            top: Number(savedScaleBarState.top),
          };
        }

        const savedLegendPanelState = initialState.legend_panel_state;
        if (
          savedLegendPanelState
          && typeof savedLegendPanelState === "object"
          && Number.isFinite(Number(savedLegendPanelState.left))
          && Number.isFinite(Number(savedLegendPanelState.top))
          && Number.isFinite(Number(savedLegendPanelState.width))
          && Number.isFinite(Number(savedLegendPanelState.height))
        ) {
          legendPanelRectState = {
            left: Number(savedLegendPanelState.left),
            top: Number(savedLegendPanelState.top),
            width: Number(savedLegendPanelState.width),
            height: Number(savedLegendPanelState.height),
          };
        }
        if (typeof initialState.legend_panel_user_sized === "boolean") {
          legendPanelUserSized = initialState.legend_panel_user_sized;
        }

        const savedLegendSectionsOpen = initialState.legend_sections_open;
        if (savedLegendSectionsOpen && typeof savedLegendSectionsOpen === "object") {
          if (typeof savedLegendSectionsOpen.traces === "boolean") {
            legendSectionOpenState.traces = savedLegendSectionsOpen.traces;
          }
          if (typeof savedLegendSectionsOpen.volumes === "boolean") {
            legendSectionOpenState.volumes = savedLegendSectionsOpen.volumes;
          }
        }

        if (initialState.active_volume_key && volumeStateKeySet.has(String(initialState.active_volume_key))) {
          activeVolumeKey = String(initialState.active_volume_key);
        }

        const savedLegendState = initialState.legend_state;
        if (savedLegendState && typeof savedLegendState === "object") {
          legendItems.forEach((item) => {
            const itemKey = String(item.key);
            if (Object.prototype.hasOwnProperty.call(savedLegendState, itemKey)) {
              legendState[itemKey] = Boolean(savedLegendState[itemKey]);
            }
          });
        }

        const savedTraceStyles = initialState.trace_style_state;
        if (savedTraceStyles && typeof savedTraceStyles === "object") {
          Object.entries(savedTraceStyles).forEach(([traceKey, state]) => {
            const target = traceStyleStateByKey[String(traceKey)];
            if (!target || !state || typeof state !== "object") {
              return;
            }
            if (typeof state.color === "string" && state.color) {
              target.color = state.color;
            }
            if (Number.isFinite(Number(state.opacity))) {
              target.opacity = clamp01(state.opacity);
            }
            if (Number.isFinite(Number(state.sizeScale))) {
              target.sizeScale = Math.max(Number(state.sizeScale), 0.05);
            }
          });
        }

        const savedVolumeState = initialState.volume_state_by_key;
        if (savedVolumeState && typeof savedVolumeState === "object") {
          Object.entries(savedVolumeState).forEach(([layerKey, state]) => {
            const target = volumeStateByKey[String(layerKey)];
            if (!target || !state || typeof state !== "object") {
              return;
            }
            if (typeof state.visible === "boolean") {
              target.visible = state.visible;
            }
            if (Number.isFinite(Number(state.vmin))) {
              target.vmin = Number(state.vmin);
            }
            if (Number.isFinite(Number(state.vmax))) {
              target.vmax = Number(state.vmax);
            }
            if (Number.isFinite(Number(state.opacity))) {
              target.opacity = Number(state.opacity);
            }
            if (Number.isFinite(Number(state.steps))) {
              target.steps = Number(state.steps);
            }
            if (Number.isFinite(Number(state.alphaCoef))) {
              target.alphaCoef = Number(state.alphaCoef);
            }
            if (Number.isFinite(Number(state.gradientStep))) {
              target.gradientStep = Number(state.gradientStep);
            }
            if (typeof state.showAllTimes === "boolean") {
              target.showAllTimes = state.showAllTimes;
            }
            if (typeof state.stretch === "string" && state.stretch) {
              target.stretch = normalizeVolumeStretch(state.stretch);
            }
            if (typeof state.colormap === "string" && state.colormap) {
              target.colormap = state.colormap;
            }
          });
        }

        const savedClusterFilterState = initialState.cluster_filter_state;
        if (savedClusterFilterState && typeof savedClusterFilterState === "object") {
          if (typeof savedClusterFilterState.parameter_key === "string" && savedClusterFilterState.parameter_key) {
            clusterFilterParameterKey = String(savedClusterFilterState.parameter_key);
          }
          const savedRanges = savedClusterFilterState.ranges_by_key;
          if (savedRanges && typeof savedRanges === "object") {
            Object.entries(savedRanges).forEach(([parameterKey, rangeState]) => {
              if (!rangeState || typeof rangeState !== "object") {
                return;
              }
              clusterFilterRangeStateByKey[String(parameterKey)] = {
                min: Number(rangeState.min),
                max: Number(rangeState.max),
              };
            });
          }
        }
        const savedDendrogramState = initialState.dendrogram_state;
        if (savedDendrogramState && typeof savedDendrogramState === "object") {
          if (typeof savedDendrogramState.trace_key === "string" && savedDendrogramState.trace_key) {
            dendrogramTraceKey = String(savedDendrogramState.trace_key);
          }
          if (typeof savedDendrogramState.connection_mode === "string" && savedDendrogramState.connection_mode) {
            dendrogramConnectionMode = String(savedDendrogramState.connection_mode);
          }
          if (typeof savedDendrogramState.threshold_mode === "string" && savedDendrogramState.threshold_mode) {
            dendrogramThresholdMode = String(savedDendrogramState.threshold_mode);
          }
          if (Number.isFinite(Number(savedDendrogramState.threshold_pc))) {
            dendrogramThresholdPc = Number(savedDendrogramState.threshold_pc);
          }
          if (Number.isFinite(Number(savedDendrogramState.threshold_age_myr))) {
            dendrogramThresholdAgeMyr = Number(savedDendrogramState.threshold_age_myr);
          }
        }
        const savedSelectionBoxState = initialState.selection_box_state;
        if (savedSelectionBoxState && typeof savedSelectionBoxState === "object") {
          const savedCenter = savedSelectionBoxState.center_local_pc || {};
          if (Number.isFinite(Number(savedCenter.x))) {
            selectionBoxState.center.x = Number(savedCenter.x);
          }
          if (Number.isFinite(Number(savedCenter.y))) {
            selectionBoxState.center.y = Number(savedCenter.y);
          }
          if (Number.isFinite(Number(savedCenter.z))) {
            selectionBoxState.center.z = Number(savedCenter.z);
          }
          if (Number.isFinite(Number(savedSelectionBoxState.half_width_pc))) {
            selectionBoxState.halfWidthPc = clampSelectionBoxHalfWidth(savedSelectionBoxState.half_width_pc);
          }
          if (savedSelectionBoxState.visible !== undefined) {
            selectionBoxState.visible = Boolean(savedSelectionBoxState.visible);
          }
          if (savedSelectionBoxState.twopcf_band_pc !== undefined) {
            selectionBoxState.twopcfBandPc = normalizeSelectionBoxBand(savedSelectionBoxState.twopcf_band_pc);
          }
        }

        currentSelections = uniqueSelections(Array.isArray(initialState.current_selections) ? initialState.current_selections : []);
        currentSelection = initialState.current_selection && typeof initialState.current_selection === "object"
          ? initialState.current_selection
          : null;
        if (Array.isArray(initialState.selected_cluster_keys) && initialState.selected_cluster_keys.length) {
          selectedClusterKeys = new Set(
            initialState.selected_cluster_keys
              .map((value) => normalizeMemberKey(value))
              .filter(Boolean)
          );
        } else {
          selectedClusterKeys = new Set(
            currentSelections
              .map((selection) => normalizedSelectionKeyFor(selection))
              .filter(Boolean)
          );
        }
        currentSelectionMode = currentSelection ? "click" : ((currentSelections.length || Boolean(initialState.lasso_selection_mask)) ? "lasso" : "none");

        const savedWidgets = initialState.widgets;
        if (savedWidgets && typeof savedWidgets === "object") {
          ["sky", "box_metrics", "age_kde", "cluster_filter", "dendrogram"].forEach((widgetKey) => {
            const panelState = savedWidgets[widgetKey];
            const panelEl = widgetPanelForKey(widgetKey);
            if (!panelEl || !panelState || typeof panelState !== "object") {
              return;
            }
            if (panelState.rect && typeof panelState.rect === "object") {
              if (Number.isFinite(Number(panelState.rect.left))) {
                panelEl.dataset.normalLeft = String(panelState.rect.left);
              }
              if (Number.isFinite(Number(panelState.rect.top))) {
                panelEl.dataset.normalTop = String(panelState.rect.top);
              }
              if (Number.isFinite(Number(panelState.rect.width))) {
                panelEl.dataset.normalWidth = String(panelState.rect.width);
              }
              if (Number.isFinite(Number(panelState.rect.height))) {
                panelEl.dataset.normalHeight = String(panelState.rect.height);
              }
            }
            if (typeof panelState.mode === "string") {
              setWidgetModeValue(widgetKey, panelState.mode);
            }
          });
        }

        const savedCamera = initialState.camera;
        if (savedCamera && typeof savedCamera === "object") {
          const target = savedCamera.target || {};
          const position = savedCamera.position || {};
          const up = savedCamera.up || {};
          const viewOffset = savedCamera.view_offset || savedCamera.viewOffset || {};
          if (Number.isFinite(Number(target.x)) && Number.isFinite(Number(target.y)) && Number.isFinite(Number(target.z))) {
            controls.target.set(Number(target.x), Number(target.y), Number(target.z));
          }
          if (Number.isFinite(Number(position.x)) && Number.isFinite(Number(position.y)) && Number.isFinite(Number(position.z))) {
            camera.position.set(Number(position.x), Number(position.y), Number(position.z));
          }
          if (Number.isFinite(Number(up.x)) && Number.isFinite(Number(up.y)) && Number.isFinite(Number(up.z))) {
            camera.up.set(Number(up.x), Number(up.y), Number(up.z));
          }
          if (typeof applyActionCameraViewOffset === "function") {
            applyActionCameraViewOffset(viewOffset);
          }
          camera.lookAt(controls.target);
          controls.update();
          initialCameraState.position.copy(camera.position);
          initialCameraState.target.copy(controls.target);
          initialCameraState.up.copy(camera.up);
          initialCameraState.fov = Number(camera.fov);
          initialCameraState.viewOffset = typeof normalizeActionViewOffset === "function"
            ? normalizeActionViewOffset(viewOffset)
            : {
              x: Number.isFinite(Number(viewOffset.x)) ? Number(viewOffset.x) : 0.0,
              y: Number.isFinite(Number(viewOffset.y)) ? Number(viewOffset.y) : 0.0,
            };
        }

        const savedDrawers = initialState.drawers;
        if (savedDrawers && typeof savedDrawers === "object") {
          if (typeof savedDrawers.tools_open === "boolean") {
            setToolsDrawerOpen(savedDrawers.tools_open);
          }
          if (typeof savedDrawers.controls_open === "boolean") {
            setControlsDrawerOpen(savedDrawers.controls_open);
          }
          if (typeof savedDrawers.sky_open === "boolean") {
            setSkyControlsDrawerOpen(savedDrawers.sky_open);
          }
        }

        if (typeof initialState.legend_open === "boolean") {
          legendPanelOpen = initialState.legend_open;
        }

        if (minimalModeEnabled) {
          clickSelectionEnabled = false;
          lassoVolumeSelectionEnabled = false;
          lassoArmed = false;
          currentSelection = null;
          currentSelections = [];
          selectedClusterKeys = new Set();
          focusSelectionKey = "";
          activeLegendEditorKey = "";
          legendPanelRectState = null;
          legendPanelUserSized = false;
          ["sky", "box_metrics", "age_kde", "cluster_filter", "dendrogram"].forEach((widgetKey) => {
            setWidgetModeValue(widgetKey, "hidden");
          });
          setToolsDrawerOpen(false);
          setControlsDrawerOpen(false);
          setSkyControlsDrawerOpen(false);
          legendPanelOpen = true;
        }

        applyGlobalControlState();
        clampClusterFilterRangeForParameter(activeClusterFilterParameterSpec());
        pruneSelectionsToActiveClusterFilter();
        applyCameraViewMode();
        applyThemePreset(activeThemeKey, { rerender: false, syncInput: false });
        renderSceneControls();
        setLegendPanelOpen(legendPanelOpen);
        setZenMode(zenModeEnabled);
      }

      function captureWidgetState(widgetKey) {
        const panelEl = widgetPanelForKey(widgetKey);
        if (!panelEl) {
          return null;
        }
        if (widgetModeForKey(widgetKey) === "normal") {
          storeWidgetRect(widgetKey);
        }
        const left = Number(panelEl.dataset.normalLeft);
        const top = Number(panelEl.dataset.normalTop);
        const width = Number(panelEl.dataset.normalWidth);
        const height = Number(panelEl.dataset.normalHeight);
        return {
          mode: widgetModeForKey(widgetKey),
          rect: {
            left: Number.isFinite(left) ? left : widgetDefaultRect(widgetKey).left,
            top: Number.isFinite(top) ? top : widgetDefaultRect(widgetKey).top,
            width: Number.isFinite(width) ? width : widgetDefaultRect(widgetKey).width,
            height: Number.isFinite(height) ? height : widgetDefaultRect(widgetKey).height,
          },
        };
      }

      function captureRuntimeState() {
        return safeJsonClone({
          current_group: currentGroup,
          current_frame_index: currentFrameIndex,
          manual_labels: manualLabels,
          active_manual_label_id: activeManualLabelId,
          legend_state: legendState,
          trace_style_state: traceStyleStateByKey,
          click_selection_enabled: clickSelectionEnabled,
          lasso_volume_selection_enabled: lassoVolumeSelectionEnabled,
          lasso_armed: lassoArmed,
          current_selection: currentSelection,
          current_selections: currentSelections,
          current_selection_mode: currentSelectionMode,
          selected_cluster_keys: Array.from(selectedClusterKeys),
          active_volume_key: activeVolumeKey,
          volume_state_by_key: volumeStateByKey,
          global_controls: {
            theme_key: activeThemeKey,
            scroll_speed: globalScrollSpeed,
            camera_fov: camera.fov,
            point_size_scale: globalPointSizeScale,
            point_opacity_scale: globalPointOpacityScale,
            point_glow_strength: globalPointGlowStrength,
            size_points_by_stars_enabled: sizePointsByStarsEnabled,
            fade_in_time_myr: fadeInTimeMyr,
            fade_in_and_out_enabled: fadeInAndOutEnabled,
            focus_trace_key: focusTraceKey,
            focus_selection_key: focusSelectionKey,
            axes_visible: axesVisible,
            galactic_reference_visible: galacticReferenceVisible,
            nearby_region_labels_visible: nearbyRegionLabelsVisible,
            sky_dome_source_key: activeSkyDomeSourceKey || defaultSkyDomeSourceKey(),
            sky_dome_enabled: Boolean(skyDomeSpec.enabled),
            sky_dome_force_visible: Boolean(skyDomeForceVisible),
            sky_dome_opacity: Number(skyDomeSpec.opacity ?? 0.55),
            sky_dome_hips_brightness: Number(skyDomeSpec.hips_brightness ?? 2.4),
            sky_dome_hips_contrast: Number(skyDomeSpec.hips_contrast ?? 1.25),
            sky_dome_hips_gamma: Number(skyDomeSpec.hips_gamma ?? 1.35),
            camera_auto_orbit_enabled: cameraAutoOrbitEnabled,
            camera_auto_orbit_speed: cameraAutoOrbitBaseSpeed,
            camera_auto_orbit_direction: cameraAutoOrbitDirection,
            camera_view_mode: cameraViewMode,
            earth_view_focus_distance_pc: earthViewFocusDistance,
          },
          camera: {
            position: { x: camera.position.x, y: camera.position.y, z: camera.position.z },
            target: { x: controls.target.x, y: controls.target.y, z: controls.target.z },
            up: { x: camera.up.x, y: camera.up.y, z: camera.up.z },
            view_offset: typeof normalizeActionViewOffset === "function"
              ? normalizeActionViewOffset(currentActionCameraViewOffset)
              : { x: 0.0, y: 0.0 },
          },
          drawers: {
            tools_open: Boolean(toolsShellEl && toolsShellEl.dataset.open === "true"),
            controls_open: Boolean(controlsShellEl && controlsShellEl.dataset.open === "true"),
            sky_open: Boolean(skyControlsShellEl && skyControlsShellEl.dataset.open === "true"),
          },
          legend_open: legendPanelOpen,
          legend_panel_state: captureLegendPanelState(),
          legend_panel_user_sized: legendPanelUserSized,
          legend_sections_open: safeJsonClone(legendSectionOpenState, { traces: true, volumes: true }),
          sky_layers: serializableSkyLayers(),
          active_sky_layer_key: activeSkyLayerKey || "",
          zen_mode_enabled: zenModeEnabled,
          widgets: {
            sky: captureWidgetState("sky"),
            box_metrics: captureWidgetState("box_metrics"),
            age_kde: captureWidgetState("age_kde"),
            cluster_filter: captureWidgetState("cluster_filter"),
            dendrogram: captureWidgetState("dendrogram"),
          },
          scale_bar_state: captureScaleBarState(),
          selection_box_state: {
            center_local_pc: {
              x: Number(selectionBoxState.center.x) || 0.0,
              y: Number(selectionBoxState.center.y) || 0.0,
              z: Number(selectionBoxState.center.z) || 0.0,
            },
            visible: selectionBoxState.visible !== false,
            half_width_pc: clampSelectionBoxHalfWidth(selectionBoxState.halfWidthPc),
            twopcf_band_pc: activeSelectionBoxBandPc(),
          },
          cluster_filter_state: {
            parameter_key: clusterFilterParameterKey,
            ranges_by_key: clusterFilterRangeStateByKey,
          },
          dendrogram_state: {
            trace_key: dendrogramTraceKey,
            connection_mode: dendrogramConnectionMode,
            threshold_mode: dendrogramThresholdMode,
            threshold_pc: dendrogramThresholdPc,
            threshold_age_myr: dendrogramThresholdAgeMyr,
          },
          lasso_selection_mask: captureLassoSelectionMaskState(currentLassoSelectionMask),
        }, {});
      }

      function slugifyFilename(text) {
        const cleaned = String(text || "oviz-threejs-state")
          .trim()
          .replace(/[^a-zA-Z0-9._-]+/g, "_")
          .replace(/^_+|_+$/g, "");
        return cleaned || "oviz-threejs-state";
      }

      function defaultExportFilename() {
        const title = slugifyFilename(sceneSpec.title || "oviz-threejs-state");
        const frame = currentFrame();
        const frameLabel = frame ? slugifyFilename(`t${frame.name}`) : "t0";
        return `${title}-${frameLabel}.html`;
      }

      async function buildExportHtml(exportSceneSpec) {
        const currentHtml = removeExistingSceneSpecPayloadHtml(
          "<!DOCTYPE html>\\n" + document.documentElement.outerHTML
        );
        const startMarker = "/*__SCENE_SPEC_START__*/";
        const endMarker = "/*__SCENE_SPEC_END__*/";
        const startIndex = currentHtml.indexOf(startMarker);
        const endIndex = currentHtml.indexOf(endMarker);
        if (startIndex < 0 || endIndex < 0 || endIndex <= startIndex) {
          throw new Error("Could not locate scene export markers in the current figure HTML.");
        }
        const payloadSceneSpec = sceneSpecForRawExport(exportSceneSpec);
        const sceneJsonText = JSON.stringify(payloadSceneSpec);
        let payloadMetadata = sceneSpecPayloadMetadataForJsonText(sceneJsonText);
        let sceneSpecMarkerText = sceneSpecMarkerWrappedJsonText(sceneJsonText);
        let payloadHtml = "";
        if (shouldCompressExportSceneSpec(sceneJsonText)) {
          try {
            const compressedPayload = await gzipBase64EncodeText(sceneJsonText);
            if (compressedPayload) {
              const payloadId = sceneSpecPayloadElementId();
              payloadMetadata = sceneSpecPayloadMetadataForCompressedJsonText(
                sceneJsonText,
                compressedPayload.compressedSizeBytes,
                compressedPayload.embeddedBase64SizeBytes
              );
              sceneSpecMarkerText = sceneSpecMarkerWrappedCompressedPayload(payloadId);
              payloadHtml = sceneSpecPayloadScriptHtml(payloadId, compressedPayload.encoded);
            } else {
              console.warn(
                "oviz Save State could not compress the scene payload because this browser does not support CompressionStream('gzip'). Exporting raw scene JSON."
              );
            }
          } catch (err) {
            console.warn(
              "oviz Save State could not compress the scene payload. Exporting raw scene JSON.",
              err
            );
          }
        }
        let exportHtml = (
          currentHtml.slice(0, startIndex)
          + sceneSpecMarkerText
          + currentHtml.slice(endIndex + endMarker.length)
        );
        exportHtml = insertSceneSpecPayloadHtml(exportHtml, payloadHtml);
        const metadataStartMarker = "/*__SCENE_SPEC_METADATA_START__*/";
        const metadataEndMarker = "/*__SCENE_SPEC_METADATA_END__*/";
        const metadataStartIndex = exportHtml.indexOf(metadataStartMarker);
        const metadataEndIndex = exportHtml.indexOf(metadataEndMarker);
        if (
          metadataStartIndex >= 0
          && metadataEndIndex > metadataStartIndex
        ) {
          exportHtml = (
            exportHtml.slice(0, metadataStartIndex)
            + sceneSpecMetadataMarkerWrappedJson(payloadMetadata)
            + exportHtml.slice(metadataEndIndex + metadataEndMarker.length)
          );
        }
        return exportHtml;
      }

      async function saveSceneStateToHtml() {
        const exportSceneSpec = safeJsonClone(sceneSpec, {});
        exportSceneSpec.initial_state = captureRuntimeState();
        exportSceneSpec.width = Math.max(root.clientWidth || sceneSpec.width || 900, 1);
        exportSceneSpec.height = Math.max(root.clientHeight || sceneSpec.height || 700, 1);
        const htmlText = await buildExportHtml(exportSceneSpec);
        const suggestedName = defaultExportFilename();

        if (typeof window.showSaveFilePicker === "function") {
          try {
            const handle = await window.showSaveFilePicker({
              suggestedName,
              types: [{
                description: "HTML file",
                accept: { "text/html": [".html"] },
              }],
            });
            const writable = await handle.createWritable();
            await writable.write(htmlText);
            await writable.close();
            return;
          } catch (err) {
            if (err && err.name === "AbortError") {
              return;
            }
          }
        }

        const requestedName = window.prompt("Save figure as", suggestedName);
        if (!requestedName) {
          return;
        }
        const filename = String(requestedName).toLowerCase().endsWith(".html")
          ? String(requestedName)
          : `${requestedName}.html`;
        const blob = new Blob([htmlText], { type: "text/html;charset=utf-8" });
        const url = URL.createObjectURL(blob);
        const link = document.createElement("a");
        link.href = url;
        link.download = filename;
        document.body.appendChild(link);
        link.click();
        link.remove();
        window.setTimeout(() => URL.revokeObjectURL(url), 1000);
      }

      function selectionKeyFor(selection) {
        if (!selection || typeof selection !== "object") {
          return "";
        }
        return String(selection.cluster_name || selection.trace_name || "").trim();
      }

      function selectionIdentityKeyFor(selection) {
        if (!selection || typeof selection !== "object") {
          return "";
        }
        const clusterName = selection.cluster_name ? String(selection.cluster_name).trim() : "";
        if (clusterName) {
          return clusterName;
        }
        const x0 = Number(selection.x0);
        const y0 = Number(selection.y0);
        const z0 = Number(selection.z0);
        if (Number.isFinite(x0) && Number.isFinite(y0) && Number.isFinite(z0)) {
          return `${x0.toFixed(6)}|${y0.toFixed(6)}|${z0.toFixed(6)}`;
        }
        const ra = Number(selection.ra_deg);
        const dec = Number(selection.dec_deg);
        if (Number.isFinite(ra) && Number.isFinite(dec)) {
          return `${ra.toFixed(6)}|${dec.toFixed(6)}`;
        }
        return String(selection.trace_name || "").trim();
      }

      function normalizedSelectionKeyFor(selection) {
        const key = selectionIdentityKeyFor(selection);
        return key ? normalizeMemberKey(key) : "";
      }

      function uniqueSelections(selections) {
        const items = Array.isArray(selections) ? selections : [selections];
        const unique = [];
        const seen = new Set();
        items.forEach((selection) => {
          const key = normalizedSelectionKeyFor(selection);
          if (!key || seen.has(key)) {
            return;
          }
          seen.add(key);
          unique.push(selection);
        });
        return unique;
      }

      function crossHoverEnabled() {
        return Boolean(skySpec.enabled && currentSelections.length && !currentSelection);
      }

      function dendrogramFocusTraceKey() {
        return dendrogramSpec.enabled && widgetModeForKey("dendrogram") !== "hidden"
          ? String(dendrogramTraceKey || "")
          : "";
      }

      function normalizeSelectionKeySet(values) {
        return new Set(
          (Array.isArray(values) ? values : Array.from(values || []))
            .map((value) => normalizeMemberKey(value))
            .filter(Boolean)
        );
      }

      function selectionKeySetsEqual(a, b) {
        const aSet = a instanceof Set ? a : new Set(a || []);
        const bSet = b instanceof Set ? b : new Set(b || []);
        return aSet.size === bSet.size && Array.from(aSet).every((value) => bSet.has(value));
      }

      function activeDendrogramSelectionKeys() {
        return dendrogramPinnedSelectionKeys.size ? dendrogramPinnedSelectionKeys : dendrogramHoveredSelectionKeys;
      }

      function setDendrogramHoveredSelectionKeys(keys, label = "", count = 0, regionKey = "") {
        const nextKeys = normalizeSelectionKeySet(keys);
        const nextRegionKey = String(regionKey || "");
        if (
          selectionKeySetsEqual(dendrogramHoveredSelectionKeys, nextKeys)
          && dendrogramHoveredBranchLabel === String(label || "")
          && dendrogramHoveredBranchCount === Number(count || 0)
          && dendrogramHoveredRegionKey === nextRegionKey
        ) {
          return;
        }
        dendrogramHoveredSelectionKeys = nextKeys;
        dendrogramHoveredBranchLabel = String(label || "");
        dendrogramHoveredBranchCount = Number(count || 0);
        dendrogramHoveredRegionKey = nextRegionKey;
        applySceneHoverState();
        if (widgetModeForKey("dendrogram") !== "hidden") {
          renderFrame(currentFrameIndex);
        }
        renderDendrogramWidget();
      }

      function setDendrogramPinnedSelectionKeys(keys, label = "", count = 0, regionKey = "") {
        const nextKeys = normalizeSelectionKeySet(keys);
        const nextRegionKey = String(regionKey || "");
        if (
          selectionKeySetsEqual(dendrogramPinnedSelectionKeys, nextKeys)
          && dendrogramPinnedBranchLabel === String(label || "")
          && dendrogramPinnedBranchCount === Number(count || 0)
          && dendrogramPinnedRegionKey === nextRegionKey
        ) {
          return;
        }
        dendrogramPinnedSelectionKeys = nextKeys;
        dendrogramPinnedBranchLabel = String(label || "");
        dendrogramPinnedBranchCount = Number(count || 0);
        dendrogramPinnedRegionKey = nextRegionKey;
        applySceneHoverState();
        if (widgetModeForKey("dendrogram") !== "hidden") {
          renderFrame(currentFrameIndex);
        }
        renderDendrogramWidget();
      }

      function clearDendrogramHoverState() {
        setDendrogramHoveredSelectionKeys([], "", 0, "");
      }

      function clearDendrogramPinnedState() {
        setDendrogramPinnedSelectionKeys([], "", 0, "");
      }

      function clearDendrogramSelectionState() {
        clearDendrogramPinnedState();
        clearDendrogramHoverState();
      }

      function activeHoveredClusterKeys() {
        const keys = new Set();
        [localHoveredClusterKey]
          .map((value) => normalizeMemberKey(value))
          .filter(Boolean)
          .forEach((value) => keys.add(value));
        if (crossHoverEnabled()) {
          [skyHoveredClusterKey]
            .map((value) => normalizeMemberKey(value))
            .filter(Boolean)
            .forEach((value) => keys.add(value));
        }
        if (widgetModeForKey("dendrogram") !== "hidden") {
          activeDendrogramSelectionKeys().forEach((value) => keys.add(normalizeMemberKey(value)));
        }
        return keys;
      }

      function applySceneHoverState() {
        updateCameraResponsivePointSprites();
      }

      function interpolateDendrogramPosition(entry, targetTimeMyr) {
        if (!entry || typeof entry !== "object") {
          return null;
        }
        const times = Array.isArray(entry.time_samples) ? entry.time_samples : [];
        const xs = Array.isArray(entry.x_samples) ? entry.x_samples : [];
        const ys = Array.isArray(entry.y_samples) ? entry.y_samples : [];
        const zs = Array.isArray(entry.z_samples) ? entry.z_samples : [];
        const sampleCount = Math.min(times.length, xs.length, ys.length, zs.length);
        if (!sampleCount) {
          return null;
        }
        const samples = [];
        for (let index = 0; index < sampleCount; index += 1) {
          const time = Number(times[index]);
          const x = Number(xs[index]);
          const y = Number(ys[index]);
          const z = Number(zs[index]);
          if (!Number.isFinite(time) || !Number.isFinite(x) || !Number.isFinite(y) || !Number.isFinite(z)) {
            continue;
          }
          samples.push({ time, x, y, z });
        }
        if (!samples.length) {
          return null;
        }
        samples.sort((a, b) => a.time - b.time);
        const targetTime = Number(targetTimeMyr);
        if (!Number.isFinite(targetTime)) {
          return null;
        }
        if (targetTime <= samples[0].time) {
          return samples[0];
        }
        if (targetTime >= samples[samples.length - 1].time) {
          return samples[samples.length - 1];
        }
        for (let index = 1; index < samples.length; index += 1) {
          const right = samples[index];
          if (targetTime > right.time) {
            continue;
          }
          const left = samples[index - 1];
          const dt = right.time - left.time;
          const fraction = Math.abs(dt) <= 1e-12 ? 0.0 : (targetTime - left.time) / dt;
          return {
            time: targetTime,
            x: left.x + (right.x - left.x) * fraction,
            y: left.y + (right.y - left.y) * fraction,
            z: left.z + (right.z - left.z) * fraction,
          };
        }
        return samples[samples.length - 1];
      }

      function applyLassoSelectionMaskUniforms(uniforms, mask) {
        const activeMask = mask && mask.maskTexture && mask.viewProjectionMatrix ? mask : null;
        uniforms.useSelectionPolygon.value = Boolean(activeMask);
        uniforms.selectionViewProjectionMatrix.value.copy(
          activeMask ? activeMask.viewProjectionMatrix : new THREE.Matrix4()
        );
        uniforms.selectionMaskTexture.value = activeMask ? activeMask.maskTexture : null;
        uniforms.selectionDimOutside.value = activeMask ? 0.0 : 1.0;
      }

      function selectionToolbarText(selections, focusSelection) {
        const activeSelections = uniqueSelections(selections);
        const clickHint = `Click select: ${clickSelectionEnabled ? "on" : "off"}`;
        const focusLabel = selectionKeyFor(focusSelection);
        const volumeHint = lassoVolumeSelectionEnabled
          ? (hasActiveLassoSelectionMask() ? "Volume lasso: on" : "Volume lasso: armed")
          : "Volume lasso: off";
        if (!activeSelections.length && !focusLabel) {
          return `Shift+drag or use Lasso to select clusters on the current frame.\n${volumeHint}\n${clickHint}`;
        }
        if (!activeSelections.length && focusLabel) {
          return `Focused: ${focusLabel}\n${volumeHint}\n${clickHint}`;
        }
        const labels = activeSelections
          .map((selection) => selectionKeyFor(selection))
          .filter(Boolean);
        const preview = labels.slice(0, 4).join(", ");
        const suffix = labels.length > 4 ? ` +${labels.length - 4} more` : "";
        const focusText = focusLabel ? `\nFocused: ${focusLabel}` : "";
        return `${labels.length} selected\n${preview}${suffix}${focusText}\n${volumeHint}\n${clickHint}`;
      }

      function skyReadoutText(selections, catalogPayload, mode = "overview", volumeOverlay = null) {
        const activeSelections = uniqueSelections(selections);
        const clusterCatalogPayload = catalogPayload || [];
        const overlayPixelCount = volumeOverlay && Number.isFinite(Number(volumeOverlay.non_zero_pixels))
          ? Number(volumeOverlay.non_zero_pixels)
          : 0;
        const overlaySampleCount = volumeOverlay && Number.isFinite(Number(volumeOverlay.sample_count))
          ? Number(volumeOverlay.sample_count)
          : 0;
        const overviewText = (
          "View: Mollweide all-sky\\n"
          + "Center: Galactic Center\\n"
          + `Survey: ${String(skySpec.survey || "P/DSS2/color")}`
        );
        if (mode !== "click") {
          if (!activeSelections.length && !overlayPixelCount) {
            return overviewText;
          }
          const labels = activeSelections
            .map((selection) => selectionKeyFor(selection))
            .filter(Boolean);
          const starCount = clusterCatalogPayload.reduce(
            (total, catalog) => total + (((catalog && catalog.points) || []).length || 0),
            0
          );
          const overlayLine = overlayPixelCount
            ? `Sky overlay pixels: ${overlayPixelCount}${overlaySampleCount ? ` from ${overlaySampleCount} samples` : ""}\n`
            : "";
          if (!activeSelections.length) {
            return overlayLine + overviewText;
          }
          return (
            `Clusters selected: ${labels.length}\n`
            + `${labels.slice(0, 6).join(", ")}${labels.length > 6 ? ` +${labels.length - 6} more` : ""}\n`
            + `Stars shown: ${starCount}\n`
            + overlayLine
            + overviewText
          );
        }

        if (!activeSelections.length) {
          return overviewText;
        }

        if (activeSelections.length === 1) {
          const selection = activeSelections[0];
          const clusterLabel = selectionKeyFor(selection);
          const clusterName = clusterLabel ? `Cluster: ${clusterLabel}\n` : "";
          const nStars = ((catalogPayload || [])[0] || {}).points ? (catalogPayload[0].points || []).length : 0;
          const starLine = nStars > 1 ? `Stars shown: ${nStars}\n` : "";
          return (
            clusterName
            + `Selected direction (t=0): l=${Number(selection.l_deg).toFixed(4)} deg, b=${Number(selection.b_deg).toFixed(4)} deg\n`
            + starLine
            + `Beam radius: ${Number(skySpec.radius_deg || 1.0).toFixed(2)} deg\n`
            + `Survey: ${String(skySpec.survey || "P/DSS2/color")}`
          );
        }
        return overviewText;
      }

      function buildEmptySkySrcdoc() {
        return buildAladinSrcdoc([], [], "overview", null);
      }

__SKY_RUNTIME_JS__

      preloadNativeHipsStartupTiles();
      if (typeof skyDomeUsesNativeHips === "function" && skyDomeUsesNativeHips()) {
        initializeSkyDomeFromSceneSpec();
        await waitForNativeHipsStartup();
      }

      function validSkyDomeImageDataUrl(value) {
        return typeof value === "string" && value.startsWith("data:image/");
      }

      function skyDomeHasLocalSources() {
        return skyDomeSourceOptions.length > 0;
      }

      function skyDomeSourceDataUrl(source) {
        if (!source || typeof source !== "object") {
          return "";
        }
        const candidates = [
          source.data_url,
          source.dataUrl,
          source.image_data_url,
          source.imageDataUrl,
          source.image,
          source.src,
          source.url,
        ];
        for (const candidate of candidates) {
          const dataUrl = String(candidate || "");
          if (validSkyDomeImageDataUrl(dataUrl)) {
            return dataUrl;
          }
        }
        return "";
      }

      function skyDomeFaceDataUrl(face) {
        if (!face || typeof face !== "object") {
          return "";
        }
        const candidates = [face.data_url, face.dataUrl, face.image_data_url, face.imageDataUrl, face.src, face.url];
        for (const candidate of candidates) {
          const dataUrl = String(candidate || "");
          if (validSkyDomeImageDataUrl(dataUrl)) {
            return dataUrl;
          }
        }
        return "";
      }

      function skyDomeSourceFaces(source) {
        const faceKeys = ["px", "nx", "py", "ny", "pz", "nz"];
        const rawFaces = source && typeof source === "object" ? source.faces : null;
        if (Array.isArray(rawFaces)) {
          const faces = rawFaces
            .map((face) => ({
              key: String((face && (face.key || face.face || face.name)) || "").toLowerCase(),
              dataUrl: skyDomeFaceDataUrl(face),
            }))
            .filter((face) => faceKeys.includes(face.key) && validSkyDomeImageDataUrl(face.dataUrl));
          return faceKeys.every((key) => faces.some((face) => face.key === key)) ? faces : [];
        }
        if (rawFaces && typeof rawFaces === "object") {
          const faces = faceKeys
            .map((key) => {
              const value = rawFaces[key];
              const face = value && typeof value === "object" ? value : { data_url: value };
              return { key, dataUrl: skyDomeFaceDataUrl(face) };
            })
            .filter((face) => validSkyDomeImageDataUrl(face.dataUrl));
          return faces.length === faceKeys.length ? faces : [];
        }
        return [];
      }

      function skyDomeSourceLabel(source, fallbackKey, index) {
        const candidates = [source.label, source.name, source.title, source.survey, fallbackKey];
        for (const candidate of candidates) {
          const label = String(candidate || "").trim();
          if (label) {
            return label;
          }
        }
        return `Sky image ${index + 1}`;
      }

      function addSkyDomeRawSource(rawSources, source, keyHint) {
        if (!source || typeof source !== "object") {
          return;
        }
        const key = String(source.key || source.id || source.name || keyHint || "").trim();
        rawSources.push(Object.assign({}, source, key ? { key } : {}));
      }

      function rawSkyDomeSources(spec) {
        const rawSources = [];
        if (!spec || typeof spec !== "object") {
          return rawSources;
        }
        if (Array.isArray(spec.sources)) {
          spec.sources.forEach((source, index) => addSkyDomeRawSource(rawSources, source, `source-${index + 1}`));
        } else if (spec.sources && typeof spec.sources === "object") {
          Object.entries(spec.sources).forEach(([key, source]) => {
            addSkyDomeRawSource(rawSources, source, key);
          });
        }
        if (Array.isArray(spec.local_sources)) {
          spec.local_sources.forEach((source, index) => addSkyDomeRawSource(rawSources, source, `local-${index + 1}`));
        }
        if (spec.source && typeof spec.source === "object") {
          addSkyDomeRawSource(rawSources, spec.source, "source");
        }
        if (skyDomeSourceDataUrl(spec)) {
          addSkyDomeRawSource(rawSources, spec, spec.source_key || spec.default_source_key || spec.source || "default");
        }
        return rawSources;
      }

      function normalizeSkyDomeSourceOptions(spec) {
        const options = [];
        const seenKeys = new Set();
        rawSkyDomeSources(spec).forEach((source, index) => {
          const faces = skyDomeSourceFaces(source);
          const dataUrl = skyDomeSourceDataUrl(source);
          if (!dataUrl) {
            return;
          }
          const baseKey = String(source.key || source.id || source.name || source.label || `source-${index + 1}`).trim()
            || `source-${index + 1}`;
          let key = baseKey;
          let suffix = 2;
          while (seenKeys.has(key)) {
            key = `${baseKey}-${suffix}`;
            suffix += 1;
          }
          seenKeys.add(key);
          options.push({
            key,
            label: skyDomeSourceLabel(source, key, index),
            projection: String(source.projection || source.projection_name || spec.projection || "CAR"),
            survey: String(source.survey || source.label || source.name || key),
            dataUrl,
            faces,
            source,
          });
        });
        return options;
      }

      function skyDomePreferredSourceKeyCandidates(includeSavedState = true) {
        const candidates = [];
        const savedGlobalControls = initialState && initialState.global_controls;
        if (
          includeSavedState
          && savedGlobalControls
          && typeof savedGlobalControls === "object"
          && typeof savedGlobalControls.sky_dome_source_key === "string"
        ) {
          candidates.push(savedGlobalControls.sky_dome_source_key);
        }
        [
          skyDomeSpec.source_key,
          skyDomeSpec.selected_source_key,
          skyDomeSpec.default_source_key,
          skyDomeSpec.main_source_key,
          skyDomeSpec.active_source_key,
          skyDomeSpec.chosen_source_key,
          typeof skyDomeSpec.source === "string" ? skyDomeSpec.source : "",
        ].forEach((value) => {
          if (value !== undefined && value !== null) {
            candidates.push(value);
          }
        });
        return candidates.map((value) => String(value || "").trim()).filter(Boolean);
      }

      function defaultSkyDomeSourceKey(options = {}) {
        if (!skyDomeSourceOptions.length) {
          return "";
        }
        for (const key of skyDomePreferredSourceKeyCandidates(options.includeSavedState !== false)) {
          if (skyDomeSourceOptionByKey.has(key)) {
            return key;
          }
        }
        const flagged = skyDomeSourceOptions.find((option) => {
          const source = option.source || {};
          return Boolean(source.selected || source.default || source.main || source.active || source.chosen);
        });
        return flagged ? flagged.key : skyDomeSourceOptions[0].key;
      }

      function syncSkyDomeSourceSelector() {
        if (!skyDomeSourceFieldEl || !skyDomeSourceSelectEl) {
          syncSkyPanelSourceSelector();
          return;
        }
        if (!skyDomeHasLocalSources()) {
          skyDomeSourceFieldEl.hidden = true;
          skyDomeSourceSelectEl.innerHTML = "";
          syncSkyPanelSourceSelector();
          return;
        }
        skyDomeSourceFieldEl.hidden = false;
        if (skyDomeSourceSelectEl.options.length !== skyDomeSourceOptions.length) {
          skyDomeSourceSelectEl.innerHTML = "";
          skyDomeSourceOptions.forEach((source) => {
            const option = document.createElement("option");
            option.value = source.key;
            option.textContent = source.label;
            skyDomeSourceSelectEl.appendChild(option);
          });
        }
        skyDomeSourceSelectEl.value = activeSkyDomeSourceKey || defaultSkyDomeSourceKey();
        syncSkyPanelSourceSelector();
      }

      function setSkyDomeSourceByKey(sourceKey, options = {}) {
        const key = String(sourceKey || defaultSkyDomeSourceKey(options) || "");
        const source = skyDomeSourceOptionByKey.get(key) || skyDomeSourceOptions[0] || null;
        if (!source) {
          skyDomeLocalSourceActive = false;
          syncSkyDomeSourceSelector();
          return false;
        }
        if (!options.force && skyDomeLocalSourceActive && activeSkyDomeSourceKey === source.key && skyDomeMesh) {
          syncSkyDomeSourceSelector();
          return true;
        }
        activeSkyDomeSourceKey = source.key;
        skyDomeLocalSourceActive = true;
        if (root && root.dataset) {
          root.dataset.skyDomeSource = source.key;
        }
        syncSkyDomeSourceSelector();
        setSkyDomeTextureFromDataUrl(source.dataUrl, source.survey, source.projection);
        if (typeof updateSkyPanel === "function") {
          updateSkyPanel();
        }
        return true;
      }

      function applyInitialSkyDomeSource() {
        syncSkyDomeSourceSelector();
        if (skyDomeHasLocalSources()) {
          setSkyDomeSourceByKey(activeSkyDomeSourceKey || defaultSkyDomeSourceKey(), { force: true });
        } else {
          initializeSkyDomeFromSceneSpec();
        }
      }

      function resetSkyDomeSourceSelection() {
        if (skyDomeHasLocalSources()) {
          setSkyDomeSourceByKey(defaultSkyDomeSourceKey({ includeSavedState: false }), { force: true, includeSavedState: false });
        }
      }

      function activeSkyDomeSourceOption() {
        return (
          skyDomeSourceOptionByKey.get(activeSkyDomeSourceKey)
          || skyDomeSourceOptionByKey.get(defaultSkyDomeSourceKey())
          || skyDomeSourceOptions[0]
          || null
        );
      }

      function syncSkyPanelSourceSelector() {
        if (!skyImageSelectEl) {
          return;
        }
        if (!skyDomeHasLocalSources()) {
          skyImageSelectEl.innerHTML = "";
          skyImageSelectEl.disabled = true;
          return;
        }
        skyImageSelectEl.disabled = false;
        if (skyImageSelectEl.options.length !== skyDomeSourceOptions.length) {
          skyImageSelectEl.innerHTML = "";
          skyDomeSourceOptions.forEach((source) => {
            const option = document.createElement("option");
            option.value = source.key;
            option.textContent = source.label;
            skyImageSelectEl.appendChild(option);
          });
        }
        skyImageSelectEl.value = (activeSkyDomeSourceKey || defaultSkyDomeSourceKey());
      }

      function cleanSkyLayerSurvey(value) {
        return String(value || "").trim();
      }

      function normalizeSkyLayerKey(value) {
        return cleanSkyLayerSurvey(value).replace(/\s+/g, " ");
      }

      function readSkyPresetField(value, fieldNames) {
        const fields = Array.isArray(fieldNames) ? fieldNames : [fieldNames];
        for (const field of fields) {
          if (!field) {
            continue;
          }
          if (value && Object.prototype.hasOwnProperty.call(value, field)) {
            const direct = cleanSkyLayerSurvey(value[field]);
            if (direct) {
              return direct;
            }
          }
          const metadata = value && (value.meta || value.metadata || value.properties || value.hips);
          if (metadata && Object.prototype.hasOwnProperty.call(metadata, field)) {
            const nested = cleanSkyLayerSurvey(metadata[field]);
            if (nested) {
              return nested;
            }
          }
        }
        return "";
      }

      function normalizeSkyLayerPreset(rawPreset, index = 0) {
        if (typeof rawPreset === "string") {
          const survey = cleanSkyLayerSurvey(rawPreset);
          return survey ? {
            key: survey,
            label: survey,
            survey,
            color: "#8fbfff",
          } : null;
        }
        if (!rawPreset || typeof rawPreset !== "object") {
          return null;
        }
        const survey = (
          readSkyPresetField(rawPreset, ["survey", "id", "ID", "obs_id", "creator_did", "hips_service_url", "url"])
          || cleanSkyLayerSurvey(rawPreset.toString && rawPreset.toString !== Object.prototype.toString ? rawPreset.toString() : "")
        );
        if (!survey || survey === "[object Object]") {
          return null;
        }
        const label = (
          readSkyPresetField(rawPreset, ["label", "name", "shortName", "short_name", "obs_title", "title", "description"])
          || survey
        );
        return {
          key: normalizeSkyLayerKey(readSkyPresetField(rawPreset, ["key"]) || survey) || survey,
          label,
          survey,
          color: readSkyPresetField(rawPreset, ["color"]) || fallbackSkyLayerPresetOptions[index % fallbackSkyLayerPresetOptions.length].color || "#8fbfff",
        };
      }

      function buildSkyLayerPresetBySurvey(presets) {
        const presetMap = new Map();
        (Array.isArray(presets) ? presets : []).forEach((rawPreset, index) => {
          const preset = normalizeSkyLayerPreset(rawPreset, index);
          if (preset && preset.survey) {
            presetMap.set(String(preset.survey), preset);
          }
        });
        return presetMap;
      }

      function mergeSkyLayerPresetOptions(primaryPresets, secondaryPresets) {
        const seen = new Set();
        const merged = [];
        [primaryPresets, secondaryPresets].forEach((presetList) => {
          (Array.isArray(presetList) ? presetList : []).forEach((rawPreset, index) => {
            const preset = normalizeSkyLayerPreset(rawPreset, index);
            if (!preset || !preset.survey || seen.has(preset.survey)) {
              return;
            }
            seen.add(preset.survey);
            merged.push(preset);
          });
        });
        return merged;
      }

      function setAladinDefaultSkyLayerPresets(rawPresets) {
        const mergedPresets = mergeSkyLayerPresetOptions(rawPresets, fallbackSkyLayerPresetOptions);
        if (!mergedPresets.length) {
          return;
        }
        const currentValue = skyLayerPresetSelectEl ? skyLayerPresetSelectEl.value : "";
        skyLayerPresetOptions = mergedPresets;
        skyLayerPresetBySurvey = buildSkyLayerPresetBySurvey(skyLayerPresetOptions);
        if (skyLayerPresetSelectEl) {
          skyLayerPresetSelectEl.innerHTML = "";
          syncSkyLayerControls();
          if (currentValue && skyLayerPresetBySurvey.has(currentValue)) {
            skyLayerPresetSelectEl.value = currentValue;
          }
        } else {
          syncSkyLayerControls();
        }
      }

      function skyLayerHipsRegistryUrl() {
        const params = new URLSearchParams();
        params.set("dataproduct_type", "image");
        params.set("fmt", "json");
        params.set("get", "record");
        params.set("MAXREC", "5000");
        params.set("fields", [
          "ID",
          "obs_id",
          "obs_title",
          "obs_regime",
          "hips_service_url",
          "hips_frame",
          "hips_tile_format",
          "hips_order",
        ].join(","));
        return `https://alasky.cds.unistra.fr/MocServer/query?${params.toString()}`;
      }

      function normalizeSkyLayerHipsRegistryRecord(rawRecord, index = 0) {
        const preset = normalizeSkyLayerPreset(rawRecord, index);
        if (!preset || !preset.survey) {
          return null;
        }
        const regime = readSkyPresetField(rawRecord, ["obs_regime", "regime"]);
        const frame = readSkyPresetField(rawRecord, ["hips_frame", "frame"]);
        const tileFormat = readSkyPresetField(rawRecord, ["hips_tile_format", "tile_format"]);
        const title = cleanSkyLayerSurvey(preset.label || preset.survey);
        return {
          ...preset,
          label: title,
          searchText: [
            preset.survey,
            title,
            regime,
            frame,
            tileFormat,
          ].filter(Boolean).join(" ").toLowerCase(),
        };
      }

      function setSkyLayerHipsRegistry(rawRecords) {
        skyLayerHipsRegistry = (Array.isArray(rawRecords) ? rawRecords : [])
          .map((record, index) => normalizeSkyLayerHipsRegistryRecord(record, index))
          .filter(Boolean);
        skyLayerHipsRegistryBySurvey = buildSkyLayerPresetBySurvey(skyLayerHipsRegistry);
        skyLayerHipsRegistryLoaded = skyLayerHipsRegistry.length > 0;
        skyLayerHipsRegistryError = "";
      }

      function fetchSkyLayerHipsRegistry() {
        if (skyLayerHipsRegistryLoaded) {
          return Promise.resolve(skyLayerHipsRegistry);
        }
        if (skyLayerHipsRegistryPromise) {
          return skyLayerHipsRegistryPromise;
        }
        if (skyLayerHipsRegistryFetchAttempted && skyLayerHipsRegistryError) {
          return Promise.resolve(skyLayerHipsRegistry);
        }
        if (typeof fetch !== "function") {
          skyLayerHipsRegistryError = "This browser cannot query the CDS HiPS registry.";
          return Promise.resolve(skyLayerHipsRegistry);
        }
        skyLayerHipsRegistryFetchAttempted = true;
        skyLayerHipsRegistryPromise = fetch(skyLayerHipsRegistryUrl(), { mode: "cors" })
          .then((response) => {
            if (!response || !response.ok) {
              throw new Error(`CDS HiPS registry request failed (${response ? response.status : "network"})`);
            }
            return response.json();
          })
          .then((records) => {
            setSkyLayerHipsRegistry(records);
            return skyLayerHipsRegistry;
          })
          .catch((error) => {
            skyLayerHipsRegistryError = error && error.message ? error.message : String(error || "Unable to load the CDS HiPS registry.");
            return skyLayerHipsRegistry;
          })
          .finally(() => {
            skyLayerHipsRegistryPromise = null;
          });
        return skyLayerHipsRegistryPromise;
      }

      function skyLayerSearchCandidatesForQuery(query) {
        const cleanQuery = cleanSkyLayerSurvey(query).toLowerCase();
        const tokens = cleanQuery.split(/\s+/).filter(Boolean);
        const seen = new Set();
        const candidates = [];
        [skyLayerPresetOptions, skyLayerHipsRegistry].forEach((presetList) => {
          (Array.isArray(presetList) ? presetList : []).forEach((rawPreset, index) => {
            const preset = normalizeSkyLayerPreset(rawPreset, index);
            if (!preset || !preset.survey || seen.has(preset.survey)) {
              return;
            }
            if (rawPreset && rawPreset.searchText) {
              preset.searchText = String(rawPreset.searchText);
            }
            seen.add(preset.survey);
            candidates.push(preset);
          });
        });
        if (!tokens.length) {
          return skyLayerPresetOptions.slice(0, 24);
        }
        return candidates
          .map((preset, index) => {
            const searchText = String(
              preset.searchText
              || [preset.label, preset.survey, preset.key].filter(Boolean).join(" ")
            ).toLowerCase();
            const matches = tokens.every((token) => searchText.includes(token));
            if (!matches) {
              return null;
            }
            const survey = String(preset.survey || "").toLowerCase();
            const label = String(preset.label || "").toLowerCase();
            let score = 0;
            if (survey === cleanQuery || label === cleanQuery) {
              score -= 50;
            } else if (survey.startsWith(cleanQuery) || label.startsWith(cleanQuery)) {
              score -= 25;
            } else if (survey.includes(`/${cleanQuery}`) || label.includes(cleanQuery)) {
              score -= 10;
            }
            if (skyLayerPresetBySurvey.has(preset.survey)) {
              score -= 5;
            }
            return { preset, score: score + index * 0.001 };
          })
          .filter(Boolean)
          .sort((a, b) => a.score - b.score)
          .slice(0, 36)
          .map((entry) => entry.preset);
      }

      function updateSkyLayerSearchOptions(query = "") {
        if (!skyLayerSearchDatalistEl) {
          return;
        }
        const candidates = skyLayerSearchCandidatesForQuery(query);
        const nextSignature = candidates
          .map((preset) => `${preset.survey}|${preset.label || ""}`)
          .join("\\n");
        if (nextSignature === skyLayerSearchRenderSignature) {
          return;
        }
        skyLayerSearchRenderSignature = nextSignature;
        skyLayerSearchDatalistEl.innerHTML = "";
        candidates.forEach((preset) => {
          const option = document.createElement("option");
          option.value = preset.survey;
          option.label = preset.label && preset.label !== preset.survey
            ? `${preset.label} (${preset.survey})`
            : preset.survey;
          option.textContent = option.label;
          skyLayerSearchDatalistEl.appendChild(option);
        });
      }

      function scheduleSkyLayerSearchOptionsUpdate(query = "") {
        if (skyLayerSearchUpdateTimer) {
          window.clearTimeout(skyLayerSearchUpdateTimer);
          skyLayerSearchUpdateTimer = 0;
        }
        skyLayerSearchUpdateTimer = window.setTimeout(() => {
          skyLayerSearchUpdateTimer = 0;
          updateSkyLayerSearchOptions(query);
        }, 140);
      }

      function handleSkyLayerSearchInput() {
        const query = skyLayerCustomInputEl ? skyLayerCustomInputEl.value : "";
        scheduleSkyLayerSearchOptionsUpdate(query);
        if (skyLayerHipsRegistryLoaded || skyLayerHipsRegistryPromise || skyLayerHipsRegistryError) {
          return;
        }
        fetchSkyLayerHipsRegistry().then(() => {
          const currentQuery = skyLayerCustomInputEl ? skyLayerCustomInputEl.value : query;
          scheduleSkyLayerSearchOptionsUpdate(currentQuery);
        });
      }

      function skyLayerPresetForSurvey(survey) {
        const cleanSurvey = cleanSkyLayerSurvey(survey);
        return skyLayerPresetBySurvey.get(cleanSurvey) || skyLayerHipsRegistryBySurvey.get(cleanSurvey) || null;
      }

      function defaultSkyLayerSurvey() {
        return (
          cleanSkyLayerSurvey(skyDomeSpec.survey)
          || cleanSkyLayerSurvey(skyDomeSpec.hips_survey)
          || cleanSkyLayerSurvey(skySpec.survey)
          || "P/DSS2/color"
        );
      }

      function clampSkyLayerOpacity(value, fallback = skyDomeDefaultOpacity) {
        const numeric = Number(value);
        return Math.min(Math.max(Number.isFinite(numeric) ? numeric : Number(fallback), 0.0), 1.0);
      }

      function normalizeSkyLayerStretch(value) {
        const normalized = String(value || "").trim().toLowerCase();
        if (normalized === "log" || normalized === "log10" || normalized === "logarithmic") {
          return "log";
        }
        if (normalized === "asinh" || normalized === "arcsinh") {
          return "asinh";
        }
        return "linear";
      }

      function normalizeSkyLayerColormap(value) {
        const normalized = String(value || "").trim().toLowerCase();
        if (normalized === "gray" || normalized === "grey") {
          return "grayscale";
        }
        return skyLayerColormapOptions.some((option) => option.value === normalized)
          ? normalized
          : "native";
      }

      function skyLayerCutValue(value) {
        if (value === "" || value === null || value === undefined) {
          return "";
        }
        const numeric = Number(value);
        return Number.isFinite(numeric) ? numeric : "";
      }

      function skyLayerCutInputValue(value) {
        const cutValue = skyLayerCutValue(value);
        return cutValue === "" ? "" : String(cutValue);
      }

      function skyLayerFromSurvey(survey, options = {}) {
        const cleanSurvey = cleanSkyLayerSurvey(survey) || defaultSkyLayerSurvey();
        const preset = skyLayerPresetForSurvey(cleanSurvey);
        const key = normalizeSkyLayerKey(options.key || cleanSurvey);
        return {
          key,
          label: cleanSkyLayerSurvey(options.label) || (preset ? preset.label : cleanSurvey),
          survey: cleanSurvey,
          color: cleanSkyLayerSurvey(options.color) || (preset ? preset.color : "#8fbfff"),
          opacity: clampSkyLayerOpacity(options.opacity, skyDomeDefaultOpacity),
          visible: options.visible === undefined ? true : Boolean(options.visible),
          stretch: normalizeSkyLayerStretch(options.stretch),
          colormap: normalizeSkyLayerColormap(options.colormap),
          cutMin: skyLayerCutValue(options.cutMin ?? options.cut_min),
          cutMax: skyLayerCutValue(options.cutMax ?? options.cut_max),
        };
      }

      function normalizeSkyLayer(rawLayer) {
        if (typeof rawLayer === "string") {
          return skyLayerFromSurvey(rawLayer);
        }
        if (!rawLayer || typeof rawLayer !== "object") {
          return null;
        }
        const survey = (
          rawLayer.survey
          || rawLayer.hips_survey
          || rawLayer.id
          || rawLayer.key
          || rawLayer.name
          || ""
        );
        if (!cleanSkyLayerSurvey(survey)) {
          return null;
        }
        return skyLayerFromSurvey(survey, rawLayer);
      }

      function uniqueSkyLayers(layers) {
        const seen = new Set();
        const unique = [];
        (Array.isArray(layers) ? layers : []).forEach((rawLayer) => {
          const layer = normalizeSkyLayer(rawLayer);
          if (!layer || !layer.key || seen.has(layer.key)) {
            return;
          }
          seen.add(layer.key);
          unique.push(layer);
        });
        return unique;
      }

      function savedSkyLayers() {
        const savedGlobalControls = initialState && initialState.global_controls;
        const candidates = [
          initialState.sky_layers,
          savedGlobalControls && savedGlobalControls.sky_layers,
          skyDomeSpec.sky_layers,
          skyDomeSpec.layers,
        ];
        for (const candidate of candidates) {
          if (Array.isArray(candidate) && candidate.length) {
            return uniqueSkyLayers(candidate);
          }
        }
        return [];
      }

      function savedActiveSkyLayerKey() {
        const savedGlobalControls = initialState && initialState.global_controls;
        return normalizeSkyLayerKey(
          initialState.active_sky_layer_key
          || (savedGlobalControls && savedGlobalControls.active_sky_layer_key)
          || skyDomeSpec.active_sky_layer_key
          || skyDomeSpec.active_layer_key
          || ""
        );
      }

      function ensureInitialSkyLayers() {
        if (skyLayerStateInitialized) {
          return;
        }
        skyLayerStateInitialized = true;
        skyLayerState = savedSkyLayers();
        if (!skyLayerState.length && skyDomeControlsAvailable()) {
          skyLayerState = [skyLayerFromSurvey(defaultSkyLayerSurvey(), {
            opacity: skyDomeDefaultOpacity,
            visible: skyDomeDefaultEnabled,
          })];
        }
        const savedActiveKey = savedActiveSkyLayerKey();
        activeSkyLayerKey = (
          savedActiveKey && skyLayerState.some((layer) => layer.key === savedActiveKey)
            ? savedActiveKey
            : (skyLayerState[0] ? skyLayerState[0].key : "")
        );
      }

      function activeSkyLayer() {
        ensureInitialSkyLayers();
        return skyLayerState.find((layer) => layer.key === activeSkyLayerKey) || skyLayerState[0] || null;
      }

      function visibleSkyLayers() {
        ensureInitialSkyLayers();
        return skyLayerState.filter((layer) => (
          layer
          && layer.visible !== false
          && clampSkyLayerOpacity(layer.opacity, skyDomeDefaultOpacity) > 0.0
        ));
      }

      function serializableSkyLayers() {
        ensureInitialSkyLayers();
        return skyLayerState.map((layer) => ({
          key: layer.key,
          label: layer.label,
          survey: layer.survey,
          color: layer.color,
          opacity: clampSkyLayerOpacity(layer.opacity, skyDomeDefaultOpacity),
          visible: layer.visible !== false,
          stretch: normalizeSkyLayerStretch(layer.stretch),
          colormap: normalizeSkyLayerColormap(layer.colormap),
          cut_min: skyLayerCutValue(layer.cutMin),
          cut_max: skyLayerCutValue(layer.cutMax),
        }));
      }

      function postSkyLayerStateToAladin() {
        if (!skyDomeFrameEl || !skyDomeFrameEl.contentWindow) {
          return;
        }
        try {
          skyDomeFrameEl.contentWindow.postMessage({
            type: "oviz-sky-layer-state",
            layers: serializableSkyLayers(),
            activeKey: activeSkyLayerKey || "",
          }, "*");
        } catch (_err) {
        }
      }

      function skyLayerActiveIndex() {
        ensureInitialSkyLayers();
        return skyLayerState.findIndex((layer) => layer.key === activeSkyLayerKey);
      }

      function renderSkyLayerList() {
        if (!skyLayerListEl) {
          return;
        }
        ensureInitialSkyLayers();
        skyLayerListEl.innerHTML = "";
        skyLayerListEl.hidden = false;
        if (!skyLayerState.length) {
          const emptyEl = document.createElement("div");
          emptyEl.className = "oviz-three-sky-layer-order";
          emptyEl.textContent = "No sky layers";
          skyLayerListEl.appendChild(emptyEl);
          return;
        }
        skyLayerState.forEach((layer, index) => {
          const row = document.createElement("details");
          row.className = "oviz-three-sky-layer-row";
          row.dataset.active = layer.key === activeSkyLayerKey ? "true" : "false";
          row.open = layer.key === activeSkyLayerKey || skyLayerState.length <= 2;
          const summaryEl = document.createElement("summary");
          summaryEl.className = "oviz-three-sky-layer-summary";
          const nameEl = document.createElement("div");
          nameEl.className = "oviz-three-sky-layer-name";
          const layerName = String(layer.label || layer.survey || layer.key);
          nameEl.textContent = layerName;
          nameEl.title = String(layer.survey || layer.key || "");
          const orderEl = document.createElement("div");
          orderEl.className = "oviz-three-sky-layer-order";
          orderEl.textContent = index === 0 ? "Top" : (index === skyLayerState.length - 1 ? "Base" : `Layer ${index + 1}`);
          summaryEl.appendChild(nameEl);
          summaryEl.appendChild(orderEl);
          row.appendChild(summaryEl);

          const bodyEl = document.createElement("div");
          bodyEl.className = "oviz-three-sky-layer-body";

          const opacityLabel = document.createElement("label");
          const opacityTextEl = document.createElement("span");
          opacityTextEl.textContent = `Opacity (${clampSkyLayerOpacity(layer.opacity, skyDomeDefaultOpacity).toFixed(2)})`;
          const opacityInput = document.createElement("input");
          opacityInput.type = "range";
          opacityInput.className = "oviz-three-sky-layer-opacity";
          opacityInput.min = "0";
          opacityInput.max = "1";
          opacityInput.step = "0.01";
          opacityInput.value = String(clampSkyLayerOpacity(layer.opacity, skyDomeDefaultOpacity));
          opacityInput.addEventListener("input", () => {
            activeSkyLayerKey = layer.key;
            layer.opacity = clampSkyLayerOpacity(opacityInput.value, skyDomeDefaultOpacity);
            layer.visible = layer.opacity > 0.0;
            opacityTextEl.textContent = `Opacity (${layer.opacity.toFixed(2)})`;
            applySkyLayerState({ forceTiles: false, renderLegend: false, syncControls: false });
          });
          opacityInput.addEventListener("change", () => {
            activeSkyLayerKey = layer.key;
            syncSkyLayerControls();
          });
          opacityLabel.appendChild(opacityTextEl);
          opacityLabel.appendChild(opacityInput);
          bodyEl.appendChild(opacityLabel);

          const stretchGridEl = document.createElement("div");
          stretchGridEl.className = "oviz-three-sky-layer-grid";
          const stretchLabel = document.createElement("label");
          stretchLabel.appendChild(document.createTextNode("Stretch"));
          const stretchSelect = document.createElement("select");
          stretchSelect.className = "oviz-three-sky-layer-stretch";
          skyLayerStretchOptions.forEach((optionDef) => {
            const option = document.createElement("option");
            option.value = optionDef.value;
            option.textContent = optionDef.label;
            stretchSelect.appendChild(option);
          });
          stretchSelect.value = normalizeSkyLayerStretch(layer.stretch);
          stretchSelect.addEventListener("change", () => {
            activeSkyLayerKey = layer.key;
            layer.stretch = normalizeSkyLayerStretch(stretchSelect.value);
            applySkyLayerState({ forceTiles: false, renderLegend: false });
          });
          stretchLabel.appendChild(stretchSelect);
          stretchGridEl.appendChild(stretchLabel);

          const colormapLabel = document.createElement("label");
          colormapLabel.appendChild(document.createTextNode("Colormap"));
          const colormapSelect = document.createElement("select");
          colormapSelect.className = "oviz-three-sky-layer-colormap";
          skyLayerColormapOptions.forEach((optionDef) => {
            const option = document.createElement("option");
            option.value = optionDef.value;
            option.textContent = optionDef.label;
            colormapSelect.appendChild(option);
          });
          colormapSelect.value = normalizeSkyLayerColormap(layer.colormap);
          colormapSelect.addEventListener("change", () => {
            activeSkyLayerKey = layer.key;
            layer.colormap = normalizeSkyLayerColormap(colormapSelect.value);
            applySkyLayerState({ forceTiles: false, renderLegend: false });
          });
          colormapLabel.appendChild(colormapSelect);
          stretchGridEl.appendChild(colormapLabel);
          bodyEl.appendChild(stretchGridEl);

          const cutGridEl = document.createElement("div");
          cutGridEl.className = "oviz-three-sky-layer-grid";
          const minLabel = document.createElement("label");
          minLabel.appendChild(document.createTextNode("Min"));
          const minInput = document.createElement("input");
          minInput.type = "number";
          minInput.className = "oviz-three-sky-layer-cut-min";
          minInput.step = "any";
          minInput.placeholder = "auto";
          minInput.value = skyLayerCutInputValue(layer.cutMin);
          minInput.addEventListener("change", () => {
            activeSkyLayerKey = layer.key;
            layer.cutMin = skyLayerCutValue(minInput.value);
            applySkyLayerState({ forceTiles: false, renderLegend: false });
          });
          minLabel.appendChild(minInput);
          cutGridEl.appendChild(minLabel);
          const maxLabel = document.createElement("label");
          maxLabel.appendChild(document.createTextNode("Max"));
          const maxInput = document.createElement("input");
          maxInput.type = "number";
          maxInput.className = "oviz-three-sky-layer-cut-max";
          maxInput.step = "any";
          maxInput.placeholder = "auto";
          maxInput.value = skyLayerCutInputValue(layer.cutMax);
          maxInput.addEventListener("change", () => {
            activeSkyLayerKey = layer.key;
            layer.cutMax = skyLayerCutValue(maxInput.value);
            applySkyLayerState({ forceTiles: false, renderLegend: false });
          });
          maxLabel.appendChild(maxInput);
          cutGridEl.appendChild(maxLabel);
          bodyEl.appendChild(cutGridEl);

          const actionEl = document.createElement("div");
          actionEl.className = "oviz-three-sky-layer-actions";
          const upButton = document.createElement("button");
          upButton.type = "button";
          upButton.textContent = "↑";
          upButton.title = "Move image up";
          upButton.setAttribute("aria-label", "Move sky image up");
          upButton.disabled = index === 0;
          upButton.addEventListener("click", () => {
            activeSkyLayerKey = layer.key;
            moveActiveSkyLayer(-1);
            focusViewer();
          });
          const downButton = document.createElement("button");
          downButton.type = "button";
          downButton.textContent = "↓";
          downButton.title = "Move image down";
          downButton.setAttribute("aria-label", "Move sky image down");
          downButton.disabled = index >= skyLayerState.length - 1;
          downButton.addEventListener("click", () => {
            activeSkyLayerKey = layer.key;
            moveActiveSkyLayer(1);
            focusViewer();
          });
          const removeButton = document.createElement("button");
          removeButton.type = "button";
          removeButton.textContent = "Remove";
          removeButton.title = "Remove image";
          removeButton.setAttribute("aria-label", `Remove ${layerName}`);
          removeButton.addEventListener("click", () => {
            activeSkyLayerKey = layer.key;
            removeActiveSkyLayer();
            focusViewer();
          });
          actionEl.appendChild(upButton);
          actionEl.appendChild(downButton);
          actionEl.appendChild(removeButton);
          bodyEl.appendChild(actionEl);
          row.appendChild(bodyEl);
          skyLayerListEl.appendChild(row);
        });
      }

      function syncSkyLayerControls() {
        ensureInitialSkyLayers();
        if (skyLayerPresetSelectEl && !skyLayerPresetSelectEl.options.length) {
          skyLayerPresetOptions.forEach((preset) => {
            const option = document.createElement("option");
            option.value = preset.survey;
            option.textContent = preset.label;
            skyLayerPresetSelectEl.appendChild(option);
          });
        }
        renderSkyLayerList();
      }

      function applySkyLayerState(options = {}) {
        ensureInitialSkyLayers();
        const visibleLayers = visibleSkyLayers();
        const topLayer = visibleLayers[0] || null;
        const baseLayer = visibleLayers[visibleLayers.length - 1] || null;
        if (!activeSkyLayerKey && skyLayerState.length) {
          activeSkyLayerKey = skyLayerState[0].key;
        }
        if (baseLayer) {
          skyDomeSpec.enabled = true;
          skyDomeSpec.opacity = 1.0;
          skyDomeSpec.survey = baseLayer.survey;
          skySpec.survey = baseLayer.survey;
          skyDomeSurvey = String((topLayer && (topLayer.label || topLayer.survey)) || baseLayer.label || baseLayer.survey || "Sky");
        } else {
          skyDomeSpec.enabled = false;
          skyDomeSpec.opacity = 0.0;
        }
        if (options.syncControls !== false) {
          syncSkyLayerControls();
        }
        postSkyLayerStateToAladin();
        if (options.update !== false && typeof updateSkyDomeFromControls === "function") {
          updateSkyDomeFromControls({
            forceTiles: options.forceTiles !== false,
            syncControls: options.syncControls !== false,
          });
        }
        if (options.renderLegend !== false && typeof renderLegend === "function") {
          renderLegend();
        }
      }

      function addOrActivateSkyLayer(survey) {
        ensureInitialSkyLayers();
        const layer = skyLayerFromSurvey(survey || defaultSkyLayerSurvey(), {
          opacity: skyDomeDefaultOpacity,
          visible: true,
        });
        const existing = skyLayerState.find((candidate) => candidate.key === layer.key);
        if (existing) {
          existing.visible = true;
          activeSkyLayerKey = existing.key;
        } else {
          skyLayerState.unshift(layer);
          activeSkyLayerKey = layer.key;
        }
        applySkyLayerState({ forceTiles: true });
      }

      function removeActiveSkyLayer() {
        ensureInitialSkyLayers();
        if (!activeSkyLayerKey) {
          return;
        }
        const index = skyLayerState.findIndex((layer) => layer.key === activeSkyLayerKey);
        if (index < 0) {
          return;
        }
        skyLayerState.splice(index, 1);
        const nextLayer = skyLayerState[Math.min(index, Math.max(skyLayerState.length - 1, 0))] || null;
        activeSkyLayerKey = nextLayer ? nextLayer.key : "";
        applySkyLayerState({ forceTiles: true });
      }

      function moveActiveSkyLayer(offset) {
        ensureInitialSkyLayers();
        const index = skyLayerActiveIndex();
        if (index < 0) {
          return;
        }
        const targetIndex = Math.min(Math.max(index + Math.round(Number(offset) || 0), 0), skyLayerState.length - 1);
        if (targetIndex === index) {
          syncSkyLayerControls();
          return;
        }
        const layer = skyLayerState[index];
        skyLayerState.splice(index, 1);
        skyLayerState.splice(targetIndex, 0, layer);
        activeSkyLayerKey = layer.key;
        applySkyLayerState({ forceTiles: true });
      }

      function skyPanelProjectionName(source) {
        const projection = String((source && source.projection) || (skyDomeSpec && skyDomeSpec.projection) || "CAR").trim().toUpperCase();
        if (projection === "MOL" || projection === "MOLLWEIDE") {
          return "MOL";
        }
        return "CAR";
      }

      function skyWrap01(value) {
        return value - Math.floor(value);
      }

      function skyPanelSolveMollweideTheta(latRad) {
        if (Math.abs(Math.abs(latRad) - (0.5 * Math.PI)) < 0.0001) {
          return Math.sign(latRad) * 0.5 * Math.PI;
        }
        let theta = latRad;
        for (let index = 0; index < 8; index += 1) {
          const f = (2.0 * theta) + Math.sin(2.0 * theta) - (Math.PI * Math.sin(latRad));
          const fp = Math.max(2.0 + (2.0 * Math.cos(2.0 * theta)), 0.0001);
          theta -= f / fp;
        }
        return theta;
      }

      function skyPanelUvForLonLat(lDeg, bDeg, projectionName) {
        const lon = wrapLongitudeDeltaDeg(Number(lDeg)) * Math.PI / 180.0;
        const lat = Math.min(Math.max(Number(bDeg), -90.0), 90.0) * Math.PI / 180.0;
        if (!Number.isFinite(lon) || !Number.isFinite(lat)) {
          return null;
        }
        if (skyPanelProjectionName({ projection: projectionName }) === "MOL") {
          const rootTwo = Math.SQRT2;
          const theta = skyPanelSolveMollweideTheta(lat);
          const x = (2.0 * rootTwo / Math.PI) * lon * Math.cos(theta);
          const y = rootTwo * Math.sin(theta);
          return {
            u: 0.5 - (x / (4.0 * rootTwo)),
            v: 0.5 - (y / (2.0 * rootTwo)),
          };
        }
        return {
          u: skyWrap01(0.5 - (normalizeSkyLongitude(lDeg) / 360.0)),
          v: Math.min(Math.max(0.5 - (Number(bDeg) / 180.0), 0.0), 1.0),
        };
      }

      function skyPanelLonLatForPoint(point) {
        const lDeg = Number(point && point.l);
        const bDeg = Number(point && point.b);
        if (Number.isFinite(lDeg) && Number.isFinite(bDeg)) {
          return { l: normalizeSkyLongitude(lDeg), b: bDeg };
        }
        const raDeg = Number(point && point.ra);
        const decDeg = Number(point && point.dec);
        if (Number.isFinite(raDeg) && Number.isFinite(decDeg)) {
          return galacticDegFromIcrsDeg(raDeg, decDeg);
        }
        return null;
      }

      function skyPanelLoadImage(source) {
        if (!source || !source.dataUrl) {
          return Promise.reject(new Error("No embedded sky image is available."));
        }
        const cacheKey = source.key;
        const cached = skyPanelImageCache.get(cacheKey);
        if (cached) {
          return cached;
        }
        const pending = new Promise((resolve, reject) => {
          const image = new Image();
          image.onload = () => resolve(image);
          image.onerror = () => reject(new Error("The embedded sky image could not be decoded."));
          image.src = source.dataUrl;
        });
        skyPanelImageCache.set(cacheKey, pending);
        return pending;
      }

      function resizeSkyPanelCanvas() {
        if (!skyCanvasEl) {
          return null;
        }
        const rect = skyCanvasEl.getBoundingClientRect();
        const width = Math.max(1, Math.floor(rect.width || 1));
        const height = Math.max(1, Math.floor(rect.height || 1));
        const dpr = Math.min(Math.max(window.devicePixelRatio || 1, 1), 2);
        const pixelWidth = Math.max(1, Math.floor(width * dpr));
        const pixelHeight = Math.max(1, Math.floor(height * dpr));
        if (skyCanvasEl.width !== pixelWidth || skyCanvasEl.height !== pixelHeight) {
          skyCanvasEl.width = pixelWidth;
          skyCanvasEl.height = pixelHeight;
        }
        const ctx = skyCanvasEl.getContext("2d");
        if (!ctx) {
          return null;
        }
        ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
        return { ctx, width, height, dpr };
      }

      function skyPanelFitImageRect(width, height, image) {
        const imageWidth = Math.max(1, Number(image.naturalWidth || image.width || 1));
        const imageHeight = Math.max(1, Number(image.naturalHeight || image.height || 1));
        const scale = Math.min(width / imageWidth, height / imageHeight);
        const drawWidth = imageWidth * scale;
        const drawHeight = imageHeight * scale;
        return {
          x: (width - drawWidth) * 0.5,
          y: (height - drawHeight) * 0.5,
          width: drawWidth,
          height: drawHeight,
          imageWidth,
          imageHeight,
        };
      }

      function skyPanelDrawWrappedImage(ctx, image, srcX, srcY, srcW, srcH, dstX, dstY, dstW, dstH) {
        const imageWidth = Math.max(1, Number(image.naturalWidth || image.width || 1));
        const sourceY = Math.min(Math.max(srcY, 0), Math.max(0, image.height - srcH));
        let remaining = Math.min(srcW, imageWidth);
        let currentSrcX = ((srcX % imageWidth) + imageWidth) % imageWidth;
        let drawnSourceW = 0.0;
        while (remaining > 0.01) {
          const sliceW = Math.min(remaining, imageWidth - currentSrcX);
          const sliceDstX = dstX + (drawnSourceW / srcW) * dstW;
          const sliceDstW = (sliceW / srcW) * dstW;
          ctx.drawImage(image, currentSrcX, sourceY, sliceW, srcH, sliceDstX, dstY, sliceDstW, dstH);
          remaining -= sliceW;
          drawnSourceW += sliceW;
          currentSrcX = 0;
        }
      }

      function skyPanelCanvasPointForLonLat(lDeg, bDeg, drawState) {
        const uv = skyPanelUvForLonLat(lDeg, bDeg, drawState.projection);
        if (!uv || uv.u < -0.05 || uv.u > 1.05 || uv.v < -0.05 || uv.v > 1.05) {
          return null;
        }
        if (drawState.mode === "click") {
          const imageX = uv.u * drawState.imageWidth;
          const imageY = uv.v * drawState.imageHeight;
          let relX = imageX - drawState.srcX;
          const halfImageWidth = drawState.imageWidth * 0.5;
          if (relX < -halfImageWidth) {
            relX += drawState.imageWidth;
          } else if (relX > halfImageWidth) {
            relX -= drawState.imageWidth;
          }
          const relY = imageY - drawState.srcY;
          return {
            x: (relX / drawState.srcW) * drawState.width,
            y: (relY / drawState.srcH) * drawState.height,
          };
        }
        return {
          x: drawState.rect.x + (uv.u * drawState.rect.width),
          y: drawState.rect.y + (uv.v * drawState.rect.height),
        };
      }

      function drawSkyPanelCatalog(ctx, payload, drawState, mode) {
        const catalogs = Array.isArray(payload) ? payload : [];
        catalogs.forEach((catalog) => {
          const points = Array.isArray(catalog.points) ? catalog.points : [];
          const color = String(catalog.color || theme.footprint || "#6ec5ff");
          points.forEach((point) => {
            const lonLat = skyPanelLonLatForPoint(point);
            if (!lonLat) {
              return;
            }
            const canvasPoint = skyPanelCanvasPointForLonLat(lonLat.l, lonLat.b, drawState);
            if (!canvasPoint) {
              return;
            }
            if (
              canvasPoint.x < -10
              || canvasPoint.x > drawState.width + 10
              || canvasPoint.y < -10
              || canvasPoint.y > drawState.height + 10
            ) {
              return;
            }
            ctx.save();
            ctx.globalAlpha = 0.92;
            ctx.fillStyle = color;
            ctx.strokeStyle = "rgba(255, 255, 255, 0.86)";
            ctx.lineWidth = mode === "click" ? 1.6 : 1.0;
            ctx.beginPath();
            ctx.arc(canvasPoint.x, canvasPoint.y, mode === "click" ? 4.2 : 2.8, 0, Math.PI * 2.0);
            ctx.fill();
            ctx.stroke();
            ctx.restore();
          });
        });
      }

      function drawSkyPanelSelectionBeams(ctx, selections, drawState, focus) {
        const radiusDeg = Math.max(Number(skySpec.radius_deg || 1.0), 0.01);
        uniqueSelections(selections).forEach((selection) => {
          const lonLat = {
            l: Number(selection.l_deg),
            b: Number(selection.b_deg),
          };
          if (!Number.isFinite(lonLat.l) || !Number.isFinite(lonLat.b)) {
            return;
          }
          const canvasPoint = skyPanelCanvasPointForLonLat(lonLat.l, lonLat.b, drawState);
          if (!canvasPoint) {
            return;
          }
          const fovDeg = focus && Number.isFinite(Number(focus.fovDeg)) ? Number(focus.fovDeg) : 360.0;
          const beamRadius = drawState.mode === "click"
            ? Math.max(7, Math.min(drawState.width, drawState.height) * radiusDeg / Math.max(fovDeg, 0.01))
            : Math.max(3, drawState.rect.width * radiusDeg / 360.0);
          ctx.save();
          ctx.strokeStyle = String(theme.footprint || "#6ec5ff");
          ctx.lineWidth = 1.5;
          ctx.globalAlpha = 0.92;
          ctx.beginPath();
          ctx.arc(canvasPoint.x, canvasPoint.y, beamRadius, 0, Math.PI * 2.0);
          ctx.stroke();
          ctx.restore();
        });
      }

      function renderLocalSkyPanel(selections, catalogPayload, mode = "overview", volumeOverlay = null) {
        syncSkyPanelSourceSelector();
        const canvasState = resizeSkyPanelCanvas();
        if (!canvasState) {
          return;
        }
        const { ctx, width, height } = canvasState;
        ctx.clearRect(0, 0, width, height);
        ctx.fillStyle = "#000000";
        ctx.fillRect(0, 0, width, height);
        const source = activeSkyDomeSourceOption();
        if (skyReadoutEl) {
          skyReadoutEl.textContent = skyReadoutText(selections, catalogPayload, mode, volumeOverlay).split("\\n")[0] || "";
        }
        if (!source) {
          if (skyStatusEl) {
            skyStatusEl.textContent = "No embedded sky image";
          }
          return;
        }
        if (skyStatusEl) {
          skyStatusEl.textContent = "Loading " + source.label;
        }
        const renderSerial = ++skyPanelRenderSerial;
        skyPanelLoadImage(source)
          .then((image) => {
            if (renderSerial !== skyPanelRenderSerial) {
              return;
            }
            ctx.clearRect(0, 0, width, height);
            ctx.fillStyle = "#000000";
            ctx.fillRect(0, 0, width, height);
            const projection = skyPanelProjectionName(source);
            const imageWidth = Math.max(1, Number(image.naturalWidth || image.width || 1));
            const imageHeight = Math.max(1, Number(image.naturalHeight || image.height || 1));
            const focus = mode === "click" ? skyFocusFromPayload(uniqueSelections(selections), catalogPayload) : null;
            let drawState = {
              mode: "overview",
              projection,
              width,
              height,
              rect: skyPanelFitImageRect(width, height, image),
              imageWidth,
              imageHeight,
            };
            if (focus) {
              const galacticFocus = galacticDegFromIcrsDeg(Number(focus.ra), Number(focus.dec));
              const uv = galacticFocus ? skyPanelUvForLonLat(galacticFocus.l, galacticFocus.b, projection) : null;
              if (uv) {
                const aspect = width / Math.max(height, 1);
                const fovDeg = Math.min(Math.max(Number(focus.fovDeg) || 12.0, 1.0), 120.0);
                let srcW = Math.max(32, imageWidth * fovDeg / 360.0);
                let srcH = Math.max(32, imageHeight * fovDeg / 180.0);
                if (srcW / Math.max(srcH, 1) < aspect) {
                  srcW = srcH * aspect;
                } else {
                  srcH = srcW / aspect;
                }
                srcW = Math.min(srcW, imageWidth);
                srcH = Math.min(srcH, imageHeight);
                const srcX = (uv.u * imageWidth) - (srcW * 0.5);
                const srcY = Math.min(Math.max((uv.v * imageHeight) - (srcH * 0.5), 0), Math.max(0, imageHeight - srcH));
                skyPanelDrawWrappedImage(ctx, image, srcX, srcY, srcW, srcH, 0, 0, width, height);
                drawState = {
                  mode: "click",
                  projection,
                  width,
                  height,
                  imageWidth,
                  imageHeight,
                  srcX,
                  srcY,
                  srcW,
                  srcH,
                  rect: { x: 0, y: 0, width, height },
                };
              } else {
                ctx.drawImage(image, drawState.rect.x, drawState.rect.y, drawState.rect.width, drawState.rect.height);
              }
            } else {
              ctx.drawImage(image, drawState.rect.x, drawState.rect.y, drawState.rect.width, drawState.rect.height);
            }
            drawSkyPanelCatalog(ctx, catalogPayload, drawState, mode);
            drawSkyPanelSelectionBeams(ctx, selections, drawState, focus);
            if (skyStatusEl) {
              const projectionLabel = projection === "MOL" ? "Mollweide" : "plate carree";
              skyStatusEl.textContent = `${source.label} | ${projectionLabel} | embedded`;
            }
          })
          .catch((err) => {
            if (renderSerial !== skyPanelRenderSerial) {
              return;
            }
            if (skyStatusEl) {
              skyStatusEl.textContent = err && err.message ? err.message : String(err || "Sky image unavailable");
            }
          });
      }

      function buildAladinSrcdoc(selections, catalogPayload, mode = "overview", volumeOverlay = null, captureOnly = false) {
        const activeSelections = uniqueSelections(selections);
        const requestedMode = mode === "click" ? "click" : "overview";
        const captureOnlyMode = Boolean(captureOnly);
        const focus = requestedMode === "click" ? skyFocusFromPayload(activeSelections, catalogPayload) : null;
        const actualMode = requestedMode === "click" && focus ? "click" : "overview";
        const bg = String(theme.panel_solid || theme.paper_bgcolor || "#121212");
        const txt = String(theme.text_color || "#d0d0d0");
        const beamColor = JSON.stringify(String(theme.footprint || "#6ec5ff"));
        const survey = JSON.stringify(String(skySpec.survey || "P/DSS2/color"));
        const cooFrame = JSON.stringify(
          actualMode === "click"
            ? (String(skySpec.frame || "galactic") === "galactic" ? "galactic" : "equatorial")
            : "galactic"
        );
        const payloadJson = JSON.stringify(catalogPayload || []);
        const volumeOverlayJson = JSON.stringify(volumeOverlay || null);
        const captureOnlyJson = JSON.stringify(captureOnlyMode);
        const skyDomeCaptureEnabledJson = JSON.stringify(captureOnlyMode && Boolean(skyDomeSpec && skyDomeSpec.enabled));
        const skyDomeBackgroundOnlyJson = JSON.stringify(
          captureOnlyMode
          && skyDomeUsesAladinBackground()
        );
        const skyDomeCaptureWidth = Math.max(512, Math.min(8192, Math.round(Number(skyDomeSpec.capture_width_px) || 4096)));
        const skyDomeCaptureHeight = Math.max(256, Math.min(4096, Math.round(Number(skyDomeSpec.capture_height_px) || 2048)));
        const skyDomeProjectionJson = JSON.stringify(String((skyDomeSpec && skyDomeSpec.projection) || "MOL").toUpperCase());
        const skyDomeCaptureFormatJson = JSON.stringify(
          String((skyDomeSpec && skyDomeSpec.capture_format) || "image/jpeg") === "image/png"
            ? "image/png"
            : "image/jpeg"
        );
        const skyDomeCaptureQuality = Math.min(
          Math.max(Number((skyDomeSpec && skyDomeSpec.capture_quality) || 0.94), 0.1),
          1.0
        );
        const beamCentersJson = JSON.stringify((actualMode === "click" ? activeSelections : [])
          .map((selection) => ({
            ra: Number(selection.ra_deg),
            dec: Number(selection.dec_deg),
          }))
          .filter((beam) => Number.isFinite(beam.ra) && Number.isFinite(beam.dec)));
        const modeJson = JSON.stringify(actualMode);
        const target = JSON.stringify(
          actualMode === "click"
            ? `${Number(focus.ra).toFixed(6)} ${Number(focus.dec).toFixed(6)}`
            : "0.000000 0.000000"
        );
        const radiusDeg = Number(skySpec.radius_deg || 1.0);
        const fovDeg = actualMode === "click" ? Number(focus.fovDeg) : 360.0;
        const ra = actualMode === "click" ? Number(focus.ra) : 0.0;
        const dec = actualMode === "click" ? Number(focus.dec) : 0.0;
        return `<!doctype html>
<html>
  <head>
    <meta charset="utf-8" />
    <meta name="viewport" content="width=device-width, initial-scale=1" />
    <style>
      html, body { margin: 0; padding: 0; width: 100%; height: 100%; background: ${bg}; color: ${txt}; overflow: hidden; }
      #oviz-wrap { position: relative; width: 100%; height: 100%; }
      #aladin-lite-div {
        width: 100%;
        height: 100%;
      }
      .aladin-stack-box {
        max-height: min(56vh, 420px) !important;
        overflow-y: auto !important;
        overscroll-behavior: contain;
      }
      .aladin-logo,
      .aladin-logo-container,
      .aladin-location,
      .aladin-coordinates,
      .aladin-fov,
      .aladin-status-bar,
      .aladin-frameChoice,
      .aladin-projectionChoice {
        display: none !important;
      }
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
        const viewMode = ${modeJson};
        const payload = ${payloadJson};
        const volumeOverlay = ${volumeOverlayJson};
        const skyDomeCaptureOnly = ${captureOnlyJson};
        const skyDomeCaptureEnabled = ${skyDomeCaptureEnabledJson};
        const skyDomeBackgroundOnly = ${skyDomeBackgroundOnlyJson};
        const skyDomeProjection = ${skyDomeProjectionJson};
        const skyDomeCaptureWidth = ${skyDomeCaptureWidth};
        const skyDomeCaptureHeight = ${skyDomeCaptureHeight};
        const skyDomeCaptureFormat = ${skyDomeCaptureFormatJson};
        const skyDomeCaptureQuality = ${skyDomeCaptureQuality};
        const beamCenters = ${beamCentersJson};
        const statusEl = document.getElementById("oviz-status");
        const catalogs = [];
        const markerCanvasCache = new Map();
        let hoveredClusterKey = "";
        let aladinInstance = null;
        let pendingSkyBackgroundView = null;
        let skyBackgroundViewAnimationFrame = 0;
        let lastAppliedSkyBackgroundSignature = "";
        let activeImageSurvey = ${survey};
        let activeSkyLayerStackSignature = "";
        const managedSkyOverlayLayerNames = new Set();
        function normalizeClusterKey(value) {
          return String(value || "")
            .trim()
            .toLowerCase()
            .replace(/[_\s]+/g, " ");
        }
        function markerCanvas(color, diameterPx, opacity) {
          const size = Math.max(3, Math.round(Number(diameterPx) || 0));
          const alpha = Math.min(Math.max(Number(opacity), 0.0), 1.0);
          const cacheKey = [String(color || "#ffffff"), size, alpha.toFixed(3)].join("|");
          if (markerCanvasCache.has(cacheKey)) {
            return markerCanvasCache.get(cacheKey);
          }
          const canvas = document.createElement("canvas");
          const pixelSize = Math.max(12, size * 4);
          canvas.width = pixelSize;
          canvas.height = pixelSize;
          const ctx = canvas.getContext("2d");
          if (ctx) {
            const radius = (pixelSize * 0.5) - 2.0;
            ctx.clearRect(0, 0, pixelSize, pixelSize);
            ctx.fillStyle = String(color || "#ffffff");
            ctx.globalAlpha = alpha;
            ctx.beginPath();
            ctx.arc(pixelSize * 0.5, pixelSize * 0.5, radius, 0, Math.PI * 2.0);
            ctx.fill();
            ctx.globalAlpha = Math.min(1.0, alpha + 0.18);
            ctx.lineWidth = Math.max(1.0, pixelSize * 0.06);
            ctx.strokeStyle = "rgba(255,255,255,0.65)";
            ctx.stroke();
          }
          markerCanvasCache.set(cacheKey, canvas);
          return canvas;
        }
        function refreshCatalogHoverShapes() {
          catalogs.forEach((cat) => {
            if (!cat || !cat.__ovizShapeFn) {
              return;
            }
            if (typeof cat.setShape === "function") {
              cat.setShape(cat.__ovizShapeFn);
            } else if (typeof cat.updateShape === "function") {
              cat.updateShape({ shape: cat.__ovizShapeFn });
            }
          });
        }
        function postHoveredClusterKey(clusterKey) {
          if (!window.parent || window.parent === window) {
            return;
          }
          window.parent.postMessage({
            type: "oviz-sky-hover-cluster",
            clusterKey: clusterKey || null,
          }, "*");
        }
        function postSkyDomeSnapshot(status, dataUrl, message, projectionOverride) {
          if (!window.parent || window.parent === window || !skyDomeCaptureEnabled || viewMode !== "overview") {
            return;
          }
          window.parent.postMessage({
            type: "oviz-aladin-sky-snapshot",
            status: String(status || ""),
            survey: ${survey},
            mode: viewMode,
            projection: String(projectionOverride || skyDomeProjection || "MOL"),
            dataUrl: dataUrl || null,
            message: message ? String(message) : "",
          }, "*");
        }
        function validImageDataUrl(value) {
          return (
            typeof value === "string"
            && value.startsWith("data:image/")
            && value.length > 1000
          );
        }
        function scheduleSkyDomeSnapshot(aladin, delayMs) {
          if (!skyDomeCaptureEnabled || !aladin || typeof aladin.getViewDataURL !== "function") {
            return;
          }
          window.setTimeout(() => captureSkyDomeSnapshot(aladin), Math.max(Number(delayMs) || 0, 0));
        }
        function captureSkyDomeSnapshot(aladin) {
          if (skyDomeBackgroundOnly || !skyDomeCaptureEnabled || !aladin || typeof aladin.getViewDataURL !== "function") {
            return;
          }
          const delaysMs = [1200, 2800, 5200, 9000];
          let attemptIndex = 0;
          const retryOrFail = (message) => {
            if (attemptIndex < delaysMs.length) {
              const delayMs = delaysMs[attemptIndex];
              attemptIndex += 1;
              window.setTimeout(tryCapture, delayMs);
            } else {
              postSkyDomeSnapshot("unavailable", null, message);
            }
          };
          const tryCapture = () => {
            try {
              const result = aladin.getViewDataURL({
                width: skyDomeCaptureWidth,
                height: skyDomeCaptureHeight,
                format: skyDomeCaptureFormat,
                quality: skyDomeCaptureQuality,
                logo: false,
              });
              Promise.resolve(result)
                .then((dataUrl) => {
                  if (validImageDataUrl(dataUrl)) {
                    postSkyDomeSnapshot("ready", dataUrl, "", skyDomeProjection);
                  } else {
                    retryOrFail("Aladin sky capture returned an empty image.");
                  }
                })
                .catch((err) => {
                  retryOrFail(err && err.message ? err.message : String(err));
                });
            } catch (err) {
              retryOrFail(err && err.message ? err.message : String(err));
            }
          };
          retryOrFail("Aladin sky capture did not become available.");
        }
        function postSkyBackgroundReady() {
          if (!window.parent || window.parent === window || !skyDomeBackgroundOnly) {
            return;
          }
          window.parent.postMessage({
            type: "oviz-aladin-sky-background-ready",
            survey: ${survey},
          }, "*");
        }
        function aladinFavoriteValue(favorite, fieldNames) {
          const fields = Array.isArray(fieldNames) ? fieldNames : [fieldNames];
          for (const field of fields) {
            if (!field) {
              continue;
            }
            if (favorite && Object.prototype.hasOwnProperty.call(favorite, field) && favorite[field] != null) {
              const value = String(favorite[field]).trim();
              if (value) {
                return value;
              }
            }
            const metadata = favorite && (favorite.meta || favorite.metadata || favorite.properties || favorite.hips);
            if (metadata && Object.prototype.hasOwnProperty.call(metadata, field) && metadata[field] != null) {
              const value = String(metadata[field]).trim();
              if (value) {
                return value;
              }
            }
          }
          return "";
        }
        function normalizeAladinFavorite(favorite) {
          if (typeof favorite === "string") {
            const survey = favorite.trim();
            return survey ? { survey, label: survey } : null;
          }
          if (!favorite || typeof favorite !== "object") {
            return null;
          }
          const survey = aladinFavoriteValue(favorite, ["survey", "id", "ID", "obs_id", "creator_did", "hips_service_url", "url"]);
          if (!survey) {
            return null;
          }
          const label = aladinFavoriteValue(favorite, ["label", "name", "shortName", "short_name", "obs_title", "title", "description"]) || survey;
          return {
            key: survey,
            label,
            survey,
            frame: aladinFavoriteValue(favorite, ["hips_frame", "hips_frame_name", "frame"]),
            order: aladinFavoriteValue(favorite, ["hips_order", "order"]),
            format: aladinFavoriteValue(favorite, ["hips_tile_format", "format"]),
          };
        }
        function postAladinHipsFavorites() {
          if (!window.parent || window.parent === window || !aladinInstance) {
            return;
          }
          const rawFavorites = Array.isArray(aladinInstance.hipsFavorites) ? aladinInstance.hipsFavorites : [];
          const favorites = rawFavorites
            .map((favorite) => normalizeAladinFavorite(favorite))
            .filter((favorite) => favorite && favorite.survey);
          if (!favorites.length) {
            return;
          }
          window.parent.postMessage({
            type: "oviz-aladin-hips-favorites",
            favorites,
          }, "*");
        }
        function applySkyBackgroundViewNow(data) {
          if (!skyDomeBackgroundOnly || !aladinInstance || !data || typeof data !== "object") {
            return;
          }
          const ra = Number(data.ra);
          const dec = Number(data.dec);
          const lDeg = Number(data.l);
          const bDeg = Number(data.b);
          const fovDeg = Number(data.fovDeg);
          const signature = [
            Number.isFinite(lDeg) ? lDeg.toFixed(5) : "",
            Number.isFinite(bDeg) ? bDeg.toFixed(5) : "",
            Number.isFinite(ra) ? ra.toFixed(5) : "",
            Number.isFinite(dec) ? dec.toFixed(5) : "",
            Number.isFinite(fovDeg) ? fovDeg.toFixed(3) : "",
          ].join("|");
          if (signature === lastAppliedSkyBackgroundSignature) {
            return;
          }
          lastAppliedSkyBackgroundSignature = signature;
          if (typeof aladinInstance.stopAnimation === "function") {
            aladinInstance.stopAnimation();
          }
          if (Number.isFinite(fovDeg) && fovDeg > 0.0 && typeof aladinInstance.setFoV === "function") {
            aladinInstance.setFoV(Math.min(Math.max(fovDeg, 0.05), 179.0));
          }
          if (Number.isFinite(lDeg) && Number.isFinite(bDeg) && typeof aladinInstance.gotoPosition === "function") {
            aladinInstance.gotoPosition(lDeg, bDeg);
          } else if (Number.isFinite(ra) && Number.isFinite(dec) && typeof aladinInstance.gotoRaDec === "function") {
            aladinInstance.gotoRaDec(ra, dec);
          }
        }
        function skyLayerNameFor(layer) {
          return "oviz-" + String((layer && layer.key) || (layer && layer.survey) || "layer")
            .replace(/[^a-zA-Z0-9_-]+/g, "-")
            .slice(0, 80);
        }
        function normalizeSkyLayerStretch(value) {
          const normalized = String(value || "").trim().toLowerCase();
          if (normalized === "log" || normalized === "log10" || normalized === "logarithmic") {
            return "log";
          }
          if (normalized === "asinh" || normalized === "arcsinh") {
            return "asinh";
          }
          return "linear";
        }
        function normalizeSkyLayerColormap(value) {
          const normalized = String(value || "").trim().toLowerCase();
          return normalized && normalized !== "native" ? normalized : "";
        }
        function numericSkyLayerCut(value) {
          if (value === "" || value === null || value === undefined) {
            return null;
          }
          const numeric = Number(value);
          return Number.isFinite(numeric) ? numeric : null;
        }
        function removeManagedSkyOverlays() {
          if (!aladinInstance) {
            managedSkyOverlayLayerNames.clear();
            return;
          }
          managedSkyOverlayLayerNames.forEach((layerName) => {
            try {
              if (typeof aladinInstance.removeImageLayer === "function") {
                aladinInstance.removeImageLayer(layerName);
              } else if (typeof aladinInstance.removeOverlayImageLayer === "function") {
                aladinInstance.removeOverlayImageLayer(layerName);
              }
            } catch (_err) {
            }
          });
          managedSkyOverlayLayerNames.clear();
        }
        function imageLayerForName(layerName, isBase) {
          if (!aladinInstance) {
            return null;
          }
          try {
            if (isBase && typeof aladinInstance.getBaseImageLayer === "function") {
              return aladinInstance.getBaseImageLayer();
            }
            if (typeof aladinInstance.getOverlayImageLayer === "function") {
              return aladinInstance.getOverlayImageLayer(layerName);
            }
            if (typeof aladinInstance.getImageLayer === "function") {
              return aladinInstance.getImageLayer(layerName);
            }
          } catch (_err) {
          }
          return null;
        }
        function applySkyImageLayerOptions(layerName, layer, isBase) {
          const imageLayer = imageLayerForName(layerName, isBase);
          if (!imageLayer) {
            return;
          }
          const opacity = Math.min(Math.max(Number(layer && layer.opacity), 0.0), 1.0);
          const visibleOpacity = layer && layer.visible === false ? 0.0 : opacity;
          if (typeof imageLayer.setOpacity === "function") {
            imageLayer.setOpacity(visibleOpacity);
          } else if (typeof imageLayer.setAlpha === "function") {
            imageLayer.setAlpha(visibleOpacity);
          } else if ("opacity" in imageLayer) {
            imageLayer.opacity = visibleOpacity;
          }
          const stretch = normalizeSkyLayerStretch(layer && layer.stretch);
          const colormap = normalizeSkyLayerColormap(layer && layer.colormap);
          if (colormap && typeof imageLayer.setColormap === "function") {
            try {
              imageLayer.setColormap(colormap, { stretch });
            } catch (_err) {
            }
          } else if (typeof imageLayer.setStretch === "function") {
            try {
              imageLayer.setStretch(stretch);
            } catch (_err) {
            }
          }
          const cutMin = numericSkyLayerCut(layer && (layer.cut_min ?? layer.cutMin));
          const cutMax = numericSkyLayerCut(layer && (layer.cut_max ?? layer.cutMax));
          if (cutMin !== null && cutMax !== null && cutMax > cutMin && typeof imageLayer.setCuts === "function") {
            try {
              imageLayer.setCuts(cutMin, cutMax);
            } catch (_err) {
            }
          }
        }
        function setBaseSkyImageLayer(survey) {
          if (!aladinInstance || !survey) {
            return;
          }
          try {
            if (typeof aladinInstance.setImageSurvey === "function") {
              aladinInstance.setImageSurvey(survey);
            } else if (typeof aladinInstance.setBaseImageLayer === "function") {
              aladinInstance.setBaseImageLayer(survey);
            }
            activeImageSurvey = survey;
          } catch (_err) {
            if (typeof aladinInstance.setImageSurvey === "function") {
              try {
                aladinInstance.setImageSurvey(survey);
                activeImageSurvey = survey;
              } catch (__err) {
              }
            }
          }
        }
        function skyImageSurveyObject(survey) {
          if (!aladinInstance || !survey) {
            return null;
          }
          try {
            if (typeof aladinInstance.newImageSurvey === "function") {
              return aladinInstance.newImageSurvey(survey);
            }
          } catch (_err) {
          }
          return survey;
        }
        function setOverlaySkyImageLayer(layerName, survey, layer, expectedSignature) {
          if (!aladinInstance || !survey || !layerName) {
            return;
          }
          const attachOverlay = (overlaySurvey) => {
            if (!overlaySurvey || (expectedSignature && activeSkyLayerStackSignature !== expectedSignature)) {
              return;
            }
            try {
              if (typeof aladinInstance.setOverlayImageLayer === "function") {
                aladinInstance.setOverlayImageLayer(overlaySurvey, layerName);
                managedSkyOverlayLayerNames.add(layerName);
                applySkyImageLayerOptions(layerName, layer, false);
              }
            } catch (_err) {
            }
          };
          const overlaySurvey = skyImageSurveyObject(survey);
          if (overlaySurvey && typeof overlaySurvey.then === "function") {
            overlaySurvey.then(attachOverlay).catch(() => {});
            return;
          }
          attachOverlay(overlaySurvey);
        }
        function applySkyLayerState(data) {
          if (!data || typeof data !== "object") {
            return;
          }
          const layers = (Array.isArray(data.layers) ? data.layers : [])
            .filter((layer) => layer && String(layer.survey || layer.key || "").trim());
          const visibleLayers = layers.filter((layer) => (
            layer.visible !== false
            && Math.min(Math.max(Number(layer.opacity), 0.0), 1.0) > 0.0
          ));
          const aladinEl = document.getElementById("aladin-lite-div");
          if (!visibleLayers.length) {
            if (aladinEl) {
              aladinEl.style.opacity = "0";
            }
            removeManagedSkyOverlays();
            return;
          }
          if (aladinEl) {
            aladinEl.style.opacity = "1";
          }
          const stackSignature = visibleLayers
            .map((layer) => String(layer.key || layer.survey) + ":" + String(layer.survey || layer.key))
            .join("|");
          const baseLayer = visibleLayers[visibleLayers.length - 1];
          const baseSurvey = String(baseLayer.survey || baseLayer.key || "").trim();
          if (stackSignature !== activeSkyLayerStackSignature) {
            activeSkyLayerStackSignature = stackSignature;
            removeManagedSkyOverlays();
            setBaseSkyImageLayer(baseSurvey);
            for (let index = visibleLayers.length - 2; index >= 0; index -= 1) {
              const layer = visibleLayers[index];
              setOverlaySkyImageLayer(
                skyLayerNameFor(layer),
                String(layer.survey || layer.key || "").trim(),
                layer,
                stackSignature
              );
            }
          } else if (baseSurvey && baseSurvey !== activeImageSurvey) {
            setBaseSkyImageLayer(baseSurvey);
          }
          applySkyImageLayerOptions("base", baseLayer, true);
          for (let index = visibleLayers.length - 2; index >= 0; index -= 1) {
            const layer = visibleLayers[index];
            applySkyImageLayerOptions(skyLayerNameFor(layer), layer, false);
          }
        }
        function scheduleSkyBackgroundView(data) {
          pendingSkyBackgroundView = data;
          if (skyBackgroundViewAnimationFrame) {
            return;
          }
          const scheduleFrame = typeof window.requestAnimationFrame === "function"
            ? (callback) => window.requestAnimationFrame(callback)
            : (callback) => window.setTimeout(callback, 0);
          skyBackgroundViewAnimationFrame = scheduleFrame(() => {
            skyBackgroundViewAnimationFrame = 0;
            const latestView = pendingSkyBackgroundView;
            pendingSkyBackgroundView = null;
            applySkyBackgroundViewNow(latestView);
            if (window.parent && window.parent !== window && latestView && latestView.seq != null) {
              window.parent.postMessage({
                type: "oviz-aladin-sky-background-view-applied",
                seq: latestView.seq,
              }, "*");
            }
          });
        }
        function setHoveredClusterKey(clusterKey, emitToParent) {
          const nextKey = normalizeClusterKey(clusterKey);
          if (nextKey === hoveredClusterKey) {
            return;
          }
          hoveredClusterKey = nextKey;
          refreshCatalogHoverShapes();
          if (emitToParent) {
            postHoveredClusterKey(nextKey);
          }
        }
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
          const aladinOptions = {
            survey: ${survey},
            fov: ${fovDeg},
            target: ${target},
            cooFrame: ${cooFrame},
            expandLayersControl: false,
            showReticle: viewMode === "click" && !skyDomeCaptureOnly,
            showLayersControl: !skyDomeCaptureOnly,
            showGotoControl: !skyDomeCaptureOnly,
            showFrame: !skyDomeCaptureOnly,
            showCooGrid: false,
            showCooGridControl: false,
            showProjectionControl: !skyDomeCaptureOnly,
            showZoomControl: !skyDomeCaptureOnly,
            showFullscreenControl: !skyDomeCaptureOnly,
            showShareControl: false,
            showContextMenu: false
          };
          if (viewMode !== "click") {
            aladinOptions.projection = skyDomeBackgroundOnly ? "TAN" : "MOL";
          }
          const aladin = A.aladin("#aladin-lite-div", aladinOptions);
          aladinInstance = aladin;
          postAladinHipsFavorites();
          if (aladin && typeof aladin.setImageSurvey === "function") {
            aladin.setImageSurvey(${survey});
          }
          if (viewMode === "click" && aladin && typeof aladin.gotoRaDec === "function") {
            aladin.gotoRaDec(${ra.toFixed(8)}, ${dec.toFixed(8)});
          }
          if (viewMode !== "click" && aladin && typeof aladin.setProjection === "function") {
            aladin.setProjection(skyDomeBackgroundOnly ? "TAN" : "MOL");
          }
          if (viewMode !== "click" && aladin && typeof aladin.setFoV === "function") {
            aladin.setFoV(skyDomeBackgroundOnly ? 90.0 : 360.0);
          }
          if (viewMode === "click") {
            const beam = A.graphicOverlay({ color: ${beamColor}, lineWidth: 2, opacity: 0.95 });
            aladin.addOverlay(beam);
            beamCenters.forEach((beamDef) => {
              beam.add(A.circle(Number(beamDef.ra), Number(beamDef.dec), ${radiusDeg.toFixed(8)}));
            });
          }
          if (
            viewMode !== "click"
            && volumeOverlay
            && volumeOverlay.data_url
            && volumeOverlay.wcs
            && typeof A.image === "function"
            && typeof aladin.setOverlayImageLayer === "function"
          ) {
            const imageLayer = A.image(String(volumeOverlay.data_url), {
              name: String(volumeOverlay.name || "Selected Volume"),
              imgFormat: "png",
              wcs: volumeOverlay.wcs,
            });
            aladin.setOverlayImageLayer(imageLayer, String(volumeOverlay.name || "Selected Volume"));
          }
          payload.forEach((catDef) => {
            const baseSize = Math.max(Number(catDef.sourceSize) || 7, 3);
            const catalogShapeFn = (source) => {
              const sourceClusterKey = source && source.data ? normalizeClusterKey(source.data.clusterKey) : "";
              const highlightSize = sourceClusterKey && sourceClusterKey === hoveredClusterKey
                ? Math.max(baseSize + 3, baseSize * 1.55)
                : baseSize;
              return markerCanvas(catDef.color, highlightSize, Number(catDef.opacity));
            };
            const cat = A.catalog({
              name: catDef.name,
              color: catDef.color,
              sourceSize: baseSize,
              shape: catalogShapeFn,
              opacity: catDef.opacity
            });
            cat.__ovizShapeFn = catalogShapeFn;
            aladin.addCatalog(cat);
            catalogs.push(cat);
            const sources = [];
            (catDef.points || []).forEach((pt) => {
              sources.push(A.source(Number(pt.ra), Number(pt.dec), {
                popupTitle: pt.label || catDef.name,
                clusterKey: pt.clusterKey || "",
              }));
            });
            if (sources.length) {
              cat.addSources(sources);
            }
          });
          if (viewMode !== "click" && aladin && typeof aladin.on === "function") {
            aladin.on("objectHovered", (object) => {
              const clusterKey = object && object.data ? object.data.clusterKey : "";
              setHoveredClusterKey(clusterKey, true);
            });
          }
          window.addEventListener("message", (event) => {
            const data = event && event.data;
            if (!data || typeof data !== "object") {
              return;
            }
            if (data.type !== "oviz-parent-hover-cluster") {
              if (data.type === "oviz-sky-background-view") {
                scheduleSkyBackgroundView(data);
              }
              if (data.type === "oviz-sky-layer-state") {
                applySkyLayerState(data);
              }
              return;
            }
            setHoveredClusterKey(data.clusterKey, false);
          });
          if (statusEl) {
            statusEl.style.display = "none";
          }
          if (skyDomeBackgroundOnly && viewMode === "overview") {
            postSkyBackgroundReady();
          } else if (skyDomeCaptureEnabled && viewMode === "overview") {
            scheduleSkyDomeSnapshot(aladin, 1400);
            scheduleSkyDomeSnapshot(aladin, 5200);
          }
        }).catch((err) => {
          fail("Aladin initialization error: " + (err && err.message ? err.message : String(err)));
        });
      })();
    <\/script>
  </body>
</html>`;
      }

      function updateSkyDomeCaptureFrame() {
        if (!skySpec.enabled || !skyDomeIsEnabled() || !skyDomeFrameEl) {
          return;
        }
        if (skyDomeUsesNativeHips() || skyDomeUsesHips2Fits()) {
          skyDomeFrameEl.removeAttribute("srcdoc");
          skyDomeFrameEl.style.width = "1px";
          skyDomeFrameEl.style.height = "1px";
          skyDomeCaptureFrameSignature = skyDomeUsesHips2Fits() ? "hips2fits" : "native-hips";
          return;
        }
        const signature = JSON.stringify({
          survey: skySpec.survey || "",
          mode: skyDomeUsesAladinBackground() ? "aladin-background" : "snapshot",
          projection: skyDomeSpec.projection || "MOL",
          width: skyDomeSpec.capture_width_px || null,
          height: skyDomeSpec.capture_height_px || null,
          format: skyDomeSpec.capture_format || null,
          quality: skyDomeSpec.capture_quality || null,
        });
        if (skyDomeCaptureFrameSignature === signature) {
          return;
        }
        skyDomeBackgroundFrameReady = false;
        skyDomeBackgroundViewSignature = "";
        if (skyDomeUsesAladinBackground()) {
          skyDomeFrameEl.style.width = "100%";
          skyDomeFrameEl.style.height = "100%";
        } else {
          const captureWidth = Math.max(512, Math.min(8192, Math.round(Number(skyDomeSpec.capture_width_px) || 4096)));
          const captureHeight = Math.max(256, Math.min(4096, Math.round(Number(skyDomeSpec.capture_height_px) || 2048)));
          skyDomeFrameEl.style.width = `${captureWidth}px`;
          skyDomeFrameEl.style.height = `${captureHeight}px`;
        }
        skyDomeFrameEl.srcdoc = buildAladinSrcdoc([], [], "overview", null, true);
        skyDomeCaptureFrameSignature = signature;
      }

      function updateSkyPanel() {
        if (!skySpec.enabled) {
          return;
        }
        if (widgetModeForKey("sky") === "hidden") {
          return;
        }
        const selectionsForPanel = currentSelection ? [currentSelection] : currentSelections;
        const mode = currentSelection ? "click" : "overview";
        const payload = buildAladinCatalogPayload(selectionsForPanel, mode);
        const volumeOverlay = buildVolumeSkyImageOverlaySpec(mode);
        setSkyHoveredClusterKey("");
        lastSentSkyHoverClusterKey = null;
        skyFrameEl.srcdoc = buildAladinSrcdoc(selectionsForPanel, payload, mode, volumeOverlay);
      }

      function updateSelectionUI() {
        if (!clearSelectionButtonEl || !clickSelectToggleEl || !volumeLassoToggleEl) {
          return;
        }
        if (minimalModeEnabled) {
          clearSelectionButtonEl.disabled = true;
          clickSelectToggleEl.checked = false;
          clickSelectToggleEl.disabled = true;
          volumeLassoToggleEl.checked = false;
          volumeLassoToggleEl.disabled = true;
          return;
        }
        clearSelectionButtonEl.disabled = currentSelections.length === 0 && !currentSelection && !hasActiveLassoSelectionMask();
        clickSelectToggleEl.checked = clickSelectionEnabled;
        volumeLassoToggleEl.checked = lassoVolumeSelectionEnabled;
      }

      function onWindowMessage(event) {
        if (!skySpec.enabled) {
          return;
        }
        const fromSkyPanel = Boolean(skyFrameEl && event.source === skyFrameEl.contentWindow);
        const fromSkyDomeCapture = Boolean(skyDomeFrameEl && event.source === skyDomeFrameEl.contentWindow);
        if (!fromSkyPanel && !fromSkyDomeCapture) {
          return;
        }
        const data = event && event.data;
        if (!data || typeof data !== "object") {
          return;
        }
        if (data.type === "oviz-sky-hover-cluster") {
          if (fromSkyDomeCapture) {
            return;
          }
          setSkyHoveredClusterKey(data.clusterKey);
        } else if (data.type === "oviz-aladin-sky-background-view-applied") {
          if (!fromSkyDomeCapture) {
            return;
          }
          markSkyDomeBackgroundViewApplied(data.seq);
        } else if (data.type === "oviz-aladin-sky-background-ready") {
          if (!fromSkyDomeCapture) {
            return;
          }
          if (String(skyDomeSpec && skyDomeSpec.source || "").toLowerCase() !== "aladin") {
            return;
          }
          skyDomeBackgroundFrameReady = true;
          skyDomeSurvey = skyDomeLocalSourceName(data.survey || (skySpec && skySpec.survey) || "aladin");
          setSkyDomeSnapshotStatus("loaded", "");
          postSkyLayerStateToAladin();
          updateSkyDomeBackgroundFrame(
            (typeof performance !== "undefined" && performance.now) ? performance.now() : Date.now()
          );
        } else if (data.type === "oviz-aladin-hips-favorites") {
          if (!fromSkyDomeCapture) {
            return;
          }
          setAladinDefaultSkyLayerPresets(data.favorites);
        } else if (data.type === "oviz-aladin-sky-snapshot") {
          if (String(skyDomeSpec && skyDomeSpec.source || "").toLowerCase() !== "aladin") {
            return;
          }
          if (data.status === "ready" && data.dataUrl) {
            setSkyDomeTextureFromDataUrl(String(data.dataUrl), data.survey, data.projection);
          } else {
            setSkyDomeSnapshotStatus(
              String(data.status || "unavailable"),
              data.message ? String(data.message) : ""
            );
          }
        }
      }

      function traceStyleStateForKey(traceKey) {
        return traceStyleStateByKey[String(traceKey)] || null;
      }

      function volumeStateKeyForLayer(layer) {
        return String((layer && (layer.state_key || layer.key)) || "");
      }

      function volumeStateNameForLayer(layer) {
        if (!layer) {
          return "";
        }
        return String(layer.state_name || layer.name || volumeStateKeyForLayer(layer));
      }

      function volumeBaseNameForLayer(layer) {
        if (!layer) {
          return "";
        }
        const explicitBaseName = String(layer.base_state_name || "").trim();
        if (explicitBaseName) {
          return explicitBaseName;
        }
        const variantLabel = String(layer.variant_label || "").trim();
        const stateName = volumeStateNameForLayer(layer);
        if (variantLabel) {
          const suffix = ` (${variantLabel})`;
          if (stateName.endsWith(suffix)) {
            return stateName.slice(0, -suffix.length);
          }
        }
        return stateName;
      }

      function frameVolumeLayers(frame = currentFrame()) {
        if (!frame || !Array.isArray(frame.decorations)) {
          return [];
        }
        return frame.decorations
          .filter((decoration) => decoration && decoration.kind === "volume_layer")
          .map((decoration) => volumeLayersByKey.get(String(decoration.key)))
          .filter(Boolean);
      }

      function frameVolumeLayerForStateKey(stateKey, frame = currentFrame()) {
        const wantedKey = String(stateKey || "");
        if (!wantedKey) {
          return null;
        }
        const frameLayers = frameVolumeLayers(frame);
        const matchedLayer = frameLayers.find((layer) => volumeStateKeyForLayer(layer) === wantedKey);
        if (matchedLayer) {
          return matchedLayer;
        }
        return volumeLayers.find((layer) => volumeStateKeyForLayer(layer) === wantedKey) || null;
      }

      function volumeLayerForKey(layerKey) {
        const exactLayer = volumeLayersByKey.get(String(layerKey));
        if (exactLayer) {
          return exactLayer;
        }
        return frameVolumeLayerForStateKey(layerKey);
      }

      function volumeVariantGroupForLayer(layer) {
        return String((layer && layer.variant_group) || "").trim();
      }

      function volumeVariantLabelForLayer(layer) {
        if (!layer) {
          return "";
        }
        return String(layer.variant_label || volumeStateNameForLayer(layer) || "").trim();
      }

      function volumeVariantLayersForGroup(variantGroup) {
        const groupKey = String(variantGroup || "").trim();
        if (!groupKey) {
          return [];
        }
        const variantLayers = [];
        const seenStateKeys = new Set();
        volumeStateKeys.forEach((stateKey) => {
          const layer = volumeLayerForKey(stateKey);
          if (!layer || volumeVariantGroupForLayer(layer) !== groupKey) {
            return;
          }
          const resolvedStateKey = volumeStateKeyForLayer(layer);
          if (!resolvedStateKey || seenStateKeys.has(resolvedStateKey)) {
            return;
          }
          seenStateKeys.add(resolvedStateKey);
          variantLayers.push(layer);
        });
        variantLayers.sort((left, right) => {
          const leftOrder = Number(left.variant_order);
          const rightOrder = Number(right.variant_order);
          const leftHasOrder = Number.isFinite(leftOrder);
          const rightHasOrder = Number.isFinite(rightOrder);
          if (leftHasOrder || rightHasOrder) {
            return (leftHasOrder ? leftOrder : Number.MAX_SAFE_INTEGER)
              - (rightHasOrder ? rightOrder : Number.MAX_SAFE_INTEGER);
          }
          return volumeVariantLabelForLayer(left).localeCompare(volumeVariantLabelForLayer(right));
        });
        return variantLayers;
      }

      function firstVolumeVariantStateKey(variantGroup) {
        const variantLayers = volumeVariantLayersForGroup(variantGroup);
        if (!variantLayers.length) {
          return "";
        }
        return volumeStateKeyForLayer(variantLayers[0]);
      }

      function activeVolumeControlKey() {
        const layer = selectedVolumeLayer() || volumeLayerForKey(activeVolumeKey);
        if (!layer) {
          return "";
        }
        const variantGroup = volumeVariantGroupForLayer(layer);
        if (variantGroup) {
          return `variant:${variantGroup}`;
        }
        const stateKey = volumeStateKeyForLayer(layer);
        return stateKey ? `state:${stateKey}` : "";
      }

      function volumeControlOptions() {
        const options = [];
        const seenControlKeys = new Set();
        volumeStateKeys.forEach((stateKey) => {
          const layer = volumeLayerForKey(stateKey);
          if (!layer) {
            return;
          }
          const variantGroup = volumeVariantGroupForLayer(layer);
          if (variantGroup) {
            const controlKey = `variant:${variantGroup}`;
            if (seenControlKeys.has(controlKey)) {
              return;
            }
            seenControlKeys.add(controlKey);
            options.push({
              controlKey,
              label: volumeBaseNameForLayer(layer),
            });
            return;
          }
          const resolvedStateKey = volumeStateKeyForLayer(layer);
          const controlKey = `state:${resolvedStateKey}`;
          if (seenControlKeys.has(controlKey)) {
            return;
          }
          seenControlKeys.add(controlKey);
          options.push({
            controlKey,
            label: volumeStateNameForLayer(layer),
          });
        });
        return options;
      }

      function setExclusiveVolumeVariantSelection(stateKey) {
        const targetLayer = volumeLayerForKey(stateKey);
        const variantGroup = volumeVariantGroupForLayer(targetLayer);
        if (!variantGroup) {
          return;
        }
        const targetStateKey = String(stateKey || "");
        const stateKeysInGroup = new Set(
          volumeLayers
            .filter((layer) => volumeVariantGroupForLayer(layer) === variantGroup)
            .map((layer) => volumeStateKeyForLayer(layer))
            .filter(Boolean)
        );
        stateKeysInGroup.forEach((groupStateKey) => {
          const state = volumeStateByKey[String(groupStateKey)];
          if (!state) {
            return;
          }
          state.visible = String(groupStateKey) === targetStateKey;
        });
      }

      function volumeLegendColorForLayer(layer) {
        if (!layer) {
          return theme.text_color || theme.axis_color;
        }
        const state = volumeStateByKey[volumeStateKeyForLayer(layer)] || {};
        const option = volumeColormapOptionFor(layer, state.colormap);
        return (option && option.legend_color) || layer.legend_color || theme.text_color || theme.axis_color;
      }

      function legendColorForItem(item) {
        const itemKey = String(item.key);
        const volumeLayer = volumeLayerForKey(itemKey);
        if (volumeLayer) {
          return volumeLegendColorForLayer(volumeLayer);
        }
        const state = traceStyleStateForKey(itemKey);
        return (state && state.color) || item.default_color || item.color || theme.text_color || theme.axis_color;
      }

      function volumeSummaryTextFor(layer, state) {
        if (!layer || !state) {
          return "";
        }
        const unitText = layer.value_unit ? ` ${layer.value_unit}` : "";
        const supportText = volumeSupported
          ? (
            Number.isFinite(Number(layer.time_myr))
              ? "Rendered as a time-resolved WebGL2 ray-marched volume matched to the current animation frame."
              : (
                state.showAllTimes
                  ? (
                    layer.co_rotate_with_frame
                      ? "Rendered across all animation frames with a co-rotating local reference frame."
                      : "Rendered across all animation frames as a static WebGL2 ray-marched volume."
                  )
                  : "Rendered at t=0 only as a WebGL2 ray-marched volume."
              )
          )
          : "Volume rendering requires WebGL2 support in the browser.";
        const resampleMethod = String(layer.downsample_method || "resampled");
        return [
          `${volumeStateNameForLayer(layer)}`,
          `${Number(layer.shape.x)} x ${Number(layer.shape.y)} x ${Number(layer.shape.z)} voxels`,
          `Source: ${Number(layer.source_shape.x)} x ${Number(layer.source_shape.y)} x ${Number(layer.source_shape.z)} | ${resampleMethod} scale ${formatVolumeNumber(Number(layer.downsample_step.x))} x ${formatVolumeNumber(Number(layer.downsample_step.y))} x ${formatVolumeNumber(Number(layer.downsample_step.z))}`,
          `Samples: ${Math.round(Number(state.steps))} | alpha_coef: ${Math.round(Number(state.alphaCoef))} | stretch: ${normalizeVolumeStretch(state.stretch)}`,
          `Data range: ${formatVolumeNumber((layer.data_range || [0])[0])} to ${formatVolumeNumber((layer.data_range || [0, 1])[1])}${unitText}`,
          supportText,
        ].join("\\n");
      }

      function createLegendField(labelText, controlEl) {
        const field = document.createElement("label");
        field.className = "oviz-three-legend-field";
        const label = document.createElement("span");
        label.textContent = labelText;
        field.appendChild(label);
        field.appendChild(controlEl);
        return { field, label };
      }

      function createLegendControlRow() {
        const row = document.createElement("div");
        row.className = "oviz-three-legend-control-row";
        return row;
      }

      function base64ToUint8Array(encoded) {
        const binary = window.atob(String(encoded || ""));
        const bytes = new Uint8Array(binary.length);
        for (let i = 0; i < binary.length; i += 1) {
          bytes[i] = binary.charCodeAt(i);
        }
        return bytes;
      }

      function base64ToUint16Array(encoded) {
        const bytes = base64ToUint8Array(encoded);
        return new Uint16Array(bytes.buffer, bytes.byteOffset, Math.floor(bytes.byteLength / 2));
      }

      function formatVolumeNumber(value) {
        const num = Number(value);
        if (!Number.isFinite(num)) {
          return "";
        }
        const abs = Math.abs(num);
        if ((abs > 0.0 && abs < 1e-3) || abs >= 1e4) {
          return num.toExponential(2);
        }
        if (abs >= 100.0) {
          return num.toFixed(2);
        }
        if (abs >= 1.0) {
          return num.toFixed(4);
        }
        return num.toPrecision(3);
      }

      function formatVolumeInputNumber(value) {
        const num = Number(value);
        if (!Number.isFinite(num)) {
          return "";
        }
        const abs = Math.abs(num);
        const text = (abs > 0.0 && abs < 1e-4) || abs >= 1e6
          ? num.toExponential(8)
          : num.toPrecision(10);
        return text
          .replace(/(\.\d*?[1-9])0+(e|$)/, "$1$2")
          .replace(/\.0+(e|$)/, "$1")
          .replace(/\.(e|$)/, "$1");
      }

      function volumeDataRangeForLayer(layer) {
        const dataRange = layer && Array.isArray(layer.data_range) ? layer.data_range : [0.0, 1.0];
        const dataMin = Number(dataRange[0]);
        const dataMax = Number(dataRange[1]);
        if (!Number.isFinite(dataMin) || !Number.isFinite(dataMax) || !(dataMax > dataMin)) {
          return null;
        }
        return { min: dataMin, max: dataMax, span: dataMax - dataMin };
      }

      function volumeWindowStepForLayer(layer) {
        const dataRange = volumeDataRangeForLayer(layer);
        if (!dataRange || !(dataRange.span > 0.0)) {
          return "any";
        }
        const rawStep = dataRange.span / 500.0;
        const magnitude = Math.pow(10.0, Math.floor(Math.log10(rawStep)));
        const normalized = rawStep / magnitude;
        const niceMultiplier = normalized <= 1.0 ? 1.0 : (normalized <= 2.0 ? 2.0 : (normalized <= 5.0 ? 5.0 : 10.0));
        return formatVolumeInputNumber(niceMultiplier * magnitude);
      }

      function syncVolumeWindowInput(inputEl, value, layer) {
        if (!inputEl) {
          return;
        }
        const dataRange = volumeDataRangeForLayer(layer);
        const stepValue = volumeWindowStepForLayer(layer);
        inputEl.step = stepValue;
        if (dataRange) {
          const numericStep = Number(stepValue);
          if (Number.isFinite(numericStep) && numericStep > 0.0) {
            inputEl.min = formatVolumeInputNumber(Math.floor(dataRange.min / numericStep) * numericStep);
            inputEl.max = formatVolumeInputNumber(Math.ceil(dataRange.max / numericStep) * numericStep);
          } else {
            inputEl.min = formatVolumeInputNumber(dataRange.min);
            inputEl.max = formatVolumeInputNumber(dataRange.max);
          }
        } else {
          inputEl.removeAttribute("min");
          inputEl.removeAttribute("max");
        }
        inputEl.value = formatVolumeInputNumber(value);
      }

      function finiteNumberInputValue(inputEl) {
        if (!inputEl || String(inputEl.value || "").trim() === "") {
          return null;
        }
        const value = Number(inputEl.value);
        return Number.isFinite(value) ? value : null;
      }

      function selectedVolumeLayer() {
        if (!activeVolumeKey) {
          return null;
        }
        return frameVolumeLayerForStateKey(activeVolumeKey);
      }

      function selectedVolumeState() {
        if (!activeVolumeKey) {
          return null;
        }
        return volumeStateByKey[String(activeVolumeKey)] || null;
      }

      function clampVolumeStateForLayer(layer, state) {
        if (!layer || !state) {
          return;
        }
        const dataRange = Array.isArray(layer.data_range) ? layer.data_range : [0.0, 1.0];
        const dataMin = Number(dataRange[0]);
        const dataMax = Number(dataRange[1]);
        if (!Number.isFinite(state.vmin)) {
          state.vmin = Number((layer.default_controls || {}).vmin);
        }
        if (!Number.isFinite(state.vmax)) {
          state.vmax = Number((layer.default_controls || {}).vmax);
        }
        if (Number.isFinite(dataMin) && Number.isFinite(dataMax)) {
          state.vmin = Math.min(Math.max(state.vmin, dataMin), dataMax);
          state.vmax = Math.min(Math.max(state.vmax, dataMin), dataMax);
        }
        if (!(state.vmax > state.vmin)) {
          if (Number.isFinite(dataMax) && dataMax > state.vmin) {
            state.vmax = dataMax;
          } else {
            state.vmax = state.vmin + 1e-6;
          }
        }
        state.opacity = Math.min(Math.max(Number(state.opacity), 0.0), 1.0);
        state.steps = Math.round(Math.min(Math.max(Number(state.steps), 24.0), 768.0));
        if (!Number.isFinite(state.steps) || state.steps < 24) {
          state.steps = Number((layer.default_controls || {}).steps || 100);
        }
        state.alphaCoef = Math.min(Math.max(Number(state.alphaCoef), 1.0), 200.0);
        if (!Number.isFinite(state.alphaCoef)) {
          state.alphaCoef = Number((layer.default_controls || {}).alpha_coef || 50.0);
        }
        state.gradientStep = Math.min(Math.max(Number(state.gradientStep), 1e-4), 0.05);
        if (!Number.isFinite(state.gradientStep)) {
          state.gradientStep = Number((layer.default_controls || {}).gradient_step || 0.005);
        }
        state.stretch = normalizeVolumeStretch(
          state.stretch !== undefined ? state.stretch : (layer.default_controls || {}).stretch
        );
      }

      function volumeColormapOptionFor(layer, colormapName) {
        const options = (layer && layer.colormap_options) || [];
        const requested = String(colormapName || "").trim().toLowerCase();
        for (const option of options) {
          if (String(option.name || "").trim().toLowerCase() === requested) {
            return option;
          }
        }
        return options.length ? options[0] : null;
      }

      function volumeStretchOptions() {
        return [
          { value: "linear", label: "Linear" },
          { value: "log10", label: "log10" },
          { value: "asinh", label: "asinh" },
        ];
      }

      function normalizeVolumeStretch(stretchName) {
        const requested = String(stretchName || "").trim().toLowerCase();
        const option = volumeStretchOptions().find((item) => item.value === requested);
        return option ? option.value : "linear";
      }

      function volumeStretchModeValue(stretchName) {
        const stretch = normalizeVolumeStretch(stretchName);
        if (stretch === "log10") {
          return 1.0;
        }
        if (stretch === "asinh") {
          return 2.0;
        }
        return 0.0;
      }

      function volumeScalarArrayFor(layer) {
        const layerKey = String(layer.key);
        if (volumeScalarDataCache.has(layerKey)) {
          return volumeScalarDataCache.get(layerKey);
        }
        const encoding = String(layer.data_encoding || "uint16_le");
        let data = null;
        if (encoding === "uint16_le") {
          data = base64ToUint16Array(layer.data_b64 || "");
        } else {
          data = base64ToUint8Array(layer.data_b64 || "");
        }
        volumeScalarDataCache.set(layerKey, data);
        return data;
      }

      function volumeSkyScalarArrayFor(layer) {
        const layerKey = `${String(layer.key)}::sky-overlay`;
        if (volumeScalarDataCache.has(layerKey)) {
          return volumeScalarDataCache.get(layerKey);
        }
        if (!layer || !layer.sky_overlay_data_b64) {
          const fallback = volumeScalarArrayFor(layer);
          volumeScalarDataCache.set(layerKey, fallback);
          return fallback;
        }
        const encoding = String(layer.sky_overlay_data_encoding || layer.data_encoding || "uint8");
        if (encoding === "png_atlas_uint8") {
          if (!volumeScalarDataPendingCache.has(layerKey)) {
            const shape = layer.sky_overlay_shape || layer.shape || {};
            const nx = Math.max(1, Math.round(Number(shape.x) || 0));
            const ny = Math.max(1, Math.round(Number(shape.y) || 0));
            const nz = Math.max(1, Math.round(Number(shape.z) || 0));
            const tiles = layer.sky_overlay_atlas_tiles || {};
            const tileCols = Math.max(1, Math.round(Number(tiles.x) || Math.ceil(Math.sqrt(nz))));
            const tileRows = Math.max(1, Math.round(Number(tiles.y) || Math.ceil(nz / tileCols)));
            const pending = new Promise((resolve) => {
              const image = new Image();
              image.onload = () => {
                const atlasWidth = Math.max(1, Number(image.naturalWidth || image.width || 0));
                const atlasHeight = Math.max(1, Number(image.naturalHeight || image.height || 0));
                const canvasEl = document.createElement("canvas");
                canvasEl.width = atlasWidth;
                canvasEl.height = atlasHeight;
                const ctx = canvasEl.getContext("2d");
                if (!ctx) {
                  volumeScalarDataPendingCache.delete(layerKey);
                  resolve(null);
                  return;
                }
                ctx.drawImage(image, 0, 0, atlasWidth, atlasHeight);
                const imageData = ctx.getImageData(0, 0, atlasWidth, atlasHeight);
                const rgba = imageData.data || [];
                const values = new Uint8Array(nx * ny * nz);
                for (let zIndex = 0; zIndex < nz; zIndex += 1) {
                  const tileRow = Math.floor(zIndex / tileCols);
                  const tileCol = zIndex % tileCols;
                  if (tileRow >= tileRows) {
                    break;
                  }
                  const xOffset = tileCol * nx;
                  const yOffset = tileRow * ny;
                  for (let yIndex = 0; yIndex < ny; yIndex += 1) {
                    for (let xIndex = 0; xIndex < nx; xIndex += 1) {
                      const atlasIndex = (((yOffset + yIndex) * atlasWidth) + (xOffset + xIndex)) * 4;
                      const voxelIndex = (zIndex * ny * nx) + (yIndex * nx) + xIndex;
                      values[voxelIndex] = Number(rgba[atlasIndex] || 0);
                    }
                  }
                }
                volumeScalarDataCache.set(layerKey, values);
                volumeScalarDataPendingCache.delete(layerKey);
                if (skySpec.enabled) {
                  updateSkyPanel();
                }
                resolve(values);
              };
              image.onerror = () => {
                volumeScalarDataPendingCache.delete(layerKey);
                resolve(null);
              };
              image.src = `data:image/png;base64,${String(layer.sky_overlay_data_b64 || "")}`;
            });
            volumeScalarDataPendingCache.set(layerKey, pending);
          }
          return null;
        }
        let data = null;
        if (encoding === "uint16_le") {
          data = base64ToUint16Array(layer.sky_overlay_data_b64 || "");
        } else {
          data = base64ToUint8Array(layer.sky_overlay_data_b64 || "");
        }
        volumeScalarDataCache.set(layerKey, data);
        return data;
      }

      function volumeColorTextureFor(option) {
        const optionKey = String(option.name || "volume-colormap");
        if (volumeColorTextureCache.has(optionKey)) {
          return volumeColorTextureCache.get(optionKey);
        }
        const bytes = base64ToUint8Array(option.lut_b64 || "");
        const width = Math.max(1, Math.floor(bytes.length / 4));
        const texture = new THREE.DataTexture(bytes, width, 1, THREE.RGBAFormat);
        texture.minFilter = THREE.NearestFilter;
        texture.magFilter = THREE.NearestFilter;
        texture.wrapS = THREE.ClampToEdgeWrapping;
        texture.wrapT = THREE.ClampToEdgeWrapping;
        texture.unpackAlignment = 1;
        texture.generateMipmaps = false;
        texture.needsUpdate = true;
        if ("colorSpace" in texture && THREE.SRGBColorSpace) {
          texture.colorSpace = THREE.SRGBColorSpace;
        } else if ("encoding" in texture && THREE.sRGBEncoding) {
          texture.encoding = THREE.sRGBEncoding;
        }
        volumeColorTextureCache.set(optionKey, texture);
        return texture;
      }

      function volumeColorBytesForOption(option) {
        const optionKey = String((option && option.name) || "volume-colormap");
        if (volumeColorBytesCache.has(optionKey)) {
          return volumeColorBytesCache.get(optionKey);
        }
        const bytes = base64ToUint8Array((option && option.lut_b64) || "");
        volumeColorBytesCache.set(optionKey, bytes);
        return bytes;
      }

      function volumeTextureFor(layer) {
        const layerKey = String(layer.key);
        if (volumeTextureCache.has(layerKey)) {
          return volumeTextureCache.get(layerKey);
        }
        const data = volumeScalarArrayFor(layer);
        const shape = layer.shape || {};
        const nx = Math.max(1, Number(shape.x || 1));
        const ny = Math.max(1, Number(shape.y || 1));
        const nz = Math.max(1, Number(shape.z || 1));
        const VolumeTextureCtor = THREE.Data3DTexture || THREE.DataTexture3D;
        if (!VolumeTextureCtor) {
          throw new Error("Three.js volume textures are unavailable in this browser build.");
        }
        const texture = new VolumeTextureCtor(data, nx, ny, nz);
        texture.format = THREE.RedFormat;
        texture.type = THREE.UnsignedByteType;
        texture.minFilter = layer.interpolation === false ? THREE.NearestFilter : THREE.LinearFilter;
        texture.magFilter = layer.interpolation === false ? THREE.NearestFilter : THREE.LinearFilter;
        texture.wrapS = THREE.ClampToEdgeWrapping;
        texture.wrapT = THREE.ClampToEdgeWrapping;
        texture.wrapR = THREE.ClampToEdgeWrapping;
        texture.unpackAlignment = 1;
        texture.generateMipmaps = false;
        texture.needsUpdate = true;
        const volumeTexture = { texture, nx, ny, nz };
        volumeTextureCache.set(layerKey, volumeTexture);
        return volumeTexture;
      }

      function volumeJitterTextureFor() {
        if (volumeJitterTexture) {
          return volumeJitterTexture;
        }
        const bytes = new Uint8Array(64 * 64);
        for (let i = 0; i < bytes.length; i += 1) {
          bytes[i] = Math.floor(Math.random() * 256.0);
        }
        volumeJitterTexture = new THREE.DataTexture(bytes, 64, 64, THREE.RedFormat, THREE.UnsignedByteType);
        volumeJitterTexture.minFilter = THREE.LinearFilter;
        volumeJitterTexture.magFilter = THREE.LinearFilter;
        volumeJitterTexture.wrapS = THREE.MirroredRepeatWrapping;
        volumeJitterTexture.wrapT = THREE.MirroredRepeatWrapping;
        volumeJitterTexture.generateMipmaps = false;
        volumeJitterTexture.unpackAlignment = 1;
        volumeJitterTexture.needsUpdate = true;
        return volumeJitterTexture;
      }

      function normalizedVolumeWindowFor(layer, state) {
        const dataRange = Array.isArray(layer.data_range) ? layer.data_range : [0.0, 1.0];
        const dataMin = Number(dataRange[0]);
        const dataMax = Number(dataRange[1]);
        if (!Number.isFinite(dataMin) || !Number.isFinite(dataMax) || !(dataMax > dataMin)) {
          return { low: 0.0, high: 1.0 };
        }
        const span = dataMax - dataMin;
        const low = Math.min(Math.max((Number(state.vmin) - dataMin) / span, 0.0), 1.0);
        let high = Math.min(Math.max((Number(state.vmax) - dataMin) / span, 0.0), 1.0);
        if (!(high > low)) {
          high = Math.min(1.0, low + 1e-6);
        }
        return { low, high };
      }

      function volumeLayerTimeMyr(layer) {
        if (!layer) {
          return null;
        }
        const rawTime = layer.time_myr;
        if (rawTime === null || rawTime === undefined || rawTime === "" || rawTime === false) {
          return null;
        }
        const timeValue = Number(rawTime);
        return Number.isFinite(timeValue) ? timeValue : null;
      }

      function volumeSupportsShowAllTimes(layer) {
        return Boolean(
          layer
          && layer.supports_show_all_times
          && volumeLayerTimeMyr(layer) === null
        );
      }

      function volumeVisibleForFrame(layer, state, frame = currentFrame()) {
        if (!layer || !state || state.visible === false) {
          return false;
        }
        const frameTime = frame ? Number(frame.time) : 0.0;
        const layerTime = volumeLayerTimeMyr(layer);
        if (layerTime !== null) {
          return approximatelyZero(frameTime - layerTime);
        }
        if (state.showAllTimes && volumeSupportsShowAllTimes(layer)) {
          return true;
        }
        if (layer.only_at_t0 === false) {
          return true;
        }
        return approximatelyZero(frameTime);
      }

      function volumeRotationAngleForFrame(layer, state, frame = currentFrame()) {
        if (
          !layer
          || !state
          || !state.showAllTimes
          || !layer.co_rotate_with_frame
          || !volumeSupportsShowAllTimes(layer)
          || !Number.isFinite(volumeCoRotationRateRadPerMyr)
        ) {
          return 0.0;
        }
        const frameTime = frame ? Number(frame.time) : 0.0;
        if (!Number.isFinite(frameTime)) {
          return 0.0;
        }
        const referenceTime = Number.isFinite(Number(layer.reference_time_myr))
          ? Number(layer.reference_time_myr)
          : 0.0;
        return volumeCoRotationRateRadPerMyr * (frameTime - referenceTime);
      }

      function volumeQuaternionForZRotation(angleRad) {
        const halfAngle = 0.5 * (Number.isFinite(angleRad) ? angleRad : 0.0);
        return new THREE.Vector4(0.0, 0.0, Math.sin(halfAngle), Math.cos(halfAngle));
      }

      function createVolumeRuntime(layer) {
        if (!volumeSupported) {
          return null;
        }
        const state = volumeStateByKey[volumeStateKeyForLayer(layer)];
        if (!state) {
          return null;
        }
        clampVolumeStateForLayer(layer, state);
        const bounds = layer.bounds || {};
        const xBounds = Array.isArray(bounds.x) ? bounds.x : [-0.5, 0.5];
        const yBounds = Array.isArray(bounds.y) ? bounds.y : [-0.5, 0.5];
        const zBounds = Array.isArray(bounds.z) ? bounds.z : [-0.5, 0.5];
        const sizeX = Math.max(1e-6, Number(xBounds[1]) - Number(xBounds[0]));
        const sizeY = Math.max(1e-6, Number(yBounds[1]) - Number(yBounds[0]));
        const sizeZ = Math.max(1e-6, Number(zBounds[1]) - Number(zBounds[0]));
        const centerX = 0.5 * (Number(xBounds[0]) + Number(xBounds[1]));
        const centerY = 0.5 * (Number(yBounds[0]) + Number(yBounds[1]));
        const centerZ = 0.5 * (Number(zBounds[0]) + Number(zBounds[1]));
        const option = volumeColormapOptionFor(layer, state.colormap);
        if (!option) {
          return null;
        }
        const volumeTexture = volumeTextureFor(layer);
        const windowState = normalizedVolumeWindowFor(layer, state);
        const uniforms = THREE.UniformsUtils.merge([
          THREE.UniformsLib.lights,
          {
            volumeTexture: { value: volumeTexture.texture },
            colormap: { value: volumeColorTextureFor(option) },
            jitterTexture: { value: volumeJitterTextureFor() },
            low: { value: Number(windowState.low) },
            high: { value: Number(windowState.high) },
            opacity: { value: Number(state.opacity) },
            samples: { value: Number(state.steps) },
            alpha_coef: { value: Number(state.alphaCoef) },
            gradient_step: { value: Number(state.gradientStep) },
            stretch_mode: { value: volumeStretchModeValue(state.stretch) },
            scale: { value: new THREE.Vector4(sizeX, sizeY, sizeZ, 1.0) },
            translation: { value: new THREE.Vector4(centerX, centerY, centerZ, 1.0) },
            rotation: { value: new THREE.Vector4(0.0, 0.0, 0.0, 1.0) },
            useSelectionPolygon: { value: false },
            selectionViewProjectionMatrix: { value: new THREE.Matrix4() },
            selectionMaskTexture: { value: null },
            selectionDimOutside: { value: 1.0 },
          },
        ]);
        applyLassoSelectionMaskUniforms(uniforms, activeVolumeLassoSelectionMask());

        const material = new THREE.ShaderMaterial({
          uniforms,
          vertexShader: VOLUME_VERTEX_SHADER,
          fragmentShader: VOLUME_FRAGMENT_SHADER,
          transparent: true,
          side: THREE.BackSide,
          depthTest: true,
          depthWrite: false,
          lights: true,
        });

        const geometry = new THREE.BoxBufferGeometry(1, 1, 1);
        const mesh = new THREE.Mesh(geometry, material);
        mesh.position.set(centerX, centerY, centerZ);
        mesh.scale.set(sizeX, sizeY, sizeZ);
        mesh.renderOrder = -30;
        mesh.frustumCulled = false;
        const runtime = { mesh, material, geometry, layer };
        applyVolumeStateToRuntime(layer, runtime);
        return runtime;
      }

      function applyVolumeStateToRuntime(layer, runtime, frame = currentFrame()) {
        if (!runtime || !runtime.material || !runtime.mesh) {
          return;
        }
        const state = volumeStateByKey[volumeStateKeyForLayer(layer)] || {};
        const option = volumeColormapOptionFor(layer, state.colormap);
        const windowState = normalizedVolumeWindowFor(layer, state);
        runtime.mesh.visible = volumeVisibleForFrame(layer, state, frame);
        runtime.material.uniforms.low.value = Number(windowState.low);
        runtime.material.uniforms.high.value = Number(windowState.high);
        runtime.material.uniforms.opacity.value = Number(state.opacity);
        runtime.material.uniforms.samples.value = Number(state.steps);
        runtime.material.uniforms.alpha_coef.value = Number(state.alphaCoef);
        runtime.material.uniforms.gradient_step.value = Number(state.gradientStep);
        runtime.material.uniforms.stretch_mode.value = volumeStretchModeValue(state.stretch);
        runtime.material.uniforms.rotation.value.copy(
          volumeQuaternionForZRotation(volumeRotationAngleForFrame(layer, state, frame))
        );
        applyLassoSelectionMaskUniforms(runtime.material.uniforms, activeVolumeLassoSelectionMask());
        if (option) {
          runtime.material.uniforms.colormap.value = volumeColorTextureFor(option);
        }
      }

      function renderVolumeControls() {
        if (
          !volumePanelEl
          || !volumeSelectEl
          || !volumeVisibleEl
          || !volumeColormapEl
          || !volumeStretchEl
          || !volumeVMinEl
          || !volumeVMaxEl
          || !volumeOpacityEl
          || !volumeAlphaEl
          || !volumeStepsEl
          || !volumeOpacityLabelEl
          || !volumeAlphaLabelEl
          || !volumeStepsLabelEl
          || !volumeSummaryEl
        ) {
          return;
        }
        const enabled = volumeStateKeys.length > 0;
        volumePanelEl.dataset.enabled = enabled ? "true" : "false";
        if (!enabled) {
          return;
        }

        if (!activeVolumeKey || !Object.prototype.hasOwnProperty.call(volumeStateByKey, String(activeVolumeKey))) {
          activeVolumeKey = String(volumeStateKeys[0]);
        }

        const layer = selectedVolumeLayer();
        const state = selectedVolumeState();
        if (!state) {
          return;
        }
        if (layer) {
          clampVolumeStateForLayer(layer, state);
        }

        volumeSelectEl.innerHTML = "";
        const controlOptions = volumeControlOptions();
        const selectedControlKey = activeVolumeControlKey();
        controlOptions.forEach((option) => {
          const optionEl = document.createElement("option");
          optionEl.value = String(option.controlKey);
          optionEl.textContent = String(option.label);
          if (String(option.controlKey) === String(selectedControlKey)) {
            optionEl.selected = true;
          }
          volumeSelectEl.appendChild(optionEl);
        });
        volumeSelectEl.disabled = !volumeSupported || controlOptions.length <= 1;

        const controlLayer = layer || volumeLayerForKey(activeVolumeKey);
        if (!controlLayer) {
          return;
        }
        const controlVariantGroup = volumeVariantGroupForLayer(controlLayer);
        const smoothingLayers = controlVariantGroup
          ? volumeVariantLayersForGroup(controlVariantGroup)
          : [];
        const showSmoothingControl = smoothingLayers.length > 1;
        if (volumeSmoothingFieldEl) {
          volumeSmoothingFieldEl.style.display = showSmoothingControl ? "" : "none";
        }
        if (volumeSmoothingEl) {
          volumeSmoothingEl.innerHTML = "";
          smoothingLayers.forEach((variantLayer) => {
            const optionEl = document.createElement("option");
            optionEl.value = volumeStateKeyForLayer(variantLayer);
            optionEl.textContent = volumeVariantLabelForLayer(variantLayer);
            if (String(optionEl.value) === String(activeVolumeKey)) {
              optionEl.selected = true;
            }
            volumeSmoothingEl.appendChild(optionEl);
          });
          volumeSmoothingEl.disabled = !volumeSupported || !showSmoothingControl;
        }
        volumeVisibleEl.checked = state.visible !== false;
        volumeColormapEl.innerHTML = "";
        ((controlLayer.colormap_options || [])).forEach((option) => {
          const optionEl = document.createElement("option");
          optionEl.value = String(option.name);
          optionEl.textContent = String(option.label || option.name);
          if (String(option.name) === String(state.colormap)) {
            optionEl.selected = true;
          }
          volumeColormapEl.appendChild(optionEl);
        });
        volumeColormapEl.value = String(state.colormap);
        volumeStretchEl.innerHTML = "";
        volumeStretchOptions().forEach((option) => {
          const optionEl = document.createElement("option");
          optionEl.value = String(option.value);
          optionEl.textContent = String(option.label);
          if (String(option.value) === String(state.stretch)) {
            optionEl.selected = true;
          }
          volumeStretchEl.appendChild(optionEl);
        });
        volumeStretchEl.value = String(state.stretch);

        syncVolumeWindowInput(volumeVMinEl, state.vmin, controlLayer);
        syncVolumeWindowInput(volumeVMaxEl, state.vmax, controlLayer);
        volumeOpacityEl.value = String(state.opacity);
        volumeAlphaEl.value = String(state.alphaCoef);
        volumeStepsEl.value = String(state.steps);
        volumeOpacityLabelEl.textContent = `Opacity (${Number(state.opacity).toFixed(2)})`;
        volumeAlphaLabelEl.textContent = `Alpha coef (${Math.round(Number(state.alphaCoef))})`;
        volumeStepsLabelEl.textContent = `Samples (${Math.round(Number(state.steps))})`;
        volumeVisibleEl.disabled = !volumeSupported;
        volumeColormapEl.disabled = !volumeSupported;
        volumeStretchEl.disabled = !volumeSupported;
        volumeVMinEl.disabled = !volumeSupported;
        volumeVMaxEl.disabled = !volumeSupported;
        volumeOpacityEl.disabled = !volumeSupported;
        volumeAlphaEl.disabled = !volumeSupported;
        volumeStepsEl.disabled = !volumeSupported;

        volumeSummaryEl.textContent = volumeSummaryTextFor(controlLayer, state);
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

      function smoothstep(edge0, edge1, value) {
        const width = Math.max(Number(edge1) - Number(edge0), 1e-9);
        const t = clampRange((Number(value) - Number(edge0)) / width, 0.0, 1.0);
        return t * t * (3.0 - 2.0 * t);
      }

      function writeStarPsfTexture(canvasEl, alphaForSample) {
        const ctx = canvasEl.getContext("2d");
        const width = canvasEl.width;
        const height = canvasEl.height;
        const cx = (width - 1) * 0.5;
        const cy = (height - 1) * 0.5;
        const maxRadius = Math.max(Math.min(width, height) * 0.5, 1.0);
        const imageData = ctx.createImageData(width, height);
        const data = imageData.data;
        for (let y = 0; y < height; y += 1) {
          const dy = (y - cy) / maxRadius;
          for (let x = 0; x < width; x += 1) {
            const dx = (x - cx) / maxRadius;
            const radius = Math.sqrt(dx * dx + dy * dy);
            const edgeTaper = 1.0 - smoothstep(0.78, 1.0, radius);
            const angle = Math.atan2(dy, dx);
            const alpha = clampRange(Number(alphaForSample(radius, dx, dy, angle)) * edgeTaper, 0.0, 1.0);
            const idx = (y * width + x) * 4;
            data[idx] = 255;
            data[idx + 1] = 255;
            data[idx + 2] = 255;
            data[idx + 3] = Math.round(alpha * 255.0);
          }
        }
        ctx.putImageData(imageData, 0, 0);
      }

      function makeSpriteTextureFromCanvas(canvasEl) {
        const texture = new THREE.CanvasTexture(canvasEl);
        texture.colorSpace = THREE.SRGBColorSpace;
        texture.generateMipmaps = true;
        texture.minFilter = THREE.LinearMipmapLinearFilter;
        texture.magFilter = THREE.LinearFilter;
        texture.needsUpdate = true;
        return texture;
      }

      function starGlowTextureFor(kind = "halo") {
        const cacheKey = String(kind || "halo");
        if (starGlowTextureCache.has(cacheKey)) {
          return starGlowTextureCache.get(cacheKey);
        }

        const canvasEl = document.createElement("canvas");
        canvasEl.width = 1024;
        canvasEl.height = 1024;

        writeStarPsfTexture(canvasEl, (radius) => {
          const airyCore = 0.82 * Math.exp(-0.5 * (radius / 0.024) ** 2.0);
          const seeingHalo = 0.34 * ((1.0 + (radius / 0.075) ** 2.0) ** -2.15);
          const aureole = 0.055 * ((1.0 + (radius / 0.25) ** 2.0) ** -2.45);
          return airyCore + seeingHalo + aureole;
        });

        const texture = makeSpriteTextureFromCanvas(canvasEl);
        starGlowTextureCache.set(cacheKey, texture);
        return texture;
      }

      function starGlowMaterialFor(color, opacity) {
        const cacheKey = [color ?? "#ffffff", Number(opacity ?? 0.0).toFixed(4)].join("|");
        if (starGlowMaterialCache.has(cacheKey)) {
          return starGlowMaterialCache.get(cacheKey);
        }
        const material = new THREE.SpriteMaterial({
          map: starGlowTextureFor("halo"),
          color: color ?? "#ffffff",
          transparent: true,
          opacity: opacity ?? 0.0,
          depthWrite: false,
          depthTest: true,
          blending: THREE.AdditiveBlending,
        });
        material.userData = { cached: true, glow: true };
        starGlowMaterialCache.set(cacheKey, material);
        return material;
      }

      function starCoreTextureFor(kind = "stellar_core") {
        const cacheKey = String(kind || "stellar_core");
        if (starCoreTextureCache.has(cacheKey)) {
          return starCoreTextureCache.get(cacheKey);
        }

        const canvasEl = document.createElement("canvas");
        canvasEl.width = 1024;
        canvasEl.height = 1024;
        writeStarPsfTexture(canvasEl, (radius) => {
          const saturatedCore = 1.00 * Math.exp(-0.5 * (radius / 0.020) ** 2.0);
          const seeingDisk = 0.28 * ((1.0 + (radius / 0.052) ** 2.0) ** -2.65);
          return saturatedCore + seeingDisk;
        });

        const texture = makeSpriteTextureFromCanvas(canvasEl);
        starCoreTextureCache.set(cacheKey, texture);
        return texture;
      }

      function starCoreMaterialFor(color, opacity) {
        const cacheKey = [color ?? "#ffffff", Number(opacity ?? 0.0).toFixed(4)].join("|");
        if (starCoreMaterialCache.has(cacheKey)) {
          return starCoreMaterialCache.get(cacheKey);
        }
        const material = new THREE.SpriteMaterial({
          map: starCoreTextureFor("stellar_core"),
          color: color ?? "#ffffff",
          transparent: true,
          opacity: opacity ?? 0.0,
          depthWrite: false,
          depthTest: true,
          blending: THREE.AdditiveBlending,
        });
        material.userData = { cached: true, glow: true, core: true };
        starCoreMaterialCache.set(cacheKey, material);
        return material;
      }

      // Glow keeps its old apparent radius while the regular marker stays on the smaller baseline.
      function glowScaleForPoint(basePointScale, position, glowStrength) {
        const strength = Math.max(Number(glowStrength) || 0.0, 0.0);
        return Math.max(Number(basePointScale) || 0.0, 0.0) * (6.45 + 1.26 * strength);
      }

      function starCoreScaleForPoint(basePointScale, position, glowStrength) {
        const strength = Math.max(Number(glowStrength) || 0.0, 0.0);
        return Math.max(Number(basePointScale) || 0.0, 0.0) * (3.00 + 0.24 * strength);
      }

      function registerCameraResponsivePointSprite(sprite, scaleKind, position, pointScaleValue, selectionKey) {
        const entry = {
          sprite,
          scaleKind,
          position: position instanceof THREE.Vector3 ? position.clone() : new THREE.Vector3(
            Number(position && position.x) || 0.0,
            Number(position && position.y) || 0.0,
            Number(position && position.z) || 0.0,
          ),
          pointScale: Math.max(Number(pointScaleValue) || 0.0, 0.0),
          baseScale: Math.max(Number(pointScaleValue) || 0.0, 0.0),
          selectionKey: selectionKey ? String(selectionKey) : "",
        };
        cameraResponsivePointEntries.push(entry);
        if (entry.selectionKey) {
          if (!selectionSpriteEntriesByKey.has(entry.selectionKey)) {
            selectionSpriteEntriesByKey.set(entry.selectionKey, []);
          }
          selectionSpriteEntriesByKey.get(entry.selectionKey).push(entry);
        }
        return entry;
      }

      function cameraResponsivePointBaseScale(entry) {
        if (!entry || !entry.sprite) {
          return 0.0;
        }
        if (entry.scaleKind === "glow") {
          return glowScaleForPoint(entry.pointScale, entry.position, globalPointGlowStrength);
        }
        if (entry.scaleKind === "core") {
          return starCoreScaleForPoint(entry.pointScale, entry.position, globalPointGlowStrength);
        }
        return Math.max(Number(entry.pointScale) || 0.0, 0.0);
      }

      function updateCameraResponsivePointSprites() {
        if (!cameraResponsivePointEntries.length) {
          return;
        }
        const activeKeys = activeHoveredClusterKeys();
        cameraResponsivePointEntries.forEach((entry) => {
          if (!entry || !entry.sprite) {
            return;
          }
          const baseScale = cameraResponsivePointBaseScale(entry);
          entry.baseScale = baseScale;
          if (entry.sprite.userData && typeof entry.sprite.userData === "object") {
            entry.sprite.userData.baseScale = baseScale;
          }
          const isActive = entry.selectionKey && activeKeys.has(entry.selectionKey);
          const scale = baseScale * (isActive ? 1.45 : 1.0);
          entry.sprite.scale.set(scale, scale, 1.0);
        });
      }

      function registerCameraResponsiveImagePlane(mesh, options = {}) {
        if (!mesh || !mesh.material) {
          return null;
        }
        const entry = {
          mesh,
          baseOpacity: clamp01(Number(options.baseOpacity) || 0.0),
          hideBelowScaleBarPc: Math.max(Number(options.hideBelowScaleBarPc) || 0.0, 0.0),
          fadeStartScaleBarPc: Math.max(Number(options.fadeStartScaleBarPc) || 0.0, 0.0),
        };
        cameraResponsiveImagePlaneEntries.push(entry);
        return entry;
      }

      function scaleBarLengthPcForCurrentView() {
        const canvasHeight = Math.max(canvas.clientHeight || root.clientHeight || 0, 1);
        const distance = cameraViewMode === "earth" && Number.isFinite(earthViewFocusDistance) && earthViewFocusDistance > 0.0
          ? earthViewFocusDistance
          : Math.max(camera.position.distanceTo(controls.target), 1e-6);
        const worldPerPixel = (2.0 * distance * Math.tan(THREE.MathUtils.degToRad(camera.fov * 0.5))) / canvasHeight;
        return worldPerPixel * 120.0;
      }

      function updateCameraResponsiveImagePlanes() {
        if (!cameraResponsiveImagePlaneEntries.length) {
          return;
        }
        cameraResponsiveImagePlaneEntries.forEach((entry) => {
          const mesh = entry && entry.mesh;
          const material = mesh && mesh.material;
          if (!mesh || !material) {
            return;
          }
          const hideBelow = Math.max(Number(entry.hideBelowScaleBarPc) || 0.0, 0.0);
          const fadeStart = Math.max(Number(entry.fadeStartScaleBarPc) || hideBelow, hideBelow);
          let fadeFactor = 1.0;
          if (cameraViewMode === "earth") {
            fadeFactor = 0.0;
          }
          if (cameraViewMode !== "earth" && hideBelow > 0.0 && Number.isFinite(currentScaleBarLengthPc)) {
            if (currentScaleBarLengthPc <= hideBelow) {
              fadeFactor = 0.0;
            } else if (fadeStart > hideBelow && currentScaleBarLengthPc < fadeStart) {
              fadeFactor = (currentScaleBarLengthPc - hideBelow) / Math.max(fadeStart - hideBelow, 1e-6);
            }
          }
          const effectiveOpacity = clamp01(Number(entry.baseOpacity) || 0.0) * clamp01(fadeFactor);
          material.opacity = effectiveOpacity;
          mesh.visible = effectiveOpacity > 1e-4;
        });
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
        if (options.screenStable) {
          sprite.userData.screenStable = {
            pixelHeight: Math.max(Number(options.screenPixels) || 12.0, 6.0),
            scaleMultiplier: Math.max(Number(options.screenScaleMultiplier) || 1.0, 0.1),
          };
        }
        return sprite;
      }

      function updateScreenStableTextSprite(sprite) {
        if (!sprite || !sprite.material || !sprite.material.map) {
          return;
        }
        const settings = sprite.userData ? sprite.userData.screenStable : null;
        if (!settings) {
          return;
        }
        const texture = sprite.material.map;
        const aspect = texture.userData.width / Math.max(texture.userData.height, 1);
        const pixelHeight = Math.max(Number(settings.pixelHeight) || 12.0, 6.0);
        const scaleMultiplier = Math.max(Number(settings.scaleMultiplier) || 1.0, 0.1);
        let worldHeight = 0.0;

        if (camera.isPerspectiveCamera) {
          const worldPosition = new THREE.Vector3();
          sprite.getWorldPosition(worldPosition);
          const distance = Math.max(camera.position.distanceTo(worldPosition), 1e-6);
          const fovRad = THREE.MathUtils.degToRad(camera.fov);
          worldHeight = 2.0 * Math.tan(fovRad * 0.5) * distance * (pixelHeight / Math.max(renderer.domElement.clientHeight, 1));
        } else if (camera.isOrthographicCamera) {
          worldHeight = ((camera.top - camera.bottom) / Math.max(camera.zoom, 1e-6)) * (pixelHeight / Math.max(renderer.domElement.clientHeight, 1));
        }

        if (!(worldHeight > 0.0) || !Number.isFinite(worldHeight)) {
          return;
        }
        const height = worldHeight * scaleMultiplier;
        sprite.scale.set(height * aspect, height, 1.0);
      }

      function updateScreenStableTextSprites() {
        for (const sprite of screenStableTextSprites) {
          updateScreenStableTextSprite(sprite);
        }
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

      function imagePlaneSpecForDecoration(decoration) {
        const key = String((decoration && decoration.key) || "");
        if (!key) {
          return null;
        }
        return imagePlaneSpecByKey.get(key) || null;
      }

      function imagePlaneTextureFor(imagePlaneSpec) {
        const dataUrl = String(imagePlaneSpec && imagePlaneSpec.image_data_url ? imagePlaneSpec.image_data_url : "");
        if (!dataUrl) {
          return null;
        }
        if (imagePlaneTextureCache.has(dataUrl)) {
          return imagePlaneTextureCache.get(dataUrl);
        }
        const texture = new THREE.TextureLoader().load(dataUrl);
        texture.minFilter = THREE.LinearFilter;
        texture.magFilter = THREE.LinearFilter;
        texture.generateMipmaps = true;
        if ("colorSpace" in texture && THREE.SRGBColorSpace) {
          texture.colorSpace = THREE.SRGBColorSpace;
        } else if ("encoding" in texture && THREE.sRGBEncoding) {
          texture.encoding = THREE.sRGBEncoding;
        }
        imagePlaneTextureCache.set(dataUrl, texture);
        return texture;
      }

      function createImagePlane(decoration) {
        const imagePlaneSpec = imagePlaneSpecForDecoration(decoration);
        const texture = imagePlaneTextureFor(imagePlaneSpec);
        if (!texture) {
          return null;
        }
        const width = Math.max(Number(imagePlaneSpec && imagePlaneSpec.width_pc) || 0.0, 1.0);
        const height = Math.max(Number(imagePlaneSpec && imagePlaneSpec.height_pc) || 0.0, 1.0);
        const opacity = clamp01(Number(decoration.opacity) || 0.0);
        if (!(opacity > 0.0)) {
          return null;
        }
        const center = decoration.center || { x: 0.0, y: 0.0, z: 0.0 };
        const geometry = new THREE.PlaneGeometry(width, height, 1, 1);
        const material = new THREE.MeshBasicMaterial({
          map: texture,
          transparent: true,
          opacity,
          side: THREE.DoubleSide,
          depthWrite: false,
          toneMapped: false,
        });
        const mesh = new THREE.Mesh(geometry, material);
        mesh.position.set(
          Number(center.x) || 0.0,
          Number(center.y) || 0.0,
          Number(center.z) || 0.0
        );
        mesh.rotation.z = Number(decoration.rotation_rad) || 0.0;
        const renderOrder = Number.isFinite(Number(decoration.render_order))
          ? Number(decoration.render_order)
          : Number(imagePlaneSpec && imagePlaneSpec.render_order);
        mesh.renderOrder = Number.isFinite(renderOrder) ? renderOrder : -20;
        registerCameraResponsiveImagePlane(mesh, {
          baseOpacity: opacity,
          hideBelowScaleBarPc: Number(imagePlaneSpec && imagePlaneSpec.hide_below_scale_bar_pc),
          fadeStartScaleBarPc: Number(imagePlaneSpec && imagePlaneSpec.fade_start_scale_bar_pc),
        });
        return mesh;
      }

      function makeLineObject(trace, materialBucket) {
        if (!trace.segments || !trace.segments.length) {
          return null;
        }
        const traceState = traceStyleStateForKey(trace.key);
        const lineColor = (traceState && traceState.color) || (trace.line || {}).color || "#ffffff";
        let lineOpacity = traceState ? clamp01(traceState.opacity) : (trace.opacity ?? 1.0);
        lineOpacity *= traceVisibilityOpacityMultiplier(trace);
        const focusedTraceKey = dendrogramFocusTraceKey();
        if (focusedTraceKey && String(trace.key) !== focusedTraceKey) {
          lineOpacity *= 0.16;
        }
        const positions = [];
        trace.segments.forEach((segment) => {
          positions.push(segment[0], segment[1], segment[2], segment[3], segment[4], segment[5]);
        });
        const geometry = new LineSegmentsGeometry();
        geometry.setPositions(positions);
        const material = new LineMaterial({
          color: lineColor,
          linewidth: Math.max((trace.line || {}).width ?? 1.0, 1.0),
          dashed: ((trace.line || {}).dash || "solid") !== "solid",
          dashSize: 8.0,
          gapSize: 5.0,
          transparent: lineOpacity < 1.0,
          opacity: lineOpacity,
          worldUnits: false,
        });
        material.resolution.set(root.clientWidth, root.clientHeight);
        const line = new LineSegments2(geometry, material);
        line.computeLineDistances();
        materialBucket.push(material);
        return line;
      }

      function makeSingleSegmentLine(startPoint, endPoint, options = {}, materialBucket = frameLineMaterials) {
        if (!(startPoint instanceof THREE.Vector3) || !(endPoint instanceof THREE.Vector3)) {
          return null;
        }
        const positions = [
          startPoint.x, startPoint.y, startPoint.z,
          endPoint.x, endPoint.y, endPoint.z,
        ];
        const geometry = new LineSegmentsGeometry();
        geometry.setPositions(positions);
        const material = new LineMaterial({
          color: String(options.color || "#ffffff"),
          linewidth: Math.max(Number(options.widthPx) || 1.0, 0.25),
          dashed: false,
          transparent: clamp01(Number(options.opacity) || 1.0) < 1.0,
          opacity: clamp01(Number(options.opacity) || 1.0),
          worldUnits: false,
        });
        material.resolution.set(root.clientWidth, root.clientHeight);
        const line = new LineSegments2(geometry, material);
        line.computeLineDistances();
        materialBucket.push(material);
        if (Number.isFinite(Number(options.renderOrder))) {
          line.renderOrder = Number(options.renderOrder);
        }
        return line;
      }

      function createTaperedPolyline(decoration) {
        const rawPoints = Array.isArray(decoration && decoration.points) ? decoration.points : [];
        if (rawPoints.length < 2) {
          return null;
        }
        const points = rawPoints
          .map((point) => new THREE.Vector3(
            Number(point && point.x) || 0.0,
            Number(point && point.y) || 0.0,
            Number(point && point.z) || 0.0
          ));
        const segmentCount = Math.max(points.length - 1, 1);
        const startWidth = Math.max(Number(decoration.start_width_px) || 2.0, 0.25);
        const endWidth = Math.max(Number(decoration.end_width_px) || 0.5, 0.1);
        const baseOpacity = clamp01(Number(decoration.opacity) || 1.0);
        const color = String(decoration.color || "#ffffff");
        const renderOrder = Number.isFinite(Number(decoration.render_order))
          ? Number(decoration.render_order)
          : -6;
        const group = new THREE.Group();
        for (let index = 0; index < segmentCount; index += 1) {
          const alpha = segmentCount <= 1 ? 0.0 : (index / (segmentCount - 1));
          const line = makeSingleSegmentLine(
            points[index],
            points[index + 1],
            {
              color,
              widthPx: startWidth + (endWidth - startWidth) * alpha,
              opacity: baseOpacity,
              renderOrder,
            },
            frameLineMaterials,
          );
          if (line) {
            group.add(line);
          }
        }
        return group;
      }

      function createSolarSystemMarker(decoration) {
        const base = decoration && decoration.base ? decoration.base : {};
        const x = Number(base.x);
        const y = Number(base.y);
        if (!Number.isFinite(x) || !Number.isFinite(y)) {
          return null;
        }
        const bottomZ = Number.isFinite(Number(decoration.bottom_z)) ? Number(decoration.bottom_z) : -400.0;
        const topZ = Number.isFinite(Number(decoration.top_z)) ? Number(decoration.top_z) : 1800.0;
        const labelZ = Number.isFinite(Number(decoration.label_z)) ? Number(decoration.label_z) : (topZ + 320.0);
        const color = String(decoration.color || "#ffe45c");
        const renderOrder = Number.isFinite(Number(decoration.render_order))
          ? Number(decoration.render_order)
          : 8;
        const hideBelowScaleBarPc = Math.max(Number(decoration.hide_below_scale_bar_pc) || 0.0, 0.0);
        const fadeStartScaleBarPc = Math.max(Number(decoration.fade_start_scale_bar_pc) || 0.0, 0.0);
        const startPoint = new THREE.Vector3(x, y, bottomZ);
        const endPoint = new THREE.Vector3(x, y, topZ);
        const group = new THREE.Group();

        const core = makeSingleSegmentLine(startPoint, endPoint, {
          color,
          widthPx: 1.55,
          opacity: 0.92,
          renderOrder: renderOrder + 0.1,
        }, frameLineMaterials);
        if (core) {
          registerCameraResponsiveImagePlane(core, {
            baseOpacity: 0.92,
            hideBelowScaleBarPc,
            fadeStartScaleBarPc,
          });
          group.add(core);
        }

        const labelText = String(decoration.label || "Solar System");
        const label = makeTextSprite(labelText, {
          color,
          size: 26,
          family: "Helvetica",
          screenStable: true,
          screenPixels: 28,
          screenScaleMultiplier: 1.0,
        });
        label.material.opacity = 0.98;
        label.position.set(x, y, labelZ);
        label.renderOrder = renderOrder + 0.2;
        registerCameraResponsiveImagePlane(label, {
          baseOpacity: 0.98,
          hideBelowScaleBarPc,
          fadeStartScaleBarPc,
        });
        screenStableTextSprites.push(label);
        group.add(label);

        return group;
      }

      function createGalacticCenterAxes(decoration) {
        const center = decoration && decoration.center ? decoration.center : {};
        const x = Number(center.x);
        const y = Number(center.y);
        const z = Number(center.z);
        if (!Number.isFinite(x) || !Number.isFinite(y) || !Number.isFinite(z)) {
          return null;
        }
        const halfLengthXY = Math.max(Number(decoration.half_length_xy) || 0.0, 1.0);
        const halfLengthZ = Math.max(Number(decoration.half_length_z) || 0.0, 1.0);
        const color = String(decoration.color || "#7a8089");
        const opacity = clamp01(Number(decoration.opacity) || 0.62);
        const widthPx = Math.max(Number(decoration.width_px) || 1.0, 0.25);
        const renderOrder = Number.isFinite(Number(decoration.render_order))
          ? Number(decoration.render_order)
          : -8;
        const group = new THREE.Group();
        const segments = [
          [new THREE.Vector3(x - halfLengthXY, y, z), new THREE.Vector3(x + halfLengthXY, y, z)],
          [new THREE.Vector3(x, y - halfLengthXY, z), new THREE.Vector3(x, y + halfLengthXY, z)],
          [new THREE.Vector3(x, y, z - halfLengthZ), new THREE.Vector3(x, y, z + halfLengthZ)],
        ];
        segments.forEach(([startPoint, endPoint]) => {
          const line = makeSingleSegmentLine(
            startPoint,
            endPoint,
            {
              color,
              widthPx,
              opacity,
              renderOrder,
            },
            frameLineMaterials,
          );
          if (line) {
            group.add(line);
          }
        });
        return group;
      }

      function buildSelectionBoxGroupForFrame(timeMyr) {
        selectionBoxHitObjects = [];
        if (!selectionBoxShouldRender()) {
          return null;
        }
        const displayState = selectionBoxDisplayStateAtTime(timeMyr);
        if (!displayState) {
          return null;
        }
        const group = new THREE.Group();
        const line = makeLineObject({
          key: "selection-box",
          segments: displayState.segments,
          line: {
            color: String(selectionBoxSpec.line_color || "#ffffff"),
            width: Math.max(Number(selectionBoxSpec.line_width) || 1.0, 1.0),
            dash: "solid",
          },
          opacity: clamp01(selectionBoxSpec.line_opacity ?? 1.0),
        }, frameLineMaterials);
        if (line) {
          group.add(line);
        }

        const cubeSidePc = 2.0 * displayState.halfWidthPc;
        const dragMesh = new THREE.Mesh(
          new THREE.BoxGeometry(cubeSidePc, cubeSidePc, cubeSidePc),
          new THREE.MeshBasicMaterial({
            color: 0xffffff,
            transparent: true,
            opacity: 0.0,
            depthWrite: false,
            side: THREE.DoubleSide,
          })
        );
        dragMesh.position.copy(displayState.centerDisplay);
        dragMesh.rotation.z = displayState.angle;
        dragMesh.userData.selectionBoxHandle = { kind: "drag" };
        group.add(dragMesh);
        selectionBoxHitObjects.push(dragMesh);

        const handleRadiusPc = clampRange(displayState.halfWidthPc * 0.075, 14.0, 40.0);
        displayState.cornersDisplay.forEach((cornerDisplay) => {
          const handleMesh = new THREE.Mesh(
            new THREE.SphereGeometry(handleRadiusPc, 12, 10),
            new THREE.MeshBasicMaterial({
              color: 0xffffff,
              transparent: true,
              opacity: 0.0,
              depthWrite: false,
            })
          );
          handleMesh.position.copy(cornerDisplay);
          handleMesh.userData.selectionBoxHandle = { kind: "resize" };
          group.add(handleMesh);
          selectionBoxHitObjects.push(handleMesh);
        });
        return group;
      }

      function smoothstep01(value) {
        const t = clampRange(Number(value) || 0.0, 0.0, 1.0);
        return t * t * (3.0 - 2.0 * t);
      }

      function fadeVisibilityFactor(timeMyr, ageNowMyr, fadeTimeMyrValue, fadeInOut) {
        const t = Number(timeMyr);
        const ageNow = Number(ageNowMyr);
        const fade = Math.max(Number(fadeTimeMyrValue) || 0.0, 0.0);
        if (!Number.isFinite(t) || !Number.isFinite(ageNow)) {
          return 1.0;
        }
        const birthTime = -ageNow;
        if (fade <= 1e-9) {
          if (t < birthTime) {
            return 0.0;
          }
          if (!fadeInOut) {
            return 1.0;
          }
          return Math.abs(t - birthTime) <= 1e-9 ? 1.0 : 0.0;
        }
        if (t < birthTime - fade) {
          return 0.0;
        }
        if (t <= birthTime) {
          return smoothstep01((t - (birthTime - fade)) / fade);
        }
        if (!fadeInOut) {
          return 1.0;
        }
        if (t <= birthTime + fade) {
          return 1.0 - smoothstep01((t - birthTime) / fade);
        }
        return 0.0;
      }

      function animatedPointState(point) {
        const motion = point && typeof point.motion === "object" ? point.motion : null;
        if (!motion) {
          return {
            size: Number(point.size ?? 0.0),
            opacity: Number(point.opacity ?? 1.0),
          };
        }
        const opacityFactor = fadeVisibilityFactor(
          motion.time_myr,
          motion.age_now_myr,
          fadeInTimeMyr,
          fadeInAndOutEnabled
        );
        return {
          size: Number(point.size ?? 0.0),
          opacity: Number(point.opacity ?? 1.0) * opacityFactor,
        };
      }

      function pointHoverText(point, trace) {
        if (point && point.hovertext) {
          return String(point.hovertext);
        }
        const motion = point && typeof point.motion === "object" ? point.motion : null;
        const selection = point && typeof point.selection === "object" ? point.selection : null;
        const rawName = String(
          (selection && selection.cluster_name)
          || (motion && motion.key)
          || (trace && trace.name)
          || "Point"
        ).trim();
        const label = rawName.replace(/_/g, " ") || "Point";
        const parts = [`<b style="font-size:16px;">${escapeHtml(label)}</b>`];
        const traceName = String((trace && trace.name) || "").trim();
        if (traceName && traceName !== rawName) {
          parts.push(escapeHtml(traceName));
        }
        const frame = typeof currentFrame === "function" ? currentFrame() : null;
        const timeMyr = Number((motion && motion.time_myr) ?? (frame && frame.time));
        const ageNowMyr = Number(motion && motion.age_now_myr);
        if (Number.isFinite(ageNowMyr)) {
          parts.push(`Age (now) = ${ageNowMyr.toFixed(1)} Myr`);
          if (Number.isFinite(timeMyr)) {
            parts.push(`Age (t) = ${(ageNowMyr + timeMyr).toFixed(1)} Myr`);
          }
        }
        const xValue = Number(point && point.x);
        const yValue = Number(point && point.y);
        const zValue = Number(point && point.z);
        if (Number.isFinite(xValue) && Number.isFinite(yValue) && Number.isFinite(zValue)) {
          parts.push(`(x,y,z) = (${xValue.toFixed(1)}, ${yValue.toFixed(1)}, ${zValue.toFixed(1)})`);
        }
        return parts.join("<br>");
      }

      function starSizeStatsForTrace(trace) {
        if (!trace || !Array.isArray(trace.points)) {
          return { hasValues: false, median: 1.0 };
        }
        if (trace._starSizeStats) {
          return trace._starSizeStats;
        }
        const values = trace.points
          .map((item) => Number(item && item.n_stars))
          .filter((value) => Number.isFinite(value) && value > 0.0)
          .sort((a, b) => a - b);
        if (!values.length) {
          trace._starSizeStats = { hasValues: false, median: 1.0 };
          return trace._starSizeStats;
        }
        const mid = Math.floor(values.length / 2);
        const median = values.length % 2
          ? values[mid]
          : 0.5 * (values[mid - 1] + values[mid]);
        trace._starSizeStats = {
          hasValues: true,
          median: Math.max(Number(median) || 1.0, 1e-3),
        };
        return trace._starSizeStats;
      }

      function sizeByStarsFactorForPoint(point, trace, traceState) {
        if (!traceState || !traceState.hasNStars) {
          return 1.0;
        }
        const nStars = Number(point && point.n_stars);
        if (!Number.isFinite(nStars) || nStars <= 0.0) {
          return 1.0;
        }
        const stats = starSizeStatsForTrace(trace);
        if (!stats.hasValues) {
          return 1.0;
        }
        const factor = clampRange(Math.sqrt(nStars / stats.median), 0.35, 2.75);
        if (sizePointsByStarsEnabled) {
          return traceState.sizeByNStarsDefault ? 1.0 : factor;
        }
        return traceState.sizeByNStarsDefault ? (1.0 / factor) : 1.0;
      }

      function addMarkerTrace(parent, trace) {
        if (!trace.points || !trace.points.length) {
          return;
        }
        const group = new THREE.Group();
        const traceState = traceStyleStateForKey(trace.key);
        const traceColor = traceState ? traceState.color : null;
        const traceOpacityMultiplier = traceState
          ? clamp01(traceState.opacity) / Math.max(clamp01(Number(trace.default_opacity ?? 1.0)), 1e-6)
          : 1.0;
        const traceVisibilityMultiplier = traceVisibilityOpacityMultiplier(trace);
        const sizeScaleFactor = traceState ? Math.max(Number(traceState.sizeScale), 0.05) : 1.0;
        trace.points.forEach((point) => {
          if (!clusterFilterPassesPoint(point)) {
            return;
          }
          const pointState = animatedPointState(point);
          if (!Number.isFinite(pointState.size) || pointState.size <= 0) {
            return;
          }
          let opacityMultiplier = 1.0;
          const pointKey = clusterFilterSelectionKeyForPoint(point) || normalizedSelectionKeyFor(point.selection);
          const focusedTraceKey = dendrogramFocusTraceKey();
          const dendrogramActiveKeys = activeDendrogramSelectionKeys();
          if (focusedTraceKey) {
            if (String(trace.key) !== focusedTraceKey) {
              opacityMultiplier *= 0.14;
            } else if (dendrogramActiveKeys.size && pointKey) {
              opacityMultiplier *= dendrogramActiveKeys.has(pointKey) ? 1.0 : 0.24;
            }
          }
          if (selectedClusterKeys.size && pointKey) {
            opacityMultiplier *= selectedClusterKeys.has(pointKey) ? 1.0 : 0.16;
          }
          const baseOpacity = Number(pointState.opacity ?? point.opacity ?? 1.0);
          const effectiveOpacity = Math.min(1.0, Math.max(0.0, baseOpacity * traceOpacityMultiplier * traceVisibilityMultiplier * opacityMultiplier * globalPointOpacityScale));
          if (effectiveOpacity <= 0.001) {
            return;
          }
          const scaleFloor = pointScale * 0.5 * Math.max(globalPointSizeScale, 0.05);
          const starsFactor = sizeByStarsFactorForPoint(point, trace, traceState);
          const scale = Math.max(
            pointState.size * sizeScaleFactor * starsFactor * globalPointSizeScale * pointScale,
            scaleFloor
          );
          const selectionKey = clusterFilterSelectionKeyForPoint(point) || normalizedSelectionKeyFor(point.selection);
          const pointPosition = new THREE.Vector3(point.x, point.y, point.z);
          const glowStrength = Math.max(globalPointGlowStrength, 0.0);
          if (glowStrength > 0.02) {
            const glowOpacity = clampRange(effectiveOpacity * (0.22 + 0.18 * glowStrength), 0.0, 0.70);
            const glowSprite = new THREE.Sprite(starGlowMaterialFor(traceColor || point.color, glowOpacity));
            const glowScale = glowScaleForPoint(
              scale,
              pointPosition,
              glowStrength
            );
            glowSprite.position.copy(pointPosition);
            glowSprite.scale.set(glowScale, glowScale, 1.0);
            glowSprite.renderOrder = -2;
            glowSprite.userData = {
              selection: point.selection || null,
              selectionKey,
              baseScale: glowScale,
              isGlow: true,
            };
            group.add(glowSprite);
            registerCameraResponsivePointSprite(glowSprite, "glow", pointPosition, scale, selectionKey);

            const coreOpacity = clampRange(effectiveOpacity * (1.00 + 0.24 * glowStrength), 0.0, 1.0);
            const coreSprite = new THREE.Sprite(starCoreMaterialFor("#ffffff", coreOpacity));
            const coreScale = starCoreScaleForPoint(
              scale,
              pointPosition,
              glowStrength
            );
            coreSprite.position.copy(pointPosition);
            coreSprite.scale.set(coreScale, coreScale, 1.0);
            coreSprite.userData = {
              hovertext: pointHoverText(point, trace),
              selection: point.selection || null,
              selectionKey,
              baseScale: coreScale,
              isGlowCore: true,
            };
            group.add(coreSprite);
            hoverTargets.push(coreSprite);
            registerCameraResponsivePointSprite(coreSprite, "core", pointPosition, scale, selectionKey);
          } else {
            const sprite = new THREE.Sprite(markerMaterialFor(point.symbol, traceColor || point.color, effectiveOpacity));
            sprite.position.copy(pointPosition);
            sprite.scale.set(scale, scale, scale);
            sprite.userData = {
              hovertext: pointHoverText(point, trace),
              selection: point.selection || null,
              selectionKey,
              baseScale: scale,
            };
            group.add(sprite);
            hoverTargets.push(sprite);
            registerCameraResponsivePointSprite(sprite, "marker", pointPosition, scale, selectionKey);
          }
        });
        parent.add(group);
      }

      function addTextTrace(parent, trace) {
        if (!trace.labels || !trace.labels.length) {
          return;
        }
        const group = new THREE.Group();
        const traceState = traceStyleStateForKey(trace.key);
        const sizeScaleFactor = traceState ? Math.max(Number(traceState.sizeScale), 0.05) : 1.0;
        const opacityValue = (traceState ? clamp01(traceState.opacity) : 1.0) * traceVisibilityOpacityMultiplier(trace);
        trace.labels.forEach((label) => {
          if (!label.text) {
            return;
          }
          const sprite = makeTextSprite(label.text, {
            color: (traceState && traceState.color) || label.color || theme.axis_color,
            size: (label.size ?? 12) * sizeScaleFactor,
            family: label.family ?? "Helvetica",
            screenStable: label.screen_stable === true,
            screenPixels: label.screen_px,
            screenScaleMultiplier: sizeScaleFactor,
          });
          sprite.material.opacity = opacityValue;
          sprite.position.set(label.x, label.y, label.z);
          if (label.screen_stable === true) {
            screenStableTextSprites.push(sprite);
          }
          group.add(sprite);
        });
        parent.add(group);
      }

      function addManualLabels(parent) {
        if (!manualLabels.length) {
          return;
        }
        const group = new THREE.Group();
        manualLabels.forEach((label, index) => {
          if (!label || !label.text) {
            return;
          }
          const sprite = makeTextSprite(label.text, {
            color: theme.axis_color || "#808080",
            size: label.size ?? DEFAULT_MANUAL_LABEL_SIZE,
            family: "Helvetica",
            screenStable: true,
            screenPixels: clampManualLabelSize(label.size ?? DEFAULT_MANUAL_LABEL_SIZE),
          });
          sprite.position.set(label.x, label.y, label.z);
          sprite.userData = {
            hovertext: `<strong>${escapeHtml(label.text)}</strong><br/>Drag to reposition`,
            manualLabelId: String(label.id || `manual-label-${index + 1}`),
          };
          screenStableTextSprites.push(sprite);
          hoverTargets.push(sprite);
          group.add(sprite);
        });
        parent.add(group);
      }

      function addDecoration(parent, decoration) {
        if (!decoration || !decoration.kind) {
          return;
        }
        if (decoration.kind === "volume_layer") {
          const layer = volumeLayersByKey.get(String(decoration.key));
          if (!layer) {
            return;
          }
          const stateKey = volumeStateKeyForLayer(layer);
          if (legendState[stateKey] === false) {
            return;
          }
          const runtime = createVolumeRuntime(layer);
          if (!runtime) {
            return;
          }
          volumeRuntimeByKey.set(String(layer.key), runtime);
          parent.add(runtime.mesh);
          return;
        }
        if (decoration.kind === "milky_way_model") {
          parent.add(createMilkyWayModel(decoration));
          return;
        }
        if (decoration.kind === "tapered_polyline") {
          const polyline = createTaperedPolyline(decoration);
          if (polyline) {
            parent.add(polyline);
          }
          return;
        }
        if (decoration.kind === "solar_system_marker") {
          const marker = createSolarSystemMarker(decoration);
          if (marker) {
            parent.add(marker);
          }
          return;
        }
        if (decoration.kind === "galactic_center_axes") {
          const axes = createGalacticCenterAxes(decoration);
          if (axes) {
            parent.add(axes);
          }
          return;
        }
        if (decoration.kind === "image_plane") {
          const imagePlane = createImagePlane(decoration);
          if (imagePlane) {
            parent.add(imagePlane);
          }
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
        if (!axesVisible) {
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
        if (isGalacticReferenceTrace(trace) && (!galacticReferenceVisible || cameraViewMode === "earth")) {
          return false;
        }
        if (isNearbyRegionLabelTrace(trace) && !nearbyRegionLabelsVisible) {
          return false;
        }
        const actionState = actionTraceVisibilityState(trace);
        if (actionState) {
          return Boolean(actionState.visible);
        }
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

      function isGalacticReferenceTrace(trace) {
        if (!trace || typeof trace !== "object") {
          return false;
        }
        const traceName = String(trace.name || "");
        if (galacticReferenceTraceNames.has(traceName)) {
          return true;
        }
        return traceName.startsWith("R = ") && traceName.endsWith("kpc");
      }

      function isNearbyRegionLabelTrace(trace) {
        if (!trace || typeof trace !== "object") {
          return false;
        }
        const traceName = String(trace.name || "");
        return nearbyRegionLabelTraceNames.has(traceName);
      }

      function medianCoordinate(values) {
        const finiteValues = values.filter((value) => Number.isFinite(value)).sort((a, b) => a - b);
        if (!finiteValues.length) {
          return NaN;
        }
        const middle = Math.floor(finiteValues.length / 2);
        if (finiteValues.length % 2 === 1) {
          return finiteValues[middle];
        }
        return 0.5 * (finiteValues[middle - 1] + finiteValues[middle]);
      }

      function focusTrackingOffsetForFrame(frame) {
        const requestedSelectionKey = normalizeMemberKey(focusSelectionKey);
        if (requestedSelectionKey && frame && Array.isArray(frame.traces)) {
          for (const trace of frame.traces) {
            if (!trace || !Array.isArray(trace.points)) {
              continue;
            }
            for (const point of trace.points) {
              const pointKey = clusterFilterSelectionKeyForPoint(point);
              if (!pointKey || pointKey !== requestedSelectionKey) {
                continue;
              }
              const xValue = Number(point.x);
              const yValue = Number(point.y);
              const zValue = Number(point.z);
              if (Number.isFinite(xValue) && Number.isFinite(yValue) && Number.isFinite(zValue)) {
                return new THREE.Vector3(xValue, yValue, zValue);
              }
            }
          }
        }

        const requestedKey = String(focusTraceKey || "");
        if (!requestedKey || !frame || !Array.isArray(frame.traces)) {
          return null;
        }
        const targetTrace = frame.traces.find((trace) => String(trace.key) === requestedKey && Array.isArray(trace.points) && trace.points.length);
        if (!targetTrace) {
          return null;
        }
        const xValues = [];
        const yValues = [];
        const zValues = [];
        targetTrace.points.forEach((point) => {
          const xValue = Number(point.x);
          const yValue = Number(point.y);
          const zValue = Number(point.z);
          if (Number.isFinite(xValue) && Number.isFinite(yValue) && Number.isFinite(zValue)) {
            xValues.push(xValue);
            yValues.push(yValue);
            zValues.push(zValue);
          }
        });
        if (!xValues.length) {
          return null;
        }
        const centerX = medianCoordinate(xValues);
        const centerY = medianCoordinate(yValues);
        const centerZ = medianCoordinate(zValues);
        if (!Number.isFinite(centerX) || !Number.isFinite(centerY) || !Number.isFinite(centerZ)) {
          return null;
        }
        return new THREE.Vector3(centerX, centerY, centerZ);
      }

      function buildTraceLegendControls(item, toggleButton) {
        const itemKey = String(item.key);
        const state = traceStyleStateForKey(itemKey);
        if (!state) {
          return null;
        }

        const controls = document.createElement("div");
        controls.className = "oviz-three-legend-controls";

        const colorInput = document.createElement("input");
        colorInput.type = "color";
        colorInput.value = cssColorToHex(state.color, "#ffffff");
        const colorField = createLegendField("Color", colorInput);

        const opacityInput = document.createElement("input");
        opacityInput.type = "range";
        opacityInput.min = "0";
        opacityInput.max = "1";
        opacityInput.step = "0.01";
        opacityInput.value = String(clamp01(state.opacity));
        const opacityField = createLegendField(`Opacity (${clamp01(state.opacity).toFixed(2)})`, opacityInput);

        const firstRow = createLegendControlRow();
        firstRow.appendChild(colorField.field);
        firstRow.appendChild(opacityField.field);
        controls.appendChild(firstRow);

        colorInput.addEventListener("input", () => {
          state.color = String(colorInput.value);
          toggleButton.style.color = state.color;
          renderFrame(currentFrameIndex);
        });
        opacityInput.addEventListener("input", () => {
          state.opacity = clamp01(opacityInput.value);
          opacityField.label.textContent = `Opacity (${state.opacity.toFixed(2)})`;
          renderFrame(currentFrameIndex);
        });

        if (state.hasPoints || state.hasLabels) {
          const sizeInput = document.createElement("input");
          sizeInput.type = "range";
          sizeInput.min = "0.25";
          sizeInput.max = "4";
          sizeInput.step = "0.05";
          sizeInput.value = String(Math.max(Number(state.sizeScale), 0.25));
          const sizeLabel = state.hasPoints && state.hasLabels
            ? "Size"
            : (state.hasLabels ? "Text size" : "Point size");
          const sizeField = createLegendField(`${sizeLabel} (${Number(state.sizeScale).toFixed(2)}x)`, sizeInput);
          controls.appendChild(sizeField.field);
          sizeInput.addEventListener("input", () => {
            state.sizeScale = Math.max(Number(sizeInput.value), 0.05);
            sizeField.label.textContent = `${sizeLabel} (${state.sizeScale.toFixed(2)}x)`;
            renderFrame(currentFrameIndex);
          });

          const summary = document.createElement("div");
          summary.className = "oviz-three-legend-summary";
          summary.textContent = state.hasPoints && state.hasLabels
            ? "Size is applied as a multiplier on the original marker and label sizes."
            : (
              state.hasLabels
                ? "Text size is applied as a multiplier on the original label sizes."
                : "Point size is applied as a multiplier on the original marker sizes."
            );
          controls.appendChild(summary);
        }

        return controls;
      }

      function buildVolumeLegendControls(item, toggleButton) {
        const layer = volumeLayerForKey(item.key);
        const stateKey = layer ? volumeStateKeyForLayer(layer) : String(item.key || "");
        const state = stateKey ? volumeStateByKey[stateKey] : null;
        if (!layer || !state) {
          return null;
        }

        const controls = document.createElement("div");
        controls.className = "oviz-three-legend-controls";

        const visibleToggle = document.createElement("label");
        visibleToggle.className = "oviz-three-legend-toggle";
        const visibleInput = document.createElement("input");
        visibleInput.type = "checkbox";
        visibleToggle.appendChild(visibleInput);
        const visibleText = document.createElement("span");
        visibleText.textContent = volumeLayerTimeMyr(layer) !== null ? "Show volume series" : "Show volume";
        visibleToggle.appendChild(visibleText);
        controls.appendChild(visibleToggle);

        let showAllTimesInput = null;
        let showAllTimesToggle = null;
        let showAllTimesText = null;
        if (volumeSupportsShowAllTimes(layer)) {
          showAllTimesToggle = document.createElement("label");
          showAllTimesToggle.className = "oviz-three-legend-toggle";
          showAllTimesInput = document.createElement("input");
          showAllTimesInput.type = "checkbox";
          showAllTimesToggle.appendChild(showAllTimesInput);
          showAllTimesText = document.createElement("span");
          showAllTimesText.textContent = layer.co_rotate_with_frame
            ? "Show at all times (co-rotating)"
            : "Show at all times";
          showAllTimesToggle.appendChild(showAllTimesText);
          controls.appendChild(showAllTimesToggle);
        }

        const colormapSelect = document.createElement("select");
        ((layer.colormap_options || [])).forEach((option) => {
          const optionEl = document.createElement("option");
          optionEl.value = String(option.name);
          optionEl.textContent = String(option.label || option.name);
          colormapSelect.appendChild(optionEl);
        });
        const colormapField = createLegendField("Colormap", colormapSelect);
        controls.appendChild(colormapField.field);

        const stretchSelect = document.createElement("select");
        volumeStretchOptions().forEach((option) => {
          const optionEl = document.createElement("option");
          optionEl.value = String(option.value);
          optionEl.textContent = String(option.label);
          stretchSelect.appendChild(optionEl);
        });
        const stretchField = createLegendField("Stretch", stretchSelect);
        controls.appendChild(stretchField.field);

        const vminInput = document.createElement("input");
        vminInput.type = "number";
        const vminField = createLegendField("vmin", vminInput);

        const vmaxInput = document.createElement("input");
        vmaxInput.type = "number";
        const vmaxField = createLegendField("vmax", vmaxInput);

        const valueRow = createLegendControlRow();
        valueRow.appendChild(vminField.field);
        valueRow.appendChild(vmaxField.field);
        controls.appendChild(valueRow);

        const opacityInput = document.createElement("input");
        opacityInput.type = "range";
        opacityInput.min = "0";
        opacityInput.max = "1";
        opacityInput.step = "0.01";
        const opacityField = createLegendField("Opacity", opacityInput);
        controls.appendChild(opacityField.field);

        const alphaInput = document.createElement("input");
        alphaInput.type = "range";
        alphaInput.min = "1";
        alphaInput.max = "200";
        alphaInput.step = "1";
        const alphaField = createLegendField("Alpha coef", alphaInput);
        controls.appendChild(alphaField.field);

        const samplesInput = document.createElement("input");
        samplesInput.type = "range";
        samplesInput.min = "24";
        samplesInput.max = "768";
        samplesInput.step = "1";
        const samplesField = createLegendField("Samples", samplesInput);
        controls.appendChild(samplesField.field);

        const summary = document.createElement("div");
        summary.className = "oviz-three-legend-summary";
        controls.appendChild(summary);

        function refreshVolumeControls(syncInputs = true) {
          const summaryLayer = frameVolumeLayerForStateKey(stateKey) || layer;
          clampVolumeStateForLayer(summaryLayer, state);
          if (syncInputs) {
            visibleInput.checked = state.visible !== false;
            if (showAllTimesInput) {
              showAllTimesInput.checked = Boolean(state.showAllTimes);
            }
            colormapSelect.value = String(state.colormap);
            stretchSelect.value = String(state.stretch);
            syncVolumeWindowInput(vminInput, state.vmin, summaryLayer);
            syncVolumeWindowInput(vmaxInput, state.vmax, summaryLayer);
            opacityInput.value = String(state.opacity);
            alphaInput.value = String(state.alphaCoef);
            samplesInput.value = String(state.steps);
          }
          opacityField.label.textContent = `Opacity (${Number(state.opacity).toFixed(2)})`;
          alphaField.label.textContent = `Alpha coef (${Math.round(Number(state.alphaCoef))})`;
          samplesField.label.textContent = `Samples (${Math.round(Number(state.steps))})`;
          summary.textContent = volumeSummaryTextFor(summaryLayer, state);
          toggleButton.style.color = volumeLegendColorForLayer(summaryLayer);
        }

        function updateVolumeFromLegend(shouldRerenderLegend = false) {
          activeVolumeKey = String(stateKey);
          if (shouldRerenderLegend) {
            renderLegend();
          }
          updateActiveVolumeRuntime();
        }

        refreshVolumeControls(true);

        visibleInput.addEventListener("change", () => {
          state.visible = Boolean(visibleInput.checked);
          legendState[String(item.key || stateKey)] = state.visible;
          legendState[stateKey] = state.visible;
          const entryEl = toggleButton.closest ? toggleButton.closest(".oviz-three-legend-entry") : null;
          if (entryEl) {
            entryEl.dataset.active = state.visible ? "true" : "false";
          }
          toggleButton.dataset.active = state.visible ? "true" : "false";
          refreshVolumeControls(false);
          updateVolumeFromLegend(false);
        });
        if (showAllTimesInput) {
          showAllTimesInput.addEventListener("change", () => {
            state.showAllTimes = Boolean(showAllTimesInput.checked);
            refreshVolumeControls(false);
            updateVolumeFromLegend(false);
          });
        }
        colormapSelect.addEventListener("change", () => {
          state.colormap = String(colormapSelect.value);
          refreshVolumeControls(false);
          updateVolumeFromLegend(false);
        });
        stretchSelect.addEventListener("change", () => {
          state.stretch = normalizeVolumeStretch(stretchSelect.value);
          refreshVolumeControls(false);
          updateVolumeFromLegend(false);
        });
        function updateVolumeWindowFromLegendInput(inputEl, key, options = {}) {
          const value = finiteNumberInputValue(inputEl);
          if (value === null) {
            return;
          }
          state[key] = value;
          if (options.syncInputs) {
            refreshVolumeControls(true);
          } else {
            const summaryLayer = frameVolumeLayerForStateKey(stateKey) || layer;
            clampVolumeStateForLayer(summaryLayer, state);
            summary.textContent = volumeSummaryTextFor(summaryLayer, state);
            toggleButton.style.color = volumeLegendColorForLayer(summaryLayer);
          }
          updateVolumeFromLegend(false);
        }
        vminInput.addEventListener("input", () => {
          updateVolumeWindowFromLegendInput(vminInput, "vmin", { syncInputs: false });
        });
        vminInput.addEventListener("change", () => {
          updateVolumeWindowFromLegendInput(vminInput, "vmin", { syncInputs: true });
        });
        vmaxInput.addEventListener("input", () => {
          updateVolumeWindowFromLegendInput(vmaxInput, "vmax", { syncInputs: false });
        });
        vmaxInput.addEventListener("change", () => {
          updateVolumeWindowFromLegendInput(vmaxInput, "vmax", { syncInputs: true });
        });
        opacityInput.addEventListener("input", () => {
          state.opacity = Number(opacityInput.value);
          refreshVolumeControls(false);
          updateVolumeFromLegend(false);
        });
        alphaInput.addEventListener("input", () => {
          state.alphaCoef = Number(alphaInput.value);
          refreshVolumeControls(false);
          updateVolumeFromLegend(false);
        });
        samplesInput.addEventListener("input", () => {
          state.steps = Number(samplesInput.value);
          refreshVolumeControls(false);
          updateVolumeFromLegend(false);
        });

        return controls;
      }

__LEGEND_RUNTIME_JS__

__WIDGET_RUNTIME_JS__

__WIDGET_CONTENT_RUNTIME_JS__

__INTERACTION_RUNTIME_JS__

__SCENE_RUNTIME_JS__

__VIEWER_RUNTIME_JS__

__ACTION_RUNTIME_JS__

      function manualLabelHitFromEvent(event) {
        const hitObject = pickSprite(event);
        if (!hitObject || !hitObject.userData || !hitObject.userData.manualLabelId) {
          return null;
        }
        const label = manualLabelById(hitObject.userData.manualLabelId);
        if (!label) {
          return null;
        }
        return { hitObject, label };
      }

      function manualLabelPlanePointFromEvent(event, plane) {
        pointerRayFromEvent(event);
        const point = new THREE.Vector3();
        return raycaster.ray.intersectPlane(plane, point) ? point : null;
      }

      function startManualLabelInteraction(event) {
        if (minimalModeEnabled || widgetPointerState || lassoState || event.button !== 0) {
          return false;
        }
        const hitInfo = manualLabelHitFromEvent(event);
        if (!hitInfo) {
          return false;
        }
        activeManualLabelId = String(hitInfo.label.id || "");
        syncManualLabelDraftFromSelection();
        renderManualLabelControls();
        const worldPosition = new THREE.Vector3();
        hitInfo.hitObject.getWorldPosition(worldPosition);
        const cameraDirection = new THREE.Vector3();
        camera.getWorldDirection(cameraDirection);
        if (cameraDirection.lengthSq() <= 1e-12) {
          cameraDirection.set(0.0, 0.0, -1.0);
        } else {
          cameraDirection.normalize();
        }
        const dragPlane = new THREE.Plane().setFromNormalAndCoplanarPoint(cameraDirection, worldPosition);
        const planePoint = manualLabelPlanePointFromEvent(event, dragPlane) || worldPosition.clone();
        manualLabelPointerState = {
          pointerId: event.pointerId,
          labelId: String(hitInfo.label.id || ""),
          plane: dragPlane.clone(),
          dragOffsetWorld: planePoint.clone().sub(worldPosition),
        };
        controls.enabled = false;
        document.body.style.userSelect = "none";
        if (typeof canvas.setPointerCapture === "function" && event.pointerId !== undefined) {
          try {
            canvas.setPointerCapture(event.pointerId);
          } catch (_err) {
          }
        }
        tooltipEl.style.display = "none";
        suppressNextCanvasClick = true;
        event.preventDefault();
        return true;
      }

      function updateManualLabelInteraction(event) {
        if (!manualLabelPointerState) {
          return false;
        }
        const selectedLabel = manualLabelById(manualLabelPointerState.labelId);
        if (!selectedLabel) {
          return finishManualLabelInteraction(event);
        }
        const planePoint = manualLabelPlanePointFromEvent(event, manualLabelPointerState.plane);
        if (!planePoint) {
          return true;
        }
        const nextWorld = planePoint.clone().sub(manualLabelPointerState.dragOffsetWorld);
        const nextLocal = nextWorld.sub(plotGroup.position);
        selectedLabel.x = Number.isFinite(nextLocal.x) ? nextLocal.x : selectedLabel.x;
        selectedLabel.y = Number.isFinite(nextLocal.y) ? nextLocal.y : selectedLabel.y;
        selectedLabel.z = Number.isFinite(nextLocal.z) ? nextLocal.z : selectedLabel.z;
        renderFrame(currentFrameIndex);
        event.preventDefault();
        return true;
      }

      function finishManualLabelInteraction(event) {
        if (!manualLabelPointerState) {
          return false;
        }
        if (typeof canvas.releasePointerCapture === "function" && manualLabelPointerState.pointerId !== undefined) {
          try {
            canvas.releasePointerCapture(manualLabelPointerState.pointerId);
          } catch (_err) {
          }
        }
        manualLabelPointerState = null;
        controls.enabled = true;
        document.body.style.userSelect = "";
        renderManualLabelControls();
        if (event) {
          event.preventDefault();
        }
        return true;
      }

      function onPointerMove(event) {
        if (skyViewDragState) {
          tooltipEl.style.display = "none";
          return;
        }
        if (selectionBoxPointerState || manualLabelPointerState) {
          tooltipEl.style.display = "none";
          return;
        }
        if (lassoState) {
          onLassoPointerMove(event);
          return;
        }
        const rect = canvas.getBoundingClientRect();
        pointer.x = ((event.clientX - rect.left) / rect.width) * 2.0 - 1.0;
        pointer.y = -((event.clientY - rect.top) / rect.height) * 2.0 + 1.0;
        raycaster.setFromCamera(pointer, camera);
        const hits = raycaster.intersectObjects(hoverTargets, false);
        const hitObject = hits.length ? hits[0].object : null;
        const hoveredSelection = hitObject && hitObject.userData ? hitObject.userData.selection : null;
        setLocalHoveredClusterKey(hoveredSelection ? normalizedSelectionKeyFor(hoveredSelection) : "");
        if (!hitObject || !hitObject.userData.hovertext) {
          tooltipEl.style.display = "none";
          return;
        }
        tooltipEl.style.display = "block";
        tooltipEl.innerHTML = hitObject.userData.hovertext;
        tooltipEl.style.left = `${event.clientX - rect.left + 14}px`;
        tooltipEl.style.top = `${event.clientY - rect.top + 14}px`;
      }

      function onPointerLeave() {
        tooltipEl.style.display = "none";
        setLocalHoveredClusterKey("");
      }

      function initControls() {
        renderActionBar();
        setToolsDrawerOpen(Boolean(toolsShellEl && toolsShellEl.dataset.open === "true"));
        setControlsDrawerOpen(Boolean(controlsShellEl && controlsShellEl.dataset.open === "true"));
        setSkyControlsDrawerOpen(Boolean(skyControlsShellEl && skyControlsShellEl.dataset.open === "true"));
        setLegendPanelOpen(legendPanelOpen);
        applyLegendPanelRect(legendPanelRectState || defaultLegendPanelRect());
        groupSelectEl.innerHTML = "";
        const groups = sceneSpec.group_order || [defaultGroup];
        groups.forEach((groupName) => {
          const option = document.createElement("option");
          option.value = groupName;
          option.textContent = groupName;
          groupSelectEl.appendChild(option);
        });
        renderLegendGroupDropdown(groups);
        syncLegendGroupChooser({ centerActive: true, smooth: false });
        if (focusGroupSelectEl) {
          focusGroupSelectEl.innerHTML = "";
          const focusNoneOption = document.createElement("option");
          focusNoneOption.value = "";
          focusNoneOption.textContent = "None";
          focusGroupSelectEl.appendChild(focusNoneOption);
          (animationSpec.focus_options || []).forEach((item) => {
            const option = document.createElement("option");
            option.value = String(item.key || "");
            option.textContent = String(item.name || item.key || "");
            focusGroupSelectEl.appendChild(option);
          });
        }
        renderWidgetMenu();
        groupSelectEl.addEventListener("change", () => {
          setLegendGroup(groupSelectEl.value, { centerActive: true, smooth: true });
        });
        if (legendPanelEl) {
          legendPanelEl.addEventListener("pointerdown", (event) => {
            const target = event.target;
            if (!target || actionInterruptsMuted()) {
              return;
            }
            if (
              target.closest(".oviz-three-legend-entry")
              || target.closest(".oviz-three-legend-section-toggle")
              || target.closest(".oviz-three-group-dropdown")
              || target.closest(".oviz-three-group-select")
            ) {
              interruptActionRun("legend", { disableOrbit: false });
            }
          });
        }
        renderSceneControls();
        ensureInitialSkyLayers();
        applySkyLayerState({ update: false, renderLegend: false });
        applyInitialSkyDomeSource();
        if (widgetSelectEl) {
          widgetSelectEl.addEventListener("change", () => {
            const widgetKey = String(widgetSelectEl.value || "");
            if (widgetKey) {
              setWidgetMode(widgetKey, "normal");
            }
            widgetSelectEl.value = "";
          });
        }
        if (zenModeButtonEl) {
          zenModeButtonEl.addEventListener("click", () => {
            setZenMode(!zenModeEnabled);
            focusViewer();
          });
        }
        if (resetViewButtonEl) {
          resetViewButtonEl.addEventListener("click", () => {
            resetCameraAndSelections();
            focusViewer();
          });
        }
        if (saveStateButtonEl) {
          saveStateButtonEl.addEventListener("click", async () => {
            saveStateButtonEl.disabled = true;
            const previousText = saveStateButtonEl.textContent;
            saveStateButtonEl.textContent = "Saving...";
            try {
              await saveSceneStateToHtml();
            } finally {
              saveStateButtonEl.disabled = false;
              saveStateButtonEl.textContent = previousText;
            }
          });
        }

        sliderEl.max = String(Math.max(frameSpecs.length - 1, 0));
        renderTimeSliderTicks();
        const handleLegendPanelToggle = (event) => {
          event.preventDefault();
          event.stopPropagation();
          setLegendPanelOpen(!legendPanelOpen);
          if (activeLegendEditorKey) {
            renderLegend();
          }
          focusViewer();
        };
        if (legendPanelToggleEl) {
          legendPanelToggleEl.addEventListener("pointerdown", (event) => {
            event.stopPropagation();
          });
          legendPanelToggleEl.addEventListener("click", handleLegendPanelToggle);
        }
        legendPanelEl.addEventListener("click", (event) => {
          const toggleButton = event.target && event.target.closest
            ? event.target.closest(".oviz-three-legend-panel-toggle")
            : null;
          if (!toggleButton) {
            return;
          }
          handleLegendPanelToggle(event);
        });
        if (legendTraceSectionToggleEl) {
          legendTraceSectionToggleEl.addEventListener("click", (event) => {
            event.preventDefault();
            event.stopPropagation();
            setLegendSectionOpen("traces", !(legendSectionOpenState.traces !== false));
            focusViewer();
          });
        }
        if (legendVolumeSectionToggleEl) {
          legendVolumeSectionToggleEl.addEventListener("click", (event) => {
            event.preventDefault();
            event.stopPropagation();
            setLegendSectionOpen("volumes", !(legendSectionOpenState.volumes !== false));
            focusViewer();
          });
        }
        if (legendDragHandleEl) {
          legendDragHandleEl.addEventListener("pointerdown", onLegendPointerStart);
        }
        legendResizeEls.forEach((handle) => handle.addEventListener("pointerdown", onLegendPointerStart));
        sliderEl.addEventListener("pointerdown", () => {
          handleManualCameraInteractionStart();
        });
        sliderEl.addEventListener("input", () => {
          handleManualCameraInteractionStart();
          if (!actionInterruptsMuted()) {
            interruptActionRun("time", { disableOrbit: false });
          }
          pause({ snap: false });
          scheduleSliderScrubRender(Number(sliderEl.value));
        });
        if (playBackwardButtonEl) {
          playBackwardButtonEl.addEventListener("click", () => {
            if (!actionInterruptsMuted()) {
              interruptActionRun("time", { disableOrbit: false });
            }
            play(-1);
          });
        }
        if (playForwardButtonEl) {
          playForwardButtonEl.addEventListener("click", () => {
            if (!actionInterruptsMuted()) {
              interruptActionRun("time", { disableOrbit: false });
            }
            play(1);
          });
        }
        if (toolsToggleEl && toolsShellEl) {
          toolsToggleEl.addEventListener("click", () => {
            setToolsDrawerOpen(toolsShellEl.dataset.open !== "true");
          });
        }
        if (controlsToggleEl && controlsShellEl) {
          controlsToggleEl.addEventListener("click", () => {
            setControlsDrawerOpen(controlsShellEl.dataset.open !== "true");
          });
        }
        if (skyControlsToggleEl && skyControlsShellEl) {
          skyControlsToggleEl.addEventListener("click", () => {
            setSkyControlsDrawerOpen(skyControlsShellEl.dataset.open !== "true");
          });
        }
        if (keyHelpButtonEl && keyHelpEl) {
          keyHelpButtonEl.addEventListener("click", () => {
            const nextOpen = keyHelpEl.dataset.open !== "true";
            setKeyHelpOpen(nextOpen);
            if (!nextOpen) {
              focusViewer();
            }
          });
        }
        if (keyHelpCloseEl) {
          keyHelpCloseEl.addEventListener("click", () => {
            setKeyHelpOpen(false);
            focusViewer();
          });
        }
        if (clearSelectionButtonEl) {
          clearSelectionButtonEl.addEventListener("click", clearClusterSelections);
        }
        if (clickSelectToggleEl) {
          clickSelectToggleEl.addEventListener("change", () => {
            clickSelectionEnabled = Boolean(clickSelectToggleEl.checked);
            updateSelectionUI();
          });
        }
        if (volumeLassoToggleEl) {
          volumeLassoToggleEl.addEventListener("change", () => {
            lassoVolumeSelectionEnabled = Boolean(volumeLassoToggleEl.checked);
            updateSelectionUI();
            const frame = currentFrame();
            if (frame && frameVolumeLayers(frame).length) {
              renderFrame(currentFrameIndex);
            }
            if (skySpec.enabled && !currentSelection) {
              updateSkyPanel();
            }
          });
        }
        if (themeSelectEl) {
          themeSelectEl.addEventListener("change", () => {
            applyThemePreset(themeSelectEl.value);
          });
        }
        if (skyDomeSourceSelectEl) {
          skyDomeSourceSelectEl.addEventListener("change", () => {
            setSkyDomeSourceByKey(skyDomeSourceSelectEl.value, { force: true });
            updateSkyPanel();
            updateSkyDomeFromControls({ forceTiles: true });
            focusViewer();
          });
        }
        if (skyLayerAddButtonEl) {
          skyLayerAddButtonEl.addEventListener("click", () => {
            const customSurvey = skyLayerCustomInputEl ? cleanSkyLayerSurvey(skyLayerCustomInputEl.value) : "";
            const presetSurvey = skyLayerPresetSelectEl ? cleanSkyLayerSurvey(skyLayerPresetSelectEl.value) : "";
            addOrActivateSkyLayer(customSurvey || presetSurvey || defaultSkyLayerSurvey());
            if (skyLayerCustomInputEl) {
              skyLayerCustomInputEl.value = "";
            }
            focusViewer();
          });
        }
        if (skyLayerCustomInputEl) {
          skyLayerCustomInputEl.addEventListener("focus", () => {
            handleSkyLayerSearchInput();
          });
          skyLayerCustomInputEl.addEventListener("input", () => {
            handleSkyLayerSearchInput();
          });
          skyLayerCustomInputEl.addEventListener("change", () => {
            updateSkyLayerSearchOptions(skyLayerCustomInputEl.value);
          });
          skyLayerCustomInputEl.addEventListener("keydown", (event) => {
            if (event.key !== "Enter") {
              return;
            }
            event.preventDefault();
            if (skyLayerAddButtonEl) {
              skyLayerAddButtonEl.click();
            }
          });
        }
        if (skyDomeVisibleToggleEl) {
          skyDomeVisibleToggleEl.addEventListener("change", () => {
            skyDomeSpec.enabled = Boolean(skyDomeVisibleToggleEl.checked);
            updateSkyDomeFromControls({ forceTiles: true });
            focusViewer();
          });
        }
        if (skyDomeForceVisibleToggleEl) {
          skyDomeForceVisibleToggleEl.addEventListener("change", () => {
            skyDomeForceVisible = Boolean(skyDomeForceVisibleToggleEl.checked);
            updateSkyDomeFromControls({ forceTiles: true });
            focusViewer();
          });
        }
        if (skyDomeOpacityEl) {
          skyDomeOpacityEl.addEventListener("input", () => {
            const opacity = Math.min(Math.max(Number(skyDomeOpacityEl.value), 0.0), 1.0);
            const layer = activeSkyLayer();
            if (layer) {
              layer.opacity = opacity;
              layer.visible = opacity > 0.0;
              activeSkyLayerKey = layer.key;
              applySkyLayerState({ forceTiles: false, syncControls: false });
              return;
            }
            skyDomeSpec.opacity = opacity;
            updateSkyDomeFromControls({ syncControls: false });
          });
          skyDomeOpacityEl.addEventListener("change", () => {
            syncSkyDomeControls();
          });
        }
        if (skyDomeBrightnessEl) {
          skyDomeBrightnessEl.addEventListener("input", () => {
            skyDomeSpec.hips_brightness = Math.min(Math.max(Number(skyDomeBrightnessEl.value), 0.1), 8.0);
            updateSkyDomeFromControls();
          });
        }
        if (skyDomeContrastEl) {
          skyDomeContrastEl.addEventListener("input", () => {
            skyDomeSpec.hips_contrast = Math.min(Math.max(Number(skyDomeContrastEl.value), 0.1), 4.0);
            updateSkyDomeFromControls();
          });
        }
        if (skyDomeGammaEl) {
          skyDomeGammaEl.addEventListener("input", () => {
            skyDomeSpec.hips_gamma = Math.min(Math.max(Number(skyDomeGammaEl.value), 0.2), 4.0);
            updateSkyDomeFromControls();
          });
        }
        if (scrollSpeedEl) {
          scrollSpeedEl.addEventListener("input", () => {
            globalScrollSpeed = Number(scrollSpeedEl.value);
            applyGlobalControlState();
            renderSceneControls();
          });
        }
        if (cameraFovEl) {
          cameraFovEl.addEventListener("input", () => {
            camera.fov = Number(cameraFovEl.value);
            applyGlobalControlState();
            renderSceneControls();
          });
        }
        if (globalPointSizeEl) {
          globalPointSizeEl.addEventListener("input", () => {
            globalPointSizeScale = Number(globalPointSizeEl.value);
            applyGlobalControlState();
            renderSceneControls();
            renderFrame(currentFrameIndex);
          });
        }
        if (globalPointOpacityEl) {
          globalPointOpacityEl.addEventListener("input", () => {
            globalPointOpacityScale = Number(globalPointOpacityEl.value);
            applyGlobalControlState();
            renderSceneControls();
            renderFrame(currentFrameIndex);
          });
        }
        if (globalPointGlowEl) {
          globalPointGlowEl.addEventListener("input", () => {
            globalPointGlowStrength = Number(globalPointGlowEl.value);
            applyGlobalControlState();
            renderSceneControls();
            renderFrame(currentFrameIndex);
          });
        }
        if (sizeByStarsToggleEl) {
          sizeByStarsToggleEl.addEventListener("change", () => {
            sizePointsByStarsEnabled = Boolean(sizeByStarsToggleEl.checked);
            renderSceneControls();
            renderFrame(currentFrameIndex);
          });
        }
        if (focusGroupSelectEl) {
          focusGroupSelectEl.addEventListener("change", () => {
            focusTraceKey = String(focusGroupSelectEl.value || "");
            focusSelectionKey = "";
            applyGlobalControlState();
            renderSceneControls();
            renderFrame(currentFrameIndex);
          });
        }
        if (fadeTimeEl) {
          fadeTimeEl.addEventListener("change", () => {
            fadeInTimeMyr = Number(fadeTimeEl.value);
            applyGlobalControlState();
            renderSceneControls();
            renderFrame(currentFrameIndex);
          });
        }
        if (fadeInOutToggleEl) {
          fadeInOutToggleEl.addEventListener("change", () => {
            fadeInAndOutEnabled = Boolean(fadeInOutToggleEl.checked);
            applyGlobalControlState();
            renderSceneControls();
            renderFrame(currentFrameIndex);
          });
        }
        if (axesVisibleToggleEl) {
          axesVisibleToggleEl.addEventListener("change", () => {
            axesVisible = Boolean(axesVisibleToggleEl.checked);
            renderSceneControls();
            buildAxes();
          });
        }
        if (galacticReferenceToggleEl) {
          galacticReferenceToggleEl.addEventListener("change", () => {
            galacticReferenceVisible = Boolean(galacticReferenceToggleEl.checked);
            renderSceneControls();
            renderFrame(currentFrameIndex);
          });
        }
        nearbyRegionLabelsToggleEls.forEach((toggleEl) => {
          toggleEl.addEventListener("change", () => {
            nearbyRegionLabelsVisible = Boolean(toggleEl.checked);
            renderSceneControls();
            renderFrame(currentFrameIndex);
          });
        });
        if (manualLabelSelectEl) {
          manualLabelSelectEl.addEventListener("change", () => {
            activeManualLabelId = String(manualLabelSelectEl.value || "");
            ensureActiveManualLabel(activeManualLabelId);
            syncManualLabelDraftFromSelection();
            renderManualLabelControls();
          });
        }
        if (manualLabelTextEl) {
          manualLabelTextEl.addEventListener("input", () => {
            manualLabelDraftText = String(manualLabelTextEl.value || "").slice(0, MAX_MANUAL_LABEL_TEXT_LENGTH);
          });
        }
        if (manualLabelSizeEl) {
          manualLabelSizeEl.addEventListener("input", () => {
            manualLabelDraftSize = clampManualLabelSize(manualLabelSizeEl.value);
          });
          manualLabelSizeEl.addEventListener("change", () => {
            manualLabelDraftSize = clampManualLabelSize(manualLabelSizeEl.value);
            renderManualLabelControls();
          });
        }
        if (manualLabelAddButtonEl) {
          manualLabelAddButtonEl.addEventListener("click", () => {
            addManualLabelFromDraft();
            focusViewer();
          });
        }
        if (manualLabelApplyButtonEl) {
          manualLabelApplyButtonEl.addEventListener("click", () => {
            applyManualLabelDraftToSelection();
            focusViewer();
          });
        }
        if (manualLabelDeleteButtonEl) {
          manualLabelDeleteButtonEl.addEventListener("click", () => {
            deleteActiveManualLabel();
            focusViewer();
          });
        }
        if (viewFromEarthButtonEl) {
          viewFromEarthButtonEl.addEventListener("click", () => {
            if (!actionInterruptsMuted()) {
              interruptActionRun("camera", { disableOrbit: true });
            }
            viewFromEarth();
          });
        }
        if (earthViewToggleButtonEl) {
          earthViewToggleButtonEl.addEventListener("click", () => {
            if (!actionInterruptsMuted()) {
              interruptActionRun("camera", { disableOrbit: true });
            }
            toggleEarthView();
            focusViewer();
          });
        }
        orbitCameraButtons.forEach((buttonEl) => {
          buttonEl.addEventListener("click", () => {
            if (!actionInterruptsMuted()) {
              interruptActionRun("camera", { disableOrbit: true });
            }
            setCameraAutoOrbitEnabled(!cameraAutoOrbitEnabled);
            focusViewer();
          });
        });
        if (resetCameraButtonEl) {
          resetCameraButtonEl.addEventListener("click", () => {
            if (!actionInterruptsMuted()) {
              interruptActionRun("camera", { disableOrbit: true });
            }
            resetCameraView();
            renderSceneControls();
          });
        }
        if (resetControlsButtonEl) {
          resetControlsButtonEl.addEventListener("click", () => {
            activeThemeKey = "default";
            globalScrollSpeed = 1.0;
            globalPointSizeScale = 1.0;
            globalPointOpacityScale = 1.0;
            globalPointGlowStrength = 0.60;
            sizePointsByStarsEnabled = false;
            fadeInTimeMyr = Number(animationSpec.fade_in_time_default);
            fadeInAndOutEnabled = Boolean(animationSpec.fade_in_and_out_default);
            focusTraceKey = String(animationSpec.focus_trace_key_default || "");
            axesVisible = Boolean(sceneSpec.show_axes);
            galacticReferenceVisible = true;
            nearbyRegionLabelsVisible = true;
            skyDomeSpec.enabled = skyDomeDefaultEnabled;
            skyDomeForceVisible = skyDomeDefaultForceVisible;
            skyDomeSpec.opacity = skyDomeDefaultOpacity;
            skyDomeSpec.hips_brightness = skyDomeDefaultHipsBrightness;
            skyDomeSpec.hips_contrast = skyDomeDefaultHipsContrast;
            skyDomeSpec.hips_gamma = skyDomeDefaultHipsGamma;
            resetSkyDomeSourceSelection();
            setCameraAutoOrbitEnabled(false);
            cameraViewMode = "free";
            earthViewFocusDistance = null;
            earthViewReturnCameraState = null;
            camera.fov = Number(initialCameraState.fov);
            if (typeof applyActionCameraViewOffset === "function") {
              applyActionCameraViewOffset(initialCameraState.viewOffset);
            }
            applyGlobalControlState();
            applyCameraViewMode();
            applyThemePreset(activeThemeKey, { rerender: false });
            renderSceneControls();
            updateSkyDomeFromControls({ forceTiles: true });
            buildAxes();
            renderLegend();
            updateSkyPanel();
            renderFrame(currentFrameIndex);
          });
        }
        widgetDragHandles.forEach((handle) => handle.addEventListener("pointerdown", onWidgetPointerStart));
        widgetResizeEls.forEach((handle) => handle.addEventListener("pointerdown", onWidgetPointerStart));
        setZenMode(zenModeEnabled);
      }

      function onCanvasPointerDown(event) {
        if (!actionInterruptsMuted()) {
          interruptActionRun("camera", { disableOrbit: true });
        }
        handleManualCameraInteractionStart();
        focusViewer();
        if (startSkyViewCameraDrag(event)) {
          return;
        }
        if (startSelectionBoxInteraction(event)) {
          return;
        }
        if (startManualLabelInteraction(event)) {
          return;
        }
        startLassoSelection(event);
      }

      function onWindowPointerMove(event) {
        if (updateSkyViewCameraDrag(event)) {
          return;
        }
        if (updateSelectionBoxInteraction(event)) {
          return;
        }
        if (updateManualLabelInteraction(event)) {
          return;
        }
        if (updateLegendPanelInteraction(event)) {
          return;
        }
        if (updateScaleBarInteraction(event)) {
          return;
        }
        onWidgetPointerMove(event);
      }

      function onWindowPointerDown(event) {
        const target = event.target;
        if (target) {
          if (toolsShellEl && toolsShellEl.dataset.open === "true" && !toolsShellEl.contains(target)) {
            setToolsDrawerOpen(false);
          }
          if (controlsShellEl && controlsShellEl.dataset.open === "true" && !controlsShellEl.contains(target)) {
            setControlsDrawerOpen(false);
          }
          if (skyControlsShellEl && skyControlsShellEl.dataset.open === "true" && !skyControlsShellEl.contains(target)) {
            setSkyControlsDrawerOpen(false);
          }
        }
        if (!activeLegendEditorKey || !target) {
          return;
        }
        if ((legendPopoverEl && legendPopoverEl.contains(target)) || (legendPanelEl && legendPanelEl.contains(target))) {
          return;
        }
        closeLegendPopover();
        renderLegend();
      }

      function onWindowPointerEnd(event) {
        if (finishSkyViewCameraDrag(event)) {
          return;
        }
        if (finishSelectionBoxInteraction(event)) {
          return;
        }
        if (finishManualLabelInteraction(event)) {
          return;
        }
        if (finishLegendPanelInteraction(event)) {
          return;
        }
        if (finishScaleBarInteraction(event)) {
          return;
        }
        onWidgetPointerEnd(event);
        finishLassoSelection(event);
      }

      function updateAnimatedFramePlayback(now) {
        if (timeActionTrack) {
          return;
        }
        if (playbackDirection !== 0 && frameSpecs.length > 1) {
          if (lastPlaybackAdvanceTimestamp === null) {
            lastPlaybackAdvanceTimestamp = now;
            return;
          }
          const elapsedMs = Math.max(0.0, now - lastPlaybackAdvanceTimestamp);
          if (elapsedMs < playbackIntervalMs) {
            return;
          }
          const steps = Math.max(1, Math.floor(elapsedMs / playbackIntervalMs));
          lastPlaybackAdvanceTimestamp += steps * playbackIntervalMs;
          const frameCount = frameSpecs.length;
          const nextIndex = ((currentFrameIndex + (playbackDirection * steps)) % frameCount + frameCount) % frameCount;
          renderFrame(nextIndex);
        }
      }

      function animate(timestamp) {
        window.requestAnimationFrame(animate);
        const now = Number(timestamp) || 0.0;
        const deltaSeconds = lastAnimationTimestamp === null
          ? 0.0
          : clampRange((now - lastAnimationTimestamp) / 1000.0, 0.0, 0.05);
        lastAnimationTimestamp = now;
        updateViewerActions(now);
        updateAnimatedFramePlayback(now);
        updateKeyboardMotion(deltaSeconds);
        updateGalacticSimpleDefaultOrbit(deltaSeconds);
        if (controls.enabled) {
          controls.update();
        } else if (cameraViewMode === "earth" && !cameraTransitionAnimationFrame && !skyViewDragState) {
          lockEarthViewCameraToTarget();
        }
        updateCameraResponsivePointSprites();
        updateScaleBar();
        updateSkyDome(now);
        updateSkyDomeBackgroundFrame(now);
        updateCameraResponsiveImagePlanes();
        updateScreenStableTextSprites();
        renderAgeKdeWidget();
        renderer.render(scene, camera);
      }

      applyInitialStateSync();
      buildAxes();
      initControls();
      initSkyPanel();
      updateSkyDomeCaptureFrame();
      initBoxMetricsPanel();
      initAgeKdePanel();
      initClusterFilterPanel();
      initDendrogramPanel();
      initVolumeControls();
      await restoreInitialLassoSelectionMask();
      renderLegend();
      updateSelectionUI();
      updatePlaybackButtons();
      renderFrame(currentFrameIndex);
      resize();
      setCameraAutoOrbitEnabled(cameraAutoOrbitEnabled);
      initialActionViewState = captureCurrentActionViewState();
      syncActionButtons();
      window.setTimeout(() => focusViewer(), 0);
      animate();

      canvas.addEventListener("pointerdown", onCanvasPointerDown, { capture: true });
      canvas.addEventListener("pointerenter", focusViewer);
      canvas.addEventListener("pointermove", onPointerMove);
      canvas.addEventListener("pointerleave", onPointerLeave);
      canvas.addEventListener("wheel", onCanvasWheel, { passive: false, capture: true });
      canvas.addEventListener("click", onCanvasClick);
      canvas.addEventListener("dblclick", onCanvasDoubleClick);
      scaleBarEl.addEventListener("pointerdown", onScaleBarPointerStart);
      window.addEventListener("keydown", onKeyDown);
      window.addEventListener("keyup", onKeyUp);
      window.addEventListener("blur", clearPressedKeys);
      window.addEventListener("pointerdown", onWindowPointerDown);
      window.addEventListener("pointermove", onWindowPointerMove);
      window.addEventListener("pointerup", onWindowPointerEnd);
      window.addEventListener("pointercancel", onWindowPointerEnd);
      window.addEventListener("resize", resize);
      window.addEventListener("message", onWindowMessage);
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

    def __init__(
        self,
        scene_spec: dict[str, Any],
        *,
        compress_scene_spec: bool | str | None = "auto",
        scene_spec_compression_threshold_bytes: int | None = None,
    ):
        self.scene_spec = scene_spec
        self._root_id = f"oviz-three-{uuid.uuid4().hex}"
        self.compress_scene_spec = compress_scene_spec
        self.scene_spec_compression_threshold_bytes = (
            DEFAULT_SCENE_SPEC_COMPRESSION_THRESHOLD_BYTES
            if scene_spec_compression_threshold_bytes is None
            else int(scene_spec_compression_threshold_bytes)
        )

    def to_dict(self) -> dict[str, Any]:
        return self.scene_spec

    def to_html(
        self,
        *,
        compress_scene_spec: bool | str | None = None,
        scene_spec_compression_threshold_bytes: int | None = None,
    ) -> str:
        return render_threejs_html(
            self.scene_spec,
            root_id=self._root_id,
            html_template=_THREEJS_HTML_TEMPLATE,
            topbar_html=_THREEJS_TOPBAR_HTML,
            minimal_topbar_html=_THREEJS_MINIMAL_TOPBAR_HTML,
            shell_html=THREEJS_SHELL_HTML,
            legend_runtime_js=THREEJS_LEGEND_RUNTIME_JS,
            widget_runtime_js=THREEJS_WIDGET_RUNTIME_JS,
            widget_content_runtime_js=THREEJS_WIDGET_CONTENT_RUNTIME_JS,
            interaction_runtime_js=THREEJS_INTERACTION_RUNTIME_JS,
            scene_runtime_js=THREEJS_SCENE_RUNTIME_JS,
            sky_runtime_js=THREEJS_SKY_RUNTIME_JS,
            viewer_runtime_js=THREEJS_VIEWER_RUNTIME_JS,
            action_runtime_js=THREEJS_ACTION_RUNTIME_JS,
            compress_scene_spec=(
                self.compress_scene_spec
                if compress_scene_spec is None
                else compress_scene_spec
            ),
            scene_spec_compression_threshold_bytes=(
                self.scene_spec_compression_threshold_bytes
                if scene_spec_compression_threshold_bytes is None
                else scene_spec_compression_threshold_bytes
            ),
        )

    def _data_url(self) -> str:
        return threejs_data_url(self.to_html())

    def _iframe_html(self) -> str:
        return threejs_iframe_html(self._data_url())

    def _repr_html_(self) -> str:
        return self._iframe_html()

    def write_html(self, file: str | Path, **kwargs: Any) -> None:
        compress_scene_spec = kwargs.pop("compress_scene_spec", None)
        scene_spec_compression_threshold_bytes = kwargs.pop(
            "scene_spec_compression_threshold_bytes",
            kwargs.pop("compress_scene_spec_threshold_bytes", None),
        )
        Path(file).write_text(
            self.to_html(
                compress_scene_spec=compress_scene_spec,
                scene_spec_compression_threshold_bytes=scene_spec_compression_threshold_bytes,
            ),
            encoding="utf-8",
        )

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
