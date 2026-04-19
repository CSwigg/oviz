from __future__ import annotations

import tempfile
import uuid
import webbrowser
from pathlib import Path
from typing import Any

from .threejs_embed import render_threejs_html, threejs_data_url, threejs_iframe_html
from .threejs_shell import THREEJS_SHELL_HTML
from .threejs_runtime_legend import THREEJS_LEGEND_RUNTIME_JS
from .threejs_runtime_interactions import THREEJS_INTERACTION_RUNTIME_JS
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
          <div class="oviz-three-tools-shell" data-open="false">
            <button class="oviz-three-tools-toggle" type="button" title="Show or hide selection controls">Selection ▸</button>
            <div class="oviz-three-tools-drawer">
              <div class="oviz-three-selection">
                <div class="oviz-three-selection-row">
                  <button class="oviz-three-selection-clear" type="button" title="Clear current cluster selection">Clear</button>
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
            </div>
          </div>
          <div class="oviz-three-controls-shell" data-open="false">
            <button class="oviz-three-controls-toggle" type="button" title="Show or hide the global scene controls">Controls ▸</button>
            <div class="oviz-three-controls-drawer">
                <div class="oviz-three-controls">
                  <div class="oviz-three-controls-title">Controls</div>
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
                    <input class="oviz-three-camera-fov" type="range" min="18" max="90" step="1" />
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
                  <button class="oviz-three-view-from-earth" type="button" title="Move the camera to the Earth position and look toward the Galactic center or active selection">View from Earth</button>
                  <button class="oviz-three-auto-orbit" type="button" title="Rotate around the current camera target" aria-pressed="false">Orbit camera</button>
                  <button class="oviz-three-reset-camera" type="button" title="Reset the camera to the initial view">Reset camera</button>
                  <button class="oviz-three-reset-controls" type="button" title="Reset the global control sliders">Reset controls</button>
                </div>
                <div class="oviz-three-controls-hint">Point size, glow, and opacity act as global multipliers on top of each trace's existing settings.</div>
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
        width: 100%;
        height: 100%;
        outline: none;
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
        display: none !important;
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
        left: 68%;
        bottom: clamp(84px, 16vh, 132px);
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
        gap: 4px;
        min-width: 0;
        color: var(--oviz-text);
        font-size: 12px;
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
        padding: 6px 8px;
        border-radius: 6px;
        border: 1px solid var(--oviz-panel-border);
        background: rgba(0, 0, 0, 0.18);
        color: var(--oviz-text);
        font: 12px Helvetica, Arial, sans-serif;
      }
      #__ROOT_ID__ .oviz-three-legend-field input[type="range"] {
        accent-color: var(--oviz-axis);
      }
      #__ROOT_ID__ .oviz-three-legend-field input[type="color"] {
        height: 34px;
        border: 1px solid var(--oviz-panel-border);
        border-radius: 6px;
        padding: 2px;
        background: rgba(0, 0, 0, 0.18);
      }
      #__ROOT_ID__ .oviz-three-legend-toggle {
        display: flex;
        align-items: center;
        gap: 8px;
        color: var(--oviz-text);
        font-size: 12px;
      }
      #__ROOT_ID__ .oviz-three-legend-toggle input {
        margin: 0;
        accent-color: var(--oviz-axis);
      }
      #__ROOT_ID__ .oviz-three-legend-summary {
        color: var(--oviz-text);
        font-size: 11px;
        line-height: 1.4;
        white-space: pre-wrap;
        opacity: 0.88;
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
      #__ROOT_ID__ .oviz-three-widget-menu {
        position: static;
        grid-column: 3;
        display: flex;
        align-items: center;
        justify-self: end;
        justify-content: flex-end;
        flex-wrap: wrap;
        gap: 4px;
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
        z-index: 6;
        width: min(220px, 19vw);
        min-width: 176px;
        min-height: 112px;
        display: flex;
        flex-direction: column;
        overflow: hidden;
      }
      #__ROOT_ID__ .oviz-three-legend-panel[data-open="false"] {
        min-height: 0;
      }
      #__ROOT_ID__ .oviz-three-legend-panel-head {
        display: flex;
        align-items: center;
        justify-content: space-between;
        gap: 10px;
        padding: 10px 12px;
        border-bottom: 1px solid rgba(255, 255, 255, 0.06);
        cursor: grab;
        user-select: none;
        touch-action: auto;
      }
      #__ROOT_ID__ .oviz-three-legend-panel[data-dragging="true"] .oviz-three-legend-panel-head {
        cursor: grabbing;
      }
      #__ROOT_ID__ .oviz-three-legend-panel-head {
        padding: 8px 10px;
        border-bottom-color: rgba(255, 255, 255, 0.05);
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
        padding: 6px 8px 8px;
        overflow: auto;
        transition: max-height 0.18s ease, padding 0.18s ease, opacity 0.18s ease;
      }
      #__ROOT_ID__ .oviz-three-legend-panel-body {
        padding: 5px 7px 7px;
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
        gap: 6px;
      }
      #__ROOT_ID__ .oviz-three-legend-section {
        display: flex;
        flex-direction: column;
        gap: 6px;
      }
      #__ROOT_ID__ .oviz-three-legend-section[data-open="false"] {
        gap: 0;
      }
      #__ROOT_ID__ .oviz-three-legend-section[data-empty="true"] {
        display: none;
      }
      #__ROOT_ID__ .oviz-three-legend-section-toggle {
        display: flex;
        align-items: center;
        justify-content: space-between;
        gap: 8px;
        width: 100%;
        padding: 1px 2px 0;
        border: 0;
        background: transparent;
        cursor: pointer;
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
        background: linear-gradient(180deg, rgba(255, 255, 255, 0.028), rgba(255, 255, 255, 0.016));
      }
      #__ROOT_ID__ .oviz-three-legend-volume-section {
        padding-top: 7px;
        border-top: 1px solid rgba(255, 255, 255, 0.05);
      }
      #__ROOT_ID__ .oviz-three-legend-section:last-child {
        padding-bottom: 2px;
      }
      #__ROOT_ID__ .oviz-three-legend-group-field {
        display: flex;
        flex-direction: column;
        gap: 4px;
        padding: 0 1px 8px;
        color: var(--oviz-muted-text);
        font: 500 10px/1.2 -apple-system, BlinkMacSystemFont, "SF Pro Text", "Helvetica Neue", sans-serif;
      }
      #__ROOT_ID__ .oviz-three-legend-group-field select {
        height: 30px;
        width: 100%;
      }
      #__ROOT_ID__ .oviz-three-group-select {
        height: 30px;
        width: 100%;
        padding: 0 28px 0 12px;
        border-radius: 4px;
        border: 1px solid rgba(255, 255, 255, 0.06);
        background: rgba(29, 32, 37, 0.98);
        color: var(--oviz-text);
        cursor: pointer;
        box-shadow: none;
        font: 500 11px/1 -apple-system, BlinkMacSystemFont, "SF Pro Text", "Helvetica Neue", sans-serif;
        -webkit-appearance: none;
        appearance: none;
      }
      #__ROOT_ID__ .oviz-three-group-select:hover {
        background: rgba(38, 42, 48, 0.98);
        border-color: rgba(255, 255, 255, 0.10);
      }
      #__ROOT_ID__ .oviz-three-legend-title {
        padding: 0 2px 6px;
        color: var(--oviz-muted-text);
        font: 500 10px/1.3 -apple-system, BlinkMacSystemFont, "SF Pro Text", "Helvetica Neue", sans-serif;
      }
      #__ROOT_ID__ .oviz-three-legend-entry {
        display: grid;
        grid-template-columns: minmax(0, 1fr) auto;
        align-items: center;
        gap: 6px;
        padding: 5px 7px;
        border: 1px solid rgba(255, 255, 255, 0.05);
        border-radius: 11px;
        background: rgba(255, 255, 255, 0.022);
        transition: background 0.14s ease, border-color 0.14s ease, opacity 0.14s ease;
      }
      #__ROOT_ID__ .oviz-three-legend-entry {
        gap: 4px;
        padding: 4px 6px;
        border-radius: 7px;
        border-color: rgba(255, 255, 255, 0.04);
        background: rgba(255, 255, 255, 0.01);
      }
      #__ROOT_ID__ .oviz-three-legend-entry[data-active="false"] {
        opacity: 0.46;
      }
      #__ROOT_ID__ .oviz-three-legend-item {
        flex: 1 1 auto;
        display: flex;
        align-items: center;
        gap: 8px;
        min-width: 0;
        border: 0;
        background: transparent;
        padding: 0;
        color: var(--oviz-text);
        text-align: left;
        font: 600 11px/1.2 -apple-system, BlinkMacSystemFont, "SF Pro Text", "Helvetica Neue", sans-serif;
      }
      #__ROOT_ID__ .oviz-three-legend-item:hover,
      #__ROOT_ID__ .oviz-three-legend-entry[data-editor-open="true"] {
        background: rgba(255, 255, 255, 0.05);
      }
      #__ROOT_ID__ .oviz-three-legend-item:hover,
      #__ROOT_ID__ .oviz-three-legend-entry[data-editor-open="true"] {
        background: rgba(255, 255, 255, 0.028);
      }
      #__ROOT_ID__ .oviz-three-legend-swatch {
        flex: 0 0 auto;
        width: 8px;
        height: 8px;
        border-radius: 999px;
        box-shadow: 0 0 0 1px rgba(255, 255, 255, 0.08);
      }
      #__ROOT_ID__ .oviz-three-legend-swatch {
        width: 7px;
        height: 7px;
        border-radius: 999px;
      }
      #__ROOT_ID__ .oviz-three-legend-meta {
        display: flex;
        flex-direction: column;
        min-width: 0;
      }
      #__ROOT_ID__ .oviz-three-legend-name {
        min-width: 0;
        overflow: hidden;
        text-overflow: ellipsis;
        white-space: nowrap;
      }
      #__ROOT_ID__ .oviz-three-legend-edit {
        flex: 0 0 auto;
        width: 22px;
        height: 22px;
        border: 0;
        border-radius: 999px;
        background: rgba(255, 255, 255, 0.05);
        color: var(--oviz-muted-text);
        font: 700 12px/1 -apple-system, BlinkMacSystemFont, "SF Pro Text", "Helvetica Neue", sans-serif;
      }
      #__ROOT_ID__ .oviz-three-legend-edit,
      #__ROOT_ID__ .oviz-three-legend-panel-toggle {
        border-radius: 6px;
        border: 1px solid rgba(255, 255, 255, 0.06);
        background: rgba(255, 255, 255, 0.024);
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
        pointer-events: auto;
        cursor: grab;
        user-select: none;
        touch-action: none;
        border-radius: 4px;
        background: rgba(18, 20, 24, 0.92);
        backdrop-filter: blur(4px) saturate(103%);
        -webkit-backdrop-filter: blur(4px) saturate(103%);
      }
      #__ROOT_ID__ .oviz-three-scale-bar[data-dragging="true"] {
        cursor: grabbing;
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
      #__ROOT_ID__[data-zen="true"] .oviz-three-legend-panel,
      #__ROOT_ID__[data-zen="true"] .oviz-three-key-help,
      #__ROOT_ID__[data-zen="true"] .oviz-three-widget-panel,
      #__ROOT_ID__[data-zen="true"] .oviz-three-note,
      #__ROOT_ID__[data-zen="true"] .oviz-three-scale-bar {
        display: none !important;
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
          justify-self: center;
          justify-content: center;
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
        #__ROOT_ID__ .oviz-three-scale-bar {
          display: none;
        }
      }
      @media (max-width: 1100px), (max-height: 760px) {
        #__ROOT_ID__[data-minimal="true"][data-galactic-simple="true"] .oviz-three-footer {
          left: 50%;
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
      __SHELL_HTML__
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

      /*__SCENE_SPEC_START__*/const sceneSpec = __SCENE_JSON__;/*__SCENE_SPEC_END__*/
      const root = document.getElementById("__ROOT_ID__");
      const initialState = sceneSpec.initial_state || {};
      const minimalModeEnabled = Boolean(initialState.minimal_mode_enabled);
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
      const titleEl = root.querySelector(".oviz-three-title");
      const zenModeButtonEl = root.querySelector(".oviz-three-zen-mode");
      const resetViewButtonEl = root.querySelector(".oviz-three-reset-view");
      const saveStateButtonEl = root.querySelector(".oviz-three-save-state");
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
      const sliderEl = root.querySelector(".oviz-three-slider");
      const sliderTrackWrapEl = root.querySelector(".oviz-three-slider-track-wrap");
      const sliderMinorTicksEl = root.querySelector(".oviz-three-slider-ticks-minor");
      const sliderMajorTicksEl = root.querySelector(".oviz-three-slider-ticks-major");
      const sliderLabelsEl = root.querySelector(".oviz-three-slider-labels");
      const timeLabelEl = root.querySelector(".oviz-three-time-label");
      const playBackwardButtonEl = root.querySelector(".oviz-three-play-backward");
      const playForwardButtonEl = root.querySelector(".oviz-three-play-forward");
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
      const legendItems = (sceneSpec.legend || {}).items || [];
      const groupVisibility = sceneSpec.group_visibility || {};
      const defaultGroup = sceneSpec.default_group || "All";
      const skySpec = sceneSpec.sky_panel || { enabled: false };
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
      const pointScale = Math.max(sceneSpec.max_span || 1, 1) / 4000.0;
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
      const starBloomMaterialCache = new Map();
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
      let lassoVolumeSelectionEnabled = false;
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
      let globalPointGlowStrength = 1.2;
      let sizePointsByStarsEnabled = false;
      let globalScrollSpeed = 1.0;
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
        alpha: false,
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
      scene.background = new THREE.Color(theme.scene_bgcolor || theme.paper_bgcolor || "#000000");

      const sceneUp = sceneSpec.camera_up || { x: 0.0, y: 0.0, z: 1.0 };
      const sceneUpVector = new THREE.Vector3(sceneUp.x ?? 0.0, sceneUp.y ?? 0.0, sceneUp.z ?? 1.0).normalize();
      const CAMERA_AUTO_ORBIT_SPEED = 1.2;
      const camera = new THREE.PerspectiveCamera(38, 1, 0.1, Math.max((sceneSpec.max_span || 1) * 20.0, 10000.0));
      camera.up.set(sceneUp.x ?? 0.0, sceneUp.y ?? 0.0, sceneUp.z ?? 1.0);
      const controls = new OrbitControls(camera, renderer.domElement);
      controls.enableDamping = true;
      controls.dampingFactor = 0.08;
      controls.rotateSpeed = 0.7;
      controls.panSpeed = 0.7;
      controls.zoomSpeed = globalScrollSpeed;
      controls.autoRotate = cameraAutoOrbitEnabled;
      controls.autoRotateSpeed = CAMERA_AUTO_ORBIT_SPEED;
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
      const initialCameraState = {
        position: camera.position.clone(),
        target: controls.target.clone(),
        up: camera.up.clone(),
        fov: Number(camera.fov),
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

      function sceneSpecMarkerWrappedJson(exportSceneSpec) {
        return `/*__SCENE_SPEC_START__*/const sceneSpec = ${JSON.stringify(exportSceneSpec)};/*__SCENE_SPEC_END__*/`;
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
        const computed = window.getComputedStyle(scaleBarEl);
        let left = Number.parseFloat(computed.left);
        let top = Number.parseFloat(computed.top);
        if (!Number.isFinite(top)) {
          const bottom = Number.parseFloat(computed.bottom);
          if (Number.isFinite(bottom)) {
            top = rootRect.height - size.height - bottom;
          }
        }
        if (!Number.isFinite(left)) {
          left = 18.0;
        }
        if (!Number.isFinite(top)) {
          top = Math.max(6.0, rootRect.height - size.height - 22.0);
        }
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
        if (minimalModeEnabled) {
          const visibleItems = visibleLegendItemsForCurrentGroup();
          const estimatedHeight = Math.max(visibleItems.length * 16 + 8, 24);
          const width = Math.min(Math.max(root.clientWidth * 0.18, 140.0), 240.0);
          const height = Math.min(estimatedHeight, Math.max(root.clientHeight - 24.0, 24.0));
          return clampLegendPanelRect(12.0, 12.0, width, height);
        }
        const width = Math.min(Math.max(root.clientWidth * 0.18, 204.0), 240.0);
        const visibleItems = visibleLegendItemsForCurrentGroup();
        const traceCount = visibleItems.filter((item) => !volumeLayerForKey(item.key)).length;
        const volumeCount = visibleItems.length - traceCount;
        const preferredTraceRows = traceCount ? clampRange(traceCount, 4, 8) : 0;
        const preferredVolumeRows = volumeCount ? clampRange(volumeCount, 1, 3) : 0;
        const estimatedHeight = 108 + preferredTraceRows * 34 + preferredVolumeRows * 32 + ((traceCount && volumeCount) ? 14 : 0);
        const minHeight = volumeCount ? 304.0 : 272.0;
        const height = Math.min(Math.max(estimatedHeight, minHeight), Math.min(root.clientHeight - 18.0, 520.0));
        return clampLegendPanelRect(14.0, 14.0, width, height);
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
          if (activeLegendEditorKey && legendEditButtonByKey.has(activeLegendEditorKey)) {
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
        const next = scaleBarPosition && Number.isFinite(Number(scaleBarPosition.left)) && Number.isFinite(Number(scaleBarPosition.top))
          ? clampScaleBarPosition(scaleBarPosition.left, scaleBarPosition.top, size.width, size.height)
          : defaultScaleBarPosition();
        scaleBarEl.style.left = `${next.left}px`;
        scaleBarEl.style.top = `${next.top}px`;
        scaleBarEl.style.right = "auto";
        scaleBarEl.style.bottom = "auto";
        if (scaleBarPosition) {
          scaleBarPosition = { left: next.left, top: next.top };
        }
      }

      function captureScaleBarState() {
        if (!scaleBarEl) {
          return null;
        }
        if (!scaleBarPosition) {
          const rootRect = root.getBoundingClientRect();
          const rect = scaleBarEl.getBoundingClientRect();
          if ((Number(rect.width) || 0.0) > 0.0 && (Number(rect.height) || 0.0) > 0.0) {
            scaleBarPosition = clampScaleBarPosition(
              rect.left - rootRect.left,
              rect.top - rootRect.top,
              rect.width,
              rect.height,
            );
          }
        }
        if (!scaleBarPosition) {
          return null;
        }
        return {
          left: Number(scaleBarPosition.left) || 0.0,
          top: Number(scaleBarPosition.top) || 0.0,
        };
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
          if (typeof savedGlobalControls.camera_auto_orbit_enabled === "boolean") {
            cameraAutoOrbitEnabled = savedGlobalControls.camera_auto_orbit_enabled;
          }
          if (typeof savedGlobalControls.camera_view_mode === "string" && savedGlobalControls.camera_view_mode) {
            cameraViewMode = String(savedGlobalControls.camera_view_mode);
          }
          if (Number.isFinite(Number(savedGlobalControls.earth_view_focus_distance_pc))) {
            earthViewFocusDistance = Number(savedGlobalControls.earth_view_focus_distance_pc);
          }
        }

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
          if (Number.isFinite(Number(target.x)) && Number.isFinite(Number(target.y)) && Number.isFinite(Number(target.z))) {
            controls.target.set(Number(target.x), Number(target.y), Number(target.z));
          }
          if (Number.isFinite(Number(position.x)) && Number.isFinite(Number(position.y)) && Number.isFinite(Number(position.z))) {
            camera.position.set(Number(position.x), Number(position.y), Number(position.z));
          }
          if (Number.isFinite(Number(up.x)) && Number.isFinite(Number(up.y)) && Number.isFinite(Number(up.z))) {
            camera.up.set(Number(up.x), Number(up.y), Number(up.z));
          }
          camera.lookAt(controls.target);
          controls.update();
          initialCameraState.position.copy(camera.position);
          initialCameraState.target.copy(controls.target);
          initialCameraState.up.copy(camera.up);
          initialCameraState.fov = Number(camera.fov);
        }

        const savedDrawers = initialState.drawers;
        if (savedDrawers && typeof savedDrawers === "object") {
          if (typeof savedDrawers.tools_open === "boolean") {
            setToolsDrawerOpen(savedDrawers.tools_open);
          }
          if (typeof savedDrawers.controls_open === "boolean") {
            setControlsDrawerOpen(savedDrawers.controls_open);
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
            camera_auto_orbit_enabled: cameraAutoOrbitEnabled,
            camera_view_mode: cameraViewMode,
            earth_view_focus_distance_pc: earthViewFocusDistance,
          },
          camera: {
            position: { x: camera.position.x, y: camera.position.y, z: camera.position.z },
            target: { x: controls.target.x, y: controls.target.y, z: controls.target.z },
            up: { x: camera.up.x, y: camera.up.y, z: camera.up.z },
          },
          drawers: {
            tools_open: Boolean(toolsShellEl && toolsShellEl.dataset.open === "true"),
            controls_open: Boolean(controlsShellEl && controlsShellEl.dataset.open === "true"),
          },
          legend_open: legendPanelOpen,
          legend_panel_state: captureLegendPanelState(),
          legend_panel_user_sized: legendPanelUserSized,
          legend_sections_open: safeJsonClone(legendSectionOpenState, { traces: true, volumes: true }),
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

      function buildExportHtml(exportSceneSpec) {
        const currentHtml = "<!DOCTYPE html>\\n" + document.documentElement.outerHTML;
        const startMarker = "/*__SCENE_SPEC_START__*/";
        const endMarker = "/*__SCENE_SPEC_END__*/";
        const startIndex = currentHtml.indexOf(startMarker);
        const endIndex = currentHtml.indexOf(endMarker);
        if (startIndex < 0 || endIndex < 0 || endIndex <= startIndex) {
          throw new Error("Could not locate scene export markers in the current figure HTML.");
        }
        return (
          currentHtml.slice(0, startIndex)
          + sceneSpecMarkerWrappedJson(exportSceneSpec)
          + currentHtml.slice(endIndex + endMarker.length)
        );
      }

      async function saveSceneStateToHtml() {
        const exportSceneSpec = safeJsonClone(sceneSpec, {});
        exportSceneSpec.initial_state = captureRuntimeState();
        exportSceneSpec.width = Math.max(root.clientWidth || sceneSpec.width || 900, 1);
        exportSceneSpec.height = Math.max(root.clientHeight || sceneSpec.height || 700, 1);
        const htmlText = buildExportHtml(exportSceneSpec);
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

      function buildAladinSrcdoc(selections, catalogPayload, mode = "overview", volumeOverlay = null) {
        const activeSelections = uniqueSelections(selections);
        const requestedMode = mode === "click" ? "click" : "overview";
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
      #aladin-lite-div { width: 100%; height: 100%; }
      .aladin-stack-box {
        max-height: min(56vh, 420px) !important;
        overflow-y: auto !important;
        overscroll-behavior: contain;
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
        const beamCenters = ${beamCentersJson};
        const statusEl = document.getElementById("oviz-status");
        const catalogs = [];
        const markerCanvasCache = new Map();
        let hoveredClusterKey = "";
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
            showReticle: viewMode === "click",
            showLayersControl: true,
            showGotoControl: true,
            showFrame: true
          };
          if (viewMode !== "click") {
            aladinOptions.projection = "MOL";
          }
          const aladin = A.aladin("#aladin-lite-div", aladinOptions);
          if (aladin && typeof aladin.setImageSurvey === "function") {
            aladin.setImageSurvey(${survey});
          }
          if (viewMode === "click" && aladin && typeof aladin.gotoRaDec === "function") {
            aladin.gotoRaDec(${ra.toFixed(8)}, ${dec.toFixed(8)});
          }
          if (viewMode !== "click" && aladin && typeof aladin.setProjection === "function") {
            aladin.setProjection("MOL");
          }
          if (viewMode !== "click" && aladin && typeof aladin.setFoV === "function") {
            aladin.setFoV(360.0);
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
              return;
            }
            setHoveredClusterKey(data.clusterKey, false);
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
        if (!skySpec.enabled || !skyFrameEl || event.source !== skyFrameEl.contentWindow) {
          return;
        }
        const data = event && event.data;
        if (!data || typeof data !== "object") {
          return;
        }
        if (data.type === "oviz-sky-hover-cluster") {
          setSkyHoveredClusterKey(data.clusterKey);
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
        state.steps = Math.round(Math.min(Math.max(Number(state.steps), 24.0), 768.0) / 8.0) * 8.0;
        if (!Number.isFinite(state.steps) || state.steps < 24) {
          state.steps = Number((layer.default_controls || {}).steps || 192);
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

      function volumeSupportsShowAllTimes(layer) {
        return Boolean(
          layer
          && layer.supports_show_all_times
          && !Number.isFinite(Number(layer.time_myr))
        );
      }

      function volumeVisibleForFrame(layer, state, frame = currentFrame()) {
        if (!layer || !state || state.visible === false) {
          return false;
        }
        const frameTime = frame ? Number(frame.time) : 0.0;
        const layerTime = Number(layer.time_myr);
        if (Number.isFinite(layerTime)) {
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

        volumeVMinEl.value = formatVolumeNumber(state.vmin);
        volumeVMaxEl.value = formatVolumeNumber(state.vmax);
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

      function starGlowTextureFor(kind = "realistic") {
        const cacheKey = String(kind || "realistic");
        if (starGlowTextureCache.has(cacheKey)) {
          return starGlowTextureCache.get(cacheKey);
        }

        const canvasEl = document.createElement("canvas");
        canvasEl.width = 256;
        canvasEl.height = 256;
        const ctx = canvasEl.getContext("2d");
        const cx = canvasEl.width * 0.5;
        const cy = canvasEl.height * 0.5;
        const size = canvasEl.width;
        ctx.clearRect(0, 0, canvasEl.width, canvasEl.height);

        if (cacheKey === "bloom") {
          const outerBloom = ctx.createRadialGradient(cx, cy, 0, cx, cy, size * 0.50);
          outerBloom.addColorStop(0.0, "rgba(255,255,255,0.38)");
          outerBloom.addColorStop(0.08, "rgba(255,255,255,0.26)");
          outerBloom.addColorStop(0.24, "rgba(255,255,255,0.12)");
          outerBloom.addColorStop(0.58, "rgba(255,255,255,0.035)");
          outerBloom.addColorStop(1.0, "rgba(255,255,255,0.0)");
          ctx.fillStyle = outerBloom;
          ctx.beginPath();
          ctx.arc(cx, cy, size * 0.50, 0, Math.PI * 2.0);
          ctx.fill();
        } else {
          const outerHalo = ctx.createRadialGradient(cx, cy, 0, cx, cy, size * 0.50);
          outerHalo.addColorStop(0.0, "rgba(255,255,255,0.92)");
          outerHalo.addColorStop(0.06, "rgba(255,255,255,0.84)");
          outerHalo.addColorStop(0.16, "rgba(255,255,255,0.50)");
          outerHalo.addColorStop(0.38, "rgba(255,255,255,0.18)");
          outerHalo.addColorStop(0.70, "rgba(255,255,255,0.05)");
          outerHalo.addColorStop(1.0, "rgba(255,255,255,0.0)");
          ctx.fillStyle = outerHalo;
          ctx.beginPath();
          ctx.arc(cx, cy, size * 0.50, 0, Math.PI * 2.0);
          ctx.fill();

          const innerHalo = ctx.createRadialGradient(cx, cy, 0, cx, cy, size * 0.24);
          innerHalo.addColorStop(0.0, "rgba(255,255,255,0.95)");
          innerHalo.addColorStop(0.10, "rgba(255,255,255,0.90)");
          innerHalo.addColorStop(0.26, "rgba(255,255,255,0.62)");
          innerHalo.addColorStop(0.54, "rgba(255,255,255,0.14)");
          innerHalo.addColorStop(1.0, "rgba(255,255,255,0.0)");
          ctx.fillStyle = innerHalo;
          ctx.beginPath();
          ctx.arc(cx, cy, size * 0.24, 0, Math.PI * 2.0);
          ctx.fill();
        }

        const texture = new THREE.CanvasTexture(canvasEl);
        texture.colorSpace = THREE.SRGBColorSpace;
        texture.generateMipmaps = false;
        texture.minFilter = THREE.LinearFilter;
        texture.magFilter = THREE.LinearFilter;
        texture.needsUpdate = true;
        starGlowTextureCache.set(cacheKey, texture);
        return texture;
      }

      function starCoreTextureFor(kind = "stellar_core") {
        const cacheKey = String(kind || "stellar_core");
        if (starCoreTextureCache.has(cacheKey)) {
          return starCoreTextureCache.get(cacheKey);
        }

        const canvasEl = document.createElement("canvas");
        canvasEl.width = 256;
        canvasEl.height = 256;
        const ctx = canvasEl.getContext("2d");
        const cx = canvasEl.width * 0.5;
        const cy = canvasEl.height * 0.5;
        const size = canvasEl.width;
        ctx.clearRect(0, 0, canvasEl.width, canvasEl.height);

        const core = ctx.createRadialGradient(cx, cy, 0, cx, cy, size * 0.22);
        core.addColorStop(0.0, "rgba(255,255,255,1.0)");
        core.addColorStop(0.08, "rgba(255,255,255,1.0)");
        core.addColorStop(0.18, "rgba(255,255,255,0.95)");
        core.addColorStop(0.38, "rgba(255,255,255,0.56)");
        core.addColorStop(0.72, "rgba(255,255,255,0.05)");
        core.addColorStop(1.0, "rgba(255,255,255,0.0)");
        ctx.fillStyle = core;
        ctx.beginPath();
        ctx.arc(cx, cy, size * 0.22, 0, Math.PI * 2.0);
        ctx.fill();

        const bloom = ctx.createRadialGradient(cx, cy, size * 0.02, cx, cy, size * 0.34);
        bloom.addColorStop(0.0, "rgba(255,255,255,0.45)");
        bloom.addColorStop(0.16, "rgba(255,255,255,0.28)");
        bloom.addColorStop(0.52, "rgba(255,255,255,0.07)");
        bloom.addColorStop(1.0, "rgba(255,255,255,0.0)");
        ctx.fillStyle = bloom;
        ctx.beginPath();
        ctx.arc(cx, cy, size * 0.34, 0, Math.PI * 2.0);
        ctx.fill();

        const texture = new THREE.CanvasTexture(canvasEl);
        texture.colorSpace = THREE.SRGBColorSpace;
        texture.generateMipmaps = false;
        texture.minFilter = THREE.LinearFilter;
        texture.magFilter = THREE.LinearFilter;
        texture.needsUpdate = true;
        starCoreTextureCache.set(cacheKey, texture);
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

      function starBloomMaterialFor(color, opacity) {
        const cacheKey = [color ?? "#ffffff", Number(opacity ?? 0.0).toFixed(4)].join("|");
        if (starBloomMaterialCache.has(cacheKey)) {
          return starBloomMaterialCache.get(cacheKey);
        }
        const material = new THREE.SpriteMaterial({
          map: starGlowTextureFor("bloom"),
          color: color ?? "#ffffff",
          transparent: true,
          opacity: opacity ?? 0.0,
          depthWrite: false,
          depthTest: true,
          blending: THREE.AdditiveBlending,
        });
        material.userData = { cached: true, glow: true, bloom: true };
        starBloomMaterialCache.set(cacheKey, material);
        return material;
      }

      function starCoreMaterialFor(color, opacity) {
        const cacheKey = [Number(opacity ?? 0.0).toFixed(4)].join("|");
        if (starCoreMaterialCache.has(cacheKey)) {
          return starCoreMaterialCache.get(cacheKey);
        }
        const material = new THREE.SpriteMaterial({
          map: starCoreTextureFor("stellar_core"),
          color: "#ffffff",
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

      function worldUnitsPerScreenPixelAt(position) {
        const point = position instanceof THREE.Vector3
          ? position
          : new THREE.Vector3(
            Number(position && position.x) || 0.0,
            Number(position && position.y) || 0.0,
            Number(position && position.z) || 0.0
          );
        const distance = Math.max(camera.position.distanceTo(point), 1e-3);
        const viewportHeight = Math.max(renderer.domElement.clientHeight || root.clientHeight || 1, 1);
        return (2.0 * distance * Math.tan(THREE.MathUtils.degToRad(camera.fov) * 0.5)) / viewportHeight;
      }

      function apparentPixelSpanForPoint(basePointScale, position) {
        const worldPerPixel = Math.max(worldUnitsPerScreenPixelAt(position), 1e-6);
        return Math.max(Number(basePointScale) || 0.0, 0.0) / worldPerPixel;
      }

      function screenPixelsToWorldScale(pixelSpan, position) {
        return Math.max(Number(pixelSpan) || 0.0, 0.0) * Math.max(worldUnitsPerScreenPixelAt(position), 1e-6);
      }

      // Optical point-spread functions are set primarily by the imaging system
      // (diffraction / seeing / detector sampling), not by the source's 3D distance.
      // So keep a persistent screen-space halo floor and only let apparent size gently
      // broaden the core and wings as a source becomes more resolved on screen.
      function glowScaleForPoint(basePointScale, position, glowStrength) {
        const strength = Math.max(Number(glowStrength) || 0.0, 0.0);
        const apparentPixels = apparentPixelSpanForPoint(basePointScale, position);
        const targetPixels = (7.0 + 4.5 * strength)
          + Math.min(apparentPixels * (0.32 + 0.05 * strength), 6.0 + 2.0 * strength);
        const clampedPixels = clampRange(targetPixels, 7.0 + 4.5 * strength, 18.0 + 8.0 * strength);
        return screenPixelsToWorldScale(clampedPixels, position);
      }

      function starCoreScaleForPoint(basePointScale, position, glowStrength) {
        const strength = Math.max(Number(glowStrength) || 0.0, 0.0);
        const apparentPixels = apparentPixelSpanForPoint(basePointScale, position);
        const targetPixels = Math.max(
          3.0 + 1.0 * strength,
          apparentPixels * (0.72 + 0.04 * strength)
        );
        const clampedPixels = clampRange(targetPixels, 3.0 + 0.8 * strength, 18.0 + 2.0 * strength);
        return screenPixelsToWorldScale(clampedPixels, position);
      }

      function starBloomScaleForPoint(basePointScale, position, glowStrength) {
        const strength = Math.max(Number(glowStrength) || 0.0, 0.0);
        const apparentPixels = apparentPixelSpanForPoint(basePointScale, position);
        const targetPixels = (14.0 + 8.0 * strength)
          + Math.min(apparentPixels * (0.48 + 0.08 * strength), 12.0 + 5.0 * strength);
        const clampedPixels = clampRange(targetPixels, 14.0 + 8.0 * strength, 34.0 + 14.0 * strength);
        return screenPixelsToWorldScale(clampedPixels, position);
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
        if (entry.scaleKind === "bloom") {
          return starBloomScaleForPoint(entry.pointScale, entry.position, globalPointGlowStrength);
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
          if (hideBelow > 0.0 && Number.isFinite(currentScaleBarLengthPc)) {
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
          const pointKey = normalizedSelectionKeyFor(point.selection);
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
          const effectiveOpacity = Math.min(1.0, Math.max(0.0, baseOpacity * traceOpacityMultiplier * opacityMultiplier * globalPointOpacityScale));
          if (effectiveOpacity <= 0.001) {
            return;
          }
          const scaleFloor = pointScale * 0.5 * Math.max(globalPointSizeScale, 0.05);
          const starsFactor = sizeByStarsFactorForPoint(point, trace, traceState);
          const scale = Math.max(
            pointState.size * sizeScaleFactor * starsFactor * globalPointSizeScale * pointScale,
            scaleFloor
          );
          const selectionKey = normalizedSelectionKeyFor(point.selection);
          const pointPosition = new THREE.Vector3(point.x, point.y, point.z);
          const glowStrength = Math.max(globalPointGlowStrength, 0.0);
          let interactionSprite = null;
          if (glowStrength > 0.02) {
            const bloomOpacity = clampRange(effectiveOpacity * (0.08 + 0.16 * glowStrength), 0.0, 0.82);
            const bloomSprite = new THREE.Sprite(starBloomMaterialFor(traceColor || point.color, bloomOpacity));
            const bloomScale = starBloomScaleForPoint(
              scale,
              pointPosition,
              glowStrength
            );
            bloomSprite.position.copy(pointPosition);
            bloomSprite.scale.set(bloomScale, bloomScale, 1.0);
            bloomSprite.renderOrder = -3;
            bloomSprite.userData = {
              selection: point.selection || null,
              selectionKey,
              baseScale: bloomScale,
              isBloom: true,
            };
            group.add(bloomSprite);
            registerCameraResponsivePointSprite(bloomSprite, "bloom", pointPosition, scale, selectionKey);

            const glowOpacity = clampRange(effectiveOpacity * (0.24 + 0.30 * glowStrength), 0.0, 0.98);
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

            const coreOpacity = clampRange(effectiveOpacity * (0.88 + 0.36 * glowStrength), 0.0, 1.0);
            const coreSprite = new THREE.Sprite(starCoreMaterialFor("#ffffff", coreOpacity));
            const coreScale = starCoreScaleForPoint(
              scale,
              pointPosition,
              glowStrength
            );
            coreSprite.position.copy(pointPosition);
            coreSprite.scale.set(coreScale, coreScale, 1.0);
            coreSprite.userData = {
              hovertext: point.hovertext || trace.name || "",
              selection: point.selection || null,
              selectionKey,
              baseScale: coreScale,
              isGlowCore: true,
            };
            group.add(coreSprite);
            hoverTargets.push(coreSprite);
            interactionSprite = coreSprite;
            registerCameraResponsivePointSprite(coreSprite, "core", pointPosition, scale, selectionKey);
          } else {
            const sprite = new THREE.Sprite(markerMaterialFor(point.symbol, traceColor || point.color, effectiveOpacity));
            sprite.position.copy(pointPosition);
            sprite.scale.set(scale, scale, scale);
            sprite.userData = {
              hovertext: point.hovertext || trace.name || "",
              selection: point.selection || null,
              selectionKey,
              baseScale: scale,
            };
            group.add(sprite);
            hoverTargets.push(sprite);
            interactionSprite = sprite;
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
        const opacityValue = traceState ? clamp01(traceState.opacity) : 1.0;
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
        if (isGalacticReferenceTrace(trace) && !galacticReferenceVisible) {
          return false;
        }
        if (isNearbyRegionLabelTrace(trace) && !nearbyRegionLabelsVisible) {
          return false;
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
        visibleText.textContent = Number.isFinite(Number(layer.time_myr)) ? "Show volume series" : "Show volume";
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
        vminInput.step = "any";
        const vminField = createLegendField("vmin", vminInput);

        const vmaxInput = document.createElement("input");
        vmaxInput.type = "number";
        vmaxInput.step = "any";
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
        samplesInput.step = "8";
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
            vminInput.value = formatVolumeNumber(state.vmin);
            vmaxInput.value = formatVolumeNumber(state.vmax);
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
        vminInput.addEventListener("change", () => {
          state.vmin = Number(vminInput.value);
          refreshVolumeControls(true);
          updateVolumeFromLegend(false);
        });
        vmaxInput.addEventListener("change", () => {
          state.vmax = Number(vmaxInput.value);
          refreshVolumeControls(true);
          updateVolumeFromLegend(false);
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
        setToolsDrawerOpen(Boolean(toolsShellEl && toolsShellEl.dataset.open === "true"));
        setControlsDrawerOpen(Boolean(controlsShellEl && controlsShellEl.dataset.open === "true"));
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
        groupSelectEl.value = currentGroup;
        groupSelectEl.style.display = groups.length > 1 ? "block" : "none";
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
          currentGroup = groupSelectEl.value;
          resetLegendState(currentGroup);
          renderLegend();
          renderFrame(currentFrameIndex);
        });
        renderSceneControls();
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
        sliderEl.addEventListener("input", () => {
          pause({ snap: false });
          scheduleSliderScrubRender(Number(sliderEl.value));
        });
        if (playBackwardButtonEl) {
          playBackwardButtonEl.addEventListener("click", () => {
            play(-1);
          });
        }
        if (playForwardButtonEl) {
          playForwardButtonEl.addEventListener("click", () => {
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
            viewFromEarth();
          });
        }
        orbitCameraButtons.forEach((buttonEl) => {
          buttonEl.addEventListener("click", () => {
            setCameraAutoOrbitEnabled(!cameraAutoOrbitEnabled);
            focusViewer();
          });
        });
        if (resetCameraButtonEl) {
          resetCameraButtonEl.addEventListener("click", () => {
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
            globalPointGlowStrength = 0.85;
            sizePointsByStarsEnabled = false;
            fadeInTimeMyr = Number(animationSpec.fade_in_time_default);
            fadeInAndOutEnabled = Boolean(animationSpec.fade_in_and_out_default);
            focusTraceKey = String(animationSpec.focus_trace_key_default || "");
            axesVisible = Boolean(sceneSpec.show_axes);
            galacticReferenceVisible = true;
            nearbyRegionLabelsVisible = true;
            setCameraAutoOrbitEnabled(false);
            cameraViewMode = "free";
            earthViewFocusDistance = null;
            camera.fov = Number(initialCameraState.fov);
            applyGlobalControlState();
            applyCameraViewMode();
            applyThemePreset(activeThemeKey, { rerender: false });
            renderSceneControls();
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
        focusViewer();
        if (galacticSimpleModeEnabled && galacticSimpleTracksOrbitTargetToSun) {
          enableGalacticSimpleOrbitTargetTracking();
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
        updateAnimatedFramePlayback(now);
        updateKeyboardMotion(deltaSeconds);
        controls.update();
        updateCameraResponsivePointSprites();
        updateScaleBar();
        updateCameraResponsiveImagePlanes();
        updateScreenStableTextSprites();
        renderAgeKdeWidget();
        renderer.render(scene, camera);
      }

      applyInitialStateSync();
      buildAxes();
      initControls();
      initSkyPanel();
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
      window.setTimeout(() => focusViewer(), 0);
      animate();

      canvas.addEventListener("pointerdown", onCanvasPointerDown);
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

    def __init__(self, scene_spec: dict[str, Any]):
        self.scene_spec = scene_spec
        self._root_id = f"oviz-three-{uuid.uuid4().hex}"

    def to_dict(self) -> dict[str, Any]:
        return self.scene_spec

    def to_html(self) -> str:
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
        )

    def _data_url(self) -> str:
        return threejs_data_url(self.to_html())

    def _iframe_html(self) -> str:
        return threejs_iframe_html(self._data_url())

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
