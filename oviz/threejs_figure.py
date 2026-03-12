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
      #__ROOT_ID__:focus {
        outline: none;
      }
      #__ROOT_ID__ .oviz-three-canvas {
        display: block;
        width: 100%;
        height: 100%;
        outline: none;
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
      #__ROOT_ID__ .oviz-three-widget-menu {
        position: absolute;
        top: 12px;
        right: 12px;
        z-index: 6;
        display: flex;
        align-items: center;
        gap: 8px;
      }
      #__ROOT_ID__ .oviz-three-tools-shell,
      #__ROOT_ID__ .oviz-three-controls-shell {
        position: relative;
        min-height: 34px;
      }
      #__ROOT_ID__ .oviz-three-tools-toggle,
      #__ROOT_ID__ .oviz-three-controls-toggle {
        border: 1px solid var(--oviz-panel-border);
        border-radius: 8px;
        background: var(--oviz-panel-bg);
        color: var(--oviz-text);
        cursor: pointer;
        font: 12px Helvetica, Arial, sans-serif;
        padding: 7px 11px;
      }
      #__ROOT_ID__ .oviz-three-tools-drawer,
      #__ROOT_ID__ .oviz-three-controls-drawer {
        position: absolute;
        top: 0;
        left: 52px;
        display: flex;
        flex-direction: column;
        gap: 10px;
        width: min(320px, 32vw);
        transform: translateX(-20px);
        opacity: 0;
        pointer-events: none;
        transition: transform 0.18s ease, opacity 0.18s ease;
      }
      #__ROOT_ID__ .oviz-three-tools-shell[data-open="true"] .oviz-three-tools-drawer,
      #__ROOT_ID__ .oviz-three-controls-shell[data-open="true"] .oviz-three-controls-drawer {
        transform: translateX(0);
        opacity: 1;
        pointer-events: auto;
      }
      #__ROOT_ID__ .oviz-three-group-select,
      #__ROOT_ID__ .oviz-three-widget-select {
        min-width: 180px;
        padding: 6px 10px;
        border-radius: 6px;
        border: 1px solid var(--oviz-panel-border);
        background: var(--oviz-panel-bg);
        color: var(--oviz-text);
        font-size: 14px;
        font-family: inherit;
      }
      #__ROOT_ID__ .oviz-three-save-state {
        padding: 7px 11px;
        border-radius: 8px;
        border: 1px solid var(--oviz-panel-border);
        background: var(--oviz-panel-bg);
        color: var(--oviz-text);
        cursor: pointer;
        font: 12px Helvetica, Arial, sans-serif;
      }
      #__ROOT_ID__ .oviz-three-zen-mode {
        padding: 7px 11px;
        border-radius: 8px;
        border: 1px solid var(--oviz-panel-border);
        background: var(--oviz-panel-bg);
        color: var(--oviz-text);
        cursor: pointer;
        font: 12px Helvetica, Arial, sans-serif;
      }
      #__ROOT_ID__ .oviz-three-reset-view {
        padding: 7px 11px;
        border-radius: 8px;
        border: 1px solid var(--oviz-panel-border);
        background: var(--oviz-panel-bg);
        color: var(--oviz-text);
        cursor: pointer;
        font: 12px Helvetica, Arial, sans-serif;
      }
      #__ROOT_ID__ .oviz-three-zen-mode[data-active="true"] {
        border-color: var(--oviz-axis);
        box-shadow: inset 0 0 0 1px rgba(255, 255, 255, 0.10);
      }
      #__ROOT_ID__ .oviz-three-save-state:disabled {
        opacity: 0.55;
        cursor: default;
      }
      #__ROOT_ID__[data-zen="true"] .oviz-three-toolbar,
      #__ROOT_ID__[data-zen="true"] .oviz-three-widget-select,
      #__ROOT_ID__[data-zen="true"] .oviz-three-reset-view,
      #__ROOT_ID__[data-zen="true"] .oviz-three-save-state,
      #__ROOT_ID__[data-zen="true"] .oviz-three-key-help,
      #__ROOT_ID__[data-zen="true"] .oviz-three-widget-panel {
        display: none !important;
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
      #__ROOT_ID__ .oviz-three-legend-entry:last-child {
        padding-bottom: 0;
        border-bottom: 0;
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
      #__ROOT_ID__ .oviz-three-controls-field select {
        width: 100%;
        min-width: 0;
        box-sizing: border-box;
      }
      #__ROOT_ID__ .oviz-three-controls-field select {
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
      #__ROOT_ID__ .oviz-three-controls-hint {
        color: var(--oviz-text);
        font-size: 11px;
        line-height: 1.4;
        white-space: pre-wrap;
        opacity: 0.88;
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
        display: none !important;
      }
      #__ROOT_ID__ .oviz-three-volume[data-enabled="true"] {
        display: none !important;
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
      #__ROOT_ID__ .oviz-three-scale-bar {
        position: absolute;
        left: 18px;
        bottom: 86px;
        z-index: 6;
        display: flex;
        flex-direction: column;
        align-items: flex-start;
        gap: 6px;
        pointer-events: none;
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
      #__ROOT_ID__ .oviz-three-widget-drag {
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
      #__ROOT_ID__ .oviz-three-widget-panel[data-mode="fullscreen"] .oviz-three-widget-drag {
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
      #__ROOT_ID__ .oviz-three-age-body {
        position: absolute;
        top: 28px;
        right: 0;
        bottom: 82px;
        left: 0;
        display: flex;
        flex-direction: column;
        gap: 8px;
        padding: 10px 10px 8px 10px;
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
      #__ROOT_ID__ .oviz-three-age-readout {
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
        pointer-events: none;
      }
      #__ROOT_ID__ .oviz-three-filter-body {
        position: absolute;
        top: 28px;
        right: 0;
        bottom: 0;
        left: 0;
        display: flex;
        flex-direction: column;
        gap: 8px;
        padding: 10px 10px 72px 10px;
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
      #__ROOT_ID__ .oviz-three-filter-readout {
        position: absolute;
        right: 0;
        bottom: 0;
        left: 0;
        height: 64px;
        margin: 0;
        padding: 8px 10px;
        overflow: auto;
        border-top: 1px solid var(--oviz-panel-border);
        background: rgba(0, 0, 0, 0.18);
        color: var(--oviz-text);
        font: 11px/1.45 Menlo, Monaco, Consolas, monospace;
        white-space: pre-wrap;
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
      #__ROOT_ID__ .oviz-three-sky-full,
      #__ROOT_ID__ .oviz-three-sky-hide,
      #__ROOT_ID__ .oviz-three-age-full,
      #__ROOT_ID__ .oviz-three-age-hide,
      #__ROOT_ID__ .oviz-three-filter-full,
      #__ROOT_ID__ .oviz-three-filter-hide {
        position: absolute;
        top: 6px;
        z-index: 3;
      }
      #__ROOT_ID__ .oviz-three-sky-full,
      #__ROOT_ID__ .oviz-three-age-full,
      #__ROOT_ID__ .oviz-three-filter-full {
        right: 52px;
      }
      #__ROOT_ID__ .oviz-three-sky-hide,
      #__ROOT_ID__ .oviz-three-age-hide,
      #__ROOT_ID__ .oviz-three-filter-hide {
        right: 6px;
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
    </style>
  </head>
  <body>
    <div id="__ROOT_ID__" tabindex="0" data-zen="false">
      <div class="oviz-three-title"></div>
      <div class="oviz-three-widget-menu">
        <button class="oviz-three-zen-mode" type="button" title="Hide interface panels and keep only the time slider visible">Zen</button>
        <button class="oviz-three-reset-view" type="button" title="Reset the camera and clear current lasso and click selections">Reset</button>
        <button class="oviz-three-save-state" type="button" title="Export an HTML copy of the figure with the current state">Save State</button>
        <select class="oviz-three-widget-select">
          <option value="">Widgets</option>
        </select>
      </div>
      <div class="oviz-three-toolbar">
        <select class="oviz-three-group-select"></select>
        <div class="oviz-three-legend"></div>
        <div class="oviz-three-tools-shell" data-open="false">
          <button class="oviz-three-tools-toggle" type="button" title="Show or hide the left control drawer">Tools ▸</button>
          <div class="oviz-three-tools-drawer">
            <div class="oviz-three-selection">
              <div class="oviz-three-selection-row">
                <button class="oviz-three-lasso-button" type="button" title="Lasso clusters on the 3D plot">Lasso</button>
                <button class="oviz-three-selection-clear" type="button" title="Clear current cluster selection">Clear</button>
              </div>
              <label class="oviz-three-selection-toggle">
                <input class="oviz-three-click-select-toggle" type="checkbox" />
                <span>Enable click select</span>
              </label>
              <label class="oviz-three-selection-toggle">
                <input class="oviz-three-volume-lasso-toggle" type="checkbox" checked />
                <span>Lasso volumetric data</span>
              </label>
              <div class="oviz-three-selection-readout">Shift+drag or use Lasso to select clusters on the current frame.</div>
            </div>
            <div class="oviz-three-volume" data-enabled="false">
              <div class="oviz-three-volume-title">Volume</div>
              <label class="oviz-three-volume-field">
                <span>Layer</span>
                <select class="oviz-three-volume-select"></select>
              </label>
              <label class="oviz-three-volume-toggle">
                <input class="oviz-three-volume-visible" type="checkbox" />
                <span>Show at t=0</span>
              </label>
              <label class="oviz-three-volume-field">
                <span>Colormap</span>
                <select class="oviz-three-volume-colormap"></select>
              </label>
              <label class="oviz-three-volume-field">
                <span>Stretch</span>
                <select class="oviz-three-volume-stretch"></select>
              </label>
              <div class="oviz-three-volume-row">
                <label class="oviz-three-volume-field">
                  <span>vmin</span>
                  <input class="oviz-three-volume-vmin" type="number" step="any" />
                </label>
                <label class="oviz-three-volume-field">
                  <span>vmax</span>
                  <input class="oviz-three-volume-vmax" type="number" step="any" />
                </label>
              </div>
              <label class="oviz-three-volume-field">
                <span class="oviz-three-volume-opacity-label">Opacity</span>
                <input class="oviz-three-volume-opacity" type="range" min="0" max="1" step="0.01" />
              </label>
              <label class="oviz-three-volume-field">
                <span class="oviz-three-volume-alpha-label">Alpha coef</span>
                <input class="oviz-three-volume-alpha" type="range" min="1" max="200" step="1" />
              </label>
              <label class="oviz-three-volume-field">
                <span class="oviz-three-volume-steps-label">Samples</span>
                <input class="oviz-three-volume-steps" type="range" min="24" max="768" step="8" />
              </label>
              <div class="oviz-three-volume-summary"></div>
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
                <input class="oviz-three-axes-visible-toggle" type="checkbox" checked />
                <span>Show axes</span>
              </label>
              <div class="oviz-three-controls-actions">
                <button class="oviz-three-key-help-button" type="button" title="Show keyboard controls">Keyboard help</button>
                <button class="oviz-three-view-from-earth" type="button" title="Move the camera to the Earth position and look toward the Galactic center or active selection">View from Earth</button>
                <button class="oviz-three-reset-camera" type="button" title="Reset the camera to the initial view">Reset camera</button>
                <button class="oviz-three-reset-controls" type="button" title="Reset the global control sliders">Reset controls</button>
              </div>
              <div class="oviz-three-controls-hint">Point size and opacity act as multiplicative global factors on top of each trace's existing settings.</div>
            </div>
          </div>
        </div>
      </div>
      <div class="oviz-three-key-help" data-open="false">
        <div class="oviz-three-key-help-head">
          <div class="oviz-three-key-help-title">Keyboard Controls</div>
          <button class="oviz-three-key-help-close" type="button" title="Close keyboard help">Close</button>
        </div>
        <p class="oviz-three-key-help-text">Keyboard controls are active as soon as the viewer loads. Shift switches WASD from orbiting to panning and makes frame stepping larger.</p>
        <div class="oviz-three-key-help-grid">
          <div class="oviz-three-key-help-keys">W / A / S / D</div>
          <div>Orbit around the current anchor point.</div>
          <div class="oviz-three-key-help-keys">Shift + W / A / S / D</div>
          <div>Pan in camera-relative space while keeping the current view direction.</div>
          <div class="oviz-three-key-help-keys">R / F</div>
          <div>Move up and down along the z-up axis.</div>
          <div class="oviz-three-key-help-keys">Q / E</div>
          <div>Zoom out and in toward the current anchor point.</div>
          <div class="oviz-three-key-help-keys">Left / Right</div>
          <div>Step the time slider backward or forward.</div>
          <div class="oviz-three-key-help-keys">Shift + Left / Right</div>
          <div>Jump five time steps backward or forward.</div>
          <div class="oviz-three-key-help-keys">1 - 9</div>
          <div>Toggle the first nine legend items in the current group.</div>
          <div class="oviz-three-key-help-keys">Shift + 1 - 9</div>
          <div>Solo the corresponding legend item.</div>
          <div class="oviz-three-key-help-keys">Space</div>
          <div>Play or pause the animation.</div>
          <div class="oviz-three-key-help-keys">L / C</div>
          <div>Toggle lasso mode or click selection.</div>
          <div class="oviz-three-key-help-keys">V</div>
          <div>Activate View from Earth.</div>
          <div class="oviz-three-key-help-keys">?</div>
          <div>Open this keyboard help panel.</div>
          <div class="oviz-three-key-help-keys">Esc</div>
          <div>Close help or clear the current selection.</div>
        </div>
      </div>
      <canvas class="oviz-three-canvas" tabindex="0"></canvas>
      <div class="oviz-three-lasso-overlay" data-active="false">
        <svg preserveAspectRatio="none" aria-hidden="true">
          <polyline points=""></polyline>
        </svg>
      </div>
      <div class="oviz-three-tooltip"></div>
      <div class="oviz-three-scale-bar">
        <div class="oviz-three-scale-label"></div>
        <div class="oviz-three-scale-line"></div>
      </div>
      <div class="oviz-three-footer">
        <button class="oviz-three-play" type="button" title="Play">▶</button>
        <button class="oviz-three-pause" type="button" title="Pause">⏸</button>
        <span class="oviz-three-time-label"></span>
        <input class="oviz-three-slider" type="range" min="0" max="0" step="1" value="0" />
      </div>
      <div class="oviz-three-note"></div>
      <div class="oviz-three-sky-panel oviz-three-widget-panel" data-widget-key="sky" data-mode="hidden">
        <div class="oviz-three-sky-drag oviz-three-widget-drag">Sky</div>
        <iframe class="oviz-three-sky-frame" loading="eager" referrerpolicy="no-referrer"></iframe>
        <pre class="oviz-three-sky-readout"></pre>
        <button class="oviz-three-sky-full" type="button" title="Toggle fullscreen sky panel">Full</button>
        <button class="oviz-three-sky-hide" type="button" title="Hide sky panel">Hide</button>
        <div class="oviz-three-sky-resize oviz-three-widget-resize" data-dir="nw"></div>
        <div class="oviz-three-sky-resize oviz-three-widget-resize" data-dir="ne"></div>
        <div class="oviz-three-sky-resize oviz-three-widget-resize" data-dir="sw"></div>
        <div class="oviz-three-sky-resize oviz-three-widget-resize" data-dir="se"></div>
      </div>
      <div class="oviz-three-age-panel oviz-three-widget-panel" data-widget-key="age_kde" data-mode="hidden">
        <div class="oviz-three-age-drag oviz-three-widget-drag">Age KDE</div>
        <div class="oviz-three-age-body">
          <div class="oviz-three-age-kde-shell">
            <canvas class="oviz-three-age-canvas"></canvas>
          </div>
          <div class="oviz-three-age-filter">
            <div class="oviz-three-age-filter-slider-shell">
              <input class="oviz-three-age-filter-range oviz-three-age-filter-range-min" type="range" min="0" max="1000" step="1" value="0" />
              <input class="oviz-three-age-filter-range oviz-three-age-filter-range-max" type="range" min="0" max="1000" step="1" value="1000" />
            </div>
            <div class="oviz-three-age-filter-labels">
              <span class="oviz-three-age-filter-range-readout-min"></span>
              <span class="oviz-three-age-filter-range-readout-max"></span>
            </div>
            <div class="oviz-three-age-filter-hint">Age filter for clusters contributing to the current KDE view</div>
          </div>
        </div>
        <pre class="oviz-three-age-readout"></pre>
        <button class="oviz-three-age-full" type="button" title="Toggle fullscreen age KDE panel">Full</button>
        <button class="oviz-three-age-hide" type="button" title="Hide age KDE panel">Hide</button>
        <div class="oviz-three-age-resize oviz-three-widget-resize" data-dir="nw"></div>
        <div class="oviz-three-age-resize oviz-three-widget-resize" data-dir="ne"></div>
        <div class="oviz-three-age-resize oviz-three-widget-resize" data-dir="sw"></div>
        <div class="oviz-three-age-resize oviz-three-widget-resize" data-dir="se"></div>
      </div>
      <div class="oviz-three-filter-panel oviz-three-widget-panel" data-widget-key="cluster_filter" data-mode="hidden">
        <div class="oviz-three-filter-drag oviz-three-widget-drag">Filter</div>
        <div class="oviz-three-filter-body">
          <div class="oviz-three-filter-row">
            <label class="oviz-three-filter-field">
              <span>Parameter</span>
              <select class="oviz-three-filter-parameter"></select>
            </label>
          </div>
          <div class="oviz-three-filter-hist">
            <canvas class="oviz-three-filter-canvas"></canvas>
          </div>
          <div class="oviz-three-filter-slider-shell">
            <input class="oviz-three-filter-range oviz-three-filter-range-min" type="range" min="0" max="1000" step="1" value="0" />
            <input class="oviz-three-filter-range oviz-three-filter-range-max" type="range" min="0" max="1000" step="1" value="1000" />
          </div>
          <div class="oviz-three-filter-labels">
            <span class="oviz-three-filter-range-readout-min"></span>
            <span class="oviz-three-filter-range-readout-max"></span>
          </div>
        </div>
        <pre class="oviz-three-filter-readout"></pre>
        <button class="oviz-three-filter-full" type="button" title="Toggle fullscreen filter panel">Full</button>
        <button class="oviz-three-filter-hide" type="button" title="Hide filter panel">Hide</button>
        <div class="oviz-three-filter-resize oviz-three-widget-resize" data-dir="nw"></div>
        <div class="oviz-three-filter-resize oviz-three-widget-resize" data-dir="ne"></div>
        <div class="oviz-three-filter-resize oviz-three-widget-resize" data-dir="sw"></div>
        <div class="oviz-three-filter-resize oviz-three-widget-resize" data-dir="se"></div>
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

      /*__SCENE_SPEC_START__*/const sceneSpec = __SCENE_JSON__;/*__SCENE_SPEC_END__*/
      const initialState = sceneSpec.initial_state || {};
      const root = document.getElementById("__ROOT_ID__");
      const canvas = root.querySelector(".oviz-three-canvas");
      const titleEl = root.querySelector(".oviz-three-title");
      const zenModeButtonEl = root.querySelector(".oviz-three-zen-mode");
      const resetViewButtonEl = root.querySelector(".oviz-three-reset-view");
      const saveStateButtonEl = root.querySelector(".oviz-three-save-state");
      const groupSelectEl = root.querySelector(".oviz-three-group-select");
      const widgetSelectEl = root.querySelector(".oviz-three-widget-select");
      const legendEl = root.querySelector(".oviz-three-legend");
      const keyHelpEl = root.querySelector(".oviz-three-key-help");
      const keyHelpButtonEl = root.querySelector(".oviz-three-key-help-button");
      const keyHelpCloseEl = root.querySelector(".oviz-three-key-help-close");
      const toolsShellEl = root.querySelector(".oviz-three-tools-shell");
      const toolsToggleEl = root.querySelector(".oviz-three-tools-toggle");
      const controlsShellEl = root.querySelector(".oviz-three-controls-shell");
      const controlsToggleEl = root.querySelector(".oviz-three-controls-toggle");
      const sliderEl = root.querySelector(".oviz-three-slider");
      const timeLabelEl = root.querySelector(".oviz-three-time-label");
      const playButtonEl = root.querySelector(".oviz-three-play");
      const pauseButtonEl = root.querySelector(".oviz-three-pause");
      const tooltipEl = root.querySelector(".oviz-three-tooltip");
      const scaleBarEl = root.querySelector(".oviz-three-scale-bar");
      const scaleLabelEl = root.querySelector(".oviz-three-scale-label");
      const noteEl = root.querySelector(".oviz-three-note");
      const selectionReadoutEl = root.querySelector(".oviz-three-selection-readout");
      const lassoButtonEl = root.querySelector(".oviz-three-lasso-button");
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
      const focusGroupSelectEl = root.querySelector(".oviz-three-focus-group-select");
      const fadeTimeEl = root.querySelector(".oviz-three-fade-time");
      const fadeInOutToggleEl = root.querySelector(".oviz-three-fade-in-out-toggle");
      const axesVisibleToggleEl = root.querySelector(".oviz-three-axes-visible-toggle");
      const viewFromEarthButtonEl = root.querySelector(".oviz-three-view-from-earth");
      const resetCameraButtonEl = root.querySelector(".oviz-three-reset-camera");
      const resetControlsButtonEl = root.querySelector(".oviz-three-reset-controls");
      const volumePanelEl = root.querySelector(".oviz-three-volume");
      const volumeSelectEl = root.querySelector(".oviz-three-volume-select");
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
      const skyReadoutEl = root.querySelector(".oviz-three-sky-readout");
      const skyFullButtonEl = root.querySelector(".oviz-three-sky-full");
      const skyHideButtonEl = root.querySelector(".oviz-three-sky-hide");
      const ageKdePanelEl = root.querySelector(".oviz-three-age-panel");
      const ageKdeCanvasEl = root.querySelector(".oviz-three-age-canvas");
      const ageKdeReadoutEl = root.querySelector(".oviz-three-age-readout");
      const ageKdeFilterRangeMinEl = root.querySelector(".oviz-three-age-filter-range-min");
      const ageKdeFilterRangeMaxEl = root.querySelector(".oviz-three-age-filter-range-max");
      const ageKdeFilterRangeReadoutMinEl = root.querySelector(".oviz-three-age-filter-range-readout-min");
      const ageKdeFilterRangeReadoutMaxEl = root.querySelector(".oviz-three-age-filter-range-readout-max");
      const ageKdeFullButtonEl = root.querySelector(".oviz-three-age-full");
      const ageKdeHideButtonEl = root.querySelector(".oviz-three-age-hide");
      const clusterFilterPanelEl = root.querySelector(".oviz-three-filter-panel");
      const clusterFilterCanvasEl = root.querySelector(".oviz-three-filter-canvas");
      const clusterFilterReadoutEl = root.querySelector(".oviz-three-filter-readout");
      const clusterFilterParameterEl = root.querySelector(".oviz-three-filter-parameter");
      const clusterFilterRangeMinEl = root.querySelector(".oviz-three-filter-range-min");
      const clusterFilterRangeMaxEl = root.querySelector(".oviz-three-filter-range-max");
      const clusterFilterRangeReadoutMinEl = root.querySelector(".oviz-three-filter-range-readout-min");
      const clusterFilterRangeReadoutMaxEl = root.querySelector(".oviz-three-filter-range-readout-max");
      const clusterFilterFullButtonEl = root.querySelector(".oviz-three-filter-full");
      const clusterFilterHideButtonEl = root.querySelector(".oviz-three-filter-hide");
      const widgetPanels = Array.from(root.querySelectorAll(".oviz-three-widget-panel"));
      const widgetDragHandles = Array.from(root.querySelectorAll(".oviz-three-widget-drag"));
      const widgetResizeEls = Array.from(root.querySelectorAll(".oviz-three-widget-resize"));

      const baseTheme = safeJsonClone(sceneSpec.theme || {}, {});
      const theme = safeJsonClone(baseTheme, {});
      titleEl.textContent = sceneSpec.title || "";

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
        };
      }

      function applyThemeCssVars() {
        root.style.setProperty("--oviz-paper-bg", theme.paper_bgcolor || "#000000");
        root.style.setProperty("--oviz-scene-bg", theme.scene_bgcolor || theme.paper_bgcolor || "#000000");
        root.style.setProperty("--oviz-text", theme.text_color || "#d0d0d0");
        root.style.setProperty("--oviz-axis", theme.axis_color || "#808080");
        root.style.setProperty("--oviz-panel-bg", theme.panel_bg || "rgba(0, 0, 0, 0.45)");
        root.style.setProperty("--oviz-panel-border", theme.panel_border || "rgba(128, 128, 128, 0.50)");
        root.style.setProperty("--oviz-panel-solid", theme.panel_solid || theme.paper_bgcolor || "#121212");
        root.style.setProperty("--oviz-footprint", theme.footprint || "#6ec5ff");
      }

      applyThemeCssVars();
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
      const ageKdeSpec = sceneSpec.age_kde || { enabled: false };
      const clusterFilterSpec = sceneSpec.cluster_filter || { enabled: false };
      const volumeLayers = ((sceneSpec.volumes || {}).layers || []);
      const volumeLayersByKey = new Map(volumeLayers.map((layer) => [String(layer.key), layer]));
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
      const themePresets = buildThemePresets(baseTheme);
      const pointScale = Math.max(sceneSpec.max_span || 1, 1) / 4000.0;
      const hoverTargets = [];
      const selectionSpriteEntriesByKey = new Map();
      const axisLineMaterials = [];
      const frameLineMaterials = [];
      const screenStableTextSprites = [];
      const markerTextureCache = new Map();
      const markerMaterialCache = new Map();
      const textTextureCache = new Map();
      const galaxyTextureCache = new Map();
      const volumeScalarDataCache = new Map();
      const volumeScalarDataPendingCache = new Map();
      const volumeTextureCache = new Map();
      const volumeColorTextureCache = new Map();
      const volumeColorBytesCache = new Map();
      const volumeRuntimeByKey = new Map();
      let volumeJitterTexture = null;
      let legendState = {};
      let currentGroup = defaultGroup;
      let currentFrameIndex = sceneSpec.initial_frame_index || 0;
      let currentSelection = null;
      let currentSelections = [];
      let currentSelectionMode = "none";
      let clickSelectionEnabled = false;
      let lassoVolumeSelectionEnabled = true;
      let playbackTimer = null;
      let skyPanelMode = "hidden";
      let ageKdePanelMode = "hidden";
      let clusterFilterPanelMode = "hidden";
      let widgetPointerState = null;
      let widgetZIndexCounter = 8;
      let lassoState = null;
      let lassoArmed = false;
      let suppressNextCanvasClick = false;
      let selectedClusterKeys = new Set();
      let localHoveredClusterKey = "";
      let skyHoveredClusterKey = "";
      let lastSentSkyHoverClusterKey = null;
      let activeThemeKey = "default";
      let axesVisible = Boolean(sceneSpec.show_axes);
      let globalPointSizeScale = 1.0;
      let globalPointOpacityScale = 1.0;
      let globalScrollSpeed = 1.0;
      let zenModeEnabled = Boolean(initialState.zen_mode_enabled);
      let fadeInTimeMyr = Number(animationSpec.fade_in_time_default);
      let fadeInAndOutEnabled = Boolean(animationSpec.fade_in_and_out_default);
      let focusTraceKey = String(animationSpec.focus_trace_key_default || "");
      let focusSelectionKey = "";
      let cameraViewMode = "free";
      let earthViewFocusDistance = null;
      let currentLassoSelectionMask = null;
      const pressedKeys = new Set();
      let lastAnimationTimestamp = null;
      let clusterFilterParameterKey = String(clusterFilterSpec.default_parameter_key || "");
      const clusterFilterRangeStateByKey = {};
      let activeVolumeKey = volumeLayers.length ? String(volumeLayers[0].key) : null;
      const volumeStateByKey = {};
      const traceStyleStateByKey = {};
      const legendPanelStateByKey = {};

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
      const volumeSupported = volumeLayers.length > 0 && Boolean(renderer.capabilities && renderer.capabilities.isWebGL2);

      const scene = new THREE.Scene();
      scene.background = new THREE.Color(theme.scene_bgcolor || theme.paper_bgcolor || "#000000");

      const sceneUp = sceneSpec.camera_up || { x: 0.0, y: 0.0, z: 1.0 };
      const sceneUpVector = new THREE.Vector3(sceneUp.x ?? 0.0, sceneUp.y ?? 0.0, sceneUp.z ?? 1.0).normalize();
      const camera = new THREE.PerspectiveCamera(38, 1, 0.1, Math.max((sceneSpec.max_span || 1) * 20.0, 10000.0));
      camera.up.set(sceneUp.x ?? 0.0, sceneUp.y ?? 0.0, sceneUp.z ?? 1.0);
      const controls = new OrbitControls(camera, renderer.domElement);
      controls.enableDamping = true;
      controls.dampingFactor = 0.08;
      controls.rotateSpeed = 0.7;
      controls.panSpeed = 0.7;
      controls.zoomSpeed = globalScrollSpeed;
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
        const defaults = layer.default_controls || {};
        volumeStateByKey[String(layer.key)] = {
          visible: layer.visible !== false,
          vmin: Number(defaults.vmin),
          vmax: Number(defaults.vmax),
          opacity: Number(defaults.opacity),
          steps: Number(defaults.steps),
          alphaCoef: Number(defaults.alpha_coef),
          gradientStep: Number(defaults.gradient_step),
          stretch: normalizeVolumeStretch(defaults.stretch),
          colormap: String(defaults.colormap || (((layer.colormap_options || [])[0] || {}).name || "inferno")),
        };
      });

      legendItems.forEach((item) => {
        const itemKey = String(item.key);
        legendPanelStateByKey[itemKey] = false;
        if (volumeLayersByKey.has(itemKey)) {
          return;
        }
        traceStyleStateByKey[itemKey] = {
          color: String(item.default_color || item.color || theme.axis_color || "#ffffff"),
          opacity: Math.min(Math.max(Number(item.default_opacity ?? 1.0), 0.0), 1.0),
          sizeScale: 1.0,
          hasPoints: Boolean(item.has_points),
          hasSegments: Boolean(item.has_segments),
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

      function widgetPanelForKey(widgetKey) {
        if (widgetKey === "sky") {
          return skyPanelEl;
        }
        if (widgetKey === "age_kde") {
          return ageKdePanelEl;
        }
        if (widgetKey === "cluster_filter") {
          return clusterFilterPanelEl;
        }
        return null;
      }

      function widgetEnabled(widgetKey) {
        if (widgetKey === "sky") {
          return Boolean(skySpec.enabled);
        }
        if (widgetKey === "age_kde") {
          return Boolean(ageKdeSpec.enabled);
        }
        if (widgetKey === "cluster_filter") {
          return Boolean(clusterFilterSpec.enabled);
        }
        return false;
      }

      function widgetModeForKey(widgetKey) {
        if (widgetKey === "sky") {
          return skyPanelMode;
        }
        if (widgetKey === "age_kde") {
          return ageKdePanelMode;
        }
        if (widgetKey === "cluster_filter") {
          return clusterFilterPanelMode;
        }
        return "hidden";
      }

      function setWidgetModeValue(widgetKey, mode) {
        if (widgetKey === "sky") {
          skyPanelMode = mode;
        } else if (widgetKey === "age_kde") {
          ageKdePanelMode = mode;
        } else if (widgetKey === "cluster_filter") {
          clusterFilterPanelMode = mode;
        }
      }

      function widgetDefaultRect(widgetKey) {
        if (widgetKey === "sky") {
          return {
            left: Math.max(12, window.innerWidth - Math.min(window.innerWidth * 0.38, 560) - 12),
            top: 72,
            width: Math.min(window.innerWidth * 0.38, 560),
            height: Math.min(window.innerHeight * 0.56, 560),
          };
        }
        if (widgetKey === "age_kde") {
          return {
            left: Math.max(12, window.innerWidth - Math.min(window.innerWidth * 0.36, 540) - 44),
            top: 96,
            width: Math.min(window.innerWidth * 0.36, 540),
            height: Math.min(window.innerHeight * 0.42, 360),
          };
        }
        if (widgetKey === "cluster_filter") {
          return {
            left: Math.max(12, window.innerWidth - Math.min(window.innerWidth * 0.34, 480) - 72),
            top: 108,
            width: Math.min(window.innerWidth * 0.34, 480),
            height: Math.min(window.innerHeight * 0.40, 340),
          };
        }
        return { left: 12, top: 72, width: 360, height: 260 };
      }

      function raiseWidget(widgetKey) {
        const panelEl = widgetPanelForKey(widgetKey);
        if (!panelEl) {
          return;
        }
        widgetZIndexCounter += 1;
        panelEl.style.zIndex = String(widgetZIndexCounter);
      }

      function clampWidgetPosition(left, top, width, height) {
        const margin = 6;
        return {
          left: Math.min(Math.max(margin, left), Math.max(margin, window.innerWidth - width - margin)),
          top: Math.min(Math.max(margin, top), Math.max(margin, window.innerHeight - height - margin)),
        };
      }

      function resizeWidgetRect(state, clientX, clientY) {
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

      function storeWidgetRect(widgetKey) {
        const panelEl = widgetPanelForKey(widgetKey);
        if (!panelEl) {
          return;
        }
        const rect = panelEl.getBoundingClientRect();
        panelEl.dataset.normalLeft = String(rect.left);
        panelEl.dataset.normalTop = String(rect.top);
        panelEl.dataset.normalWidth = String(rect.width);
        panelEl.dataset.normalHeight = String(rect.height);
      }

      function restoreWidgetRect(widgetKey) {
        const panelEl = widgetPanelForKey(widgetKey);
        if (!panelEl) {
          return;
        }
        const left = Number(panelEl.dataset.normalLeft);
        const top = Number(panelEl.dataset.normalTop);
        const width = Number(panelEl.dataset.normalWidth);
        const height = Number(panelEl.dataset.normalHeight);
        const defaults = widgetDefaultRect(widgetKey);
        const next = [left, top, width, height].every(Number.isFinite)
          ? clampWidgetPosition(left, top, width, height)
          : clampWidgetPosition(defaults.left, defaults.top, defaults.width, defaults.height);
        const nextWidth = Number.isFinite(width) ? width : defaults.width;
        const nextHeight = Number.isFinite(height) ? height : defaults.height;
        panelEl.style.left = `${next.left}px`;
        panelEl.style.top = `${next.top}px`;
        panelEl.style.right = "auto";
        panelEl.style.bottom = "auto";
        panelEl.style.width = `${nextWidth}px`;
        panelEl.style.height = `${nextHeight}px`;
      }

      function restoreInitialLassoSelectionMask() {
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
          if (typeof savedGlobalControls.camera_view_mode === "string" && savedGlobalControls.camera_view_mode) {
            cameraViewMode = String(savedGlobalControls.camera_view_mode);
          }
          if (Number.isFinite(Number(savedGlobalControls.earth_view_focus_distance_pc))) {
            earthViewFocusDistance = Number(savedGlobalControls.earth_view_focus_distance_pc);
          }
        }

        if (initialState.active_volume_key && volumeLayersByKey.has(String(initialState.active_volume_key))) {
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

        const savedLegendPanels = initialState.legend_panel_state;
        if (savedLegendPanels && typeof savedLegendPanels === "object") {
          Object.entries(savedLegendPanels).forEach(([itemKey, isOpen]) => {
            if (Object.prototype.hasOwnProperty.call(legendPanelStateByKey, String(itemKey))) {
              legendPanelStateByKey[String(itemKey)] = Boolean(isOpen);
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
          ["sky", "age_kde", "cluster_filter"].forEach((widgetKey) => {
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

        applyGlobalControlState();
        clampClusterFilterRangeForParameter(activeClusterFilterParameterSpec());
        pruneSelectionsToActiveClusterFilter();
        applyCameraViewMode();
        applyThemePreset(activeThemeKey, { rerender: false, syncInput: false });
        renderSceneControls();
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

      function captureLassoSelectionMaskState(mask) {
        if (!mask || !mask.maskTexture || !mask.viewProjectionMatrix) {
          return null;
        }
        const image = mask.maskTexture.image;
        let dataUrl = null;
        if (image && typeof image.toDataURL === "function") {
          dataUrl = image.toDataURL("image/png");
        } else if (image && typeof image.src === "string" && image.src) {
          dataUrl = image.src;
        }
        if (!dataUrl) {
          return null;
        }
        return {
          data_url: dataUrl,
          view_projection_matrix: Array.from(mask.viewProjectionMatrix.elements || []).map((value) => Number(value)),
          polygon_ndc: Array.isArray(mask.polygonNdc)
            ? mask.polygonNdc
              .map((point) => ({
                x: Number(point && point.x),
                y: Number(point && point.y),
              }))
              .filter((point) => Number.isFinite(point.x) && Number.isFinite(point.y))
            : [],
        };
      }

      function captureRuntimeState() {
        return safeJsonClone({
          current_group: currentGroup,
          current_frame_index: currentFrameIndex,
          legend_state: legendState,
          trace_style_state: traceStyleStateByKey,
          legend_panel_state: legendPanelStateByKey,
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
            fade_in_time_myr: fadeInTimeMyr,
            fade_in_and_out_enabled: fadeInAndOutEnabled,
            focus_trace_key: focusTraceKey,
            focus_selection_key: focusSelectionKey,
            axes_visible: axesVisible,
            camera_view_mode: cameraViewMode,
            earth_view_focus_distance_pc: earthViewFocusDistance,
          },
          camera: {
            position: { x: camera.position.x, y: camera.position.y, z: camera.position.z },
            target: { x: controls.target.x, y: controls.target.y, z: controls.target.z },
            up: { x: camera.up.x, y: camera.up.y, z: camera.up.z },
          },
          drawers: {
            tools_open: toolsShellEl.dataset.open === "true",
            controls_open: controlsShellEl.dataset.open === "true",
          },
          zen_mode_enabled: zenModeEnabled,
          widgets: {
            sky: captureWidgetState("sky"),
            age_kde: captureWidgetState("age_kde"),
            cluster_filter: captureWidgetState("cluster_filter"),
          },
          cluster_filter_state: {
            parameter_key: clusterFilterParameterKey,
            ranges_by_key: clusterFilterRangeStateByKey,
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

      function activeHoveredClusterKeys() {
        if (!crossHoverEnabled()) {
          return new Set();
        }
        return new Set(
          [localHoveredClusterKey, skyHoveredClusterKey]
            .map((value) => normalizeMemberKey(value))
            .filter(Boolean)
        );
      }

      function applySceneHoverState() {
        const activeKeys = activeHoveredClusterKeys();
        selectionSpriteEntriesByKey.forEach((entries, key) => {
          const isActive = activeKeys.has(key);
          entries.forEach((entry) => {
            const scale = entry.baseScale * (isActive ? 1.45 : 1.0);
            entry.sprite.scale.set(scale, scale, scale);
          });
        });
      }

      function postParentHoverToSkyFrame() {
        if (!skySpec.enabled || !skyFrameEl || !skyFrameEl.contentWindow) {
          return;
        }
        const clusterKey = crossHoverEnabled() ? normalizeMemberKey(localHoveredClusterKey) : "";
        if (clusterKey === lastSentSkyHoverClusterKey) {
          return;
        }
        lastSentSkyHoverClusterKey = clusterKey;
        skyFrameEl.contentWindow.postMessage({
          type: "oviz-parent-hover-cluster",
          clusterKey: clusterKey || null,
        }, "*");
      }

      function setLocalHoveredClusterKey(clusterKey) {
        const nextKey = crossHoverEnabled() ? normalizeMemberKey(clusterKey) : "";
        if (nextKey === localHoveredClusterKey) {
          return;
        }
        localHoveredClusterKey = nextKey;
        applySceneHoverState();
        postParentHoverToSkyFrame();
      }

      function setSkyHoveredClusterKey(clusterKey) {
        const nextKey = crossHoverEnabled() ? normalizeMemberKey(clusterKey) : "";
        if (nextKey === skyHoveredClusterKey) {
          return;
        }
        skyHoveredClusterKey = nextKey;
        applySceneHoverState();
      }

      function clearCrossHoverState() {
        localHoveredClusterKey = "";
        skyHoveredClusterKey = "";
        lastSentSkyHoverClusterKey = null;
        applySceneHoverState();
        postParentHoverToSkyFrame();
      }

      function hasActiveLassoSelectionMask() {
        return Boolean(
          currentLassoSelectionMask
          && currentLassoSelectionMask.maskTexture
          && currentLassoSelectionMask.viewProjectionMatrix
        );
      }

      function activeVolumeLassoSelectionMask() {
        return lassoVolumeSelectionEnabled ? currentLassoSelectionMask : null;
      }

      function disposeLassoSelectionMask(mask) {
        if (!mask || !mask.maskTexture) {
          return;
        }
        mask.maskTexture.dispose();
      }

      function canvasPointToNdc(point) {
        const width = Math.max(canvas.clientWidth || 0, 1);
        const height = Math.max(canvas.clientHeight || 0, 1);
        return {
          x: (Number(point.x) / width) * 2.0 - 1.0,
          y: 1.0 - (Number(point.y) / height) * 2.0,
        };
      }

      function captureLassoSelectionMask(points) {
        const source = Array.isArray(points) ? points : [];
        if (source.length < 3) {
          return null;
        }
        const pointsNdc = source
          .map((point) => canvasPointToNdc(point))
          .filter((point) => Number.isFinite(point.x) && Number.isFinite(point.y));
        if (pointsNdc.length < 3) {
          return null;
        }

        const maskCanvasSize = 512;
        const maskCanvas = document.createElement("canvas");
        maskCanvas.width = maskCanvasSize;
        maskCanvas.height = maskCanvasSize;
        const maskCtx = maskCanvas.getContext("2d");
        maskCtx.clearRect(0, 0, maskCanvasSize, maskCanvasSize);
        maskCtx.fillStyle = "#ffffff";
        maskCtx.beginPath();
        pointsNdc.forEach((point, index) => {
          const x = (point.x * 0.5 + 0.5) * maskCanvasSize;
          const y = (0.5 - point.y * 0.5) * maskCanvasSize;
          if (index === 0) {
            maskCtx.moveTo(x, y);
          } else {
            maskCtx.lineTo(x, y);
          }
        });
        maskCtx.closePath();
        maskCtx.fill();

        const maskTexture = new THREE.CanvasTexture(maskCanvas);
        maskTexture.minFilter = THREE.NearestFilter;
        maskTexture.magFilter = THREE.NearestFilter;
        maskTexture.wrapS = THREE.ClampToEdgeWrapping;
        maskTexture.wrapT = THREE.ClampToEdgeWrapping;
        maskTexture.flipY = false;
        maskTexture.generateMipmaps = false;
        maskTexture.needsUpdate = true;
        const maskImageData = maskCtx.getImageData(0, 0, maskCanvasSize, maskCanvasSize);

        camera.updateMatrixWorld(true);
        camera.updateProjectionMatrix();
        const viewProjectionMatrix = new THREE.Matrix4()
          .multiplyMatrices(camera.projectionMatrix, camera.matrixWorldInverse);
        return {
          maskTexture,
          viewProjectionMatrix,
          polygonNdc: pointsNdc.map((point) => ({ x: Number(point.x), y: Number(point.y) })),
          maskSize: maskCanvasSize,
          maskAlphaData: maskImageData ? maskImageData.data : null,
        };
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
        const dustHint = lassoVolumeSelectionEnabled
          ? (hasActiveLassoSelectionMask() ? "Dust lasso: on" : "Dust lasso: armed")
          : "Dust lasso: off";
        if (!activeSelections.length && !focusLabel) {
          return `Shift+drag or use Lasso to select clusters on the current frame.\n${dustHint}\n${clickHint}`;
        }
        if (!activeSelections.length && focusLabel) {
          return `Focused: ${focusLabel}\n${dustHint}\n${clickHint}`;
        }
        const labels = activeSelections
          .map((selection) => selectionKeyFor(selection))
          .filter(Boolean);
        const preview = labels.slice(0, 4).join(", ");
        const suffix = labels.length > 4 ? ` +${labels.length - 4} more` : "";
        const focusText = focusLabel ? `\nFocused: ${focusLabel}` : "";
        return `${labels.length} selected\n${preview}${suffix}${focusText}\n${dustHint}\n${clickHint}`;
      }

      function skyReadoutText(selections, catalogPayload, mode = "overview", volumeOverlay = null) {
        const activeSelections = uniqueSelections(selections);
        const clusterCatalogPayload = catalogPayload || [];
        const dustPixelCount = volumeOverlay && Number.isFinite(Number(volumeOverlay.non_zero_pixels))
          ? Number(volumeOverlay.non_zero_pixels)
          : 0;
        const dustSampleCount = volumeOverlay && Number.isFinite(Number(volumeOverlay.sample_count))
          ? Number(volumeOverlay.sample_count)
          : 0;
        const overviewText = (
          "View: Mollweide all-sky\\n"
          + "Center: Galactic Center\\n"
          + `Survey: ${String(skySpec.survey || "P/DSS2/color")}`
        );
        if (mode !== "click") {
          if (!activeSelections.length && !dustPixelCount) {
            return overviewText;
          }
          const labels = activeSelections
            .map((selection) => selectionKeyFor(selection))
            .filter(Boolean);
          const starCount = clusterCatalogPayload.reduce(
            (total, catalog) => total + (((catalog && catalog.points) || []).length || 0),
            0
          );
          const dustLine = dustPixelCount
            ? `Dust sky pixels: ${dustPixelCount}${dustSampleCount ? ` from ${dustSampleCount} samples` : ""}\n`
            : "";
          if (!activeSelections.length) {
            return dustLine + overviewText;
          }
          return (
            `Clusters selected: ${labels.length}\n`
            + `${labels.slice(0, 6).join(", ")}${labels.length > 6 ? ` +${labels.length - 6} more` : ""}\n`
            + `Stars shown: ${starCount}\n`
            + dustLine
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

      function normalizeMemberKey(value) {
        return String(value || "")
          .trim()
          .toLowerCase()
          .replace(/[_\s]+/g, " ");
      }

      const GALACTIC_TO_ICRS_MATRIX = [
        -0.0548755604162154, 0.4941094278755837, -0.8676661490190047,
        -0.8734370902348850, -0.4448296299600112, -0.1980763734312015,
        -0.4838350155487132, 0.7469822444972189, 0.4559837761750669,
      ];

      function normalizeSkyLongitude(value) {
        let lon = Number(value);
        if (!Number.isFinite(lon)) {
          return NaN;
        }
        lon %= 360.0;
        if (lon < 0.0) {
          lon += 360.0;
        }
        return lon;
      }

      function wrapLongitudeDeltaDeg(value) {
        let delta = Number(value);
        if (!Number.isFinite(delta)) {
          return NaN;
        }
        while (delta <= -180.0) {
          delta += 360.0;
        }
        while (delta > 180.0) {
          delta -= 360.0;
        }
        return delta;
      }

      function galacticLonLatDegFromCartesian(x, y, z) {
        const xx = Number(x);
        const yy = Number(y);
        const zz = Number(z);
        if (!Number.isFinite(xx) || !Number.isFinite(yy) || !Number.isFinite(zz)) {
          return null;
        }
        const distance = Math.sqrt(xx * xx + yy * yy + zz * zz);
        if (!(distance > 1e-9)) {
          return null;
        }
        return {
          l: normalizeSkyLongitude(Math.atan2(yy, xx) * 180.0 / Math.PI),
          b: Math.asin(Math.min(1.0, Math.max(-1.0, zz / distance))) * 180.0 / Math.PI,
          distance,
        };
      }

      function icrsDegFromGalacticDeg(lDeg, bDeg) {
        const lon = Number(lDeg) * Math.PI / 180.0;
        const lat = Number(bDeg) * Math.PI / 180.0;
        if (!Number.isFinite(lon) || !Number.isFinite(lat)) {
          return null;
        }
        const cosLat = Math.cos(lat);
        const xGal = cosLat * Math.cos(lon);
        const yGal = cosLat * Math.sin(lon);
        const zGal = Math.sin(lat);
        const xEq = (
          GALACTIC_TO_ICRS_MATRIX[0] * xGal
          + GALACTIC_TO_ICRS_MATRIX[1] * yGal
          + GALACTIC_TO_ICRS_MATRIX[2] * zGal
        );
        const yEq = (
          GALACTIC_TO_ICRS_MATRIX[3] * xGal
          + GALACTIC_TO_ICRS_MATRIX[4] * yGal
          + GALACTIC_TO_ICRS_MATRIX[5] * zGal
        );
        const zEq = (
          GALACTIC_TO_ICRS_MATRIX[6] * xGal
          + GALACTIC_TO_ICRS_MATRIX[7] * yGal
          + GALACTIC_TO_ICRS_MATRIX[8] * zGal
        );
        const dec = Math.asin(Math.min(1.0, Math.max(-1.0, zEq))) * 180.0 / Math.PI;
        const ra = normalizeSkyLongitude(Math.atan2(yEq, xEq) * 180.0 / Math.PI);
        if (!Number.isFinite(ra) || !Number.isFinite(dec)) {
          return null;
        }
        return { ra, dec };
      }

      function volumeOverlayStretchValue(value, stretchName) {
        const clamped = Math.min(Math.max(Number(value), 0.0), 1.0);
        const stretch = normalizeVolumeStretch(stretchName);
        if (stretch === "log10") {
          const strength = 999.0;
          return Math.log(1.0 + strength * clamped) / Math.log(1.0 + strength);
        }
        if (stretch === "asinh") {
          const strength = 10.0;
          const numer = Math.log(strength * clamped + Math.sqrt((strength * clamped) * (strength * clamped) + 1.0));
          const denom = Math.log(strength + Math.sqrt(strength * strength + 1.0));
          return denom > 0.0 ? numer / denom : clamped;
        }
        return clamped;
      }

      function galacticDirectionVectorFromLonLatDeg(lDeg, bDeg) {
        const lonRad = Number(lDeg) * Math.PI / 180.0;
        const latRad = Number(bDeg) * Math.PI / 180.0;
        if (!Number.isFinite(lonRad) || !Number.isFinite(latRad)) {
          return null;
        }
        const cosLat = Math.cos(latRad);
        return {
          x: cosLat * Math.cos(lonRad),
          y: cosLat * Math.sin(lonRad),
          z: Math.sin(latRad),
        };
      }

      function intersectRayWithBounds(direction, xBounds, yBounds, zBounds) {
        if (!direction) {
          return null;
        }
        const dir = [
          Number(direction.x) || 0.0,
          Number(direction.y) || 0.0,
          Number(direction.z) || 0.0,
        ];
        const bounds = [
          Array.isArray(xBounds) ? xBounds : [-0.5, 0.5],
          Array.isArray(yBounds) ? yBounds : [-0.5, 0.5],
          Array.isArray(zBounds) ? zBounds : [-0.5, 0.5],
        ];
        let tMin = -Infinity;
        let tMax = Infinity;
        for (let axis = 0; axis < 3; axis += 1) {
          const d = dir[axis];
          const low = Number(bounds[axis][0]);
          const high = Number(bounds[axis][1]);
          if (!Number.isFinite(low) || !Number.isFinite(high)) {
            return null;
          }
          if (Math.abs(d) < 1e-8) {
            if (0.0 < Math.min(low, high) || 0.0 > Math.max(low, high)) {
              return null;
            }
            continue;
          }
          let t0 = low / d;
          let t1 = high / d;
          if (t0 > t1) {
            const swap = t0;
            t0 = t1;
            t1 = swap;
          }
          tMin = Math.max(tMin, t0);
          tMax = Math.min(tMax, t1);
          if (!(tMax > tMin)) {
            return null;
          }
        }
        const near = Math.max(tMin, 0.0);
        if (!(tMax > near)) {
          return null;
        }
        return { tMin: near, tMax };
      }

      function sampleVolumeScalarTrilinear(scalarData, nx, ny, nz, xNorm, yNorm, zNorm) {
        if (!scalarData || !scalarData.length) {
          return 0.0;
        }
        const x = Math.min(Math.max(Number(xNorm), 0.0), 1.0) * Math.max(nx - 1, 0);
        const y = Math.min(Math.max(Number(yNorm), 0.0), 1.0) * Math.max(ny - 1, 0);
        const z = Math.min(Math.max(Number(zNorm), 0.0), 1.0) * Math.max(nz - 1, 0);
        const x0 = Math.floor(x);
        const y0 = Math.floor(y);
        const z0 = Math.floor(z);
        const x1 = Math.min(x0 + 1, nx - 1);
        const y1 = Math.min(y0 + 1, ny - 1);
        const z1 = Math.min(z0 + 1, nz - 1);
        const tx = x - x0;
        const ty = y - y0;
        const tz = z - z0;
        const strideY = nx;
        const strideZ = nx * ny;
        const index000 = z0 * strideZ + y0 * strideY + x0;
        const index100 = z0 * strideZ + y0 * strideY + x1;
        const index010 = z0 * strideZ + y1 * strideY + x0;
        const index110 = z0 * strideZ + y1 * strideY + x1;
        const index001 = z1 * strideZ + y0 * strideY + x0;
        const index101 = z1 * strideZ + y0 * strideY + x1;
        const index011 = z1 * strideZ + y1 * strideY + x0;
        const index111 = z1 * strideZ + y1 * strideY + x1;
        const c000 = Number(scalarData[index000] || 0) / 255.0;
        const c100 = Number(scalarData[index100] || 0) / 255.0;
        const c010 = Number(scalarData[index010] || 0) / 255.0;
        const c110 = Number(scalarData[index110] || 0) / 255.0;
        const c001 = Number(scalarData[index001] || 0) / 255.0;
        const c101 = Number(scalarData[index101] || 0) / 255.0;
        const c011 = Number(scalarData[index011] || 0) / 255.0;
        const c111 = Number(scalarData[index111] || 0) / 255.0;
        const c00 = c000 * (1.0 - tx) + c100 * tx;
        const c10 = c010 * (1.0 - tx) + c110 * tx;
        const c01 = c001 * (1.0 - tx) + c101 * tx;
        const c11 = c011 * (1.0 - tx) + c111 * tx;
        const c0 = c00 * (1.0 - ty) + c10 * ty;
        const c1 = c01 * (1.0 - ty) + c11 * ty;
        return c0 * (1.0 - tz) + c1 * tz;
      }

      function pointInsideProjectedLassoMask(worldX, worldY, worldZ, mask) {
        if (
          !mask
          || !mask.viewProjectionMatrix
          || !Array.isArray(mask.polygonNdc)
          || mask.polygonNdc.length < 3
        ) {
          return false;
        }
        const e = mask.viewProjectionMatrix.elements || [];
        if (e.length !== 16) {
          return false;
        }
        const clipX = e[0] * worldX + e[4] * worldY + e[8] * worldZ + e[12];
        const clipY = e[1] * worldX + e[5] * worldY + e[9] * worldZ + e[13];
        const clipZ = e[2] * worldX + e[6] * worldY + e[10] * worldZ + e[14];
        const clipW = e[3] * worldX + e[7] * worldY + e[11] * worldZ + e[15];
        if (!(clipW > 0.0)) {
          return false;
        }
        const ndcPoint = {
          x: clipX / clipW,
          y: clipY / clipW,
          z: clipZ / clipW,
        };
        if (!Number.isFinite(ndcPoint.x) || !Number.isFinite(ndcPoint.y) || !Number.isFinite(ndcPoint.z)) {
          return false;
        }
        if (ndcPoint.z < -1.0 || ndcPoint.z > 1.0) {
          return false;
        }
        if (ndcPoint.x < -1.0 || ndcPoint.x > 1.0 || ndcPoint.y < -1.0 || ndcPoint.y > 1.0) {
          return false;
        }
        const maskSize = Math.max(0, Math.round(Number(mask.maskSize) || 0));
        const maskAlphaData = mask.maskAlphaData;
        if (maskSize > 0 && maskAlphaData && maskAlphaData.length >= maskSize * maskSize * 4) {
          const maskX = Math.max(0, Math.min(maskSize - 1, Math.round((ndcPoint.x * 0.5 + 0.5) * (maskSize - 1))));
          const maskY = Math.max(0, Math.min(maskSize - 1, Math.round((0.5 - ndcPoint.y * 0.5) * (maskSize - 1))));
          const alphaIndex = ((maskY * maskSize) + maskX) * 4 + 3;
          return Number(maskAlphaData[alphaIndex] || 0) > 0;
        }
        return pointInPolygon(ndcPoint, mask.polygonNdc);
      }

      function buildVolumeSkyImageOverlaySpec(mode = "overview") {
        if (mode === "click") {
          return null;
        }
        const activeMask = activeVolumeLassoSelectionMask();
        if (
          !activeMask
          || !activeMask.viewProjectionMatrix
          || !Array.isArray(activeMask.polygonNdc)
          || activeMask.polygonNdc.length < 3
        ) {
          return null;
        }

        const frame = currentFrame();
        const frameTime = frame ? Number(frame.time) : NaN;
        const displayOffsetX = Number(plotGroup.position.x) || 0.0;
        const displayOffsetY = Number(plotGroup.position.y) || 0.0;
        const displayOffsetZ = Number(plotGroup.position.z) || 0.0;
        const collectedSamples = [];
        const renderLayers = [];
        let usedLayerCount = 0;

        volumeLayers.forEach((layer) => {
          const layerKey = String(layer.key);
          const state = volumeStateByKey[layerKey];
          if (!state || state.visible === false || legendState[layerKey] === false) {
            return;
          }
          if (layer.only_at_t0 && !approximatelyZero(frameTime)) {
            return;
          }

          const shape = layer.sky_overlay_shape || layer.shape || {};
          const nx = Math.max(1, Math.round(Number(shape.x) || 0));
          const ny = Math.max(1, Math.round(Number(shape.y) || 0));
          const nz = Math.max(1, Math.round(Number(shape.z) || 0));
          const totalVoxels = nx * ny * nz;
          if (!(totalVoxels > 0)) {
            return;
          }

          const scalarData = volumeSkyScalarArrayFor(layer);
          if (!scalarData || !scalarData.length) {
            return;
          }

          const option = volumeColormapOptionFor(layer, state.colormap);
          const colorBytes = option ? volumeColorBytesForOption(option) : null;
          if (!colorBytes || colorBytes.length < 4) {
            return;
          }

          const bounds = layer.bounds || {};
          const xBounds = Array.isArray(bounds.x) ? bounds.x : [-0.5, 0.5];
          const yBounds = Array.isArray(bounds.y) ? bounds.y : [-0.5, 0.5];
          const zBounds = Array.isArray(bounds.z) ? bounds.z : [-0.5, 0.5];
          const xSpan = Number(xBounds[1]) - Number(xBounds[0]);
          const ySpan = Number(yBounds[1]) - Number(yBounds[0]);
          const zSpan = Number(zBounds[1]) - Number(zBounds[0]);
          if (!(Number.isFinite(xSpan) && Number.isFinite(ySpan) && Number.isFinite(zSpan))) {
            return;
          }

          const windowState = normalizedVolumeWindowFor(layer, state);
          const low = Number(windowState.low);
          const high = Number(windowState.high);
          const span = Math.max(high - low, 1e-6);
          const desiredScanCount = 2500000;
          const stride = Math.max(1, Math.ceil(Math.cbrt(totalVoxels / desiredScanCount)));
          const minScaledThreshold = 0.02;
          const layerOpacityWeight = Math.max(0.05, Number(state.opacity) || 0.0);
          const lutSamples = Math.max(1, Math.floor(colorBytes.length / 4));
          const cellSizeX = Math.abs(xSpan / Math.max(nx, 1));
          const cellSizeY = Math.abs(ySpan / Math.max(ny, 1));
          const cellSizeZ = Math.abs(zSpan / Math.max(nz, 1));
          const cellSizeMin = Math.max(1e-6, Math.min(cellSizeX, cellSizeY, cellSizeZ));
          let layerUsed = false;

          for (let iz = 0; iz < nz; iz += stride) {
            const localZ = Number(zBounds[0]) + ((iz + 0.5) / nz) * zSpan;
            const displayedZ = localZ + displayOffsetZ;
            for (let iy = 0; iy < ny; iy += stride) {
              const localY = Number(yBounds[0]) + ((iy + 0.5) / ny) * ySpan;
              const displayedY = localY + displayOffsetY;
              for (let ix = 0; ix < nx; ix += stride) {
                const dataIndex = iz * nx * ny + iy * nx + ix;
                const normalizedValue = Number(scalarData[dataIndex] || 0) / 255.0;
                if (!(normalizedValue > low)) {
                  continue;
                }
                const scaledValue = Math.min(Math.max((normalizedValue - low) / span, 0.0), 1.0);
                const stretchedValue = volumeOverlayStretchValue(scaledValue, state.stretch);
                if (!(stretchedValue > minScaledThreshold)) {
                  continue;
                }

                const localX = Number(xBounds[0]) + ((ix + 0.5) / nx) * xSpan;
                const displayedX = localX + displayOffsetX;
                if (!pointInsideProjectedLassoMask(displayedX, displayedY, displayedZ, activeMask)) {
                  continue;
                }

                const galactic = galacticLonLatDegFromCartesian(localX, localY, localZ);
                if (!galactic) {
                  continue;
                }

                const lDeg = normalizeSkyLongitude(galactic.l);
                const bDeg = Number(galactic.b);
                if (!Number.isFinite(lDeg) || !Number.isFinite(bDeg)) {
                  continue;
                }

                const weight = stretchedValue * layerOpacityWeight;
                collectedSamples.push({
                  l: lDeg,
                  b: bDeg,
                  weight,
                });
                layerUsed = true;
              }
            }
          }

          if (layerUsed) {
            usedLayerCount += 1;
            renderLayers.push({
              layer,
              state,
              scalarData,
              colorBytes,
              lutSamples,
              low,
              high,
              span,
              opacityWeight: layerOpacityWeight,
              xBounds,
              yBounds,
              zBounds,
              xSpan,
              ySpan,
              zSpan,
              cellSizeMin,
              nx,
              ny,
              nz,
              baseRaySamples: Math.max(96, Math.min(320, Math.round((Number(state.samples) || 96) * 4.0))),
            });
          }
        });

        if (!collectedSamples.length || !renderLayers.length) {
          return null;
        }

        let sumSin = 0.0;
        let sumCos = 0.0;
        let sumLat = 0.0;
        let sumWeight = 0.0;
        collectedSamples.forEach((sample) => {
          const weight = Math.max(Number(sample.weight) || 0.0, 1e-6);
          const lonRad = Number(sample.l) * Math.PI / 180.0;
          sumSin += Math.sin(lonRad) * weight;
          sumCos += Math.cos(lonRad) * weight;
          sumLat += Number(sample.b) * weight;
          sumWeight += weight;
        });

        const baseCenterLon = normalizeSkyLongitude(Math.atan2(sumSin, sumCos) * 180.0 / Math.PI);
        const baseCenterLat = sumWeight > 0.0 ? (sumLat / sumWeight) : 0.0;
        let minRelLon = Infinity;
        let maxRelLon = -Infinity;
        let minLat = Infinity;
        let maxLat = -Infinity;
        collectedSamples.forEach((sample) => {
          const relLon = wrapLongitudeDeltaDeg(Number(sample.l) - baseCenterLon);
          sample.relLonBase = relLon;
          minRelLon = Math.min(minRelLon, relLon);
          maxRelLon = Math.max(maxRelLon, relLon);
          minLat = Math.min(minLat, Number(sample.b));
          maxLat = Math.max(maxLat, Number(sample.b));
        });

        const lonSpanRaw = Math.max(maxRelLon - minRelLon, 0.25);
        const latSpanRaw = Math.max(maxLat - minLat, 0.25);
        const lonSpan = Math.min(Math.max(lonSpanRaw * 1.12, 2.0), 180.0);
        const latSpan = Math.min(Math.max(latSpanRaw * 1.12, 2.0), 180.0);
        const lonMidOffset = 0.5 * (minRelLon + maxRelLon);
        const overlayCenterLon = normalizeSkyLongitude(baseCenterLon + lonMidOffset);
        const overlayCenterLat = Math.min(Math.max(baseCenterLat + 0.5 * ((minLat + maxLat) - (2.0 * baseCenterLat)), -89.0), 89.0);
        const aspect = Math.max(latSpan / lonSpan, 0.35);
        const width = lonSpan <= 10.0 ? 2048 : (lonSpan <= 28.0 ? 1536 : (lonSpan <= 70.0 ? 1024 : 768));
        const height = Math.max(384, Math.min(1536, Math.round(width * aspect)));
        const pixelCount = width * height;
        const weightGrid = new Float32Array(pixelCount);
        const redGrid = new Float32Array(pixelCount);
        const greenGrid = new Float32Array(pixelCount);
        const blueGrid = new Float32Array(pixelCount);
        const halfLonSpan = lonSpan * 0.5;
        const halfLatSpan = latSpan * 0.5;
        for (let yIndex = 0; yIndex < height; yIndex += 1) {
          const lat = overlayCenterLat + halfLatSpan - (((yIndex + 0.5) / height) * latSpan);
          const latRad = lat * Math.PI / 180.0;
          const sinLat = Math.sin(latRad);
          const cosLat = Math.cos(latRad);
          for (let xIndex = 0; xIndex < width; xIndex += 1) {
            const relLon = halfLonSpan - (((xIndex + 0.5) / width) * lonSpan);
            const lon = normalizeSkyLongitude(overlayCenterLon + relLon);
            const lonRad = lon * Math.PI / 180.0;
            const rayDirection = {
              x: cosLat * Math.cos(lonRad),
              y: cosLat * Math.sin(lonRad),
              z: sinLat,
            };
            const pixelIndex = yIndex * width + xIndex;

            for (let layerIndex = 0; layerIndex < renderLayers.length; layerIndex += 1) {
              const layerInfo = renderLayers[layerIndex];
              const hit = intersectRayWithBounds(
                rayDirection,
                layerInfo.xBounds,
                layerInfo.yBounds,
                layerInfo.zBounds,
              );
              if (!hit) {
                continue;
              }
              const rayLength = hit.tMax - hit.tMin;
              if (!(rayLength > 0.0)) {
                continue;
              }
              const sampleCount = Math.max(
                48,
                Math.min(
                  layerInfo.baseRaySamples,
                  Math.round(rayLength / Math.max(layerInfo.cellSizeMin * 0.8, 1e-6))
                )
              );
              if (!(sampleCount > 0)) {
                continue;
              }
              const dt = rayLength / sampleCount;
              const contributionScale = dt / Math.max(layerInfo.cellSizeMin, 1e-6);

              for (let sampleIndex = 0; sampleIndex < sampleCount; sampleIndex += 1) {
                const t = hit.tMin + ((sampleIndex + 0.5) / sampleCount) * rayLength;
                const localX = rayDirection.x * t;
                const localY = rayDirection.y * t;
                const localZ = rayDirection.z * t;
                if (
                  !pointInsideProjectedLassoMask(
                    localX + displayOffsetX,
                    localY + displayOffsetY,
                    localZ + displayOffsetZ,
                    activeMask,
                  )
                ) {
                  continue;
                }
                const xNorm = (localX - Number(layerInfo.xBounds[0])) / Math.max(layerInfo.xSpan, 1e-6);
                const yNorm = (localY - Number(layerInfo.yBounds[0])) / Math.max(layerInfo.ySpan, 1e-6);
                const zNorm = (localZ - Number(layerInfo.zBounds[0])) / Math.max(layerInfo.zSpan, 1e-6);
                if (
                  xNorm < 0.0 || xNorm > 1.0
                  || yNorm < 0.0 || yNorm > 1.0
                  || zNorm < 0.0 || zNorm > 1.0
                ) {
                  continue;
                }
                const normalizedValue = sampleVolumeScalarTrilinear(
                  layerInfo.scalarData,
                  layerInfo.nx,
                  layerInfo.ny,
                  layerInfo.nz,
                  xNorm,
                  yNorm,
                  zNorm,
                );
                if (!(normalizedValue > layerInfo.low)) {
                  continue;
                }
                const scaledValue = Math.min(
                  Math.max((normalizedValue - layerInfo.low) / layerInfo.span, 0.0),
                  1.0
                );
                const stretchedValue = volumeOverlayStretchValue(scaledValue, layerInfo.state.stretch);
                if (!(stretchedValue > 0.005)) {
                  continue;
                }
                const contribution = stretchedValue * layerInfo.opacityWeight * contributionScale;
                if (!(contribution > 0.0)) {
                  continue;
                }
                const colorIndex = Math.max(
                  0,
                  Math.min(
                    layerInfo.lutSamples - 1,
                    Math.round(stretchedValue * (layerInfo.lutSamples - 1))
                  )
                ) * 4;
                weightGrid[pixelIndex] += contribution;
                redGrid[pixelIndex] += (Number(layerInfo.colorBytes[colorIndex]) || 0) * contribution;
                greenGrid[pixelIndex] += (Number(layerInfo.colorBytes[colorIndex + 1]) || 0) * contribution;
                blueGrid[pixelIndex] += (Number(layerInfo.colorBytes[colorIndex + 2]) || 0) * contribution;
              }
            }
          }
        }

        let maxWeight = 0.0;
        let nonZeroPixels = 0;
        for (let index = 0; index < pixelCount; index += 1) {
          const weight = Number(weightGrid[index]);
          if (weight > 0.0) {
            nonZeroPixels += 1;
            if (weight > maxWeight) {
              maxWeight = weight;
            }
          }
        }

        if (!(maxWeight > 0.0) || !nonZeroPixels) {
          return null;
        }

        const overlayCanvas = document.createElement("canvas");
        overlayCanvas.width = width;
        overlayCanvas.height = height;
        const overlayCtx = overlayCanvas.getContext("2d");
        if (!overlayCtx) {
          return null;
        }
        const imageData = overlayCtx.createImageData(width, height);
        const rgba = imageData.data;
        const normalizer = maxWeight > 0.0 ? maxWeight : 1.0;

        for (let index = 0; index < pixelCount; index += 1) {
          const weight = Number(weightGrid[index]);
          if (!(weight > 0.0)) {
            continue;
          }
          const normalized = Math.min(weight / normalizer, 1.0);
          const intensity = Math.pow(normalized, 0.72);
          const outIndex = index * 4;
          const invWeight = 1.0 / weight;
          rgba[outIndex] = Math.max(
            0,
            Math.min(255, Math.round((Number(redGrid[index]) * invWeight) * (0.30 + 0.70 * intensity)))
          );
          rgba[outIndex + 1] = Math.max(
            0,
            Math.min(255, Math.round((Number(greenGrid[index]) * invWeight) * (0.30 + 0.70 * intensity)))
          );
          rgba[outIndex + 2] = Math.max(
            0,
            Math.min(255, Math.round((Number(blueGrid[index]) * invWeight) * (0.30 + 0.70 * intensity)))
          );
          rgba[outIndex + 3] = Math.max(0, Math.min(255, Math.round(255.0 * Math.pow(normalized, 0.78))));
        }

        overlayCtx.putImageData(imageData, 0, 0);
        return {
          kind: "volume_image",
          name: usedLayerCount > 1 ? "Selected Dust Intensity" : "Selected Dust",
          data_url: overlayCanvas.toDataURL("image/png"),
          width,
          height,
          sample_count: collectedSamples.length,
          non_zero_pixels: nonZeroPixels,
          wcs: {
            NAXIS: 2,
            NAXIS1: width,
            NAXIS2: height,
            CTYPE1: "GLON-CAR",
            CTYPE2: "GLAT-CAR",
            CUNIT1: "deg",
            CUNIT2: "deg",
            CRPIX1: (width / 2.0) + 0.5,
            CRPIX2: (height / 2.0) + 0.5,
            CRVAL1: overlayCenterLon,
            CRVAL2: overlayCenterLat,
            CD1_1: -lonSpan / width,
            CD1_2: 0.0,
            CD2_1: 0.0,
            CD2_2: -latSpan / height,
            LONPOLE: 0.0,
          },
        };
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

      function buildAladinCatalogPayload(selections, mode = "overview") {
        const activeSelections = uniqueSelections(selections);
        const payload = activeSelections.map((selection) => {
          const resolvedMembers = resolveMemberPoints(selection);
          const lookupName = resolvedMembers.key || selectionKeyFor(selection) || "Selection";
          const traceName = selection && selection.trace_name ? String(selection.trace_name) : lookupName;
          const clusterKey = normalizedSelectionKeyFor(selection);
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
              clusterKey,
            }));
          } else if (Number.isFinite(Number(selection.ra_deg)) && Number.isFinite(Number(selection.dec_deg))) {
            points = [{
              l: Number(selection.l_deg),
              b: Number(selection.b_deg),
              ra: Number(selection.ra_deg),
              dec: Number(selection.dec_deg),
              label: lookupName,
              clusterKey,
            }];
          }
          return {
            name: lookupName,
            traceName,
            color: clusterColor,
            opacity: 1.0,
            sourceSize: points.length > 1 ? 4 : 7,
            points,
          };
        }).filter((catalog) => (catalog.points || []).length);

        if (mode === "click" || payload.length <= 1) {
          return payload;
        }

        const grouped = new Map();
        payload.forEach((catalog) => {
          const groupKey = String(catalog.traceName || catalog.name || "Selected trace");
          if (!grouped.has(groupKey)) {
            grouped.set(groupKey, {
              name: groupKey,
              color: catalog.color,
              opacity: 1.0,
              sourceSize: 4,
              points: [],
              seen: new Set(),
            });
          }
          const group = grouped.get(groupKey);
          (catalog.points || []).forEach((point) => {
            const ra = Number(point.ra);
            const dec = Number(point.dec);
            const label = String(point.label || catalog.name || "Selection");
            const clusterKey = String(point.clusterKey || "");
            const key = `${ra.toFixed(8)}|${dec.toFixed(8)}|${label}|${clusterKey}`;
            if (!Number.isFinite(ra) || !Number.isFinite(dec) || group.seen.has(key)) {
              return;
            }
            group.seen.add(key);
            group.points.push({
              l: Number(point.l),
              b: Number(point.b),
              ra,
              dec,
              label,
              clusterKey,
            });
          });
        });

        return Array.from(grouped.values())
          .map((group) => ({
            name: group.name,
            color: group.color,
            opacity: group.opacity,
            sourceSize: group.sourceSize,
            points: group.points,
          }))
          .filter((group) => group.points.length);
      }

      function angularSeparationDeg(ra1Deg, dec1Deg, ra2Deg, dec2Deg) {
        const rad = Math.PI / 180.0;
        const sin1 = Math.sin(dec1Deg * rad);
        const sin2 = Math.sin(dec2Deg * rad);
        const cos1 = Math.cos(dec1Deg * rad);
        const cos2 = Math.cos(dec2Deg * rad);
        const deltaRa = (ra1Deg - ra2Deg) * rad;
        const cosSep = Math.min(1.0, Math.max(-1.0, sin1 * sin2 + cos1 * cos2 * Math.cos(deltaRa)));
        return Math.acos(cosSep) * 180.0 / Math.PI;
      }

      function skyFocusFromPayload(selections, catalogPayload) {
        const focusPoints = [];
        (catalogPayload || []).forEach((catalog) => {
          (catalog.points || []).forEach((point) => {
            const ra = Number(point.ra);
            const dec = Number(point.dec);
            if (Number.isFinite(ra) && Number.isFinite(dec)) {
              focusPoints.push({ ra, dec });
            }
          });
        });
        if (!focusPoints.length) {
          uniqueSelections(selections).forEach((selection) => {
            const ra = Number(selection.ra_deg);
            const dec = Number(selection.dec_deg);
            if (Number.isFinite(ra) && Number.isFinite(dec)) {
              focusPoints.push({ ra, dec });
            }
          });
        }
        if (!focusPoints.length) {
          return null;
        }

        const rad = Math.PI / 180.0;
        let sx = 0.0;
        let sy = 0.0;
        let sz = 0.0;
        focusPoints.forEach((point) => {
          const ra = point.ra * rad;
          const dec = point.dec * rad;
          const cosDec = Math.cos(dec);
          sx += cosDec * Math.cos(ra);
          sy += cosDec * Math.sin(ra);
          sz += Math.sin(dec);
        });

        const norm = Math.sqrt(sx * sx + sy * sy + sz * sz);
        let centerRa = focusPoints[0].ra;
        let centerDec = focusPoints[0].dec;
        if (norm > 1e-9) {
          centerRa = Math.atan2(sy, sx) * 180.0 / Math.PI;
          if (centerRa < 0.0) {
            centerRa += 360.0;
          }
          centerDec = Math.asin(Math.min(1.0, Math.max(-1.0, sz / norm))) * 180.0 / Math.PI;
        }

        let maxSep = 0.0;
        focusPoints.forEach((point) => {
          maxSep = Math.max(maxSep, angularSeparationDeg(centerRa, centerDec, point.ra, point.dec));
        });

        const radiusDeg = Number(skySpec.radius_deg || 1.0);
        return {
          ra: centerRa,
          dec: centerDec,
          fovDeg: Math.min(Math.max(radiusDeg * 2.4, maxSep * 2.8, 1.2), 180.0),
        };
      }

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
              name: String(volumeOverlay.name || "Selected Dust"),
              imgFormat: "png",
              wcs: volumeOverlay.wcs,
            });
            aladin.setOverlayImageLayer(imageLayer, String(volumeOverlay.name || "Selected Dust"));
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
        skyReadoutEl.textContent = skyReadoutText(selectionsForPanel, payload, mode, volumeOverlay);
      }

      function updateSelectionUI() {
        selectionReadoutEl.textContent = selectionToolbarText(currentSelections, currentSelection);
        clearSelectionButtonEl.disabled = currentSelections.length === 0 && !currentSelection && !hasActiveLassoSelectionMask();
        lassoButtonEl.dataset.active = lassoArmed ? "true" : "false";
        lassoButtonEl.setAttribute("aria-pressed", lassoArmed ? "true" : "false");
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

      function clamp01(value) {
        return Math.min(Math.max(Number(value), 0.0), 1.0);
      }

      function clampRange(value, minValue, maxValue) {
        const numeric = Number(value);
        if (!Number.isFinite(numeric)) {
          return minValue;
        }
        return Math.min(Math.max(numeric, minValue), maxValue);
      }

      function themePresetForKey(themeKey) {
        const requestedKey = String(themeKey || "default");
        return themePresets[requestedKey] || themePresets.default || {};
      }

      function applyThemePreset(themeKey, options = {}) {
        activeThemeKey = Object.prototype.hasOwnProperty.call(themePresets, String(themeKey))
          ? String(themeKey)
          : "default";
        const nextTheme = safeJsonClone(themePresetForKey(activeThemeKey), {});
        Object.keys(theme).forEach((key) => {
          delete theme[key];
        });
        Object.assign(theme, nextTheme);
        applyThemeCssVars();
        scene.background = new THREE.Color(theme.scene_bgcolor || theme.paper_bgcolor || "#000000");
        if (options.syncInput !== false && themeSelectEl) {
          themeSelectEl.value = activeThemeKey;
        }
        if (options.rerender !== false) {
          buildAxes();
          renderLegend();
          updateSkyPanel();
          renderFrame(currentFrameIndex);
        }
      }

      function applyGlobalControlState() {
        globalScrollSpeed = clampRange(globalScrollSpeed, 0.2, 4.0);
        globalPointSizeScale = clampRange(globalPointSizeScale, 0.25, 4.0);
        globalPointOpacityScale = clampRange(globalPointOpacityScale, 0.0, 2.0);
        fadeInTimeMyr = Math.max(Number.isFinite(Number(fadeInTimeMyr)) ? Number(fadeInTimeMyr) : 0.0, 0.0);
        focusTraceKey = String(focusTraceKey || "");
        camera.fov = clampRange(camera.fov, 18.0, 90.0);
        controls.zoomSpeed = globalScrollSpeed;
        camera.updateProjectionMatrix();
      }

      function applyCameraViewMode() {
        const isEarthView = cameraViewMode === "earth";
        controls.enableRotate = true;
        controls.enableZoom = !isEarthView;
        controls.enablePan = !isEarthView;
      }

      function resetCameraView() {
        cameraViewMode = "free";
        earthViewFocusDistance = null;
        camera.position.copy(initialCameraState.position);
        controls.target.copy(initialCameraState.target);
        camera.up.copy(initialCameraState.up);
        camera.fov = Number(initialCameraState.fov);
        applyGlobalControlState();
        applyCameraViewMode();
        controls.update();
        renderFrame(currentFrameIndex);
      }

      function resetCameraAndSelections() {
        focusSelectionKey = "";
        lassoArmed = false;
        clearClusterSelections();
        resetCameraView();
      }

      function centroidFromSpriteEntries(entries) {
        if (!Array.isArray(entries) || !entries.length) {
          return null;
        }
        let count = 0;
        const centroid = new THREE.Vector3();
        entries.forEach((entry) => {
          if (!entry || !entry.sprite) {
            return;
          }
          const worldPoint = new THREE.Vector3();
          entry.sprite.getWorldPosition(worldPoint);
          if (!Number.isFinite(worldPoint.x) || !Number.isFinite(worldPoint.y) || !Number.isFinite(worldPoint.z)) {
            return;
          }
          centroid.add(worldPoint);
          count += 1;
        });
        if (!count) {
          return null;
        }
        return centroid.multiplyScalar(1.0 / count);
      }

      function selectionFallbackPoint(selection) {
        if (!selection || typeof selection !== "object") {
          return null;
        }
        const x = Number(selection.x0);
        const y = Number(selection.y0);
        const z = Number(selection.z0);
        if (Number.isFinite(x) && Number.isFinite(y) && Number.isFinite(z)) {
          return new THREE.Vector3(x, y, z);
        }
        return null;
      }

      function currentSelectionCentroidWorldPoint() {
        if (currentSelection) {
          const key = normalizedSelectionKeyFor(currentSelection);
          const centroid = key ? centroidFromSpriteEntries(selectionSpriteEntriesByKey.get(key)) : null;
          if (centroid) {
            return centroid;
          }
          return selectionFallbackPoint(currentSelection);
        }

        if (currentSelections.length) {
          let count = 0;
          const centroid = new THREE.Vector3();
          currentSelections.forEach((selection) => {
            const key = normalizedSelectionKeyFor(selection);
            const selectionCentroid = key ? centroidFromSpriteEntries(selectionSpriteEntriesByKey.get(key)) : null;
            const fallbackPoint = selectionCentroid || selectionFallbackPoint(selection);
            if (!fallbackPoint) {
              return;
            }
            centroid.add(fallbackPoint);
            count += 1;
          });
          if (count) {
            return centroid.multiplyScalar(1.0 / count);
          }
        }
        return null;
      }

      function earthViewTargetPoint() {
        const selectedPoint = currentSelectionCentroidWorldPoint();
        if (selectedPoint) {
          return selectedPoint;
        }
        const gcX = Math.max(
          Number((sceneSpec.ranges || {}).x ? sceneSpec.ranges.x[1] : NaN) || 0.0,
          8122.0
        );
        return new THREE.Vector3(gcX, 0.0, 0.0);
      }

      function viewFromEarth() {
        const earthPoint = new THREE.Vector3(0.0, 0.0, 20.8);
        let targetPoint = earthViewTargetPoint();
        if (!targetPoint || !Number.isFinite(targetPoint.x) || !Number.isFinite(targetPoint.y) || !Number.isFinite(targetPoint.z)) {
          targetPoint = new THREE.Vector3(8122.0, 0.0, 0.0);
        }
        const direction = new THREE.Vector3().subVectors(targetPoint, earthPoint);
        if (direction.lengthSq() <= 1e-12) {
          direction.set(1.0, 0.0, 0.0);
          targetPoint = earthPoint.clone().add(direction);
        }
        earthViewFocusDistance = Math.max(direction.length(), 1e-6);
        direction.normalize();
        const orbitRadius = Math.max(1e-3, Math.min(0.05, earthViewFocusDistance * 1e-6));
        cameraViewMode = "earth";
        controls.target.copy(earthPoint);
        camera.position.copy(earthPoint.clone().sub(direction.clone().multiplyScalar(orbitRadius)));
        camera.up.set(sceneUp.x ?? 0.0, sceneUp.y ?? 0.0, sceneUp.z ?? 1.0);
        camera.fov = 90.0;
        applyGlobalControlState();
        applyCameraViewMode();
        controls.update();
        renderSceneControls();
        updateScaleBar();
      }

      function renderSceneControls() {
        if (themeSelectEl) {
          themeSelectEl.value = activeThemeKey;
        }
        if (scrollSpeedEl) {
          scrollSpeedEl.value = String(globalScrollSpeed);
        }
        if (scrollSpeedLabelEl) {
          scrollSpeedLabelEl.textContent = `Scroll speed (${globalScrollSpeed.toFixed(2)}x)`;
        }
        if (cameraFovEl) {
          cameraFovEl.value = String(camera.fov);
        }
        if (cameraFovLabelEl) {
          cameraFovLabelEl.textContent = `Camera FOV (${Math.round(Number(camera.fov))} deg)`;
        }
        if (globalPointSizeEl) {
          globalPointSizeEl.value = String(globalPointSizeScale);
        }
        if (globalPointSizeLabelEl) {
          globalPointSizeLabelEl.textContent = `Point size (${globalPointSizeScale.toFixed(2)}x)`;
        }
        if (globalPointOpacityEl) {
          globalPointOpacityEl.value = String(globalPointOpacityScale);
        }
        if (globalPointOpacityLabelEl) {
          globalPointOpacityLabelEl.textContent = `Point opacity (${globalPointOpacityScale.toFixed(2)}x)`;
        }
        if (focusGroupSelectEl) {
          focusGroupSelectEl.value = focusTraceKey;
        }
        if (fadeTimeEl) {
          fadeTimeEl.value = Number(fadeInTimeMyr).toFixed(1).replace(/\.0$/, "");
        }
        if (fadeInOutToggleEl) {
          fadeInOutToggleEl.checked = fadeInAndOutEnabled;
        }
        if (axesVisibleToggleEl) {
          axesVisibleToggleEl.checked = axesVisible;
        }
      }

      function setZenMode(enabled) {
        zenModeEnabled = Boolean(enabled);
        root.dataset.zen = zenModeEnabled ? "true" : "false";
        if (zenModeButtonEl) {
          zenModeButtonEl.dataset.active = zenModeEnabled ? "true" : "false";
          zenModeButtonEl.textContent = zenModeEnabled ? "Exit Zen" : "Zen";
          zenModeButtonEl.title = zenModeEnabled
            ? "Restore the interface panels and controls"
            : "Hide interface panels and keep only the time slider visible";
        }
        if (zenModeEnabled) {
          tooltipEl.style.display = "none";
        }
      }

      function setKeyHelpOpen(isOpen) {
        if (!keyHelpEl) {
          return;
        }
        keyHelpEl.dataset.open = isOpen ? "true" : "false";
      }

      function normalizedKeyboardKey(key) {
        const text = String(key || "");
        return text.length === 1 ? text.toLowerCase() : text.toLowerCase();
      }

      function focusViewer() {
        const focusTarget = canvas || root;
        try {
          window.focus();
        } catch (_err) {}
        try {
          focusTarget.focus({ preventScroll: true });
        } catch (_err) {
          focusTarget.focus();
        }
      }

      function clearPressedKeys() {
        pressedKeys.clear();
      }

      function keyboardTargetIsEditable(target) {
        if (!target || target === document.body || target === root || target === canvas) {
          return false;
        }
        if (typeof target.closest === "function" && target.closest(".oviz-three-key-help")) {
          return true;
        }
        if (target.isContentEditable) {
          return true;
        }
        const tagName = String(target.tagName || "").toLowerCase();
        return ["input", "select", "textarea", "button"].includes(tagName);
      }

      function keyboardLegendItems() {
        const defaults = groupDefaults(currentGroup);
        return legendItems.filter((item) => {
          const mode = defaults[item.key];
          return !(mode === false || mode === undefined);
        });
      }

      function toggleLegendItemByIndex(index, solo = false) {
        const items = keyboardLegendItems();
        if (index < 0 || index >= items.length) {
          return false;
        }
        if (solo) {
          items.forEach((item, itemIndex) => {
            legendState[item.key] = itemIndex === index;
          });
        } else {
          const item = items[index];
          legendState[item.key] = !legendState[item.key];
        }
        renderLegend();
        renderFrame(currentFrameIndex);
        return true;
      }

      function soloTraceLegendItem(itemKey) {
        const targetKey = String(itemKey || "");
        if (!targetKey || volumeLayerForKey(targetKey)) {
          return false;
        }
        const defaults = groupDefaults(currentGroup);
        let foundTarget = false;
        legendItems.forEach((item) => {
          const key = String(item.key || "");
          const mode = defaults[key];
          if (mode === false || mode === undefined || volumeLayerForKey(key)) {
            return;
          }
          legendState[key] = key === targetKey;
          if (key === targetKey) {
            foundTarget = true;
          }
        });
        if (!foundTarget) {
          return false;
        }
        renderLegend();
        renderFrame(currentFrameIndex);
        return true;
      }

      function cameraTravelDistance(fast = false) {
        const referenceDistance = cameraViewMode === "earth" && Number.isFinite(earthViewFocusDistance) && earthViewFocusDistance > 0.0
          ? earthViewFocusDistance
          : Math.max(camera.position.distanceTo(controls.target), 1.0);
        const step = clampRange(referenceDistance * 0.01, 1.0, Math.max((sceneSpec.max_span || 1) * 0.03, 1.0));
        return fast ? step * 4.0 : step;
      }

      function enterFreeCameraMode() {
        if (cameraViewMode !== "free") {
          cameraViewMode = "free";
          earthViewFocusDistance = null;
          applyCameraViewMode();
        }
      }

      function translateCameraAndTarget(delta) {
        if (!delta || delta.lengthSq() <= 1e-18) {
          return;
        }
        enterFreeCameraMode();
        camera.position.add(delta);
        controls.target.add(delta);
        controls.update();
        updateScaleBar();
      }

      function rotateCameraYaw(sign, fast = false) {
        const angle = (fast ? 0.12 : 0.04) * sign;
        const offset = camera.position.clone().sub(controls.target);
        if (offset.lengthSq() <= 1e-18) {
          return;
        }
        offset.applyAxisAngle(sceneUpVector, angle);
        camera.position.copy(controls.target.clone().add(offset));
        camera.up.set(sceneUp.x ?? 0.0, sceneUp.y ?? 0.0, sceneUp.z ?? 1.0);
        controls.update();
      }

      function stepFrame(delta) {
        if (!delta) {
          return;
        }
        pause();
        const nextIndex = Math.max(0, Math.min(currentFrameIndex + delta, frameSpecs.length - 1));
        if (nextIndex !== currentFrameIndex) {
          renderFrame(nextIndex);
        }
      }

      function movementVectors() {
        const forward = new THREE.Vector3();
        camera.getWorldDirection(forward);
        if (forward.lengthSq() <= 1e-18) {
          forward.set(1.0, 0.0, 0.0);
        } else {
          forward.normalize();
        }
        const right = new THREE.Vector3().crossVectors(forward, sceneUpVector).normalize();
        if (right.lengthSq() <= 1e-18) {
          right.set(0.0, 1.0, 0.0);
        }
        const up = sceneUpVector.clone();
        return { forward, right, up };
      }

      function orbitCameraByKeyboard(deltaTheta, deltaPhi) {
        if ((!Number.isFinite(deltaTheta) || Math.abs(deltaTheta) <= 1e-12) && (!Number.isFinite(deltaPhi) || Math.abs(deltaPhi) <= 1e-12)) {
          return;
        }
        const offset = camera.position.clone().sub(controls.target);
        if (offset.lengthSq() <= 1e-18) {
          return;
        }
        const quat = new THREE.Quaternion().setFromUnitVectors(sceneUpVector, new THREE.Vector3(0.0, 1.0, 0.0));
        const quatInverse = quat.clone().invert();
        offset.applyQuaternion(quat);
        const spherical = new THREE.Spherical().setFromVector3(offset);
        spherical.theta += Number.isFinite(deltaTheta) ? deltaTheta : 0.0;
        spherical.phi += Number.isFinite(deltaPhi) ? deltaPhi : 0.0;
        const minPolar = Number.isFinite(controls.minPolarAngle) ? controls.minPolarAngle : 0.02;
        const maxPolar = Number.isFinite(controls.maxPolarAngle) ? controls.maxPolarAngle : (Math.PI - 0.02);
        spherical.phi = clampRange(spherical.phi, minPolar, maxPolar);
        if (typeof spherical.makeSafe === "function") {
          spherical.makeSafe();
        }
        offset.setFromSpherical(spherical);
        offset.applyQuaternion(quatInverse);
        camera.position.copy(controls.target.clone().add(offset));
        camera.up.set(sceneUp.x ?? 0.0, sceneUp.y ?? 0.0, sceneUp.z ?? 1.0);
      }

      function zoomCameraByKeyboard(sign, deltaSeconds) {
        if (!Number.isFinite(sign) || Math.abs(sign) <= 1e-12 || !Number.isFinite(deltaSeconds) || deltaSeconds <= 0.0) {
          return;
        }
        const offset = camera.position.clone().sub(controls.target);
        const distance = offset.length();
        if (!(distance > 1e-12)) {
          return;
        }
        const speedScale = Math.max(globalScrollSpeed, 0.2);
        const zoomFactor = Math.exp(sign * deltaSeconds * 1.8 * speedScale);
        const maxSpan = Math.max(sceneSpec.max_span || 1, 1);
        const minDistance = cameraViewMode === "earth"
          ? Math.max(1e-6, Math.min(0.005, distance * 0.2))
          : Math.max(maxSpan * 1e-5, 0.05);
        const maxDistance = Math.max(maxSpan * 12.0, 10.0);
        const nextDistance = clampRange(distance * zoomFactor, minDistance, maxDistance);
        offset.setLength(nextDistance);
        camera.position.copy(controls.target.clone().add(offset));
        camera.up.set(sceneUp.x ?? 0.0, sceneUp.y ?? 0.0, sceneUp.z ?? 1.0);
        updateScaleBar();
      }

      function updateKeyboardMotion(deltaSeconds) {
        if (!Number.isFinite(deltaSeconds) || deltaSeconds <= 0.0 || pressedKeys.size === 0) {
          return;
        }
        if (keyboardTargetIsEditable(document.activeElement) || (keyHelpEl && keyHelpEl.dataset.open === "true")) {
          return;
        }

        const shiftHeld = pressedKeys.has("shift");
        const speedScale = Math.max(globalScrollSpeed, 0.2);
        const panDistance = cameraTravelDistance(shiftHeld) * deltaSeconds * 6.0 * speedScale;
        const verticalDistance = cameraTravelDistance(shiftHeld) * deltaSeconds * 5.0 * speedScale;
        const orbitAzimuthSpeed = deltaSeconds * 1.4 * speedScale;
        const orbitPolarSpeed = deltaSeconds * 1.15 * speedScale;
        const { forward, right, up } = movementVectors();

        if (shiftHeld) {
          const translation = new THREE.Vector3();
          if (pressedKeys.has("w")) {
            translation.add(forward.clone().multiplyScalar(panDistance));
          }
          if (pressedKeys.has("s")) {
            translation.add(forward.clone().multiplyScalar(-panDistance));
          }
          if (pressedKeys.has("a")) {
            translation.add(right.clone().multiplyScalar(-panDistance));
          }
          if (pressedKeys.has("d")) {
            translation.add(right.clone().multiplyScalar(panDistance));
          }
          if (translation.lengthSq() > 1e-18) {
            translateCameraAndTarget(translation);
          }
        } else {
          let deltaTheta = 0.0;
          let deltaPhi = 0.0;
          if (pressedKeys.has("a")) {
            deltaTheta += orbitAzimuthSpeed;
          }
          if (pressedKeys.has("d")) {
            deltaTheta -= orbitAzimuthSpeed;
          }
          if (pressedKeys.has("w")) {
            deltaPhi -= orbitPolarSpeed;
          }
          if (pressedKeys.has("s")) {
            deltaPhi += orbitPolarSpeed;
          }
          if (Math.abs(deltaTheta) > 1e-12 || Math.abs(deltaPhi) > 1e-12) {
            orbitCameraByKeyboard(deltaTheta, deltaPhi);
          }
        }

        if (pressedKeys.has("q")) {
          zoomCameraByKeyboard(1.0, deltaSeconds);
        }
        if (pressedKeys.has("e")) {
          zoomCameraByKeyboard(-1.0, deltaSeconds);
        }

        let verticalTranslation = null;
        if (pressedKeys.has("r")) {
          verticalTranslation = up.clone().multiplyScalar(verticalDistance);
        } else if (pressedKeys.has("f")) {
          verticalTranslation = up.clone().multiplyScalar(-verticalDistance);
        }
        if (verticalTranslation) {
          translateCameraAndTarget(verticalTranslation);
        }
      }

      function onKeyDown(event) {
        if (keyboardTargetIsEditable(event.target)) {
          return;
        }
        if (event.metaKey || event.ctrlKey || event.altKey) {
          return;
        }

        const key = String(event.key || "");
        const lowerKey = normalizedKeyboardKey(key);
        const fast = Boolean(event.shiftKey);
        const isMovementKey = ["w", "a", "s", "d", "q", "e", "r", "f", "shift"].includes(lowerKey);

        if (keyHelpEl && keyHelpEl.dataset.open === "true") {
          if (key === "Escape" || key === "?" || (key === "/" && event.shiftKey)) {
            clearPressedKeys();
            setKeyHelpOpen(false);
            focusViewer();
            event.preventDefault();
          }
          return;
        }

        if (isMovementKey) {
          pressedKeys.add(lowerKey);
          event.preventDefault();
          return;
        }

        if (event.repeat && (key === " " || key === "Escape" || /^[1-9]$/.test(key) || lowerKey === "l" || lowerKey === "c" || lowerKey === "v" || key === "?" || (key === "/" && event.shiftKey))) {
          event.preventDefault();
          return;
        }

        if (key === "?" || (key === "/" && event.shiftKey)) {
          clearPressedKeys();
          setKeyHelpOpen(true);
          event.preventDefault();
          return;
        }

        if (key === "Escape") {
          clearPressedKeys();
          clearClusterSelections();
          lassoArmed = false;
          updateSelectionUI();
          event.preventDefault();
          return;
        }

        if (key === " ") {
          if (playbackTimer) {
            pause();
          } else {
            play();
          }
          event.preventDefault();
          return;
        }

        if (key === "ArrowLeft") {
          stepFrame(fast ? -5 : -1);
          event.preventDefault();
          return;
        }
        if (key === "ArrowRight") {
          stepFrame(fast ? 5 : 1);
          event.preventDefault();
          return;
        }

        if (/^[1-9]$/.test(key)) {
          const handled = toggleLegendItemByIndex(Number(key) - 1, fast);
          if (handled) {
            event.preventDefault();
          }
          return;
        }

        if (lowerKey === "l") {
          lassoArmed = !lassoArmed;
          updateSelectionUI();
          event.preventDefault();
          return;
        }
        if (lowerKey === "c") {
          clickSelectionEnabled = !clickSelectionEnabled;
          updateSelectionUI();
          event.preventDefault();
          return;
        }
        if (lowerKey === "v") {
          viewFromEarth();
          event.preventDefault();
          return;
        }
      }

      function onKeyUp(event) {
        if (keyboardTargetIsEditable(event.target)) {
          return;
        }
        if (event.metaKey || event.ctrlKey || event.altKey) {
          return;
        }
        const lowerKey = normalizedKeyboardKey(event.key || "");
        if (pressedKeys.has(lowerKey)) {
          pressedKeys.delete(lowerKey);
          event.preventDefault();
        }
      }

      function cssColorToHex(value, fallback = "#ffffff") {
        try {
          return `#${new THREE.Color(value || fallback).getHexString()}`;
        } catch (_err) {
          return fallback;
        }
      }

      function cssColorWithAlpha(value, alpha, fallback = "#ffffff") {
        try {
          const color = new THREE.Color(value || fallback);
          const r = Math.round(color.r * 255.0);
          const g = Math.round(color.g * 255.0);
          const b = Math.round(color.b * 255.0);
          return `rgba(${r}, ${g}, ${b}, ${clamp01(alpha)})`;
        } catch (_err) {
          return fallback;
        }
      }

      function formatCompactNumber(value) {
        const num = Number(value);
        if (!Number.isFinite(num)) {
          return "";
        }
        const abs = Math.abs(num);
        if (abs >= 100.0) {
          return String(Math.round(num));
        }
        if (abs >= 10.0) {
          return num.toFixed(1).replace(/\.0$/, "");
        }
        if (abs >= 1.0) {
          return num.toFixed(2).replace(/0$/, "").replace(/\.$/, "");
        }
        return num.toPrecision(2).replace(/\.0+e/, "e");
      }

      function formatDistanceLabelPc(distancePc) {
        const value = Number(distancePc);
        if (!Number.isFinite(value) || value <= 0.0) {
          return "";
        }
        if (value >= 1000.0) {
          return `${formatCompactNumber(value / 1000.0)} kpc`;
        }
        return `${formatCompactNumber(value)} pc`;
      }

      function updateScaleBar() {
        if (!scaleBarEl || !scaleLabelEl) {
          return;
        }
        const canvasHeight = Math.max(canvas.clientHeight || root.clientHeight || 0, 1);
        const distance = cameraViewMode === "earth" && Number.isFinite(earthViewFocusDistance) && earthViewFocusDistance > 0.0
          ? earthViewFocusDistance
          : Math.max(camera.position.distanceTo(controls.target), 1e-6);
        const worldPerPixel = (2.0 * distance * Math.tan(THREE.MathUtils.degToRad(camera.fov * 0.5))) / canvasHeight;
        const barLengthPc = worldPerPixel * 120.0;
        scaleLabelEl.textContent = formatDistanceLabelPc(barLengthPc);
        scaleBarEl.style.display = Number.isFinite(barLengthPc) ? "flex" : "none";
      }

      function clusterFilterParameterSpecForKey(parameterKey) {
        return clusterFilterParameters.find((parameter) => String(parameter.key) === String(parameterKey)) || null;
      }

      function activeClusterFilterParameterSpec() {
        return clusterFilterParameterSpecForKey(clusterFilterParameterKey)
          || (clusterFilterParameters.length ? clusterFilterParameters[0] : null);
      }

      function clampClusterFilterRangeForParameter(parameter) {
        if (!parameter) {
          return null;
        }
        const key = String(parameter.key || "");
        const parameterMin = Number(parameter.min);
        const parameterMax = Number(parameter.max);
        const rangeState = clusterFilterRangeStateByKey[key] || {
          min: parameterMin,
          max: parameterMax,
        };
        let minValue = Number(rangeState.min);
        let maxValue = Number(rangeState.max);
        if (!Number.isFinite(minValue)) {
          minValue = parameterMin;
        }
        if (!Number.isFinite(maxValue)) {
          maxValue = parameterMax;
        }
        minValue = clampRange(minValue, parameterMin, parameterMax);
        maxValue = clampRange(maxValue, parameterMin, parameterMax);
        if (minValue > maxValue) {
          const middle = 0.5 * (minValue + maxValue);
          minValue = middle;
          maxValue = middle;
        }
        rangeState.min = minValue;
        rangeState.max = maxValue;
        clusterFilterRangeStateByKey[key] = rangeState;
        return rangeState;
      }

      function formatClusterFilterValue(value, parameter) {
        const numericValue = Number(value);
        if (!Number.isFinite(numericValue)) {
          return "";
        }
        const unit = String((parameter && parameter.unit) || "").trim();
        const valueText = numericValue >= 1000.0
          ? formatCompactNumber(numericValue)
          : formatCompactNumber(numericValue);
        return unit ? `${valueText} ${unit}` : valueText;
      }

      function clusterFilterEntryValue(entry, parameterKey) {
        if (!entry || typeof entry !== "object") {
          return NaN;
        }
        return Number(entry[String(parameterKey || "")]);
      }

      function clusterFilterSelectionKeyForPoint(point) {
        if (point && point.motion && point.motion.key) {
          return normalizeMemberKey(point.motion.key);
        }
        if (point && point.selection) {
          return normalizedSelectionKeyFor(point.selection);
        }
        return "";
      }

      function clusterFilterPassesSelectionKey(selectionKey) {
        if (!clusterFilterSpec.enabled) {
          return true;
        }
        const key = normalizeMemberKey(selectionKey);
        if (!key) {
          return true;
        }
        const parameter = activeClusterFilterParameterSpec();
        if (!parameter) {
          return true;
        }
        const rangeState = clampClusterFilterRangeForParameter(parameter);
        const entry = clusterFilterEntryByKey.get(key);
        if (!entry || !rangeState) {
          return true;
        }
        const value = clusterFilterEntryValue(entry, parameter.key);
        if (!Number.isFinite(value)) {
          return true;
        }
        return value >= Number(rangeState.min) && value <= Number(rangeState.max);
      }

      function clusterFilterPassesSelectionKeyExcludingParameter(selectionKey, excludedParameterKey) {
        if (!clusterFilterSpec.enabled) {
          return true;
        }
        const key = normalizeMemberKey(selectionKey);
        if (!key) {
          return true;
        }
        const parameter = activeClusterFilterParameterSpec();
        if (!parameter) {
          return true;
        }
        if (excludedParameterKey && String(parameter.key) === String(excludedParameterKey)) {
          return true;
        }
        const rangeState = clampClusterFilterRangeForParameter(parameter);
        const entry = clusterFilterEntryByKey.get(key);
        if (!entry || !rangeState) {
          return true;
        }
        const value = clusterFilterEntryValue(entry, parameter.key);
        if (!Number.isFinite(value)) {
          return true;
        }
        return value >= Number(rangeState.min) && value <= Number(rangeState.max);
      }

      function clusterFilterPassesSelection(selection) {
        return clusterFilterPassesSelectionKey(normalizedSelectionKeyFor(selection));
      }

      function clusterFilterPassesPoint(point) {
        return clusterFilterPassesSelectionKey(clusterFilterSelectionKeyForPoint(point));
      }

      function pruneSelectionsToActiveClusterFilter() {
        currentSelections = uniqueSelections(currentSelections.filter((selection) => clusterFilterPassesSelection(selection)));
        if (currentSelection && !clusterFilterPassesSelection(currentSelection)) {
          currentSelection = null;
        }
        const nextSelectedKeys = new Set();
        currentSelections.forEach((selection) => {
          const key = normalizedSelectionKeyFor(selection);
          if (key) {
            nextSelectedKeys.add(key);
          }
        });
        if (currentSelection) {
          const focusKey = normalizedSelectionKeyFor(currentSelection);
          if (focusKey) {
            nextSelectedKeys.add(focusKey);
          }
        }
        selectedClusterKeys = nextSelectedKeys;
        currentSelectionMode = currentSelection ? "click" : ((currentSelections.length || hasActiveLassoSelectionMask()) ? "lasso" : "none");
        if (!clusterFilterPassesSelectionKey(localHoveredClusterKey)) {
          setLocalHoveredClusterKey("");
        }
        if (!clusterFilterPassesSelectionKey(skyHoveredClusterKey)) {
          skyHoveredClusterKey = "";
          lastSentSkyHoverClusterKey = null;
        }
      }

      function ageKdeGridValues() {
        const traces = Array.isArray(ageKdeSpec.traces) ? ageKdeSpec.traces : [];
        if (traces.length && Array.isArray(traces[0].x) && traces[0].x.length >= 2) {
          return traces[0].x.map((value) => Number(value)).filter(Number.isFinite);
        }
        const xRange = Array.isArray(ageKdeSpec.x_range) ? ageKdeSpec.x_range : [-1.0, 0.0];
        const xMin = Number(xRange[0]);
        const xMax = Number(xRange[1]);
        const count = 300;
        return Array.from({ length: count }, (_, index) => {
          const frac = count <= 1 ? 0.0 : index / (count - 1);
          return xMin + frac * (xMax - xMin);
        });
      }

      function gaussianKde(values, xGrid, bandwidth) {
        const finiteValues = (Array.isArray(values) ? values : []).map((value) => Number(value)).filter(Number.isFinite);
        if (!finiteValues.length) {
          return xGrid.map(() => 0.0);
        }
        const bw = Math.max(Number(bandwidth) || 1.0, 1e-6);
        const norm = finiteValues.length * bw * Math.sqrt(2.0 * Math.PI);
        return xGrid.map((xValue) => {
          let sum = 0.0;
          for (const value of finiteValues) {
            const u = (Number(xValue) - value) / bw;
            sum += Math.exp(-0.5 * u * u);
          }
          return sum / norm;
        });
      }

      function ageKdeFilterParameterSpec() {
        return clusterFilterParameterSpecForKey("age_now_myr");
      }

      function ageKdeAxisRange() {
        const xRange = Array.isArray(ageKdeSpec.x_range) ? ageKdeSpec.x_range : [-1.0, 0.0];
        let minValue = Number(xRange[0]);
        let maxValue = Number(xRange[1]);
        if (!Number.isFinite(minValue)) {
          minValue = -1.0;
        }
        if (!Number.isFinite(maxValue)) {
          maxValue = 0.0;
        }
        if (minValue > maxValue) {
          const swap = minValue;
          minValue = maxValue;
          maxValue = swap;
        }
        return { min: minValue, max: maxValue };
      }

      function ageNowToKdeAxisValue(ageNowMyr) {
        const value = Number(ageNowMyr);
        if (!Number.isFinite(value)) {
          return NaN;
        }
        return -Math.abs(value);
      }

      function kdeAxisValueToAgeNow(axisValue) {
        const value = Number(axisValue);
        if (!Number.isFinite(value)) {
          return NaN;
        }
        return Math.abs(value);
      }

      function ageKdeAxisValueToSlider(axisValue) {
        const axisRange = ageKdeAxisRange();
        const denom = Math.max(axisRange.max - axisRange.min, 1e-9);
        return Math.round(1000.0 * clampRange((Number(axisValue) - axisRange.min) / denom, 0.0, 1.0));
      }

      function ageKdeSliderValueToAxisValue(sliderValue) {
        const axisRange = ageKdeAxisRange();
        const normalized = clampRange(Number(sliderValue) / 1000.0, 0.0, 1.0);
        return axisRange.min + normalized * (axisRange.max - axisRange.min);
      }

      function ageKdeAxisFilterRange() {
        const parameter = ageKdeFilterParameterSpec();
        if (!parameter) {
          return null;
        }
        const rangeState = clampClusterFilterRangeForParameter(parameter);
        const axisRange = ageKdeAxisRange();
        let minAxis = clampRange(ageNowToKdeAxisValue(rangeState.max), axisRange.min, axisRange.max);
        let maxAxis = clampRange(ageNowToKdeAxisValue(rangeState.min), axisRange.min, axisRange.max);
        if (minAxis > maxAxis) {
          const swap = minAxis;
          minAxis = maxAxis;
          maxAxis = swap;
        }
        return { min: minAxis, max: maxAxis };
      }

      function formatAgeKdeAxisValue(axisValue) {
        const value = Number(axisValue);
        if (!Number.isFinite(value)) {
          return "";
        }
        return `${formatCompactNumber(value)} Myr`;
      }

      function setClusterAgeFilterFromKdeAxisRange(minAxisValue, maxAxisValue) {
        const parameter = ageKdeFilterParameterSpec();
        if (!parameter) {
          return;
        }
        const axisRange = ageKdeAxisRange();
        const clampedMinAxis = clampRange(Math.min(Number(minAxisValue), Number(maxAxisValue)), axisRange.min, axisRange.max);
        const clampedMaxAxis = clampRange(Math.max(Number(minAxisValue), Number(maxAxisValue)), axisRange.min, axisRange.max);
        const youngerAge = clampRange(
          kdeAxisValueToAgeNow(clampedMaxAxis),
          Number(parameter.min),
          Number(parameter.max),
        );
        const olderAge = clampRange(
          kdeAxisValueToAgeNow(clampedMinAxis),
          Number(parameter.min),
          Number(parameter.max),
        );
        clusterFilterParameterKey = String(parameter.key || "");
        clusterFilterRangeStateByKey[String(parameter.key)] = {
          min: Math.min(youngerAge, olderAge),
          max: Math.max(youngerAge, olderAge),
        };
        applyClusterFilterState();
      }

      function ageKdeBaseEntries() {
        const ageParameter = ageKdeFilterParameterSpec();
        if (!ageParameter) {
          return [];
        }
        const selectedKeys = selectedClusterKeys.size ? selectedClusterKeys : null;
        return clusterFilterEntries.filter((entry) => {
          if (!entry || typeof entry !== "object") {
            return false;
          }
          const selectionKey = normalizeMemberKey(entry.selection_key || "");
          if (!selectionKey) {
            return false;
          }
          if (!clusterFilterEntryVisibleInScene(entry)) {
            return false;
          }
          if (selectedKeys && selectedKeys.size && !selectedKeys.has(selectionKey)) {
            return false;
          }
          if (!clusterFilterPassesSelectionKeyExcludingParameter(selectionKey, ageParameter.key)) {
            return false;
          }
          return Number.isFinite(clusterFilterEntryValue(entry, ageParameter.key));
        });
      }

      function filteredAgeKdeSeries() {
        const defaults = groupDefaults(currentGroup);
        const selectedKeys = selectedClusterKeys.size ? selectedClusterKeys : null;
        const traceMetaByKey = new Map();
        const traceMetaByName = new Map();
        (ageKdeSpec.traces || []).forEach((traceSpec) => {
          if (traceSpec.trace_key) {
            traceMetaByKey.set(String(traceSpec.trace_key), traceSpec);
          }
          if (traceSpec.trace_name) {
            traceMetaByName.set(String(traceSpec.trace_name), traceSpec);
          }
        });

        const groupedAges = new Map();
        (ageKdeSpec.cluster_points || []).forEach((point) => {
          const traceKey = point && point.trace_key ? String(point.trace_key) : null;
          const traceName = point && point.trace_name ? String(point.trace_name) : "";
          const traceLookupKey = traceKey || traceName;
          if (!traceLookupKey) {
            return;
          }
          if (traceKey) {
            const mode = defaults[traceKey];
            if (mode !== true || legendState[traceKey] === false) {
              return;
            }
          }
          if (selectedKeys && selectedKeys.size) {
            const selectionLike = {
              cluster_name: point.cluster_name,
              trace_name: traceName,
            };
            const pointKey = normalizedSelectionKeyFor(selectionLike);
            if (!pointKey || !selectedKeys.has(pointKey)) {
              return;
            }
          }
          const pointSelectionKey = normalizedSelectionKeyFor({
            cluster_name: point.cluster_name,
            trace_name: traceName,
          });
          if (!clusterFilterPassesSelectionKey(pointSelectionKey)) {
            return;
          }
          const ageNow = Number(point.age_now_myr);
          if (!Number.isFinite(ageNow)) {
            return;
          }
          if (!groupedAges.has(traceLookupKey)) {
            groupedAges.set(traceLookupKey, []);
          }
          groupedAges.get(traceLookupKey).push(-Math.abs(ageNow));
        });

        const xGrid = ageKdeGridValues();
        const series = [];
        groupedAges.forEach((lookbackValues, traceLookupKey) => {
          if (!lookbackValues.length) {
            return;
          }
          const meta = traceMetaByKey.get(String(traceLookupKey)) || traceMetaByName.get(String(traceLookupKey)) || {};
          const traceKey = meta.trace_key ? String(meta.trace_key) : (String(traceLookupKey).startsWith("trace-") ? String(traceLookupKey) : null);
          const styleState = traceKey ? traceStyleStateForKey(traceKey) : null;
          const densityRaw = gaussianKde(lookbackValues, xGrid, ageKdeSpec.bandwidth_myr);
          const densityMax = Math.max(...densityRaw, 0.0);
          const density = densityMax > 0.0 ? densityRaw.map((value) => value / densityMax) : densityRaw;
          series.push({
            traceName: String(meta.trace_name || traceLookupKey),
            traceKey,
            color: styleState && styleState.color ? styleState.color : (meta.color || axisSpec.x?.linecolor || theme.axis_color || "#808080"),
            opacity: clamp01((styleState && Number.isFinite(styleState.opacity) ? styleState.opacity : 1.0) * Number(meta.opacity ?? 1.0)),
            x: xGrid,
            y: density,
            count: lookbackValues.length,
          });
        });

        series.sort((a, b) => String(a.traceName).localeCompare(String(b.traceName)));
        return series;
      }

      function renderAgeKdeWidget() {
        if (!ageKdeSpec.enabled || !ageKdeCanvasEl || widgetModeForKey("age_kde") === "hidden") {
          return;
        }

        const rect = ageKdeCanvasEl.getBoundingClientRect();
        const cssWidth = Math.max(1, Math.round(rect.width));
        const cssHeight = Math.max(1, Math.round(rect.height));
        const dpr = Math.max(window.devicePixelRatio || 1, 1);
        if (ageKdeCanvasEl.width !== Math.round(cssWidth * dpr) || ageKdeCanvasEl.height !== Math.round(cssHeight * dpr)) {
          ageKdeCanvasEl.width = Math.round(cssWidth * dpr);
          ageKdeCanvasEl.height = Math.round(cssHeight * dpr);
        }

        const ctx = ageKdeCanvasEl.getContext("2d");
        if (!ctx) {
          return;
        }
        ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
        ctx.clearRect(0, 0, cssWidth, cssHeight);
        ctx.fillStyle = theme.scene_bgcolor || theme.paper_bgcolor || "#000000";
        ctx.fillRect(0, 0, cssWidth, cssHeight);

        const margin = { left: 44, right: 16, top: 18, bottom: 28 };
        const plotWidth = Math.max(40, cssWidth - margin.left - margin.right);
        const plotHeight = Math.max(40, cssHeight - margin.top - margin.bottom);
        const xRange = Array.isArray(ageKdeSpec.x_range) ? ageKdeSpec.x_range : [-1.0, 0.0];
        const xMin = Number(xRange[0]);
        const xMax = Number(xRange[1]);
        const axisColor = String(ageKdeSpec.axis_color || theme.axis_color || "#808080");
        const filteredSeries = filteredAgeKdeSeries();
        const ageFilterParameter = ageKdeFilterParameterSpec();
        const ageRangeState = ageFilterParameter ? clampClusterFilterRangeForParameter(ageFilterParameter) : null;
        const ageAxisRangeState = ageFilterParameter ? ageKdeAxisFilterRange() : null;
        const ageEntries = ageKdeBaseEntries();
        const yMax = Math.max(
          ...filteredSeries.map((series) => Math.max(...series.y, 0.0)),
          Number((Array.isArray(ageKdeSpec.y_range) ? ageKdeSpec.y_range[1] : 1.0)) || 1.0,
          1e-6,
        );

        function xToPx(value) {
          const denom = Math.max(xMax - xMin, 1e-6);
          return margin.left + ((Number(value) - xMin) / denom) * plotWidth;
        }

        function yToPx(value) {
          return margin.top + plotHeight - (Math.max(Number(value), 0.0) / yMax) * plotHeight;
        }

        ctx.strokeStyle = axisColor;
        ctx.lineWidth = 1.5;
        ctx.beginPath();
        ctx.moveTo(margin.left, margin.top);
        ctx.lineTo(margin.left, margin.top + plotHeight);
        ctx.lineTo(margin.left + plotWidth, margin.top + plotHeight);
        ctx.stroke();

        ctx.fillStyle = axisColor;
        ctx.font = "11px Helvetica, Arial, sans-serif";
        ctx.textBaseline = "middle";
        ctx.textAlign = "right";
        [0.0, 0.5, 1.0].forEach((frac) => {
          const yValue = frac * yMax;
          const yPx = yToPx(yValue);
          ctx.globalAlpha = 0.16;
          ctx.beginPath();
          ctx.moveTo(margin.left, yPx);
          ctx.lineTo(margin.left + plotWidth, yPx);
          ctx.stroke();
          ctx.globalAlpha = 1.0;
          ctx.fillText(formatCompactNumber(yValue), margin.left - 8, yPx);
        });

        ctx.textAlign = "center";
        ctx.textBaseline = "top";
        [xMin, 0.5 * (xMin + xMax), xMax].forEach((xValue) => {
          const xPx = xToPx(xValue);
          ctx.fillText(formatCompactNumber(xValue), xPx, margin.top + plotHeight + 8);
        });

        const visibleTraceNames = [];
        filteredSeries.forEach((traceSpec) => {
          const lineColor = traceSpec.color || axisColor;
          const lineOpacity = clamp01(traceSpec.opacity);
          const xValues = Array.isArray(traceSpec.x) ? traceSpec.x : [];
          const yValues = Array.isArray(traceSpec.y) ? traceSpec.y : [];
          if (xValues.length < 2 || yValues.length !== xValues.length) {
            return;
          }
          visibleTraceNames.push(String(traceSpec.traceName || traceSpec.traceKey || "trace"));

          ctx.beginPath();
          ctx.moveTo(xToPx(xValues[0]), margin.top + plotHeight);
          for (let i = 0; i < xValues.length; i += 1) {
            ctx.lineTo(xToPx(xValues[i]), yToPx(yValues[i]));
          }
          ctx.lineTo(xToPx(xValues[xValues.length - 1]), margin.top + plotHeight);
          ctx.closePath();
          ctx.fillStyle = cssColorWithAlpha(lineColor, Math.min(0.28, 0.14 + 0.18 * lineOpacity), axisColor);
          ctx.fill();

          ctx.beginPath();
          for (let i = 0; i < xValues.length; i += 1) {
            const xPx = xToPx(xValues[i]);
            const yPx = yToPx(yValues[i]);
            if (i === 0) {
              ctx.moveTo(xPx, yPx);
            } else {
              ctx.lineTo(xPx, yPx);
            }
          }
          ctx.strokeStyle = cssColorWithAlpha(lineColor, lineOpacity, axisColor);
          ctx.lineWidth = 2.0;
          ctx.stroke();
        });

        const frame = currentFrame();
        const timeValue = frame ? Number(frame.time) : 0.0;
        const markerX = xToPx(Math.min(Math.max(timeValue, xMin), xMax));
        ctx.save();
        ctx.setLineDash([6, 6]);
        ctx.strokeStyle = axisColor;
        ctx.lineWidth = 1.5;
        ctx.beginPath();
        ctx.moveTo(markerX, margin.top);
        ctx.lineTo(markerX, margin.top + plotHeight);
        ctx.stroke();
        ctx.restore();

        if (ageFilterParameter && ageRangeState && ageAxisRangeState) {
          if (ageKdeFilterRangeMinEl) {
            ageKdeFilterRangeMinEl.value = String(ageKdeAxisValueToSlider(ageAxisRangeState.min));
          }
          if (ageKdeFilterRangeMaxEl) {
            ageKdeFilterRangeMaxEl.value = String(ageKdeAxisValueToSlider(ageAxisRangeState.max));
          }
          if (ageKdeFilterRangeReadoutMinEl) {
            ageKdeFilterRangeReadoutMinEl.textContent = formatAgeKdeAxisValue(ageAxisRangeState.min);
          }
          if (ageKdeFilterRangeReadoutMaxEl) {
            ageKdeFilterRangeReadoutMaxEl.textContent = formatAgeKdeAxisValue(ageAxisRangeState.max);
          }
        }

        ageKdeReadoutEl.textContent = [
          `${String(ageKdeSpec.title || "Relative SFH")}`,
          `Time: ${frame ? frame.name : "0"} Myr`,
          `Bandwidth: ${formatCompactNumber(ageKdeSpec.bandwidth_myr)} Myr`,
          ageFilterParameter && ageAxisRangeState
            ? `Age filter: ${formatAgeKdeAxisValue(ageAxisRangeState.min)} to ${formatAgeKdeAxisValue(ageAxisRangeState.max)}`
            : "Age filter: unavailable",
          `Traces shown: ${visibleTraceNames.length}`,
          `Clusters contributing: ${filteredSeries.reduce((total, series) => total + Number(series.count || 0), 0)}`,
          ageFilterParameter
            ? `Visible age candidates: ${ageEntries.length}`
            : null,
        ].filter(Boolean).join("\\n");
      }

      function clusterFilterSliderValueToActual(sliderValue, parameter) {
        const spec = parameter || activeClusterFilterParameterSpec();
        if (!spec) {
          return 0.0;
        }
        const minValue = Number(spec.min);
        const maxValue = Number(spec.max);
        if (!Number.isFinite(minValue) || !Number.isFinite(maxValue) || Math.abs(maxValue - minValue) <= 1e-12) {
          return minValue;
        }
        const normalized = clampRange(Number(sliderValue) / 1000.0, 0.0, 1.0);
        return minValue + normalized * (maxValue - minValue);
      }

      function clusterFilterActualValueToSlider(actualValue, parameter) {
        const spec = parameter || activeClusterFilterParameterSpec();
        if (!spec) {
          return 0;
        }
        const minValue = Number(spec.min);
        const maxValue = Number(spec.max);
        if (!Number.isFinite(minValue) || !Number.isFinite(maxValue) || Math.abs(maxValue - minValue) <= 1e-12) {
          return 0;
        }
        return Math.round(1000.0 * clampRange((Number(actualValue) - minValue) / (maxValue - minValue), 0.0, 1.0));
      }

      function clusterFilterEntryVisibleInScene(entry) {
        if (!entry || typeof entry !== "object") {
          return false;
        }
        const traceKey = entry.trace_key ? String(entry.trace_key) : "";
        if (!traceKey) {
          return true;
        }
        const defaults = groupDefaults(currentGroup);
        const mode = defaults[traceKey];
        if (mode === false || mode === undefined) {
          return false;
        }
        return legendState[traceKey] !== false;
      }

      function applyClusterFilterState() {
        clampClusterFilterRangeForParameter(activeClusterFilterParameterSpec());
        pruneSelectionsToActiveClusterFilter();
        updateSelectionUI();
        updateSkyPanel();
        renderFrame(currentFrameIndex);
      }

      function renderClusterFilterWidget() {
        if (!clusterFilterSpec.enabled || !clusterFilterCanvasEl || widgetModeForKey("cluster_filter") === "hidden") {
          return;
        }

        const parameter = activeClusterFilterParameterSpec();
        if (!parameter) {
          clusterFilterReadoutEl.textContent = "No filterable cluster parameters available.";
          return;
        }
        const rangeState = clampClusterFilterRangeForParameter(parameter);
        clusterFilterParameterKey = String(parameter.key);
        if (clusterFilterParameterEl) {
          clusterFilterParameterEl.value = clusterFilterParameterKey;
        }
        if (clusterFilterRangeMinEl) {
          clusterFilterRangeMinEl.value = String(clusterFilterActualValueToSlider(rangeState.min, parameter));
        }
        if (clusterFilterRangeMaxEl) {
          clusterFilterRangeMaxEl.value = String(clusterFilterActualValueToSlider(rangeState.max, parameter));
        }
        if (clusterFilterRangeReadoutMinEl) {
          clusterFilterRangeReadoutMinEl.textContent = formatClusterFilterValue(rangeState.min, parameter);
        }
        if (clusterFilterRangeReadoutMaxEl) {
          clusterFilterRangeReadoutMaxEl.textContent = formatClusterFilterValue(rangeState.max, parameter);
        }

        const finiteEntries = clusterFilterEntries.filter((entry) => Number.isFinite(clusterFilterEntryValue(entry, parameter.key)));
        const rect = clusterFilterCanvasEl.getBoundingClientRect();
        const cssWidth = Math.max(1, Math.round(rect.width));
        const cssHeight = Math.max(1, Math.round(rect.height));
        const dpr = Math.max(window.devicePixelRatio || 1, 1);
        if (clusterFilterCanvasEl.width !== Math.round(cssWidth * dpr) || clusterFilterCanvasEl.height !== Math.round(cssHeight * dpr)) {
          clusterFilterCanvasEl.width = Math.round(cssWidth * dpr);
          clusterFilterCanvasEl.height = Math.round(cssHeight * dpr);
        }

        const ctx = clusterFilterCanvasEl.getContext("2d");
        if (!ctx) {
          return;
        }
        ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
        ctx.clearRect(0, 0, cssWidth, cssHeight);
        ctx.fillStyle = theme.scene_bgcolor || theme.paper_bgcolor || "#000000";
        ctx.fillRect(0, 0, cssWidth, cssHeight);

        const margin = { left: 34, right: 12, top: 10, bottom: 22 };
        const plotWidth = Math.max(40, cssWidth - margin.left - margin.right);
        const plotHeight = Math.max(40, cssHeight - margin.top - margin.bottom);
        const xMin = Number(parameter.min);
        const xMax = Number(parameter.max);
        const axisColor = String(theme.axis_color || "#808080");
        const bins = Math.max(12, Math.min(28, Math.round(plotWidth / 18)));
        const counts = new Array(bins).fill(0);
        const highlightedCounts = new Array(bins).fill(0);
        const denom = Math.max(xMax - xMin, 1e-9);

        finiteEntries.forEach((entry) => {
          const value = clusterFilterEntryValue(entry, parameter.key);
          const frac = clampRange((value - xMin) / denom, 0.0, 0.999999);
          const binIndex = Math.max(0, Math.min(bins - 1, Math.floor(frac * bins)));
          counts[binIndex] += 1;
          if (value >= rangeState.min && value <= rangeState.max) {
            highlightedCounts[binIndex] += 1;
          }
        });

        const yMax = Math.max(...counts, 1);
        const binWidth = plotWidth / bins;
        ctx.strokeStyle = axisColor;
        ctx.lineWidth = 1.0;
        ctx.beginPath();
        ctx.moveTo(margin.left, margin.top);
        ctx.lineTo(margin.left, margin.top + plotHeight);
        ctx.lineTo(margin.left + plotWidth, margin.top + plotHeight);
        ctx.stroke();

        for (let i = 0; i < bins; i += 1) {
          const totalCount = counts[i];
          if (!totalCount) {
            continue;
          }
          const x = margin.left + i * binWidth + 1;
          const fullHeight = (totalCount / yMax) * plotHeight;
          const selectedHeight = (highlightedCounts[i] / yMax) * plotHeight;
          ctx.fillStyle = cssColorWithAlpha(theme.axis_color || "#6f7f8f", 0.22, axisColor);
          ctx.fillRect(x, margin.top + plotHeight - fullHeight, Math.max(1, binWidth - 2), fullHeight);
          if (selectedHeight > 0) {
            ctx.fillStyle = cssColorWithAlpha(theme.text_color || "#ffffff", 0.80, "#ffffff");
            ctx.fillRect(x, margin.top + plotHeight - selectedHeight, Math.max(1, binWidth - 2), selectedHeight);
          }
        }

        ctx.fillStyle = axisColor;
        ctx.font = "10px Menlo, Monaco, Consolas, monospace";
        ctx.textAlign = "left";
        ctx.textBaseline = "top";
        ctx.fillText(formatCompactNumber(xMin), margin.left, margin.top + plotHeight + 6);
        ctx.textAlign = "right";
        ctx.fillText(formatCompactNumber(xMax), margin.left + plotWidth, margin.top + plotHeight + 6);

        const inRangeCount = finiteEntries.filter((entry) => {
          const value = clusterFilterEntryValue(entry, parameter.key);
          return value >= rangeState.min && value <= rangeState.max;
        }).length;
        const visibleInRangeCount = finiteEntries.filter((entry) => {
          if (!clusterFilterEntryVisibleInScene(entry)) {
            return false;
          }
          const value = clusterFilterEntryValue(entry, parameter.key);
          return value >= rangeState.min && value <= rangeState.max;
        }).length;
        clusterFilterReadoutEl.textContent = [
          `Filter: ${String(parameter.label || parameter.key)}`,
          `Range: ${formatClusterFilterValue(rangeState.min, parameter)} to ${formatClusterFilterValue(rangeState.max, parameter)}`,
          `Clusters in range: ${inRangeCount} / ${finiteEntries.length}`,
          `Currently visible: ${visibleInRangeCount}`,
        ].join("\\n");
      }

      function traceStyleStateForKey(traceKey) {
        return traceStyleStateByKey[String(traceKey)] || null;
      }

      function volumeLayerForKey(layerKey) {
        return volumeLayersByKey.get(String(layerKey)) || null;
      }

      function volumeLegendColorForLayer(layer) {
        if (!layer) {
          return theme.text_color || theme.axis_color;
        }
        const state = volumeStateByKey[String(layer.key)] || {};
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
        const unitText = layer.value_unit ? ` ${layer.value_unit}` : "";
        const supportText = volumeSupported
          ? "Rendered at t=0 only as a WebGL2 ray-marched volume."
          : "Volume rendering requires WebGL2 support in the browser.";
        const resampleMethod = String(layer.downsample_method || "resampled");
        return [
          `${layer.name}`,
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

      function setClusterSelections(selections, mode = "click") {
        const nextSelections = uniqueSelections(selections);
        const normalizedMode = String(mode || "click");

        if (normalizedMode === "lasso") {
          currentSelections = nextSelections;
          if (!currentSelections.length && !hasActiveLassoSelectionMask()) {
            currentSelection = null;
          }
          selectedClusterKeys = new Set(
            currentSelections
              .map((selection) => normalizedSelectionKeyFor(selection))
              .filter(Boolean)
          );
          if (currentSelection) {
            const focusKey = normalizedSelectionKeyFor(currentSelection);
            if (!focusKey || !selectedClusterKeys.has(focusKey)) {
              currentSelection = null;
            }
          }
        } else {
          const nextFocus = nextSelections.length ? nextSelections[0] : null;
          if (currentSelections.length) {
            const focusKey = normalizedSelectionKeyFor(nextFocus);
            if (!nextFocus || !focusKey || !selectedClusterKeys.has(focusKey)) {
              return;
            }
          }
          currentSelection = nextFocus;
        }

        currentSelectionMode = currentSelection ? "click" : ((currentSelections.length || hasActiveLassoSelectionMask()) ? "lasso" : "none");
        clearCrossHoverState();
        updateSelectionUI();
        updateSkyPanel();
        renderFrame(currentFrameIndex);
      }

      function clearClusterSelections() {
        currentSelections = [];
        currentSelection = null;
        disposeLassoSelectionMask(currentLassoSelectionMask);
        currentLassoSelectionMask = null;
        currentSelectionMode = "none";
        selectedClusterKeys = new Set();
        clearCrossHoverState();
        updateSelectionUI();
        updateSkyPanel();
        renderFrame(currentFrameIndex);
      }

      function applySkyPanelMode() {
        applyWidgetMode("sky");
        skyFullButtonEl.textContent = skyPanelMode === "fullscreen" ? "Window" : "Full";
      }

      function applyAgeKdePanelMode() {
        applyWidgetMode("age_kde");
        ageKdeFullButtonEl.textContent = ageKdePanelMode === "fullscreen" ? "Window" : "Full";
      }

      function applyClusterFilterPanelMode() {
        applyWidgetMode("cluster_filter");
        clusterFilterFullButtonEl.textContent = clusterFilterPanelMode === "fullscreen" ? "Window" : "Full";
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
        return volumeLayersByKey.get(String(activeVolumeKey)) || null;
      }

      function selectedVolumeState() {
        const layer = selectedVolumeLayer();
        if (!layer) {
          return null;
        }
        return volumeStateByKey[String(layer.key)] || null;
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

      const VOLUME_VERTEX_SHADER = `
        uniform vec4 rotation;
        uniform vec4 translation;

        varying vec3 localPosition;
        varying vec3 transformedCameraPosition;
        varying vec3 transformedWorldPosition;

        vec3 rotate_vertex_position(vec3 pos, vec3 t, vec4 q) {
          vec3 p = pos.xyz - t.xyz;
          return p.xyz + 2.0 * cross(cross(p.xyz, q.xyz) + q.w * p.xyz, q.xyz) + t.xyz;
        }

        void main() {
          vec3 transformed = position;
          localPosition = position;
          vec4 worldPosition = modelMatrix * vec4(transformed, 1.0);
          transformedCameraPosition = rotate_vertex_position(cameraPosition.xyz, translation.xyz, rotation);
          transformedWorldPosition = rotate_vertex_position(worldPosition.xyz, translation.xyz, rotation);
          gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0);
        }
      `;

      const VOLUME_FRAGMENT_SHADER = `
        #include <common>
        #include <lights_pars_begin>

        precision highp sampler3D;

        uniform sampler3D volumeTexture;
        uniform sampler2D colormap;
        uniform sampler2D jitterTexture;
        uniform sampler2D selectionMaskTexture;
        uniform float low;
        uniform float high;
        uniform float opacity;
        uniform float samples;
        uniform float alpha_coef;
        uniform float gradient_step;
        uniform float stretch_mode;
        uniform vec4 scale;
        uniform vec4 translation;
        uniform bool useSelectionPolygon;
        uniform mat4 selectionViewProjectionMatrix;
        uniform float selectionDimOutside;

        varying vec3 localPosition;
        varying vec3 transformedCameraPosition;
        varying vec3 transformedWorldPosition;

        float inv_range;

        struct Ray {
          vec3 origin;
          vec3 direction;
          vec3 inv_direction;
          int sign[3];
        };

        vec3 aabb[2] = vec3[2](
          vec3(-0.5, -0.5, -0.5),
          vec3(0.5, 0.5, 0.5)
        );

        Ray makeRay(vec3 origin, vec3 direction) {
          vec3 inv_direction = vec3(1.0) / direction;
          return Ray(
            origin,
            direction,
            inv_direction,
            int[3](
              ((inv_direction.x < 0.0) ? 1 : 0),
              ((inv_direction.y < 0.0) ? 1 : 0),
              ((inv_direction.z < 0.0) ? 1 : 0)
            )
          );
        }

        void intersect(in Ray ray, in vec3 bounds[2], out float tmin, out float tmax) {
          float tymin;
          float tymax;
          float tzmin;
          float tzmax;
          tmin = (bounds[ray.sign[0]].x - ray.origin.x) * ray.inv_direction.x;
          tmax = (bounds[1 - ray.sign[0]].x - ray.origin.x) * ray.inv_direction.x;
          tymin = (bounds[ray.sign[1]].y - ray.origin.y) * ray.inv_direction.y;
          tymax = (bounds[1 - ray.sign[1]].y - ray.origin.y) * ray.inv_direction.y;
          tzmin = (bounds[ray.sign[2]].z - ray.origin.z) * ray.inv_direction.z;
          tzmax = (bounds[1 - ray.sign[2]].z - ray.origin.z) * ray.inv_direction.z;
          tmin = max(max(tmin, tymin), tzmin);
          tmax = min(min(tmax, tymax), tzmax);
        }

        float sampleVolume(vec3 pos) {
          return texture(volumeTexture, clamp(pos, vec3(0.0), vec3(1.0))).x;
        }

        float asinhStretch(float value) {
          return log(value + sqrt(value * value + 1.0));
        }

        float applyStretch(float value) {
          float clampedValue = clamp(value, 0.0, 1.0);
          if (stretch_mode < 0.5) {
            return clampedValue;
          }
          if (stretch_mode < 1.5) {
            float strength = 999.0;
            return log(1.0 + strength * clampedValue) / log(1.0 + strength);
          }
          float asinhStrength = 10.0;
          return asinhStretch(asinhStrength * clampedValue) / asinhStretch(asinhStrength);
        }

        bool sampleInsideSelectionPolygon(vec3 worldPos) {
          if (!useSelectionPolygon) {
            return true;
          }
          vec4 clip = selectionViewProjectionMatrix * vec4(worldPos, 1.0);
          if (clip.w <= 0.0) {
            return false;
          }
          vec3 ndc = clip.xyz / clip.w;
          if (ndc.z < -1.0 || ndc.z > 1.0) {
            return false;
          }
          vec2 uv = vec2(ndc.x * 0.5 + 0.5, 0.5 - ndc.y * 0.5);
          if (uv.x <= 0.0 || uv.x >= 1.0 || uv.y <= 0.0 || uv.y >= 1.0) {
            return false;
          }
          return texture2D(selectionMaskTexture, uv).r > 0.5;
        }

        vec3 worldGetNormal(in float px, in vec3 pos) {
          return normalize(vec3(
            px - sampleVolume(pos + vec3(gradient_step, 0.0, 0.0)),
            px - sampleVolume(pos + vec3(0.0, gradient_step, 0.0)),
            px - sampleVolume(pos + vec3(0.0, 0.0, gradient_step))
          ));
        }

        void main() {
          float jitter = texture2D(jitterTexture, gl_FragCoord.xy / 64.0).r;
          float tmin = 0.0;
          float tmax = 0.0;
          float px = 0.0;
          vec4 pxColor = vec4(0.0);
          vec4 value = vec4(0.0);
          vec3 direction = normalize(transformedWorldPosition - transformedCameraPosition);

          inv_range = 1.0 / max(high - low, 1e-6);
          aabb[0] = aabb[0] * scale.xyz + translation.xyz;
          aabb[1] = aabb[1] * scale.xyz + translation.xyz;
          intersect(makeRay(transformedCameraPosition, direction), aabb, tmin, tmax);

          if (tmax <= max(0.0, tmin)) {
            discard;
          }

          vec3 textcoord_end = localPosition + vec3(0.5);
          vec3 textcoord_start = textcoord_end - (tmax - max(0.0, tmin)) * direction / scale.xyz;
          vec3 textcoord_delta = textcoord_end - textcoord_start;
          int sampleCount = min(int(length(textcoord_delta) * samples), int(samples * 1.8));
          if (sampleCount <= 0) {
            discard;
          }

          textcoord_delta = textcoord_delta / float(sampleCount);
          textcoord_start = textcoord_start - textcoord_delta * (0.01 + 0.98 * jitter);
          vec3 textcoord = textcoord_start - textcoord_delta;
          float step = length(textcoord_delta);

          for (int count = 0; count < 2048; count++) {
            if (count >= sampleCount) {
              break;
            }

            textcoord += textcoord_delta;
            px = texture(volumeTexture, textcoord).x;
            float scaled_px = (px - low) * inv_range;

            if (scaled_px > 0.0) {
              scaled_px = applyStretch(min(scaled_px, 0.999));
              pxColor = texture(colormap, vec2(scaled_px, 0.5));
              vec3 worldPos = (textcoord - vec3(0.5)) * scale.xyz + translation.xyz;
              bool insideSelection = sampleInsideSelectionPolygon(worldPos);
              if (!insideSelection && selectionDimOutside <= 0.001) {
                continue;
              }
              pxColor.a = 1.0 - pow(1.0 - clamp(pxColor.a * opacity, 0.0, 0.999), step * alpha_coef);
              if (!insideSelection) {
                pxColor.a *= selectionDimOutside;
              }
              pxColor.a *= (1.0 - value.a);
              pxColor.rgb *= pxColor.a;

              #if NUM_DIR_LIGHTS > 0
              if (pxColor.a > 0.0) {
                vec4 addedLights = vec4(ambientLightColor / PI, 1.0);
                vec3 specularColor = vec3(0.0);
                vec3 normal = worldGetNormal(px, textcoord);
                vec3 lightDirection;
                float lightingIntensity;
                vec3 lightReflect;
                float specularFactor;

                #pragma unroll_loop_start
                for (int i = 0; i < NUM_DIR_LIGHTS; i++) {
                  lightDirection = directionalLights[i].direction;
                  lightingIntensity = clamp(dot(lightDirection, normal), 0.0, 1.0);
                  addedLights.rgb += directionalLights[i].color / PI * (0.2 + 0.8 * lightingIntensity);
                  lightReflect = normalize(reflect(lightDirection, normal));
                  specularFactor = dot(direction, lightReflect);
                  if (specularFactor > 0.0) {
                    specularColor += 0.002 * scaled_px * (1.0 / max(step, 1e-6))
                      * directionalLights[i].color / PI * pow(specularFactor, 250.0) * pxColor.a;
                  }
                }
                #pragma unroll_loop_end

                pxColor.rgb = pxColor.rgb * addedLights.xyz + specularColor;
              }
              #endif

              value += pxColor;
              if (value.a >= 0.99) {
                value.a = 1.0;
                break;
              }
            }
          }

          if (value.a <= 0.001) {
            discard;
          }
          gl_FragColor = value;
        }
      `;

      function createVolumeRuntime(layer) {
        if (!volumeSupported) {
          return null;
        }
        const state = volumeStateByKey[String(layer.key)];
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
          depthTest: false,
          depthWrite: false,
          lights: true,
        });

        const geometry = new THREE.BoxBufferGeometry(1, 1, 1);
        const mesh = new THREE.Mesh(geometry, material);
        mesh.position.set(centerX, centerY, centerZ);
        mesh.scale.set(sizeX, sizeY, sizeZ);
        mesh.visible = state.visible !== false;
        mesh.renderOrder = -30;
        mesh.frustumCulled = false;
        const runtime = { mesh, material, geometry, layer };
        return runtime;
      }

      function applyVolumeStateToRuntime(layer, runtime) {
        if (!runtime || !runtime.material || !runtime.mesh) {
          return;
        }
        const state = volumeStateByKey[String(layer.key)] || {};
        const option = volumeColormapOptionFor(layer, state.colormap);
        const windowState = normalizedVolumeWindowFor(layer, state);
        runtime.mesh.visible = state.visible !== false;
        runtime.material.uniforms.low.value = Number(windowState.low);
        runtime.material.uniforms.high.value = Number(windowState.high);
        runtime.material.uniforms.opacity.value = Number(state.opacity);
        runtime.material.uniforms.samples.value = Number(state.steps);
        runtime.material.uniforms.alpha_coef.value = Number(state.alphaCoef);
        runtime.material.uniforms.gradient_step.value = Number(state.gradientStep);
        runtime.material.uniforms.stretch_mode.value = volumeStretchModeValue(state.stretch);
        applyLassoSelectionMaskUniforms(runtime.material.uniforms, activeVolumeLassoSelectionMask());
        if (option) {
          runtime.material.uniforms.colormap.value = volumeColorTextureFor(option);
        }
      }

      function renderVolumeControls() {
        const enabled = volumeLayers.length > 0;
        volumePanelEl.dataset.enabled = enabled ? "true" : "false";
        if (!enabled) {
          return;
        }

        if (!activeVolumeKey || !volumeLayersByKey.has(String(activeVolumeKey))) {
          activeVolumeKey = String(volumeLayers[0].key);
        }

        const layer = selectedVolumeLayer();
        const state = selectedVolumeState();
        if (!layer || !state) {
          return;
        }
        clampVolumeStateForLayer(layer, state);

        volumeSelectEl.innerHTML = "";
        volumeLayers.forEach((volumeLayer) => {
          const optionEl = document.createElement("option");
          optionEl.value = String(volumeLayer.key);
          optionEl.textContent = volumeLayer.name || String(volumeLayer.key);
          if (String(volumeLayer.key) === String(activeVolumeKey)) {
            optionEl.selected = true;
          }
          volumeSelectEl.appendChild(optionEl);
        });
        volumeSelectEl.disabled = !volumeSupported || volumeLayers.length <= 1;

        volumeVisibleEl.checked = state.visible !== false;
        volumeColormapEl.innerHTML = "";
        ((layer.colormap_options || [])).forEach((option) => {
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

        const unitText = layer.value_unit ? ` ${layer.value_unit}` : "";
        const supportText = volumeSupported
          ? "Rendered at t=0 only as a WebGL2 ray-marched volume."
          : "Volume rendering requires WebGL2 support in the browser.";
        const resampleMethod = String(layer.downsample_method || "resampled");
        volumeSummaryEl.textContent = [
          `${layer.name}`,
          `${Number(layer.shape.x)} x ${Number(layer.shape.y)} x ${Number(layer.shape.z)} voxels`,
          `Source: ${Number(layer.source_shape.x)} x ${Number(layer.source_shape.y)} x ${Number(layer.source_shape.z)} | ${resampleMethod} scale ${formatVolumeNumber(Number(layer.downsample_step.x))} x ${formatVolumeNumber(Number(layer.downsample_step.y))} x ${formatVolumeNumber(Number(layer.downsample_step.z))}`,
          `Samples: ${Math.round(Number(state.steps))} | alpha_coef: ${Math.round(Number(state.alphaCoef))} | stretch: ${normalizeVolumeStretch(state.stretch)}`,
          `Data range: ${formatVolumeNumber((layer.data_range || [0])[0])} to ${formatVolumeNumber((layer.data_range || [0, 1])[1])}${unitText}`,
          supportText,
        ].join("\\n");
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

      function makeLineObject(trace, materialBucket) {
        if (!trace.segments || !trace.segments.length) {
          return null;
        }
        const traceState = traceStyleStateForKey(trace.key);
        const lineColor = (traceState && traceState.color) || (trace.line || {}).color || "#ffffff";
        const lineOpacity = traceState ? clamp01(traceState.opacity) : (trace.opacity ?? 1.0);
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
          if (selectedClusterKeys.size && pointKey) {
            opacityMultiplier = selectedClusterKeys.has(pointKey) ? 1.0 : 0.16;
          }
          const baseOpacity = Number(pointState.opacity ?? point.opacity ?? 1.0);
          const effectiveOpacity = Math.min(1.0, Math.max(0.0, baseOpacity * traceOpacityMultiplier * opacityMultiplier * globalPointOpacityScale));
          if (effectiveOpacity <= 0.001) {
            return;
          }
          const sprite = new THREE.Sprite(markerMaterialFor(point.symbol, traceColor || point.color, effectiveOpacity));
          const scaleFloor = pointScale * 0.5 * Math.max(globalPointSizeScale, 0.05);
          const scale = Math.max(pointState.size * sizeScaleFactor * globalPointSizeScale * pointScale, scaleFloor);
          const selectionKey = normalizedSelectionKeyFor(point.selection);
          sprite.position.set(point.x, point.y, point.z);
          sprite.scale.set(scale, scale, scale);
          sprite.userData = {
            hovertext: point.hovertext || trace.name || "",
            selection: point.selection || null,
            selectionKey,
            baseScale: scale,
          };
          group.add(sprite);
          hoverTargets.push(sprite);
          if (selectionKey) {
            if (!selectionSpriteEntriesByKey.has(selectionKey)) {
              selectionSpriteEntriesByKey.set(selectionKey, []);
            }
            selectionSpriteEntriesByKey.get(selectionKey).push({
              sprite,
              baseScale: scale,
            });
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
        trace.labels.forEach((label) => {
          if (!label.text) {
            return;
          }
          const sprite = makeTextSprite(label.text, {
            color: (traceState && traceState.color) || label.color || theme.axis_color,
            size: label.size ?? 12,
            family: label.family ?? "Helvetica",
            screenStable: label.screen_stable === true,
            screenPixels: label.screen_px,
          });
          sprite.position.set(label.x, label.y, label.z);
          if (label.screen_stable === true) {
            screenStableTextSprites.push(sprite);
          }
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
          if (legendState[String(layer.key)] === false) {
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

        if (state.hasPoints) {
          const sizeInput = document.createElement("input");
          sizeInput.type = "range";
          sizeInput.min = "0.25";
          sizeInput.max = "4";
          sizeInput.step = "0.05";
          sizeInput.value = String(Math.max(Number(state.sizeScale), 0.25));
          const sizeField = createLegendField(`Point size (${Number(state.sizeScale).toFixed(2)}x)`, sizeInput);
          controls.appendChild(sizeField.field);
          sizeInput.addEventListener("input", () => {
            state.sizeScale = Math.max(Number(sizeInput.value), 0.05);
            sizeField.label.textContent = `Point size (${state.sizeScale.toFixed(2)}x)`;
            renderFrame(currentFrameIndex);
          });

          const summary = document.createElement("div");
          summary.className = "oviz-three-legend-summary";
          summary.textContent = "Point size is applied as a multiplier on the original marker sizes.";
          controls.appendChild(summary);
        }

        return controls;
      }

      function buildVolumeLegendControls(item, toggleButton) {
        const layer = volumeLayerForKey(item.key);
        const state = layer ? volumeStateByKey[String(layer.key)] : null;
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
        visibleText.textContent = "Show at t=0";
        visibleToggle.appendChild(visibleText);
        controls.appendChild(visibleToggle);

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
          clampVolumeStateForLayer(layer, state);
          if (syncInputs) {
            visibleInput.checked = state.visible !== false;
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
          summary.textContent = volumeSummaryTextFor(layer, state);
          toggleButton.style.color = volumeLegendColorForLayer(layer);
        }

        function updateVolumeFromLegend(shouldRerenderLegend = false) {
          activeVolumeKey = String(layer.key);
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

      function renderLegend() {
        legendEl.innerHTML = "";
        const title = document.createElement("div");
        title.className = "oviz-three-legend-title";
        title.textContent = "Click names to toggle. Use arrows for controls.";
        legendEl.appendChild(title);

        const defaults = groupDefaults(currentGroup);
        legendItems.forEach((item) => {
          const mode = defaults[item.key];
          if (mode === false || mode === undefined) {
            return;
          }

          const entry = document.createElement("div");
          entry.className = "oviz-three-legend-entry";
          entry.dataset.open = legendPanelStateByKey[String(item.key)] ? "true" : "false";

          const row = document.createElement("div");
          row.className = "oviz-three-legend-row";
          const button = document.createElement("button");
          button.type = "button";
          button.className = "oviz-three-legend-item";
          button.dataset.active = legendState[item.key] ? "true" : "false";
          button.textContent = item.name;
          button.style.color = legendColorForItem(item);
          button.addEventListener("click", () => {
            legendState[item.key] = !legendState[item.key];
            renderLegend();
            renderFrame(currentFrameIndex);
          });
          if (!volumeLayerForKey(item.key)) {
            button.addEventListener("dblclick", (event) => {
              event.preventDefault();
              event.stopPropagation();
              soloTraceLegendItem(item.key);
            });
          }
          row.appendChild(button);

          const disclosure = document.createElement("button");
          disclosure.type = "button";
          disclosure.className = "oviz-three-legend-disclosure";
          disclosure.dataset.open = legendPanelStateByKey[String(item.key)] ? "true" : "false";
          disclosure.textContent = legendPanelStateByKey[String(item.key)] ? "▾" : "▸";
          disclosure.addEventListener("click", () => {
            const itemKey = String(item.key);
            legendPanelStateByKey[itemKey] = !legendPanelStateByKey[itemKey];
            renderLegend();
          });
          row.appendChild(disclosure);

          entry.appendChild(row);

          const controls = volumeLayerForKey(item.key)
            ? buildVolumeLegendControls(item, button)
            : buildTraceLegendControls(item, button);
          if (controls) {
            entry.appendChild(controls);
          } else {
            disclosure.disabled = true;
          }

          legendEl.appendChild(entry);
        });
      }

      function setToolsDrawerOpen(isOpen) {
        const nextOpen = Boolean(isOpen);
        toolsShellEl.dataset.open = nextOpen ? "true" : "false";
        toolsToggleEl.textContent = nextOpen ? "Tools ◂" : "Tools ▸";
      }

      function setControlsDrawerOpen(isOpen) {
        const nextOpen = Boolean(isOpen);
        controlsShellEl.dataset.open = nextOpen ? "true" : "false";
        controlsToggleEl.textContent = nextOpen ? "Controls ◂" : "Controls ▸";
      }

      function renderFrame(index) {
        currentFrameIndex = Math.max(0, Math.min(index, frameSpecs.length - 1));
        sliderEl.value = String(currentFrameIndex);
        const frame = frameSpecs[currentFrameIndex];
        timeLabelEl.textContent = `Time (Myr): ${frame.name}`;
        tooltipEl.style.display = "none";
        hoverTargets.length = 0;
        selectionSpriteEntriesByKey.clear();
        screenStableTextSprites.length = 0;
        clearGroup(plotGroup);
        plotGroup.position.set(0.0, 0.0, 0.0);
        frameLineMaterials.length = 0;
        volumeRuntimeByKey.clear();

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

        if (currentSelectionMode === "click" && currentSelection && approximatelyZero(Number(frame.time))) {
          const footprint = buildSelectionFootprint(currentSelection, frameLineMaterials);
          if (footprint) {
            plotGroup.add(footprint);
          }
        }

        (frame.decorations || []).forEach((decoration) => {
          addDecoration(plotGroup, decoration);
        });

        const focusOffset = focusTrackingOffsetForFrame(frame);
        if (focusOffset) {
          plotGroup.position.copy(focusOffset.multiplyScalar(-1.0));
        }

        applySceneHoverState();
        resize();
        renderAgeKdeWidget();
        renderClusterFilterWidget();
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
        updateScaleBar();
        renderAgeKdeWidget();
        renderClusterFilterWidget();
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

      function applyWidgetMode(widgetKey) {
        const panelEl = widgetPanelForKey(widgetKey);
        if (!panelEl) {
          return;
        }
        const mode = widgetEnabled(widgetKey) ? widgetModeForKey(widgetKey) : "hidden";
        panelEl.dataset.mode = mode;
      }

      function setWidgetMode(widgetKey, mode) {
        if (!widgetEnabled(widgetKey)) {
          return;
        }
        const panelEl = widgetPanelForKey(widgetKey);
        if (!panelEl) {
          return;
        }
        const currentMode = widgetModeForKey(widgetKey);
        const nextMode = ["normal", "fullscreen", "hidden"].includes(mode) ? mode : "normal";
        if (nextMode === "fullscreen" && currentMode === "normal") {
          storeWidgetRect(widgetKey);
        }
        setWidgetModeValue(widgetKey, nextMode);
        if (widgetKey === "sky") {
          applySkyPanelMode();
        } else if (widgetKey === "age_kde") {
          applyAgeKdePanelMode();
        } else if (widgetKey === "cluster_filter") {
          applyClusterFilterPanelMode();
        } else {
          applyWidgetMode(widgetKey);
        }
        if (nextMode === "normal") {
          restoreWidgetRect(widgetKey);
          raiseWidget(widgetKey);
        } else if (nextMode === "fullscreen") {
          raiseWidget(widgetKey);
          panelEl.style.left = "0px";
          panelEl.style.top = "0px";
          panelEl.style.right = "0px";
          panelEl.style.bottom = "0px";
          panelEl.style.width = "auto";
          panelEl.style.height = "auto";
        }
        resize();
        if (widgetKey === "sky" && nextMode !== "hidden") {
          updateSkyPanel();
        }
        renderAgeKdeWidget();
        renderClusterFilterWidget();
      }

      function renderWidgetMenu() {
        widgetSelectEl.innerHTML = "";
        const placeholder = document.createElement("option");
        placeholder.value = "";
        placeholder.textContent = "Widgets";
        placeholder.selected = true;
        widgetSelectEl.appendChild(placeholder);

        const items = [];
        if (skySpec.enabled) {
          items.push({ key: "sky", label: "Sky View" });
        }
        if (ageKdeSpec.enabled) {
          items.push({ key: "age_kde", label: "Age KDE" });
        }
        if (clusterFilterSpec.enabled) {
          items.push({ key: "cluster_filter", label: "Cluster Filter" });
        }

        items.forEach((item) => {
          const option = document.createElement("option");
          option.value = item.key;
          option.textContent = item.label;
          widgetSelectEl.appendChild(option);
        });
        widgetSelectEl.style.display = items.length ? "block" : "none";
        widgetSelectEl.value = "";
      }

      function currentFrameAllowsSelection() {
        const frame = currentFrame();
        return Boolean(frame);
      }

      function canvasPointFromEvent(event) {
        const rect = canvas.getBoundingClientRect();
        return {
          x: event.clientX - rect.left,
          y: event.clientY - rect.top,
        };
      }

      function updateLassoOverlay() {
        if (!lassoState || !Array.isArray(lassoState.points) || !lassoState.points.length) {
          lassoOverlayEl.dataset.active = "false";
          lassoPolylineEl.setAttribute("points", "");
          return;
        }
        const points = lassoState.points.slice();
        if (points.length > 2) {
          points.push(points[0]);
        }
        lassoOverlayEl.dataset.active = "true";
        lassoPolylineEl.setAttribute("points", points.map((point) => `${point.x},${point.y}`).join(" "));
      }

      function pointInPolygon(point, polygon) {
        let inside = false;
        for (let i = 0, j = polygon.length - 1; i < polygon.length; j = i, i += 1) {
          const xi = polygon[i].x;
          const yi = polygon[i].y;
          const xj = polygon[j].x;
          const yj = polygon[j].y;
          const denom = yj - yi;
          if (Math.abs(denom) < 1e-9) {
            continue;
          }
          const intersect = ((yi > point.y) !== (yj > point.y))
            && (point.x < ((xj - xi) * (point.y - yi)) / denom + xi);
          if (intersect) {
            inside = !inside;
          }
        }
        return inside;
      }

      function spriteScreenPoint(sprite) {
        const projected = new THREE.Vector3();
        sprite.getWorldPosition(projected);
        projected.project(camera);
        if (
          !Number.isFinite(projected.x)
          || !Number.isFinite(projected.y)
          || !Number.isFinite(projected.z)
          || projected.z < -1.0
          || projected.z > 1.0
        ) {
          return null;
        }
        return {
          x: (projected.x * 0.5 + 0.5) * canvas.clientWidth,
          y: (-projected.y * 0.5 + 0.5) * canvas.clientHeight,
        };
      }

      function startLassoSelection(event) {
        if (!currentFrameAllowsSelection() || widgetPointerState || event.button !== 0) {
          return false;
        }
        if (!(event.shiftKey || lassoArmed)) {
          return false;
        }
        lassoState = {
          pointerId: event.pointerId,
          points: [canvasPointFromEvent(event)],
          moved: false,
        };
        setLocalHoveredClusterKey("");
        controls.enabled = false;
        document.body.style.userSelect = "none";
        if (typeof canvas.setPointerCapture === "function" && event.pointerId !== undefined) {
          try {
            canvas.setPointerCapture(event.pointerId);
          } catch (_err) {
          }
        }
        updateLassoOverlay();
        tooltipEl.style.display = "none";
        event.preventDefault();
        return true;
      }

      function onLassoPointerMove(event) {
        if (!lassoState) {
          return;
        }
        const point = canvasPointFromEvent(event);
        const lastPoint = lassoState.points[lassoState.points.length - 1];
        const dx = point.x - lastPoint.x;
        const dy = point.y - lastPoint.y;
        if ((dx * dx + dy * dy) < 4.0) {
          return;
        }
        lassoState.points.push(point);
        lassoState.moved = true;
        updateLassoOverlay();
        tooltipEl.style.display = "none";
        event.preventDefault();
      }

      function finishLassoSelection(event) {
        if (!lassoState) {
          return;
        }
        if (typeof canvas.releasePointerCapture === "function" && lassoState.pointerId !== undefined) {
          try {
            canvas.releasePointerCapture(lassoState.pointerId);
          } catch (_err) {
          }
        }
        controls.enabled = true;
        document.body.style.userSelect = "";
        const polygon = Array.isArray(lassoState.points) ? lassoState.points.slice() : [];
        const shouldSuppressClick = Boolean(lassoState.moved || lassoArmed);
        lassoState = null;
        updateLassoOverlay();
        if (shouldSuppressClick) {
          suppressNextCanvasClick = true;
        }
        if (polygon.length < 3) {
          return;
        }
        disposeLassoSelectionMask(currentLassoSelectionMask);
        currentLassoSelectionMask = captureLassoSelectionMask(polygon);
        const selected = [];
        hoverTargets.forEach((sprite) => {
          const selection = sprite && sprite.userData ? sprite.userData.selection : null;
          if (!selection) {
            return;
          }
          const screenPoint = spriteScreenPoint(sprite);
          if (screenPoint && pointInPolygon(screenPoint, polygon)) {
            selected.push(selection);
          }
        });
        setClusterSelections(selected, "lasso");
        if (event) {
          event.preventDefault();
        }
      }

      function pickSprite(event) {
        const rect = canvas.getBoundingClientRect();
        pointer.x = ((event.clientX - rect.left) / rect.width) * 2.0 - 1.0;
        pointer.y = -((event.clientY - rect.top) / rect.height) * 2.0 + 1.0;
        raycaster.setFromCamera(pointer, camera);
        const hits = raycaster.intersectObjects(hoverTargets, false);
        return hits.length ? hits[0].object : null;
      }

      function pointerRayFromEvent(event) {
        const rect = canvas.getBoundingClientRect();
        pointer.x = ((event.clientX - rect.left) / rect.width) * 2.0 - 1.0;
        pointer.y = -((event.clientY - rect.top) / rect.height) * 2.0 + 1.0;
        raycaster.setFromCamera(pointer, camera);
        return raycaster.ray;
      }

      function doubleClickTargetFromEvent(event) {
        pointerRayFromEvent(event);
        const spriteHits = raycaster.intersectObjects(hoverTargets, false);
        if (spriteHits.length && spriteHits[0].object) {
          const worldPoint = new THREE.Vector3();
          const hitObject = spriteHits[0].object;
          hitObject.getWorldPosition(worldPoint);
          return {
            worldPoint,
            selection: hitObject.userData ? (hitObject.userData.selection || null) : null,
            selectionKey: hitObject.userData ? normalizeMemberKey(hitObject.userData.selectionKey || "") : "",
          };
        }

        const plotHits = raycaster.intersectObjects(plotGroup.children, true);
        for (const hit of plotHits) {
          if (!hit || !hit.object) {
            continue;
          }
          if (hoverTargets.includes(hit.object)) {
            continue;
          }
          if (hit.point && Number.isFinite(hit.point.x) && Number.isFinite(hit.point.y) && Number.isFinite(hit.point.z)) {
            return {
              worldPoint: hit.point.clone(),
              selection: null,
              selectionKey: "",
            };
          }
        }

        const cameraDirection = new THREE.Vector3();
        camera.getWorldDirection(cameraDirection);
        if (cameraDirection.lengthSq() > 1e-12) {
          const focusPlane = new THREE.Plane().setFromNormalAndCoplanarPoint(
            cameraDirection.normalize(),
            controls.target.clone()
          );
          const planePoint = new THREE.Vector3();
          if (raycaster.ray.intersectPlane(focusPlane, planePoint)) {
            return {
              worldPoint: planePoint,
              selection: null,
              selectionKey: "",
            };
          }
        }
        return null;
      }

      function recenterCameraTarget(worldPoint) {
        if (!worldPoint || !Number.isFinite(worldPoint.x) || !Number.isFinite(worldPoint.y) || !Number.isFinite(worldPoint.z)) {
          return;
        }
        cameraViewMode = "free";
        earthViewFocusDistance = null;
        const delta = new THREE.Vector3().subVectors(worldPoint, controls.target);
        if (delta.lengthSq() <= 1e-18) {
          return;
        }
        controls.target.add(delta);
        camera.position.add(delta);
        applyCameraViewMode();
        controls.update();
        updateScaleBar();
      }

      function onCanvasClick(event) {
        if (suppressNextCanvasClick) {
          suppressNextCanvasClick = false;
          return;
        }
        if (widgetPointerState || !clickSelectionEnabled) {
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
        setClusterSelections([selection], "click");
      }

      function onCanvasDoubleClick(event) {
        if (widgetPointerState || lassoState || event.button !== 0) {
          return;
        }
        const target = doubleClickTargetFromEvent(event);
        if (!target || !target.worldPoint) {
          return;
        }
        if (target.selectionKey) {
          focusSelectionKey = target.selectionKey;
          renderFrame(currentFrameIndex);
          recenterCameraTarget(new THREE.Vector3(0.0, 0.0, 0.0));
        } else {
          focusSelectionKey = "";
          recenterCameraTarget(target.worldPoint);
        }
        event.preventDefault();
      }

      function onWidgetPointerStart(event) {
        const panelEl = event.target.closest(".oviz-three-widget-panel");
        if (!panelEl) {
          return;
        }
        const widgetKey = String(panelEl.dataset.widgetKey || "");
        if (widgetModeForKey(widgetKey) !== "normal") {
          return;
        }
        const resizeHandle = event.target.closest(".oviz-three-widget-resize");
        const dragHandle = event.target.closest(".oviz-three-widget-drag");
        if (!resizeHandle && !dragHandle) {
          return;
        }
        const rect = panelEl.getBoundingClientRect();
        widgetPointerState = {
          widgetKey,
          panelEl,
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
        panelEl.style.left = `${rect.left}px`;
        panelEl.style.top = `${rect.top}px`;
        panelEl.style.right = "auto";
        panelEl.style.bottom = "auto";
        panelEl.style.width = `${rect.width}px`;
        panelEl.style.height = `${rect.height}px`;
        raiseWidget(widgetKey);
        controls.enabled = false;
        if (widgetPointerState.handle && typeof widgetPointerState.handle.setPointerCapture === "function" && event.pointerId !== undefined) {
          try {
            widgetPointerState.handle.setPointerCapture(event.pointerId);
          } catch (_err) {
          }
        }
        document.body.style.userSelect = "none";
        event.preventDefault();
        event.stopPropagation();
      }

      function onWidgetPointerMove(event) {
        if (!widgetPointerState) {
          return;
        }
        if (widgetPointerState.mode === "drag") {
          const left = widgetPointerState.startLeft + (event.clientX - widgetPointerState.startX);
          const top = widgetPointerState.startTop + (event.clientY - widgetPointerState.startY);
          const next = clampWidgetPosition(left, top, widgetPointerState.startWidth, widgetPointerState.startHeight);
          widgetPointerState.panelEl.style.left = `${next.left}px`;
          widgetPointerState.panelEl.style.top = `${next.top}px`;
        } else {
          const next = resizeWidgetRect(widgetPointerState, event.clientX, event.clientY);
          widgetPointerState.panelEl.style.left = `${next.left}px`;
          widgetPointerState.panelEl.style.top = `${next.top}px`;
          widgetPointerState.panelEl.style.width = `${next.width}px`;
          widgetPointerState.panelEl.style.height = `${next.height}px`;
        }
        resize();
        renderAgeKdeWidget();
        event.preventDefault();
        event.stopPropagation();
      }

      function onWidgetPointerEnd(event) {
        if (!widgetPointerState) {
          return;
        }
        if (widgetPointerState.handle && typeof widgetPointerState.handle.releasePointerCapture === "function" && event.pointerId !== undefined) {
          try {
            widgetPointerState.handle.releasePointerCapture(event.pointerId);
          } catch (_err) {
          }
        }
        controls.enabled = true;
        document.body.style.userSelect = "";
        storeWidgetRect(widgetPointerState.widgetKey);
        widgetPointerState = null;
      }

      function initSkyPanel() {
        if (!skySpec.enabled) {
          applySkyPanelMode();
          return;
        }
        skyFrameEl.addEventListener("load", () => {
          lastSentSkyHoverClusterKey = null;
          postParentHoverToSkyFrame();
        });
        skyReadoutEl.textContent = skyReadoutText([], [], "overview");
        skyFrameEl.srcdoc = buildEmptySkySrcdoc();
        skyHideButtonEl.addEventListener("click", () => setWidgetMode("sky", "hidden"));
        skyFullButtonEl.addEventListener("click", () => {
          setWidgetMode("sky", skyPanelMode === "fullscreen" ? "normal" : "fullscreen");
        });
        applySkyPanelMode();
      }

      function initAgeKdePanel() {
        if (!ageKdeSpec.enabled) {
          applyAgeKdePanelMode();
          return;
        }
        ageKdeHideButtonEl.addEventListener("click", () => setWidgetMode("age_kde", "hidden"));
        ageKdeFullButtonEl.addEventListener("click", () => {
          setWidgetMode("age_kde", ageKdePanelMode === "fullscreen" ? "normal" : "fullscreen");
        });
        if (ageKdeFilterRangeMinEl && ageKdeFilterRangeMaxEl) {
          ageKdeFilterRangeMinEl.addEventListener("input", () => {
            const parameter = ageKdeFilterParameterSpec();
            if (!parameter) {
              return;
            }
            const currentAxisRange = ageKdeAxisFilterRange();
            const minAxisValue = ageKdeSliderValueToAxisValue(ageKdeFilterRangeMinEl.value);
            const maxAxisValue = currentAxisRange ? Number(currentAxisRange.max) : ageKdeSliderValueToAxisValue(ageKdeFilterRangeMaxEl.value);
            setClusterAgeFilterFromKdeAxisRange(Math.min(minAxisValue, maxAxisValue), maxAxisValue);
          });
          ageKdeFilterRangeMaxEl.addEventListener("input", () => {
            const parameter = ageKdeFilterParameterSpec();
            if (!parameter) {
              return;
            }
            const currentAxisRange = ageKdeAxisFilterRange();
            const minAxisValue = currentAxisRange ? Number(currentAxisRange.min) : ageKdeSliderValueToAxisValue(ageKdeFilterRangeMinEl.value);
            const maxAxisValue = ageKdeSliderValueToAxisValue(ageKdeFilterRangeMaxEl.value);
            setClusterAgeFilterFromKdeAxisRange(minAxisValue, Math.max(minAxisValue, maxAxisValue));
          });
        }
        renderAgeKdeWidget();
        applyAgeKdePanelMode();
      }

      function initClusterFilterPanel() {
        if (!clusterFilterSpec.enabled) {
          applyClusterFilterPanelMode();
          return;
        }
        clusterFilterParameterEl.innerHTML = "";
        clusterFilterParameters.forEach((parameter) => {
          const option = document.createElement("option");
          option.value = String(parameter.key || "");
          option.textContent = String(parameter.label || parameter.key || "");
          clusterFilterParameterEl.appendChild(option);
        });
        if (!clusterFilterParameterKey && clusterFilterParameters.length) {
          clusterFilterParameterKey = String(clusterFilterParameters[0].key || "");
        }
        clusterFilterHideButtonEl.addEventListener("click", () => setWidgetMode("cluster_filter", "hidden"));
        clusterFilterFullButtonEl.addEventListener("click", () => {
          setWidgetMode("cluster_filter", clusterFilterPanelMode === "fullscreen" ? "normal" : "fullscreen");
        });
        clusterFilterParameterEl.addEventListener("change", () => {
          clusterFilterParameterKey = String(clusterFilterParameterEl.value || "");
          clampClusterFilterRangeForParameter(activeClusterFilterParameterSpec());
          applyClusterFilterState();
        });
        clusterFilterRangeMinEl.addEventListener("input", () => {
          const parameter = activeClusterFilterParameterSpec();
          if (!parameter) {
            return;
          }
          const rangeState = clampClusterFilterRangeForParameter(parameter);
          rangeState.min = Math.min(clusterFilterSliderValueToActual(clusterFilterRangeMinEl.value, parameter), Number(rangeState.max));
          clusterFilterRangeStateByKey[String(parameter.key)] = rangeState;
          applyClusterFilterState();
        });
        clusterFilterRangeMaxEl.addEventListener("input", () => {
          const parameter = activeClusterFilterParameterSpec();
          if (!parameter) {
            return;
          }
          const rangeState = clampClusterFilterRangeForParameter(parameter);
          rangeState.max = Math.max(clusterFilterSliderValueToActual(clusterFilterRangeMaxEl.value, parameter), Number(rangeState.min));
          clusterFilterRangeStateByKey[String(parameter.key)] = rangeState;
          applyClusterFilterState();
        });
        renderClusterFilterWidget();
        applyClusterFilterPanelMode();
      }

      function updateActiveVolumeRuntime() {
        const layer = selectedVolumeLayer();
        const state = selectedVolumeState();
        if (!layer || !state) {
          renderVolumeControls();
          return;
        }
        clampVolumeStateForLayer(layer, state);
        renderVolumeControls();
        const runtime = volumeRuntimeByKey.get(String(layer.key));
        if (runtime) {
          applyVolumeStateToRuntime(layer, runtime);
          if (skySpec.enabled && !currentSelection) {
            updateSkyPanel();
          }
          return;
        }
        const frame = currentFrame();
        if (frame && approximatelyZero(Number(frame.time))) {
          renderFrame(currentFrameIndex);
          if (skySpec.enabled && !currentSelection) {
            updateSkyPanel();
          }
        }
      }

      function initVolumeControls() {
        renderVolumeControls();
        if (!volumeLayers.length) {
          return;
        }

        volumeSelectEl.addEventListener("change", () => {
          activeVolumeKey = String(volumeSelectEl.value);
          renderVolumeControls();
        });
        volumeVisibleEl.addEventListener("change", () => {
          const state = selectedVolumeState();
          if (!state) {
            return;
          }
          state.visible = Boolean(volumeVisibleEl.checked);
          updateActiveVolumeRuntime();
        });
        volumeColormapEl.addEventListener("change", () => {
          const state = selectedVolumeState();
          if (!state) {
            return;
          }
          state.colormap = String(volumeColormapEl.value);
          updateActiveVolumeRuntime();
        });
        volumeStretchEl.addEventListener("change", () => {
          const state = selectedVolumeState();
          if (!state) {
            return;
          }
          state.stretch = normalizeVolumeStretch(volumeStretchEl.value);
          updateActiveVolumeRuntime();
        });
        volumeOpacityEl.addEventListener("input", () => {
          const state = selectedVolumeState();
          if (!state) {
            return;
          }
          state.opacity = Number(volumeOpacityEl.value);
          updateActiveVolumeRuntime();
        });
        volumeAlphaEl.addEventListener("input", () => {
          const state = selectedVolumeState();
          if (!state) {
            return;
          }
          state.alphaCoef = Number(volumeAlphaEl.value);
          updateActiveVolumeRuntime();
        });
        volumeStepsEl.addEventListener("input", () => {
          const state = selectedVolumeState();
          if (!state) {
            return;
          }
          state.steps = Number(volumeStepsEl.value);
          updateActiveVolumeRuntime();
        });
        volumeVMinEl.addEventListener("change", () => {
          const state = selectedVolumeState();
          if (!state) {
            return;
          }
          state.vmin = Number(volumeVMinEl.value);
          updateActiveVolumeRuntime();
        });
        volumeVMaxEl.addEventListener("change", () => {
          const state = selectedVolumeState();
          if (!state) {
            return;
          }
          state.vmax = Number(volumeVMaxEl.value);
          updateActiveVolumeRuntime();
        });
      }

      function onPointerMove(event) {
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
        setToolsDrawerOpen(toolsShellEl.dataset.open === "true");
        setControlsDrawerOpen(controlsShellEl.dataset.open === "true");
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
        renderWidgetMenu();
        groupSelectEl.addEventListener("change", () => {
          currentGroup = groupSelectEl.value;
          resetLegendState(currentGroup);
          renderLegend();
          renderFrame(currentFrameIndex);
        });
        renderSceneControls();
        widgetSelectEl.addEventListener("change", () => {
          const widgetKey = String(widgetSelectEl.value || "");
          if (widgetKey) {
            setWidgetMode(widgetKey, "normal");
          }
          widgetSelectEl.value = "";
        });
        zenModeButtonEl.addEventListener("click", () => {
          setZenMode(!zenModeEnabled);
          focusViewer();
        });
        resetViewButtonEl.addEventListener("click", () => {
          resetCameraAndSelections();
          focusViewer();
        });
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

        sliderEl.max = String(Math.max(frameSpecs.length - 1, 0));
        sliderEl.addEventListener("input", () => {
          pause();
          renderFrame(Number(sliderEl.value));
        });
        playButtonEl.addEventListener("click", play);
        pauseButtonEl.addEventListener("click", pause);
        toolsToggleEl.addEventListener("click", () => {
          setToolsDrawerOpen(toolsShellEl.dataset.open !== "true");
        });
        controlsToggleEl.addEventListener("click", () => {
          setControlsDrawerOpen(controlsShellEl.dataset.open !== "true");
        });
        keyHelpButtonEl.addEventListener("click", () => {
          const nextOpen = keyHelpEl.dataset.open !== "true";
          setKeyHelpOpen(nextOpen);
          if (!nextOpen) {
            focusViewer();
          }
        });
        keyHelpCloseEl.addEventListener("click", () => {
          setKeyHelpOpen(false);
          focusViewer();
        });
        lassoButtonEl.addEventListener("click", () => {
          lassoArmed = !lassoArmed;
          updateSelectionUI();
        });
        clearSelectionButtonEl.addEventListener("click", clearClusterSelections);
        clickSelectToggleEl.addEventListener("change", () => {
          clickSelectionEnabled = Boolean(clickSelectToggleEl.checked);
          updateSelectionUI();
        });
        volumeLassoToggleEl.addEventListener("change", () => {
          lassoVolumeSelectionEnabled = Boolean(volumeLassoToggleEl.checked);
          updateSelectionUI();
          const frame = currentFrame();
          if (frame && approximatelyZero(Number(frame.time)) && volumeLayers.length) {
            renderFrame(currentFrameIndex);
          }
          if (skySpec.enabled && !currentSelection) {
            updateSkyPanel();
          }
        });
        themeSelectEl.addEventListener("change", () => {
          applyThemePreset(themeSelectEl.value);
        });
        scrollSpeedEl.addEventListener("input", () => {
          globalScrollSpeed = Number(scrollSpeedEl.value);
          applyGlobalControlState();
          renderSceneControls();
        });
        cameraFovEl.addEventListener("input", () => {
          camera.fov = Number(cameraFovEl.value);
          applyGlobalControlState();
          renderSceneControls();
        });
        globalPointSizeEl.addEventListener("input", () => {
          globalPointSizeScale = Number(globalPointSizeEl.value);
          applyGlobalControlState();
          renderSceneControls();
          renderFrame(currentFrameIndex);
        });
        globalPointOpacityEl.addEventListener("input", () => {
          globalPointOpacityScale = Number(globalPointOpacityEl.value);
          applyGlobalControlState();
          renderSceneControls();
          renderFrame(currentFrameIndex);
        });
        focusGroupSelectEl.addEventListener("change", () => {
          focusTraceKey = String(focusGroupSelectEl.value || "");
          focusSelectionKey = "";
          applyGlobalControlState();
          renderSceneControls();
          renderFrame(currentFrameIndex);
        });
        fadeTimeEl.addEventListener("change", () => {
          fadeInTimeMyr = Number(fadeTimeEl.value);
          applyGlobalControlState();
          renderSceneControls();
          renderFrame(currentFrameIndex);
        });
        fadeInOutToggleEl.addEventListener("change", () => {
          fadeInAndOutEnabled = Boolean(fadeInOutToggleEl.checked);
          applyGlobalControlState();
          renderSceneControls();
          renderFrame(currentFrameIndex);
        });
        axesVisibleToggleEl.addEventListener("change", () => {
          axesVisible = Boolean(axesVisibleToggleEl.checked);
          renderSceneControls();
          buildAxes();
        });
        viewFromEarthButtonEl.addEventListener("click", () => {
          viewFromEarth();
        });
        resetCameraButtonEl.addEventListener("click", () => {
          resetCameraView();
          renderSceneControls();
        });
        resetControlsButtonEl.addEventListener("click", () => {
          activeThemeKey = "default";
          globalScrollSpeed = 1.0;
          globalPointSizeScale = 1.0;
          globalPointOpacityScale = 1.0;
          fadeInTimeMyr = Number(animationSpec.fade_in_time_default);
          fadeInAndOutEnabled = Boolean(animationSpec.fade_in_and_out_default);
          focusTraceKey = String(animationSpec.focus_trace_key_default || "");
          axesVisible = Boolean(sceneSpec.show_axes);
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
        widgetDragHandles.forEach((handle) => handle.addEventListener("pointerdown", onWidgetPointerStart));
        widgetResizeEls.forEach((handle) => handle.addEventListener("pointerdown", onWidgetPointerStart));
        setZenMode(zenModeEnabled);
      }

      function onCanvasPointerDown(event) {
        focusViewer();
        startLassoSelection(event);
      }

      function onWindowPointerMove(event) {
        onWidgetPointerMove(event);
      }

      function onWindowPointerEnd(event) {
        onWidgetPointerEnd(event);
        finishLassoSelection(event);
      }

      function animate(timestamp) {
        window.requestAnimationFrame(animate);
        const now = Number(timestamp) || 0.0;
        const deltaSeconds = lastAnimationTimestamp === null
          ? 0.0
          : clampRange((now - lastAnimationTimestamp) / 1000.0, 0.0, 0.05);
        lastAnimationTimestamp = now;
        updateKeyboardMotion(deltaSeconds);
        controls.update();
        updateScreenStableTextSprites();
        updateScaleBar();
        renderAgeKdeWidget();
        renderer.render(scene, camera);
      }

      applyInitialStateSync();
      buildAxes();
      initControls();
      initSkyPanel();
      initAgeKdePanel();
      initClusterFilterPanel();
      initVolumeControls();
      await restoreInitialLassoSelectionMask();
      renderLegend();
      updateSelectionUI();
      renderFrame(currentFrameIndex);
      resize();
      window.setTimeout(() => focusViewer(), 0);
      animate();

      canvas.addEventListener("pointerdown", onCanvasPointerDown);
      canvas.addEventListener("pointerenter", focusViewer);
      canvas.addEventListener("pointermove", onPointerMove);
      canvas.addEventListener("pointerleave", onPointerLeave);
      canvas.addEventListener("click", onCanvasClick);
      canvas.addEventListener("dblclick", onCanvasDoubleClick);
      window.addEventListener("keydown", onKeyDown);
      window.addEventListener("keyup", onKeyUp);
      window.addEventListener("blur", clearPressedKeys);
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
