#!/usr/bin/env python3
"""Build a local, deterministic Actions/States audit copy of the July figure.

This helper never modifies the canonical July HTML.  Its output is intentionally
an untracked local diagnostic artifact and is not part of the upload workflow.
"""

from __future__ import annotations

import argparse
import base64
import copy
import gzip
import json
import re
from pathlib import Path
from typing import Any

from oviz.threejs_figure import ThreeJSFigure


HERE = Path(__file__).resolve().parent
DEFAULT_INPUT = HERE / "main_figure_chronos_july4.html"
DEFAULT_OUTPUT = HERE / "main_figure_chronos_july4_actions_audit.html"

RUNTIME_SNAPSHOT_KEYS = {
    "current_group",
    "current_frame_index",
    "current_frame_value",
    "playback_state",
    "manual_labels",
    "active_manual_label_id",
    "legend_state",
    "trace_style_state",
    "click_selection_enabled",
    "lasso_volume_selection_enabled",
    "lasso_selection_filter_enabled",
    "lasso_armed",
    "current_selection",
    "current_selections",
    "current_selection_mode",
    "selected_cluster_keys",
    "active_volume_key",
    "volume_state_by_key",
    "global_controls",
    "camera",
    "drawers",
    "legend_open",
    "legend_panel_state",
    "legend_panel_user_sized",
    "legend_sections_open",
    "sky_layers",
    "active_sky_layer_key",
    "zen_mode_enabled",
    "widgets",
    "scale_bar_state",
    "selection_box_state",
    "cluster_filter_state",
    "dendrogram_state",
    "lasso_selection_mask",
}


def read_embedded_scene_spec(path: Path) -> dict[str, Any]:
    html = Path(path).read_text(encoding="utf-8")
    marker = re.search(
        r"/\*__SCENE_SPEC_START__\*/const sceneSpec = (.*?);/\*__SCENE_SPEC_END__\*/",
        html,
        re.DOTALL,
    )
    if marker is None:
        raise ValueError(f"Could not find the scene payload marker in {path}")
    expression = marker.group(1).strip()
    if expression.startswith("{"):
        return json.loads(expression)

    payload_match = re.fullmatch(
        r"await\s+inflateOvizGzipBase64SceneSpec\(readOvizSceneSpecPayload\((.*?)\)\)",
        expression,
        re.DOTALL,
    )
    if payload_match is None:
        raise ValueError("Unsupported scene payload expression")
    payload_id = json.loads(payload_match.group(1))
    chunks = re.findall(
        (
            r"<script\b(?=[^>]*type=[\"']application/octet-stream[\"'])"
            r"(?=[^>]*data-oviz-payload-id=[\"']"
            + re.escape(payload_id)
            + r"[\"'])[^>]*data-oviz-payload-index=[\"'](\d+)[\"'][^>]*>(.*?)</script>"
        ),
        html,
        re.DOTALL,
    )
    if not chunks:
        raise ValueError(f"Could not find compressed payload chunks for {payload_id}")
    encoded = "".join(
        re.sub(r"\s+", "", chunk)
        for _, chunk in sorted(chunks, key=lambda item: int(item[0]))
    )
    return json.loads(gzip.decompress(base64.b64decode(encoded)))


def _camera(scene: dict[str, Any], *, azimuth: float, elevation: float, fov: float) -> dict[str, Any]:
    center = scene.get("center") or {"x": 0.0, "y": 0.0, "z": 0.0}
    span = max(float(scene.get("max_span") or 10000.0), 1000.0)
    return {
        "position": {
            "x": float(center.get("x", 0.0)) + azimuth * span,
            "y": float(center.get("y", 0.0)) - 0.82 * span,
            "z": float(center.get("z", 0.0)) + elevation * span,
        },
        "target": {axis: float(center.get(axis, 0.0)) for axis in ("x", "y", "z")},
        "up": {"x": 0.0, "y": 0.0, "z": 1.0},
        "view_offset": {"x": 0.0, "y": 0.0},
        "fov": fov,
    }


def _trace_key(frame: dict[str, Any], predicate) -> str:
    for trace in frame.get("traces", []):
        if predicate(trace):
            return str(trace.get("key") or "")
    return ""


def _snapshot(
    scene: dict[str, Any],
    frame_index: int,
    *,
    camera: dict[str, Any],
    red_only: bool = False,
    sky: bool = False,
    lasso_seed: int | None = None,
) -> dict[str, Any]:
    frames = scene["frames"]
    frame_index = max(0, min(int(frame_index), len(frames) - 1))
    frame = frames[frame_index]
    initial_state = scene.get("initial_state") or {}
    snapshot = {
        key: copy.deepcopy(value)
        for key, value in initial_state.items()
        if key in RUNTIME_SNAPSHOT_KEYS
    }
    snapshot.update(
        {
            "current_group": "Clusters" if "Clusters" in scene.get("group_order", []) else scene.get("default_group", "All"),
            "current_frame_index": frame_index,
            "current_frame_value": float(frame_index),
            "playback_state": {"direction": 0, "interval_ms": int(scene.get("playback_interval_ms") or 160)},
            "camera": {key: value for key, value in camera.items() if key != "fov"},
            "trace_style_state": {},
            "click_selection_enabled": False,
            "lasso_volume_selection_enabled": False,
            "current_selection": None,
            "current_selections": [],
            "current_selection_mode": "none",
            "selected_cluster_keys": [],
            "lasso_selection_mask": None,
            "lasso_selection_filter_enabled": True,
            "zen_mode_enabled": True,
        }
    )
    controls = copy.deepcopy(snapshot.get("global_controls") or {})
    controls.update(
        {
            "camera_fov": float(camera["fov"]),
            "camera_view_mode": "earth" if sky else "free",
            "camera_auto_orbit_enabled": False,
            "point_opacity_scale": 1.0,
            "point_size_scale": 1.0,
        }
    )
    snapshot["global_controls"] = controls

    if lasso_seed is not None:
        selections: list[dict[str, Any]] = []
        seen: set[str] = set()
        for trace in frame.get("traces", []):
            for point in trace.get("points", []):
                selection = point.get("selection")
                if not isinstance(selection, dict):
                    continue
                identity = str(
                    selection.get("cluster_name")
                    or selection.get("trace_name")
                    or ""
                ).strip()
                key = " ".join(identity.lower().replace("_", " ").split())
                if not key or key == "sun" or key in seen:
                    continue
                seen.add(key)
                selections.append(copy.deepcopy(selection))
        if selections:
            offset = max(int(lasso_seed), 0) % len(selections)
            chosen = (selections[offset:] + selections[:offset])[: min(12, len(selections))]
            selected_keys = [
                " ".join(
                    str(item.get("cluster_name") or item.get("trace_name") or "")
                    .lower()
                    .replace("_", " ")
                    .split()
                )
                for item in chosen
            ]
            snapshot.update(
                {
                    "current_selection": copy.deepcopy(chosen[0]),
                    "current_selections": copy.deepcopy(chosen),
                    "current_selection_mode": "lasso",
                    "selected_cluster_keys": selected_keys,
                    "lasso_selection_filter_enabled": True,
                    "lasso_selection_mask": {
                        "data_url": (
                            "data:image/png;base64,"
                            "iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAQAAAC1HAwC"
                            "AAAAC0lEQVR42mP8/x8AAusB9Wl6w7sAAAAASUVORK5CYII="
                        ),
                        "view_projection_matrix": [
                            1.0, 0.0, 0.0, 0.0,
                            0.0, 1.0, 0.0, 0.0,
                            0.0, 0.0, 1.0, 0.0,
                            0.0, 0.0, 0.0, 1.0,
                        ],
                        "polygon_ndc": [
                            {"x": -0.8, "y": -0.8},
                            {"x": 0.8, "y": -0.8},
                            {"x": 0.0, "y": 0.8},
                        ],
                    },
                }
            )

    if red_only:
        young_key = _trace_key(frame, lambda trace: "< 15" in str(trace.get("name") or ""))
        legend = {}
        for trace in frame.get("traces", []):
            key = str(trace.get("key") or "")
            if key and trace.get("showlegend", True):
                legend[key] = key == young_key
        snapshot["legend_state"] = legend
        if young_key:
            snapshot["trace_style_state"] = {
                young_key: {"color": "#ff3038", "opacity": 1.0, "sizeScale": 1.15}
            }
    return snapshot


def build_audit_scene(scene_spec: dict[str, Any]) -> dict[str, Any]:
    scene = copy.deepcopy(scene_spec)
    scene["initial_state"] = copy.deepcopy(scene.get("initial_state") or {})
    scene["initial_state"]["zen_mode_enabled"] = True
    scene["initial_state"]["minimal_mode_enabled"] = True
    frames = scene.get("frames") or []
    if len(frames) < 2:
        raise ValueError("The Actions audit requires at least two timeline frames")
    present = len(frames) - 1
    lookback = max(0, present - 5)
    past_time = max(0, present - 20)
    topology_boundary = max(0, present - 1)
    state_items = [
        {
            "id": "audit-lookback",
            "name": "Audit All Domains",
            "transition": {"duration_ms": 2400, "easing": "easeInOutCubic"},
            "snapshot": _snapshot(
                scene,
                lookback,
                camera=_camera(scene, azimuth=-0.22, elevation=0.38, fov=44.0),
                red_only=True,
                lasso_seed=1,
            ),
        },
        {
            "id": "audit-topology-boundary",
            "name": "Audit t=-1 Boundary",
            "transition": {"duration_ms": 2400, "easing": "linear"},
            "snapshot": _snapshot(
                scene,
                topology_boundary,
                camera=_camera(scene, azimuth=0.18, elevation=0.18, fov=48.0),
                lasso_seed=7,
            ),
        },
        {
            "id": "audit-past-time",
            "name": "Audit Past Time",
            "transition": {"duration_ms": 2400, "easing": "easeInOutCubic"},
            "snapshot": _snapshot(
                scene,
                past_time,
                camera=_camera(scene, azimuth=-0.05, elevation=0.30, fov=46.0),
            ),
        },
        {
            "id": "audit-present",
            "name": "Audit Present Day",
            "transition": {"duration_ms": 2400, "easing": "easeInOutCubic"},
            "snapshot": _snapshot(
                scene,
                present,
                camera=_camera(scene, azimuth=0.42, elevation=0.24, fov=52.0),
                red_only=True,
            ),
        },
    ]
    scene["states"] = {
        "schema_version": 1,
        "project_id": "project-july4-actions-audit",
        "revision": 1,
        "synchronized_revision": 1,
        "default_mode": "present",
        "default_transition": {"duration_ms": 2400, "easing": "easeInOutCubic"},
        "items": state_items,
        "assets": {},
    }
    group = "Clusters" if "Clusters" in scene.get("group_order", []) else scene.get("default_group", "All")
    scene["actions"] = {
        "enabled": True,
        "items": [
            {
                "key": "audit-state-present",
                "label": "Audit: State Present",
                "description": "State-backed transition to the present-day topology.",
                "steps": [{"type": "state", "start": "after_previous", "delay_ms": 0, "state": "audit-present"}],
            },
            {
                "key": "audit-state-lookback",
                "label": "Audit: State Lookback",
                "description": "State-backed time, camera, and appearance transition.",
                "steps": [{"type": "state", "start": "after_previous", "delay_ms": 0, "state": "audit-lookback"}],
            },
            {
                "key": "audit-state-past-time",
                "label": "Audit: State Past Time",
                "description": "Unfiltered state-backed transition to t = -20 Myr.",
                "steps": [{"type": "state", "start": "after_previous", "delay_ms": 0, "state": "audit-past-time"}],
            },
            {
                "key": "audit-sequential",
                "label": "Audit: Sequential",
                "description": "Legend, camera, then fractional time using actual completion ordering.",
                "steps": [
                    {
                        "type": "legend_group",
                        "start": "after_previous",
                        "delay_ms": 0,
                        "group": group,
                        "duration_ms": 800,
                        "easing": "easeInOutCubic",
                    },
                    {
                        "type": "camera",
                        "start": "after_previous",
                        "delay_ms": 40,
                        "target": {"kind": "group", "name": group},
                        "duration_ms": 800,
                        "easing": "easeInOutCubic",
                        "fit_padding": 1.15,
                        "distance_scale": 0.9,
                        "orbit": {"enabled": False, "speed_multiplier": 1.0, "persist": True, "direction": 1.0},
                    },
                    {
                        "type": "time",
                        "start": "after_previous",
                        "delay_ms": 40,
                        "direction": "forward",
                        "interval_ms": 160,
                        "stop_at_time_myr": 0.0,
                    },
                ],
            },
            {
                "key": "audit-concurrent",
                "label": "Audit: Concurrent",
                "description": "Disjoint camera and time channels sharing the same actual start.",
                "steps": [
                    {
                        "type": "camera",
                        "start": "after_previous",
                        "delay_ms": 0,
                        "target": {"kind": "group", "name": group},
                        "duration_ms": 1200,
                        "easing": "linear",
                        "fit_padding": 1.15,
                        "distance_scale": 0.82,
                        "orbit": {"enabled": False, "speed_multiplier": 1.0, "persist": True, "direction": 1.0},
                    },
                    {
                        "type": "time",
                        "start": "with_previous",
                        "delay_ms": 0,
                        "direction": "backward",
                        "interval_ms": 160,
                        "stop_after_frames": 5,
                    },
                ],
            },
        ],
    }
    scene["debug_transitions"] = True
    return scene


def _diagnostics_script(root_id: str, *, auto_run: bool) -> str:
    auto_run_js = "true" if auto_run else "false"
    return f"""
    <style>
      #oviz-actions-audit-panel {{ position: fixed; right: 12px; bottom: 12px; z-index: 100000;
        width: min(430px, calc(100vw - 24px)); max-height: 44vh; overflow: auto;
        padding: 10px; border: 1px solid #5d7188; border-radius: 8px;
        background: rgba(4, 10, 18, .94); color: #dbe9f6; font: 12px/1.35 ui-monospace, monospace; }}
      #oviz-actions-audit-panel button {{ margin-bottom: 7px; }}
      #oviz-actions-audit-log {{ white-space: pre-wrap; }}
    </style>
    <aside id="oviz-actions-audit-panel">
      <button id="oviz-actions-audit-run" type="button">Run deterministic audit</button>
      <div id="oviz-actions-audit-status">Waiting for Oviz…</div>
      <pre id="oviz-actions-audit-log"></pre>
    </aside>
    <script>
    (() => {{
      const root = document.getElementById({json.dumps(root_id)});
      const status = document.getElementById("oviz-actions-audit-status");
      const logEl = document.getElementById("oviz-actions-audit-log");
      const runButton = document.getElementById("oviz-actions-audit-run");
      const eventNames = ["states-ready", "transition-start", "transition-progress",
        "transition-complete", "transition-error", "action-start", "action-step-start",
        "action-progress", "action-complete", "action-cancel", "action-error"];
      const events = [];
      const errors = [];
      const frameSamples = [];
      let running = false;
      function log(value) {{
        logEl.textContent += `${{typeof value === "string" ? value : JSON.stringify(value)}}\n`;
        logEl.scrollTop = logEl.scrollHeight;
      }}
      eventNames.forEach((name) => root.addEventListener(name, (event) => {{
        events.push({{ name, detail: event.detail || null, at: performance.now() }});
        if (name.endsWith("error")) errors.push({{ name, detail: event.detail || null }});
      }}));
      window.addEventListener("error", (event) => errors.push({{ name: "window-error", detail: event.message }}));
      window.addEventListener("unhandledrejection", (event) => errors.push({{ name: "unhandled-rejection", detail: String(event.reason) }}));
      function parseDataset(name) {{
        try {{ return root.dataset[name] ? JSON.parse(root.dataset[name]) : null; }}
        catch (_error) {{ return null; }}
      }}
      function require(condition, message) {{
        if (!condition) throw new Error(message);
      }}
      function waitFrames(count = 1) {{
        return new Promise((resolve) => {{
          const step = (remaining) => requestAnimationFrame(() => remaining <= 1 ? resolve() : step(remaining - 1));
          step(Math.max(1, Number(count) || 1));
        }});
      }}
      function sampleFrame() {{
        if (!running) return;
        const canvas = root.querySelector("canvas");
        frameSamples.push({{
          at: performance.now(),
          canvasOpacity: canvas ? Number(getComputedStyle(canvas).opacity) : NaN,
          diagnostics: parseDataset("transitionDiagnostics"),
          retained: parseDataset("retainedSceneDebug"),
          owner: root.dataset.transitionOwner || "",
          runId: root.dataset.transitionRunId || "",
          actionRunId: root.dataset.actionRunId || "",
          actionPhase: root.dataset.actionPhase || "",
        }});
        requestAnimationFrame(sampleFrame);
      }}
      function clickAction(label, options = {{}}) {{
        const button = Array.from(root.querySelectorAll(".oviz-three-action-button"))
          .find((item) => item.textContent.trim() === label);
        if (!button) throw new Error(`Missing Action button: ${{label}}`);
        return new Promise((resolve, reject) => {{
          const timeout = setTimeout(() => reject(new Error(`Timed out: ${{label}}`)), 20000);
          const cleanup = () => {{
            clearTimeout(timeout);
            root.removeEventListener("action-complete", finish);
            root.removeEventListener("action-cancel", cancel);
            root.removeEventListener("action-error", fail);
          }};
          const matches = (event) => {{
            const detail = event.detail || {{}};
            return !detail.actionKey || detail.actionKey === button.dataset.actionKey;
          }};
          const finish = (event) => {{
            if (!matches(event)) return;
            cleanup(); resolve({{ outcome: "complete", detail: event.detail || {{}} }});
          }};
          const cancel = (event) => {{
            if (!matches(event)) return;
            cleanup();
            if (options.allowCancel) resolve({{ outcome: "cancel", detail: event.detail || {{}} }});
            else reject(new Error(`Unexpected Action cancellation: ${{label}}`));
          }};
          const fail = (event) => {{
            if (!matches(event)) return;
            cleanup(); reject(new Error(JSON.stringify(event.detail || {{}})));
          }};
          root.addEventListener("action-complete", finish);
          root.addEventListener("action-cancel", cancel);
          root.addEventListener("action-error", fail);
          button.click();
        }});
      }}
      function orderedPhases(eventSlice) {{
        const result = [];
        eventSlice.forEach((event) => {{
          if (!event.name.startsWith("transition-")) return;
          const phase = event.detail && event.detail.phase;
          if (phase && result.at(-1) !== phase) result.push(phase);
        }});
        return result;
      }}
      function requirePhaseOrder(phases, expected, label) {{
        let cursor = -1;
        expected.forEach((phase) => {{
          const index = phases.indexOf(phase, cursor + 1);
          require(index > cursor, `${{label}} missing ordered phase ${{phase}}: ${{phases.join(",")}}`);
          cursor = index;
        }});
      }}
      function requireMonotonic(values, direction, label) {{
        const finite = values.filter((value) => Number.isFinite(Number(value))).map(Number);
        require(finite.length >= 2, `${{label}} did not produce enough samples`);
        for (let index = 1; index < finite.length; index += 1) {{
          const delta = finite[index] - finite[index - 1];
          if (direction < 0) require(delta <= 1e-5, `${{label}} reversed at ${{index}}`);
          else require(delta >= -1e-5, `${{label}} reversed at ${{index}}`);
        }}
      }}
      function requireStateFidelity(label) {{
        const logical = parseDataset("stateFidelity");
        const rendered = parseDataset("renderedTraceFidelity");
        require(logical && logical.exact === true, `${{label}} logical fidelity failed: ${{JSON.stringify(logical)}}`);
        require(rendered && rendered.exact === true, `${{label}} rendered fidelity failed: ${{JSON.stringify(rendered)}}`);
      }}
      async function run() {{
        if (running) return;
        running = true; runButton.disabled = true; errors.length = 0; events.length = 0;
        frameSamples.length = 0; logEl.textContent = ""; requestAnimationFrame(sampleFrame);
        status.textContent = "Running…";
        try {{
          const states = window.Oviz.get({json.dumps(root_id)}).states;
          await states.original();
          await waitFrames(2);

          const allDomainEventStart = events.length;
          const allDomainSampleStart = frameSamples.length;
          const allDomainResult = await states.goTo("audit-lookback");
          require(!allDomainResult.cancelled, "Original to all-domains State was cancelled");
          await waitFrames(2);
          requireStateFidelity("Original to all-domains");
          const allDomainPhases = orderedPhases(events.slice(allDomainEventStart));
          requirePhaseOrder(allDomainPhases, ["camera+time", "appearance"], "Original to all-domains");
          const allDomainSamples = frameSamples.slice(allDomainSampleStart).filter(
            (sample) => sample.diagnostics && sample.diagnostics.targetId === "audit-lookback"
          );
          requireMonotonic(
            allDomainSamples.filter((sample) => sample.diagnostics.phase === "time")
              .map((sample) => sample.diagnostics.frameValue),
            -1,
            "all-domains timeline",
          );
          const lassoProgress = allDomainSamples
            .filter((sample) => sample.diagnostics.phase === "appearance" && sample.retained)
            .map((sample) => sample.retained.lassoProgress)
            .filter((value) => value !== null && value !== undefined);
          requireMonotonic(lassoProgress, 1, "all-domains lasso membership");

          const reverseResult = await states.previous();
          require(!reverseResult.cancelled, "Previous back to Original was cancelled");
          await waitFrames(2);
          requireStateFidelity("reverse to Original");

          require((await clickAction("Audit: State Present")).outcome === "complete", "State-backed present Action failed");
          requireStateFidelity("state-backed present Action");
          require((await clickAction("Audit: State Lookback")).outcome === "complete", "State-backed lookback Action failed");
          requireStateFidelity("state-backed lookback Action");
          require((await clickAction("Audit: Sequential")).outcome === "complete", "Sequential Action failed");

          const actionToState = clickAction("Audit: Concurrent", {{ allowCancel: true }});
          await waitFrames(8);
          const actionToStateResult = await states.goTo("audit-topology-boundary");
          const actionToStateAction = await actionToState;
          require(actionToStateAction.outcome === "cancel", "Action to State did not cancel the Action");
          require(!actionToStateResult.cancelled, "Action to State target was cancelled");
          requireStateFidelity("Action to State");

          const stateToActionState = states.goTo("audit-lookback");
          await waitFrames(8);
          const stateToActionAction = clickAction("Audit: Sequential");
          const stateToActionStateResult = await stateToActionState;
          require(stateToActionStateResult.cancelled === true, "State to Action did not cancel the State");
          require((await stateToActionAction).outcome === "complete", "State to Action legacy Action failed");

          require((await clickAction("Audit: Concurrent")).outcome === "complete", "Repeated concurrent Action failed");
          require((await clickAction("Audit: Sequential")).outcome === "complete", "Repeated sequential Action failed");

          const rapidOne = states.goTo("audit-lookback");
          await waitFrames(8);
          const rapidTwo = states.goTo("audit-topology-boundary");
          await waitFrames(8);
          const rapidThree = states.goTo("audit-present");
          const rapidResults = await Promise.all([rapidOne, rapidTwo, rapidThree]);
          require(rapidResults[0].cancelled === true, "Rapid State 1 was not cancelled");
          require(rapidResults[1].cancelled === true, "Rapid State 2 was not cancelled");
          require(rapidResults[2].cancelled !== true, "Rapid State 3 did not complete");
          await waitFrames(2);
          requireStateFidelity("rapid final State");

          const canvas = root.querySelector("canvas");
          const canvasOpacity = canvas ? Number(getComputedStyle(canvas).opacity) : NaN;
          require(Math.abs(canvasOpacity - 1) < 1e-6, "Final canvas opacity was not one");
          require(frameSamples.every((sample) => Math.abs(sample.canvasOpacity - 1) < 1e-6), "Canvas opacity dropped during the audit");
          require(errors.length === 0, `Runtime errors: ${{JSON.stringify(errors)}}`);
          const activeStates = states.list().filter((item) => item.active);
          require(activeStates.length === 1 && activeStates[0].id === "audit-present", "Final active State was not audit-present");
          const sequentialCompletes = events.filter(
            (event) => event.name === "action-complete" && event.detail && event.detail.actionKey === "audit-sequential"
          );
          require(sequentialCompletes.length >= 3, "Sequential Action was not exercised repeatedly");
          const actionDiagnostics = events.filter((event) => event.name.startsWith("action-") && event.detail);
          require(actionDiagnostics.every((event) => event.detail.owner === "action" && event.detail.runId), "Action diagnostics lack owner/run ID");
          status.textContent = "PASS";
          log({{ canvasOpacity, allDomainPhases, rapidResults, errors, eventCount: events.length,
            sampleCount: frameSamples.length,
            fidelity: parseDataset("stateFidelity"),
            renderedFidelity: parseDataset("renderedTraceFidelity"),
            transitionMetrics: root.dataset.stateTransitionMetrics || null,
            retainedMetrics: root.dataset.retainedSceneMetrics || null,
            action: {{ phase: root.dataset.actionPhase || null,
              progress: root.dataset.actionProgress || null,
              effectiveDurationMs: root.dataset.actionEffectiveDurationMs || null }} }});
        }} catch (error) {{
          errors.push({{ name: "audit", detail: String(error && error.stack || error) }});
          status.textContent = "FAIL"; log(errors.at(-1));
        }} finally {{ running = false; runButton.disabled = false; }}
      }}
      runButton.addEventListener("click", run);
      root.addEventListener("states-ready", () => {{ status.textContent = "Ready"; if ({auto_run_js}) setTimeout(run, 0); }}, {{ once: true }});
    }})();
    </script>
    """


def render_audit_html(scene_spec: dict[str, Any], *, auto_run: bool = True) -> str:
    figure = ThreeJSFigure(build_audit_scene(scene_spec), compress_scene_spec=True)
    html = figure.to_html(compress_scene_spec=True)
    root_match = re.search(r'<div id="(oviz-three-[^"]+)"', html)
    if root_match is None:
        raise ValueError("Could not find the Oviz root in generated HTML")
    body_end = html.rfind("</body>")
    if body_end < 0:
        raise ValueError("Could not find the outer HTML body terminator")
    return (
        html[:body_end]
        + _diagnostics_script(root_match.group(1), auto_run=auto_run)
        + "\n"
        + html[body_end:]
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--no-auto-run", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    source = args.input.expanduser().resolve()
    output = args.output.expanduser().resolve()
    if source == output:
        raise ValueError("The audit output must not overwrite the canonical input")
    scene_spec = read_embedded_scene_spec(source)
    html = render_audit_html(scene_spec, auto_run=not args.no_auto_run)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(html, encoding="utf-8")
    print(f"Wrote local audit figure: {output}")


if __name__ == "__main__":
    main()
