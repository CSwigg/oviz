from __future__ import annotations

import json
import shutil
import subprocess
import textwrap

import pytest

from oviz.threejs_runtime_actions import THREEJS_ACTION_RUNTIME_JS


def _run_action_runtime(actions: list[dict], body: str) -> dict:
    """Execute the production Actions runtime against a deterministic viewer stub."""

    if shutil.which("node") is None:
        pytest.skip("node is not available")
    scene_spec = {
        "width": 900,
        "height": 600,
        "max_span": 1000,
        "actions": {"enabled": True, "items": actions},
    }
    script = f"""
    "use strict";

    let clockNow = 0;
    Object.defineProperty(globalThis, "performance", {{
      configurable: true,
      value: {{ now: () => clockNow }},
    }});
    globalThis.window = globalThis;
    globalThis.CustomEvent = class CustomEvent {{
      constructor(type, options = {{}}) {{
        this.type = type;
        this.detail = options.detail || {{}};
      }}
    }};

    const emittedEvents = [];
    const capturedErrors = [];
    console.error = (...parts) => capturedErrors.push(parts.map(String).join(" "));
    const root = {{
      id: "viewer-under-test",
      dataset: {{}},
      clientWidth: 900,
      clientHeight: 600,
      dispatchEvent: (event) => emittedEvents.push({{ name: event.type, detail: event.detail }}),
    }};
    const sceneSpec = {json.dumps(scene_spec)};
    const actionBarEl = null;
    const groupSelectEl = null;
    const legendItems = [
      {{ key: "young", showlegend: true }},
      {{ key: "old", showlegend: true }},
    ];
    const frameSpecs = [
      {{ time: -1, traces: [{{ key: "young", showlegend: true }}] }},
      {{ time: 0, traces: [{{ key: "young", showlegend: true }}, {{ key: "old", showlegend: true }}] }},
    ];
    let currentGroup = "Default";
    let legendState = {{ young: true, old: false }};
    let currentFrameIndex = 0;
    let displayedFrameValue = 0;
    let playbackDirection = 0;
    let lastPlaybackAdvanceTimestamp = null;
    let playbackIntervalMs = 100;
    let minimalModeEnabled = false;
    let cameraAutoOrbitEnabled = false;
    let cameraAutoOrbitSpeedMultiplier = 1;
    let cameraAutoOrbitDirection = 1;
    let nearbyRegionLabelsVisible = true;
    let cameraViewMode = "free";
    let milkyWayViewOpacityScale = 1.0;
    let skyDomeViewOpacityScale = 0.0;
    let ovizHeldSelectionTransition = null;
    let ovizStateControllerReady = true;
    let ovizStateTransition = null;
    let ovizActiveStateId = null;
    const ovizStatesProject = {{ items: [] }};
    const renderCalls = [];
    const cancelledStateTransitions = [];
    let captureSnapshot = {{ source: "captureRuntimeState" }};
    let stateBeginHook = (target) => Promise.resolve({{ cancelled: false, target }});

    function clampRange(value, minimum, maximum) {{
      return Math.min(Math.max(Number(value), Number(minimum)), Number(maximum));
    }}
    function clampFrameIndex(value) {{
      return Math.min(Math.max(Math.round(Number(value) || 0), 0), frameSpecs.length - 1);
    }}
    function clampFrameValue(value) {{
      return Math.min(Math.max(Number(value) || 0, 0), frameSpecs.length - 1);
    }}
    function frameTimeForValue(value) {{
      const index = clampFrameIndex(value);
      return Number(frameSpecs[index].time);
    }}
    function currentFrame() {{ return frameSpecs[clampFrameIndex(currentFrameIndex)]; }}
    function groupDefaults(name) {{
      return name === "All" ? {{ young: true, old: true }} : {{ young: true, old: "legendonly" }};
    }}
    function isGalacticReferenceTrace() {{ return false; }}
    function galacticReferenceMotionVisible() {{ return true; }}
    function isNearbyRegionLabelTrace() {{ return false; }}
    function traceVisible(trace) {{
      return Object.prototype.hasOwnProperty.call(legendState, trace.key)
        ? Boolean(legendState[trace.key])
        : groupDefaults(currentGroup)[trace.key] === true;
    }}
    function stateTraceVisibilityState() {{ return null; }}
    function ovizTraceCandidatesForSnapshots(fromSnapshot, toSnapshot) {{
      const candidates = [];
      [fromSnapshot, toSnapshot].forEach((snapshot) => {{
        const frameValue = snapshot && snapshot.timeline
          ? snapshot.timeline.frame_value
          : displayedFrameValue;
        const frame = frameSpecs[clampFrameIndex(frameValue)];
        (frame && frame.traces || []).forEach((trace) => candidates.push(trace));
      }});
      return candidates;
    }}
    function renderLegend() {{ renderCalls.push({{ kind: "legend", at: clockNow }}); }}
    function renderFrame(index) {{ renderCalls.push({{ kind: "frame", index, at: clockNow }}); }}
    function renderInterpolatedFrameValue(value, options = {{}}) {{
      renderCalls.push({{ kind: "interpolated", value, options, at: clockNow }});
    }}
    function updateTimelineMotionOpacity() {{}}
    function updatePlaybackButtons() {{}}
    function syncLegendGroupChooser() {{}}
    function focusViewer() {{}}
    function setCameraAutoOrbitEnabled(enabled) {{ cameraAutoOrbitEnabled = Boolean(enabled); }}
    function setMilkyWayModelOpacityScale(value) {{ milkyWayViewOpacityScale = Number(value); }}
    function setSkyDomeViewOpacityScale(value) {{ skyDomeViewOpacityScale = Number(value); }}
    function releaseOvizHeldSelectionTransition() {{ ovizHeldSelectionTransition = null; return true; }}
    function galacticSimpleDefaultOrbitActive() {{ return false; }}
    function ovizStateEvent(name, detail) {{ emittedEvents.push({{ name, detail }}); }}
    function safeJsonClone(value, fallback) {{
      try {{ return JSON.parse(JSON.stringify(value)); }} catch (_error) {{ return fallback; }}
    }}
    function captureRuntimeState() {{ return safeJsonClone(captureSnapshot, {{}}); }}
    function ovizStateTargetFor(state) {{
      return {{ id: String(state), index: 0, name: String(state), snapshot: {{ state }} }};
    }}
    function ovizBeginStateTransition(target, options = {{}}) {{
      return stateBeginHook(target, options);
    }}
    function ovizCancelStateTransitionWithoutSnap(reason, options = {{}}) {{
      cancelledStateTransitions.push({{ reason, options }});
      ovizStateTransition = null;
      return true;
    }}

    {THREEJS_ACTION_RUNTIME_JS}

    (async () => {{
      {body}
    }})().catch((error) => {{
      process.stderr.write(String(error && error.stack || error));
      process.exit(1);
    }});
    """
    result = subprocess.run(
        ["node"],
        input=textwrap.dedent(script),
        text=True,
        capture_output=True,
        check=False,
    )
    if result.returncode:
        raise AssertionError(result.stderr or result.stdout)
    return json.loads(result.stdout)


def test_restore_delegates_a_complete_runtime_snapshot_without_crashing():
    payload = _run_action_runtime(
        [],
        r"""
        const fullSnapshot = {
          camera: { position: { x: 1, y: 2, z: 3 }, target: { x: 4, y: 5, z: 6 } },
          global_controls: { camera_fov: 42, camera_view_mode: "free", zen_mode: true },
          timeline: { frame_value: 0.375, playback_direction: -1, playback_interval_ms: 77 },
          traces: { young: { visible: true, opacity: 0.65, color: "#ff0000" } },
          legend: { group: "Default", state: { young: true, old: false } },
          selections: { points: ["cluster-1"] },
          lasso: { selected_keys: ["cluster-1"], inverted: false },
          volumes: { dust: { opacity: 0.2, lasso_weights: [1, 0.5, 0] } },
          panels: { legend: { left: 12, top: 24, width: 220, height: 160 } },
          aladin: { layers: [{ id: "dust", opacity: 0.7 }] },
        };
        initialActionViewState = fullSnapshot;
        clockNow = 42;
        let submitted = null;
        let resolveState = null;
        stateBeginHook = (target, options) => {
          submitted = { target, options };
          ovizStateTransition = { phasePlan: { effectiveDurationMs: 2400 }, startedAt: clockNow };
          return new Promise((resolve) => {
            resolveState = (result) => {
              ovizStateTransition = null;
              resolve(result);
            };
          });
        };

        const started = startRestoreInitialView();
        updateViewerActions(clockNow, { schedulerOnly: true });
        const runId = activeActionRun && activeActionRun.id;
        const submittedBeforeResolve = safeJsonClone(submitted, null);
        const sameReference = submitted && submitted.target.snapshot === fullSnapshot;
        resolveState({ cancelled: false });
        await Promise.resolve();
        await Promise.resolve();
        updateViewerActions(clockNow, { schedulerOnly: true });

        process.stdout.write(JSON.stringify({
          started,
          runId,
          submitted: submittedBeforeResolve,
          sameReference,
          activeAfterResolve: activeActionRun && activeActionRun.id,
          completed: emittedEvents.some((event) => event.name === "action-complete"),
          errors: capturedErrors,
        }));
        """,
    )

    assert payload["started"] is True
    assert payload["runId"].startswith("action-")
    assert payload["sameReference"] is False
    assert payload["submitted"]["options"] == {"preserveActionRun": True}
    snapshot = payload["submitted"]["target"]["snapshot"]
    assert snapshot["timeline"]["frame_value"] == pytest.approx(0.375)
    assert snapshot["camera"]["position"] == {"x": 1, "y": 2, "z": 3}
    assert snapshot["traces"]["young"]["opacity"] == pytest.approx(0.65)
    assert snapshot["lasso"]["selected_keys"] == ["cluster-1"]
    assert snapshot["volumes"]["dust"]["lasso_weights"] == [1, 0.5, 0]
    assert snapshot["panels"]["legend"]["width"] == 220
    assert payload["activeAfterResolve"] is None
    assert payload["completed"] is True
    assert payload["errors"] == []


def test_scheduler_uses_actual_predecessor_start_and_completion_times():
    payload = _run_action_runtime(
        [],
        r"""
        const starts = [];
        startCameraAction = (_step, index, now) => starts.push({ type: "camera", index, now });
        startTimeAction = (_step, index, now) => starts.push({ type: "time", index, now });
        startLegendGroupAction = (_step, index, now) => starts.push({ type: "legend", index, now });

        const run = beginActionRun({
          key: "dependency-test",
          steps: [
            { type: "camera", start: "after_previous", delay_ms: 30 },
            { type: "time", start: "with_previous", delay_ms: 20 },
            { type: "legend_group", start: "after_previous", delay_ms: 15 },
          ],
        }, 100);

        updateViewerActions(100, { schedulerOnly: true });
        updateViewerActions(129, { schedulerOnly: true });
        updateViewerActions(130, { schedulerOnly: true });
        updateViewerActions(149, { schedulerOnly: true });
        updateViewerActions(150, { schedulerOnly: true });
        markActionStepFinished(1, 195, run.id);
        updateViewerActions(209, { schedulerOnly: true });
        updateViewerActions(210, { schedulerOnly: true });

        process.stdout.write(JSON.stringify({
          starts,
          records: run.steps.map((record) => ({
            index: record.index,
            startedAt: record.startedAt,
            finishedAt: record.finishedAt,
          })),
        }));
        """,
    )

    assert payload["starts"] == [
        {"type": "camera", "index": 0, "now": 130},
        {"type": "time", "index": 1, "now": 150},
        {"type": "legend", "index": 2, "now": 210},
    ]
    assert payload["records"][0]["startedAt"] == 130
    assert payload["records"][1]["startedAt"] == 150
    assert payload["records"][1]["finishedAt"] == 195
    assert payload["records"][2]["startedAt"] == 210


def test_domain_locks_allow_disjoint_concurrency_and_block_collisions():
    payload = _run_action_runtime(
        [],
        r"""
        const starts = [];
        startCameraAction = (_step, index, now) => starts.push({ type: "camera", index, now });
        startTimeAction = (_step, index, now) => starts.push({ type: "time", index, now });

        const run = beginActionRun({
          key: "domain-test",
          steps: [
            { type: "camera", delay_ms: 0 },
            { type: "time", start: "with_previous", delay_ms: 0 },
            { type: "camera", start: "with_previous", delay_ms: 0 },
          ],
        }, 100);
        updateViewerActions(100, { schedulerOnly: true });
        const ownersWhileBlocked = Array.from(run.domainOwners.entries());
        const thirdStartedWhileBlocked = run.steps[2].started;
        markActionStepFinished(0, 120, run.id);
        updateViewerActions(120, { schedulerOnly: true });

        process.stdout.write(JSON.stringify({
          starts,
          ownersWhileBlocked,
          thirdStartedWhileBlocked,
          thirdStartedAt: run.steps[2].startedAt,
        }));
        """,
    )

    assert payload["starts"][:2] == [
        {"type": "camera", "index": 0, "now": 100},
        {"type": "time", "index": 1, "now": 100},
    ]
    assert payload["thirdStartedWhileBlocked"] is False
    assert sorted(payload["ownersWhileBlocked"]) == [["camera", 0], ["time", 1]]
    assert payload["starts"][2] == {"type": "camera", "index": 2, "now": 120}
    assert payload["thirdStartedAt"] == 120


def test_replaced_action_ignores_stale_state_completion_callbacks():
    actions = [
        {
            "key": "state-a",
            "label": "State A",
            "steps": [{"type": "state", "runtime_snapshot": {"marker": "A"}}],
        },
        {
            "key": "state-b",
            "label": "State B",
            "steps": [{"type": "state", "runtime_snapshot": {"marker": "B"}}],
        },
    ]
    payload = _run_action_runtime(
        actions,
        r"""
        const deferred = [];
        stateBeginHook = (target) => {
          let resolve = null;
          const promise = new Promise((resolver) => { resolve = resolver; });
          deferred.push({ target, resolve });
          ovizStateTransition = { phasePlan: { effectiveDurationMs: 800 }, startedAt: clockNow };
          return promise;
        };

        clockNow = 0;
        startViewerAction("state-a");
        updateViewerActions(clockNow, { schedulerOnly: true });
        const firstRunId = activeActionRun.id;

        clockNow = 10;
        startViewerAction("state-b");
        updateViewerActions(clockNow, { schedulerOnly: true });
        const secondRunId = activeActionRun.id;
        deferred[0].resolve({ cancelled: false });
        await Promise.resolve();
        await Promise.resolve();
        const secondFinishedAfterStaleCallback = activeActionRun.steps[0].finished;

        ovizStateTransition = null;
        deferred[1].resolve({ cancelled: false });
        await Promise.resolve();
        await Promise.resolve();
        updateViewerActions(clockNow, { schedulerOnly: true });

        process.stdout.write(JSON.stringify({
          firstRunId,
          secondRunId,
          targets: deferred.map((item) => item.target.snapshot.marker),
          secondFinishedAfterStaleCallback,
          activeAfterSecondCompletion: activeActionRun && activeActionRun.id,
          cancels: emittedEvents.filter((event) => event.name === "action-cancel"),
          completes: emittedEvents.filter((event) => event.name === "action-complete"),
        }));
        """,
    )

    assert payload["firstRunId"] != payload["secondRunId"]
    assert payload["targets"] == ["A", "B"]
    assert payload["secondFinishedAfterStaleCallback"] is False
    assert payload["activeAfterSecondCompletion"] is None
    assert len(payload["cancels"]) == 1
    assert payload["cancels"][0]["detail"]["runId"] == payload["firstRunId"]
    assert payload["cancels"][0]["detail"]["reason"] == "replaced"
    assert len(payload["completes"]) == 1
    assert payload["completes"][0]["detail"]["runId"] == payload["secondRunId"]


def test_scheduler_contains_step_errors_and_can_run_the_next_action():
    actions = [
        {"key": "bad", "steps": [{"type": "time", "fail_for_test": True}]},
        {"key": "good", "steps": [{"type": "time", "fail_for_test": False}]},
    ]
    payload = _run_action_runtime(
        actions,
        r"""
        startTimeAction = (step, index, now) => {
          if (step.fail_for_test) throw new Error("intentional scheduler failure");
          markActionStepFinished(index, now, activeActionRun.id);
        };

        clockNow = 5;
        const badStarted = startViewerAction("bad");
        const badUpdateReturned = updateViewerActions(clockNow, { schedulerOnly: true });
        const activeAfterFailure = activeActionRun && activeActionRun.id;
        const statusAfterFailure = root.dataset.actionStatus;

        clockNow = 10;
        const goodStarted = startViewerAction("good");
        updateViewerActions(clockNow, { schedulerOnly: true });
        const activeAfterRecovery = activeActionRun && activeActionRun.id;

        process.stdout.write(JSON.stringify({
          badStarted,
          badUpdateReturned,
          activeAfterFailure,
          statusAfterFailure,
          goodStarted,
          activeAfterRecovery,
          errorEvents: emittedEvents.filter((event) => event.name === "action-error"),
          completeEvents: emittedEvents.filter((event) => event.name === "action-complete"),
          capturedErrors,
        }));
        """,
    )

    assert payload["badStarted"] is True
    assert payload["badUpdateReturned"] is False
    assert payload["activeAfterFailure"] is None
    assert payload["statusAfterFailure"] == "error"
    assert len(payload["errorEvents"]) == 1
    assert "intentional scheduler failure" in payload["errorEvents"][0]["detail"]["error"]
    assert payload["goodStarted"] is True
    assert payload["activeAfterRecovery"] is None
    assert len(payload["completeEvents"]) == 1
    assert payload["capturedErrors"]


def test_legacy_action_retargets_from_the_live_state_trace_opacity():
    actions = [
        {
            "key": "continue-time",
            "steps": [
                {
                    "type": "time",
                    "direction": "forward",
                    "stop_after_ms": 100,
                    "interval_ms": 100,
                }
            ],
        }
    ]
    payload = _run_action_runtime(
        actions,
        r"""
        stateTraceVisibilityState = (trace) => ({
          visible: true,
          opacity: trace.key === "young" ? 0.35 : 0.7,
        });
        ovizStateTransition = {
          fromSnapshot: { timeline: { frame_value: 0 } },
          targetSnapshot: { timeline: { frame_value: 1 } },
          phasePlan: { effectiveDurationMs: 800 },
        };

        clockNow = 25;
        const started = startViewerAction("continue-time");
        updateViewerActions(clockNow, { schedulerOnly: true });

        process.stdout.write(JSON.stringify({
          started,
          fromOpacityByKey: legendTransitionState && legendTransitionState.fromOpacityByKey,
          cancelledStateTransitions,
        }));
        """,
    )

    assert payload["started"] is True
    assert payload["fromOpacityByKey"]["young"] == pytest.approx(0.35)
    assert payload["fromOpacityByKey"]["old"] == pytest.approx(0.7)
    assert payload["cancelledStateTransitions"] == [
        {
            "reason": "action-start",
            "options": {
                "restorePresentation": False,
                "preserveRenderedSelection": True,
            },
        }
    ]


def test_legend_action_preserves_fractional_timeline_position_while_fading():
    actions = [
        {
            "key": "show-all",
            "steps": [
                {
                    "type": "legend_group",
                    "group": "All",
                    "duration_ms": 100,
                    "easing": "linear",
                }
            ],
        }
    ]
    payload = _run_action_runtime(
        actions,
        r"""
        displayedFrameValue = 0.4;
        currentFrameIndex = 0;
        clockNow = 0;
        const started = startViewerAction("show-all");
        updateViewerActions(clockNow, { schedulerOnly: false });
        clockNow = 50;
        updateViewerActions(clockNow, { schedulerOnly: false });
        clockNow = 100;
        updateViewerActions(clockNow, { schedulerOnly: false });

        process.stdout.write(JSON.stringify({
          started,
          displayedFrameValue,
          currentFrameIndex,
          renderCalls,
        }));
        """,
    )

    assert payload["started"] is True
    assert payload["displayedFrameValue"] == pytest.approx(0.4)
    assert payload["currentFrameIndex"] == 0
    assert not [call for call in payload["renderCalls"] if call["kind"] == "frame"]
    interpolated = [call for call in payload["renderCalls"] if call["kind"] == "interpolated"]
    assert interpolated
    assert all(call["value"] == pytest.approx(0.4) for call in interpolated)


def test_sky_camera_action_exclusively_owns_domains_before_with_previous_time():
    actions = [
        {
            "key": "sky-camera-and-time",
            "steps": [
                {"type": "camera", "duration_ms": 800},
                {
                    "type": "time",
                    "start": "with_previous",
                    "direction": "forward",
                    "stop_after_ms": 100,
                    "interval_ms": 100,
                },
            ],
        }
    ]
    payload = _run_action_runtime(
        actions,
        r"""
        cameraViewMode = "earth";
        const starts = [];
        startCameraAction = (_step, index) => starts.push({ type: "camera", index });
        startTimeAction = (_step, index, now) => {
          starts.push({ type: "time", index, now });
          markActionStepFinished(index, now, activeActionRun.id);
        };

        clockNow = 0;
        startViewerAction("sky-camera-and-time");
        updateViewerActions(clockNow, { schedulerOnly: true });
        const afterFirstFrame = {
          starts: safeJsonClone(starts, []),
          domains: activeActionRun.steps.map((record) => record.domains),
          timeStarted: activeActionRun.steps[1].started,
        };
        clockNow = 100;
        markActionStepFinished(0, clockNow, activeActionRun.id);
        updateViewerActions(clockNow, { schedulerOnly: true });

        process.stdout.write(JSON.stringify({
          afterFirstFrame,
          starts,
          activeAfterCompletion: activeActionRun && activeActionRun.id,
        }));
        """,
    )

    assert payload["afterFirstFrame"]["starts"] == [{"type": "camera", "index": 0}]
    assert payload["afterFirstFrame"]["domains"][0] == [
        "camera",
        "appearance",
        "selection",
        "time",
    ]
    assert payload["afterFirstFrame"]["timeStarted"] is False
    assert payload["starts"][1]["type"] == "time"
    assert payload["activeAfterCompletion"] is None


def test_concurrent_legend_and_time_action_renders_scene_once_per_update():
    actions = [
        {
            "key": "appearance-and-time",
            "steps": [
                {
                    "type": "legend_group",
                    "group": "All",
                    "duration_ms": 100,
                    "easing": "linear",
                },
                {
                    "type": "time",
                    "start": "with_previous",
                    "direction": "forward",
                    "stop_after_ms": 100,
                    "interval_ms": 100,
                },
            ],
        }
    ]
    payload = _run_action_runtime(
        actions,
        r"""
        displayedFrameValue = 0;
        currentFrameIndex = 0;
        clockNow = 0;
        startViewerAction("appearance-and-time");
        updateViewerActions(clockNow, { schedulerOnly: false });
        const firstFrameCalls = safeJsonClone(renderCalls, []);
        renderCalls.length = 0;
        clockNow = 50;
        updateViewerActions(clockNow, { schedulerOnly: false });
        const secondFrameCalls = safeJsonClone(renderCalls, []);

        process.stdout.write(JSON.stringify({ firstFrameCalls, secondFrameCalls }));
        """,
    )

    assert len([call for call in payload["firstFrameCalls"] if call["kind"] == "interpolated"]) == 1
    assert len([call for call in payload["secondFrameCalls"] if call["kind"] == "interpolated"]) == 1


def test_state_interruption_rolls_back_held_presentation_and_selection_then_releases_it():
    actions = [
        {
            "key": "long-time",
            "steps": [
                {
                    "type": "time",
                    "direction": "forward",
                    "stop_after_ms": 500,
                    "interval_ms": 500,
                }
            ],
        }
    ]
    payload = _run_action_runtime(
        actions,
        r"""
        actionHeldAppearanceRollback = {
          fromSnapshot: {},
          targetSnapshot: {},
          frozenAppearanceProgress: 0.6,
          currentAppearanceProgress: 0.6,
          startedAt: 0,
          durationMs: 240,
          fromMilkyWayOpacity: 0.4,
          fromSkyOpacity: 0.6,
          toMilkyWayOpacity: 1.0,
          toSkyOpacity: 0.0,
        };
        ovizHeldSelectionTransition = { progress: 0.0 };

        clockNow = 0;
        startViewerAction("long-time");
        updateViewerActions(clockNow, { schedulerOnly: false });
        clockNow = 120;
        updateViewerActions(clockNow, { schedulerOnly: false });
        const midpoint = {
          appearance: actionHeldAppearanceRollback.currentAppearanceProgress,
          selection: ovizHeldSelectionTransition.progress,
          milkyWay: milkyWayViewOpacityScale,
          sky: skyDomeViewOpacityScale,
        };
        clockNow = 240;
        updateViewerActions(clockNow, { schedulerOnly: false });

        process.stdout.write(JSON.stringify({
          midpoint,
          heldAppearance: actionHeldAppearanceRollback,
          heldSelection: ovizHeldSelectionTransition,
          milkyWayViewOpacityScale,
          skyDomeViewOpacityScale,
          activeRun: activeActionRun && activeActionRun.id,
        }));
        """,
    )

    assert payload["midpoint"]["appearance"] == pytest.approx(0.3)
    assert payload["midpoint"]["selection"] == pytest.approx(0.5)
    assert payload["midpoint"]["milkyWay"] == pytest.approx(0.7)
    assert payload["midpoint"]["sky"] == pytest.approx(0.3)
    assert payload["heldAppearance"] is None
    assert payload["heldSelection"] is None
    assert payload["milkyWayViewOpacityScale"] == pytest.approx(1.0)
    assert payload["skyDomeViewOpacityScale"] == pytest.approx(0.0)
    assert payload["activeRun"] is not None
