from __future__ import annotations

import importlib
import json
import shutil
import subprocess
import textwrap

import pytest


def _transition_runtime_js() -> str:
    """Return the production DOM-free runtime, or skip only these Node tests."""

    try:
        module = importlib.import_module("oviz.threejs_transition_runtime")
    except ModuleNotFoundError:
        pytest.skip("the unified transition runtime has not been integrated yet")
    runtime_js = getattr(module, "THREEJS_TRANSITION_RUNTIME_JS", "")
    if not isinstance(runtime_js, str) or not runtime_js.strip():
        pytest.skip("the unified transition runtime JavaScript is not available")
    return runtime_js


def _run_node(body: str) -> dict:
    if shutil.which("node") is None:
        pytest.skip("node is not available")
    runtime_js = _transition_runtime_js()
    script = f"""
    "use strict";
    {runtime_js}
    const runtime = globalThis.OvizTransitionRuntime;
    if (!runtime || typeof runtime !== "object") {{
      throw new Error("THREEJS_TRANSITION_RUNTIME_JS must expose globalThis.OvizTransitionRuntime");
    }}
    (async () => {{
      {body}
    }})().catch((error) => {{
      console.error(error && error.stack ? error.stack : error);
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


def test_trace_union_and_point_matching_use_stable_keys():
    payload = _run_node(
        r"""
        const source = [
          { key: "a", points: [
            { motion: { key: "motion-1" }, selection: { cluster_name: "wrong-source" } },
            { selection: { cluster_name: "selection-2" } },
          ] },
          { key: "b", points: [] },
        ];
        const target = [
          { key: "c", points: [] },
          { key: "a", points: [
            { selection: { cluster_name: "selection-2" } },
            { motion: { key: "motion-1" }, selection: { cluster_name: "wrong-target" } },
          ] },
        ];
        const union = runtime.buildTraceUnion(source, target);
        const keys = union.map((entry) => String(entry.key || entry.traceKey));
        const pairs = runtime.matchTracePoints("a", source[0].points, target[1].points);
        const normalizedPairs = pairs.map((pair) => ({
          source: Number(pair.sourceIndex ?? pair.fromIndex),
          target: Number(pair.targetIndex ?? pair.toIndex),
          kind: String(pair.kind || pair.matchKind || pair.matchedBy || ""),
          key: String(pair.key || ""),
        }));
        process.stdout.write(JSON.stringify({ keys, normalizedPairs }));
        """
    )

    assert payload["keys"] == ["a", "b", "c"]
    assert {tuple((item["source"], item["target"])) for item in payload["normalizedPairs"]} == {
        (0, 1),
        (1, 0),
    }
    assert any(item["key"].startswith("motion:") for item in payload["normalizedPairs"])
    assert any(item["key"].startswith("selection:") for item in payload["normalizedPairs"])
    assert all(item["kind"] == "stable" for item in payload["normalizedPairs"])


def test_trace_union_keeps_destination_only_traces_at_t0():
    payload = _run_node(
        r"""
        const source = Array.from({ length: 19 }, (_, index) => ({ key: `shared-${index}`, points: [] }));
        const target = [
          ...source.slice().reverse(),
          ...Array.from({ length: 4 }, (_, index) => ({ key: `present-only-${index}`, points: [] })),
        ];
        const union = runtime.buildTraceUnion(source, target);
        process.stdout.write(JSON.stringify({
          count: union.length,
          keys: union.map((entry) => String(entry.key || entry.traceKey)),
        }));
        """
    )

    assert payload["count"] == 23
    assert payload["keys"][:19] == [f"shared-{index}" for index in range(19)]
    assert payload["keys"][-4:] == [f"present-only-{index}" for index in range(4)]


def test_visibility_defaults_and_membership_weights_are_continuous():
    payload = _run_node(
        r"""
        const trace = { key: "young", showlegend: true };
        const groupVisibility = { Clusters: { young: true, old: false } };
        const omitted = runtime.resolveTraceVisibility(trace, groupVisibility.Clusters, {});
        const explicitOff = runtime.resolveTraceVisibility(trace, groupVisibility.Clusters, { young: false });
        const samples = [0, 0.25, 0.5, 0.75, 1].map((progress) => ({
          progress,
          leaving: runtime.membershipWeight("a", null, new Set(["a"]), new Set(), progress),
          entering: runtime.membershipWeight("b", null, new Set(), new Set(["b"]), progress),
          shared: runtime.membershipWeight("c", null, new Set(["c"]), new Set(["c"]), progress),
          fromUnfiltered: runtime.membershipWeight("d", null, null, new Set(), progress),
          toUnfiltered: runtime.membershipWeight("e", null, new Set(), null, progress),
        }));
        process.stdout.write(JSON.stringify({ omitted, explicitOff, samples }));
        """
    )

    assert bool(payload["omitted"]) is True
    assert bool(payload["explicitOff"]) is False
    samples = payload["samples"]
    assert [item["leaving"] for item in samples] == pytest.approx([1, 0.75, 0.5, 0.25, 0])
    assert [item["entering"] for item in samples] == pytest.approx([0, 0.25, 0.5, 0.75, 1])
    assert [item["shared"] for item in samples] == pytest.approx([1, 1, 1, 1, 1])
    assert [item["fromUnfiltered"] for item in samples] == pytest.approx([1, 0.75, 0.5, 0.25, 0])
    assert [item["toUnfiltered"] for item in samples] == pytest.approx([0, 0.25, 0.5, 0.75, 1])


def test_session_runs_one_ordered_phase_per_fake_animation_frame():
    payload = _run_node(
        r"""
        let now = 0;
        let nextFrameId = 1;
        const pending = new Map();
        const calls = [];
        const events = [];
        const requestFrame = (callback) => {
          const id = nextFrameId++;
          pending.set(id, callback);
          return id;
        };
        const cancelFrame = (id) => pending.delete(id);
        async function advance(milliseconds) {
          now += milliseconds;
          const callbacks = Array.from(pending.values());
          pending.clear();
          callbacks.forEach((callback) => callback(now));
          await Promise.resolve();
        }
        const adapter = {
          now: () => now,
          requestFrame,
          cancelFrame,
          capture: () => ({ frame: 0, camera: { x: 0 } }),
          prepare: (context) => calls.push({ type: "prepare", phase: context.phase || "" }),
          applyCamera: (progress, context) => calls.push({ type: "camera", progress, frame: now, runId: context.runId }),
          applyAppearance: (progress, context) => calls.push({ type: "appearance", progress, frame: now, runId: context.runId }),
          applyTime: (progress, context) => calls.push({ type: "time", progress, frame: now, runId: context.runId }),
          finalize: (_target, context) => calls.push({ type: "finalize", frame: now, runId: context.runId }),
          validate: (_target, context) => { calls.push({ type: "validate", frame: now, runId: context.runId }); return []; },
          emit: (name, detail) => events.push({ name, detail }),
        };
        const session = runtime.createSession(adapter);
        const completion = session.start({
          owner: "state",
          runId: "ordered",
          target: { frame: 5 },
          easing: "linear",
          phases: [
            { name: "camera", durationMs: 80 },
            { name: "appearance", durationMs: 80 },
            { name: "time", durationMs: 80 },
          ],
        });
        for (let index = 0; index < 40 && pending.size; index += 1) await advance(20);
        await completion;
        // Exact settle validation is intentionally required on the frame after finalize.
        if (pending.size) await advance(20);
        const visualCalls = calls.filter((item) => ["camera", "appearance", "time"].includes(item.type));
        process.stdout.write(JSON.stringify({ calls, events, visualCalls }));
        """
    )

    calls = payload["calls"]
    camera = [item for item in calls if item["type"] == "camera"]
    appearance = [item for item in calls if item["type"] == "appearance"]
    time = [item for item in calls if item["type"] == "time"]
    assert camera and appearance and time
    assert max(item["frame"] for item in camera) <= min(item["frame"] for item in appearance)
    assert max(item["frame"] for item in appearance) <= min(item["frame"] for item in time)
    assert [item["progress"] for item in camera] == sorted(item["progress"] for item in camera)
    assert [item["progress"] for item in appearance] == sorted(item["progress"] for item in appearance)
    assert [item["progress"] for item in time] == sorted(item["progress"] for item in time)
    assert sum(item["type"] == "finalize" for item in calls) == 1
    assert sum(item["type"] == "validate" for item in calls) == 1
    assert all(item["runId"] == "ordered" for item in payload["visualCalls"])
    progress_events = [item for item in payload["events"] if item["name"].endswith("progress")]
    assert progress_events
    assert all("phase" in item["detail"] for item in progress_events)
    assert all("phaseProgress" in item["detail"] for item in progress_events)
    assert all("effectiveDurationMs" in item["detail"] for item in progress_events)


def test_session_finishes_every_crossed_phase_and_settle_has_no_time_callback():
    payload = _run_node(
        r"""
        let now = 0;
        let nextFrameId = 1;
        const pending = new Map();
        const calls = [];
        const requestFrame = (callback) => { const id = nextFrameId++; pending.set(id, callback); return id; };
        const cancelFrame = (id) => pending.delete(id);
        async function advance(milliseconds) {
          now += milliseconds;
          const callbacks = Array.from(pending.values());
          pending.clear();
          callbacks.forEach((callback) => callback(now));
          await Promise.resolve();
          await Promise.resolve();
        }
        const session = runtime.createSession({
          now: () => now,
          requestFrame,
          cancelFrame,
          applyCamera: (progress) => calls.push({ type: "camera", progress, now }),
          applyAppearance: (progress) => calls.push({ type: "appearance", progress, now }),
          applyTime: (progress) => calls.push({ type: "time", progress, now }),
          finalize: () => calls.push({ type: "finalize", now }),
          validate: () => [],
        });
        const completion = session.start({
          owner: "state",
          runId: "long-frame",
          easing: "linear",
          phases: [
            { name: "camera", durationMs: 50 },
            { name: "appearance", durationMs: 50 },
            { name: "time", durationMs: 50 },
          ],
        });
        await advance(10);
        // Jump across both phase boundaries in a single animation frame.
        await advance(115);
        await advance(25);
        await completion;

        const beforeSettle = calls.slice();
        const settleCompletion = session.start({
          owner: "state",
          runId: "settle-only",
          phases: [],
        });
        await advance(0);
        await settleCompletion;
        process.stdout.write(JSON.stringify({ calls, beforeSettle }));
        """
    )

    crossed = [item for item in payload["beforeSettle"] if item["now"] == 125]
    assert [(item["type"], item["progress"]) for item in crossed] == [
        ("camera", 1),
        ("appearance", 1),
        ("time", 0.5),
    ]
    assert any(item["type"] == "time" and item["progress"] == 1 for item in payload["beforeSettle"])
    assert [item for item in payload["calls"][len(payload["beforeSettle"]):] if item["type"] == "time"] == []


def test_session_contains_prepare_and_cancel_cleanup_exceptions():
    payload = _run_node(
        r"""
        const prepareEvents = [];
        const prepareSession = runtime.createSession({
          now: () => 0,
          requestFrame: () => 1,
          prepare: () => { throw new Error("prepare-failed"); },
          emit: (name, detail) => prepareEvents.push({ name, detail }),
        });
        const prepareError = await prepareSession.start({
          owner: "state", runId: "prepare", phases: [{ name: "camera", durationMs: 20 }],
        }).then(() => "", (error) => String(error && error.message || error));

        const cancelEvents = [];
        const cancelSession = runtime.createSession({
          now: () => 0,
          requestFrame: () => 2,
          cancelFrame: () => {},
          cancel: () => { throw new Error("cancel-cleanup-failed"); },
          emit: (name, detail) => cancelEvents.push({ name, detail }),
        });
        const cancelled = cancelSession.start({
          owner: "action", runId: "cancel", phases: [{ name: "time", durationMs: 20 }],
        });
        const cancelReturned = cancelSession.cancel("manual");
        const cancelResult = await cancelled;
        process.stdout.write(JSON.stringify({
          prepareError,
          prepareCurrent: prepareSession.current(),
          prepareEvents,
          cancelReturned,
          cancelCurrent: cancelSession.current(),
          cancelResult,
          cancelEvents,
        }));
        """
    )

    assert payload["prepareError"] == "prepare-failed"
    assert payload["prepareCurrent"] is None
    assert any(
        item["name"] == "transition-error" and item["detail"].get("stage") == "prepare"
        for item in payload["prepareEvents"]
    )
    assert payload["cancelReturned"] is True
    assert payload["cancelCurrent"] is None
    assert payload["cancelResult"]["cancelled"] is True
    assert payload["cancelResult"]["cleanupError"] == "cancel-cleanup-failed"
    assert any(
        item["name"] == "transition-error" and item["detail"].get("stage") == "cancel"
        for item in payload["cancelEvents"]
    )
    assert any(item["name"] == "transition-cancel" for item in payload["cancelEvents"])


def test_async_finalize_uses_run_token_to_block_abandoned_commits():
    payload = _run_node(
        r"""
        let now = 0;
        let nextFrameId = 1;
        const pendingFrames = new Map();
        const pendingFinalizers = new Map();
        const mutations = [];
        const events = [];
        const requestFrame = (callback) => { const id = nextFrameId++; pendingFrames.set(id, callback); return id; };
        const cancelFrame = (id) => pendingFrames.delete(id);
        async function advance(milliseconds) {
          now += milliseconds;
          const callbacks = Array.from(pendingFrames.values());
          pendingFrames.clear();
          callbacks.forEach((callback) => callback(now));
          await Promise.resolve();
          await Promise.resolve();
        }
        const session = runtime.createSession({
          now: () => now,
          requestFrame,
          cancelFrame,
          applyCamera: () => {},
          finalize: (_target, context) => new Promise((resolve) => {
            pendingFinalizers.set(context.runId, () => {
              const committed = context.commit(() => mutations.push(context.runId));
              resolve(committed);
            });
          }),
          validate: (_target, context) => mutations.push(`validated:${context.runId}`),
          emit: (name, detail) => events.push({ name, detail }),
        });
        const first = session.start({
          owner: "state", runId: "one", phases: [{ name: "camera", durationMs: 10 }],
        });
        await advance(10);
        const second = session.retarget({
          owner: "state", runId: "two", phases: [{ name: "camera", durationMs: 10 }],
        });
        pendingFinalizers.get("one")();
        await Promise.resolve();
        await Promise.resolve();
        await advance(10);
        pendingFinalizers.get("two")();
        await Promise.resolve();
        await Promise.resolve();
        const results = await Promise.all([first, second]);
        process.stdout.write(JSON.stringify({ mutations, events, results }));
        """
    )

    assert "one" not in payload["mutations"]
    assert "validated:one" not in payload["mutations"]
    assert payload["mutations"] == ["two", "validated:two"]
    assert payload["results"][0]["cancelled"] is True
    assert payload["results"][1]["cancelled"] is False
    assert not any(
        item["name"] == "transition-complete" and item["detail"].get("runId") == "one"
        for item in payload["events"]
    )


def test_rapid_retarget_invalidates_the_abandoned_run():
    payload = _run_node(
        r"""
        let now = 0;
        let nextFrameId = 1;
        const pending = new Map();
        const applied = [];
        const events = [];
        const requestFrame = (callback) => { const id = nextFrameId++; pending.set(id, callback); return id; };
        const cancelFrame = (id) => pending.delete(id);
        async function advance(milliseconds) {
          now += milliseconds;
          const callbacks = Array.from(pending.values());
          pending.clear();
          callbacks.forEach((callback) => callback(now));
          await Promise.resolve();
        }
        const session = runtime.createSession({
          now: () => now,
          requestFrame,
          cancelFrame,
          capture: () => ({ live: applied.at(-1) || null }),
          prepare: () => {},
          applyCamera: (progress, context) => applied.push({ progress, runId: context.runId }),
          applyAppearance: () => {},
          applyTime: () => {},
          finalize: (_target, context) => applied.push({ progress: 1, runId: context.runId, final: true }),
          validate: () => [],
          emit: (name, detail) => events.push({ name, detail }),
        });
        const first = session.start({
          owner: "state", runId: "one", target: { x: 1 }, easing: "linear",
          phases: [{ name: "camera", durationMs: 200 }],
        });
        await advance(40);
        const second = session.retarget({
          owner: "state", runId: "two", target: { x: 2 }, easing: "linear",
          phases: [{ name: "camera", durationMs: 100 }],
        });
        // A stale callback captured before retargeting must be harmless.
        for (let index = 0; index < 20 && pending.size; index += 1) await advance(20);
        const results = await Promise.all([first, second]);
        process.stdout.write(JSON.stringify({ applied, events, results }));
        """
    )

    assert not any(item.get("final") and item.get("runId") == "one" for item in payload["applied"])
    assert sum(bool(item.get("final")) and item.get("runId") == "two" for item in payload["applied"]) == 1
    completes = [item for item in payload["events"] if item["name"].endswith("complete")]
    assert not any(item["detail"].get("runId") == "one" for item in completes)
    assert any(item["detail"].get("runId") == "two" for item in completes)
    assert payload["results"][0].get("cancelled") is True


def test_session_contains_errors_and_can_restore_a_complete_snapshot_afterward():
    payload = _run_node(
        r"""
        let now = 0;
        let nextFrameId = 1;
        const pending = new Map();
        const events = [];
        const finalized = [];
        let failNextTimeApply = true;
        const requestFrame = (callback) => { const id = nextFrameId++; pending.set(id, callback); return id; };
        const cancelFrame = (id) => pending.delete(id);
        async function advance(milliseconds) {
          now += milliseconds;
          const callbacks = Array.from(pending.values());
          pending.clear();
          callbacks.forEach((callback) => callback(now));
          await Promise.resolve();
        }
        const session = runtime.createSession({
          now: () => now,
          requestFrame,
          cancelFrame,
          capture: () => ({ current_frame_value: 2.5, camera: { position: { x: 1, y: 2, z: 3 } } }),
          prepare: () => {},
          applyCamera: () => {},
          applyAppearance: () => {},
          applyTime: () => {
            if (failNextTimeApply) { failNextTimeApply = false; throw new Error("intentional-frame-error"); }
          },
          finalize: (target, context) => finalized.push({ target, runId: context.runId }),
          validate: () => [],
          emit: (name, detail) => events.push({ name, detail }),
        });
        const broken = session.start({
          owner: "action", runId: "broken", target: { current_frame_value: 7 },
          phases: [{ name: "time", durationMs: 20 }],
        }).then(() => null, (error) => String(error && error.message || error));
        await advance(20);
        const brokenError = await broken;

        const restoreSnapshot = {
          current_group: "Clusters",
          current_frame_index: 5,
          current_frame_value: 5.25,
          playback_state: { direction: 0, interval_ms: 160 },
          camera: {
            position: { x: 11, y: 12, z: 13 },
            target: { x: 1, y: 2, z: 3 },
            up: { x: 0, y: 0, z: 1 },
            view_offset: { x: 0.2, y: -0.1 },
          },
          global_controls: { camera_fov: 47, camera_view_mode: "free" },
          legend_state: { young: true },
          trace_style_state: { young: { opacity: 0.8 } },
          selected_cluster_keys: ["cluster-a"],
          lasso_selection_mask: null,
        };
        const restored = session.start({
          owner: "action", runId: "restore", target: restoreSnapshot,
          phases: [{ name: "time", durationMs: 20 }],
        });
        await advance(20);
        const restoredResult = await restored;
        process.stdout.write(JSON.stringify({ brokenError, restoredResult, finalized, events }));
        """
    )

    assert payload["brokenError"] == "intentional-frame-error"
    assert any(
        item["name"] == "transition-error" and item["detail"].get("runId") == "broken"
        for item in payload["events"]
    )
    assert payload["restoredResult"]["cancelled"] is False
    assert payload["finalized"] == [
        {
            "runId": "restore",
            "target": {
                "current_group": "Clusters",
                "current_frame_index": 5,
                "current_frame_value": 5.25,
                "playback_state": {"direction": 0, "interval_ms": 160},
                "camera": {
                    "position": {"x": 11, "y": 12, "z": 13},
                    "target": {"x": 1, "y": 2, "z": 3},
                    "up": {"x": 0, "y": 0, "z": 1},
                    "view_offset": {"x": 0.2, "y": -0.1},
                },
                "global_controls": {"camera_fov": 47, "camera_view_mode": "free"},
                "legend_state": {"young": True},
                "trace_style_state": {"young": {"opacity": 0.8}},
                "selected_cluster_keys": ["cluster-a"],
                "lasso_selection_mask": None,
            },
        }
    ]


def test_frame_coordinator_updates_only_once_for_each_animation_frame():
    payload = _run_node(
        r"""
        const calls = [];
        const errors = [];
        const coordinator = runtime.createFrameCoordinator({
          onError: (error, context) => errors.push({ message: error.message, context }),
        });
        const first = coordinator.update(10, { owner: "action", runId: "a-1" }, (context) => {
          calls.push(context);
          return true;
        });
        const duplicate = coordinator.update(10, { owner: "state", runId: "s-stale" }, () => {
          calls.push({ stale: true });
          return false;
        });
        const contained = coordinator.update(11, { owner: "state", runId: "s-1" }, () => {
          throw new Error("contained-frame-error");
        });
        process.stdout.write(JSON.stringify({ first, duplicate, contained, calls, errors, snapshot: coordinator.snapshot() }));
        """
    )

    assert payload["first"] is True
    assert payload["duplicate"] is True
    assert payload["contained"] is False
    assert len(payload["calls"]) == 1
    assert payload["calls"][0]["owner"] == "action"
    assert payload["errors"][0]["message"] == "contained-frame-error"
    assert payload["snapshot"]["owner"] == "state"
    assert payload["snapshot"]["runId"] == "s-1"
    assert payload["snapshot"]["updateCount"] == 2
