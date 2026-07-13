from __future__ import annotations


THREEJS_ACTION_RUNTIME_JS = """
      const actionsSpec = sceneSpec.actions && typeof sceneSpec.actions === "object"
        ? sceneSpec.actions
        : { enabled: false, items: [] };
      const actionDefinitions = Array.isArray(actionsSpec.items) ? actionsSpec.items : [];
      const actionButtonByKey = new Map();
      let activeActionKey = "";
      let selectedActionKey = "";
      let activeActionRun = null;
      let observedActionRun = null;
      let actionRunSerial = 0;
      let legendTransitionState = null;
      let actionHeldTraceOpacityByKey = null;
      let actionHeldAppearanceRollback = null;
      let cameraActionTrack = null;
      let timeActionTrack = null;
      let actionInterruptGuardDepth = 0;
      let actionCameraOrbitOwned = false;
      let actionCameraOrbitShouldPersist = true;
      let actionSceneRenderSerial = 0;
      let initialActionViewState = null;
      let currentActionCameraViewOffset = { x: 0.0, y: 0.0 };

      function actionTimestamp(value = null) {
        if (value !== null && value !== undefined) {
          const supplied = Number(value);
          if (Number.isFinite(supplied)) {
            return supplied;
          }
        }
        return window.performance ? window.performance.now() : Date.now();
      }

      function renderActionInterpolatedFrameValue(frameValue, options = {}) {
        actionSceneRenderSerial += 1;
        if (root && root.dataset) {
          root.dataset.actionSceneRenderCount = String(actionSceneRenderSerial);
        }
        return renderInterpolatedFrameValue(frameValue, options);
      }

      function actionEvent(name, run = activeActionRun, detail = {}) {
        const payload = Object.assign({
          owner: "action",
          runId: run ? String(run.id || "") : "",
          actionKey: run ? String(run.key || "") : "",
          actionLabel: run ? String(run.label || run.key || "") : "",
        }, detail || {});
        if (typeof ovizStateEvent === "function") {
          ovizStateEvent(name, payload);
        } else {
          root.dispatchEvent(new CustomEvent(name, { detail: Object.assign({ rootId: root.id }, payload) }));
        }
      }

      function actionDebugState(run, record = null, detail = {}) {
        if (!root || !root.dataset) {
          return;
        }
        root.dataset.actionOwner = run ? "action" : "";
        root.dataset.actionRunId = run ? String(run.id || "") : "";
        root.dataset.actionKey = run ? String(run.key || "") : "";
        root.dataset.actionStep = record ? String(record.index) : "";
        root.dataset.actionPhase = String(detail.phase || (record && record.domain) || "");
        if (detail.progress !== undefined) {
          root.dataset.actionProgress = String(clampRange(Number(detail.progress) || 0.0, 0.0, 1.0));
        }
        if (detail.effectiveDurationMs !== undefined) {
          root.dataset.actionEffectiveDurationMs = String(Math.max(Number(detail.effectiveDurationMs) || 0.0, 0.0));
        }
        if (detail.status) {
          root.dataset.actionStatus = String(detail.status);
        }
      }

      function actionDomainsForStep(step) {
        const type = String(step && step.type || "");
        if (type === "state") {
          return ["camera", "appearance", "selection", "time"];
        }
        if (type === "legend_group") {
          return ["appearance"];
        }
        if (type === "camera") {
          return ["camera"];
        }
        if (type === "time") {
          return ["time"];
        }
        return [];
      }

      function actionDomainsForStepAtStart(step) {
        if (
          step
          && String(step.type || "") === "camera"
          && String(cameraViewMode || "free") === "earth"
        ) {
          // Sky camera Actions are implemented by the State controller so the
          // registered Aladin path remains the sole spatial writer.  Treat that
          // routed step as an exclusive State session; otherwise a nominally
          // disjoint with_previous time step can run behind it and jump later.
          return ["camera", "appearance", "selection", "time"];
        }
        return actionDomainsForStep(step);
      }

      function actionRecordOwnsDomains(run, record) {
        if (!run || !record || !(run.domainOwners instanceof Map)) {
          return false;
        }
        return record.domains.every((domain) => run.domainOwners.get(domain) === record.index);
      }

      function actionDomainsAvailable(run, record) {
        if (!run || !record || !(run.domainOwners instanceof Map)) {
          return false;
        }
        return record.domains.every((domain) => !run.domainOwners.has(domain));
      }

      function claimActionDomains(run, record) {
        if (!actionDomainsAvailable(run, record)) {
          return false;
        }
        record.domains.forEach((domain) => run.domainOwners.set(domain, record.index));
        return true;
      }

      function releaseActionDomains(run, record) {
        if (!run || !record || !(run.domainOwners instanceof Map)) {
          return;
        }
        record.domains.forEach((domain) => {
          if (run.domainOwners.get(domain) === record.index) {
            run.domainOwners.delete(domain);
          }
        });
      }

      function actionInterruptsMuted() {
        return actionInterruptGuardDepth > 0;
      }

      function withActionGuard(callback) {
        actionInterruptGuardDepth += 1;
        try {
          return callback();
        } finally {
          actionInterruptGuardDepth = Math.max(actionInterruptGuardDepth - 1, 0);
        }
      }

      function actionEasingValue(name, t) {
        const clamped = clampRange(t, 0.0, 1.0);
        const easingName = String(name || "easeInOutCubic");
        if (easingName === "linear") {
          return clamped;
        }
        if (easingName === "easeOutCubic") {
          return 1.0 - Math.pow(1.0 - clamped, 3.0);
        }
        if (easingName === "easeInCubic") {
          return clamped * clamped * clamped;
        }
        return clamped < 0.5
          ? 4.0 * clamped * clamped * clamped
          : 1.0 - Math.pow(-2.0 * clamped + 2.0, 3.0) / 2.0;
      }

      function cloneLegendStateMap(source) {
        const next = {};
        if (!source || typeof source !== "object") {
          return next;
        }
        Object.entries(source).forEach(([key, value]) => {
          next[String(key)] = Boolean(value);
        });
        return next;
      }

      function legendStateForGroup(groupName) {
        const nextState = {};
        const defaults = groupDefaults(groupName);
        legendItems.forEach((item) => {
          const mode = defaults[item.key];
          if (mode === true) {
            nextState[item.key] = true;
          } else if (mode === "legendonly") {
            nextState[item.key] = false;
          }
        });
        return nextState;
      }

      function traceVisibleForGroupState(trace, groupName, stateOverride) {
        if (isGalacticReferenceTrace(trace) && !galacticReferenceMotionVisible()) {
          return false;
        }
        if (isNearbyRegionLabelTrace(trace) && !nearbyRegionLabelsVisible) {
          return false;
        }
        const defaults = groupDefaults(groupName);
        const mode = defaults[trace.key];
        if (mode === false || mode === undefined) {
          return false;
        }
        if (trace.showlegend) {
          const stateMap = stateOverride && typeof stateOverride === "object"
            ? stateOverride
            : legendStateForGroup(groupName);
          if (Object.prototype.hasOwnProperty.call(stateMap, trace.key)) {
            return Boolean(stateMap[trace.key]);
          }
          return mode === true;
        }
        return mode === true;
      }

      function actionTraceVisibilityState(trace) {
        if (!trace || !trace.key) {
          return null;
        }
        if (!legendTransitionState) {
          if (
            actionHeldTraceOpacityByKey
            && Object.prototype.hasOwnProperty.call(actionHeldTraceOpacityByKey, trace.key)
          ) {
            const opacity = clampRange(Number(actionHeldTraceOpacityByKey[trace.key]) || 0.0, 0.0, 1.0);
            return { visible: opacity > 0.0001, opacity };
          }
          return null;
        }
        const fromVisible = traceVisibleForGroupState(
          trace,
          legendTransitionState.fromGroup,
          legendTransitionState.fromLegendState,
        );
        const toVisible = traceVisibleForGroupState(
          trace,
          legendTransitionState.toGroup,
          legendTransitionState.toLegendState,
        );
        const fromOpacity = legendTransitionState.fromOpacityByKey
          && Object.prototype.hasOwnProperty.call(legendTransitionState.fromOpacityByKey, trace.key)
          ? clampRange(Number(legendTransitionState.fromOpacityByKey[trace.key]) || 0.0, 0.0, 1.0)
          : (fromVisible ? 1.0 : 0.0);
        const toOpacity = toVisible ? 1.0 : 0.0;
        if (fromOpacity <= 0.0001 && toOpacity <= 0.0001) {
          return { visible: false, opacity: 0.0 };
        }
        const progress = clampRange(Number(legendTransitionState.progress) || 0.0, 0.0, 1.0);
        const opacity = fromOpacity + (toOpacity - fromOpacity) * progress;
        return { visible: opacity > 0.0001, opacity };
      }

      function captureActionTraceOpacityMap() {
        const result = {};
        const candidates = new Map();
        const remember = (trace) => {
          if (trace && trace.key && !candidates.has(trace.key)) candidates.set(trace.key, trace);
        };
        (currentFrame() && currentFrame().traces || []).forEach(remember);
        if (
          typeof ovizStateTransition !== "undefined"
          && ovizStateTransition
          && typeof ovizTraceCandidatesForSnapshots === "function"
        ) {
          ovizTraceCandidatesForSnapshots(
            ovizStateTransition.fromSnapshot,
            ovizStateTransition.targetSnapshot,
          ).forEach(remember);
        }
        candidates.forEach((trace) => {
          const state = actionTraceVisibilityState(trace);
          const stateTransition = stateTraceVisibilityState(trace);
          result[trace.key] = state
            ? clampRange(Number(state.opacity) || 0.0, 0.0, 1.0)
            : (
              stateTransition
                ? clampRange(Number(stateTransition.opacity) || 0.0, 0.0, 1.0)
                : (traceVisible(trace) ? 1.0 : 0.0)
            );
        });
        return result;
      }

      function traceVisibilityOpacityMultiplier(trace) {
        const state = actionTraceVisibilityState(trace);
        const stateTransition = stateTraceVisibilityState(trace);
        const actionOpacity = state ? clampRange(state.opacity, 0.0, 1.0) : 1.0;
        const transitionOpacity = stateTransition ? clampRange(stateTransition.opacity, 0.0, 1.0) : 1.0;
        return actionOpacity * transitionOpacity;
      }

      function renderActionBar() {
        if (!actionBarEl) {
          return;
        }
        actionBarEl.innerHTML = "";
        actionButtonByKey.clear();
        if (!minimalModeEnabled) {
          actionBarEl.dataset.visible = "false";
          actionBarEl.style.display = "none";
          return;
        }
        const visibleActions = actionDefinitions.filter((action) => action && Array.isArray(action.steps) && action.steps.length);
        actionBarEl.dataset.visible = visibleActions.length ? "true" : "false";
        if (!visibleActions.length) {
          actionBarEl.style.display = "none";
          return;
        }
        actionBarEl.style.display = "flex";
        visibleActions.forEach((action) => {
          const button = document.createElement("button");
          button.type = "button";
          button.className = "oviz-three-action-button";
          button.dataset.actionKey = String(action.key || "");
          button.textContent = String(action.label || action.key || "Action");
          if (action.description) {
            button.title = String(action.description);
            button.setAttribute("aria-label", String(action.description));
          }
          button.addEventListener("click", () => {
            const actionKey = String(action.key || "");
            if (selectedActionKey && selectedActionKey === actionKey) {
              startRestoreInitialView();
            } else {
              startViewerAction(actionKey);
            }
            focusViewer();
          });
          actionButtonByKey.set(String(action.key || ""), button);
          actionBarEl.appendChild(button);
        });
      }

      function syncActionButtons() {
        actionButtonByKey.forEach((button, key) => {
          const active = key === selectedActionKey;
          button.dataset.active = active ? "true" : "false";
          button.setAttribute("aria-pressed", active ? "true" : "false");
        });
      }

      function normalizeActionViewOffset(value) {
        if (!value || typeof value !== "object") {
          return { x: 0.0, y: 0.0 };
        }
        const x = Number(value.x);
        const y = Number(value.y);
        return {
          x: Number.isFinite(x) ? x : 0.0,
          y: Number.isFinite(y) ? y : 0.0,
        };
      }

      function applyActionCameraViewOffset(viewOffset) {
        const nextOffset = normalizeActionViewOffset(viewOffset);
        currentActionCameraViewOffset = nextOffset;
        if (Math.abs(nextOffset.x) <= 1e-6 && Math.abs(nextOffset.y) <= 1e-6) {
          if (typeof camera.clearViewOffset === "function") {
            camera.clearViewOffset();
          }
          camera.updateProjectionMatrix();
          return;
        }
        const width = Math.max(root.clientWidth || sceneSpec.width || 1, 1);
        const height = Math.max(root.clientHeight || sceneSpec.height || 1, 1);
        camera.setViewOffset(
          width,
          height,
          Math.round(width * nextOffset.x),
          Math.round(height * nextOffset.y),
          width,
          height,
        );
        camera.updateProjectionMatrix();
      }

      function captureCurrentActionViewState() {
        return captureRuntimeState();
      }

      function beginLegendTransition(targetGroup, targetLegendState, stepIndex, now, durationMs, easing) {
        const timestamp = actionTimestamp(now);
        const runId = activeActionRun ? activeActionRun.id : null;
        const previousGroup = currentGroup;
        const previousLegendState = cloneLegendStateMap(legendState);
        const nextGroup = String(targetGroup || currentGroup);
        const nextLegendState = cloneLegendStateMap(targetLegendState);
        const fromOpacityByKey = actionHeldTraceOpacityByKey;
        actionHeldTraceOpacityByKey = null;
        const nextStateMatchesCurrent = !fromOpacityByKey && nextGroup === currentGroup
          && JSON.stringify(nextLegendState) === JSON.stringify(previousLegendState);
        currentGroup = nextGroup;
        if (groupSelectEl) {
          groupSelectEl.value = currentGroup;
        }
        syncLegendGroupChooser();
        renderLegend();
        const nextDurationMs = Math.max(Number(durationMs) || 0, 0);
        const stepRecord = actionStepRecord(stepIndex, runId);
        if (stepRecord) {
          stepRecord.effectiveDurationMs = nextDurationMs;
        }
        if (nextDurationMs <= 0 || nextStateMatchesCurrent) {
          legendState = nextLegendState;
          renderLegend();
          renderActionInterpolatedFrameValue(displayedFrameValue, {
            updateWidgets: false,
            preserveCamera: true,
          });
          markActionStepFinished(stepIndex, timestamp, runId);
          return;
        }
        legendTransitionState = {
          stepIndex,
          runId,
          fromGroup: previousGroup,
          fromLegendState: previousLegendState,
          toGroup: nextGroup,
          toLegendState: nextLegendState,
          startTime: timestamp,
          durationMs: nextDurationMs,
          easing: String(easing || "easeInOutCubic"),
          progress: 0.0,
          fromOpacityByKey,
        };
      }

      function targetFrameIndexForTimeAction(step) {
        const currentIndex = clampFrameIndex(currentFrameIndex);
        if (step.stop_after_frames !== undefined) {
          const frames = Math.max(Number(step.stop_after_frames) || 0, 0);
          const direction = String(step.direction) === "backward" ? -1 : 1;
          return clampFrameIndex(currentIndex + direction * frames);
        }
        if (step.stop_at_time_myr !== undefined) {
          const targetTime = Number(step.stop_at_time_myr);
          let bestIndex = currentIndex;
          let bestDistance = Infinity;
          frameSpecs.forEach((frame, index) => {
            const frameTime = Number(frame && frame.time);
            if (!Number.isFinite(frameTime)) {
              return;
            }
            const direction = String(step.direction) === "backward" ? -1 : 1;
            if (direction < 0 && frameTime > frameTimeForValue(currentFrameIndex)) {
              return;
            }
            if (direction > 0 && frameTime < frameTimeForValue(currentFrameIndex)) {
              return;
            }
            const distance = Math.abs(frameTime - targetTime);
            if (distance < bestDistance) {
              bestDistance = distance;
              bestIndex = index;
            }
          });
          return bestIndex;
        }
        return String(step.direction) === "backward" ? 0 : Math.max(frameSpecs.length - 1, 0);
      }

      function buildScheduledActionSteps(action) {
        return (action.steps || []).map((step, index) => {
          const domains = actionDomainsForStep(step);
          return {
            index,
            step,
            domains,
            domain: domains.join("+") || String(step.type || "action"),
            started: false,
            finished: false,
            startedAt: null,
            finishedAt: null,
            effectiveDurationMs: 0.0,
            lastProgressEventAt: -Infinity,
            lastProgress: 0.0,
          };
        });
      }

      function actionDefinitionByKey(actionKey) {
        if (!actionKey) {
          return null;
        }
        return actionDefinitions.find((item) => String(item.key) === String(actionKey)) || null;
      }

      function restoreTimeIntervalMsForAction(action) {
        const fallbackIntervalMs = Math.max(Number(playbackIntervalMs) || 160, 1);
        if (!action || !Array.isArray(action.steps)) {
          return fallbackIntervalMs;
        }
        const timeStep = action.steps.find((step) => step && step.type === "time");
        if (!timeStep) {
          return fallbackIntervalMs;
        }
        return Math.max(Number(timeStep.interval_ms) || fallbackIntervalMs, 1);
      }

      function actionStepRecord(stepIndex, runId = null) {
        if (
          !activeActionRun
          || activeActionRun.status !== "running"
          || !Array.isArray(activeActionRun.steps)
          || (runId !== null && String(activeActionRun.id) !== String(runId))
        ) {
          return null;
        }
        return activeActionRun.steps.find((step) => step.index === stepIndex) || null;
      }

      function emitActionProgress(record, now, rawProgress, easedProgress, phase = null) {
        const run = activeActionRun;
        if (!run || !record || String(run.id) !== String(record.runId || run.id)) {
          return;
        }
        const raw = clampRange(Number(rawProgress) || 0.0, 0.0, 1.0);
        const eased = clampRange(Number(easedProgress) || 0.0, 0.0, 1.0);
        record.lastProgress = raw;
        const timestamp = actionTimestamp(now);
        if (raw < 1.0 && timestamp - Number(record.lastProgressEventAt) < 100.0) {
          return;
        }
        record.lastProgressEventAt = timestamp;
        const detail = {
          step: record.index,
          stepType: String(record.step && record.step.type || ""),
          phase: String(phase || record.domain || ""),
          progress: raw,
          easedProgress: eased,
          effectiveDurationMs: Math.max(Number(record.effectiveDurationMs) || 0.0, 0.0),
        };
        actionDebugState(run, record, Object.assign({ status: "running" }, detail));
        actionEvent("action-progress", run, detail);
      }

      function markActionStepFinished(stepIndex, now = null, runId = null) {
        const record = actionStepRecord(stepIndex, runId);
        if (record && !record.finished) {
          const timestamp = actionTimestamp(now);
          emitActionProgress(record, timestamp, 1.0, 1.0, record.domain);
          record.finished = true;
          record.finishedAt = timestamp;
          releaseActionDomains(activeActionRun, record);
        }
      }

      function actionStepReadyAt(run, record) {
        const delayMs = Math.max(Number(record.step && record.step.delay_ms) || 0.0, 0.0);
        if (record.index <= 0) {
          return Number(run.startedAt) + delayMs;
        }
        const previous = run.steps[record.index - 1];
        if (!previous) {
          return Number(run.startedAt) + delayMs;
        }
        if (String(record.step.start || "after_previous") === "with_previous") {
          return previous.startedAt === null ? null : Number(previous.startedAt) + delayMs;
        }
        return previous.finishedAt === null ? null : Number(previous.finishedAt) + delayMs;
      }

      function finishLegendTransition(commitTarget = true, now = null, options = {}) {
        if (!legendTransitionState) {
          return false;
        }
        const stepIndex = legendTransitionState.stepIndex;
        const runId = legendTransitionState.runId;
        const nextGroup = legendTransitionState.toGroup;
        const nextLegendState = cloneLegendStateMap(legendTransitionState.toLegendState);
        legendTransitionState = null;
        actionHeldTraceOpacityByKey = null;
        if (commitTarget) {
          currentGroup = nextGroup;
          legendState = nextLegendState;
          if (groupSelectEl) {
            groupSelectEl.value = currentGroup;
          }
          syncLegendGroupChooser();
        }
        renderLegend();
        if (options.render !== false) {
          renderActionInterpolatedFrameValue(displayedFrameValue, {
            updateWidgets: false,
            preserveCamera: true,
          });
        }
        markActionStepFinished(stepIndex, now, runId);
        return true;
      }

      function stopCameraActionTrack(options = {}) {
        if (!cameraActionTrack) {
          return false;
        }
        const track = cameraActionTrack;
        cameraActionTrack = null;
        if (options.complete !== false) {
          camera.position.copy(track.toPosition);
          controls.target.copy(track.toTarget);
          camera.up.copy(track.toUp);
          camera.fov = Number(track.toFov);
          camera.updateProjectionMatrix();
          applyActionCameraViewOffset(track.toViewOffset);
          controls.update();
        }
        if (track.orbitEnabled && options.activateOrbit !== false) {
          cameraAutoOrbitSpeedMultiplier = Math.max(Number(track.orbitSpeedMultiplier) || 1.0, 0.05);
          cameraAutoOrbitDirection = Number(track.orbitDirection) < 0.0 ? -1.0 : 1.0;
          actionCameraOrbitOwned = track.orbitOwned !== false;
          setCameraAutoOrbitEnabled(true);
          actionCameraOrbitShouldPersist = track.orbitPersist !== false;
        }
        if (options.markFinished !== false) {
          markActionStepFinished(track.stepIndex, options.now, track.runId);
        }
        return true;
      }

      function stopTimeActionTrack(options = {}) {
        if (!timeActionTrack) {
          return false;
        }
        const track = timeActionTrack;
        timeActionTrack = null;
        if (options.complete === true) {
          displayedFrameValue = clampFrameValue(track.targetFrameValue);
          currentFrameIndex = clampFrameIndex(displayedFrameValue);
          renderActionInterpolatedFrameValue(displayedFrameValue, {
            updateWidgets: false,
            preserveCamera: true,
          });
        }
        playbackDirection = 0;
        lastPlaybackAdvanceTimestamp = null;
        updateTimelineMotionOpacity();
        updatePlaybackButtons();
        if (options.markFinished !== false) {
          markActionStepFinished(track.stepIndex, options.now, track.runId);
        }
        return true;
      }

      function finishActionHeldAppearanceRollback() {
        const held = actionHeldAppearanceRollback;
        if (held) {
          held.currentAppearanceProgress = 0.0;
          setMilkyWayModelOpacityScale(held.toMilkyWayOpacity);
          setSkyDomeViewOpacityScale(held.toSkyOpacity, { force: false });
          actionHeldAppearanceRollback = null;
        }
        if (typeof ovizHeldSelectionTransition !== "undefined" && ovizHeldSelectionTransition) {
          ovizHeldSelectionTransition.progress = 1.0;
          if (typeof releaseOvizHeldSelectionTransition === "function") {
            releaseOvizHeldSelectionTransition();
          } else {
            ovizHeldSelectionTransition = null;
          }
        }
        return Boolean(held);
      }

      function updateActionHeldAppearanceRollback(now) {
        const held = actionHeldAppearanceRollback;
        if (!held) return false;
        const durationMs = Math.max(Number(held.durationMs) || 240.0, 1.0);
        const rawProgress = clampRange(
          (actionTimestamp(now) - Number(held.startedAt)) / durationMs,
          0.0,
          1.0,
        );
        const progress = actionEasingValue("easeInOutCubic", rawProgress);
        held.currentAppearanceProgress = clampRange(
          Number(held.frozenAppearanceProgress) * (1.0 - progress),
          0.0,
          1.0,
        );
        if (typeof ovizHeldSelectionTransition !== "undefined" && ovizHeldSelectionTransition) {
          ovizHeldSelectionTransition.progress = progress;
        }
        setMilkyWayModelOpacityScale(
          Number(held.fromMilkyWayOpacity)
          + (Number(held.toMilkyWayOpacity) - Number(held.fromMilkyWayOpacity)) * progress
        );
        setSkyDomeViewOpacityScale(
          Number(held.fromSkyOpacity)
          + (Number(held.toSkyOpacity) - Number(held.fromSkyOpacity)) * progress,
          { force: false },
        );
        if (rawProgress >= 1.0) {
          finishActionHeldAppearanceRollback();
        }
        return true;
      }

      function cancelActionRun(reason = "cancelled", options = {}) {
        const run = activeActionRun;
        const timestamp = actionTimestamp(options.now);
        if (run) {
          run.status = "cancelled";
          run.cancelReason = String(reason || "cancelled");
        }
        if (
          run
          && options.cancelStateTransition !== false
          && typeof ovizStateTransition !== "undefined"
          && ovizStateTransition
          && typeof ovizCancelStateTransitionWithoutSnap === "function"
        ) {
          try {
            ovizCancelStateTransitionWithoutSnap(`action-${String(reason || "cancelled")}`, {
              restorePresentation: options.restorePresentation !== false,
            });
          } catch (err) {
            console.error("Oviz Action could not cancel its State transition", err);
          }
        }
        if (legendTransitionState) {
          try {
            if (options.preserveLegendTransitionFrame === true) {
              actionHeldTraceOpacityByKey = captureActionTraceOpacityMap();
            } else if (options.commitLegendTransition !== false) {
              currentGroup = legendTransitionState.toGroup;
              legendState = cloneLegendStateMap(legendTransitionState.toLegendState);
              if (groupSelectEl) groupSelectEl.value = currentGroup;
              syncLegendGroupChooser();
              renderLegend();
            }
          } catch (err) {
            actionHeldTraceOpacityByKey = null;
            console.error("Oviz Action legend cleanup failed", err);
          }
          legendTransitionState = null;
        }
        cameraActionTrack = null;
        timeActionTrack = null;
        if (options.preserveHeldAppearance !== true) {
          finishActionHeldAppearanceRollback();
        }
        playbackDirection = 0;
        lastPlaybackAdvanceTimestamp = null;
        try {
          updateTimelineMotionOpacity();
          updatePlaybackButtons();
        } catch (err) {
          console.error("Oviz Action timeline cleanup failed", err);
        }
        if (options.disableOrbit && actionCameraOrbitOwned) {
          cameraAutoOrbitSpeedMultiplier = 1.0;
          try {
            setCameraAutoOrbitEnabled(false);
          } catch (err) {
            console.error("Oviz Action orbit cleanup failed", err);
          }
        }
        actionCameraOrbitOwned = false;
        actionCameraOrbitShouldPersist = true;
        activeActionRun = null;
        activeActionKey = "";
        if (options.clearSelectedAction) {
          selectedActionKey = "";
        }
        try {
          syncActionButtons();
        } catch (err) {
          console.error("Oviz Action button cleanup failed", err);
        }
        if (run && options.emitEvent !== false) {
          actionDebugState(run, null, { status: "cancelled", phase: "cancelled", progress: 0.0 });
          actionEvent("action-cancel", run, {
            reason: String(reason || "cancelled"),
            elapsedMs: Math.max(0.0, timestamp - Number(run.startedAt || timestamp)),
          });
        }
        if (observedActionRun === run) {
          observedActionRun = null;
        }
        return Boolean(run);
      }

      function clearActiveActionState(options = {}) {
        return cancelActionRun(String(options.reason || "cleared"), options);
      }

      function failActionRun(err, record = null, now = null) {
        const run = activeActionRun;
        if (!run) {
          return false;
        }
        const timestamp = actionTimestamp(now);
        const message = String(err && err.message || err || "Unknown Action error");
        try {
          cancelActionRun("error", {
            now: timestamp,
            preserveLegendTransitionFrame: true,
            disableOrbit: true,
            cancelStateTransition: true,
            emitEvent: false,
          });
        } catch (cleanupErr) {
          activeActionRun = null;
          if (observedActionRun === run) observedActionRun = null;
          activeActionKey = "";
          legendTransitionState = null;
          cameraActionTrack = null;
          timeActionTrack = null;
          playbackDirection = 0;
          lastPlaybackAdvanceTimestamp = null;
          actionCameraOrbitOwned = false;
          actionCameraOrbitShouldPersist = true;
          console.error("Oviz Action cleanup failed", cleanupErr);
        }
        if (root && root.dataset) {
          root.dataset.actionStatus = "error";
          root.dataset.actionError = message;
        }
        try {
          actionEvent("action-error", run, {
            step: record ? record.index : null,
            stepType: record && record.step ? String(record.step.type || "") : "",
            error: message,
            elapsedMs: Math.max(0.0, timestamp - Number(run.startedAt || timestamp)),
          });
        } catch (eventErr) {
          console.error("Oviz Action error event failed", eventErr);
        }
        console.error("Oviz Action failed", err);
        return true;
      }

      function completeActionRun(now = null) {
        const run = activeActionRun;
        if (!run || run.status !== "running") {
          return false;
        }
        const timestamp = actionTimestamp(now);
        finishActionHeldAppearanceRollback();
        run.status = "complete";
        activeActionRun = null;
        activeActionKey = "";
        actionDebugState(run, null, { status: "complete", phase: "complete", progress: 1.0 });
        actionEvent("action-complete", run, {
          elapsedMs: Math.max(0.0, timestamp - Number(run.startedAt || timestamp)),
        });
        if (observedActionRun === run) {
          observedActionRun = null;
        }
        syncActionButtons();
        return true;
      }

      function interruptActionRun(reason = "user", options = {}) {
        if (!activeActionRun || actionInterruptsMuted()) {
          return false;
        }
        cancelActionRun(reason, {
          commitLegendTransition: true,
          completeCamera: false,
          preserveLegendTransitionFrame: true,
          disableOrbit: options.disableOrbit !== false && (reason === "camera" || reason === "keyboard" || reason === "user"),
          clearSelectedAction: options.clearSelectedAction === true,
        });
        return true;
      }

      function handleManualCameraInteractionStart() {
        if (actionInterruptsMuted()) {
          return false;
        }
        if (activeActionRun) {
          return interruptActionRun("camera", { disableOrbit: true });
        }
        if (minimalModeEnabled && cameraAutoOrbitEnabled) {
          cameraAutoOrbitSpeedMultiplier = 1.0;
          setCameraAutoOrbitEnabled(false);
          actionCameraOrbitOwned = false;
          actionCameraOrbitShouldPersist = true;
          return true;
        }
        if (cameraAutoOrbitEnabled && typeof galacticSimpleDefaultOrbitActive === "function" && galacticSimpleDefaultOrbitActive()) {
          cameraAutoOrbitSpeedMultiplier = 1.0;
          setCameraAutoOrbitEnabled(false);
          actionCameraOrbitOwned = false;
          actionCameraOrbitShouldPersist = true;
          return true;
        }
        if (actionCameraOrbitOwned && cameraAutoOrbitEnabled) {
          cameraAutoOrbitSpeedMultiplier = 1.0;
          setCameraAutoOrbitEnabled(false);
          actionCameraOrbitOwned = false;
          actionCameraOrbitShouldPersist = true;
          return true;
        }
        return false;
      }

      function currentFrameGeometryPointsForTrace(trace) {
        const points = [];
        if (!trace || typeof trace !== "object") {
          return points;
        }
        (trace.points || []).forEach((point) => {
          const x = Number(point && point.x);
          const y = Number(point && point.y);
          const z = Number(point && point.z);
          if (Number.isFinite(x) && Number.isFinite(y) && Number.isFinite(z)) {
            points.push(new THREE.Vector3(x, y, z));
          }
        });
        (trace.labels || []).forEach((label) => {
          const x = Number(label && label.x);
          const y = Number(label && label.y);
          const z = Number(label && label.z);
          if (Number.isFinite(x) && Number.isFinite(y) && Number.isFinite(z)) {
            points.push(new THREE.Vector3(x, y, z));
          }
        });
        (trace.segments || []).forEach((segment) => {
          if (!Array.isArray(segment) || segment.length < 6) {
            return;
          }
          const coords = segment.map((value) => Number(value));
          if (coords.every((value) => Number.isFinite(value))) {
            points.push(new THREE.Vector3(coords[0], coords[1], coords[2]));
            points.push(new THREE.Vector3(coords[3], coords[4], coords[5]));
          }
        });
        return points;
      }

      function actionTargetPointsForSpec(targetSpec) {
        const frame = currentFrame();
        if (!frame || !targetSpec || typeof targetSpec !== "object") {
          return [];
        }
        if (String(targetSpec.kind) === "point") {
          return [new THREE.Vector3(Number(targetSpec.x) || 0.0, Number(targetSpec.y) || 0.0, Number(targetSpec.z) || 0.0)];
        }
        const traces = Array.isArray(frame.traces) ? frame.traces : [];
        if (String(targetSpec.kind) === "trace") {
          const match = traces.find((trace) => String(trace.key) === String(targetSpec.key));
          return currentFrameGeometryPointsForTrace(match);
        }
        if (String(targetSpec.kind) === "group") {
          const groupName = String(targetSpec.name || "");
          const groupLegendState = legendStateForGroup(groupName);
          const groupPoints = [];
          traces.forEach((trace) => {
            if (!traceVisibleForGroupState(trace, groupName, groupLegendState)) {
              return;
            }
            currentFrameGeometryPointsForTrace(trace).forEach((point) => groupPoints.push(point));
          });
          return groupPoints;
        }
        return [];
      }

      function actionAnchorPointForSpec(targetSpec, fallbackPoint) {
        const targetPoints = actionTargetPointsForSpec(targetSpec);
        if (!targetPoints.length) {
          return fallbackPoint ? fallbackPoint.clone() : null;
        }
        const bounds = new THREE.Box3();
        targetPoints.forEach((point) => bounds.expandByPoint(point));
        return bounds.getCenter(new THREE.Vector3());
      }

      function resolveCameraActionState(step) {
        const targetPoints = actionTargetPointsForSpec(step.target);
        if (!targetPoints.length) {
          return null;
        }
        const bounds = new THREE.Box3();
        targetPoints.forEach((point) => bounds.expandByPoint(point));
        const center = bounds.getCenter(new THREE.Vector3());
        const anchorCenter = actionAnchorPointForSpec(step.anchor_target, center) || center.clone();
        const size = bounds.getSize(new THREE.Vector3());
        const maxExtent = Math.max(size.x, size.y, size.z, sceneSpec.max_span * 0.02, 24.0);
        const fitPadding = Math.max(Number(step.fit_padding) || 1.15, 0.1);
        const currentViewDirection = camera.position.clone().sub(controls.target);
        if (currentViewDirection.lengthSq() <= 1e-10) {
          currentViewDirection.set(0.0, -1.0, 1.0);
        }
        currentViewDirection.normalize();

        const overrides = step.camera && typeof step.camera === "object" ? step.camera : {};
        const nextTarget = overrides.target
          ? new THREE.Vector3(Number(overrides.target.x) || 0.0, Number(overrides.target.y) || 0.0, Number(overrides.target.z) || 0.0)
          : anchorCenter.clone();
        const nextFov = clampRange(Number(overrides.fov) || camera.fov, 0.05, 120.0);
        const baseFitDistance = Math.max(
          (maxExtent * fitPadding) / Math.tan((nextFov * Math.PI / 180.0) * 0.5),
          sceneSpec.max_span * 0.02,
        );
        const distanceScale = Math.max(Number(step.distance_scale) || 1.0, 0.01);
        const fitDistance = Math.max(baseFitDistance * distanceScale, 12.0);
        const nextPosition = overrides.position
          ? new THREE.Vector3(Number(overrides.position.x) || 0.0, Number(overrides.position.y) || 0.0, Number(overrides.position.z) || 0.0)
          : nextTarget.clone().add(currentViewDirection.clone().multiplyScalar(fitDistance));
        const nextUp = overrides.up
          ? new THREE.Vector3(Number(overrides.up.x) || 0.0, Number(overrides.up.y) || 0.0, Number(overrides.up.z) || 1.0)
          : camera.up.clone();
        const nextViewOffset = overrides.view_offset
          ? normalizeActionViewOffset(overrides.view_offset)
          : normalizeActionViewOffset(currentActionCameraViewOffset);
        return {
          toPosition: nextPosition,
          toTarget: nextTarget,
          toUp: nextUp,
          toFov: nextFov,
          toViewOffset: nextViewOffset,
        };
      }

      function startLegendGroupAction(step, stepIndex, now) {
        const nextGroup = String(step.group || currentGroup);
        const nextLegendState = legendStateForGroup(nextGroup);
        beginLegendTransition(
          nextGroup,
          nextLegendState,
          stepIndex,
          now,
          Number(step.duration_ms) || 0,
          step.easing,
        );
      }

      function startCameraAction(step, stepIndex, now) {
        const timestamp = actionTimestamp(now);
        const ownerRun = activeActionRun;
        const runId = ownerRun ? ownerRun.id : null;
        const resolved = resolveCameraActionState(step);
        if (!resolved) {
          throw new Error(`Camera Action target could not be resolved for step ${stepIndex}.`);
        }
        if (cameraViewMode === "earth" && ovizStateControllerReady) {
          const destination = captureRuntimeState();
          destination.camera = {
            position: { x: resolved.toPosition.x, y: resolved.toPosition.y, z: resolved.toPosition.z },
            target: { x: resolved.toTarget.x, y: resolved.toTarget.y, z: resolved.toTarget.z },
            up: { x: resolved.toUp.x, y: resolved.toUp.y, z: resolved.toUp.z },
            view_offset: normalizeActionViewOffset(resolved.toViewOffset),
          };
          destination.global_controls.camera_view_mode = "free";
          destination.global_controls.camera_fov = Number(resolved.toFov);
          destination.global_controls.camera_auto_orbit_enabled = Boolean(step.orbit && step.orbit.enabled);
          destination.global_controls.camera_auto_orbit_direction = step.orbit && Number(step.orbit.direction) < 0.0 ? -1.0 : 1.0;
          const activeStateIndex = ovizActiveStateId === null
            ? -1
            : ovizStatesProject.items.findIndex((item) => item.id === ovizActiveStateId);
          const statePromise = ovizBeginStateTransition({
            id: ovizActiveStateId,
            index: activeStateIndex,
            name: "Camera action",
            snapshot: destination,
            transition: {
              duration_ms: Math.max(Number(step.duration_ms) || 0, 0),
              easing: String(step.easing || "easeInOutCubic"),
            },
          }, { preserveActionRun: true });
          const record = actionStepRecord(stepIndex, runId);
          if (record && ovizStateTransition && ovizStateTransition.phasePlan) {
            record.effectiveDurationMs = Math.max(Number(ovizStateTransition.phasePlan.effectiveDurationMs) || 0.0, 0.0);
          }
          statePromise.then((result) => {
            if (!activeActionRun || String(activeActionRun.id) !== String(runId)) return;
            if (result && result.cancelled) {
              cancelActionRun(String(result.reason || "state-transition-cancelled"), {
                preserveLegendTransitionFrame: true,
                disableOrbit: true,
                cancelStateTransition: false,
              });
              return;
            }
            if (step.orbit && step.orbit.enabled && step.orbit.persist === false) {
              setCameraAutoOrbitEnabled(false);
            }
            markActionStepFinished(stepIndex, null, runId);
          }).catch((err) => {
            if (activeActionRun && String(activeActionRun.id) === String(runId)) {
              failActionRun(err, actionStepRecord(stepIndex, runId));
            }
          });
          return;
        }
        setCameraAutoOrbitEnabled(false);
        actionCameraOrbitOwned = false;
        actionCameraOrbitShouldPersist = true;
        cameraAutoOrbitSpeedMultiplier = 1.0;
        const durationMs = Math.max(Number(step.duration_ms) || 0, 0);
        const record = actionStepRecord(stepIndex, runId);
        if (record) record.effectiveDurationMs = durationMs;
        if (durationMs <= 0) {
          camera.position.copy(resolved.toPosition);
          controls.target.copy(resolved.toTarget);
          camera.up.copy(resolved.toUp);
          camera.fov = Number(resolved.toFov);
          camera.updateProjectionMatrix();
          applyActionCameraViewOffset(resolved.toViewOffset);
          controls.update();
          if (step.orbit && step.orbit.enabled) {
            cameraAutoOrbitSpeedMultiplier = Math.max(Number(step.orbit.speed_multiplier) || 1.0, 0.05);
            cameraAutoOrbitDirection = Number(step.orbit.direction) < 0.0 ? -1.0 : 1.0;
            actionCameraOrbitOwned = true;
            setCameraAutoOrbitEnabled(true);
            actionCameraOrbitShouldPersist = step.orbit.persist !== false;
          }
          markActionStepFinished(stepIndex, timestamp, runId);
          return;
        }
        cameraActionTrack = {
          stepIndex,
          runId,
          startTime: timestamp,
          durationMs,
          easing: String(step.easing || "easeInOutCubic"),
          fromPosition: camera.position.clone(),
          fromTarget: controls.target.clone(),
          fromUp: camera.up.clone(),
          fromFov: Number(camera.fov),
          fromViewOffset: normalizeActionViewOffset(currentActionCameraViewOffset),
          toPosition: resolved.toPosition,
          toTarget: resolved.toTarget,
          toUp: resolved.toUp,
          toFov: Number(resolved.toFov),
          toViewOffset: normalizeActionViewOffset(resolved.toViewOffset),
          orbitEnabled: Boolean(step.orbit && step.orbit.enabled),
          orbitSpeedMultiplier: step.orbit ? Math.max(Number(step.orbit.speed_multiplier) || 1.0, 0.05) : 1.0,
          orbitDirection: step.orbit && Number(step.orbit.direction) < 0.0 ? -1.0 : 1.0,
          orbitPersist: step.orbit ? step.orbit.persist !== false : true,
          orbitOwned: Boolean(step.orbit && step.orbit.enabled),
        };
      }

      function startTimeAction(step, stepIndex, now) {
        const timestamp = actionTimestamp(now);
        const runId = activeActionRun ? activeActionRun.id : null;
        const direction = String(step.direction) === "backward" ? -1 : 1;
        const intervalMs = Math.max(Number(step.interval_ms) || playbackIntervalMs, 1);
        const fromFrameValue = clampFrameValue(displayedFrameValue);
        let targetFrameValue = Number(targetFrameIndexForTimeAction(step));
        let durationMs = Math.abs(targetFrameValue - fromFrameValue) * intervalMs;
        if (step.stop_after_ms !== undefined) {
          durationMs = Math.max(Number(step.stop_after_ms) || 0.0, 0.0);
          targetFrameValue = clampFrameValue(
            fromFrameValue + direction * (durationMs / intervalMs)
          );
        }
        const record = actionStepRecord(stepIndex, runId);
        if (record) record.effectiveDurationMs = durationMs;
        if (durationMs <= 0.0 || Math.abs(targetFrameValue - fromFrameValue) <= 1e-9) {
          displayedFrameValue = clampFrameValue(targetFrameValue);
          currentFrameIndex = clampFrameIndex(displayedFrameValue);
          renderActionInterpolatedFrameValue(displayedFrameValue, {
            updateWidgets: false,
            preserveCamera: true,
          });
          markActionStepFinished(stepIndex, timestamp, runId);
          return;
        }
        timeActionTrack = {
          stepIndex,
          runId,
          direction,
          intervalMs,
          startTime: timestamp,
          durationMs,
          fromFrameValue,
          targetFrameValue,
          easing: "easeInOutCubic",
        };
        playbackDirection = direction;
        lastPlaybackAdvanceTimestamp = null;
        updatePlaybackButtons();
      }

      function startStateAction(step, stepIndex, now) {
        if (!ovizStateControllerReady) {
          throw new Error("The States controller is not ready for this Action.");
        }
        let target;
        const runId = activeActionRun ? activeActionRun.id : null;
        try {
          target = step.runtime_snapshot
            ? {
                id: null,
                index: -1,
                name: String(step.runtime_name || "Original action view"),
                snapshot: safeJsonClone(step.runtime_snapshot, {}),
                transition: {
                  duration_ms: Math.max(Number(step.duration_ms) || 1200, 0.0),
                  easing: String(step.easing || "easeInOutCubic"),
                },
              }
            : ovizStateTargetFor(step.state);
        } catch (err) {
          throw new Error(`State Action target could not be resolved: ${String(err && err.message || err)}`);
        }
        const statePromise = ovizBeginStateTransition(target, { preserveActionRun: true });
        const record = actionStepRecord(stepIndex, runId);
        if (record && ovizStateTransition && ovizStateTransition.phasePlan) {
          record.effectiveDurationMs = Math.max(Number(ovizStateTransition.phasePlan.effectiveDurationMs) || 0.0, 0.0);
        }
        statePromise.then((result) => {
          if (!activeActionRun || String(activeActionRun.id) !== String(runId)) return;
          if (result && result.cancelled) {
            cancelActionRun(String(result.reason || "state-transition-cancelled"), {
              preserveLegendTransitionFrame: true,
              disableOrbit: true,
              cancelStateTransition: false,
            });
            return;
          }
          markActionStepFinished(stepIndex, null, runId);
        }).catch((err) => {
          if (activeActionRun && String(activeActionRun.id) === String(runId)) {
            failActionRun(err, actionStepRecord(stepIndex, runId));
          }
        });
      }

      function startActionStep(stepRecord, now) {
        const run = activeActionRun;
        if (!run || !stepRecord || stepRecord.started || !stepRecord.step) {
          return;
        }
        stepRecord.domains = actionDomainsForStepAtStart(stepRecord.step);
        stepRecord.domain = stepRecord.domains.join("+") || String(stepRecord.step.type || "action");
        if (!claimActionDomains(run, stepRecord)) {
          return;
        }
        const timestamp = actionTimestamp(now);
        stepRecord.started = true;
        stepRecord.startedAt = timestamp;
        stepRecord.runId = run.id;
        const step = stepRecord.step;
        if (step.type === "legend_group" || step.type === "camera") {
          stepRecord.effectiveDurationMs = Math.max(Number(step.duration_ms) || 0.0, 0.0);
        } else if (step.type === "time" && step.stop_after_ms !== undefined) {
          stepRecord.effectiveDurationMs = Math.max(Number(step.stop_after_ms) || 0.0, 0.0);
        }
        actionDebugState(run, stepRecord, {
          status: "running",
          phase: stepRecord.domain,
          progress: 0.0,
          effectiveDurationMs: stepRecord.effectiveDurationMs,
        });
        actionEvent("action-step-start", run, {
          step: stepRecord.index,
          stepType: String(step.type || ""),
          phase: stepRecord.domain,
          effectiveDurationMs: stepRecord.effectiveDurationMs,
        });
        if (step.type === "state") {
          startStateAction(step, stepRecord.index, timestamp);
          return;
        }
        if (step.type === "legend_group") {
          startLegendGroupAction(step, stepRecord.index, timestamp);
          return;
        }
        if (step.type === "camera") {
          startCameraAction(step, stepRecord.index, timestamp);
          return;
        }
        startTimeAction(step, stepRecord.index, timestamp);
      }

      function beginActionRun(action, now) {
        const timestamp = actionTimestamp(now);
        actionRunSerial += 1;
        const run = {
          id: `action-${actionRunSerial}`,
          key: String(action.key || ""),
          label: String(action.label || action.key || "Action"),
          steps: buildScheduledActionSteps(action),
          domainOwners: new Map(),
          startedAt: timestamp,
          status: "running",
        };
        activeActionRun = run;
        observedActionRun = run;
        activeActionKey = run.key;
        actionDebugState(run, null, { status: "running", phase: "scheduled", progress: 0.0 });
        actionEvent("action-start", run, { stepCount: run.steps.length });
        return run;
      }

      function startViewerAction(actionKey) {
        const action = actionDefinitionByKey(actionKey);
        if (!action) {
          return false;
        }
        try {
          const now = actionTimestamp();
          const stateBackedAction = Array.isArray(action.steps)
            && action.steps.length === 1
            && action.steps[0].type === "state";
          if (ovizStateTransition && !stateBackedAction) {
            actionHeldTraceOpacityByKey = captureActionTraceOpacityMap();
            actionHeldAppearanceRollback = {
              fromSnapshot: safeJsonClone(ovizStateTransition.fromSnapshot, {}),
              targetSnapshot: safeJsonClone(ovizStateTransition.targetSnapshot, {}),
              frozenAppearanceProgress: clampRange(
                Number(ovizStateTransition.currentAppearanceProgress) || 0.0,
                0.0,
                1.0,
              ),
              currentAppearanceProgress: clampRange(
                Number(ovizStateTransition.currentAppearanceProgress) || 0.0,
                0.0,
                1.0,
              ),
              startedAt: now,
              durationMs: 240.0,
              fromMilkyWayOpacity: clampRange(Number(milkyWayViewOpacityScale) || 0.0, 0.0, 1.0),
              fromSkyOpacity: clampRange(Number(skyDomeViewOpacityScale) || 0.0, 0.0, 1.0),
              toMilkyWayOpacity: String(cameraViewMode || "free") === "earth" ? 0.0 : 1.0,
              toSkyOpacity: String(cameraViewMode || "free") === "earth" ? 1.0 : 0.0,
            };
            ovizCancelStateTransitionWithoutSnap("action-start", {
              restorePresentation: false,
              preserveRenderedSelection: true,
            });
          }
          withActionGuard(() => {
            cancelActionRun("replaced", {
              now,
              commitLegendTransition: true,
              preserveLegendTransitionFrame: true,
              preserveHeldAppearance: true,
              disableOrbit: true,
              cancelStateTransition: false,
            });
            const run = beginActionRun(action, now);
            selectedActionKey = String(action.key);
            syncActionButtons();
            if (
              actionHeldTraceOpacityByKey
              && !action.steps.some((step) => step && step.type === "legend_group")
            ) {
              beginLegendTransition(
                currentGroup,
                legendState,
                -1,
                now,
                240,
                "easeInOutCubic",
              );
            }
            return run;
          });
          return true;
        } catch (err) {
          if (activeActionRun) {
            failActionRun(err, null);
          } else {
            console.error("Oviz Action could not start", err);
          }
          return false;
        }
      }

      function startRestoreInitialView() {
        if (!initialActionViewState || !ovizStateControllerReady) {
          return false;
        }
        try {
          const now = actionTimestamp();
          withActionGuard(() => {
            cancelActionRun("restore", {
              now,
              commitLegendTransition: true,
              preserveLegendTransitionFrame: true,
              disableOrbit: true,
              clearSelectedAction: true,
              cancelStateTransition: false,
            });
            selectedActionKey = "";
            beginActionRun({
              key: "__restore__",
              label: "Restore original view",
              steps: [{
                type: "state",
                start: "after_previous",
                delay_ms: 0,
                runtime_snapshot: safeJsonClone(initialActionViewState, {}),
                runtime_name: "Original action view",
                duration_ms: 1200,
                easing: "easeInOutCubic",
              }],
            }, now);
            syncActionButtons();
          });
          return true;
        } catch (err) {
          if (activeActionRun) {
            failActionRun(err, null);
          } else {
            console.error("Oviz Action restore could not start", err);
          }
          return false;
        }
      }

      function updateLegendTransition(now, options = {}) {
        if (!legendTransitionState) {
          return false;
        }
        if (
          legendTransitionState.runId
          && (!activeActionRun || String(activeActionRun.id) !== String(legendTransitionState.runId))
        ) {
          legendTransitionState = null;
          return false;
        }
        const durationMs = Math.max(Number(legendTransitionState.durationMs) || 0, 1);
        const rawProgress = (now - Number(legendTransitionState.startTime)) / durationMs;
        const easedProgress = actionEasingValue(legendTransitionState.easing, rawProgress);
        legendTransitionState.progress = easedProgress;
        const record = actionStepRecord(legendTransitionState.stepIndex, legendTransitionState.runId);
        if (record) emitActionProgress(record, now, rawProgress, easedProgress, "appearance");
        if (options.render !== false && rawProgress < 1.0) {
          renderActionInterpolatedFrameValue(displayedFrameValue, {
            updateWidgets: false,
            preserveCamera: true,
            transitionOwnerToken: `action:${String(legendTransitionState.runId || "")}`,
          });
        }
        if (rawProgress >= 1.0) {
          finishLegendTransition(true, now, options);
        }
        return true;
      }

      function updateCameraAction(now) {
        if (!cameraActionTrack) {
          return false;
        }
        if (!activeActionRun || String(activeActionRun.id) !== String(cameraActionTrack.runId)) {
          cameraActionTrack = null;
          return false;
        }
        const durationMs = Math.max(Number(cameraActionTrack.durationMs) || 0, 1);
        const rawProgress = (now - Number(cameraActionTrack.startTime)) / durationMs;
        const easedProgress = actionEasingValue(cameraActionTrack.easing, rawProgress);
        const record = actionStepRecord(cameraActionTrack.stepIndex, cameraActionTrack.runId);
        if (record) emitActionProgress(record, now, rawProgress, easedProgress, "camera");
        camera.position.lerpVectors(cameraActionTrack.fromPosition, cameraActionTrack.toPosition, easedProgress);
        controls.target.lerpVectors(cameraActionTrack.fromTarget, cameraActionTrack.toTarget, easedProgress);
        camera.up.lerpVectors(cameraActionTrack.fromUp, cameraActionTrack.toUp, easedProgress).normalize();
        camera.fov = cameraActionTrack.fromFov + (cameraActionTrack.toFov - cameraActionTrack.fromFov) * easedProgress;
        camera.updateProjectionMatrix();
        applyActionCameraViewOffset({
          x: cameraActionTrack.fromViewOffset.x + (cameraActionTrack.toViewOffset.x - cameraActionTrack.fromViewOffset.x) * easedProgress,
          y: cameraActionTrack.fromViewOffset.y + (cameraActionTrack.toViewOffset.y - cameraActionTrack.fromViewOffset.y) * easedProgress,
        });
        controls.update();
        if (rawProgress >= 1.0) {
          stopCameraActionTrack({ complete: true, now });
        }
        return true;
      }

      function timeActionReachedStop(track, now) {
        return !track || actionTimestamp(now) - Number(track.startTime || 0.0)
          >= Math.max(Number(track.durationMs) || 0.0, 0.0);
      }

      function updateTimeAction(now) {
        if (!timeActionTrack) {
          return false;
        }
        if (!activeActionRun || String(activeActionRun.id) !== String(timeActionTrack.runId)) {
          timeActionTrack = null;
          return false;
        }
        if (timeActionReachedStop(timeActionTrack, now)) {
          stopTimeActionTrack({ complete: true, now });
          return true;
        }
        const rawProgress = clampRange(
          (now - Number(timeActionTrack.startTime)) / Math.max(Number(timeActionTrack.durationMs) || 1.0, 1.0),
          0.0,
          1.0,
        );
        const progress = actionEasingValue(timeActionTrack.easing, rawProgress);
        const record = actionStepRecord(timeActionTrack.stepIndex, timeActionTrack.runId);
        if (record) emitActionProgress(record, now, rawProgress, progress, "time");
        displayedFrameValue = clampFrameValue(
          timeActionTrack.fromFrameValue
          + (timeActionTrack.targetFrameValue - timeActionTrack.fromFrameValue) * progress
        );
        currentFrameIndex = clampFrameIndex(displayedFrameValue);
        renderActionInterpolatedFrameValue(displayedFrameValue, {
          updateWidgets: false,
          preserveCamera: true,
        });
        if (rawProgress >= 1.0) {
          stopTimeActionTrack({ complete: true, now });
        }
        return true;
      }

      function updateStateBackedActionProgress(now) {
        const run = activeActionRun;
        if (!run || typeof ovizStateTransition === "undefined" || !ovizStateTransition) {
          return;
        }
        const record = run.steps.find((item) => item.started && !item.finished && item.step.type === "state");
        if (!record) {
          return;
        }
        const effectiveDurationMs = Math.max(
          Number(ovizStateTransition.phasePlan && ovizStateTransition.phasePlan.effectiveDurationMs)
            || Number(record.effectiveDurationMs)
            || 0.0,
          0.0,
        );
        record.effectiveDurationMs = effectiveDurationMs;
        const rawProgress = effectiveDurationMs > 0.0
          ? clampRange((actionTimestamp(now) - Number(ovizStateTransition.startedAt || now)) / effectiveDurationMs, 0.0, 1.0)
          : 1.0;
        const phase = String(ovizStateTransition.currentPhase || record.domain || "state");
        emitActionProgress(record, now, rawProgress, rawProgress, phase);
      }

      function updateViewerActions(now, options = {}) {
        let rerendered = false;
        const renderSerialAtStart = actionSceneRenderSerial;
        const timestamp = actionTimestamp(now);
        try {
          if (!activeActionRun && observedActionRun && observedActionRun.status === "running") {
            const interruptedRun = observedActionRun;
            interruptedRun.status = "cancelled";
            observedActionRun = null;
            actionDebugState(interruptedRun, null, {
              status: "cancelled",
              phase: "cancelled",
              progress: 0.0,
            });
            actionEvent("action-cancel", interruptedRun, {
              reason: "state-navigation",
              elapsedMs: Math.max(0.0, timestamp - Number(interruptedRun.startedAt || timestamp)),
            });
          }
          const run = activeActionRun;
          if (run && run.status === "running" && Array.isArray(run.steps)) {
            for (let stepIndex = 0; stepIndex < run.steps.length; stepIndex += 1) {
              if (activeActionRun !== run || run.status !== "running") break;
              const stepRecord = run.steps[stepIndex];
              if (stepRecord.started) continue;
              const readyAt = actionStepReadyAt(run, stepRecord);
              if (readyAt !== null && timestamp >= readyAt && actionDomainsAvailable(run, stepRecord)) {
                startActionStep(stepRecord, timestamp);
              }
            }
          }
          updateStateBackedActionProgress(timestamp);
          if (!options.schedulerOnly) {
            const heldAppearanceUpdated = updateActionHeldAppearanceRollback(timestamp);
            const timeWillRender = Boolean(timeActionTrack);
            updateLegendTransition(timestamp, { render: !timeWillRender });
            updateCameraAction(timestamp);
            updateTimeAction(timestamp);
            rerendered = actionSceneRenderSerial !== renderSerialAtStart;
            if (heldAppearanceUpdated && !rerendered) {
              renderActionInterpolatedFrameValue(displayedFrameValue, {
                updateWidgets: false,
                preserveCamera: true,
                transitionOwnerToken: activeActionRun
                  ? `action:${String(activeActionRun.id || "")}`
                  : "action-held-appearance",
              });
              rerendered = true;
            }
          }
          if (
            activeActionRun
            && activeActionRun.status === "running"
            && Array.isArray(activeActionRun.steps)
            && activeActionRun.steps.every((step) => step.finished)
            && !legendTransitionState
            && !cameraActionTrack
            && !timeActionTrack
            && !actionHeldAppearanceRollback
            && !(typeof ovizHeldSelectionTransition !== "undefined" && ovizHeldSelectionTransition)
          ) {
            if (actionCameraOrbitOwned && !actionCameraOrbitShouldPersist) {
              cameraAutoOrbitSpeedMultiplier = 1.0;
              setCameraAutoOrbitEnabled(false);
              actionCameraOrbitOwned = false;
              actionCameraOrbitShouldPersist = true;
            }
            completeActionRun(timestamp);
          }
        } catch (err) {
          const record = activeActionRun && Array.isArray(activeActionRun.steps)
            ? activeActionRun.steps.find((item) => item.started && !item.finished) || null
            : null;
          failActionRun(err, record, timestamp);
        }
        return rerendered;
      }
""".strip()
