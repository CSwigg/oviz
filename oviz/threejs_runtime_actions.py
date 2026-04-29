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
      let legendTransitionState = null;
      let cameraActionTrack = null;
      let timeActionTrack = null;
      let actionInterruptGuardDepth = 0;
      let actionCameraOrbitOwned = false;
      let actionCameraOrbitShouldPersist = true;
      let initialActionViewState = null;
      let currentActionCameraViewOffset = { x: 0.0, y: 0.0 };

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
        if (isGalacticReferenceTrace(trace) && !galacticReferenceVisible) {
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
          return Boolean(stateMap[trace.key]);
        }
        return mode === true;
      }

      function actionTraceVisibilityState(trace) {
        if (!legendTransitionState || !trace || !trace.key) {
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
        if (!fromVisible && !toVisible) {
          return { visible: false, opacity: 0.0 };
        }
        const progress = clampRange(Number(legendTransitionState.progress) || 0.0, 0.0, 1.0);
        if (fromVisible && toVisible) {
          return { visible: true, opacity: 1.0 };
        }
        if (fromVisible) {
          return { visible: true, opacity: 1.0 - progress };
        }
        return { visible: true, opacity: progress };
      }

      function traceVisibilityOpacityMultiplier(trace) {
        const state = actionTraceVisibilityState(trace);
        return state ? clampRange(state.opacity, 0.0, 1.0) : 1.0;
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
        return {
          group: String(currentGroup || defaultGroup || ""),
          legendState: cloneLegendStateMap(legendState),
          frameIndex: clampFrameIndex(currentFrameIndex),
          camera: {
            position: camera.position.clone(),
            target: controls.target.clone(),
            up: camera.up.clone(),
            fov: Number(camera.fov),
            viewOffset: normalizeActionViewOffset(currentActionCameraViewOffset),
          },
        };
      }

      function beginLegendTransition(targetGroup, targetLegendState, stepIndex, now, durationMs, easing) {
        const previousGroup = currentGroup;
        const previousLegendState = cloneLegendStateMap(legendState);
        const nextGroup = String(targetGroup || currentGroup);
        const nextLegendState = cloneLegendStateMap(targetLegendState);
        const nextStateMatchesCurrent = nextGroup === currentGroup
          && JSON.stringify(nextLegendState) === JSON.stringify(previousLegendState);
        currentGroup = nextGroup;
        if (groupSelectEl) {
          groupSelectEl.value = currentGroup;
        }
        syncLegendGroupChooser();
        renderLegend();
        const nextDurationMs = Math.max(Number(durationMs) || 0, 0);
        if (nextDurationMs <= 0 || nextStateMatchesCurrent) {
          legendState = nextLegendState;
          renderLegend();
          renderFrame(currentFrameIndex);
          markActionStepFinished(stepIndex);
          return;
        }
        legendTransitionState = {
          stepIndex,
          fromGroup: previousGroup,
          fromLegendState: previousLegendState,
          toGroup: nextGroup,
          toLegendState: nextLegendState,
          startTime: now,
          durationMs: nextDurationMs,
          easing: String(easing || "easeInOutCubic"),
          progress: 0.0,
        };
        renderFrame(currentFrameIndex);
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

      function estimateActionStepDurationMs(step) {
        if (!step || typeof step !== "object") {
          return 0;
        }
        if (step.type === "legend_group" || step.type === "camera") {
          return Math.max(Number(step.duration_ms) || 0, 0);
        }
        if (step.type !== "time") {
          return 0;
        }
        if (step.stop_after_ms !== undefined) {
          return Math.max(Number(step.stop_after_ms) || 0, 0);
        }
        const intervalMs = Math.max(Number(step.interval_ms) || playbackIntervalMs, 1);
        const targetIndex = targetFrameIndexForTimeAction(step);
        const frameDistance = Math.abs(targetIndex - clampFrameIndex(currentFrameIndex));
        return frameDistance * intervalMs;
      }

      function buildScheduledActionSteps(action) {
        const scheduled = [];
        let previousStartMs = 0.0;
        let previousEndMs = 0.0;
        (action.steps || []).forEach((step, index) => {
          const delayMs = Math.max(Number(step.delay_ms) || 0, 0);
          const startMs = (String(step.start || "after_previous") === "with_previous" ? previousStartMs : previousEndMs) + delayMs;
          const durationMs = estimateActionStepDurationMs(step);
          const endMs = startMs + durationMs;
          scheduled.push({
            index,
            step,
            startMs,
            endMs,
            started: false,
            finished: false,
          });
          previousStartMs = startMs;
          previousEndMs = endMs;
        });
        return scheduled;
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

      function actionStepRecord(stepIndex) {
        if (!activeActionRun || !Array.isArray(activeActionRun.steps)) {
          return null;
        }
        return activeActionRun.steps.find((step) => step.index === stepIndex) || null;
      }

      function markActionStepFinished(stepIndex) {
        const record = actionStepRecord(stepIndex);
        if (record) {
          record.finished = true;
        }
      }

      function finishLegendTransition(commitTarget = true) {
        if (!legendTransitionState) {
          return false;
        }
        const stepIndex = legendTransitionState.stepIndex;
        const nextGroup = legendTransitionState.toGroup;
        const nextLegendState = cloneLegendStateMap(legendTransitionState.toLegendState);
        legendTransitionState = null;
        if (commitTarget) {
          currentGroup = nextGroup;
          legendState = nextLegendState;
          if (groupSelectEl) {
            groupSelectEl.value = currentGroup;
          }
          syncLegendGroupChooser();
        }
        renderLegend();
        renderFrame(currentFrameIndex);
        markActionStepFinished(stepIndex);
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
        markActionStepFinished(track.stepIndex);
        return true;
      }

      function stopTimeActionTrack() {
        if (!timeActionTrack) {
          return false;
        }
        const track = timeActionTrack;
        timeActionTrack = null;
        playbackDirection = 0;
        lastPlaybackAdvanceTimestamp = null;
        updatePlaybackButtons();
        markActionStepFinished(track.stepIndex);
        return true;
      }

      function clearActiveActionState(options = {}) {
        if (legendTransitionState) {
          finishLegendTransition(options.commitLegendTransition !== false);
        }
        if (cameraActionTrack) {
          stopCameraActionTrack({
            complete: options.completeCamera !== false,
            activateOrbit: false,
          });
        }
        if (timeActionTrack) {
          stopTimeActionTrack();
        }
        if (options.disableOrbit && actionCameraOrbitOwned) {
          cameraAutoOrbitSpeedMultiplier = 1.0;
          setCameraAutoOrbitEnabled(false);
          actionCameraOrbitOwned = false;
          actionCameraOrbitShouldPersist = true;
        }
        activeActionRun = null;
        activeActionKey = "";
        if (options.clearSelectedAction) {
          selectedActionKey = "";
        }
        syncActionButtons();
      }

      function interruptActionRun(reason = "user", options = {}) {
        if (!activeActionRun || actionInterruptsMuted()) {
          return false;
        }
        clearActiveActionState({
          commitLegendTransition: true,
          completeCamera: true,
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
        const nextFov = Math.max(Number(overrides.fov) || camera.fov, 10.0);
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
        const resolved = resolveCameraActionState(step);
        if (!resolved) {
          markActionStepFinished(stepIndex);
          return;
        }
        setCameraAutoOrbitEnabled(false);
        actionCameraOrbitOwned = false;
        actionCameraOrbitShouldPersist = true;
        cameraAutoOrbitSpeedMultiplier = 1.0;
        const durationMs = Math.max(Number(step.duration_ms) || 0, 0);
        if (durationMs <= 0) {
          camera.position.copy(resolved.toPosition);
          controls.target.copy(resolved.toTarget);
          camera.up.copy(resolved.toUp);
          camera.fov = Number(resolved.toFov);
          applyActionCameraViewOffset(resolved.toViewOffset);
          controls.update();
          if (step.orbit && step.orbit.enabled) {
            cameraAutoOrbitSpeedMultiplier = Math.max(Number(step.orbit.speed_multiplier) || 1.0, 0.05);
            cameraAutoOrbitDirection = Number(step.orbit.direction) < 0.0 ? -1.0 : 1.0;
            actionCameraOrbitOwned = true;
            setCameraAutoOrbitEnabled(true);
            actionCameraOrbitShouldPersist = step.orbit.persist !== false;
          }
          markActionStepFinished(stepIndex);
          return;
        }
        cameraActionTrack = {
          stepIndex,
          startTime: now,
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
        const direction = String(step.direction) === "backward" ? -1 : 1;
        timeActionTrack = {
          stepIndex,
          direction,
          intervalMs: Math.max(Number(step.interval_ms) || playbackIntervalMs, 1),
          startTime: now,
          lastAdvanceTimestamp: null,
          advancedFrames: 0,
          stopAtTimeMyr: step.stop_at_time_myr === undefined ? null : Number(step.stop_at_time_myr),
          stopAfterFrames: step.stop_after_frames === undefined ? null : Math.max(Number(step.stop_after_frames) || 0, 0),
          stopAfterMs: step.stop_after_ms === undefined ? null : Math.max(Number(step.stop_after_ms) || 0, 0),
          targetIndex: targetFrameIndexForTimeAction(step),
        };
        playbackDirection = direction;
        lastPlaybackAdvanceTimestamp = null;
        updatePlaybackButtons();
      }

      function startActionStep(stepRecord, now) {
        if (!stepRecord || stepRecord.started || !stepRecord.step) {
          return;
        }
        stepRecord.started = true;
        const step = stepRecord.step;
        if (step.type === "legend_group") {
          startLegendGroupAction(step, stepRecord.index, now);
          return;
        }
        if (step.type === "camera") {
          startCameraAction(step, stepRecord.index, now);
          return;
        }
        startTimeAction(step, stepRecord.index, now);
      }

      function startViewerAction(actionKey) {
        const action = actionDefinitionByKey(actionKey);
        if (!action) {
          return false;
        }
        withActionGuard(() => {
          clearActiveActionState({
            commitLegendTransition: true,
            completeCamera: true,
            disableOrbit: true,
          });
          activeActionRun = {
            key: String(action.key),
            steps: buildScheduledActionSteps(action),
            startedAt: window.performance ? window.performance.now() : Date.now(),
          };
          activeActionKey = String(action.key);
          selectedActionKey = String(action.key);
          syncActionButtons();
        });
        return true;
      }

      function startRestoreInitialView() {
        if (!initialActionViewState) {
          return false;
        }
        const now = window.performance ? window.performance.now() : Date.now();
        const restoreAction = actionDefinitionByKey(selectedActionKey || activeActionKey);
        const restoreTimeIntervalMs = restoreTimeIntervalMsForAction(restoreAction);
        withActionGuard(() => {
          clearActiveActionState({
            commitLegendTransition: true,
            completeCamera: true,
            disableOrbit: true,
            clearSelectedAction: true,
          });
          const restoreSteps = [];
          let nextStepIndex = 0;
          if (
            String(currentGroup || "") !== String(initialActionViewState.group || "")
            || JSON.stringify(cloneLegendStateMap(legendState)) !== JSON.stringify(cloneLegendStateMap(initialActionViewState.legendState))
          ) {
            restoreSteps.push({ index: nextStepIndex, finished: false });
            beginLegendTransition(
              initialActionViewState.group,
              initialActionViewState.legendState,
              nextStepIndex,
              now,
              420,
              "easeInOutCubic",
            );
            nextStepIndex += 1;
          }
          if (clampFrameIndex(currentFrameIndex) !== clampFrameIndex(initialActionViewState.frameIndex)) {
            restoreSteps.push({ index: nextStepIndex, finished: false });
            timeActionTrack = {
              stepIndex: nextStepIndex,
              direction: clampFrameIndex(initialActionViewState.frameIndex) < clampFrameIndex(currentFrameIndex) ? -1 : 1,
              intervalMs: restoreTimeIntervalMs,
              startTime: now,
              lastAdvanceTimestamp: null,
              advancedFrames: 0,
              stopAtTimeMyr: null,
              stopAfterFrames: null,
              stopAfterMs: null,
              targetIndex: clampFrameIndex(initialActionViewState.frameIndex),
            };
            playbackDirection = timeActionTrack.direction;
            lastPlaybackAdvanceTimestamp = null;
            updatePlaybackButtons();
            nextStepIndex += 1;
          }
          restoreSteps.push({ index: nextStepIndex, finished: false });
          setCameraAutoOrbitEnabled(false);
          actionCameraOrbitOwned = false;
          actionCameraOrbitShouldPersist = true;
          cameraAutoOrbitSpeedMultiplier = 1.0;
          cameraAutoOrbitDirection = Number(
            initialState.global_controls && initialState.global_controls.camera_auto_orbit_direction
          ) < 0.0 ? -1.0 : 1.0;
          cameraActionTrack = {
            stepIndex: nextStepIndex,
            startTime: now,
            durationMs: 1050,
            easing: "easeInOutCubic",
            fromPosition: camera.position.clone(),
            fromTarget: controls.target.clone(),
            fromUp: camera.up.clone(),
            fromFov: Number(camera.fov),
            fromViewOffset: normalizeActionViewOffset(currentActionCameraViewOffset),
            toPosition: initialActionViewState.camera.position.clone(),
            toTarget: initialActionViewState.camera.target.clone(),
            toUp: initialActionViewState.camera.up.clone(),
            toFov: Number(initialActionViewState.camera.fov),
            toViewOffset: normalizeActionViewOffset(initialActionViewState.camera.viewOffset),
            orbitEnabled: Boolean(initialState.global_controls && initialState.global_controls.camera_auto_orbit_enabled),
            orbitSpeedMultiplier: 1.0,
            orbitDirection: Number(
              initialState.global_controls && initialState.global_controls.camera_auto_orbit_direction
            ) < 0.0 ? -1.0 : 1.0,
            orbitPersist: true,
            orbitOwned: false,
          };
          activeActionRun = {
            key: "__restore__",
            steps: restoreSteps,
            startedAt: now,
          };
          activeActionKey = "__restore__";
          syncActionButtons();
        });
        return true;
      }

      function updateLegendTransition(now) {
        if (!legendTransitionState) {
          return false;
        }
        const durationMs = Math.max(Number(legendTransitionState.durationMs) || 0, 1);
        const rawProgress = (now - Number(legendTransitionState.startTime)) / durationMs;
        const easedProgress = actionEasingValue(legendTransitionState.easing, rawProgress);
        legendTransitionState.progress = easedProgress;
        renderFrameScene(currentFrame(), frameTimeForValue(displayedFrameValue), { updateWidgets: false });
        if (rawProgress >= 1.0) {
          finishLegendTransition(true);
        }
        return true;
      }

      function updateCameraAction(now) {
        if (!cameraActionTrack) {
          return false;
        }
        const durationMs = Math.max(Number(cameraActionTrack.durationMs) || 0, 1);
        const rawProgress = (now - Number(cameraActionTrack.startTime)) / durationMs;
        const easedProgress = actionEasingValue(cameraActionTrack.easing, rawProgress);
        camera.position.lerpVectors(cameraActionTrack.fromPosition, cameraActionTrack.toPosition, easedProgress);
        controls.target.lerpVectors(cameraActionTrack.fromTarget, cameraActionTrack.toTarget, easedProgress);
        camera.up.lerpVectors(cameraActionTrack.fromUp, cameraActionTrack.toUp, easedProgress).normalize();
        camera.fov = cameraActionTrack.fromFov + (cameraActionTrack.toFov - cameraActionTrack.fromFov) * easedProgress;
        applyActionCameraViewOffset({
          x: cameraActionTrack.fromViewOffset.x + (cameraActionTrack.toViewOffset.x - cameraActionTrack.fromViewOffset.x) * easedProgress,
          y: cameraActionTrack.fromViewOffset.y + (cameraActionTrack.toViewOffset.y - cameraActionTrack.fromViewOffset.y) * easedProgress,
        });
        controls.update();
        if (rawProgress >= 1.0) {
          stopCameraActionTrack({ complete: true });
        }
        return true;
      }

      function timeActionReachedStop(track) {
        if (!track) {
          return true;
        }
        if (track.stopAfterFrames !== null && track.advancedFrames >= track.stopAfterFrames) {
          return true;
        }
        if (track.stopAfterMs !== null && (window.performance ? window.performance.now() : Date.now()) - track.startTime >= track.stopAfterMs) {
          return true;
        }
        if (track.stopAtTimeMyr !== null) {
          const currentTime = frameTimeForValue(currentFrameIndex);
          if (track.direction < 0) {
            return currentTime <= track.stopAtTimeMyr + 1e-9;
          }
          return currentTime >= track.stopAtTimeMyr - 1e-9;
        }
        return clampFrameIndex(currentFrameIndex) === clampFrameIndex(track.targetIndex);
      }

      function updateTimeAction(now) {
        if (!timeActionTrack) {
          return false;
        }
        if (timeActionReachedStop(timeActionTrack)) {
          stopTimeActionTrack();
          return true;
        }
        if (timeActionTrack.lastAdvanceTimestamp === null) {
          timeActionTrack.lastAdvanceTimestamp = now;
          return false;
        }
        const elapsedMs = Math.max(0.0, now - timeActionTrack.lastAdvanceTimestamp);
        if (elapsedMs < timeActionTrack.intervalMs) {
          return false;
        }
        const steps = Math.max(1, Math.floor(elapsedMs / timeActionTrack.intervalMs));
        timeActionTrack.lastAdvanceTimestamp += steps * timeActionTrack.intervalMs;
        let moved = false;
        for (let idx = 0; idx < steps; idx += 1) {
          const nextIndex = clampFrameIndex(currentFrameIndex + timeActionTrack.direction);
          if (nextIndex === currentFrameIndex) {
            stopTimeActionTrack();
            break;
          }
          renderFrame(nextIndex);
          timeActionTrack.advancedFrames += 1;
          moved = true;
          if (timeActionReachedStop(timeActionTrack)) {
            stopTimeActionTrack();
            break;
          }
        }
        return moved;
      }

      function updateViewerActions(now) {
        let rerendered = false;
        if (activeActionRun && Array.isArray(activeActionRun.steps)) {
          const elapsedMs = Math.max(0.0, now - Number(activeActionRun.startedAt || now));
          activeActionRun.steps.forEach((stepRecord) => {
            if (!stepRecord.started && elapsedMs >= stepRecord.startMs) {
              startActionStep(stepRecord, now);
            }
          });
        }
        if (updateLegendTransition(now)) {
          rerendered = true;
        }
        if (updateCameraAction(now)) {
          rerendered = true;
        }
        if (updateTimeAction(now)) {
          rerendered = true;
        }
        if (activeActionRun && Array.isArray(activeActionRun.steps) && activeActionRun.steps.every((step) => step.finished)) {
          if (actionCameraOrbitOwned && !actionCameraOrbitShouldPersist) {
            cameraAutoOrbitSpeedMultiplier = 1.0;
            setCameraAutoOrbitEnabled(false);
            actionCameraOrbitOwned = false;
            actionCameraOrbitShouldPersist = true;
          }
          activeActionRun = null;
          activeActionKey = "";
          syncActionButtons();
        }
        return rerendered;
      }
""".strip()
