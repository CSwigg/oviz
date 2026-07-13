"""Ordered viewer-state authoring, transitions, persistence, and automation API."""

THREEJS_STATE_RUNTIME_JS = r"""
      const OVIZ_STATES_VERSION = 1;
      const OVIZ_STATES_DEFAULT_TRANSITION = Object.freeze({
        duration_ms: 1200,
        easing: "easeInOutCubic",
      });
      const OVIZ_STATES_DB_NAME = "oviz-viewer-states";
      const OVIZ_STATES_DB_STORE = "drafts";
      let ovizStatesProject = null;
      let ovizOriginalRuntimeState = null;
      let ovizOriginalSceneInitialState = null;
      let ovizOriginalCameraResetState = null;
      let ovizStatesMode = "edit";
      let ovizActiveStateId = null;
      let ovizStateDirty = false;
      let ovizStateControllerReady = false;
      let ovizStateTransition = null;
      let ovizStateTransitionSerial = 0;
      let ovizStateTransitionTraceOpacity = null;
      let ovizStateSelectionTransition = null;
      let ovizHeldSelectionTransition = null;
      let ovizStatesShellEl = null;
      let ovizStatesDrawerEl = null;
      let ovizStatesRowsEl = null;
      let ovizStatesStatusEl = null;
      let ovizStatesRecoveryEl = null;
      let ovizStatesDraftTimer = null;
      let ovizStatesPreloadStatus = { loaded: 0, total: 0, failed: [] };
      let ovizResidentSkyBaseLayerKey = "";
      const ovizSkyLayerTransitionWaiters = new Map();

      function ovizTransitionDebugEnabled() {
        if (sceneSpec && sceneSpec.debug_transitions === true) return true;
        try {
          const query = new URLSearchParams(window.location.search || "");
          return ["1", "true", "yes", "on"].includes(
            String(query.get("ovizDebugTransitions") || "").toLowerCase()
          );
        } catch (_err) {
          return false;
        }
      }

      function ovizWriteTransitionDiagnostics(transition, values = {}) {
        if (!root || !root.dataset) return;
        root.dataset.transitionOwner = transition ? "state" : "";
        root.dataset.transitionRunId = transition ? String(transition.transitionId || "") : "";
        if (!ovizTransitionDebugEnabled()) return;
        const diagnostics = Object.assign({
          owner: transition ? "state" : "",
          runId: transition ? String(transition.transitionId || "") : "",
          targetId: transition ? transition.targetId : null,
          frameValue: Number(displayedFrameValue),
          canvasOpacity: Number(renderer && renderer.domElement && renderer.domElement.style.opacity || 1),
        }, values || {});
        root.dataset.transitionDiagnostics = JSON.stringify(diagnostics);
      }

      function ovizStatesClone(value, fallback = null) {
        try {
          return JSON.parse(JSON.stringify(value));
        } catch (_err) {
          return fallback;
        }
      }

      const OVIZ_MAX_VOLUME_MASK_SOURCE_COMPONENTS = 4;

      function ovizNormalizeVolumeMaskSourceComponents(components, maxComponents = OVIZ_MAX_VOLUME_MASK_SOURCE_COMPONENTS) {
        const merged = [];
        (Array.isArray(components) ? components : []).forEach((component) => {
          if (!component || !Number.isFinite(Number(component.weight))) return;
          const weight = Math.max(0.0, Number(component.weight));
          if (weight <= 1e-8) return;
          // A null mask is a real, unfiltered endpoint whose membership is one
          // everywhere.  Merge by runtime identity without cloning masks or
          // their GPU textures.
          const existing = merged.find((entry) => entry.mask === (component.mask || null));
          if (existing) {
            existing.weight += weight;
          } else {
            merged.push({ mask: component.mask || null, weight });
          }
        });
        const limit = Math.max(1, Math.floor(Number(maxComponents) || OVIZ_MAX_VOLUME_MASK_SOURCE_COMPONENTS));
        let retained = merged;
        if (retained.length > limit) {
          const strongest = new Set(
            retained
              .slice()
              .sort((left, right) => right.weight - left.weight)
              .slice(0, limit)
          );
          retained = retained.filter((component) => strongest.has(component));
        }
        const total = retained.reduce((sum, component) => sum + component.weight, 0.0);
        if (total <= 1e-8) return [{ mask: null, weight: 1.0 }];
        return retained.map((component) => ({
          mask: component.mask,
          weight: component.weight / total,
        }));
      }

      function ovizVolumeMaskTransitionSourceComponents(transition, fallbackMask = null) {
        if (transition && Array.isArray(transition.sourceMaskComponents)) {
          return ovizNormalizeVolumeMaskSourceComponents(transition.sourceMaskComponents);
        }
        if (transition) {
          const blend = clampRange(Number(transition.fromMaskBlend) || 0.0, 0.0, 1.0);
          return ovizNormalizeVolumeMaskSourceComponents([
            { mask: transition.fromMask || null, weight: 1.0 - blend },
            { mask: transition.fromSecondaryMask || null, weight: blend },
          ]);
        }
        return ovizNormalizeVolumeMaskSourceComponents([
          { mask: fallbackMask || null, weight: 1.0 },
        ]);
      }

      function ovizFreezeVolumeMaskSourceComponents(transition, fallbackMask = null) {
        const source = ovizVolumeMaskTransitionSourceComponents(transition, fallbackMask);
        if (!transition) return source;
        const progress = clampRange(Number(transition.progress) || 0.0, 0.0, 1.0);
        return ovizNormalizeVolumeMaskSourceComponents([
          ...source.map((component) => ({
            mask: component.mask,
            weight: component.weight * (1.0 - progress),
          })),
          // Null deliberately means the unfiltered, all-points endpoint.
          { mask: transition.toMask || null, weight: progress },
        ]);
      }

      function ovizStatesUuid(prefix = "state") {
        if (window.crypto && typeof window.crypto.randomUUID === "function") {
          return `${prefix}-${window.crypto.randomUUID()}`;
        }
        return `${prefix}-${Date.now().toString(36)}-${Math.random().toString(36).slice(2, 10)}`;
      }

      function ovizNormalizeTransition(value, fallback = OVIZ_STATES_DEFAULT_TRANSITION) {
        const source = value && typeof value === "object" ? value : {};
        const duration = Number(source.duration_ms ?? source.duration ?? fallback.duration_ms);
        const easing = String(source.easing || fallback.easing || "easeInOutCubic");
        return {
          duration_ms: Number.isFinite(duration) ? clampRange(duration, 0, 60000) : fallback.duration_ms,
          easing: ["linear", "easeInOutCubic", "easeInOutQuad", "easeOutCubic"].includes(easing)
            ? easing
            : fallback.easing,
        };
      }

      function ovizNormalizeStateRecord(value, index, seenIds) {
        const source = value && typeof value === "object" ? value : {};
        let id = String(source.id || "").trim() || ovizStatesUuid("state");
        while (seenIds.has(id)) {
          id = ovizStatesUuid("state");
        }
        seenIds.add(id);
        return {
          id,
          name: String(source.name || `State ${index + 1}`).trim() || `State ${index + 1}`,
          transition: source.transition ? ovizNormalizeTransition(source.transition) : null,
          snapshot: ovizStatesClone(source.snapshot || source.state || {}, {}),
          degraded: Boolean(source.degraded),
        };
      }

      function ovizNormalizeStatesProject(raw) {
        const embedded = raw && typeof raw === "object";
        const source = embedded ? raw : {};
        const rawItems = Array.isArray(source.items)
          ? source.items
          : (Array.isArray(source.ordered_states) ? source.ordered_states : []);
        const seenIds = new Set();
        return {
          schema_version: OVIZ_STATES_VERSION,
          project_id: String(source.project_id || "").trim() || ovizStatesUuid("project"),
          revision: Math.max(0, Math.floor(Number(source.revision) || 0)),
          default_mode: source.default_mode === "present" ? "present" : "edit",
          default_transition: ovizNormalizeTransition(source.default_transition),
          items: rawItems.map((item, index) => ovizNormalizeStateRecord(item, index, seenIds)),
          assets: source.assets && typeof source.assets === "object" ? ovizStatesClone(source.assets, {}) : {},
          embedded,
          synchronized_revision: Math.max(0, Math.floor(Number(source.synchronized_revision) || 0)),
        };
      }

      function ovizStatesPublicProject() {
        return {
          schema_version: OVIZ_STATES_VERSION,
          project_id: ovizStatesProject.project_id,
          revision: ovizStatesProject.revision,
          synchronized_revision: ovizStatesProject.synchronized_revision,
          default_mode: ovizStatesProject.default_mode,
          default_transition: ovizStatesClone(ovizStatesProject.default_transition, {}),
          items: ovizStatesClone(ovizStatesProject.items, []),
          assets: ovizStatesClone(ovizStatesProject.assets, {}),
        };
      }

      function ovizStateEvent(name, detail = {}) {
        const payload = Object.assign({
          rootId: root.id,
          projectId: ovizStatesProject ? ovizStatesProject.project_id : "",
        }, detail || {});
        root.dispatchEvent(new CustomEvent(name, { detail: payload }));
        try {
          if (window.parent && window.parent !== window) {
            window.parent.postMessage({ source: "oviz", type: name, detail: payload }, "*");
          }
        } catch (_err) {
          // Parent messaging is a convenience API; viewer events remain authoritative.
        }
      }

      function ovizStatesChanged(reason) {
        ovizStatesProject.revision += 1;
        ovizStateDirty = true;
        ovizScheduleDraftSave();
        ovizRenderStatesDrawer();
        if (typeof postSkyLayerStateToAladin === "function") {
          postSkyLayerStateToAladin();
        }
        ovizStateEvent("states-changed", {
          reason,
          revision: ovizStatesProject.revision,
          states: ovizStatesList(),
        });
      }

      function ovizResidentSkyLayers(currentLayers = []) {
        const current = Array.isArray(currentLayers) ? ovizStatesClone(currentLayers, []) : [];
        const layerKey = (layer) => String(layer && (layer.key || layer.survey) || "").trim();
        const currentByKey = new Map();
        current.forEach((layer) => {
          const key = layerKey(layer);
          if (key) currentByKey.set(key, layer);
        });
        if (!ovizResidentSkyBaseLayerKey && current.length) {
          ovizResidentSkyBaseLayerKey = layerKey(current[current.length - 1]);
        }
        const residentByKey = new Map();
        const residentOrder = [];
        const remember = (layers) => {
          (Array.isArray(layers) ? layers : []).forEach((layer) => {
            const key = layerKey(layer);
            if (!key) return;
            if (!residentByKey.has(key)) residentOrder.push(key);
            residentByKey.set(key, layer);
          });
        };
        remember(current);
        if (ovizStatesProject && Array.isArray(ovizStatesProject.items)) {
          ovizStatesProject.items.forEach((item) => remember(item && item.snapshot && item.snapshot.sky_layers));
        }
        const materialize = (key) => {
          const active = currentByKey.get(key);
          if (active) return ovizStatesClone(active, {});
          const resident = ovizStatesClone(residentByKey.get(key), {});
          resident.visible = false;
          resident.opacity = 0.0;
          return resident;
        };
        const layers = residentOrder
          .filter((key) => key !== ovizResidentSkyBaseLayerKey)
          .map(materialize);
        if (ovizResidentSkyBaseLayerKey && residentByKey.has(ovizResidentSkyBaseLayerKey)) {
          layers.push(materialize(ovizResidentSkyBaseLayerKey));
        }
        return layers;
      }

      function ovizStatesList() {
        return ovizStatesProject.items.map((item, index) => ({
          id: item.id,
          index: index + 1,
          name: item.name,
          transition: ovizStatesClone(item.transition, null),
          active: item.id === ovizActiveStateId,
          degraded: Boolean(item.degraded),
        }));
      }

      function ovizStateIndexFor(idOrIndex) {
        if (typeof idOrIndex === "number" && Number.isFinite(idOrIndex)) {
          const oneBased = Math.floor(idOrIndex);
          return oneBased === 0 ? -1 : oneBased - 1;
        }
        const text = String(idOrIndex ?? "");
        return ovizStatesProject.items.findIndex((item) => item.id === text);
      }

      function ovizStateTargetFor(idOrIndex) {
        if (idOrIndex === null || idOrIndex === undefined || idOrIndex === 0 || idOrIndex === "original") {
          return { id: null, index: -1, name: "Original", snapshot: ovizOriginalRuntimeState };
        }
        const index = ovizStateIndexFor(idOrIndex);
        if (index < 0 || index >= ovizStatesProject.items.length) {
          throw new Error(`Unknown Oviz state: ${String(idOrIndex)}`);
        }
        const item = ovizStatesProject.items[index];
        return Object.assign({ index }, item);
      }

      function ovizHydrateAssets(value) {
        if (Array.isArray(value)) {
          return value.map(ovizHydrateAssets);
        }
        if (!value || typeof value !== "object") {
          return value;
        }
        if (typeof value.__oviz_asset_ref__ === "string") {
          return ovizStatesProject.assets[value.__oviz_asset_ref__] || null;
        }
        const out = {};
        Object.entries(value).forEach(([key, nested]) => {
          out[key] = ovizHydrateAssets(nested);
        });
        return out;
      }

      async function ovizHashText(text) {
        if (window.crypto && window.crypto.subtle && typeof TextEncoder !== "undefined") {
          const bytes = new TextEncoder().encode(text);
          const digest = await window.crypto.subtle.digest("SHA-256", bytes);
          return Array.from(new Uint8Array(digest)).map((byte) => byte.toString(16).padStart(2, "0")).join("");
        }
        let hash = 2166136261;
        for (let index = 0; index < text.length; index += 1) {
          hash ^= text.charCodeAt(index);
          hash = Math.imul(hash, 16777619);
        }
        return `fnv1a-${(hash >>> 0).toString(16)}`;
      }

      async function ovizDeduplicateAssets(value, assets) {
        if (typeof value === "string" && value.startsWith("data:") && value.length >= 4096) {
          const hash = await ovizHashText(value);
          assets[hash] = value;
          return { __oviz_asset_ref__: hash };
        }
        if (Array.isArray(value)) {
          return Promise.all(value.map((nested) => ovizDeduplicateAssets(nested, assets)));
        }
        if (!value || typeof value !== "object") {
          return value;
        }
        const out = {};
        for (const [key, nested] of Object.entries(value)) {
          out[key] = await ovizDeduplicateAssets(nested, assets);
        }
        return out;
      }

      async function ovizCompactProjectForStorage(project) {
        const assets = {};
        const compact = ovizStatesClone(project, {});
        compact.items = [];
        for (const item of project.items) {
          compact.items.push(Object.assign({}, item, {
            snapshot: await ovizDeduplicateAssets(item.snapshot, assets),
          }));
        }
        compact.assets = assets;
        return compact;
      }

      function ovizRestoreOriginalCameraResetState() {
        if (!ovizOriginalCameraResetState) {
          return;
        }
        initialCameraState.position.copy(ovizOriginalCameraResetState.position);
        initialCameraState.target.copy(ovizOriginalCameraResetState.target);
        initialCameraState.up.copy(ovizOriginalCameraResetState.up);
        initialCameraState.fov = ovizOriginalCameraResetState.fov;
        initialCameraState.viewOffset = ovizStatesClone(ovizOriginalCameraResetState.viewOffset, { x: 0, y: 0 });
      }

      function ovizApplyCapturedCameraState(snapshot) {
        const savedCamera = snapshot && snapshot.camera || {};
        const savedGlobal = snapshot && snapshot.global_controls || {};
        const position = savedCamera.position || {};
        const target = savedCamera.target || {};
        const up = savedCamera.up || {};
        if ([position.x, position.y, position.z].every((value) => Number.isFinite(Number(value)))) {
          camera.position.set(Number(position.x), Number(position.y), Number(position.z));
        }
        if ([target.x, target.y, target.z].every((value) => Number.isFinite(Number(value)))) {
          controls.target.set(Number(target.x), Number(target.y), Number(target.z));
        }
        if ([up.x, up.y, up.z].every((value) => Number.isFinite(Number(value)))) {
          camera.up.set(Number(up.x), Number(up.y), Number(up.z)).normalize();
        }
        if (Number.isFinite(Number(savedGlobal.camera_fov))) {
          camera.fov = Number(savedGlobal.camera_fov);
        }
        const viewOffset = savedCamera.view_offset || savedCamera.viewOffset || { x: 0.0, y: 0.0 };
        ovizApplyInterpolatedViewOffset(viewOffset, viewOffset, 1.0);
        camera.updateProjectionMatrix();
        camera.lookAt(controls.target);
        camera.updateMatrixWorld(true);
        controls.update();
      }

      function ovizStateFidelityDifferences(expected, actual, path = "state", differences = []) {
        if (differences.length >= 64) {
          return differences;
        }
        if (typeof expected === "number" && typeof actual === "number") {
          const tolerance = 1e-6 * Math.max(1.0, Math.abs(expected));
          if (!Number.isFinite(actual) || Math.abs(expected - actual) > tolerance) {
            differences.push(path);
          }
          return differences;
        }
        if (Array.isArray(expected)) {
          if (!Array.isArray(actual) || expected.length !== actual.length) {
            differences.push(path);
            return differences;
          }
          expected.forEach((item, index) => {
            ovizStateFidelityDifferences(item, actual[index], `${path}[${index}]`, differences);
          });
          return differences;
        }
        if (expected && typeof expected === "object") {
          if (!actual || typeof actual !== "object" || Array.isArray(actual)) {
            differences.push(path);
            return differences;
          }
          Object.keys(expected).forEach((key) => {
            ovizStateFidelityDifferences(expected[key], actual[key], `${path}.${key}`, differences);
          });
          return differences;
        }
        if (expected !== actual) {
          differences.push(path);
        }
        return differences;
      }

      function ovizApplyStateImmediately(snapshot, options = {}) {
        const hydrated = ovizHydrateAssets(snapshot || {});
        const playback = hydrated.playback_state || {};
        const targetPlaybackDirection = Number(playback.direction) || 0;
        const targetOrbitEnabled = Boolean(
          hydrated.global_controls && hydrated.global_controls.camera_auto_orbit_enabled
        );
        applyViewerStateSyncInternal(hydrated, options);
        playbackDirection = 0;
        lastPlaybackAdvanceTimestamp = null;
        setCameraAutoOrbitEnabled(false);
        const frameValue = Number(hydrated.current_frame_value);
        const renderTargetFrame = () => {
          if (Number.isFinite(frameValue)) {
            renderInterpolatedFrameValue(frameValue);
          } else {
            renderFrame(Number(hydrated.current_frame_index) || 0);
          }
        };
        renderTargetFrame();
        ovizApplyCapturedCameraState(hydrated);
        if (Number.isFinite(Number(playback.interval_ms))) {
          playbackIntervalMs = Math.max(80, Number(playback.interval_ms));
        }
        ovizRestoreOriginalCameraResetState();
        const lassoRestoreOptions = Object.prototype.hasOwnProperty.call(options, "runtimeLassoMask")
          ? { runtimeMask: options.runtimeLassoMask }
          : {};
        return restoreLassoSelectionMask(hydrated.lasso_selection_mask, lassoRestoreOptions).then(() => {
          renderTargetFrame();
          ovizApplyCapturedCameraState(hydrated);
          updateSelectionUI();
          updatePlaybackButtons();
          renderLegend();
          resize();
          ovizApplyCapturedCameraState(hydrated);
          playbackDirection = targetPlaybackDirection;
          lastPlaybackAdvanceTimestamp = null;
          setCameraAutoOrbitEnabled(targetOrbitEnabled);
          updatePlaybackButtons();
          updateTimelineMotionOpacity();
          return captureRuntimeState();
        });
      }

      function ovizPointFrom(snapshot, key, fallback) {
        const cameraState = snapshot && snapshot.camera && snapshot.camera[key];
        return {
          x: Number.isFinite(Number(cameraState && cameraState.x)) ? Number(cameraState.x) : fallback.x,
          y: Number.isFinite(Number(cameraState && cameraState.y)) ? Number(cameraState.y) : fallback.y,
          z: Number.isFinite(Number(cameraState && cameraState.z)) ? Number(cameraState.z) : fallback.z,
        };
      }

      function ovizLerp(a, b, t) {
        return Number(a) + (Number(b) - Number(a)) * t;
      }

      function ovizEasing(name, progress) {
        if (typeof actionEasingValue === "function") {
          return actionEasingValue(name, progress);
        }
        const t = clampRange(progress, 0, 1);
        return t < 0.5 ? 4 * t * t * t : 1 - Math.pow(-2 * t + 2, 3) / 2;
      }

      function ovizTraceTargetVisible(trace, snapshot) {
        const groupName = String(snapshot.current_group || currentGroup);
        const legend = snapshot.legend_state || {};
        const defaults = groupVisibility[groupName] || {};
        if (
          typeof OvizTransitionRuntime !== "undefined"
          && OvizTransitionRuntime
          && typeof OvizTransitionRuntime.resolveTraceVisibility === "function"
        ) {
          return OvizTransitionRuntime.resolveTraceVisibility(trace, defaults, legend);
        }
        const mode = defaults[trace.key];
        if (mode === false || mode === undefined) return false;
        if (trace.showlegend) {
          return Object.prototype.hasOwnProperty.call(legend, trace.key)
            ? Boolean(legend[trace.key])
            : mode === true;
        }
        return mode === true;
      }

      function stateTraceVisibilityState(trace) {
        if (!ovizStateTransitionTraceOpacity || !trace) {
          return null;
        }
        const opacity = ovizStateTransitionTraceOpacity.get(trace.key);
        return opacity === undefined ? null : { visible: opacity > 0.0001, opacity };
      }

      function ovizCancelActionWithoutSnap(reason = "state-navigation") {
        if (typeof cancelActionRun === "function" && activeActionRun) {
          const cancelled = cancelActionRun(reason, {
            preserveLegendTransitionFrame: true,
            disableOrbit: true,
            cancelStateTransition: false,
            restorePresentation: false,
          });
          // ovizGoToStateRecord captured the exact live trace opacities before
          // cancelling this Action and copied them into the new State tracks.
          // Do not leave the Action's held opacity map installed: trace
          // rendering multiplies that map with the State tracks, so the
          // animated destination can look correct and then lose points when
          // the exact destination frame is rebuilt at completion.
          actionHeldTraceOpacityByKey = null;
          return cancelled;
        }
        activeActionRun = null;
        activeActionKey = "";
        cameraActionTrack = null;
        timeActionTrack = null;
        legendTransitionState = null;
        // State navigation owns the complete rendered result.  A legacy
        // Action may leave this transient opacity map behind when it is
        // interrupted; keeping it would continue to multiply the State
        // transition (and even its exact final restoration) afterward.
        actionHeldTraceOpacityByKey = null;
        actionCameraOrbitOwned = false;
        actionCameraOrbitShouldPersist = true;
        syncActionButtons();
        return false;
      }

      function ovizCancelStateTransitionWithoutSnap(reason = "cancelled", options = {}) {
        const transition = ovizStateTransition;
        if (!transition) {
          return false;
        }
        const selectionTransition = ovizStateSelectionTransition;
        if (options.preserveRenderedSelection === true && selectionTransition) {
          const traceCandidates = ovizTraceCandidatesForSnapshots(
            transition.fromSnapshot,
            transition.targetSnapshot,
          );
          const liveWeights = ovizCaptureSelectionMembershipOpacityMap(traceCandidates);
          const liveDefaultWeight = ovizSelectionMembershipOpacity("", null);
          const sourceMaskComponents = ovizFreezeVolumeMaskSourceComponents(
            selectionTransition,
            activeVolumeLassoSelectionMask(),
          );
          ovizHeldSelectionTransition = {
            fromWeightByKey: liveWeights,
            fromDefaultWeight: liveDefaultWeight,
            toEndpoint: {
              keys: new Set(Array.from(selectedClusterKeys || [])),
              filtered: lassoSelectionFilterActive(),
              hasMask: Boolean(currentLassoSelectionMask),
            },
            sourceMaskComponents,
            fromMask: sourceMaskComponents[0] ? sourceMaskComponents[0].mask : null,
            fromSecondaryMask: sourceMaskComponents[1] ? sourceMaskComponents[1].mask : null,
            fromMaskBlend: sourceMaskComponents[1] ? sourceMaskComponents[1].weight : 0.0,
            toMask: currentLassoSelectionMask,
            progress: 0.0,
            rollback: true,
          };
        }
        ovizStateTransition = null;
        ovizStateTransitionTraceOpacity = null;
        ovizStateSelectionTransition = null;
        ovizStateTimelineMotionActive = false;
        if (typeof cancelSkyViewTransitionAnimations === "function") {
          cancelSkyViewTransitionAnimations({
            reason,
            // A State retarget sends a replacement layer transition below.
            // Keeping the iframe's live resident layer stack until that message
            // arrives lets it retarget from the pixels actually on screen.
            cancelLayerTransition: options.preserveSkyLayerTransition !== true,
          });
        }
        const resolveLayerTransition = ovizSkyLayerTransitionWaiters.get(String(transition.transitionId));
        if (resolveLayerTransition) {
          ovizSkyLayerTransitionWaiters.delete(String(transition.transitionId));
          resolveLayerTransition({ cancelled: true, reason });
        }
        if (options.restorePresentation !== false) {
          renderer.domElement.style.opacity = "1";
          setMilkyWayModelOpacityScale(cameraViewMode === "earth" ? 0.0 : 1.0);
          setSkyDomeViewOpacityScale(cameraViewMode === "earth" ? 1.0 : 0.0, { force: false });
          updateTimelineMotionOpacity();
        }
        if (
          options.preserveSkyLayerTransition !== true
          && skyDomeFrameEl
          && skyDomeFrameEl.contentWindow
          && transition.skyLayerTransitionStarted
        ) {
          try {
            skyDomeFrameEl.contentWindow.postMessage({
              type: "oviz-sky-layer-transition-cancel",
              transitionId: transition.transitionId,
              reason,
            }, "*");
          } catch (_err) {
          }
        }
        if (typeof transition.resolve === "function") {
          transition.resolve({ cancelled: true, id: transition.targetId, reason });
        }
        ovizWriteTransitionDiagnostics(null, {
          cancelledRunId: transition.transitionId,
          reason,
        });
        ovizStateEvent("transition-cancel", {
          owner: "state",
          runId: transition.transitionId,
          id: transition.targetId,
          index: transition.targetIndex + 1,
          reason,
        });
        const preservedMasks = options.preserveSelectionMasks instanceof Set
          ? options.preserveSelectionMasks
          : new Set();
        if (options.preserveRenderedSelection === true && ovizHeldSelectionTransition) {
          (ovizHeldSelectionTransition.sourceMaskComponents || []).forEach((component) => {
            if (component && component.mask) preservedMasks.add(component.mask);
          });
          [
            ovizHeldSelectionTransition.fromMask,
            ovizHeldSelectionTransition.fromSecondaryMask,
            ovizHeldSelectionTransition.toMask,
          ].forEach((mask) => {
            if (mask) preservedMasks.add(mask);
          });
        }
        const cancelledMasks = new Set([transition.targetRuntimeLassoMask]);
        if (selectionTransition) {
          (selectionTransition.sourceMaskComponents || []).forEach((component) => {
            if (component && component.mask) cancelledMasks.add(component.mask);
          });
          [
            selectionTransition.fromMask,
            selectionTransition.fromSecondaryMask,
            selectionTransition.toMask,
          ].forEach((mask) => cancelledMasks.add(mask));
        }
        cancelledMasks.forEach((mask) => {
          if (mask && !preservedMasks.has(mask)) disposeSelectionMaskIfUnused(mask);
        });
        ovizRenderStatesDrawer();
        return true;
      }

      function ovizCurrentTransitionSnapshot() {
        const snapshot = captureRuntimeState();
        if (ovizStateTransition) {
          snapshot.current_frame_value = displayedFrameValue;
          snapshot.camera = {
            position: { x: camera.position.x, y: camera.position.y, z: camera.position.z },
            target: { x: controls.target.x, y: controls.target.y, z: controls.target.z },
            up: { x: camera.up.x, y: camera.up.y, z: camera.up.z },
            view_offset: ovizStatesClone(currentActionCameraViewOffset, { x: 0, y: 0 }),
          };
        }
        return snapshot;
      }

      function ovizTransitionVisualSignature(snapshot) {
        const state = snapshot || {};
        return JSON.stringify({
          frame: Number.isFinite(Number(state.current_frame_value))
            ? Number(state.current_frame_value)
            : Number(state.current_frame_index) || 0,
          group: state.current_group || "",
          legend: state.legend_state || {},
          traces: state.trace_style_state || {},
          volumes: state.volume_state_by_key || {},
          selections: {
            current: state.current_selection || null,
            currentMany: state.current_selections || [],
            selectedKeys: state.selected_cluster_keys || [],
            lassoMask: state.lasso_selection_mask || null,
          },
          sky: {
            layers: state.sky_layers || [],
            activeLayer: state.active_sky_layer_key || "",
            source: (state.global_controls || {}).sky_dome_source_key || "",
            hidden: Boolean((state.global_controls || {}).sky_background_hidden),
          },
          controls: state.global_controls || {},
          labels: state.manual_labels || [],
          panels: {
            widgets: state.widgets || {},
            legend: state.legend_panel_state || null,
            scaleBar: state.scale_bar_state || null,
          },
          selectionBox: state.selection_box_state || null,
          filters: state.cluster_filter_state || null,
          zen: Boolean(state.zen_mode_enabled),
        });
      }

      function ovizStateFrameValue(snapshot) {
        return Number.isFinite(Number(snapshot && snapshot.current_frame_value))
          ? Number(snapshot.current_frame_value)
          : Number(snapshot && snapshot.current_frame_index) || 0.0;
      }

      function ovizStateCameraSignature(snapshot) {
        const state = snapshot || {};
        const controlsState = state.global_controls || {};
        return JSON.stringify({
          camera: state.camera || {},
          viewMode: controlsState.camera_view_mode || "free",
          fov: controlsState.camera_fov,
        });
      }

      function ovizStateAppearanceSignature(snapshot) {
        const state = snapshot || {};
        const controlsState = ovizStatesClone(state.global_controls || {}, {});
        [
          "camera_view_mode",
          "camera_fov",
          "camera_auto_orbit_enabled",
          "camera_auto_orbit_direction",
          "camera_auto_orbit_speed_multiplier",
        ].forEach((key) => delete controlsState[key]);
        return JSON.stringify({
          group: state.current_group || "",
          legend: state.legend_state || {},
          traceStyles: state.trace_style_state || {},
          volumes: state.volume_state_by_key || {},
          selections: {
            current: state.current_selection || null,
            currentMany: state.current_selections || [],
            mode: state.current_selection_mode || "none",
            selectedKeys: state.selected_cluster_keys || [],
            lassoMask: state.lasso_selection_mask || null,
            lassoFilter: state.lasso_selection_filter_enabled !== false,
          },
          skyLayers: state.sky_layers || [],
          controls: controlsState,
          labels: state.manual_labels || [],
          widgets: state.widgets || {},
          legendPanel: state.legend_panel_state || null,
          scaleBar: state.scale_bar_state || null,
          selectionBox: state.selection_box_state || null,
          clusterFilter: state.cluster_filter_state || null,
          zen: Boolean(state.zen_mode_enabled),
          renderedTraceOpacity: state.__oviz_rendered_trace_opacity || null,
          renderedSelectionWeights: state.__oviz_rendered_selection_weights || null,
          renderedVolumeMaskBlend: state.__oviz_rendered_volume_mask_blend || null,
        });
      }

      function ovizBuildTransitionPhases(from, to, transitionSpec) {
        const changed = {
          camera: ovizStateCameraSignature(from) !== ovizStateCameraSignature(to),
          appearance: ovizStateAppearanceSignature(from) !== ovizStateAppearanceSignature(to),
          time: Math.abs(ovizStateFrameValue(from) - ovizStateFrameValue(to)) > 1e-9,
        };
        const names = ["camera", "appearance", "time"].filter((name) => changed[name]);
        const requestedDurationMs = Math.max(Number(transitionSpec && transitionSpec.duration_ms) || 0.0, 0.0);
        if (!names.length) {
          return {
            changed,
            phases: [{ name: "settle", startMs: 0.0, endMs: 0.0 }],
            effectiveDurationMs: 0.0,
            requestedDurationMs,
          };
        }
        const phaseMinimumDurationMs = 800.0;
        const minimumDurationMs = names.length * phaseMinimumDurationMs;
        const effectiveDurationMs = Math.max(requestedDurationMs, minimumDurationMs);
        const baseWeight = { camera: 0.35, appearance: 0.30, time: 0.35, settle: 1.0 };
        const totalWeight = names.reduce((sum, name) => sum + baseWeight[name], 0.0) || 1.0;
        const extraDurationMs = Math.max(effectiveDurationMs - minimumDurationMs, 0.0);
        let cursorMs = 0.0;
        const phases = names.map((name, index) => {
          const durationMs = phaseMinimumDurationMs + extraDurationMs * (baseWeight[name] / totalWeight);
          const phase = {
            name,
            startMs: cursorMs,
            endMs: index === names.length - 1 ? effectiveDurationMs : cursorMs + durationMs,
          };
          cursorMs = phase.endMs;
          return phase;
        });
        return { changed, phases, effectiveDurationMs, requestedDurationMs };
      }

      function ovizTransitionPhaseState(transition, elapsedMs) {
        const phases = transition.phasePlan.phases;
        const elapsed = clampRange(Number(elapsedMs) || 0.0, 0.0, transition.phasePlan.effectiveDurationMs);
        const phase = phases.find((item) => elapsed <= item.endMs + 1e-9) || phases[phases.length - 1];
        const durationMs = Math.max(phase.endMs - phase.startMs, 1e-9);
        return {
          name: phase.name,
          progress: clampRange((elapsed - phase.startMs) / durationMs, 0.0, 1.0),
          elapsedMs: elapsed,
        };
      }

      function ovizTransitionPhaseProgress(transition, phaseName, elapsedMs) {
        const phase = transition.phasePlan.phases.find((item) => item.name === phaseName);
        if (!phase) return 1.0;
        if (elapsedMs <= phase.startMs) return 0.0;
        if (elapsedMs >= phase.endMs) return 1.0;
        return clampRange((elapsedMs - phase.startMs) / Math.max(phase.endMs - phase.startMs, 1e-9), 0.0, 1.0);
      }

      function ovizTraceCandidatesForFrameValue(frameValue) {
        const candidates = new Map();
        const clampedValue = clampFrameValue(frameValue);
        const indices = new Set([
          clampFrameIndex(Math.floor(clampedValue)),
          clampFrameIndex(Math.ceil(clampedValue)),
        ]);
        indices.forEach((index) => {
          const frame = frameSpecs[index];
          (frame && Array.isArray(frame.traces) ? frame.traces : []).forEach((trace) => {
            if (trace && trace.key && !candidates.has(trace.key)) candidates.set(trace.key, trace);
          });
        });
        return candidates;
      }

      function ovizTraceCandidatesForSnapshots(from, to) {
        const candidates = new Map();
        const remember = (trace) => {
          if (trace && trace.key && !candidates.has(trace.key)) candidates.set(trace.key, trace);
        };
        (currentFrame() && currentFrame().traces || []).forEach(remember);
        const fromFrameValue = clampFrameValue(ovizStateFrameValue(from));
        const toFrameValue = clampFrameValue(ovizStateFrameValue(to));
        const firstIndex = clampFrameIndex(Math.floor(Math.min(fromFrameValue, toFrameValue)));
        const lastIndex = clampFrameIndex(Math.ceil(Math.max(fromFrameValue, toFrameValue)));
        for (let index = firstIndex; index <= lastIndex; index += 1) {
          const frame = frameSpecs[index];
          (frame && Array.isArray(frame.traces) ? frame.traces : []).forEach(remember);
        }
        return candidates;
      }

      function ovizSelectionEndpoint(snapshot) {
        const state = snapshot || {};
        const keys = new Set(
          (Array.isArray(state.selected_cluster_keys) ? state.selected_cluster_keys : [])
            .map((key) => normalizeMemberKey(key))
            .filter(Boolean)
        );
        const hasMask = Boolean(state.lasso_selection_mask);
        const filtered = state.lasso_selection_filter_enabled !== false
          && String(state.current_selection_mode || "none") === "lasso"
          && Boolean(keys.size || hasMask);
        return { keys, filtered, hasMask };
      }

      function ovizSelectionMembershipOpacity(pointKey, point = null) {
        const transition = ovizStateSelectionTransition || ovizHeldSelectionTransition;
        const key = normalizeMemberKey(pointKey);
        if (!transition) {
          if (!lassoSelectionFilterActive()) return 1.0;
          if (key && selectedClusterKeys.has(key)) return 1.0;
          if (
            !selectedClusterKeys.size
            && currentLassoSelectionMask
            && point
            && typeof pointInsideProjectedLassoMask === "function"
          ) {
            return pointInsideProjectedLassoMask(
              Number(point.x), Number(point.y), Number(point.z), currentLassoSelectionMask
            ) ? 1.0 : 0.0;
          }
          return 0.0;
        }
        const fromWeight = key && transition.fromWeightByKey.has(key)
          ? transition.fromWeightByKey.get(key)
          : transition.fromDefaultWeight;
        let toWeight = transition.toEndpoint.filtered ? 0.0 : 1.0;
        if (transition.toEndpoint.filtered && key && transition.toEndpoint.keys.has(key)) {
          toWeight = 1.0;
        } else if (
          transition.toEndpoint.filtered
          && !transition.toEndpoint.keys.size
          && transition.toMask
          && point
          && typeof pointInsideProjectedLassoMask === "function"
        ) {
          toWeight = pointInsideProjectedLassoMask(
            Number(point.x), Number(point.y), Number(point.z), transition.toMask
          ) ? 1.0 : 0.0;
        }
        return ovizLerp(fromWeight, toWeight, clampRange(Number(transition.progress) || 0.0, 0.0, 1.0));
      }

      function ovizCaptureSelectionMembershipOpacityMap(traceCandidates) {
        const result = new Map();
        traceCandidates.forEach((trace) => {
          (trace && Array.isArray(trace.points) ? trace.points : []).forEach((point) => {
            const key = normalizeMemberKey(
              clusterFilterSelectionKeyForPoint(point) || normalizedSelectionKeyFor(point.selection)
            );
            if (key && !result.has(key)) {
              result.set(key, ovizSelectionMembershipOpacity(key, point));
            }
          });
        });
        return result;
      }

      function ovizStateDestinationCameraState(snapshot) {
        const savedCamera = snapshot && snapshot.camera || {};
        const savedGlobal = snapshot && snapshot.global_controls || {};
        return {
          position: savedCamera.position,
          target: savedCamera.target,
          up: savedCamera.up,
          fov: savedGlobal.camera_fov,
          viewOffset: savedCamera.view_offset || savedCamera.viewOffset || { x: 0.0, y: 0.0 },
        };
      }

      function ovizTransitionEpochMs(performanceMs) {
        const origin = (
          typeof performance !== "undefined"
          && Number.isFinite(Number(performance.timeOrigin))
        )
          ? Number(performance.timeOrigin)
          : (Date.now() - Number(performanceMs || 0.0));
        return origin + Number(performanceMs || 0.0);
      }

      function ovizInterpolatedViewOffset(startValue, endValue, progress) {
        const start = typeof normalizeActionViewOffset === "function"
          ? normalizeActionViewOffset(startValue)
          : (startValue || { x: 0.0, y: 0.0 });
        const end = typeof normalizeActionViewOffset === "function"
          ? normalizeActionViewOffset(endValue)
          : (endValue || { x: 0.0, y: 0.0 });
        return {
          x: ovizLerp(Number(start.x) || 0.0, Number(end.x) || 0.0, progress),
          y: ovizLerp(Number(start.y) || 0.0, Number(end.y) || 0.0, progress),
        };
      }

      function ovizApplyInterpolatedViewOffset(startValue, endValue, progress) {
        if (typeof applyActionCameraViewOffset !== "function") {
          return;
        }
        applyActionCameraViewOffset(ovizInterpolatedViewOffset(startValue, endValue, progress));
      }

      function ovizPrepareDestinationSkyPresentation(snapshot, options = {}) {
        const global = snapshot && snapshot.global_controls || {};
        if (options.prepareOnly === true) {
          if (typeof updateSkyDomeCaptureFrame === "function") {
            updateSkyDomeCaptureFrame();
          }
          return;
        }
        restoreSkyLayerStateFromSnapshot(snapshot, {
          postToAladin: options.postLayersToAladin !== false,
        });
        if (typeof global.sky_dome_enabled === "boolean") {
          skyDomeSpec.enabled = global.sky_dome_enabled;
        }
        if (typeof global.sky_dome_force_visible === "boolean") {
          skyDomeForceVisible = global.sky_dome_force_visible;
        }
        if (Number.isFinite(Number(global.sky_dome_opacity))) {
          skyDomeSpec.opacity = clampRange(Number(global.sky_dome_opacity), 0.0, 1.0);
        }
        if (typeof global.sky_background_hidden === "boolean") {
          skyBackgroundHidden = global.sky_background_hidden;
        }
        if (
          typeof global.sky_dome_source_key === "string"
          && skyDomeSourceOptionByKey.has(global.sky_dome_source_key)
          && typeof setSkyDomeSourceByKey === "function"
        ) {
          setSkyDomeSourceByKey(global.sky_dome_source_key);
        }
        if (root && root.dataset) {
          root.dataset.skyBackgroundHidden = skyBackgroundHidden ? "true" : "false";
        }
        if (typeof updateSkyDomeCaptureFrame === "function") {
          updateSkyDomeCaptureFrame();
        }
      }

      function ovizStartSkyLayerTransition(transition, sourceLayers) {
        if (
          !transition
          || !skyDomeFrameEl
          || !skyDomeFrameEl.contentWindow
        ) {
          return false;
        }
        const fromLayers = ovizResidentSkyLayers(sourceLayers || []);
        const toLayers = ovizResidentSkyLayers(
          Array.isArray(transition.targetSnapshot && transition.targetSnapshot.sky_layers)
            ? transition.targetSnapshot.sky_layers
            : [],
        );
        if (
          !transition.forceSkyLayerRetarget
          && JSON.stringify(fromLayers) === JSON.stringify(toLayers)
        ) {
          transition.skyLayerPromise = Promise.resolve({ unchanged: true });
          return false;
        }
        try {
          const appearancePhase = transition.phasePlan.phases.find((phase) => phase.name === "appearance");
          const phaseDurationMs = appearancePhase
            ? Math.max(appearancePhase.endMs - appearancePhase.startMs, 0.0)
            : 0.0;
          const phaseStartedAtEpochMs = transition.startedAtEpochMs
            + (appearancePhase ? appearancePhase.startMs : 0.0);
          transition.skyLayerPromise = new Promise((resolve) => {
            ovizSkyLayerTransitionWaiters.set(String(transition.transitionId), resolve);
          });
          skyDomeFrameEl.contentWindow.postMessage({
            type: "oviz-sky-layer-transition",
            transitionId: transition.transitionId,
            durationMs: phaseDurationMs,
            easing: transition.transitionSpec.easing,
            startedAtEpochMs: phaseStartedAtEpochMs,
            residentStack: true,
            fromLayers,
            toLayers,
          }, "*");
          return true;
        } catch (_err) {
          ovizSkyLayerTransitionWaiters.delete(String(transition.transitionId));
          transition.skyLayerPromise = Promise.resolve({ failed: true });
          postSkyLayerStateToAladin();
          return false;
        }
      }

      function ovizCreateEarthCameraTrack(destination) {
        const destinationCamera = ovizStateDestinationCameraState(destination);
        const endPositionValue = destinationCamera.position || {};
        const endTargetValue = destinationCamera.target || {};
        const endPosition = new THREE.Vector3(
          Number(endPositionValue.x),
          Number(endPositionValue.y),
          Number(endPositionValue.z)
        );
        const endTarget = new THREE.Vector3(
          Number(endTargetValue.x),
          Number(endTargetValue.y),
          Number(endTargetValue.z)
        );
        if (
          !Number.isFinite(endPosition.x)
          || !Number.isFinite(endPosition.y)
          || !Number.isFinite(endPosition.z)
          || !Number.isFinite(endTarget.x)
          || !Number.isFinite(endTarget.y)
          || !Number.isFinite(endTarget.z)
        ) {
          return null;
        }
        const earthPoint = earthViewPoint();
        const startDirection = new THREE.Vector3().subVectors(controls.target, camera.position);
        if (startDirection.lengthSq() <= 1e-12) {
          camera.getWorldDirection(startDirection);
        }
        const endDirection = new THREE.Vector3().subVectors(endTarget, endPosition);
        if (startDirection.lengthSq() <= 1e-12 || endDirection.lengthSq() <= 1e-12) {
          return null;
        }
        startDirection.normalize();
        endDirection.normalize();
        const startRadius = Math.max(camera.position.distanceTo(earthPoint), 1e-9);
        const endRadius = Math.max(endPosition.distanceTo(endTarget), 1e-9);
        const endFov = Number(destinationCamera.fov);
        return {
          earthPoint,
          startDirection,
          endDirection,
          startRadius,
          endRadius,
          startFov: clampRange(Number(camera.fov) || 60.0, 0.05, 120.0),
          endFov: Number.isFinite(endFov) ? clampRange(endFov, 0.05, 120.0) : clampRange(Number(camera.fov) || 60.0, 0.05, 120.0),
          startViewOffset: ovizStatesClone(currentActionCameraViewOffset, { x: 0.0, y: 0.0 }),
          endViewOffset: ovizStatesClone(destinationCamera.viewOffset, { x: 0.0, y: 0.0 }),
        };
      }

      function ovizApplyEarthCameraTrack(track, progress) {
        if (!track) {
          return;
        }
        const direction = slerpSkyViewDirection(track.startDirection, track.endDirection, progress);
        const logRadius = ovizLerp(Math.log(track.startRadius), Math.log(track.endRadius), progress);
        const radius = Math.exp(logRadius);
        camera.position.copy(track.earthPoint).sub(direction.clone().multiplyScalar(radius));
        controls.target.copy(track.earthPoint);
        camera.up.copy(skyViewUpVectorForDirection(direction));
        camera.fov = interpolateCameraFovTangent(track.startFov, track.endFov, progress);
        ovizApplyInterpolatedViewOffset(track.startViewOffset, track.endViewOffset, progress);
        camera.updateProjectionMatrix();
        camera.lookAt(track.earthPoint);
        camera.updateMatrixWorld(true);
      }

      function ovizCreateNativeViewCameraTrack(destination) {
        const destinationCamera = ovizStateDestinationCameraState(destination);
        const endPositionValue = destinationCamera.position || {};
        const endTargetValue = destinationCamera.target || {};
        const endUpValue = destinationCamera.up || {};
        const endPosition = new THREE.Vector3(
          Number(endPositionValue.x),
          Number(endPositionValue.y),
          Number(endPositionValue.z)
        );
        const endTarget = new THREE.Vector3(
          Number(endTargetValue.x),
          Number(endTargetValue.y),
          Number(endTargetValue.z)
        );
        const endUp = new THREE.Vector3(
          Number(endUpValue.x),
          Number(endUpValue.y),
          Number(endUpValue.z)
        );
        if (
          !Number.isFinite(endPosition.x)
          || !Number.isFinite(endPosition.y)
          || !Number.isFinite(endPosition.z)
          || !Number.isFinite(endTarget.x)
          || !Number.isFinite(endTarget.y)
          || !Number.isFinite(endTarget.z)
        ) {
          return null;
        }
        const startPosition = camera.position.clone();
        const startTarget = controls.target.clone();
        const startUp = camera.up.clone().normalize();
        const startDirection = new THREE.Vector3().subVectors(startTarget, startPosition);
        const endDirection = new THREE.Vector3().subVectors(endTarget, endPosition);
        if (startDirection.lengthSq() <= 1e-12) camera.getWorldDirection(startDirection);
        if (endDirection.lengthSq() <= 1e-12) endDirection.copy(startDirection);
        startDirection.normalize();
        endDirection.normalize();
        if (endUp.lengthSq() <= 1e-12) endUp.copy(startUp);
        endUp.normalize();
        const endFovValue = Number(destinationCamera.fov);
        return {
          startPosition,
          endPosition,
          startTarget,
          endTarget,
          startUp,
          endUp,
          startDirection,
          endDirection,
          startDistance: Math.max(startPosition.distanceTo(startTarget), 1e-9),
          endDistance: Math.max(endPosition.distanceTo(endTarget), 1e-9),
          startFov: clampRange(Number(camera.fov) || 60.0, 0.05, 120.0),
          endFov: Number.isFinite(endFovValue)
            ? clampRange(endFovValue, 0.05, 120.0)
            : clampRange(Number(camera.fov) || 60.0, 0.05, 120.0),
          startViewOffset: ovizStatesClone(currentActionCameraViewOffset, { x: 0.0, y: 0.0 }),
          endViewOffset: ovizStatesClone(destinationCamera.viewOffset, { x: 0.0, y: 0.0 }),
        };
      }

      function ovizApplyNativeViewCameraTrack(track, progress) {
        if (!track) {
          return;
        }
        const t = clampRange(progress, 0.0, 1.0);
        if (t >= 1.0 - 1e-9) {
          camera.position.copy(track.endPosition);
          controls.target.copy(track.endTarget);
          camera.up.copy(track.endUp);
          camera.fov = track.endFov;
          ovizApplyInterpolatedViewOffset(track.endViewOffset, track.endViewOffset, 1.0);
        } else {
          const direction = slerpSkyViewDirection(track.startDirection, track.endDirection, t);
          const distance = ovizLerp(track.startDistance, track.endDistance, t);
          camera.position.copy(track.startPosition).lerp(track.endPosition, t);
          controls.target.copy(camera.position).add(direction.multiplyScalar(distance));
          camera.up.copy(track.startUp).lerp(track.endUp, t);
          if (camera.up.lengthSq() <= 1e-12) camera.up.copy(track.endUp);
          camera.up.normalize();
          camera.fov = interpolateCameraFovTangent(track.startFov, track.endFov, t);
          ovizApplyInterpolatedViewOffset(track.startViewOffset, track.endViewOffset, t);
        }
        camera.updateProjectionMatrix();
        camera.lookAt(controls.target);
        camera.updateMatrixWorld(true);
      }

      function ovizPrepareEarthViewModeForStateTransition() {
        if (cameraViewMode !== "earth") {
          earthViewReturnCameraState = captureEarthViewReturnCameraState();
        }
        const earthPoint = earthViewPoint();
        const targetPoint = earthViewTargetPoint();
        earthViewFocusDistance = targetPoint
          ? Math.max(targetPoint.distanceTo(earthPoint), 1e-6)
          : Math.max(Number(earthViewFocusDistance) || 8122.0, 1e-6);
        cameraViewMode = "earth";
        setSkyDomeViewOpacityScale(0.0, { force: false });
        applyGlobalControlState();
        applyCameraViewMode();
        buildAxes();
        renderFrame(currentFrameIndex);
      }

      function ovizBeginStateTransition(target, options = {}) {
        const destination = ovizHydrateAssets(target.snapshot || {});
        const transitionSpec = ovizNormalizeTransition(target.transition, ovizStatesProject.default_transition);
        const from = ovizCurrentTransitionSnapshot();
        const sourceAppearanceComposite = (
          typeof actionHeldAppearanceRollback !== "undefined"
          && actionHeldAppearanceRollback
          && Number(actionHeldAppearanceRollback.currentAppearanceProgress) > 0.0
        ) ? {
            fromSnapshot: ovizStatesClone(actionHeldAppearanceRollback.fromSnapshot, {}),
            targetSnapshot: ovizStatesClone(actionHeldAppearanceRollback.targetSnapshot, {}),
            targetWeight: clampRange(
              Number(actionHeldAppearanceRollback.currentAppearanceProgress) || 0.0,
              0.0,
              1.0,
            ),
          }
        : null;
        const traceCandidates = ovizTraceCandidatesForSnapshots(from, destination);
        const sourceTraceKeys = new Set(
          Array.from(ovizTraceCandidatesForFrameValue(ovizStateFrameValue(from)).keys())
        );
        const destinationTraceKeys = new Set(
          Array.from(ovizTraceCandidatesForFrameValue(ovizStateFrameValue(destination)).keys())
        );
        const liveTraceOpacity = new Map();
        traceCandidates.forEach((trace) => {
          liveTraceOpacity.set(
            trace.key,
            sourceTraceKeys.has(trace.key) && traceVisible(trace)
              ? clampRange(traceVisibilityOpacityMultiplier(trace), 0.0, 1.0)
              : 0.0,
          );
        });
        const liveSelectionWeights = ovizCaptureSelectionMembershipOpacityMap(traceCandidates);
        const liveSelectionDefaultWeight = ovizSelectionMembershipOpacity("", null);
        const previousSelectionTransition = ovizStateSelectionTransition || ovizHeldSelectionTransition;
        const previousStateTransition = ovizStateTransition;
        const destinationViewMode = String(
          (destination.global_controls || {}).camera_view_mode || "free"
        );
        const forceSkyLayerRetarget = Boolean(
          previousStateTransition
          && previousStateTransition.skyLayerTransitionStarted
          && destinationViewMode === "earth"
        );
        // Freeze the exact composed source membership before cancelling the
        // old transition.  Keeping four weighted masks prevents a rapid third
        // (or fourth) retarget from dropping an earlier live contribution.
        const sourceVolumeMaskComponents = ovizFreezeVolumeMaskSourceComponents(
          previousSelectionTransition,
          activeVolumeLassoSelectionMask()
        );
        const sourceVolumeMask = sourceVolumeMaskComponents[0]
          ? sourceVolumeMaskComponents[0].mask
          : null;
        const sourceVolumeSecondaryMask = sourceVolumeMaskComponents[1]
          ? sourceVolumeMaskComponents[1].mask
          : null;
        const sourceVolumeMaskBlend = sourceVolumeMaskComponents[1]
          ? sourceVolumeMaskComponents[1].weight
          : 0.0;
        if (
          ovizStateTransition
          || previousSelectionTransition
          || legendTransitionState
          || actionHeldTraceOpacityByKey
        ) {
          // Retarget phase detection must describe the frame that is actually
          // on screen, not only the last committed logical snapshot.  These
          // fields are transition-local and are never exported as authored state.
          from.__oviz_rendered_trace_opacity = Object.fromEntries(liveTraceOpacity);
          from.__oviz_rendered_selection_weights = Object.fromEntries(liveSelectionWeights);
          from.__oviz_rendered_volume_mask_blend = {
            active: Boolean(previousSelectionTransition),
            progress: previousSelectionTransition
              ? clampRange(Number(previousSelectionTransition.progress) || 0.0, 0.0, 1.0)
              : 0.0,
            has_source: Boolean(sourceVolumeMask),
            has_secondary: Boolean(sourceVolumeSecondaryMask),
            component_weights: sourceVolumeMaskComponents.map((component) => component.weight),
          };
        }
        const preservedSelectionMasks = new Set(
          sourceVolumeMaskComponents.map((component) => component.mask).filter(Boolean)
        );
        ovizCancelStateTransitionWithoutSnap("state-retarget", {
          restorePresentation: false,
          preserveSelectionMasks: preservedSelectionMasks,
          preserveSkyLayerTransition: forceSkyLayerRetarget,
        });
        ovizHeldSelectionTransition = null;
        const fromViewMode = String((from.global_controls || {}).camera_view_mode || "free");
        const toViewMode = destinationViewMode;
        const viewTransitionKind = fromViewMode === "earth" && toViewMode === "earth"
          ? "earth-to-earth"
          : (
            fromViewMode !== toViewMode && toViewMode === "earth"
              ? "enter-earth"
              : (
                fromViewMode !== toViewMode && fromViewMode === "earth"
                  ? "exit-earth"
                  : "generic"
              )
          );
        const nativeViewTransition = viewTransitionKind !== "generic";
        if (viewTransitionKind === "enter-earth" || viewTransitionKind === "exit-earth") {
          transitionSpec.duration_ms = Math.max(
            transitionSpec.duration_ms,
            toViewMode === "earth" ? 1320 : 1180
          );
        }
        const phasePlan = ovizBuildTransitionPhases(from, destination, transitionSpec);
        if (options.preserveActionRun !== true) {
          ovizCancelActionWithoutSnap("state-navigation");
        }
        playbackDirection = 0;
        lastPlaybackAdvanceTimestamp = null;
        setCameraAutoOrbitEnabled(false);
        const transitionFromFrame = ovizStateFrameValue(from);
        const transitionToFrame = ovizStateFrameValue(destination);
        ovizStateTimelineMotionActive = Math.abs(transitionToFrame - transitionFromFrame) > 1e-9;
        if (ovizStateTimelineMotionActive) {
          displayedFrameValue = clampFrameValue(transitionFromFrame);
          currentFrameIndex = clampFrameIndex(displayedFrameValue);
          updateTimelineUi(displayedFrameValue, frameTimeForValue(displayedFrameValue));
          updateTimelineMotionOpacity();
        }
        const now = performance.now();
        const sourceSkyLayers = typeof serializableSkyLayers === "function"
          ? serializableSkyLayers()
          : [];
        ovizStateTransitionSerial += 1;
        const transitionId = `state-view-${ovizStateTransitionSerial}`;
        const traceOpacity = new Map();
        traceCandidates.forEach((trace) => {
          const previous = liveTraceOpacity.has(trace.key)
            ? liveTraceOpacity.get(trace.key)
            : (traceVisible(trace) ? 1 : 0);
          traceOpacity.set(trace.key, {
            from: previous,
            to: destinationTraceKeys.has(trace.key) && ovizTraceTargetVisible(trace, destination) ? 1 : 0,
          });
        });
        ovizStateTransitionTraceOpacity = new Map();
        traceOpacity.forEach((value, key) => ovizStateTransitionTraceOpacity.set(key, value.from));
        let resolvePromise;
        let rejectPromise;
        const promise = new Promise((resolve, reject) => {
          resolvePromise = resolve;
          rejectPromise = reject;
        });
        ovizStateTransition = {
          targetId: target.id,
          targetIndex: target.index,
          targetName: target.name,
          targetSnapshot: destination,
          fromSnapshot: from,
          sourceAppearanceComposite,
          transitionSpec,
          phasePlan,
          startedAt: now,
          startedAtEpochMs: ovizTransitionEpochMs(now),
          traceOpacity,
          visualChanges: ovizTransitionVisualSignature(from) !== ovizTransitionVisualSignature(destination),
          nativeViewTransition,
          viewTransitionKind,
          transitionId,
          earthCameraTrack: viewTransitionKind === "earth-to-earth"
            ? ovizCreateEarthCameraTrack(destination)
            : null,
          nativeCameraTrack: (
            viewTransitionKind === "enter-earth" || viewTransitionKind === "exit-earth"
          ) ? ovizCreateNativeViewCameraTrack(destination) : null,
          skyBackgroundPromise: null,
          skyLayerPromise: Promise.resolve({ unchanged: true }),
          forceSkyLayerRetarget,
          fromViewMode,
          toViewMode,
          lastProgressEventAt: -Infinity,
          sceneRenderCount: 0,
          appearancePrepared: false,
          lassoReady: !destination.lasso_selection_mask,
          lassoLoadPromise: null,
          targetRuntimeLassoMask: null,
          selectionTransition: null,
          currentPhase: phasePlan.phases[0].name,
          currentAppearanceProgress: 0.0,
          animationFrameCount: 0,
          lastAnimationAt: now,
          maxFrameGapMs: 0,
          longFrameCount: 0,
          resolve: resolvePromise,
          reject: rejectPromise,
        };
        if (root && root.dataset) {
          root.dataset.stateTransitionKind = viewTransitionKind;
          root.dataset.stateTransitionProgress = "0";
          root.dataset.stateTransitionId = transitionId;
          root.dataset.stateTransitionPhase = phasePlan.phases[0].name;
          root.dataset.stateTransitionPhaseProgress = "0";
          root.dataset.stateTransitionEffectiveDurationMs = String(phasePlan.effectiveDurationMs);
        }
        ovizWriteTransitionDiagnostics(ovizStateTransition, {
          phase: phasePlan.phases[0].name,
          phaseProgress: 0.0,
          effectiveDurationMs: phasePlan.effectiveDurationMs,
          topologyPrepareCount: 0,
          retainedUpdateCount: 0,
        });
        const fromSelectionSignature = JSON.stringify({
          mode: from.current_selection_mode || "none",
          keys: from.selected_cluster_keys || [],
          mask: from.lasso_selection_mask || null,
          filter: from.lasso_selection_filter_enabled !== false,
        });
        const toSelectionSignature = JSON.stringify({
          mode: destination.current_selection_mode || "none",
          keys: destination.selected_cluster_keys || [],
          mask: destination.lasso_selection_mask || null,
          filter: destination.lasso_selection_filter_enabled !== false,
        });
        if (fromSelectionSignature !== toSelectionSignature || previousSelectionTransition) {
          const selectionTransition = {
            fromWeightByKey: liveSelectionWeights,
            fromDefaultWeight: liveSelectionDefaultWeight,
            toEndpoint: ovizSelectionEndpoint(destination),
            sourceMaskComponents: sourceVolumeMaskComponents,
            fromMask: sourceVolumeMask,
            fromSecondaryMask: sourceVolumeSecondaryMask,
            fromMaskBlend: sourceVolumeMaskBlend,
            toMask: null,
            progress: 0.0,
          };
          ovizStateTransition.selectionTransition = selectionTransition;
          ovizStateSelectionTransition = selectionTransition;
        }
        if (destination.lasso_selection_mask) {
          const transition = ovizStateTransition;
          transition.lassoLoadPromise = loadLassoSelectionMask(destination.lasso_selection_mask).then((mask) => {
            if (transition !== ovizStateTransition) {
              disposeSelectionMaskIfUnused(mask);
              return null;
            }
            transition.targetRuntimeLassoMask = mask;
            transition.lassoReady = true;
            if (
              transition.selectionTransition
              && destination.lasso_selection_filter_enabled !== false
            ) {
              transition.selectionTransition.toMask = mask;
            }
            return mask;
          });
        }
        if (nativeViewTransition) {
          const destinationCameraState = ovizStateDestinationCameraState(destination);
          if (viewTransitionKind === "enter-earth") {
            setMilkyWayModelOpacityScale(1.0);
            setSkyDomeViewOpacityScale(0.0, { force: false });
          } else if (viewTransitionKind === "exit-earth") {
            setMilkyWayModelOpacityScale(0.0);
            setSkyDomeViewOpacityScale(1.0, { force: false });
          }
          if (toViewMode === "earth") {
            if (viewTransitionKind === "enter-earth") {
              resetToSunReferenceFrameForSkyView();
              ovizStateTransition.nativeCameraTrack = ovizCreateNativeViewCameraTrack(destination);
              ovizPrepareEarthViewModeForStateTransition();
            }
            ovizPrepareDestinationSkyPresentation(destination, {
              postLayersToAladin: false,
              prepareOnly: true,
            });
            if (viewTransitionKind === "enter-earth") {
              const synchronizedStart = performance.now();
              ovizStateTransition.startedAt = synchronizedStart;
              ovizStateTransition.startedAtEpochMs = ovizTransitionEpochMs(synchronizedStart);
              ovizStateTransition.lastAnimationAt = synchronizedStart;
            }
            // The rendered Three.js camera is the only spatial source of truth.
            // updateSkyDomeBackgroundFrame() forwards that exact pose to Aladin
            // every frame, matching the smooth View-button path and preventing
            // a second iframe-local camera animation from drifting off the traces.
            ovizStateTransition.skyBackgroundPromise = null;
            // Layers are an appearance domain.  Keep the current registered
            // Aladin stack through the camera/Sky handoff and begin its own
            // crossfade only when the appearance phase starts.
            ovizStateTransition.sourceSkyLayers = sourceSkyLayers;
          }
        }
        ovizRenderStatesDrawer();
        ovizStateEvent("transition-start", {
          owner: "state",
          runId: transitionId,
          id: target.id,
          index: target.index + 1,
          name: target.name,
          durationMs: transitionSpec.duration_ms,
          effectiveDurationMs: phasePlan.effectiveDurationMs,
          phase: phasePlan.phases[0].name,
          easing: transitionSpec.easing,
        });
        return promise;
      }

      function ovizApplyTransitionNumericControls(from, to, progress) {
        const fromGlobal = from.global_controls || {};
        const toGlobal = to.global_controls || {};
        const pairs = [
          ["point_size_scale", (value) => { globalPointSizeScale = value; }],
          ["point_opacity_scale", (value) => { globalPointOpacityScale = value; }],
          ["point_glow_strength", (value) => { globalPointGlowStrength = value; }],
          ["sky_dome_opacity", (value) => { skyDomeSpec.opacity = value; }],
          ["sky_dome_hips_brightness", (value) => { skyDomeSpec.hips_brightness = value; }],
          ["sky_dome_hips_contrast", (value) => { skyDomeSpec.hips_contrast = value; }],
          ["sky_dome_hips_gamma", (value) => { skyDomeSpec.hips_gamma = value; }],
        ];
        pairs.forEach(([key, setter]) => {
          const a = Number(fromGlobal[key]);
          const b = Number(toGlobal[key]);
          if (Number.isFinite(a) && Number.isFinite(b)) {
            setter(ovizLerp(a, b, progress));
          }
        });
      }

      function ovizInterpolateColor(fromColor, toColor, progress) {
        try {
          const from = new THREE.Color(String(fromColor));
          const to = new THREE.Color(String(toColor));
          from.lerp(to, progress);
          return `#${from.getHexString()}`;
        } catch (_err) {
          return progress < 0.5 ? fromColor : toColor;
        }
      }

      function ovizApplyTransitionTraceStyles(from, to, progress) {
        const fromStyles = from.trace_style_state || {};
        const toStyles = to.trace_style_state || {};
        Object.keys(traceStyleStateByKey).forEach((key) => {
          const live = traceStyleStateByKey[key];
          const a = fromStyles[key] || live;
          const b = toStyles[key] || a;
          ["opacity", "sizeScale"].forEach((field) => {
            if (Number.isFinite(Number(a[field])) && Number.isFinite(Number(b[field]))) {
              live[field] = ovizLerp(a[field], b[field], progress);
            }
          });
          if (a.color && b.color) live.color = ovizInterpolateColor(a.color, b.color, progress);
          if (progress >= 0.5) {
            if (b.colorMode) live.colorMode = b.colorMode;
            if (b.colormap) live.colormap = b.colormap;
          }
        });
      }

      function ovizApplyTransitionVolumes(from, to, progress) {
        const fromVolumes = from.volume_state_by_key || {};
        const toVolumes = to.volume_state_by_key || {};
        Object.keys(volumeStateByKey).forEach((key) => {
          const live = volumeStateByKey[key];
          const a = fromVolumes[key] || live;
          const b = toVolumes[key] || a;
          const fromVisible = a.visible !== false;
          const toVisible = b.visible !== false;
          live.visible = fromVisible || toVisible;
          ["vmin", "vmax", "steps", "alphaCoef", "gradientStep"].forEach((field) => {
            if (Number.isFinite(Number(a[field])) && Number.isFinite(Number(b[field]))) {
              live[field] = ovizLerp(a[field], b[field], progress);
            }
          });
          const aOpacity = fromVisible ? Number(a.opacity) || 0 : 0;
          const bOpacity = toVisible ? Number(b.opacity) || 0 : 0;
          live.opacity = ovizLerp(aOpacity, bOpacity, progress);
          // Keep a volume resident for the whole crossfade.  The exact target
          // showAllTimes value is restored once all visual phases finish.
          live.showAllTimes = Boolean(a.showAllTimes || b.showAllTimes);
          if (progress >= 0.5) {
            if (b.stretch) live.stretch = b.stretch;
            if (b.colormap) live.colormap = b.colormap;
          }
        });
      }

      function ovizApplyTransitionPanelGeometry(from, to, progress) {
        const fromWidgets = from.widgets || {};
        const toWidgets = to.widgets || {};
        ["sky", "box_metrics", "age_kde", "cluster_filter", "dendrogram"].forEach((key) => {
          const panel = widgetPanelForKey(key);
          const a = fromWidgets[key] && fromWidgets[key].rect;
          const b = toWidgets[key] && toWidgets[key].rect;
          if (!panel || !a || !b) return;
          ["left", "top", "width", "height"].forEach((field) => {
            if (Number.isFinite(Number(a[field])) && Number.isFinite(Number(b[field]))) {
              panel.style[field] = `${ovizLerp(a[field], b[field], progress)}px`;
            }
          });
        });
      }

      function updateOvizStateTransition(now) {
        const transition = ovizStateTransition;
        if (!transition) {
          return false;
        }
        if (transition.finishing) {
          return true;
        }
        try {
          const frameGapMs = Math.max(0, now - transition.lastAnimationAt);
          transition.lastAnimationAt = now;
          transition.animationFrameCount += 1;
          transition.maxFrameGapMs = Math.max(transition.maxFrameGapMs, frameGapMs);
          if (frameGapMs > 50) transition.longFrameCount += 1;
          let elapsedMs = Math.max(0.0, now - transition.startedAt);
          let phaseState = ovizTransitionPhaseState(transition, elapsedMs);
          if (phaseState.name === "appearance" && !transition.lassoReady) {
            transition.startedAt += frameGapMs;
            transition.startedAtEpochMs += frameGapMs;
            elapsedMs = Math.max(0.0, now - transition.startedAt);
            phaseState = ovizTransitionPhaseState(transition, elapsedMs);
          }
          const raw = transition.phasePlan.effectiveDurationMs <= 0
            ? 1.0
            : clampRange(elapsedMs / transition.phasePlan.effectiveDurationMs, 0.0, 1.0);
          const cameraRaw = ovizTransitionPhaseProgress(transition, "camera", elapsedMs);
          const appearanceRaw = ovizTransitionPhaseProgress(transition, "appearance", elapsedMs);
          const timeRaw = ovizTransitionPhaseProgress(transition, "time", elapsedMs);
          const cameraProgress = ovizEasing(transition.transitionSpec.easing, cameraRaw);
          const appearanceProgress = ovizEasing(transition.transitionSpec.easing, appearanceRaw);
          const timeProgress = ovizEasing(transition.transitionSpec.easing, timeRaw);
          transition.currentPhase = phaseState.name;
          transition.currentAppearanceProgress = appearanceProgress;
          if (
            phaseState.name === "appearance"
            && !transition.skyLayerTransitionStarted
            && transition.sourceSkyLayers
          ) {
            transition.skyLayerTransitionStarted = true;
            ovizStartSkyLayerTransition(transition, transition.sourceSkyLayers);
          }
          const modeCameraRaw = transition.nativeViewTransition
            ? clampRange((cameraRaw - 0.20) / 0.60, 0, 1)
            : cameraRaw;
          const modeCameraProgress = ovizEasing(transition.transitionSpec.easing, modeCameraRaw);
          const from = transition.fromSnapshot;
          const to = transition.targetSnapshot;
          const fromPosition = ovizPointFrom(from, "position", camera.position);
          const toPosition = ovizPointFrom(to, "position", fromPosition);
          const fromTarget = ovizPointFrom(from, "target", controls.target);
          const toTarget = ovizPointFrom(to, "target", fromTarget);
          const fromUp = ovizPointFrom(from, "up", camera.up);
          const toUp = ovizPointFrom(to, "up", fromUp);
          if (transition.viewTransitionKind === "earth-to-earth") {
            ovizApplyEarthCameraTrack(transition.earthCameraTrack, cameraProgress);
          } else if (
            transition.viewTransitionKind === "enter-earth"
            || transition.viewTransitionKind === "exit-earth"
          ) {
            ovizApplyNativeViewCameraTrack(transition.nativeCameraTrack, modeCameraProgress);
          } else if (!transition.nativeViewTransition) {
            camera.position.set(
              ovizLerp(fromPosition.x, toPosition.x, cameraProgress),
              ovizLerp(fromPosition.y, toPosition.y, cameraProgress),
              ovizLerp(fromPosition.z, toPosition.z, cameraProgress)
            );
            controls.target.set(
              ovizLerp(fromTarget.x, toTarget.x, cameraProgress),
              ovizLerp(fromTarget.y, toTarget.y, cameraProgress),
              ovizLerp(fromTarget.z, toTarget.z, cameraProgress)
            );
            camera.up.set(
              ovizLerp(fromUp.x, toUp.x, cameraProgress),
              ovizLerp(fromUp.y, toUp.y, cameraProgress),
              ovizLerp(fromUp.z, toUp.z, cameraProgress)
            ).normalize();
          }
          const fromFov = Number((from.global_controls || {}).camera_fov);
          const toFov = Number((to.global_controls || {}).camera_fov);
          if (!transition.nativeViewTransition && Number.isFinite(fromFov) && Number.isFinite(toFov)) {
            camera.fov = ovizLerp(fromFov, toFov, cameraProgress);
            camera.updateProjectionMatrix();
          }
          if (!transition.nativeViewTransition) {
            const fromViewOffset = from.camera && (from.camera.view_offset || from.camera.viewOffset);
            const toViewOffset = to.camera && (to.camera.view_offset || to.camera.viewOffset);
            ovizApplyInterpolatedViewOffset(fromViewOffset, toViewOffset, cameraProgress);
          }
          camera.updateProjectionMatrix();
          camera.lookAt(controls.target);
          camera.updateMatrixWorld(true);
          if (root && root.dataset) {
            const liveDirection = new THREE.Vector3();
            camera.getWorldDirection(liveDirection);
            root.dataset.stateTransitionProgress = raw.toFixed(6);
            root.dataset.stateTransitionPhase = phaseState.name;
            root.dataset.stateTransitionPhaseProgress = phaseState.progress.toFixed(6);
            root.dataset.stateTransitionCamera = JSON.stringify({
              x: camera.position.x,
              y: camera.position.y,
              z: camera.position.z,
              tx: controls.target.x,
              ty: controls.target.y,
              tz: controls.target.z,
              dx: liveDirection.x,
              dy: liveDirection.y,
              dz: liveDirection.z,
              fov: Number(camera.fov),
            });
            if (typeof skyDomeBackgroundViewForCamera === "function") {
              const liveSkyView = skyDomeBackgroundViewForCamera();
              if (liveSkyView) {
                root.dataset.stateTransitionSkyView = JSON.stringify({
                  l: liveSkyView.l,
                  b: liveSkyView.b,
                  ra: liveSkyView.ra,
                  dec: liveSkyView.dec,
                  fovDeg: liveSkyView.fovDeg,
                  cameraFovDeg: liveSkyView.cameraFovDeg,
                });
              }
            }
          }
          ovizWriteTransitionDiagnostics(transition, {
            phase: phaseState.name,
            phaseProgress: phaseState.progress,
            effectiveDurationMs: transition.phasePlan.effectiveDurationMs,
            camera: {
              x: camera.position.x,
              y: camera.position.y,
              z: camera.position.z,
              tx: controls.target.x,
              ty: controls.target.y,
              tz: controls.target.z,
              fov: Number(camera.fov),
            },
            maxFrameGapMs: transition.maxFrameGapMs,
            longFrameCount: transition.longFrameCount,
            sceneRenderCount: transition.sceneRenderCount,
          });
          const fromFrame = ovizStateFrameValue(from);
          const toFrame = ovizStateFrameValue(to);
          const interpolatedFrameValue = ovizLerp(fromFrame, toFrame, timeProgress);
          ovizStateTimelineMotionActive = transition.phasePlan.changed.time && timeRaw > 0.0 && timeRaw < 1.0;
          ovizApplyTransitionNumericControls(from, to, appearanceProgress);
          ovizApplyTransitionTraceStyles(from, to, appearanceProgress);
          ovizApplyTransitionVolumes(from, to, appearanceProgress);
          ovizApplyTransitionPanelGeometry(from, to, appearanceProgress);
          transition.traceOpacity.forEach((value, key) => {
            ovizStateTransitionTraceOpacity.set(
              key,
              ovizLerp(value.from, value.to, appearanceProgress)
            );
          });
          if (transition.selectionTransition) {
            transition.selectionTransition.progress = appearanceProgress;
          }
          if (transition.viewTransitionKind === "enter-earth") {
            const milkyWayFade = 1.0 - ovizEasing(
              transition.transitionSpec.easing,
              clampRange(cameraRaw / 0.20, 0, 1)
            );
            const skyFade = ovizEasing(
              transition.transitionSpec.easing,
              clampRange((cameraRaw - 0.80) / 0.20, 0, 1)
            );
            setMilkyWayModelOpacityScale(milkyWayFade);
            setSkyDomeViewOpacityScale(skyFade, { force: false });
          } else if (transition.viewTransitionKind === "exit-earth") {
            const skyFade = 1.0 - ovizEasing(
              transition.transitionSpec.easing,
              clampRange(cameraRaw / 0.20, 0, 1)
            );
            const milkyWayFade = ovizEasing(
              transition.transitionSpec.easing,
              clampRange((cameraRaw - 0.80) / 0.20, 0, 1)
            );
            setSkyDomeViewOpacityScale(skyFade, { force: false });
            setMilkyWayModelOpacityScale(milkyWayFade);
          }
          displayedFrameValue = clampFrameValue(interpolatedFrameValue);
          currentFrameIndex = clampFrameIndex(displayedFrameValue);
          updateTimelineUi(displayedFrameValue, frameTimeForValue(displayedFrameValue));
          updateTimelineMotionOpacity();
          renderer.domElement.style.opacity = "1";
          if (appearanceRaw > 0.0 || timeRaw > 0.0) {
            renderInterpolatedFrameValue(displayedFrameValue, {
              updateWidgets: false,
              preserveCamera: true,
            });
            transition.sceneRenderCount += 1;
          }
          if (now - transition.lastProgressEventAt >= 100 || raw >= 1) {
            transition.lastProgressEventAt = now;
            ovizStateEvent("transition-progress", {
              owner: "state",
              runId: transition.transitionId,
              id: transition.targetId,
              index: transition.targetIndex + 1,
              progress: raw,
              phase: phaseState.name,
              phaseProgress: phaseState.progress,
              effectiveDurationMs: transition.phasePlan.effectiveDurationMs,
            });
          }
          if (raw >= 1) {
            ovizFinishStateTransition(transition);
          }
          return true;
        } catch (err) {
          ovizFailStateTransition(transition, err);
          return false;
        }
      }

      async function ovizFinishStateTransition(transition) {
        if (transition !== ovizStateTransition || transition.finishing) {
          return;
        }
        transition.finishing = true;
        const targetPosition = ovizPointFrom(transition.targetSnapshot, "position", camera.position);
        const targetControlsTarget = ovizPointFrom(transition.targetSnapshot, "target", controls.target);
        const preApplyCameraError = Math.sqrt(
          Math.pow(camera.position.x - targetPosition.x, 2)
          + Math.pow(camera.position.y - targetPosition.y, 2)
          + Math.pow(camera.position.z - targetPosition.z, 2)
          + Math.pow(controls.target.x - targetControlsTarget.x, 2)
          + Math.pow(controls.target.y - targetControlsTarget.y, 2)
          + Math.pow(controls.target.z - targetControlsTarget.z, 2)
        );
        const targetFov = Number((transition.targetSnapshot.global_controls || {}).camera_fov);
        const preApplyFovError = Number.isFinite(targetFov) ? Math.abs(Number(camera.fov) - targetFov) : 0;
        renderer.domElement.style.opacity = "1";
        try {
          if (transition.lassoLoadPromise) {
            await transition.lassoLoadPromise;
          }
          if (transition !== ovizStateTransition) {
            return;
          }
          if (transition.skyBackgroundPromise) {
            await Promise.race([
              transition.skyBackgroundPromise,
              new Promise((resolve) => window.setTimeout(() => resolve({ timedOut: true }), 500.0)),
            ]);
          }
          if (transition.skyLayerPromise) {
            await Promise.race([
              transition.skyLayerPromise,
              new Promise((resolve) => window.setTimeout(() => resolve({ timedOut: true }), 500.0)),
            ]);
          }
          if (transition !== ovizStateTransition) {
            return;
          }
          if (typeof finishSkyDomeBackgroundProgrammaticTransition === "function") {
            finishSkyDomeBackgroundProgrammaticTransition(transition.transitionId);
          }
          if (typeof cancelSkyViewTransitionAnimations === "function") {
            cancelSkyViewTransitionAnimations({
              cancelBackground: false,
              reason: "state-transition-complete",
            });
          }
          const finishedSelectionTransition = transition.selectionTransition;
          // The exact target render must not be multiplied by temporary
          // appearance/lasso tracks.  Previously these were cleared only
          // after ovizApplyStateImmediately(), leaving the already-rendered
          // destination frame hidden even though state fidelity was exact.
          ovizStateTransitionTraceOpacity = null;
          ovizStateSelectionTransition = null;
          ovizHeldSelectionTransition = null;
          transition.sceneRenderCount += 1;
          // The animated camera has already delivered the exact target view to
          // Aladin. Do not force the same center/FOV again here: a redundant
          // setFoV clears Aladin's canvas and can expose a black repaint frame.
          let appliedSnapshot = await ovizApplyStateImmediately(transition.targetSnapshot, {
            forceSkyBackground: false,
            postSkyLayersToAladin: false,
            runtimeLassoMask: transition.targetRuntimeLassoMask,
          });
          if (transition !== ovizStateTransition) {
            return;
          }
          // Let the exact destination survive a renderer frame before testing
          // fidelity.  Logical state alone cannot reveal a stale or empty plot.
          await new Promise((resolve) => window.requestAnimationFrame(() => resolve()));
          if (transition !== ovizStateTransition) {
            return;
          }
          appliedSnapshot = captureRuntimeState();
          let fidelityDifferences = ovizStateFidelityDifferences(
            transition.targetSnapshot,
            appliedSnapshot,
          );
          if (typeof ovizRenderedSceneFidelityDifferences === "function") {
            fidelityDifferences.push(
              ...ovizRenderedSceneFidelityDifferences(transition.targetSnapshot)
            );
          }
          if (fidelityDifferences.length) {
            appliedSnapshot = await ovizApplyStateImmediately(transition.targetSnapshot, {
              forceSkyBackground: false,
              postSkyLayersToAladin: false,
              runtimeLassoMask: transition.targetRuntimeLassoMask,
            });
            if (transition !== ovizStateTransition) {
              return;
            }
            await new Promise((resolve) => window.requestAnimationFrame(() => resolve()));
            if (transition !== ovizStateTransition) {
              return;
            }
            appliedSnapshot = captureRuntimeState();
            fidelityDifferences = ovizStateFidelityDifferences(
              transition.targetSnapshot,
              appliedSnapshot,
            );
            if (typeof ovizRenderedSceneFidelityDifferences === "function") {
              fidelityDifferences.push(
                ...ovizRenderedSceneFidelityDifferences(transition.targetSnapshot)
              );
            }
          }
          postSkyLayerStateToAladin();
          root.dataset.stateFidelity = JSON.stringify({
            exact: fidelityDifferences.length === 0,
            differences: fidelityDifferences,
          });
          if (fidelityDifferences.length) {
            throw new Error(`State fidelity check failed: ${fidelityDifferences.join(", ")}`);
          }
          ovizStateTimelineMotionActive = false;
          updateTimelineMotionOpacity();
          setMilkyWayModelOpacityScale(transition.toViewMode === "earth" ? 0.0 : 1.0);
          setSkyDomeViewOpacityScale(transition.toViewMode === "earth" ? 1.0 : 0.0, { force: false });
          if (transition !== ovizStateTransition) {
            return;
          }
          ovizStateTransition = null;
          if (finishedSelectionTransition) {
            const finishedMasks = new Set([
              finishedSelectionTransition.fromMask,
              finishedSelectionTransition.fromSecondaryMask,
              finishedSelectionTransition.toMask,
              transition.targetRuntimeLassoMask,
            ]);
            (finishedSelectionTransition.sourceMaskComponents || []).forEach((component) => {
              if (component && component.mask) finishedMasks.add(component.mask);
            });
            finishedMasks.forEach((mask) => disposeSelectionMaskIfUnused(mask));
          }
          if (root && root.dataset) {
            root.dataset.stateTransitionProgress = "1";
            root.dataset.stateTransitionId = "";
            root.dataset.stateTransitionPhase = "complete";
            root.dataset.stateTransitionPhaseProgress = "1";
          }
          ovizActiveStateId = transition.targetId;
          ovizStateDirty = false;
          const performanceMetrics = {
            durationMs: Math.max(0, performance.now() - transition.startedAt),
            effectiveDurationMs: transition.phasePlan.effectiveDurationMs,
            animationFrames: transition.animationFrameCount,
            sceneRenders: transition.sceneRenderCount,
            maxFrameGapMs: transition.maxFrameGapMs,
            longFrames: transition.longFrameCount,
            preApplyCameraError,
            preApplyFovError,
          };
          root.dataset.stateTransitionMetrics = JSON.stringify(performanceMetrics);
          ovizWriteTransitionDiagnostics(null, {
            completedRunId: transition.transitionId,
            performance: performanceMetrics,
          });
          ovizRenderStatesDrawer();
          ovizStateEvent("transition-complete", {
            owner: "state",
            runId: transition.transitionId,
            id: transition.targetId,
            index: transition.targetIndex + 1,
            name: transition.targetName,
            performance: performanceMetrics,
          });
          transition.resolve({ cancelled: false, id: transition.targetId, index: transition.targetIndex + 1 });
        } catch (err) {
          ovizFailStateTransition(transition, err);
        }
      }

      function ovizFailStateTransition(transition, err) {
        const failedSelectionTransition = transition && transition.selectionTransition;
        const failedMasks = new Set([transition && transition.targetRuntimeLassoMask]);
        if (failedSelectionTransition) {
          (failedSelectionTransition.sourceMaskComponents || []).forEach((component) => {
            if (component && component.mask) failedMasks.add(component.mask);
          });
          [
            failedSelectionTransition.fromMask,
            failedSelectionTransition.fromSecondaryMask,
            failedSelectionTransition.toMask,
          ].forEach((mask) => failedMasks.add(mask));
        }
        // Async completion from an abandoned run must never reset the active
        // camera/Sky presentation or emit an error for its replacement.
        if (transition !== ovizStateTransition) {
          failedMasks.forEach((mask) => disposeSelectionMaskIfUnused(mask));
          return false;
        }
        ovizStateTransition = null;
        ovizStateTransitionTraceOpacity = null;
        ovizStateSelectionTransition = null;
        ovizHeldSelectionTransition = null;
        failedMasks.forEach((mask) => disposeSelectionMaskIfUnused(mask));
        if (typeof cancelSkyViewTransitionAnimations === "function") {
          cancelSkyViewTransitionAnimations({ reason: "state-transition-error" });
        }
        ovizStateTimelineMotionActive = false;
        updateTimelineMotionOpacity();
        renderer.domElement.style.opacity = "1";
        setMilkyWayModelOpacityScale(cameraViewMode === "earth" ? 0.0 : 1.0);
        setSkyDomeViewOpacityScale(cameraViewMode === "earth" ? 1.0 : 0.0, { force: false });
        ovizWriteTransitionDiagnostics(null, {
          failedRunId: transition && transition.transitionId,
          error: String(err && err.message || err),
        });
        ovizRenderStatesDrawer();
        ovizStateEvent("transition-error", {
          owner: "state",
          runId: transition && transition.transitionId,
          id: transition && transition.targetId,
          message: String(err && err.message || err),
        });
        if (transition && typeof transition.reject === "function") {
          transition.reject(err);
        }
        return true;
      }

      function ovizGoToState(idOrIndex) {
        if (!ovizStateControllerReady) {
          return Promise.reject(new Error("Oviz States controller is not ready."));
        }
        try {
          return ovizBeginStateTransition(ovizStateTargetFor(idOrIndex));
        } catch (err) {
          return Promise.reject(err);
        }
      }

      function ovizStatesNext() {
        const activeIndex = ovizStateTransition
          ? Number(ovizStateTransition.targetIndex)
          : (
            ovizActiveStateId === null
              ? -1
              : ovizStatesProject.items.findIndex((item) => item.id === ovizActiveStateId)
          );
        if (activeIndex >= ovizStatesProject.items.length - 1) {
          return Promise.resolve({ boundary: true, id: ovizActiveStateId });
        }
        return ovizGoToState(activeIndex + 2);
      }

      function ovizStatesPrevious() {
        const activeIndex = ovizStateTransition
          ? Number(ovizStateTransition.targetIndex)
          : (
            ovizActiveStateId === null
              ? -1
              : ovizStatesProject.items.findIndex((item) => item.id === ovizActiveStateId)
          );
        if (activeIndex < 0) {
          return Promise.resolve({ boundary: true, id: null });
        }
        return activeIndex === 0 ? ovizGoToState("original") : ovizGoToState(activeIndex);
      }

      function ovizAssertEditable() {
        if (ovizStatesMode !== "edit") {
          throw new Error("Switch Oviz States to Edit mode first.");
        }
        if (ovizStateTransition) {
          throw new Error("Add and Update are unavailable during a transition.");
        }
      }

      function ovizAddState(options = {}) {
        ovizAssertEditable();
        const item = ovizNormalizeStateRecord({
          id: options.id,
          name: options.name || `State ${ovizStatesProject.items.length + 1}`,
          transition: options.transition || null,
          snapshot: options.snapshot || captureRuntimeState(),
        }, ovizStatesProject.items.length, new Set(ovizStatesProject.items.map((state) => state.id)));
        ovizStatesProject.items.push(item);
        ovizActiveStateId = item.id;
        ovizStatesChanged("add");
        return ovizStatesClone(item, {});
      }

      function ovizUpdateState(idOrIndex, options = {}) {
        ovizAssertEditable();
        const index = ovizStateIndexFor(idOrIndex ?? ovizActiveStateId);
        if (index < 0) throw new Error("Select a saved state to update.");
        const item = ovizStatesProject.items[index];
        item.snapshot = ovizStatesClone(options.snapshot || captureRuntimeState(), {});
        if (options.transition !== undefined) item.transition = options.transition ? ovizNormalizeTransition(options.transition) : null;
        if (options.name !== undefined) item.name = String(options.name).trim() || item.name;
        ovizActiveStateId = item.id;
        ovizStatesChanged("update");
        return ovizStatesClone(item, {});
      }

      function ovizRenameState(idOrIndex, name) {
        ovizAssertEditable();
        const index = ovizStateIndexFor(idOrIndex);
        if (index < 0) throw new Error("Unknown state.");
        ovizStatesProject.items[index].name = String(name || "").trim() || ovizStatesProject.items[index].name;
        ovizStatesChanged("rename");
        return ovizStatesClone(ovizStatesProject.items[index], {});
      }

      function ovizDuplicateState(idOrIndex, options = {}) {
        ovizAssertEditable();
        const index = ovizStateIndexFor(idOrIndex);
        if (index < 0) throw new Error("Unknown state.");
        const source = ovizStatesProject.items[index];
        const copy = ovizNormalizeStateRecord({
          name: options.name || `${source.name} copy`,
          transition: source.transition,
          snapshot: source.snapshot,
        }, index + 1, new Set(ovizStatesProject.items.map((item) => item.id)));
        ovizStatesProject.items.splice(index + 1, 0, copy);
        ovizStatesChanged("duplicate");
        return ovizStatesClone(copy, {});
      }

      function ovizMoveState(idOrIndex, destinationIndex) {
        ovizAssertEditable();
        const index = ovizStateIndexFor(idOrIndex);
        if (index < 0) throw new Error("Unknown state.");
        const destination = clampRange(Math.floor(Number(destinationIndex)) - 1, 0, Math.max(ovizStatesProject.items.length - 1, 0));
        const [item] = ovizStatesProject.items.splice(index, 1);
        ovizStatesProject.items.splice(destination, 0, item);
        ovizStatesChanged("move");
        return ovizStatesList();
      }

      function ovizRemoveState(idOrIndex) {
        ovizAssertEditable();
        const index = ovizStateIndexFor(idOrIndex);
        if (index < 0) throw new Error("Unknown state.");
        const [removed] = ovizStatesProject.items.splice(index, 1);
        if (ovizActiveStateId === removed.id) ovizActiveStateId = null;
        ovizStatesChanged("remove");
        return ovizStatesClone(removed, {});
      }

      function ovizSetStatesMode(mode) {
        const normalized = String(mode).toLowerCase();
        if (!['edit', 'present'].includes(normalized)) throw new Error("Mode must be 'edit' or 'present'.");
        ovizStatesMode = normalized;
        ovizRenderStatesDrawer();
        return ovizStatesMode;
      }

      async function ovizWriteHtmlFile(htmlText, suggestedName) {
        if (typeof window.showSaveFilePicker === "function") {
          try {
            const handle = await window.showSaveFilePicker({
              suggestedName,
              types: [{ description: "HTML file", accept: { "text/html": [".html"] } }],
            });
            const writable = await handle.createWritable();
            await writable.write(htmlText);
            await writable.close();
            return { saved: true, filename: suggestedName };
          } catch (err) {
            if (err && err.name === "AbortError") return { saved: false, cancelled: true };
          }
        }
        const blob = new Blob([htmlText], { type: "text/html;charset=utf-8" });
        const url = URL.createObjectURL(blob);
        const link = document.createElement("a");
        link.href = url;
        link.download = suggestedName;
        document.body.appendChild(link);
        link.click();
        link.remove();
        window.setTimeout(() => URL.revokeObjectURL(url), 1000);
        return { saved: true, filename: suggestedName };
      }

      async function ovizExportStatesHtml(options = {}) {
        const compact = await ovizCompactProjectForStorage(ovizStatesPublicProject());
        compact.default_mode = "present";
        compact.synchronized_revision = compact.revision;
        const exportSceneSpec = safeJsonClone(sceneSpec, {});
        exportSceneSpec.initial_state = ovizStatesClone(ovizOriginalSceneInitialState, {});
        exportSceneSpec.states = compact;
        exportSceneSpec.width = Math.max(root.clientWidth || sceneSpec.width || 900, 1);
        exportSceneSpec.height = Math.max(root.clientHeight || sceneSpec.height || 700, 1);
        const html = await buildExportHtml(exportSceneSpec);
        if (options.download === false) return html;
        const result = await ovizWriteHtmlFile(html, options.filename || `${slugifyFilename(sceneSpec.title || "oviz")}-states.html`);
        if (result.saved) {
          ovizStatesProject.synchronized_revision = ovizStatesProject.revision;
          ovizStateDirty = false;
          await ovizSaveDraftNow();
          ovizRenderStatesDrawer();
        }
        return result;
      }

      function ovizOpenDraftDb() {
        return new Promise((resolve, reject) => {
          if (!window.indexedDB) {
            resolve(null);
            return;
          }
          const request = window.indexedDB.open(OVIZ_STATES_DB_NAME, 1);
          request.onupgradeneeded = () => {
            const db = request.result;
            if (!db.objectStoreNames.contains(OVIZ_STATES_DB_STORE)) {
              db.createObjectStore(OVIZ_STATES_DB_STORE, { keyPath: "project_id" });
            }
          };
          request.onsuccess = () => resolve(request.result);
          request.onerror = () => reject(request.error);
        });
      }

      async function ovizReadDraft() {
        const db = await ovizOpenDraftDb();
        if (!db) return null;
        return new Promise((resolve) => {
          const transaction = db.transaction(OVIZ_STATES_DB_STORE, "readonly");
          const request = transaction.objectStore(OVIZ_STATES_DB_STORE).get(ovizStatesProject.project_id);
          request.onsuccess = () => resolve(request.result || null);
          request.onerror = () => resolve(null);
          transaction.oncomplete = () => db.close();
        });
      }

      async function ovizSaveDraftNow() {
        if (!ovizStatesProject) return;
        const db = await ovizOpenDraftDb();
        if (!db) return;
        const compact = await ovizCompactProjectForStorage(ovizStatesPublicProject());
        compact.saved_at = Date.now();
        return new Promise((resolve) => {
          const transaction = db.transaction(OVIZ_STATES_DB_STORE, "readwrite");
          transaction.objectStore(OVIZ_STATES_DB_STORE).put(compact);
          transaction.oncomplete = () => { db.close(); resolve(); };
          transaction.onerror = () => { db.close(); resolve(); };
        });
      }

      function ovizScheduleDraftSave() {
        if (ovizStatesDraftTimer) window.clearTimeout(ovizStatesDraftTimer);
        ovizStatesDraftTimer = window.setTimeout(() => {
          ovizStatesDraftTimer = null;
          ovizSaveDraftNow();
        }, 250);
      }

      function ovizCollectAladinSources() {
        const sources = new Set();
        ovizStatesProject.items.forEach((item) => {
          const snapshot = ovizHydrateAssets(item.snapshot);
          (Array.isArray(snapshot.sky_layers) ? snapshot.sky_layers : []).forEach((layer) => {
            const source = layer && (layer.survey || layer.hips || layer.url || layer.source || layer.key);
            if (source) sources.add(String(source));
          });
        });
        return Array.from(sources).filter(Boolean);
      }

      function ovizPreloadAladinSource(source) {
        return new Promise((resolve) => {
          const iframe = document.createElement("iframe");
          iframe.setAttribute("aria-hidden", "true");
          iframe.style.cssText = "position:absolute;width:1px;height:1px;opacity:0;pointer-events:none;left:-10000px";
          let settled = false;
          const finish = (ok) => {
            if (settled) return;
            settled = true;
            window.clearTimeout(timer);
            iframe.remove();
            resolve(ok);
          };
          iframe.onload = () => window.setTimeout(() => finish(true), 350);
          iframe.onerror = () => finish(false);
          const timer = window.setTimeout(() => finish(false), 8000);
          try {
            iframe.srcdoc = buildAladinSrcdoc([], [], "overview", null, false, { survey: source });
            root.appendChild(iframe);
          } catch (_err) {
            finish(false);
          }
        });
      }

      async function ovizPreloadAllAladinSources() {
        const sources = ovizCollectAladinSources();
        ovizStatesPreloadStatus = { loaded: 0, total: sources.length, failed: [] };
        ovizRenderStatesDrawer();
        await Promise.all(sources.map(async (source) => {
          const ok = await ovizPreloadAladinSource(source);
          ovizStatesPreloadStatus.loaded += 1;
          if (!ok) ovizStatesPreloadStatus.failed.push(source);
          ovizRenderStatesDrawer();
        }));
        if (ovizStatesPreloadStatus.failed.length) {
          ovizStatesProject.items.forEach((item) => {
            const serialized = JSON.stringify(item.snapshot || {});
            item.degraded = ovizStatesPreloadStatus.failed.some((source) => serialized.includes(source));
          });
        }
      }

      function ovizInstallStatesStyles() {
        if (document.getElementById("oviz-states-styles")) return;
        const style = document.createElement("style");
        style.id = "oviz-states-styles";
        style.textContent = `
          .oviz-states-shell{position:relative;display:flex;align-items:center}
          .oviz-states-toggle{white-space:nowrap}
          .oviz-states-drawer{position:absolute;right:0;top:calc(100% + 8px);width:min(390px,92vw);max-height:min(72vh,680px);overflow:auto;padding:10px;background:var(--oviz-panel-bg,#fff);color:var(--oviz-text,#222);border:1px solid rgba(127,127,127,.35);border-radius:10px;box-shadow:0 12px 35px rgba(0,0,0,.22);z-index:80;display:none}
          .oviz-states-shell[data-open=true] .oviz-states-drawer{display:block}
          .oviz-states-head,.oviz-states-nav,.oviz-states-editbar{display:flex;align-items:center;gap:6px;margin-bottom:8px}
          .oviz-states-head strong{flex:1}.oviz-states-head small{opacity:.7}
          .oviz-states-drawer button,.oviz-states-drawer select,.oviz-states-drawer input{font:inherit}
          .oviz-states-nav button{flex:1}.oviz-states-editbar{flex-wrap:wrap}
          .oviz-states-row{display:grid;grid-template-columns:28px minmax(0,1fr) auto;gap:6px;align-items:center;padding:5px 4px;border-radius:6px}
          .oviz-states-row:hover{background:rgba(127,127,127,.1)}.oviz-states-row[data-active=true]{background:rgba(70,130,255,.16)}
          .oviz-states-row-name{border:0;background:transparent;text-align:left;overflow:hidden;text-overflow:ellipsis;white-space:nowrap;color:inherit}
          .oviz-states-row-actions{display:flex;gap:3px}.oviz-states-row-actions button{padding:2px 5px}
          .oviz-states-status{min-height:1.2em;font-size:11px;opacity:.75}.oviz-states-recovery{padding:7px;margin-bottom:8px;background:rgba(255,174,0,.18);border-radius:6px}
          .oviz-states-mode-edit .oviz-states-present-only{display:none}.oviz-states-mode-present .oviz-states-edit-only{display:none}
          @media(max-width:700px){.oviz-states-drawer{position:fixed;inset:52px 8px 8px auto;width:min(390px,calc(100vw - 16px));max-height:none}}
        `;
        document.head.appendChild(style);
      }

      function ovizMakeButton(label, action, className = "") {
        const button = document.createElement("button");
        button.type = "button";
        button.textContent = label;
        button.className = className;
        button.addEventListener("click", (event) => {
          event.stopPropagation();
          try {
            const result = action();
            if (result && typeof result.catch === "function") result.catch((err) => ovizSetStatesStatus(err.message, true));
          } catch (err) {
            ovizSetStatesStatus(err.message, true);
          }
        });
        return button;
      }

      function ovizSetStatesStatus(message, error = false) {
        if (!ovizStatesStatusEl) return;
        ovizStatesStatusEl.textContent = String(message || "");
        ovizStatesStatusEl.style.color = error ? "#d54b4b" : "";
      }

      function ovizBuildStatesDrawer() {
        if (!widgetMenuEl || minimalModeEnabled) return;
        ovizInstallStatesStyles();
        ovizStatesShellEl = document.createElement("div");
        ovizStatesShellEl.className = "oviz-states-shell";
        ovizStatesShellEl.dataset.open = "false";
        const toggle = ovizMakeButton("States ▸", () => {
          const open = ovizStatesShellEl.dataset.open !== "true";
          ovizStatesShellEl.dataset.open = open ? "true" : "false";
          toggle.textContent = open ? "States ▾" : "States ▸";
        }, "oviz-states-toggle");
        ovizStatesDrawerEl = document.createElement("div");
        ovizStatesDrawerEl.className = "oviz-states-drawer";
        ovizStatesDrawerEl.addEventListener("pointerdown", (event) => event.stopPropagation());
        ovizStatesShellEl.append(toggle, ovizStatesDrawerEl);
        widgetMenuEl.prepend(ovizStatesShellEl);
        ovizRenderStatesDrawer();
      }

      function ovizRenderStatesDrawer() {
        if (!ovizStatesDrawerEl || !ovizStatesProject) return;
        ovizStatesDrawerEl.className = `oviz-states-drawer oviz-states-mode-${ovizStatesMode}`;
        ovizStatesDrawerEl.innerHTML = "";
        const head = document.createElement("div"); head.className = "oviz-states-head";
        const title = document.createElement("strong"); title.textContent = "States";
        const count = document.createElement("small"); count.textContent = `${ovizStatesProject.items.length} saved`;
        const modeButton = ovizMakeButton(ovizStatesMode === "edit" ? "Present" : "Edit", () => ovizSetStatesMode(ovizStatesMode === "edit" ? "present" : "edit"));
        head.append(title, count, modeButton); ovizStatesDrawerEl.append(head);
        const nav = document.createElement("div"); nav.className = "oviz-states-nav";
        nav.append(
          ovizMakeButton("Previous", ovizStatesPrevious),
          ovizMakeButton("Original", () => ovizGoToState("original")),
          ovizMakeButton("Next", ovizStatesNext)
        );
        ovizStatesDrawerEl.append(nav);
        ovizStatesRecoveryEl = document.createElement("div");
        ovizStatesRecoveryEl.className = "oviz-states-recovery";
        ovizStatesRecoveryEl.hidden = true;
        ovizStatesDrawerEl.append(ovizStatesRecoveryEl);
        const editbar = document.createElement("div"); editbar.className = "oviz-states-editbar edit-only oviz-states-edit-only";
        editbar.append(
          ovizMakeButton("Add", () => ovizAddState()),
          ovizMakeButton("Update", () => ovizUpdateState(ovizActiveStateId)),
          ovizMakeButton("Export HTML", () => ovizExportStatesHtml()),
          ovizMakeButton("Current view only", () => saveSceneStateToHtml())
        );
        if (ovizStateTransition) editbar.querySelectorAll("button").forEach((button, index) => { if (index < 2) button.disabled = true; });
        ovizStatesDrawerEl.append(editbar);
        ovizStatesRowsEl = document.createElement("div");
        ovizStatesProject.items.forEach((item, index) => {
          const row = document.createElement("div"); row.className = "oviz-states-row"; row.dataset.active = item.id === ovizActiveStateId ? "true" : "false";
          const number = document.createElement("span"); number.textContent = String(index + 1);
          const name = ovizMakeButton(`${item.name}${item.degraded ? " ⚠" : ""}`, () => ovizGoToState(item.id), "oviz-states-row-name");
          const actions = document.createElement("span"); actions.className = "oviz-states-row-actions oviz-states-edit-only";
          actions.append(
            ovizMakeButton("↑", () => ovizMoveState(item.id, Math.max(1, index))),
            ovizMakeButton("↓", () => ovizMoveState(item.id, Math.min(ovizStatesProject.items.length, index + 2))),
            ovizMakeButton("✎", () => { const value = window.prompt("State name", item.name); if (value !== null) ovizRenameState(item.id, value); }),
            ovizMakeButton("⧉", () => ovizDuplicateState(item.id)),
            ovizMakeButton("×", () => ovizRemoveState(item.id))
          );
          row.append(number, name, actions); ovizStatesRowsEl.append(row);
        });
        ovizStatesDrawerEl.append(ovizStatesRowsEl);
        ovizStatesStatusEl = document.createElement("div"); ovizStatesStatusEl.className = "oviz-states-status";
        const preload = ovizStatesPreloadStatus.total
          ? ` · assets ${ovizStatesPreloadStatus.loaded}/${ovizStatesPreloadStatus.total}${ovizStatesPreloadStatus.failed.length ? " degraded" : ""}`
          : "";
        ovizStatesStatusEl.textContent = `${ovizStateTransition ? "Transitioning" : (ovizStateDirty ? "Unsaved changes" : "Ready")}${preload}`;
        ovizStatesDrawerEl.append(ovizStatesStatusEl);
      }

      function ovizShowDraftRecovery(draft) {
        if (!ovizStatesRecoveryEl) return;
        ovizStatesRecoveryEl.hidden = false;
        ovizStatesRecoveryEl.textContent = "A newer local draft is available. ";
        ovizStatesRecoveryEl.append(
          ovizMakeButton("Restore", () => {
            ovizStatesProject = ovizNormalizeStatesProject(draft);
            ovizStatesMode = "edit";
            ovizStateDirty = true;
            ovizRenderStatesDrawer();
            ovizStateEvent("states-changed", { reason: "restore-draft", states: ovizStatesList() });
          }),
          ovizMakeButton("Discard", async () => {
            draft.revision = -1;
            const db = await ovizOpenDraftDb();
            if (db) {
              const transaction = db.transaction(OVIZ_STATES_DB_STORE, "readwrite");
              transaction.objectStore(OVIZ_STATES_DB_STORE).delete(ovizStatesProject.project_id);
              transaction.oncomplete = () => db.close();
            }
            ovizRenderStatesDrawer();
          })
        );
      }

      function ovizInstallPublicApi() {
        const api = {
          list: ovizStatesList,
          goTo: ovizGoToState,
          next: ovizStatesNext,
          previous: ovizStatesPrevious,
          original: () => ovizGoToState("original"),
          add: ovizAddState,
          update: ovizUpdateState,
          rename: ovizRenameState,
          duplicate: ovizDuplicateState,
          move: ovizMoveState,
          remove: ovizRemoveState,
          setMode: ovizSetStatesMode,
          capture: () => captureRuntimeState(),
          exportHtml: ovizExportStatesHtml,
        };
        window.Oviz = window.Oviz || {};
        const registry = window.Oviz.__viewers || new Map();
        window.Oviz.__viewers = registry;
        registry.set(root.id, { root, states: api });
        window.Oviz.get = (rootId) => registry.get(String(rootId));
      }

      function ovizPostMessageHandler(event) {
        const data = event.data;
        if (data && data.type === "oviz-aladin-sky-layer-transition-complete") {
          const transitionId = String(data.transitionId || "");
          const resolve = ovizSkyLayerTransitionWaiters.get(transitionId);
          if (resolve) {
            ovizSkyLayerTransitionWaiters.delete(transitionId);
            resolve(data);
          }
          return;
        }
        if (!data || data.source !== "oviz-command" || (data.rootId && data.rootId !== root.id)) return;
        const requestId = data.requestId || null;
        const command = String(data.command || "");
        const api = window.Oviz && window.Oviz.get(root.id) && window.Oviz.get(root.id).states;
        if (!api || typeof api[command] !== "function") return;
        Promise.resolve().then(() => api[command](...(Array.isArray(data.args) ? data.args : []))).then((result) => {
          event.source && event.source.postMessage({ source: "oviz-result", rootId: root.id, requestId, command, ok: true, result }, "*");
        }).catch((err) => {
          event.source && event.source.postMessage({ source: "oviz-result", rootId: root.id, requestId, command, ok: false, error: String(err && err.message || err) }, "*");
        });
      }

      async function initializeOvizStates() {
        ovizOriginalSceneInitialState = ovizStatesClone(sceneSpec.initial_state, {});
        ovizOriginalRuntimeState = captureRuntimeState();
        ovizOriginalCameraResetState = {
          position: initialCameraState.position.clone(),
          target: initialCameraState.target.clone(),
          up: initialCameraState.up.clone(),
          fov: Number(initialCameraState.fov),
          viewOffset: ovizStatesClone(initialCameraState.viewOffset, { x: 0, y: 0 }),
        };
        ovizStatesProject = ovizNormalizeStatesProject(sceneSpec.states);
        ovizStatesMode = ovizStatesProject.embedded ? ovizStatesProject.default_mode : "edit";
        ovizBuildStatesDrawer();
        postSkyLayerStateToAladin();
        ovizInstallPublicApi();
        window.addEventListener("message", ovizPostMessageHandler);
        root.addEventListener("pointerdown", (event) => {
          if (ovizStateTransition && (!ovizStatesShellEl || !ovizStatesShellEl.contains(event.target))) {
            ovizCancelStateTransitionWithoutSnap("user-interaction", { restorePresentation: true });
          }
        }, { capture: true });
        const draft = await ovizReadDraft();
        if (draft && Number(draft.revision) > Number(ovizStatesProject.revision)) {
          ovizShowDraftRecovery(draft);
        }
        await ovizPreloadAllAladinSources();
        ovizStateControllerReady = true;
        ovizRenderStatesDrawer();
        ovizStateEvent("states-ready", {
          mode: ovizStatesMode,
          states: ovizStatesList(),
          degraded: ovizStatesPreloadStatus.failed.length > 0,
          preload: ovizStatesClone(ovizStatesPreloadStatus, {}),
        });
      }
""".strip()
