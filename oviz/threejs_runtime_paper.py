"""Embedded JavaScript for the scroll-synchronized paper reader."""

from __future__ import annotations


THREEJS_PAPER_RUNTIME_JS = """
      let ovizPaperProject = null;
      let ovizPaperPanelEl = null;
      let ovizPaperScrollEl = null;
      let ovizPaperBodyEl = null;
      let ovizPaperRailEl = null;
      let ovizPaperSyncButtonEl = null;
      let ovizPaperToggleButtonEl = null;
      let ovizPaperStatusEl = null;
      let ovizPaperLightboxEl = null;
      let ovizPaperLightboxRestoreFocusEl = null;
      let ovizPaperHeadResizeObserver = null;
      let ovizPaperBodyResizeObserver = null;
      let ovizPaperModeObserver = null;
      let ovizPaperRailLayoutFrame = 0;
      let ovizPaperResizerEl = null;
      let ovizPaperCollapseButtonEl = null;
      let ovizPaperExpandTabEl = null;
      let ovizPaperWidthFraction = 0.42;
      let ovizPaperOffsetTweenFrame = 0;
      let ovizPaperAnchorEntries = [];
      let ovizPaperActiveAnchorId = "";
      let ovizPaperSyncMode = "on";
      let ovizPaperPausedAnchorId = "";
      let ovizPaperScrollFrame = 0;
      let ovizPaperPendingJumpTarget = null;
      let ovizPaperJumpDebounceTimer = 0;
      let ovizPaperJumpGeneration = 0;
      let ovizPaperScrollSyncSuppressed = false;
      let ovizPaperScrollSettleTimer = 0;
      let ovizPaperKatexRequested = false;

      function ovizPaperAssetValue(value) {
        if (
          value
          && typeof value === "object"
          && typeof value.__oviz_asset_ref__ === "string"
          && ovizPaperProject
          && ovizPaperProject.assets
        ) {
          return String(ovizPaperProject.assets[value.__oviz_asset_ref__] || "");
        }
        return typeof value === "string" ? value : "";
      }

      function ovizPaperSetDataset(key, value) {
        if (root && root.dataset) {
          root.dataset[key] = String(value);
        }
      }

      function ovizPaperSetStatus(status, message = "") {
        const nextStatus = String(status || "ready");
        ovizPaperSetDataset("paperStatus", nextStatus);
        if (ovizPaperStatusEl) {
          ovizPaperStatusEl.textContent = String(message || "");
          ovizPaperStatusEl.dataset.status = nextStatus;
        }
      }

      async function ovizPaperFlushStateJump(generation) {
        if (generation !== ovizPaperJumpGeneration) {
          return null;
        }
        if (typeof ovizWaitForStatesControllerReady === "function") {
          try {
            await ovizWaitForStatesControllerReady();
          } catch (_readyError) {
            if (generation === ovizPaperJumpGeneration) {
              ovizPaperPendingJumpTarget = null;
              ovizPaperSetDataset("paperJumpPending", "false");
              ovizPaperSetStatus("error", "Scene states are unavailable");
            }
            return null;
          }
        }
        // Waiting for readiness may outlive the scroll position that queued
        // this request.  Only the latest generation may touch the viewer.
        if (generation !== ovizPaperJumpGeneration) {
          return null;
        }
        const next = ovizPaperPendingJumpTarget;
        ovizPaperPendingJumpTarget = null;
        if (next === null || next === undefined || next === "") {
          ovizPaperSetDataset("paperJumpPending", "false");
          ovizPaperSetStatus("ready", "Scene ready");
          return null;
        }
        const activeTransition = (
          typeof ovizStateTransition !== "undefined" ? ovizStateTransition : null
        );
        if (
          !activeTransition
          && typeof ovizActiveStateId !== "undefined"
          && next !== "original"
          && String(ovizActiveStateId) === String(next)
        ) {
          ovizPaperSetDataset("paperJumpPending", "false");
          ovizPaperSetStatus("ready", "Scene ready");
          return null;
        }
        ovizPaperClearHoverCue();
        try {
          // ovizGoToState already retargets from the exact live rendered
          // frame.  Calling it now (rather than awaiting the old promise)
          // makes fast scroll input last-wins without a queued transition lag.
          const result = await ovizGoToState(next);
          if (generation === ovizPaperJumpGeneration) {
            ovizPaperSetDataset("paperJumpPending", "false");
            ovizPaperSetStatus("ready", "Scene updated");
          }
          return result;
        } catch (_jumpError) {
          if (generation === ovizPaperJumpGeneration) {
            ovizPaperSetDataset("paperJumpPending", "false");
            ovizPaperSetStatus("error", "Scene update failed");
          }
          return null;
        }
      }

      function ovizPaperQueueStateJump(target) {
        if (target === null || target === undefined || target === "") {
          return null;
        }
        ovizPaperPendingJumpTarget = target;
        ovizPaperJumpGeneration += 1;
        const generation = ovizPaperJumpGeneration;
        if (ovizPaperJumpDebounceTimer) {
          window.clearTimeout(ovizPaperJumpDebounceTimer);
        }
        ovizPaperSetDataset("paperJumpPending", "true");
        ovizPaperSetStatus("syncing", "Updating scene…");
        ovizPaperJumpDebounceTimer = window.setTimeout(() => {
          ovizPaperJumpDebounceTimer = 0;
          void ovizPaperFlushStateJump(generation);
        }, 110);
        return null;
      }

      let ovizPaperHoverCueEl = null;
      let ovizPaperHoverHeldOpacity = null;
      let ovizPaperHoverTracesActive = false;
      let ovizPaperHoverVolumeRestore = null;

      function ovizPaperSampleTraces() {
        const frame = sceneSpec && Array.isArray(sceneSpec.frames) ? sceneSpec.frames[0] : null;
        return ((frame && frame.traces) || []).filter((trace) => (
          trace
          && trace.key
          && trace.showlegend
          && Array.isArray(trace.points)
          && trace.points.length > 1
        ));
      }

      function ovizPaperVolumeLayerFor(nameFragment) {
        const layers = (sceneSpec.volumes && sceneSpec.volumes.layers) || [];
        const fragment = String(nameFragment || "").toLowerCase();
        return layers.find((layer) => (
          layer && String(layer.name || "").toLowerCase().includes(fragment)
        )) || null;
      }

      function ovizPaperRefreshVolume(layer) {
        if (!layer || typeof volumeRuntimeByKey === "undefined") {
          return;
        }
        const runtime = volumeRuntimeByKey.get(String(layer.key));
        if (runtime && typeof applyVolumeStateToRuntime === "function") {
          applyVolumeStateToRuntime(layer, runtime);
        }
      }

      function ovizPaperApplyHoverCue(cueEl) {
        if (!cueEl || cueEl === ovizPaperHoverCueEl) {
          return;
        }
        // Never fight an in-flight state transition; the cue re-applies on the
        // next hover once the scene has settled.
        if (typeof ovizStateTransition !== "undefined" && ovizStateTransition) {
          return;
        }
        ovizPaperClearHoverCue();
        ovizPaperHoverCueEl = cueEl;
        const targetNames = String(cueEl.dataset.cueTraces || "")
          .split("|").map((name) => name.trim()).filter(Boolean);
        const dimOpacity = Math.min(
          Math.max(Number(cueEl.dataset.cueDim ?? 0.12), 0.0), 1.0
        );
        if (targetNames.length && typeof actionHeldTraceOpacityByKey !== "undefined") {
          const wanted = new Set(targetNames);
          const overrides = {};
          let matched = false;
          ovizPaperSampleTraces().forEach((trace) => {
            const isTarget = wanted.has(String(trace.name));
            if (!isTarget) {
              // Never resurrect group-hidden traces just to dim them; only
              // targets may be summoned into view by a hover.
              const visibleNow = typeof traceVisibleForGroupState !== "function"
                || traceVisibleForGroupState(trace, currentGroup);
              if (!visibleNow) {
                return;
              }
            }
            overrides[String(trace.key)] = isTarget ? 1.0 : dimOpacity;
            matched = matched || isTarget;
          });
          if (matched) {
            ovizPaperHoverHeldOpacity = actionHeldTraceOpacityByKey;
            actionHeldTraceOpacityByKey = overrides;
            ovizPaperHoverTracesActive = true;
          }
        }
        const volumeFragment = String(cueEl.dataset.cueVolume || "").trim();
        if (volumeFragment && typeof volumeStateByKey !== "undefined") {
          const layer = ovizPaperVolumeLayerFor(volumeFragment);
          const stateKey = layer ? String(layer.state_key || layer.key || "") : "";
          const state = stateKey ? volumeStateByKey[stateKey] : null;
          if (state) {
            ovizPaperHoverVolumeRestore = {
              layer,
              stateKey,
              alphaCoef: state.alphaCoef,
              galacticLightIntensity: state.galacticLightIntensity,
              galacticAmbient: state.galacticAmbient,
            };
            state.alphaCoef = Math.min(Number(state.alphaCoef || 105.0) * 1.5, 400.0);
            state.galacticLightIntensity = Math.min(
              Number(state.galacticLightIntensity || 1.35) * 1.4, 4.0
            );
            state.galacticAmbient = Math.min(
              Number(state.galacticAmbient || 0.22) + 0.10, 1.0
            );
            ovizPaperRefreshVolume(layer);
          }
        }
        if (ovizPaperHoverTracesActive || ovizPaperHoverVolumeRestore) {
          ovizPaperSetDataset("paperHoverCue", targetNames.join("|") || volumeFragment);
          if (typeof renderDisplayedFrameForViewHandoff === "function") {
            renderDisplayedFrameForViewHandoff();
          }
          if (typeof ovizInvalidateRender === "function") {
            ovizInvalidateRender();
          }
        } else {
          ovizPaperHoverCueEl = null;
        }
      }

      function ovizPaperClearHoverCue() {
        if (!ovizPaperHoverCueEl) {
          return;
        }
        ovizPaperHoverCueEl = null;
        if (ovizPaperHoverTracesActive && typeof actionHeldTraceOpacityByKey !== "undefined") {
          actionHeldTraceOpacityByKey = ovizPaperHoverHeldOpacity;
          ovizPaperHoverHeldOpacity = null;
          ovizPaperHoverTracesActive = false;
        }
        if (ovizPaperHoverVolumeRestore && typeof volumeStateByKey !== "undefined") {
          const restore = ovizPaperHoverVolumeRestore;
          ovizPaperHoverVolumeRestore = null;
          const state = volumeStateByKey[restore.stateKey];
          if (state) {
            state.alphaCoef = restore.alphaCoef;
            state.galacticLightIntensity = restore.galacticLightIntensity;
            state.galacticAmbient = restore.galacticAmbient;
            ovizPaperRefreshVolume(restore.layer);
          }
        }
        ovizPaperSetDataset("paperHoverCue", "");
        if (typeof renderDisplayedFrameForViewHandoff === "function") {
          renderDisplayedFrameForViewHandoff();
        }
        if (typeof ovizInvalidateRender === "function") {
          ovizInvalidateRender();
        }
      }

      function ovizPaperWireHoverCues(scrollEl) {
        scrollEl.addEventListener("mouseover", (event) => {
          const cue = event.target && event.target.closest
            ? event.target.closest(".oviz-paper-cue")
            : null;
          if (cue && scrollEl.contains(cue)) {
            ovizPaperApplyHoverCue(cue);
          }
        });
        scrollEl.addEventListener("mouseout", (event) => {
          if (!ovizPaperHoverCueEl) {
            return;
          }
          const next = event.relatedTarget && event.relatedTarget.closest
            ? event.relatedTarget.closest(".oviz-paper-cue")
            : null;
          if (next !== ovizPaperHoverCueEl) {
            ovizPaperClearHoverCue();
          }
        });
        scrollEl.addEventListener("scroll", () => {
          if (ovizPaperHoverCueEl) {
            ovizPaperClearHoverCue();
          }
        }, { passive: true });
      }

      function ovizPaperSyncLabel() {
        if (ovizPaperSyncMode === "on") return "Sync: on";
        if (ovizPaperSyncMode === "paused") return "Sync: paused — resume";
        return "Sync: off";
      }

      function ovizPaperRenderSyncButton() {
        if (!ovizPaperSyncButtonEl) {
          return;
        }
        ovizPaperSyncButtonEl.textContent = ovizPaperSyncLabel();
        ovizPaperSyncButtonEl.dataset.mode = ovizPaperSyncMode;
        ovizPaperSetDataset("paperSync", ovizPaperSyncMode);
      }

      function ovizPaperSetSyncMode(mode, options = {}) {
        const nextMode = mode === "paused" || mode === "off" ? mode : "on";
        if (nextMode === ovizPaperSyncMode) {
          ovizPaperRenderSyncButton();
          return;
        }
        ovizPaperSyncMode = nextMode;
        if (nextMode === "paused") {
          ovizPaperPausedAnchorId = ovizPaperActiveAnchorId;
        }
        ovizPaperRenderSyncButton();
        if (nextMode === "on" && options.applyActive !== false) {
          const entry = ovizPaperAnchorEntries.find(
            (candidate) => candidate.anchorId === ovizPaperActiveAnchorId
          );
          if (entry && entry.stateId) {
            ovizPaperQueueStateJump(entry.stateId);
          }
        }
      }

      function ovizPaperPauseFromSceneInteraction() {
        if (!ovizPaperProject || !ovizPaperIsVisiblyOpen()) {
          return;
        }
        if (ovizPaperSyncMode === "on") {
          ovizPaperSetSyncMode("paused");
        }
      }

      function ovizPaperIsOpen() {
        return Boolean(ovizPaperPanelEl && ovizPaperPanelEl.dataset.open === "true");
      }

      function ovizPaperIsSuppressedByViewerMode() {
        return Boolean(
          root
          && root.dataset
          && (
            root.dataset.zen === "true"
            || root.dataset.presentationMode === "true"
            || root.dataset.mobile === "true"
          )
        );
      }

      function ovizPaperIsVisiblyOpen() {
        return Boolean(
          ovizPaperIsOpen()
          && !ovizPaperIsCollapsed()
          && !ovizPaperIsSuppressedByViewerMode()
        );
      }

      function ovizPaperSyncPanelAccessibility() {
        if (!ovizPaperPanelEl) {
          return;
        }
        const open = ovizPaperIsOpen();
        const collapsed = ovizPaperIsCollapsed();
        const suppressed = ovizPaperIsSuppressedByViewerMode();
        const interactive = ovizPaperIsVisiblyOpen();
        ovizPaperSetDataset("paperVisible", interactive ? "true" : "false");
        ovizPaperPanelEl.setAttribute("aria-hidden", interactive ? "false" : "true");
        if (interactive) {
          ovizPaperPanelEl.removeAttribute("inert");
        } else {
          ovizPaperPanelEl.setAttribute("inert", "");
        }
        if (ovizPaperToggleButtonEl) {
          ovizPaperToggleButtonEl.setAttribute("aria-pressed", open ? "true" : "false");
          ovizPaperToggleButtonEl.setAttribute("aria-expanded", interactive ? "true" : "false");
          ovizPaperToggleButtonEl.textContent = interactive ? "Paper ▾" : "Paper ▸";
          ovizPaperToggleButtonEl.setAttribute(
            "title",
            open ? "Close the paper reader" : "Open the paper reader"
          );
        }
        if (ovizPaperCollapseButtonEl) {
          ovizPaperCollapseButtonEl.setAttribute("aria-expanded", interactive ? "true" : "false");
        }
        if (ovizPaperExpandTabEl) {
          ovizPaperExpandTabEl.setAttribute("aria-expanded", interactive ? "true" : "false");
        }
        const activeElement = document.activeElement;
        const focusIsInPaperUi = Boolean(
          ovizPaperPanelEl.contains(activeElement)
          || activeElement === ovizPaperExpandTabEl
          || activeElement === ovizPaperToggleButtonEl
        );
        if (!interactive && focusIsInPaperUi) {
          const fallback = suppressed
            ? canvas
            : (open && collapsed ? ovizPaperExpandTabEl : ovizPaperToggleButtonEl);
          if (fallback && typeof fallback.focus === "function") {
            window.requestAnimationFrame(() => fallback.focus());
          }
        }
      }

      function ovizPaperSetOpen(open) {
        if (!ovizPaperPanelEl || !ovizPaperProject) {
          return;
        }
        const next = Boolean(open);
        ovizPaperPanelEl.dataset.open = next ? "true" : "false";
        ovizPaperSetDataset("paperOpen", next ? "true" : "false");
        if (next) {
          ovizPaperScheduleScrollSync();
          ovizPaperScheduleRailLayout();
          ovizPaperRequestKatex();
        }
        ovizPaperSyncExpandTab();
        ovizPaperSyncPanelAccessibility();
        ovizPaperTweenViewOffsetToTarget();
        if (typeof ovizInvalidateRender === "function") {
          ovizInvalidateRender();
        }
      }

      function ovizPaperIsCollapsed() {
        return Boolean(ovizPaperPanelEl && ovizPaperPanelEl.dataset.collapsed === "true");
      }

      function ovizPaperSyncExpandTab() {
        const showTab = Boolean(
          ovizPaperProject
          && ovizPaperIsOpen()
          && ovizPaperIsCollapsed()
          && !ovizPaperIsSuppressedByViewerMode()
        );
        ovizPaperSetDataset("paperCollapsed", ovizPaperIsCollapsed() ? "true" : "false");
        if (ovizPaperExpandTabEl) {
          ovizPaperExpandTabEl.hidden = !showTab;
        }
      }

      function ovizPaperSetCollapsed(collapsed) {
        if (!ovizPaperPanelEl) {
          return;
        }
        ovizPaperPanelEl.dataset.collapsed = collapsed ? "true" : "false";
        ovizPaperSyncExpandTab();
        ovizPaperSyncPanelAccessibility();
        ovizPaperTweenViewOffsetToTarget();
        if (!collapsed) {
          ovizPaperScheduleScrollSync();
        }
      }

      function ovizPaperApplyWidthFraction(fraction, options = {}) {
        const clamped = Math.min(Math.max(Number(fraction) || 0.42, 0.26), 0.62);
        ovizPaperWidthFraction = clamped;
        if (ovizPaperPanelEl) {
          ovizPaperPanelEl.style.width = `${(clamped * 100).toFixed(2)}%`;
        }
        ovizPaperScheduleRailLayout();
        if (root && root.style) {
          root.style.setProperty(
            "--oviz-paper-panel-width",
            `${(clamped * 100).toFixed(2)}%`,
          );
        }
        if (options.tweenOffset !== false) {
          ovizPaperTweenViewOffsetToTarget(options.instantOffset === true);
        }
      }

      function ovizPaperTargetViewOffset() {
        // Content recentres into the region left of the panel: the offset is
        // half the covered fraction, and zero when the panel is out of view.
        const covering = ovizPaperIsVisiblyOpen();
        return { x: covering ? ovizPaperWidthFraction / 2.0 : 0.0, y: 0.0 };
      }

      function ovizPaperTweenViewOffsetToTarget(instant = false) {
        if (typeof applyPaperCameraViewOffset !== "function") {
          return;
        }
        if (ovizPaperOffsetTweenFrame) {
          window.cancelAnimationFrame(ovizPaperOffsetTweenFrame);
          ovizPaperOffsetTweenFrame = 0;
        }
        const target = ovizPaperTargetViewOffset();
        const current = (
          typeof currentPaperCameraViewOffset !== "undefined"
          && currentPaperCameraViewOffset
        ) ? currentPaperCameraViewOffset : { x: 0.0, y: 0.0 };
        const start = { x: Number(current.x) || 0.0, y: Number(current.y) || 0.0 };
        if (
          instant
          || (Math.abs(start.x - target.x) <= 1e-4 && Math.abs(start.y - target.y) <= 1e-4)
        ) {
          applyPaperCameraViewOffset(target);
          if (typeof ovizInvalidateRender === "function") {
            ovizInvalidateRender();
          }
          return;
        }
        const durationMs = 240.0;
        const startMs = performance.now();
        const step = (nowMs) => {
          ovizPaperOffsetTweenFrame = 0;
          const linear = Math.min(Math.max((nowMs - startMs) / durationMs, 0.0), 1.0);
          const eased = linear * linear * (3.0 - (2.0 * linear));
          applyPaperCameraViewOffset({
            x: start.x + ((target.x - start.x) * eased),
            y: start.y + ((target.y - start.y) * eased),
          });
          if (typeof ovizInvalidateRender === "function") {
            ovizInvalidateRender();
          }
          if (linear < 1.0) {
            ovizPaperOffsetTweenFrame = window.requestAnimationFrame(step);
          }
        };
        ovizPaperOffsetTweenFrame = window.requestAnimationFrame(step);
      }

      function ovizPaperWireResizer() {
        if (!ovizPaperResizerEl || !ovizPaperPanelEl) {
          return;
        }
        let dragPointerId = null;
        const onMove = (event) => {
          if (dragPointerId === null || event.pointerId !== dragPointerId) {
            return;
          }
          const rootRect = root.getBoundingClientRect();
          const panelRect = ovizPaperPanelEl.getBoundingClientRect();
          const fraction = (panelRect.right - event.clientX) / Math.max(rootRect.width, 1.0);
          ovizPaperApplyWidthFraction(fraction, { instantOffset: true });
        };
        const endDrag = (event) => {
          if (dragPointerId === null || (event && event.pointerId !== dragPointerId)) {
            return;
          }
          dragPointerId = null;
          root.dataset.paperResizing = "false";
          ovizPaperRenderRail();
          ovizPaperScheduleScrollSync();
        };
        ovizPaperResizerEl.addEventListener("pointerdown", (event) => {
          dragPointerId = event.pointerId;
          root.dataset.paperResizing = "true";
          ovizPaperResizerEl.setPointerCapture(event.pointerId);
          event.preventDefault();
          event.stopPropagation();
        });
        ovizPaperResizerEl.addEventListener("pointermove", onMove);
        ovizPaperResizerEl.addEventListener("pointerup", endDrag);
        ovizPaperResizerEl.addEventListener("pointercancel", endDrag);
      }

      function ovizPaperEntryForStateId(stateId) {
        const wanted = String(stateId === undefined || stateId === null ? "" : stateId);
        if (!wanted) {
          return null;
        }
        return ovizPaperAnchorEntries.find(
          (entry) => entry.stateId !== null && String(entry.stateId) === wanted
        ) || null;
      }

      function ovizPaperHandleTransitionLifecycle(event) {
        // Some state restoration paths rewrite the camera projection after
        // the final transition frame.  Reapply only the Paper overlay; the
        // logical State/Action view offset remains untouched.
        ovizPaperTweenViewOffsetToTarget(true);
        const eventName = String(event && event.type || "");
        const detail = event && event.detail && typeof event.detail === "object"
          ? event.detail
          : {};
        if (eventName === "transition-start") {
          ovizPaperSetStatus("syncing", "Updating scene…");
          return;
        }
        if (eventName === "transition-error") {
          ovizPaperSetStatus("error", "Scene update failed");
          return;
        }
        if (eventName === "transition-cancel") {
          if (root.dataset.paperJumpPending !== "true") {
            ovizPaperSetStatus("ready", "Scene update retargeted");
          }
          return;
        }
        if (eventName !== "transition-complete") {
          return;
        }
        if (root.dataset.paperJumpPending !== "true") {
          ovizPaperSetStatus("ready", "Scene updated");
        }
        const entry = ovizPaperEntryForStateId(detail.id);
        if (
          entry
          && ovizPaperIsVisiblyOpen()
          && entry.anchorId !== ovizPaperActiveAnchorId
        ) {
          // State navigation can originate from States, Actions, keyboard, or
          // the public API.  Follow it in the document without feeding the
          // same State back into the controller.
          ovizPaperScrollToAnchor(entry.anchorId, {
            applyState: false,
            behavior: "smooth",
          });
        }
      }

      function ovizPaperWatchTransitionsForOffset() {
        if (!root) {
          return;
        }
        root.addEventListener("transition-start", ovizPaperHandleTransitionLifecycle);
        root.addEventListener("transition-complete", ovizPaperHandleTransitionLifecycle);
        root.addEventListener("transition-cancel", ovizPaperHandleTransitionLifecycle);
        root.addEventListener("transition-error", ovizPaperHandleTransitionLifecycle);
      }

      function ovizPaperSyncViewerModeVisibility() {
        ovizPaperSyncExpandTab();
        ovizPaperSyncPanelAccessibility();
        ovizPaperTweenViewOffsetToTarget();
        if (ovizPaperIsVisiblyOpen()) {
          ovizPaperScheduleScrollSync();
        }
      }

      function ovizPaperWatchViewerModes() {
        if (!root || typeof MutationObserver !== "function") {
          return;
        }
        if (ovizPaperModeObserver) {
          ovizPaperModeObserver.disconnect();
        }
        ovizPaperModeObserver = new MutationObserver(() => {
          ovizPaperSyncViewerModeVisibility();
        });
        ovizPaperModeObserver.observe(root, {
          attributes: true,
          attributeFilter: ["data-zen", "data-presentation-mode", "data-mobile"],
        });
      }

      function ovizPaperWatchHeaderGeometry() {
        if (!ovizPaperPanelEl || !root || !root.style) {
          return;
        }
        const headEl = ovizPaperPanelEl.querySelector(".oviz-three-paper-head");
        if (!headEl) {
          return;
        }
        const update = () => {
          const height = Math.max(Math.ceil(headEl.getBoundingClientRect().height), 0);
          root.style.setProperty("--oviz-paper-head-height", `${height}px`);
        };
        update();
        if (typeof ResizeObserver === "function") {
          if (ovizPaperHeadResizeObserver) {
            ovizPaperHeadResizeObserver.disconnect();
          }
          ovizPaperHeadResizeObserver = new ResizeObserver(update);
          ovizPaperHeadResizeObserver.observe(headEl);
        }
      }

      function ovizPaperActivateAnchor(entry, options = {}) {
        if (!entry) {
          return;
        }
        const changed = entry.anchorId !== ovizPaperActiveAnchorId;
        ovizPaperActiveAnchorId = entry.anchorId;
        ovizPaperSetDataset("paperActiveAnchor", entry.anchorId);
        ovizPaperAnchorEntries.forEach((candidate) => {
          const isActive = candidate.anchorId === entry.anchorId;
          if (candidate.blockEl) {
            candidate.blockEl.dataset.active = isActive ? "true" : "false";
          }
          if (candidate.railEl) {
            candidate.railEl.dataset.active = isActive ? "true" : "false";
          }
        });
        if (!changed && options.force !== true) {
          return;
        }
        if (options.applyState === false) {
          return;
        }
        if (ovizPaperSyncMode === "paused") {
          if (
            ovizPaperProject.sync.resume_policy === "next-anchor"
            && changed
            && entry.anchorId !== ovizPaperPausedAnchorId
          ) {
            ovizPaperSetSyncMode("on", { applyActive: false });
          } else {
            return;
          }
        }
        if (ovizPaperSyncMode !== "on") {
          return;
        }
        if (entry.stateId) {
          ovizPaperQueueStateJump(entry.stateId);
        }
      }

      function ovizPaperRunScrollSync() {
        ovizPaperScrollFrame = 0;
        if (
          !ovizPaperScrollEl
          || !ovizPaperAnchorEntries.length
          || !ovizPaperIsVisiblyOpen()
        ) {
          return;
        }
        if (ovizPaperScrollSyncSuppressed) {
          return;
        }
        if (ovizPaperProject.sync.mode !== "scroll") {
          return;
        }
        const panelRect = ovizPaperScrollEl.getBoundingClientRect();
        const readingLineY = panelRect.top
          + panelRect.height * Number(ovizPaperProject.sync.reading_line);
        let active = null;
        for (const entry of ovizPaperAnchorEntries) {
          if (!entry.blockEl) {
            continue;
          }
          const rect = entry.blockEl.getBoundingClientRect();
          if (rect.top <= readingLineY) {
            active = entry;
          } else {
            break;
          }
        }
        // The closing anchor sits at the very end of the document, where no
        // amount of scrolling can carry it above the reading line — treat
        // reaching the bottom as activating the last anchor.
        const atBottom = (
          ovizPaperScrollEl.scrollTop + ovizPaperScrollEl.clientHeight
          >= ovizPaperScrollEl.scrollHeight - 4.0
        );
        if (atBottom) {
          active = ovizPaperAnchorEntries[ovizPaperAnchorEntries.length - 1];
        }
        if (!active) {
          active = ovizPaperAnchorEntries[0];
        }
        ovizPaperActivateAnchor(active);
      }

      function ovizPaperScheduleScrollSync() {
        if (ovizPaperScrollFrame) {
          return;
        }
        ovizPaperScrollFrame = window.requestAnimationFrame(() => ovizPaperRunScrollSync());
      }

      function ovizPaperReleaseScrollSyncAfterSettle(delayMs = 180) {
        if (ovizPaperScrollSettleTimer) {
          window.clearTimeout(ovizPaperScrollSettleTimer);
        }
        ovizPaperScrollSettleTimer = window.setTimeout(() => {
          ovizPaperScrollSettleTimer = 0;
          ovizPaperScrollSyncSuppressed = false;
          ovizPaperScheduleScrollSync();
        }, Math.max(Number(delayMs) || 0, 0));
      }

      function ovizPaperHandleScroll() {
        if (ovizPaperScrollSyncSuppressed) {
          // Smooth programmatic scrolling emits many intermediate events.
          // Hold scene sync until the destination has actually settled.
          ovizPaperReleaseScrollSyncAfterSettle();
        }
        ovizPaperScheduleScrollSync();
      }

      function ovizPaperRequestKatex() {
        if (
          ovizPaperKatexRequested
          || !ovizPaperProject
          || !ovizPaperProject.math
          || ovizPaperProject.math.enabled === false
          || !ovizPaperBodyEl
        ) {
          return;
        }
        const hasMathDelimiters = /\\$[^$]/.test(ovizPaperBodyEl.textContent || "");
        if (!hasMathDelimiters) {
          ovizPaperKatexRequested = true;
          return;
        }
        ovizPaperKatexRequested = true;
        const version = String(ovizPaperProject.math.katex_version || "0.16.11");
        const base = `https://cdn.jsdelivr.net/npm/katex@${version}/dist/`;
        const link = document.createElement("link");
        link.rel = "stylesheet";
        link.href = `${base}katex.min.css`;
        link.dataset.ovizPaperKatex = "true";
        document.head.appendChild(link);
        const script = document.createElement("script");
        script.src = `${base}katex.min.js`;
        script.dataset.ovizPaperKatex = "true";
        script.onload = () => {
          const auto = document.createElement("script");
          auto.src = `${base}contrib/auto-render.min.js`;
          auto.dataset.ovizPaperKatex = "true";
          auto.onload = () => {
            try {
              if (typeof window.renderMathInElement === "function") {
                window.renderMathInElement(ovizPaperBodyEl, {
                  delimiters: [
                    { left: "$$", right: "$$", display: true },
                    { left: "$", right: "$", display: false },
                  ],
                  throwOnError: false,
                });
              }
            } catch (_katexError) {
              // Fail soft: raw TeX stays readable.
            }
          };
          document.head.appendChild(auto);
        };
        script.onerror = () => {
          // CDN unreachable: leave the raw TeX in place.
        };
        document.head.appendChild(script);
      }

      function ovizPaperLightboxFocusableElements() {
        if (!ovizPaperLightboxEl) {
          return [];
        }
        return Array.from(ovizPaperLightboxEl.querySelectorAll(
          'a[href], button:not([disabled]), input:not([disabled]), select:not([disabled]), '
          + 'textarea:not([disabled]), [tabindex]:not([tabindex="-1"])'
        )).filter((element) => (
          !element.hasAttribute("inert")
          && element.getAttribute("aria-hidden") !== "true"
        ));
      }

      function ovizPaperHandleLightboxKeydown(event) {
        if (!ovizPaperLightboxEl || ovizPaperLightboxEl.dataset.open !== "true") {
          return;
        }
        if (event.key === "Escape") {
          event.preventDefault();
          event.stopPropagation();
          ovizPaperCloseLightbox();
          return;
        }
        if (event.key !== "Tab") {
          return;
        }
        const focusable = ovizPaperLightboxFocusableElements();
        if (!focusable.length) {
          event.preventDefault();
          ovizPaperLightboxEl.focus();
          return;
        }
        const first = focusable[0];
        const last = focusable[focusable.length - 1];
        const active = document.activeElement;
        if (event.shiftKey && (active === first || !ovizPaperLightboxEl.contains(active))) {
          event.preventDefault();
          last.focus();
        } else if (!event.shiftKey && (
          active === last || !ovizPaperLightboxEl.contains(active)
        )) {
          event.preventDefault();
          first.focus();
        }
      }

      function ovizPaperOpenLightbox(dataUrl, captionHtml) {
        if (!ovizPaperLightboxEl) {
          return;
        }
        ovizPaperLightboxRestoreFocusEl = (
          document.activeElement && typeof document.activeElement.focus === "function"
        ) ? document.activeElement : null;
        const image = ovizPaperLightboxEl.querySelector("img");
        const caption = ovizPaperLightboxEl.querySelector(".oviz-three-paper-lightbox-caption");
        if (image) {
          image.src = dataUrl;
        }
        if (caption) {
          caption.innerHTML = captionHtml || "";
        }
        ovizPaperLightboxEl.dataset.open = "true";
        ovizPaperLightboxEl.setAttribute("aria-hidden", "false");
        ovizPaperLightboxEl.removeAttribute("inert");
        const closeButton = ovizPaperLightboxEl.querySelector(
          ".oviz-three-paper-lightbox-close"
        );
        window.requestAnimationFrame(() => {
          if (closeButton && typeof closeButton.focus === "function") {
            closeButton.focus();
          } else if (typeof ovizPaperLightboxEl.focus === "function") {
            ovizPaperLightboxEl.focus();
          }
        });
      }

      function ovizPaperCloseLightbox() {
        if (ovizPaperLightboxEl) {
          ovizPaperLightboxEl.dataset.open = "false";
          ovizPaperLightboxEl.setAttribute("aria-hidden", "true");
          ovizPaperLightboxEl.setAttribute("inert", "");
          const restore = ovizPaperLightboxRestoreFocusEl;
          ovizPaperLightboxRestoreFocusEl = null;
          if (restore && typeof restore.focus === "function") {
            window.requestAnimationFrame(() => restore.focus());
          }
        }
      }

      function ovizPaperRenderBlock(block, sectionEl) {
        const blockEl = document.createElement("div");
        blockEl.className = "oviz-three-paper-block";
        blockEl.dataset.blockId = block.id;
        blockEl.dataset.blockType = block.type;
        if (block.type === "figure") {
          const figure = document.createElement("figure");
          const dataUrl = ovizPaperAssetValue(block.image_data_url || block.asset);
          if (dataUrl) {
            const image = document.createElement("img");
            image.src = dataUrl;
            image.loading = "lazy";
            image.alt = "Paper figure";
            image.tabIndex = 0;
            image.setAttribute("role", "button");
            image.setAttribute("aria-label", "Enlarge paper figure");
            image.addEventListener("click", () => {
              ovizPaperOpenLightbox(dataUrl, block.caption_html || "");
            });
            image.addEventListener("keydown", (event) => {
              if (event.key === "Enter" || event.key === " ") {
                event.preventDefault();
                ovizPaperOpenLightbox(dataUrl, block.caption_html || "");
              }
            });
            // Decoded figures change every offset below them: refresh the
            // progress rail and the active-anchor computation.
            image.addEventListener("load", () => {
              ovizPaperRenderRail();
              ovizPaperScheduleScrollSync();
            });
            figure.appendChild(image);
          }
          if (block.live) {
            const liveBadge = document.createElement("div");
            liveBadge.className = "oviz-three-paper-live-badge";
            liveBadge.textContent = "Shown live in the scene behind";
            figure.appendChild(liveBadge);
          }
          if (block.caption_html) {
            const caption = document.createElement("figcaption");
            caption.innerHTML = block.caption_html;
            figure.appendChild(caption);
          }
          blockEl.appendChild(figure);
        } else {
          blockEl.innerHTML = block.html || "";
        }
        if (block.anchor) {
          blockEl.dataset.anchorId = block.anchor.id;
          blockEl.classList.add("oviz-three-paper-anchored");
          const marker = document.createElement("button");
          marker.type = "button";
          marker.className = "oviz-three-paper-anchor-marker";
          marker.title = block.anchor.label || "Show this view in the scene";
          marker.setAttribute("aria-label", marker.title);
          marker.addEventListener("click", () => {
            const entry = ovizPaperAnchorEntries.find(
              (candidate) => candidate.anchorId === block.anchor.id
            );
            if (entry) {
              ovizPaperSetSyncMode("on", { applyActive: false });
              ovizPaperActivateAnchor(entry, { force: true });
            }
          });
          blockEl.appendChild(marker);
        }
        sectionEl.appendChild(blockEl);
        return blockEl;
      }

      function ovizPaperRenderRail() {
        if (!ovizPaperRailEl) {
          return;
        }
        ovizPaperRailEl.innerHTML = "";
        const scrollHeight = Math.max(ovizPaperScrollEl.scrollHeight, 1);
        ovizPaperAnchorEntries.forEach((entry) => {
          if (!entry.blockEl) {
            return;
          }
          const tick = document.createElement("button");
          tick.type = "button";
          tick.className = "oviz-three-paper-rail-tick";
          tick.title = entry.label || "Anchor";
          tick.setAttribute("aria-label", entry.label || "Show paper view");
          tick.style.top = `${(100.0 * entry.blockEl.offsetTop / scrollHeight).toFixed(3)}%`;
          tick.addEventListener("click", () => {
            ovizPaperSetSyncMode("on", { applyActive: false });
            ovizPaperScrollToAnchor(entry.anchorId);
          });
          entry.railEl = tick;
          ovizPaperRailEl.appendChild(tick);
        });
      }

      function ovizPaperScheduleRailLayout() {
        if (ovizPaperRailLayoutFrame) {
          return;
        }
        ovizPaperRailLayoutFrame = window.requestAnimationFrame(() => {
          ovizPaperRailLayoutFrame = 0;
          if (!ovizPaperRailEl || !ovizPaperScrollEl || !ovizPaperBodyEl) {
            return;
          }
          ovizPaperRenderRail();
        });
      }

      function ovizPaperWatchBodyGeometry() {
        if (!ovizPaperBodyEl || typeof ResizeObserver !== "function") {
          return;
        }
        if (ovizPaperBodyResizeObserver) {
          ovizPaperBodyResizeObserver.disconnect();
        }
        ovizPaperBodyResizeObserver = new ResizeObserver(() => {
          ovizPaperScheduleRailLayout();
        });
        ovizPaperBodyResizeObserver.observe(ovizPaperBodyEl);
      }

      function ovizPaperScrollToAnchor(anchorId, options = {}) {
        const entry = ovizPaperAnchorEntries.find((candidate) => candidate.anchorId === anchorId);
        if (!entry || !entry.blockEl || !ovizPaperScrollEl) {
          return;
        }
        const readingLine = Number(ovizPaperProject.sync.reading_line);
        const target = entry.blockEl.offsetTop
          - ovizPaperScrollEl.clientHeight * readingLine + 8.0;
        ovizPaperScrollSyncSuppressed = true;
        ovizPaperReleaseScrollSyncAfterSettle(options.behavior === "auto" ? 40 : 180);
        ovizPaperScrollEl.scrollTo({
          top: Math.max(target, 0),
          behavior: options.behavior === "auto" ? "auto" : "smooth",
        });
        ovizPaperActivateAnchor(entry, {
          force: true,
          applyState: options.applyState !== false,
        });
      }

      function ovizPaperRenderBody() {
        ovizPaperBodyEl.innerHTML = "";
        ovizPaperAnchorEntries = [];
        (ovizPaperProject.sections || []).forEach((section) => {
          const sectionEl = document.createElement("section");
          sectionEl.className = "oviz-three-paper-section";
          sectionEl.dataset.sectionId = section.id;
          if (section.title_html) {
            const heading = document.createElement(
              section.level <= 1 ? "h2" : (section.level === 2 ? "h3" : "h4")
            );
            heading.innerHTML = section.title_html;
            sectionEl.appendChild(heading);
          }
          (section.blocks || []).forEach((block) => {
            const blockEl = ovizPaperRenderBlock(block, sectionEl);
            if (block.anchor) {
              ovizPaperAnchorEntries.push({
                anchorId: block.anchor.id,
                stateId: block.anchor.state_id || null,
                label: block.anchor.label || "",
                blockEl,
                railEl: null,
              });
            }
          });
          ovizPaperBodyEl.appendChild(sectionEl);
        });
        ovizPaperSetDataset("paperAnchorCount", ovizPaperAnchorEntries.length);
        // The panel is still display:none during initial construction, so its
        // block offsets are zero in this synchronous turn.  Measure on the
        // next frame after ovizPaperSetOpen() has made the reader visible.
        ovizPaperScheduleRailLayout();
      }

      function initializeOvizPaper() {
        const spec = sceneSpec.paper;
        ovizPaperPanelEl = root.querySelector(".oviz-three-paper-panel");
        ovizPaperToggleButtonEl = root.querySelector(".oviz-three-paper-toggle");
        if (!spec || !spec.available || !ovizPaperPanelEl) {
          if (ovizPaperPanelEl) {
            ovizPaperPanelEl.remove();
          }
          if (ovizPaperToggleButtonEl) {
            ovizPaperToggleButtonEl.remove();
          }
          const strayExpandTab = root.querySelector(".oviz-three-paper-expand-tab");
          if (strayExpandTab) {
            strayExpandTab.remove();
          }
          ovizPaperSetDataset("paperReady", "unavailable");
          return;
        }
        ovizPaperProject = spec;
        ovizPaperScrollEl = ovizPaperPanelEl.querySelector(".oviz-three-paper-scroll");
        ovizPaperBodyEl = ovizPaperPanelEl.querySelector(".oviz-three-paper-body");
        ovizPaperRailEl = ovizPaperPanelEl.querySelector(".oviz-three-paper-rail");
        ovizPaperSyncButtonEl = ovizPaperPanelEl.querySelector(".oviz-three-paper-sync");
        ovizPaperStatusEl = ovizPaperPanelEl.querySelector(".oviz-three-paper-status");
        ovizPaperLightboxEl = root.querySelector(".oviz-three-paper-lightbox");
        ovizPaperResizerEl = ovizPaperPanelEl.querySelector(".oviz-three-paper-resizer");
        ovizPaperCollapseButtonEl = ovizPaperPanelEl.querySelector(".oviz-three-paper-collapse");
        ovizPaperExpandTabEl = root.querySelector(".oviz-three-paper-expand-tab");
        ovizPaperPanelEl.style.setProperty(
          "--oviz-paper-panel-width",
          `${(100.0 * Number(spec.panel.width_fraction)).toFixed(2)}%`
        );
        ovizPaperApplyWidthFraction(Number(spec.panel.width_fraction), { tweenOffset: false });

        const titleEl = ovizPaperPanelEl.querySelector(".oviz-three-paper-doc-title");
        if (titleEl) {
          if (spec.link_url) {
            const link = document.createElement("a");
            link.href = spec.link_url;
            link.target = "_blank";
            link.rel = "noopener";
            link.textContent = spec.title || "Paper";
            titleEl.appendChild(link);
          } else {
            titleEl.textContent = spec.title || "Paper";
          }
        }
        const authorsEl = ovizPaperPanelEl.querySelector(".oviz-three-paper-doc-authors");
        if (authorsEl) {
          authorsEl.textContent = (spec.authors || []).join(", ");
        }
        const venueEl = ovizPaperPanelEl.querySelector(".oviz-three-paper-doc-venue");
        if (venueEl) {
          venueEl.innerHTML = spec.venue_html || "";
        }

        ovizPaperRenderBody();
        ovizPaperWatchBodyGeometry();
        ovizPaperRenderSyncButton();
        ovizPaperSetStatus("ready", "Scene ready");
        ovizPaperWatchHeaderGeometry();

        ovizPaperScrollEl.addEventListener("scroll", ovizPaperHandleScroll, { passive: true });
        ovizPaperWireHoverCues(ovizPaperScrollEl);
        ovizPaperPanelEl.addEventListener("pointerdown", (event) => event.stopPropagation());
        const closeButton = ovizPaperPanelEl.querySelector(".oviz-three-paper-close");
        if (closeButton) {
          closeButton.addEventListener("click", () => ovizPaperSetOpen(false));
        }
        if (ovizPaperSyncButtonEl) {
          ovizPaperSyncButtonEl.addEventListener("click", () => {
            ovizPaperSetSyncMode(ovizPaperSyncMode === "on" ? "off" : "on");
          });
        }
        if (ovizPaperToggleButtonEl) {
          ovizPaperToggleButtonEl.addEventListener("click", () => {
            if (ovizPaperIsOpen() && ovizPaperIsCollapsed()) {
              ovizPaperSetCollapsed(false);
              return;
            }
            ovizPaperSetOpen(!ovizPaperIsOpen());
          });
        }
        if (ovizPaperCollapseButtonEl) {
          ovizPaperCollapseButtonEl.addEventListener("click", () => ovizPaperSetCollapsed(true));
        }
        if (ovizPaperExpandTabEl) {
          ovizPaperExpandTabEl.addEventListener("click", () => ovizPaperSetCollapsed(false));
        }
        ovizPaperWireResizer();
        ovizPaperWatchTransitionsForOffset();
        ovizPaperWatchViewerModes();
        ovizPaperSyncExpandTab();
        if (ovizPaperLightboxEl) {
          ovizPaperLightboxEl.addEventListener("click", (event) => {
            if (event.target === ovizPaperLightboxEl) {
              ovizPaperCloseLightbox();
            }
          });
          const lightboxClose = ovizPaperLightboxEl.querySelector(
            ".oviz-three-paper-lightbox-close"
          );
          if (lightboxClose) {
            lightboxClose.addEventListener("click", () => ovizPaperCloseLightbox());
          }
          root.addEventListener("keydown", ovizPaperHandleLightboxKeydown);
        }
        // Reading is scroll-driven; grabbing the scene pauses sync so the
        // camera is never fought over.
        canvas.addEventListener("pointerdown", ovizPaperPauseFromSceneInteraction, { capture: true });
        canvas.addEventListener("wheel", ovizPaperPauseFromSceneInteraction, {
          capture: true,
          passive: true,
        });
        if (ovizPaperProject.sync.mode !== "scroll") {
          ovizPaperSyncMode = "off";
          ovizPaperRenderSyncButton();
        }

        try {
          if (window.Oviz && window.Oviz.__viewers && window.Oviz.__viewers.get) {
            const viewerRecord = window.Oviz.__viewers.get(root.id);
            if (viewerRecord) {
              viewerRecord.paper = {
                open: () => ovizPaperSetOpen(true),
                close: () => ovizPaperSetOpen(false),
                goToAnchor: (anchorId) => ovizPaperScrollToAnchor(String(anchorId)),
                setSync: (mode) => ovizPaperSetSyncMode(mode),
                anchors: () => ovizPaperAnchorEntries.map((entry) => ({
                  id: entry.anchorId,
                  stateId: entry.stateId,
                  label: entry.label,
                })),
              };
            }
          }
        } catch (_apiError) {
        }

        ovizPaperSetDataset("paperReady", "true");
        ovizPaperSetOpen(Boolean(spec.enabled));
        if (spec.enabled) {
          ovizPaperScheduleScrollSync();
        }
      }
"""
