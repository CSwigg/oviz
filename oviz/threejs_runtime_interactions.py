from __future__ import annotations


THREEJS_INTERACTION_RUNTIME_JS = """
      function renderTimeSliderTicks() {
        if (!sliderMinorTicksEl || !sliderMajorTicksEl || !sliderLabelsEl) {
          return;
        }
        sliderMinorTicksEl.innerHTML = "";
        sliderMajorTicksEl.innerHTML = "";
        sliderLabelsEl.innerHTML = "";
      }

      function updatePlaybackButtons() {
        const isReverse = playbackDirection < 0;
        const isForward = playbackDirection > 0;
        if (playBackwardButtonEl) {
          playBackwardButtonEl.dataset.active = isReverse ? "true" : "false";
          playBackwardButtonEl.textContent = isReverse ? "⏸" : "◀";
          playBackwardButtonEl.setAttribute("title", isReverse ? "Pause reverse playback" : "Play backward");
          playBackwardButtonEl.setAttribute("aria-label", isReverse ? "Pause reverse playback" : "Play backward");
        }
        if (playForwardButtonEl) {
          playForwardButtonEl.dataset.active = isForward ? "true" : "false";
          playForwardButtonEl.textContent = isForward ? "⏸" : "▶";
          playForwardButtonEl.setAttribute("title", isForward ? "Pause forward playback" : "Play forward");
          playForwardButtonEl.setAttribute("aria-label", isForward ? "Pause forward playback" : "Play forward");
        }
      }

      function updateTimeSliderTickState(activeFrameIndex = currentFrameIndex) {
        if (!sliderTrackWrapEl) {
          return;
        }
        const activeIndex = String(clampFrameIndex(activeFrameIndex));
        sliderTrackWrapEl.querySelectorAll(".oviz-three-slider-tick, .oviz-three-slider-tick-label").forEach((el) => {
          el.dataset.active = el.dataset.frameIndex === activeIndex ? "true" : "false";
        });
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
        const nextKey = normalizeMemberKey(clusterKey);
        if (nextKey === localHoveredClusterKey) {
          return;
        }
        localHoveredClusterKey = nextKey;
        applySceneHoverState();
        if (widgetModeForKey("dendrogram") !== "hidden") {
          renderDendrogramWidget();
        }
        postParentHoverToSkyFrame();
      }

      function setSkyHoveredClusterKey(clusterKey) {
        const nextKey = crossHoverEnabled() ? normalizeMemberKey(clusterKey) : "";
        if (nextKey === skyHoveredClusterKey) {
          return;
        }
        skyHoveredClusterKey = nextKey;
        applySceneHoverState();
        if (widgetModeForKey("dendrogram") !== "hidden") {
          renderDendrogramWidget();
        }
      }

      function clearCrossHoverState() {
        localHoveredClusterKey = "";
        skyHoveredClusterKey = "";
        lastSentSkyHoverClusterKey = null;
        applySceneHoverState();
        if (widgetModeForKey("dendrogram") !== "hidden") {
          renderDendrogramWidget();
        }
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
        return lassoVolumeSelectionEnabled && lassoSelectionFilterEnabled ? currentLassoSelectionMask : null;
      }

      function disposeLassoSelectionMask(mask) {
        if (!mask || !mask.maskTexture) {
          return;
        }
        mask.maskTexture.dispose();
      }

      function lassoSelectionFilterActive() {
        return Boolean(
          lassoSelectionFilterEnabled
          && currentSelectionMode === "lasso"
          && (selectedClusterKeys.size || hasActiveLassoSelectionMask())
        );
      }

      function selectionStateHasContent() {
        return Boolean(currentSelection || currentSelections.length || selectedClusterKeys.size || hasActiveLassoSelectionMask());
      }

      function selectionMaskInUse(mask, ignoredSnapshot = null) {
        if (!mask) {
          return false;
        }
        if (currentLassoSelectionMask === mask) {
          return true;
        }
        return selectionUndoStack.some((snapshot) => (
          snapshot
          && snapshot !== ignoredSnapshot
          && snapshot.currentLassoSelectionMask === mask
        ));
      }

      function disposeSelectionMaskIfUnused(mask, ignoredSnapshot = null) {
        if (!mask || selectionMaskInUse(mask, ignoredSnapshot)) {
          return;
        }
        disposeLassoSelectionMask(mask);
      }

      function selectionStateSnapshot() {
        return {
          currentSelection: safeJsonClone(currentSelection, null),
          currentSelections: safeJsonClone(currentSelections, []),
          selectedClusterKeys: Array.from(selectedClusterKeys),
          currentSelectionMode,
          lassoSelectionFilterEnabled,
          currentLassoSelectionMask,
        };
      }

      function releaseSelectionSnapshot(snapshot) {
        if (!snapshot || !snapshot.currentLassoSelectionMask) {
          return;
        }
        disposeSelectionMaskIfUnused(snapshot.currentLassoSelectionMask, snapshot);
      }

      function pushSelectionUndoState() {
        if (minimalModeEnabled) {
          return;
        }
        selectionUndoStack.push(selectionStateSnapshot());
        while (selectionUndoStack.length > MAX_SELECTION_UNDO_STATES) {
          releaseSelectionSnapshot(selectionUndoStack.shift());
        }
      }

      function restoreSelectionSnapshot(snapshot) {
        if (!snapshot) {
          return false;
        }
        const previousMask = currentLassoSelectionMask;
        currentSelection = null;
        currentSelections = uniqueSelections(Array.isArray(snapshot.currentSelections) ? snapshot.currentSelections : []);
        selectedClusterKeys = normalizeSelectionKeySet(snapshot.selectedClusterKeys || []);
        currentSelectionMode = currentSelections.length || selectedClusterKeys.size || snapshot.currentLassoSelectionMask
          ? "lasso"
          : "none";
        lassoSelectionFilterEnabled = snapshot.lassoSelectionFilterEnabled !== false;
        currentLassoSelectionMask = snapshot.currentLassoSelectionMask || null;
        disposeSelectionMaskIfUnused(previousMask);
        clearCrossHoverState();
        updateSelectionUI();
        updateSkyPanel();
        renderFrame(currentFrameIndex);
        return true;
      }

      function undoSelectionState() {
        const snapshot = selectionUndoStack.pop();
        if (!snapshot) {
          return false;
        }
        return restoreSelectionSnapshot(snapshot);
      }

      function toggleLassoSelectionFilter() {
        if (!selectionStateHasContent()) {
          return false;
        }
        pushSelectionUndoState();
        lassoSelectionFilterEnabled = !lassoSelectionFilterEnabled;
        currentSelectionMode = currentSelections.length || selectedClusterKeys.size || hasActiveLassoSelectionMask()
          ? "lasso"
          : "none";
        updateSelectionUI();
        updateSkyPanel();
        renderFrame(currentFrameIndex);
        return true;
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

      function setClusterSelections(selections, mode = "lasso", options = {}) {
        if (minimalModeEnabled) {
          currentSelection = null;
          currentSelections = [];
          selectedClusterKeys = new Set();
          currentSelectionMode = "none";
          updateSelectionUI();
          return;
        }
        const nextSelections = uniqueSelections(selections);
        const normalizedMode = String(mode || "lasso");
        if (normalizedMode !== "lasso") {
          return;
        }
        const hasReplacementMask = Boolean(options && Object.prototype.hasOwnProperty.call(options, "lassoMask"));
        const replacementMask = hasReplacementMask ? (options.lassoMask || null) : currentLassoSelectionMask;
        pushSelectionUndoState();

        const previousMask = currentLassoSelectionMask;
        if (hasReplacementMask) {
          currentLassoSelectionMask = replacementMask;
          disposeSelectionMaskIfUnused(previousMask);
        }
        currentSelections = nextSelections;
        currentSelection = null;
        selectedClusterKeys = new Set(
          currentSelections
            .map((selection) => normalizedSelectionKeyFor(selection))
            .filter(Boolean)
        );
        lassoSelectionFilterEnabled = true;

        currentSelectionMode = currentSelections.length || hasActiveLassoSelectionMask() ? "lasso" : "none";
        clearCrossHoverState();
        updateSelectionUI();
        updateSkyPanel();
        renderFrame(currentFrameIndex);
      }

      function clearClusterSelections() {
        if (!selectionStateHasContent()) {
          return;
        }
        pushSelectionUndoState();
        const previousMask = currentLassoSelectionMask;
        currentSelections = [];
        currentSelection = null;
        currentLassoSelectionMask = null;
        currentSelectionMode = "none";
        lassoSelectionFilterEnabled = true;
        selectedClusterKeys = new Set();
        disposeSelectionMaskIfUnused(previousMask);
        clearCrossHoverState();
        updateSelectionUI();
        updateSkyPanel();
        renderFrame(currentFrameIndex);
      }

      function setDrawerAccessibility(shellEl, toggleEl, drawerSelector, isOpen) {
        if (!shellEl || !toggleEl) {
          return;
        }
        const drawerEl = shellEl.querySelector(drawerSelector);
        toggleEl.setAttribute("aria-expanded", isOpen ? "true" : "false");
        if (drawerEl) {
          drawerEl.setAttribute("aria-hidden", isOpen ? "false" : "true");
          if (isOpen) {
            drawerEl.removeAttribute("inert");
          } else {
            drawerEl.setAttribute("inert", "");
          }
        }
      }

      function setToolsDrawerOpen(isOpen) {
        if (!toolsShellEl || !toolsToggleEl) {
          return;
        }
        const nextOpen = minimalModeEnabled ? false : Boolean(isOpen);
        toolsShellEl.dataset.open = nextOpen ? "true" : "false";
        toolsToggleEl.textContent = nextOpen ? "Selection ▾" : "Selection ▸";
        setDrawerAccessibility(toolsShellEl, toolsToggleEl, ".oviz-three-tools-drawer", nextOpen);
        if (nextOpen && controlsShellEl && controlsShellEl.dataset.open === "true") {
          controlsShellEl.dataset.open = "false";
          if (controlsToggleEl) {
            controlsToggleEl.textContent = "Controls ▸";
            setDrawerAccessibility(controlsShellEl, controlsToggleEl, ".oviz-three-controls-drawer", false);
          }
        }
        if (nextOpen && skyControlsShellEl && skyControlsShellEl.dataset.open === "true") {
          skyControlsShellEl.dataset.open = "false";
          if (skyControlsToggleEl) {
            skyControlsToggleEl.textContent = "Sky ▸";
            setDrawerAccessibility(skyControlsShellEl, skyControlsToggleEl, ".oviz-three-sky-controls-drawer", false);
          }
        }
      }

      function setControlsDrawerOpen(isOpen) {
        if (!controlsShellEl || !controlsToggleEl) {
          return;
        }
        const nextOpen = minimalModeEnabled ? false : Boolean(isOpen);
        controlsShellEl.dataset.open = nextOpen ? "true" : "false";
        controlsToggleEl.textContent = mobileModeEnabled ? "Controls" : (nextOpen ? "Controls ▾" : "Controls ▸");
        setDrawerAccessibility(controlsShellEl, controlsToggleEl, ".oviz-three-controls-drawer", nextOpen);
        if (mobileModeEnabled && nextOpen) {
          setLegendPanelOpen(false);
        }
        if (nextOpen && toolsShellEl && toolsShellEl.dataset.open === "true") {
          toolsShellEl.dataset.open = "false";
          if (toolsToggleEl) {
            toolsToggleEl.textContent = "Selection ▸";
            setDrawerAccessibility(toolsShellEl, toolsToggleEl, ".oviz-three-tools-drawer", false);
          }
        }
        if (nextOpen && skyControlsShellEl && skyControlsShellEl.dataset.open === "true") {
          skyControlsShellEl.dataset.open = "false";
          if (skyControlsToggleEl) {
            skyControlsToggleEl.textContent = "Sky ▸";
            setDrawerAccessibility(skyControlsShellEl, skyControlsToggleEl, ".oviz-three-sky-controls-drawer", false);
          }
        }
      }

      function setSkyControlsDrawerOpen(isOpen) {
        if (!skyControlsShellEl || !skyControlsToggleEl) {
          return;
        }
        const nextOpen = minimalModeEnabled ? false : Boolean(isOpen);
        skyControlsShellEl.dataset.open = nextOpen ? "true" : "false";
        skyControlsToggleEl.textContent = mobileModeEnabled ? "Layers" : (nextOpen ? "Sky ▾" : "Sky ▸");
        setDrawerAccessibility(skyControlsShellEl, skyControlsToggleEl, ".oviz-three-sky-controls-drawer", nextOpen);
        if (mobileModeEnabled && nextOpen) {
          setLegendPanelOpen(false);
        }
        if (nextOpen && toolsShellEl && toolsShellEl.dataset.open === "true") {
          toolsShellEl.dataset.open = "false";
          if (toolsToggleEl) {
            toolsToggleEl.textContent = "Selection ▸";
            setDrawerAccessibility(toolsShellEl, toolsToggleEl, ".oviz-three-tools-drawer", false);
          }
        }
        if (nextOpen && controlsShellEl && controlsShellEl.dataset.open === "true") {
          controlsShellEl.dataset.open = "false";
          if (controlsToggleEl) {
            controlsToggleEl.textContent = "Controls ▸";
            setDrawerAccessibility(controlsShellEl, controlsToggleEl, ".oviz-three-controls-drawer", false);
          }
        }
      }

      function formatDisplayedTimeLabel(timeValue) {
        return formatTick(timeValue).replace(/^-/, "−");
      }

      function updateTimelineUi(frameValue, timeValue = frameTimeForValue(frameValue)) {
        const clampedValue = clampFrameValue(frameValue);
        if (sliderEl) {
          sliderEl.value = String(clampedValue);
        }
        if (timeLabelEl) {
          timeLabelEl.textContent = `Time (Myr): ${formatDisplayedTimeLabel(timeValue)}`;
        }
        updateTimeSliderTickState(clampFrameIndex(clampedValue));
      }

      function galacticReferenceMotionVisible() {
        return Boolean(
          galacticReferenceVisible
          && cameraViewMode !== "earth"
          && (
            timelineScrubMotionActive
            || playbackDirection !== 0
            || (typeof timeActionTrack !== "undefined" && timeActionTrack)
            || ovizStateTimelineMotionActive
          )
        );
      }

      function galacticReferenceTimeOpacity() {
        const timeMyr = Number(frameTimeForValue(displayedFrameValue));
        if (!Number.isFinite(timeMyr)) {
          return 0.0;
        }
        return Math.abs(timeMyr) <= 1e-9 ? 0.0 : 1.0;
      }

      function milkyWayTimelineOpacity() {
        const timeMyr = Number(frameTimeForValue(displayedFrameValue));
        if (!Number.isFinite(timeMyr)) {
          return 1.0;
        }
        return Math.abs(timeMyr) <= 1e-9 ? 1.0 : 0.0;
      }

      function updateGalacticReferenceMotionOpacity() {
        const scale = galacticReferenceMotionVisible() ? galacticReferenceTimeOpacity() : 0.0;
        galacticReferenceOpacityGroups.forEach((group) => {
          if (!group) return;
          group.visible = scale > 0.001;
          group.traverse((object) => {
            const material = object && object.material;
            if (!material) return;
            const materials = Array.isArray(material) ? material : [material];
            materials.forEach((item) => {
              const baseOpacity = Number(item && item.userData && item.userData.ovizTimelineBaseOpacity);
              if (!item || !Number.isFinite(baseOpacity)) return;
              item.transparent = true;
              item.opacity = baseOpacity * scale;
            });
          });
        });
      }

      function updateMilkyWayTimelineOpacity() {
        const scale = milkyWayTimelineOpacity();
        cameraResponsiveImagePlaneEntries.forEach((entry) => {
          if (entry && entry.milkyWayImage) {
            entry.timeOpacityScale = scale;
          }
        });
        plotGroup.traverse((object) => {
          if (object && object.userData && object.userData.ovizDecorationKind === "milky_way_model") {
            object.userData.ovizTimeOpacityScale = scale;
          }
        });
        setMilkyWayModelOpacityScale(milkyWayViewOpacityScale);
        updateCameraResponsiveImagePlanes();
      }

      function updateTimelineMotionOpacity() {
        updateGalacticReferenceMotionOpacity();
        updateMilkyWayTimelineOpacity();
      }

      function setTimelineScrubMotionActive(active, options = {}) {
        if (timelineMotionHideTimer) {
          window.clearTimeout(timelineMotionHideTimer);
          timelineMotionHideTimer = 0;
        }
        timelineScrubMotionActive = Boolean(active);
        if (root && root.dataset) {
          root.dataset.timelineMotionActive = timelineScrubMotionActive ? "true" : "false";
        }
        const settleDelayMs = Math.max(Number(options.settleDelayMs) || 0.0, 0.0);
        if (timelineScrubMotionActive && settleDelayMs > 0.0) {
          timelineMotionHideTimer = window.setTimeout(() => {
            timelineMotionHideTimer = 0;
            timelineScrubMotionActive = false;
            if (root && root.dataset) {
              root.dataset.timelineMotionActive = "false";
            }
            if (playbackDirection === 0 && !(typeof timeActionTrack !== "undefined" && timeActionTrack)) {
              updateGalacticReferenceMotionOpacity();
            }
          }, settleDelayMs);
        }
      }

      function pause(options = {}) {
        const wasPlaying = playbackDirection !== 0;
        if (!actionInterruptsMuted()) {
          interruptActionRun("time", { disableOrbit: false });
        }
        if (timeActionTrack) {
          stopTimeActionTrack();
        }
        const snapToFrame = options.snap !== false;
        playbackDirection = 0;
        lastPlaybackAdvanceTimestamp = null;
        if (sliderScrubRenderHandle !== null && !snapToFrame) {
          window.cancelAnimationFrame(sliderScrubRenderHandle);
          sliderScrubRenderHandle = null;
        }
        if (frameTransitionState) {
          const snappedIndex = clampFrameIndex(displayedFrameValue);
          frameTransitionState = null;
          if (snapToFrame) {
            renderFrame(snappedIndex);
          }
        } else if (snapToFrame) {
          const snappedIndex = clampFrameIndex(displayedFrameValue);
          if (Math.abs(displayedFrameValue - snappedIndex) > 1e-6) {
            renderFrame(snappedIndex);
          }
        }
        updatePlaybackButtons();
        if (wasPlaying && !timelineScrubMotionActive) {
          renderFrame(currentFrameIndex);
        }
      }

      function play(direction = 1) {
        if (!actionInterruptsMuted()) {
          interruptActionRun("time", { disableOrbit: false });
        }
        if (timeActionTrack) {
          stopTimeActionTrack();
        }
        if (frameSpecs.length <= 1) {
          return;
        }
        const nextDirection = direction < 0 ? -1 : 1;
        if (playbackDirection === nextDirection) {
          pause();
          return;
        }
        frameTransitionState = null;
        pendingSliderFrameIndex = null;
        if (sliderScrubRenderHandle !== null) {
          window.cancelAnimationFrame(sliderScrubRenderHandle);
          sliderScrubRenderHandle = null;
        }
        playbackDirection = nextDirection;
        const now = window.performance ? window.performance.now() : Date.now();
        lastPlaybackAdvanceTimestamp = now;
        updatePlaybackButtons();
        const nextIndex = (currentFrameIndex + playbackDirection + frameSpecs.length) % frameSpecs.length;
        renderFrame(nextIndex);
      }

      function scheduleSliderScrubRender(index) {
        const clampedValue = clampFrameValue(index);
        pendingSliderFrameIndex = clampedValue;
        displayedFrameValue = clampedValue;
        updateTimelineUi(clampedValue, frameTimeForValue(clampedValue));
        if (sliderScrubRenderHandle !== null) {
          return;
        }
        sliderScrubRenderHandle = window.requestAnimationFrame(() => {
          sliderScrubRenderHandle = null;
          const targetIndex = pendingSliderFrameIndex;
          pendingSliderFrameIndex = null;
          if (targetIndex !== null && targetIndex !== undefined) {
            const stableFrameIndex = clampFrameIndex(targetIndex);
            if (stableFrameIndex !== currentFrameIndex || !galacticReferenceOpacityGroups.length) {
              renderFrame(stableFrameIndex);
            }
            displayedFrameValue = clampFrameValue(targetIndex);
            currentFrameIndex = stableFrameIndex;
            updateTimelineUi(displayedFrameValue, frameTimeForValue(displayedFrameValue));
            updateTimelineMotionOpacity();
          }
        });
      }

      function currentFrameAllowsSelection() {
        const frame = currentFrame();
        return !minimalModeEnabled && Boolean(frame);
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
        const nextLassoSelectionMask = captureLassoSelectionMask(polygon);
        const selected = [];
        hoverTargets.forEach((sprite) => {
          const selection = selectionForSprite(sprite);
          if (!selection) {
            return;
          }
          const screenPoint = spriteScreenPoint(sprite);
          if (screenPoint && pointInPolygon(screenPoint, polygon)) {
            selected.push(selection);
          }
        });
        setClusterSelections(selected, "lasso", { lassoMask: nextLassoSelectionMask });
        if (event) {
          event.preventDefault();
        }
      }

      function isClusterInfoSprite(sprite) {
        return Boolean(
          sprite
          && sprite.userData
          && (
            normalizeMemberKey(sprite.userData.selectionKey || "")
            || normalizedSelectionKeyFor(sprite.userData.selection)
          )
          && sprite.userData.hovertext
        );
      }

      function clusterInfoSpritePickRadiusPx(sprite) {
        if (!sprite || !sprite.userData) {
          return 0.0;
        }
        const rect = canvas.getBoundingClientRect();
        const worldPosition = new THREE.Vector3();
        sprite.getWorldPosition(worldPosition);
        const cameraDirection = new THREE.Vector3();
        camera.getWorldDirection(cameraDirection);
        const depth = worldPosition.clone().sub(camera.position).dot(cameraDirection);
        if (!Number.isFinite(depth) || depth <= 1e-6) {
          return 0.0;
        }
        const viewportHeight = Math.max(rect.height || canvas.clientHeight || 1, 1);
        const worldHeight = 2.0 * depth * Math.tan(THREE.MathUtils.degToRad(camera.fov * 0.5));
        const pickWorldScale = Math.max(
          Number(sprite.userData.pickWorldScale ?? sprite.userData.baseScale ?? sprite.scale.x) || 0.0,
          0.0
        );
        if (!Number.isFinite(worldHeight) || worldHeight <= 1e-9 || pickWorldScale <= 0.0) {
          return 0.0;
        }
        const radiusPx = 0.5 * pickWorldScale * viewportHeight / worldHeight;
        const minRadiusPx = cameraViewMode === "earth" ? 5.0 : 3.5;
        const maxRadiusPx = cameraViewMode === "earth" ? 18.0 : 11.0;
        return clampRange(radiusPx, minRadiusPx, maxRadiusPx);
      }

      function pickClusterInfoSpriteByScreenDistance(event) {
        const clickPoint = canvasPointFromEvent(event);
        let best = null;
        hoverTargets.forEach((sprite) => {
          if (!isClusterInfoSprite(sprite) || sprite.visible === false) {
            return;
          }
          const screenPoint = spriteScreenPoint(sprite);
          if (!screenPoint) {
            return;
          }
          const radius = clusterInfoSpritePickRadiusPx(sprite);
          if (!(radius > 0.0)) {
            return;
          }
          const dx = Number(clickPoint.x) - Number(screenPoint.x);
          const dy = Number(clickPoint.y) - Number(screenPoint.y);
          const distanceSq = dx * dx + dy * dy;
          const radiusSq = radius * radius;
          if (distanceSq > radiusSq) {
            return;
          }
          const score = distanceSq / Math.max(radiusSq, 1e-9);
          if (!best || score < best.score) {
            best = { sprite, score };
          }
        });
        return best ? best.sprite : null;
      }

      function pickSprite(event, options = {}) {
        const clusterHit = pickClusterInfoSpriteByScreenDistance(event);
        if (clusterHit) {
          return clusterHit;
        }
        if (options && options.clusterOnly) {
          return null;
        }
        const rect = canvas.getBoundingClientRect();
        pointer.x = ((event.clientX - rect.left) / rect.width) * 2.0 - 1.0;
        pointer.y = -((event.clientY - rect.top) / rect.height) * 2.0 + 1.0;
        raycaster.setFromCamera(pointer, camera);
        const hits = raycaster.intersectObjects(hoverTargets, false);
        for (const hit of hits) {
          const object = hit && hit.object;
          if (!object || isClusterInfoSprite(object)) {
            continue;
          }
          return object;
        }
        return null;
      }

      function pointerRayFromEvent(event) {
        const rect = canvas.getBoundingClientRect();
        pointer.x = ((event.clientX - rect.left) / rect.width) * 2.0 - 1.0;
        pointer.y = -((event.clientY - rect.top) / rect.height) * 2.0 + 1.0;
        raycaster.setFromCamera(pointer, camera);
        return raycaster.ray;
      }

      function doubleClickTargetFromEvent(event) {
        const spriteHitObject = pickSprite(event);
        if (spriteHitObject) {
          const worldPoint = new THREE.Vector3();
          const hitObject = spriteHitObject;
          hitObject.getWorldPosition(worldPoint);
          return {
            worldPoint,
            selection: hitObject.userData ? (hitObject.userData.selection || null) : null,
            selectionKey: hitObject.userData ? normalizeMemberKey(hitObject.userData.selectionKey || "") : "",
            referenceFrameKey: hitObject.userData ? normalizeMemberKey(hitObject.userData.referenceFrameKey || "") : "",
            manualLabelId: hitObject.userData ? String(hitObject.userData.manualLabelId || "") : "",
          };
        }

        pointerRayFromEvent(event);
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
              referenceFrameKey: "",
              manualLabelId: "",
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
              referenceFrameKey: "",
              manualLabelId: "",
              fallbackPlane: true,
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
        setZoomAnchorToCameraTarget();
        applyCameraViewMode();
        controls.update();
        updateScaleBar();
      }

      function hideClusterInfoTooltip() {
        activeClusterInfoSelectionKey = "";
        activeClusterInfoSprite = null;
        tooltipEl.style.display = "none";
        tooltipEl.innerHTML = "";
        tooltipEl.style.borderColor = "";
        tooltipEl.style.color = "";
        tooltipEl.style.boxShadow = "";
      }

      function clusterInfoSpriteForKey(selectionKey) {
        const key = normalizeMemberKey(selectionKey);
        const candidates = key ? (selectionSpriteEntriesByKey.get(key) || []) : [];
        const candidateEntry = candidates.find((entry) => (
          entry
          && entry.sprite
          && entry.sprite.userData
          && entry.sprite.userData.hovertext
        ));
        return candidateEntry ? candidateEntry.sprite : null;
      }

      function activeClusterInfoSpriteObject() {
        const key = normalizeMemberKey(activeClusterInfoSelectionKey);
        if (!key) {
          return null;
        }
        const currentSprite = clusterInfoSpriteForKey(key);
        if (currentSprite) {
          activeClusterInfoSprite = currentSprite;
          return currentSprite;
        }
        if (
          activeClusterInfoSprite
          && activeClusterInfoSprite.userData
          && normalizeMemberKey(activeClusterInfoSprite.userData.selectionKey || "") === key
        ) {
          return activeClusterInfoSprite;
        }
        return null;
      }

      function updateClusterInfoTooltipPosition() {
        if (!activeClusterInfoSelectionKey || tooltipEl.style.display === "none") {
          return;
        }
        const hitObject = activeClusterInfoSpriteObject();
        if (!hitObject) {
          hideClusterInfoTooltip();
          return;
        }
        const rect = canvas.getBoundingClientRect();
        const anchor = spriteScreenPoint(hitObject);
        if (!anchor) {
          hideClusterInfoTooltip();
          return;
        }
        const tooltipWidth = tooltipEl.offsetWidth || 0;
        const tooltipHeight = tooltipEl.offsetHeight || 0;
        const horizontalOffset = 16;
        let left = anchor.x + horizontalOffset;
        let top = anchor.y - Math.max(tooltipHeight * 0.5, 8);
        if (left + tooltipWidth > rect.width - 4) {
          left = anchor.x - tooltipWidth - horizontalOffset;
        }
        left = clampRange(left, 4, Math.max(rect.width - tooltipWidth - 4, 4));
        top = clampRange(top, 4, Math.max(rect.height - tooltipHeight - 4, 4));
        tooltipEl.style.left = `${left}px`;
        tooltipEl.style.top = `${top}px`;
      }

      function showClusterInfoTooltip(hitObject) {
        if (!hitObject || !hitObject.userData || !hitObject.userData.hovertext) {
          return false;
        }
        const tooltipColor = String(
          hitObject.userData.tooltipColor
          || (hitObject.userData.selection && hitObject.userData.selection.cluster_color)
          || "#ffffff"
        );
        const selectionKey = normalizeMemberKey(hitObject.userData.selectionKey || normalizedSelectionKeyFor(hitObject.userData.selection));
        if (!selectionKey) {
          return false;
        }
        activeClusterInfoSelectionKey = selectionKey;
        activeClusterInfoSprite = hitObject;
        tooltipEl.style.display = "block";
        tooltipEl.innerHTML = hitObject.userData.hovertext;
        tooltipEl.style.borderColor = tooltipColor;
        tooltipEl.style.color = tooltipColor;
        tooltipEl.style.boxShadow = `0 0 0 1px ${tooltipColor}`;
        updateClusterInfoTooltipPosition();
        return true;
      }

      function canvasClickWasDrag(event) {
        if (!canvasPointerDownInfo || event.button !== canvasPointerDownInfo.button) {
          return false;
        }
        const dx = Number(event.clientX) - Number(canvasPointerDownInfo.x);
        const dy = Number(event.clientY) - Number(canvasPointerDownInfo.y);
        const dragThresholdSquared = canvasPointerDownInfo.pointerType === "touch" ? 144.0 : 16.0;
        return (dx * dx + dy * dy) > dragThresholdSquared;
      }

      function showClusterInfoTooltipFromPointerEvent(event) {
        const hitObject = pickSprite(event);
        const clickedSelectionKey = hitObject && hitObject.userData
          ? normalizeMemberKey(hitObject.userData.selectionKey || normalizedSelectionKeyFor(hitObject.userData.selection))
          : "";
        if (!clickedSelectionKey || !showClusterInfoTooltip(hitObject)) {
          hideClusterInfoTooltip();
          setLocalHoveredClusterKey("");
          return false;
        }
        setLocalHoveredClusterKey(clickedSelectionKey);
        return true;
      }

      function handleCanvasTouchTap(event) {
        if (!canvasPointerDownInfo || canvasPointerDownInfo.pointerType !== "touch") {
          return false;
        }
        if (!event || event.type === "pointercancel" || event.pointerType !== "touch") {
          canvasPointerDownInfo = null;
          return false;
        }
        if (
          canvasPointerDownInfo.pointerId !== undefined
          && event.pointerId !== undefined
          && canvasPointerDownInfo.pointerId !== event.pointerId
        ) {
          return false;
        }
        if (canvasClickWasDrag(event)) {
          canvasPointerDownInfo = null;
          return false;
        }
        canvasPointerDownInfo = null;
        if (
          minimalModeEnabled
          || widgetPointerState
          || lassoState
          || lassoArmed
          || selectionBoxPointerState
          || manualLabelPointerState
          || skyViewDragState
          || skyViewPinchState
        ) {
          return false;
        }
        showClusterInfoTooltipFromPointerEvent(event);
        suppressNextCanvasClick = true;
        stopPointerEvent(event);
        return true;
      }

      function onCanvasClick(event) {
        if (suppressNextCanvasClick) {
          suppressNextCanvasClick = false;
          canvasPointerDownInfo = null;
          return;
        }
        if (canvasClickWasDrag(event)) {
          canvasPointerDownInfo = null;
          return;
        }
        canvasPointerDownInfo = null;
        if (minimalModeEnabled || widgetPointerState || lassoState || event.button !== 0) {
          hideClusterInfoTooltip();
          setLocalHoveredClusterKey("");
          return;
        }
        showClusterInfoTooltipFromPointerEvent(event);
      }

      function onCanvasDoubleClick(event) {
        if (cameraViewMode === "earth") {
          event.preventDefault();
          return;
        }
        if (minimalModeEnabled || widgetPointerState || lassoState || event.button !== 0) {
          return;
        }
        if (selectionBoxRayHitFromEvent(event)) {
          return;
        }
        const target = doubleClickTargetFromEvent(event);
        if (!target || !target.worldPoint || target.fallbackPlane) {
          resetCameraView();
          event.preventDefault();
          return;
        }
        if (target.manualLabelId) {
          return;
        }
        if (target.referenceFrameKey) {
          focusSelectionKey = target.referenceFrameKey;
          renderFrame(currentFrameIndex);
          recenterCameraTarget(new THREE.Vector3(0.0, 0.0, 0.0));
        } else {
          focusSelectionKey = "";
          recenterCameraTarget(target.worldPoint);
        }
        event.preventDefault();
      }

      function onCanvasWheel(event) {
        const deltaY = Number(event.deltaY);
        if (!Number.isFinite(deltaY) || Math.abs(deltaY) <= 1e-6) {
          return;
        }
        if (cameraViewMode === "earth") {
          if (!zoomEarthViewByWheelDelta(deltaY)) {
            return;
          }
          event.preventDefault();
          event.stopPropagation();
          if (typeof event.stopImmediatePropagation === "function") {
            event.stopImmediatePropagation();
          }
          return;
        }
        if (!initialZoomAnchorActive()) {
          return;
        }
        const speedScale = Math.max(globalScrollSpeed, 0.2);
        const clampedDelta = clampRange(deltaY, -240.0, 240.0);
        const zoomFactor = Math.exp(clampedDelta * 0.0015 * speedScale);
        if (!zoomCameraTowardPoint(currentZoomAnchorPoint, zoomFactor)) {
          return;
        }
        event.preventDefault();
        event.stopPropagation();
        if (typeof event.stopImmediatePropagation === "function") {
          event.stopImmediatePropagation();
        }
      }

      function onScaleBarPointerStart(event) {
        return false;
      }

      function updateScaleBarInteraction(event) {
        if (!scaleBarPointerState || !scaleBarEl) {
          return false;
        }
        const left = scaleBarPointerState.startLeft + (event.clientX - scaleBarPointerState.startX);
        const top = scaleBarPointerState.startTop + (event.clientY - scaleBarPointerState.startY);
        const next = clampScaleBarPosition(left, top, scaleBarPointerState.width, scaleBarPointerState.height);
        scaleBarPosition = { left: next.left, top: next.top };
        scaleBarEl.style.left = `${next.left}px`;
        scaleBarEl.style.top = `${next.top}px`;
        scaleBarEl.style.right = "auto";
        scaleBarEl.style.bottom = "auto";
        event.preventDefault();
        event.stopPropagation();
        return true;
      }

      function finishScaleBarInteraction(event) {
        if (!scaleBarPointerState || !scaleBarEl) {
          return false;
        }
        if (typeof scaleBarEl.releasePointerCapture === "function" && event.pointerId !== undefined) {
          try {
            scaleBarEl.releasePointerCapture(event.pointerId);
          } catch (_err) {
          }
        }
        scaleBarEl.dataset.dragging = "false";
        controls.enabled = true;
        document.body.style.userSelect = "";
        scaleBarPointerState = null;
        return true;
      }
""".strip()
