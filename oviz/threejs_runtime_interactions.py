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

      function setClusterSelections(selections, mode = "click") {
        if (minimalModeEnabled) {
          currentSelection = null;
          currentSelections = [];
          selectedClusterKeys = new Set();
          currentSelectionMode = "none";
          updateSelectionUI();
          return;
        }
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

      function setToolsDrawerOpen(isOpen) {
        if (!toolsShellEl || !toolsToggleEl) {
          return;
        }
        const nextOpen = minimalModeEnabled ? false : Boolean(isOpen);
        toolsShellEl.dataset.open = nextOpen ? "true" : "false";
        toolsToggleEl.textContent = nextOpen ? "Selection ▾" : "Selection ▸";
        if (nextOpen && controlsShellEl && controlsShellEl.dataset.open === "true") {
          controlsShellEl.dataset.open = "false";
          if (controlsToggleEl) {
            controlsToggleEl.textContent = "Controls ▸";
          }
        }
      }

      function setControlsDrawerOpen(isOpen) {
        if (!controlsShellEl || !controlsToggleEl) {
          return;
        }
        const nextOpen = minimalModeEnabled ? false : Boolean(isOpen);
        controlsShellEl.dataset.open = nextOpen ? "true" : "false";
        controlsToggleEl.textContent = nextOpen ? "Controls ▾" : "Controls ▸";
        if (nextOpen && toolsShellEl && toolsShellEl.dataset.open === "true") {
          toolsShellEl.dataset.open = "false";
          if (toolsToggleEl) {
            toolsToggleEl.textContent = "Selection ▸";
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

      function pause(options = {}) {
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
      }

      function play(direction = 1) {
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
        const clampedIndex = clampFrameIndex(index);
        pendingSliderFrameIndex = clampedIndex;
        displayedFrameValue = clampedIndex;
        currentFrameIndex = clampedIndex;
        updateTimelineUi(clampedIndex, frameTimeForValue(clampedIndex));
        if (sliderScrubRenderHandle !== null) {
          return;
        }
        sliderScrubRenderHandle = window.requestAnimationFrame(() => {
          sliderScrubRenderHandle = null;
          const targetIndex = pendingSliderFrameIndex;
          pendingSliderFrameIndex = null;
          if (targetIndex !== null && targetIndex !== undefined) {
            renderFrame(targetIndex);
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
            manualLabelId: hitObject.userData ? String(hitObject.userData.manualLabelId || "") : "",
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
              manualLabelId: "",
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
        if (minimalModeEnabled || widgetPointerState || !clickSelectionEnabled) {
          return;
        }
        if (selectionBoxRayHitFromEvent(event)) {
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
        if (minimalModeEnabled || widgetPointerState || lassoState || event.button !== 0) {
          return;
        }
        if (selectionBoxRayHitFromEvent(event)) {
          return;
        }
        const target = doubleClickTargetFromEvent(event);
        if (!target || !target.worldPoint) {
          return;
        }
        if (target.manualLabelId) {
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

      function onCanvasWheel(event) {
        if (galacticSimpleModeEnabled && galacticSimpleTracksOrbitTargetToSun) {
          enableGalacticSimpleOrbitTargetTracking();
        }
        if (!initialZoomAnchorActive()) {
          return;
        }
        const deltaY = Number(event.deltaY);
        if (!Number.isFinite(deltaY) || Math.abs(deltaY) <= 1e-6) {
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
        if (!scaleBarEl || event.button !== 0) {
          return;
        }
        const rect = scaleBarEl.getBoundingClientRect();
        const rootRect = root.getBoundingClientRect();
        scaleBarPointerState = {
          startX: event.clientX,
          startY: event.clientY,
          startLeft: rect.left - rootRect.left,
          startTop: rect.top - rootRect.top,
          width: rect.width,
          height: rect.height,
        };
        scaleBarEl.dataset.dragging = "true";
        controls.enabled = false;
        if (typeof scaleBarEl.setPointerCapture === "function" && event.pointerId !== undefined) {
          try {
            scaleBarEl.setPointerCapture(event.pointerId);
          } catch (_err) {
          }
        }
        document.body.style.userSelect = "none";
        focusViewer();
        event.preventDefault();
        event.stopPropagation();
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
