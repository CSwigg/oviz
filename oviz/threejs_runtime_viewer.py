from __future__ import annotations


THREEJS_VIEWER_RUNTIME_JS = """
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

      function clampManualLabelSize(value) {
        return clampRange(value, MIN_MANUAL_LABEL_SIZE, MAX_MANUAL_LABEL_SIZE);
      }

      function sanitizeManualLabelText(value) {
        return String(value ?? "")
          .replace(/\s+/g, " ")
          .trim()
          .slice(0, MAX_MANUAL_LABEL_TEXT_LENGTH);
      }

      function manualLabelFallbackText(index) {
        return `Label ${Math.max(1, Number(index) || 1)}`;
      }

      function trackManualLabelId(id, fallbackIndex = 0) {
        const numericFallback = Math.max(0, Number(fallbackIndex) || 0);
        const text = String(id || "");
        const match = text.match(/manual-label-(\d+)$/);
        if (match) {
          manualLabelIdCounter = Math.max(manualLabelIdCounter, Number(match[1]) || numericFallback);
        } else {
          manualLabelIdCounter = Math.max(manualLabelIdCounter, numericFallback);
        }
      }

      function nextManualLabelId() {
        manualLabelIdCounter += 1;
        return `manual-label-${manualLabelIdCounter}`;
      }

      function normalizeManualLabel(rawLabel, index) {
        if (!rawLabel || typeof rawLabel !== "object") {
          return null;
        }
        const fallbackIndex = index + 1;
        const rawId = String(rawLabel.id || "");
        const id = rawId || nextManualLabelId();
        trackManualLabelId(id, fallbackIndex);
        return {
          id,
          text: sanitizeManualLabelText(rawLabel.text) || manualLabelFallbackText(fallbackIndex),
          x: Number.isFinite(Number(rawLabel.x)) ? Number(rawLabel.x) : 0.0,
          y: Number.isFinite(Number(rawLabel.y)) ? Number(rawLabel.y) : 0.0,
          z: Number.isFinite(Number(rawLabel.z)) ? Number(rawLabel.z) : 0.0,
          size: clampManualLabelSize(rawLabel.size ?? DEFAULT_MANUAL_LABEL_SIZE),
        };
      }

      function manualLabelById(labelId) {
        const requestedId = String(labelId || "");
        return manualLabels.find((label) => String(label.id) === requestedId) || null;
      }

      function activeManualLabel() {
        return manualLabelById(activeManualLabelId);
      }

      function ensureActiveManualLabel(preferredId = "") {
        const requestedId = String(preferredId || activeManualLabelId || "");
        if (requestedId && manualLabelById(requestedId)) {
          activeManualLabelId = requestedId;
          return;
        }
        activeManualLabelId = manualLabels.length ? String(manualLabels[manualLabels.length - 1].id || "") : "";
      }

      function syncManualLabelDraftFromSelection() {
        const selectedLabel = activeManualLabel();
        if (selectedLabel) {
          manualLabelDraftText = String(selectedLabel.text || "");
          manualLabelDraftSize = clampManualLabelSize(selectedLabel.size);
          return;
        }
        manualLabelDraftText = "";
        manualLabelDraftSize = DEFAULT_MANUAL_LABEL_SIZE;
      }

      function loadManualLabels(labelState, preferredId = "") {
        manualLabelIdCounter = 0;
        manualLabels = Array.isArray(labelState)
          ? labelState
            .map((label, index) => normalizeManualLabel(label, index))
            .filter(Boolean)
          : [];
        ensureActiveManualLabel(preferredId);
        syncManualLabelDraftFromSelection();
      }

      function manualLabelScenePositionFromCurrentTarget() {
        return controls.target.clone().sub(plotGroup.position);
      }

      function currentManualLabelDraftText(fallbackText) {
        return sanitizeManualLabelText(manualLabelDraftText) || String(fallbackText || "");
      }

      function currentManualLabelDraftSize() {
        return clampManualLabelSize(manualLabelDraftSize);
      }

      function renderManualLabelControls() {
        if (!manualLabelSelectEl || !manualLabelTextEl || !manualLabelSizeEl) {
          return;
        }
        ensureActiveManualLabel();
        const selectedLabel = activeManualLabel();

        manualLabelSelectEl.innerHTML = "";
        const placeholder = document.createElement("option");
        placeholder.value = "";
        placeholder.textContent = manualLabels.length ? "Select label" : "No manual labels";
        manualLabelSelectEl.appendChild(placeholder);
        manualLabels.forEach((label, index) => {
          const option = document.createElement("option");
          option.value = String(label.id || "");
          option.textContent = `${index + 1}. ${String(label.text || manualLabelFallbackText(index + 1))}`;
          manualLabelSelectEl.appendChild(option);
        });
        manualLabelSelectEl.disabled = manualLabels.length === 0;
        manualLabelSelectEl.value = selectedLabel ? String(selectedLabel.id || "") : "";
        manualLabelTextEl.value = String(manualLabelDraftText || "");
        manualLabelSizeEl.value = String(Math.round(currentManualLabelDraftSize()));
        if (manualLabelApplyButtonEl) {
          manualLabelApplyButtonEl.disabled = !selectedLabel;
        }
        if (manualLabelDeleteButtonEl) {
          manualLabelDeleteButtonEl.disabled = !selectedLabel;
        }
        if (manualLabelReadoutEl) {
          if (selectedLabel) {
            manualLabelReadoutEl.textContent = `Selected: ${selectedLabel.text} @ (${formatTick(selectedLabel.x)}, ${formatTick(selectedLabel.y)}, ${formatTick(selectedLabel.z)}) pc. Click and drag in the scene to reposition.`;
          } else {
            manualLabelReadoutEl.textContent = "Add a label at the current camera target, then drag it in the scene.";
          }
        }
      }

      function addManualLabelFromDraft() {
        const position = manualLabelScenePositionFromCurrentTarget();
        const nextLabel = {
          id: nextManualLabelId(),
          text: currentManualLabelDraftText(manualLabelFallbackText(manualLabels.length + 1)),
          x: Number.isFinite(position.x) ? position.x : 0.0,
          y: Number.isFinite(position.y) ? position.y : 0.0,
          z: Number.isFinite(position.z) ? position.z : 0.0,
          size: currentManualLabelDraftSize(),
        };
        manualLabels.push(nextLabel);
        activeManualLabelId = nextLabel.id;
        syncManualLabelDraftFromSelection();
        renderManualLabelControls();
        renderFrame(currentFrameIndex);
      }

      function applyManualLabelDraftToSelection() {
        const selectedLabel = activeManualLabel();
        if (!selectedLabel) {
          return;
        }
        selectedLabel.text = currentManualLabelDraftText(selectedLabel.text || manualLabelFallbackText(1));
        selectedLabel.size = currentManualLabelDraftSize();
        syncManualLabelDraftFromSelection();
        renderManualLabelControls();
        renderFrame(currentFrameIndex);
      }

      function deleteActiveManualLabel() {
        if (!activeManualLabel()) {
          return;
        }
        manualLabels = manualLabels.filter((label) => String(label.id || "") !== String(activeManualLabelId || ""));
        ensureActiveManualLabel();
        syncManualLabelDraftFromSelection();
        renderManualLabelControls();
        renderFrame(currentFrameIndex);
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
        root.dataset.themeKey = activeThemeKey;
        applyThemeCssVars();
        applySceneBackground();
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
        globalPointGlowStrength = clampRange(globalPointGlowStrength, 0.0, 4.0);
        fadeInTimeMyr = Math.max(Number.isFinite(Number(fadeInTimeMyr)) ? Number(fadeInTimeMyr) : 0.0, 0.0);
        focusTraceKey = String(focusTraceKey || "");
        camera.fov = clampRange(camera.fov, 0.05, 120.0);
        controls.zoomSpeed = globalScrollSpeed;
        updateControlSensitivityForView();
        syncCameraAutoOrbitMode();
        camera.updateProjectionMatrix();
      }

      function updateControlSensitivityForView() {
        if (!controls || !camera) {
          return;
        }
        if (cameraViewMode === "earth") {
          const fovScale = clampRange((Number(camera.fov) || 90.0) / 90.0, 0.06, 1.35);
          controls.rotateSpeed = 0.0;
          controls.panSpeed = 0.55 * fovScale;
          return;
        }
        controls.rotateSpeed = 0.7;
        controls.panSpeed = 0.7;
      }

      function normalizeCameraAutoOrbitSpeed(value) {
        const speed = Number(value);
        if (!Number.isFinite(speed)) {
          return 1.0;
        }
        return Math.max(speed, 0.0);
      }

      function normalizeCameraAutoOrbitDirection(value) {
        if (typeof value === "string") {
          const normalized = String(value).trim().toLowerCase();
          if (["reverse", "clockwise", "cw", "-1", "negative"].includes(normalized)) {
            return -1.0;
          }
          if (["forward", "counterclockwise", "anticlockwise", "ccw", "1", "positive"].includes(normalized)) {
            return 1.0;
          }
        }
        const numeric = Number(value);
        if (Number.isFinite(numeric) && numeric < 0.0) {
          return -1.0;
        }
        return 1.0;
      }

      function cameraAutoOrbitAngularSpeed() {
        return CAMERA_AUTO_ORBIT_SPEED
          * normalizeCameraAutoOrbitSpeed(cameraAutoOrbitBaseSpeed)
          * Math.max(Number(cameraAutoOrbitSpeedMultiplier) || 1.0, 0.05)
          * normalizeCameraAutoOrbitDirection(cameraAutoOrbitDirection);
      }

      function syncCameraAutoOrbitUi() {
        orbitCameraButtons.forEach((buttonEl) => {
          buttonEl.dataset.active = cameraAutoOrbitEnabled ? "true" : "false";
          buttonEl.setAttribute("aria-pressed", cameraAutoOrbitEnabled ? "true" : "false");
          buttonEl.textContent = cameraAutoOrbitEnabled ? "Stop orbit" : "Orbit camera";
          buttonEl.title = cameraAutoOrbitEnabled
            ? "Stop rotating around the current camera target"
            : "Rotate around the current camera target";
        });
      }

      function galacticSimpleDefaultOrbitActive() {
        const initialViewOffset = initialCameraState && initialCameraState.viewOffset
          ? normalizeActionViewOffset(initialCameraState.viewOffset)
          : { x: 0.0, y: 0.0 };
        const framedInitialView = Math.abs(initialViewOffset.x) > 1e-6 || Math.abs(initialViewOffset.y) > 1e-6;
        return Boolean(
          cameraAutoOrbitEnabled
          && galacticSimpleModeEnabled
          && galacticSimpleTracksOrbitTargetToSun
          && !actionCameraOrbitOwned
          && !framedInitialView
        );
      }

      function syncCameraAutoOrbitMode() {
        const useControlsAutoRotate = cameraAutoOrbitEnabled && !galacticSimpleDefaultOrbitActive();
        controls.autoRotate = useControlsAutoRotate;
        controls.autoRotateSpeed = cameraAutoOrbitAngularSpeed();
      }

      function setCameraAutoOrbitEnabled(enabled) {
        cameraAutoOrbitEnabled = Boolean(enabled);
        syncCameraAutoOrbitMode();
        syncCameraAutoOrbitUi();
      }

      function syncEarthViewToggleUi() {
        const isEarthView = cameraViewMode === "earth";
        if (root && root.dataset) {
          root.dataset.cameraViewMode = cameraViewMode;
        }
        if (earthViewToggleButtonEl) {
          earthViewToggleButtonEl.textContent = isEarthView ? "Sky View" : "3D View";
          earthViewToggleButtonEl.dataset.active = isEarthView ? "true" : "false";
          earthViewToggleButtonEl.setAttribute("aria-pressed", isEarthView ? "true" : "false");
          earthViewToggleButtonEl.title = isEarthView
            ? "Return to the previous 3D view"
            : "Move to the sky view from the observer position";
        }
        if (viewFromEarthButtonEl) {
          viewFromEarthButtonEl.textContent = isEarthView ? "Sky View" : "3D View";
          viewFromEarthButtonEl.dataset.active = isEarthView ? "true" : "false";
          viewFromEarthButtonEl.setAttribute("aria-pressed", isEarthView ? "true" : "false");
          viewFromEarthButtonEl.title = isEarthView
            ? "Return to the previous 3D view"
            : "Move to the sky view from the observer position";
        }
      }

      function applyCameraViewMode() {
        const isEarthView = cameraViewMode === "earth";
        controls.enabled = !isEarthView;
        controls.enableRotate = !isEarthView;
        controls.enableZoom = !isEarthView;
        controls.enablePan = !isEarthView;
        updateControlSensitivityForView();
        if (skyDomeFrameEl && !isEarthView) {
          skyDomeFrameEl.style.opacity = "0";
        }
        if (typeof updateSkyDomeBackgroundFrame === "function") {
          updateSkyDomeBackgroundFrame(
            (typeof performance !== "undefined" && performance.now) ? performance.now() : Date.now(),
            { force: true }
          );
        }
        syncEarthViewToggleUi();
      }

      function trackedZoomAnchorPointForFrame(frame) {
        const frameTraces = frame && Array.isArray(frame.traces) ? frame.traces : [];
        for (const trace of frameTraces) {
          if (!trace || String(trace.name || "") !== "Sun") {
            continue;
          }
          const points = Array.isArray(trace.points) ? trace.points : [];
          if (!points.length) {
            continue;
          }
          const point = points[0] || {};
          const x = Number(point.x);
          const y = Number(point.y);
          const z = Number(point.z);
          if (Number.isFinite(x) && Number.isFinite(y) && Number.isFinite(z)) {
            return new THREE.Vector3(x, y, z);
          }
        }
        return initialZoomAnchorPoint ? initialZoomAnchorPoint.clone() : null;
      }

      function initialZoomAnchorActive() {
        return Boolean(
          currentZoomAnchorPoint
          && cameraViewMode === "free"
          && controls
          && controls.target
        );
      }

      function zoomCameraTowardPoint(anchorPoint, zoomFactor) {
        if (!(anchorPoint instanceof THREE.Vector3) || !Number.isFinite(Number(zoomFactor)) || zoomFactor <= 0.0) {
          return false;
        }
        const cameraOffset = camera.position.clone().sub(anchorPoint);
        const targetOffset = controls.target.clone().sub(anchorPoint);
        const viewOffset = controls.target.clone().sub(camera.position);
        const viewDistance = viewOffset.length();
        if (!(viewDistance > 1e-9)) {
          return false;
        }
        const maxSpan = Math.max(sceneSpec.max_span || 1, 1);
        const minDistance = Math.max(maxSpan * 1e-5, 0.05);
        const maxDistance = Math.max(maxSpan * 12.0, 10.0);
        const nextViewDistance = clampRange(viewDistance * zoomFactor, minDistance, maxDistance);
        const appliedFactor = nextViewDistance / Math.max(viewDistance, 1e-9);
        const preserveTargetOffset = Boolean(
          galacticSimpleModeEnabled
          && galacticSimpleTracksOrbitTargetToSun
        );
        camera.position.copy(anchorPoint.clone().add(cameraOffset.multiplyScalar(appliedFactor)));
        controls.target.copy(
          preserveTargetOffset
            ? anchorPoint.clone().add(targetOffset)
            : anchorPoint.clone().add(targetOffset.multiplyScalar(appliedFactor))
        );
        camera.up.set(sceneUp.x ?? 0.0, sceneUp.y ?? 0.0, sceneUp.z ?? 1.0);
        controls.update();
        updateScaleBar();
        return true;
      }

      function enableGalacticSimpleOrbitTargetTracking() {
        if (!galacticSimpleModeEnabled || !galacticSimpleTracksOrbitTargetToSun) {
          return false;
        }
        if (!(currentZoomAnchorPoint instanceof THREE.Vector3) || !controls || !controls.target) {
          return false;
        }
        galacticSimpleOrbitTargetTrackingActive = true;
        controls.target.copy(currentZoomAnchorPoint);
        return true;
      }

      function updateGalacticSimpleDefaultOrbit(deltaSeconds) {
        if (!galacticSimpleDefaultOrbitActive()) {
          return false;
        }
        if (!Number.isFinite(deltaSeconds) || deltaSeconds <= 0.0) {
          return false;
        }
        const anchorPoint = trackedZoomAnchorPointForFrame(currentFrame()) || currentZoomAnchorPoint || initialZoomAnchorPoint;
        if (!(anchorPoint instanceof THREE.Vector3)) {
          return false;
        }
        const angle = deltaSeconds
          * ((2.0 * Math.PI) / 60.0)
          * cameraAutoOrbitAngularSpeed();
        if (Math.abs(angle) <= 1e-12) {
          return false;
        }
        const cameraOffset = camera.position.clone().sub(anchorPoint);
        const targetOffset = controls.target.clone().sub(anchorPoint);
        if (cameraOffset.lengthSq() <= 1e-18) {
          return false;
        }
        cameraOffset.applyAxisAngle(sceneUpVector, angle);
        targetOffset.applyAxisAngle(sceneUpVector, angle);
        camera.position.copy(anchorPoint.clone().add(cameraOffset));
        controls.target.copy(anchorPoint.clone().add(targetOffset));
        camera.up.set(sceneUp.x ?? 0.0, sceneUp.y ?? 0.0, sceneUp.z ?? 1.0);
        return true;
      }

      function resetCameraView() {
        setCameraAutoOrbitEnabled(false);
        cameraViewMode = "free";
        earthViewFocusDistance = null;
        earthViewReturnCameraState = null;
        galacticSimpleOrbitTargetTrackingActive = false;
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

      function earthViewPoint() {
        return new THREE.Vector3(0.0, 0.0, 0.0);
      }

      function cameraDirectionForEarthView(targetPoint = null) {
        const direction = new THREE.Vector3();
        camera.getWorldDirection(direction);
        if (direction.lengthSq() > 1e-12) {
          return direction.normalize();
        }
        const earthPoint = earthViewPoint();
        const fallbackTarget = targetPoint || earthViewTargetPoint();
        direction.subVectors(fallbackTarget, earthPoint);
        if (direction.lengthSq() <= 1e-12) {
          direction.set(1.0, 0.0, 0.0);
        }
        return direction.normalize();
      }

      function leveledSkyViewDirectionForEarthView(targetPoint = null) {
        const direction = cameraDirectionForEarthView(targetPoint);
        direction.z = 0.0;
        if (direction.lengthSq() <= 1e-12) {
          const earthPoint = earthViewPoint();
          const fallbackTarget = targetPoint || earthViewTargetPoint();
          direction.subVectors(fallbackTarget, earthPoint);
          direction.z = 0.0;
        }
        if (direction.lengthSq() <= 1e-12) {
          direction.set(1.0, 0.0, 0.0);
        }
        return direction.normalize();
      }

      function skyViewUpVectorForDirection(direction) {
        const up = sceneUpVector && sceneUpVector.lengthSq && sceneUpVector.lengthSq() > 1e-12
          ? sceneUpVector.clone()
          : new THREE.Vector3(0.0, 0.0, 1.0);
        const safeDirection = direction && direction.lengthSq && direction.lengthSq() > 1e-12
          ? direction.clone().normalize()
          : new THREE.Vector3(1.0, 0.0, 0.0);
        up.sub(safeDirection.clone().multiplyScalar(up.dot(safeDirection)));
        if (up.lengthSq() <= 1e-12) {
          up.set(0.0, 0.0, 1.0);
        }
        return up.normalize();
      }

      function lockEarthViewCameraToTarget() {
        if (cameraViewMode !== "earth") {
          return;
        }
        const earthPoint = earthViewPoint();
        const direction = new THREE.Vector3().subVectors(earthPoint, camera.position);
        if (direction.lengthSq() <= 1e-12) {
          camera.getWorldDirection(direction);
        }
        if (direction.lengthSq() <= 1e-12) {
          direction.set(1.0, 0.0, 0.0);
        }
        direction.normalize();
        controls.target.copy(earthPoint);
        camera.up.copy(skyViewUpVectorForDirection(direction));
        camera.lookAt(earthPoint);
        camera.updateMatrixWorld(true);
      }

      function cameraFovForScaleBarLength(scaleBarPc, referenceDistancePc) {
        const targetScaleBarPc = Number(scaleBarPc);
        const referenceDistance = Number(referenceDistancePc);
        if (!Number.isFinite(targetScaleBarPc) || !(targetScaleBarPc > 0.0) || !Number.isFinite(referenceDistance) || !(referenceDistance > 0.0)) {
          return Number(camera.fov) || 90.0;
        }
        const canvasHeight = Math.max(canvas.clientHeight || root.clientHeight || 0, 1);
        const tangent = targetScaleBarPc * canvasHeight / Math.max(240.0 * referenceDistance, 1e-9);
        const fovDeg = 2.0 * Math.atan(Math.max(tangent, 1e-9)) * 180.0 / Math.PI;
        return clampRange(fovDeg, 0.05, 120.0);
      }

      function cameraDistanceForScaleBarLength(scaleBarPc, fovDeg) {
        const targetScaleBarPc = Number(scaleBarPc);
        const fov = Number(fovDeg);
        if (!Number.isFinite(targetScaleBarPc) || !(targetScaleBarPc > 0.0) || !Number.isFinite(fov) || !(fov > 0.0)) {
          return Math.max(camera.position.distanceTo(controls.target), 1.0);
        }
        const canvasHeight = Math.max(canvas.clientHeight || root.clientHeight || 0, 1);
        const tangent = Math.tan(THREE.MathUtils.degToRad(fov * 0.5));
        if (!(tangent > 1e-12)) {
          return Math.max(camera.position.distanceTo(controls.target), 1.0);
        }
        return Math.max(targetScaleBarPc * canvasHeight / (240.0 * 2.0 * tangent), 1e-6);
      }

      function cancelCameraTransition() {
        if (cameraTransitionAnimationFrame) {
          window.cancelAnimationFrame(cameraTransitionAnimationFrame);
          cameraTransitionAnimationFrame = 0;
        }
      }

      function cancelSkyDomeOpacityAnimation() {
        if (skyDomeOpacityAnimationFrame) {
          window.cancelAnimationFrame(skyDomeOpacityAnimationFrame);
          skyDomeOpacityAnimationFrame = 0;
        }
      }

      function setSkyDomeViewOpacityScale(value, options = {}) {
        skyDomeViewOpacityScale = Math.min(Math.max(Number(value) || 0.0, 0.0), 1.0);
        if (skyDomeFrameEl && skyDomeViewOpacityScale <= 0.002) {
          skyDomeFrameEl.style.opacity = "0";
        }
        const now = (typeof performance !== "undefined" && performance.now) ? performance.now() : Date.now();
        if (typeof updateSkyDome === "function") {
          updateSkyDome(now);
        }
        if (typeof updateSkyDomeBackgroundFrame === "function") {
          updateSkyDomeBackgroundFrame(now, { force: options.force !== false });
        }
        if (typeof refreshSkyDomeControlStatus === "function") {
          refreshSkyDomeControlStatus();
        }
      }

      function animateSkyDomeViewOpacity(targetOpacity, options = {}, onComplete = null) {
        cancelSkyDomeOpacityAnimation();
        const startOpacity = Math.min(Math.max(Number(skyDomeViewOpacityScale) || 0.0, 0.0), 1.0);
        const endOpacity = Math.min(Math.max(Number(targetOpacity) || 0.0, 0.0), 1.0);
        const durationMs = Math.max(Number(options.durationMs) || 420.0, 1.0);
        const startMs = (typeof performance !== "undefined" && performance.now) ? performance.now() : Date.now();
        const step = (timestampMs) => {
          const now = Number(timestampMs) || ((typeof performance !== "undefined" && performance.now) ? performance.now() : Date.now());
          const linear = clampRange((now - startMs) / durationMs, 0.0, 1.0);
          const eased = linear * linear * linear * (linear * (linear * 6.0 - 15.0) + 10.0);
          setSkyDomeViewOpacityScale(startOpacity + ((endOpacity - startOpacity) * eased), { force: true });
          if (linear < 1.0) {
            skyDomeOpacityAnimationFrame = window.requestAnimationFrame(step);
            return;
          }
          skyDomeOpacityAnimationFrame = 0;
          setSkyDomeViewOpacityScale(endOpacity, { force: true });
          if (typeof onComplete === "function") {
            onComplete();
          }
        };
        skyDomeOpacityAnimationFrame = window.requestAnimationFrame(step);
      }

      function animateCameraTransition(targetPosition, targetControlsTarget, targetFov, onComplete = null, options = {}) {
        cancelCameraTransition();
        const startPosition = camera.position.clone();
        const startTarget = controls.target.clone();
        const startUp = camera.up.clone();
        const startFov = Number(camera.fov) || 60.0;
        const endPosition = targetPosition.clone();
        const endTarget = targetControlsTarget.clone();
        const endUp = options && options.targetUp ? options.targetUp.clone() : startUp.clone();
        const endFov = clampRange(targetFov, 0.05, 120.0);
        const lockDirection = Boolean(options && options.lockDirection && options.direction);
        const lockedDirection = lockDirection
          ? options.direction.clone().normalize()
          : null;
        const startLockedDirection = lockDirection && options && options.startDirection
          ? options.startDirection.clone().normalize()
          : null;
        const startDistance = Math.max(startPosition.distanceTo(startTarget), 1e-6);
        const rawEndDistance = Number(options && options.endDistance);
        const endDistance = Number.isFinite(rawEndDistance) && rawEndDistance > 0.0
          ? rawEndDistance
          : Math.max(endPosition.distanceTo(endTarget), 1e-6);
        const startMs = (typeof performance !== "undefined" && performance.now) ? performance.now() : Date.now();
        const durationMs = Math.max(Number(options && options.durationMs) || 720.0, 1.0);
        const step = (timestampMs) => {
          const now = Number(timestampMs) || ((typeof performance !== "undefined" && performance.now) ? performance.now() : Date.now());
          const linear = clampRange((now - startMs) / durationMs, 0.0, 1.0);
          const eased = linear * linear * linear * (linear * (linear * 6.0 - 15.0) + 10.0);
          camera.position.copy(startPosition).lerp(endPosition, eased);
          if (lockDirection && lockedDirection) {
            const distance = startDistance + ((endDistance - startDistance) * eased);
            let stepDirection = lockedDirection;
            if (startLockedDirection) {
              stepDirection = startLockedDirection.clone().lerp(lockedDirection, eased);
              if (stepDirection.lengthSq() <= 1e-12) {
                stepDirection = lockedDirection;
              } else {
                stepDirection.normalize();
              }
            }
            controls.target.copy(camera.position).add(stepDirection.clone().multiplyScalar(distance));
          } else {
            controls.target.copy(startTarget).lerp(endTarget, eased);
          }
          camera.up.copy(startUp).lerp(endUp, eased);
          if (camera.up.lengthSq() <= 1e-12) {
            camera.up.copy(endUp);
          }
          camera.up.normalize();
          camera.fov = startFov + (endFov - startFov) * eased;
          camera.updateProjectionMatrix();
          updateControlSensitivityForView();
          controls.update();
          renderSceneControls();
          updateScaleBar();
          updateCameraResponsiveImagePlanes();
          updateSkyDomeBackgroundFrame(now, { force: true });
          if (linear < 1.0) {
            cameraTransitionAnimationFrame = window.requestAnimationFrame(step);
            return;
          }
          cameraTransitionAnimationFrame = 0;
          if (typeof onComplete === "function") {
            onComplete();
          }
        };
        cameraTransitionAnimationFrame = window.requestAnimationFrame(step);
      }

      function setCameraFovFromZoomFactor(zoomFactor) {
        const factor = Number(zoomFactor);
        if (!Number.isFinite(factor) || !(factor > 0.0)) {
          return false;
        }
        const previousFov = Number(camera.fov) || 90.0;
        const nextFov = clampRange(previousFov * factor, 0.05, 120.0);
        if (Math.abs(nextFov - previousFov) <= 1e-6) {
          return false;
        }
        cancelCameraTransition();
        camera.fov = nextFov;
        camera.updateProjectionMatrix();
        updateControlSensitivityForView();
        if (cameraViewMode === "earth") {
          lockEarthViewCameraToTarget();
        } else {
          controls.update();
        }
        renderSceneControls();
        updateScaleBar();
        updateCameraResponsiveImagePlanes();
        updateSkyDomeBackgroundFrame(
          (typeof performance !== "undefined" && performance.now) ? performance.now() : Date.now(),
          { force: true }
        );
        return true;
      }

      function zoomEarthViewByWheelDelta(deltaY) {
        const numericDelta = Number(deltaY);
        if (!Number.isFinite(numericDelta) || Math.abs(numericDelta) <= 1e-6) {
          return false;
        }
        const speedScale = Math.max(globalScrollSpeed, 0.2);
        const clampedDelta = clampRange(numericDelta, -240.0, 240.0);
        const zoomFactor = Math.exp(clampedDelta * 0.0015 * speedScale);
        return setCameraFovFromZoomFactor(zoomFactor);
      }

      function stopPointerEvent(event) {
        if (!event) {
          return;
        }
        event.preventDefault();
        event.stopPropagation();
        if (typeof event.stopImmediatePropagation === "function") {
          event.stopImmediatePropagation();
        }
      }

      function rotateSkyViewCameraByPixels(deltaX, deltaY) {
        const dx = Number(deltaX);
        const dy = Number(deltaY);
        if (!Number.isFinite(dx) || !Number.isFinite(dy) || (Math.abs(dx) <= 1e-6 && Math.abs(dy) <= 1e-6)) {
          return false;
        }
        const height = Math.max(Number(canvas && canvas.clientHeight) || 1.0, 1.0);
        const width = Math.max(Number(canvas && canvas.clientWidth) || height, 1.0);
        const verticalFovRad = THREE.MathUtils.degToRad(clampRange(Number(camera.fov) || 90.0, 0.05, 120.0));
        const horizontalFovRad = 2.0 * Math.atan(Math.tan(verticalFovRad * 0.5) * (width / height));
        const direction = new THREE.Vector3();
        camera.getWorldDirection(direction);
        if (direction.lengthSq() <= 1e-12) {
          return false;
        }
        direction.normalize();
        const up = camera.up.clone();
        if (up.lengthSq() <= 1e-12) {
          up.copy(skyViewUpVectorForDirection(direction));
        }
        up.normalize();
        const right = new THREE.Vector3().crossVectors(up, direction);
        if (right.lengthSq() <= 1e-12) {
          return false;
        }
        right.normalize();
        const xOffset = (2.0 * dx / width) * Math.tan(horizontalFovRad * 0.5);
        const yOffset = (2.0 * dy / height) * Math.tan(verticalFovRad * 0.5);
        const nextDirection = direction.clone()
          .addScaledVector(right, xOffset)
          .addScaledVector(up, yOffset);
        if (nextDirection.lengthSq() <= 1e-12) {
          return false;
        }
        nextDirection.normalize();
        const earthPoint = earthViewPoint();
        const orbitDistance = Math.max(camera.position.distanceTo(earthPoint), 1e-6);
        controls.target.copy(earthPoint);
        camera.position.copy(earthPoint).sub(nextDirection.clone().multiplyScalar(orbitDistance));
        camera.up.copy(skyViewUpVectorForDirection(nextDirection));
        camera.lookAt(earthPoint);
        camera.updateMatrixWorld(true);
        updateScaleBar();
        updateCameraResponsiveImagePlanes();
        updateSkyDomeBackgroundFrame(
          (typeof performance !== "undefined" && performance.now) ? performance.now() : Date.now(),
          { force: true }
        );
        return true;
      }

      function startSkyViewCameraDrag(event) {
        if (
          cameraViewMode !== "earth"
          || widgetPointerState
          || selectionBoxPointerState
          || manualLabelPointerState
          || lassoState
          || !event
          || event.button !== 0
          || event.shiftKey
          || lassoArmed
        ) {
          return false;
        }
        skyViewDragState = {
          pointerId: event.pointerId,
          lastX: Number(event.clientX) || 0.0,
          lastY: Number(event.clientY) || 0.0,
          moved: false,
        };
        controls.enabled = false;
        setSkyDomeBackgroundCameraActive(true);
        document.body.style.userSelect = "none";
        if (typeof canvas.setPointerCapture === "function" && event.pointerId !== undefined) {
          try {
            canvas.setPointerCapture(event.pointerId);
          } catch (_err) {
          }
        }
        tooltipEl.style.display = "none";
        stopPointerEvent(event);
        return true;
      }

      function updateSkyViewCameraDrag(event) {
        if (!skyViewDragState || !event) {
          return false;
        }
        if (skyViewDragState.pointerId !== undefined && event.pointerId !== undefined && skyViewDragState.pointerId !== event.pointerId) {
          return true;
        }
        const nextX = Number(event.clientX) || skyViewDragState.lastX;
        const nextY = Number(event.clientY) || skyViewDragState.lastY;
        const dx = nextX - skyViewDragState.lastX;
        const dy = nextY - skyViewDragState.lastY;
        skyViewDragState.lastX = nextX;
        skyViewDragState.lastY = nextY;
        if (rotateSkyViewCameraByPixels(dx, dy)) {
          skyViewDragState.moved = true;
        }
        tooltipEl.style.display = "none";
        stopPointerEvent(event);
        return true;
      }

      function finishSkyViewCameraDrag(event) {
        if (!skyViewDragState) {
          return false;
        }
        if (typeof canvas.releasePointerCapture === "function" && skyViewDragState.pointerId !== undefined) {
          try {
            canvas.releasePointerCapture(skyViewDragState.pointerId);
          } catch (_err) {
          }
        }
        const moved = Boolean(skyViewDragState.moved);
        skyViewDragState = null;
        controls.enabled = cameraViewMode !== "earth";
        setSkyDomeBackgroundCameraActive(false);
        document.body.style.userSelect = "";
        if (moved) {
          suppressNextCanvasClick = true;
        }
        if (event) {
          stopPointerEvent(event);
        }
        return true;
      }

      function captureEarthViewReturnCameraState() {
        return {
          position: camera.position.clone(),
          target: controls.target.clone(),
          up: camera.up.clone(),
          fov: Number(camera.fov) || 60.0,
        };
      }

      function resetToSunReferenceFrameForSkyView() {
        const hasSelectionFocus = Boolean(normalizeMemberKey(focusSelectionKey));
        const hasTraceFocus = Boolean(String(focusTraceKey || ""));
        if (!hasSelectionFocus && !hasTraceFocus) {
          return false;
        }
        const previousPlotOffset = plotGroup.position.clone();
        focusSelectionKey = "";
        focusTraceKey = "";
        renderFrame(currentFrameIndex);
        const nextPlotOffset = plotGroup.position.clone();
        const worldDelta = nextPlotOffset.clone().sub(previousPlotOffset);
        if (worldDelta.lengthSq() > 1e-18) {
          camera.position.add(worldDelta);
          controls.target.add(worldDelta);
          controls.update();
          camera.updateMatrixWorld(true);
        }
        renderSceneControls();
        updateScaleBar();
        updateCameraResponsiveImagePlanes();
        return true;
      }

      function fallbackEarthViewReturnCameraState() {
        return {
          position: initialCameraState.position.clone(),
          target: initialCameraState.target.clone(),
          up: initialCameraState.up.clone(),
          fov: Number(initialCameraState.fov) || 60.0,
        };
      }

      function enterEarthViewFromCurrentCamera(options = {}) {
        const wasEarthView = cameraViewMode === "earth";
        const transitionSerial = ++skyViewTransitionSerial;
        cancelSkyDomeOpacityAnimation();
        setSkyDomeViewOpacityScale(0.0, { force: true });
        if (!wasEarthView && (!options || options.storeReturnState !== false)) {
          earthViewReturnCameraState = captureEarthViewReturnCameraState();
        }
        const earthPoint = earthViewPoint();
        let targetPoint = earthViewTargetPoint();
        if (!targetPoint || !Number.isFinite(targetPoint.x) || !Number.isFinite(targetPoint.y) || !Number.isFinite(targetPoint.z)) {
          targetPoint = new THREE.Vector3(8122.0, 0.0, 0.0);
        }
        const previousScaleBarPc = scaleBarLengthPcForCurrentView();
        const preserveScaleBar = Boolean(options && options.preserveScaleBar);
        const preserveDirection = !options || options.preserveDirection !== false;
        const preserveFov = Boolean(options && options.preserveFov);
        const startDirection = preserveDirection
          ? cameraDirectionForEarthView(targetPoint)
          : null;
        const direction = preserveDirection
          ? leveledSkyViewDirectionForEarthView(targetPoint)
          : new THREE.Vector3().subVectors(targetPoint, earthPoint);
        if (!preserveDirection) {
          direction.z = 0.0;
        }
        if (direction.lengthSq() <= 1e-12) {
          direction.set(1.0, 0.0, 0.0);
        }
        direction.normalize();
        const focusDistance = Math.max(targetPoint.distanceTo(earthPoint), 1e-6);
        const orbitRadius = Math.max(1e-3, Math.min(0.05, focusDistance * 1e-6));
        const targetPosition = earthPoint.clone().sub(direction.clone().multiplyScalar(orbitRadius));
        const targetFov = preserveFov
          ? clampRange(Number(camera.fov) || 60.0, 0.05, 120.0)
          : (
            preserveScaleBar
              ? cameraFovForScaleBarLength(previousScaleBarPc, focusDistance)
              : 90.0
          );
        const targetUp = skyViewUpVectorForDirection(direction);
        cameraViewMode = "earth";
        earthViewFocusDistance = focusDistance;
        applyGlobalControlState();
        applyCameraViewMode();
        animateCameraTransition(
          targetPosition,
          earthPoint,
          targetFov,
          () => {
            controls.target.copy(earthPoint);
            camera.position.copy(targetPosition);
            camera.up.copy(targetUp);
            camera.fov = targetFov;
            renderFrame(currentFrameIndex);
            applyGlobalControlState();
            applyCameraViewMode();
            controls.update();
            renderSceneControls();
            updateScaleBar();
            updateCameraResponsiveImagePlanes();
            updateSkyDomeBackgroundFrame(
              (typeof performance !== "undefined" && performance.now) ? performance.now() : Date.now(),
              { force: true }
            );
            window.setTimeout(() => {
              if (transitionSerial !== skyViewTransitionSerial || cameraViewMode !== "earth") {
                return;
              }
              animateSkyDomeViewOpacity(1.0, { durationMs: 420.0 });
            }, 70);
          },
          {
            lockDirection: preserveDirection,
            direction,
            startDirection,
            endDistance: orbitRadius,
            durationMs: 820.0,
            targetUp,
          }
        );
      }

      function exitEarthView() {
        if (cameraViewMode !== "earth") {
          return false;
        }
        const transitionSerial = ++skyViewTransitionSerial;
        const returnState = earthViewReturnCameraState || fallbackEarthViewReturnCameraState();
        const targetPosition = returnState.position.clone();
        const targetControlsTarget = returnState.target.clone();
        const targetUp = returnState.up.clone();
        const targetFov = Number(returnState.fov) || 60.0;
        const exitTargetDirection = new THREE.Vector3().subVectors(targetControlsTarget, targetPosition);
        const exitTargetDistance = Math.max(exitTargetDirection.length(), 1e-6);
        if (exitTargetDirection.lengthSq() <= 1e-12) {
          camera.getWorldDirection(exitTargetDirection);
        }
        if (exitTargetDirection.lengthSq() <= 1e-12) {
          exitTargetDirection.set(1.0, 0.0, 0.0);
        }
        exitTargetDirection.normalize();
        animateSkyDomeViewOpacity(0.0, { durationMs: 360.0 }, () => {
          if (transitionSerial !== skyViewTransitionSerial) {
            return;
          }
          const exitStartDirection = new THREE.Vector3();
          camera.getWorldDirection(exitStartDirection);
          if (exitStartDirection.lengthSq() <= 1e-12) {
            exitStartDirection.subVectors(controls.target, camera.position);
          }
          if (exitStartDirection.lengthSq() <= 1e-12) {
            exitStartDirection.copy(exitTargetDirection);
          }
          exitStartDirection.normalize();
          setSkyDomeViewOpacityScale(0.0, { force: true });
          animateCameraTransition(
            targetPosition,
            targetControlsTarget,
            targetFov,
            () => {
              cameraViewMode = "free";
              earthViewFocusDistance = null;
              controls.target.copy(targetControlsTarget);
              camera.position.copy(targetPosition);
              camera.up.copy(targetUp);
              camera.fov = targetFov;
              earthViewReturnCameraState = null;
              setSkyDomeViewOpacityScale(0.0, { force: true });
              applyGlobalControlState();
              applyCameraViewMode();
              controls.update();
              renderSceneControls();
              updateScaleBar();
              updateCameraResponsiveImagePlanes();
            },
            {
              durationMs: 820.0,
              targetUp,
              lockDirection: true,
              direction: exitTargetDirection,
              startDirection: exitStartDirection,
              endDistance: exitTargetDistance,
            }
          );
        });
        return true;
      }

      function viewFromEarth() {
        resetToSunReferenceFrameForSkyView();
        enterEarthViewFromCurrentCamera({ preserveScaleBar: false, preserveDirection: true, preserveFov: false });
      }

      function toggleEarthView() {
        if (cameraViewMode === "earth") {
          return exitEarthView();
        }
        viewFromEarth();
        return true;
      }

      function skyDomeControlsUseNativeHips() {
        const mode = String(
          (skyDomeSpec && (
            skyDomeSpec.background_mode
            || skyDomeSpec.mode
            || skyDomeSpec.render_mode
          ))
          || ""
        ).toLowerCase();
        const source = String(skyDomeSpec && skyDomeSpec.source || "").toLowerCase();
        return Boolean(
          skyDomeSpec
          && (
            source === "hips"
            || source === "native_hips"
            || source === "native-hips"
            || mode === "native_hips"
            || mode === "native-hips"
            || mode === "hips"
          )
        );
      }

      function skyDomeControlsAvailable() {
        return Boolean(
          skyDomeSpec
          && typeof skyDomeSpec === "object"
          && (
            Object.prototype.hasOwnProperty.call(skyDomeSpec, "enabled")
            || skyDomeSpec.source
            || skyDomeSpec.background_mode
            || skyDomeSpec.hips_base_url
            || skyDomeHasLocalSources()
          )
        );
      }

      function refreshSkyDomeControlStatus() {
        if (!skyDomeStatusEl) {
          return;
        }
        if (!skyDomeControlsAvailable()) {
          skyDomeStatusEl.textContent = "";
          skyDomeStatusEl.hidden = true;
          return;
        }
        const mode = root && root.dataset ? String(root.dataset.skyDomeMode || "") : "";
        const status = root && root.dataset ? String(root.dataset.skyDomeSnapshotStatus || skyDomeSnapshotStatus || "idle") : String(skyDomeSnapshotStatus || "idle");
        const message = root && root.dataset ? String(root.dataset.skyDomeSnapshotMessage || "") : "";
        const liveOpacity = root && root.dataset && root.dataset.skyDomeOpacity
          ? Number(root.dataset.skyDomeOpacity)
          : skyDomeOpacityForCurrentView();
        const parts = [
          `Status: ${status}`,
          `Opacity: ${Number.isFinite(liveOpacity) ? liveOpacity.toFixed(2) : "0.00"}`,
        ];
        if (skyDomeControlsUseNativeHips()) {
          parts.push(`Source: ${skyDomeHipsSurveyName()} HiPS`);
          const loadedTiles = root && root.dataset ? String(root.dataset.skyDomeHipsLoadedTiles || "0") : "0";
          const activeTiles = root && root.dataset ? String(root.dataset.skyDomeHipsActiveTiles || "0") : "0";
          const loadingTiles = root && root.dataset ? String(root.dataset.skyDomeHipsLoadingTiles || "0") : "0";
          const pendingTiles = root && root.dataset ? String(root.dataset.skyDomeHipsPendingTiles || "0") : "0";
          parts.push(`Tiles: ${loadedTiles}/${activeTiles} loaded, ${loadingTiles} loading, ${pendingTiles} pending`);
        } else if (mode) {
          parts.push(`Mode: ${mode}`);
        }
        const lowerStatus = status.toLowerCase();
        const quietStatus = lowerStatus === "loaded" || lowerStatus === "idle" || lowerStatus === "disabled";
        if (!message && quietStatus && !skyDomeControlsUseNativeHips()) {
          skyDomeStatusEl.textContent = "";
          skyDomeStatusEl.hidden = true;
          return;
        }
        skyDomeStatusEl.hidden = false;
        if (message) {
          skyDomeStatusEl.textContent = message;
        } else if (lowerStatus.includes("error") || lowerStatus.includes("unavailable")) {
          skyDomeStatusEl.textContent = "Sky image unavailable.";
        } else if (lowerStatus.includes("loading")) {
          skyDomeStatusEl.textContent = "Loading sky image...";
        } else {
          skyDomeStatusEl.textContent = skyDomeControlsUseNativeHips() ? parts.join(" | ") : "";
          skyDomeStatusEl.hidden = !skyDomeStatusEl.textContent;
        }
      }

      function syncSkyDomeControls() {
        if (!skyDomeControlsEl) {
          return;
        }
        const available = skyDomeControlsAvailable();
        skyDomeControlsEl.hidden = !available;
        if (!available) {
          return;
        }
        const nativeHips = skyDomeControlsUseNativeHips();
        const stretchControlsAvailable = (
          nativeHips
          || (typeof skyDomeUsesAladinBackground === "function" && skyDomeUsesAladinBackground())
          || (typeof skyDomeUsesHips2Fits === "function" && skyDomeUsesHips2Fits())
        );
        skyDomeHipsControlEls.forEach((controlEl) => {
          controlEl.hidden = !stretchControlsAvailable;
        });
        if (skyDomeVisibleToggleEl) {
          skyDomeVisibleToggleEl.checked = skyDomeIsEnabled();
        }
        if (skyDomeForceVisibleToggleEl) {
          skyDomeForceVisibleToggleEl.checked = Boolean(skyDomeForceVisible);
        }
        if (skyDomeOpacityEl) {
          skyDomeOpacityEl.value = String(Math.min(Math.max(Number(skyDomeSpec.opacity ?? 0.55), 0.0), 1.0));
        }
        if (skyDomeOpacityLabelEl) {
          skyDomeOpacityLabelEl.textContent = `Opacity (${Math.min(Math.max(Number(skyDomeSpec.opacity ?? 0.55), 0.0), 1.0).toFixed(2)})`;
        }
        if (skyDomeBrightnessEl) {
          skyDomeBrightnessEl.value = String(skyDomeHipsBrightness());
        }
        if (skyDomeBrightnessLabelEl) {
          skyDomeBrightnessLabelEl.textContent = `Brightness (${skyDomeHipsBrightness().toFixed(2)}x)`;
        }
        if (skyDomeContrastEl) {
          skyDomeContrastEl.value = String(skyDomeHipsContrast());
        }
        if (skyDomeContrastLabelEl) {
          skyDomeContrastLabelEl.textContent = `Contrast (${skyDomeHipsContrast().toFixed(2)}x)`;
        }
        if (skyDomeGammaEl) {
          skyDomeGammaEl.value = String(skyDomeHipsGamma());
        }
        if (skyDomeGammaLabelEl) {
          skyDomeGammaLabelEl.textContent = `Gamma (${skyDomeHipsGamma().toFixed(2)})`;
        }
        if (typeof syncSkyLayerControls === "function") {
          syncSkyLayerControls();
        }
        refreshSkyDomeControlStatus();
      }

      function updateSkyDomeFromControls(options = {}) {
        if (skyDomeSpec.enabled && !skyDomeMesh) {
          initializeSkyDomeFromSceneSpec();
        }
        if (skyDomeUsesNativeHips() && skyDomeHipsState) {
          if (skyDomeHipsState.allskyMesh) {
            applyNativeHipsMaterialSettings(
              skyDomeHipsState.allskyMesh.material,
              skyDomeHipsState.allskyTexture,
              skyDomeOpacityForCurrentView()
            );
          }
          skyDomeHipsState.tileCache.forEach((entry) => {
            if (entry && entry.mesh) {
              applyNativeHipsMaterialSettings(entry.mesh.material, entry.texture, skyDomeOpacityForCurrentView());
            }
          });
          updateNativeHipsTiles(
            (typeof performance !== "undefined" && performance.now) ? performance.now() : Date.now(),
            Boolean(options.forceTiles)
          );
        }
        if (typeof applySkyDomeFrameVisualSettings === "function") {
          applySkyDomeFrameVisualSettings();
        }
        updateSkyDome((typeof performance !== "undefined" && performance.now) ? performance.now() : Date.now());
        if (options.syncControls !== false) {
          syncSkyDomeControls();
        }
      }

      function renderSceneControls() {
        syncEarthViewToggleUi();
        const cameraFovValue = Number.isFinite(Number(camera.fov)) ? Number(camera.fov) : 60.0;
        const cameraFovLabel = cameraFovValue < 1.0
          ? cameraFovValue.toFixed(2).replace(/0+$/, "").replace(/\.$/, "")
          : String(Math.round(cameraFovValue));
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
          cameraFovLabelEl.textContent = `Camera FOV (${cameraFovLabel} deg)`;
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
        if (globalPointGlowEl) {
          globalPointGlowEl.value = String(globalPointGlowStrength);
        }
        if (globalPointGlowLabelEl) {
          globalPointGlowLabelEl.textContent = `Star glow (${globalPointGlowStrength.toFixed(2)}x)`;
        }
        if (sizeByStarsToggleEl) {
          sizeByStarsToggleEl.checked = sizePointsByStarsEnabled;
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
        if (galacticReferenceToggleEl) {
          galacticReferenceToggleEl.checked = galacticReferenceVisible;
        }
        nearbyRegionLabelsToggleEls.forEach((toggleEl) => {
          toggleEl.checked = nearbyRegionLabelsVisible;
        });
        syncSkyDomeControls();
        renderManualLabelControls();
        syncCameraAutoOrbitUi();
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
          if (typeof setSkyControlsDrawerOpen === "function") {
            setSkyControlsDrawerOpen(false);
          }
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
        if (minimalModeEnabled) {
          return false;
        }
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
          earthViewReturnCameraState = null;
          applyCameraViewMode();
          renderFrame(currentFrameIndex);
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
        if (galacticSimpleModeEnabled && galacticSimpleTracksOrbitTargetToSun) {
          enableGalacticSimpleOrbitTargetTracking();
        }
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
        if (galacticSimpleModeEnabled && galacticSimpleTracksOrbitTargetToSun) {
          enableGalacticSimpleOrbitTargetTracking();
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
        if (cameraViewMode === "earth") {
          const speedScale = Math.max(globalScrollSpeed, 0.2);
          setCameraFovFromZoomFactor(Math.exp(sign * deltaSeconds * 1.8 * speedScale));
          return;
        }
        if (initialZoomAnchorActive()) {
          const speedScale = Math.max(globalScrollSpeed, 0.2);
          const zoomFactor = Math.exp(sign * deltaSeconds * 1.8 * speedScale);
          if (zoomCameraTowardPoint(currentZoomAnchorPoint, zoomFactor)) {
            return;
          }
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
        const interruptsAction = isMovementKey
          || key === " "
          || key === "ArrowLeft"
          || key === "ArrowRight"
          || /^[1-9]$/.test(key)
          || ["l", "c", "v", "o"].includes(lowerKey);

        if (interruptsAction && !actionInterruptsMuted()) {
          interruptActionRun("keyboard", { disableOrbit: true });
        }

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

        if (event.repeat && (key === " " || key === "Escape" || /^[1-9]$/.test(key) || lowerKey === "l" || lowerKey === "c" || lowerKey === "v" || lowerKey === "o" || key === "?" || (key === "/" && event.shiftKey))) {
          event.preventDefault();
          return;
        }

        if (!minimalModeEnabled && (key === "?" || (key === "/" && event.shiftKey))) {
          clearPressedKeys();
          setKeyHelpOpen(true);
          event.preventDefault();
          return;
        }

        if (key === "Escape") {
          clearPressedKeys();
          if (activeLegendEditorKey) {
            closeLegendPopover();
            renderLegend();
            event.preventDefault();
            return;
          }
          clearClusterSelections();
          lassoArmed = false;
          updateSelectionUI();
          event.preventDefault();
          return;
        }

        if (key === " ") {
          if (playbackDirection !== 0) {
            pause();
          } else {
            play(1);
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

        if (minimalModeEnabled) {
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
        if (lowerKey === "o") {
          setCameraAutoOrbitEnabled(!cameraAutoOrbitEnabled);
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
        const barLengthPc = scaleBarLengthPcForCurrentView();
        currentScaleBarLengthPc = barLengthPc;
        scaleLabelEl.textContent = formatDistanceLabelPc(barLengthPc);
        scaleBarEl.style.display = Number.isFinite(barLengthPc) ? "flex" : "none";
        if (Number.isFinite(barLengthPc)) {
          applyScaleBarPosition();
        }
      }
""".strip()
