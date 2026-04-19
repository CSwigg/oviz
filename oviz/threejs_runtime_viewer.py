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
        scene.background = new THREE.Color(theme.scene_bgcolor || theme.paper_bgcolor || "#000000");
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
        camera.fov = clampRange(camera.fov, 18.0, 90.0);
        controls.zoomSpeed = globalScrollSpeed;
        controls.autoRotate = cameraAutoOrbitEnabled;
        controls.autoRotateSpeed = CAMERA_AUTO_ORBIT_SPEED;
        camera.updateProjectionMatrix();
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

      function setCameraAutoOrbitEnabled(enabled) {
        cameraAutoOrbitEnabled = Boolean(enabled);
        controls.autoRotate = cameraAutoOrbitEnabled;
        controls.autoRotateSpeed = CAMERA_AUTO_ORBIT_SPEED;
        syncCameraAutoOrbitUi();
      }

      function applyCameraViewMode() {
        const isEarthView = cameraViewMode === "earth";
        controls.enableRotate = true;
        controls.enableZoom = !isEarthView;
        controls.enablePan = !isEarthView;
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
        camera.position.copy(anchorPoint.clone().add(cameraOffset.multiplyScalar(appliedFactor)));
        controls.target.copy(anchorPoint.clone().add(targetOffset.multiplyScalar(appliedFactor)));
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

      function resetCameraView() {
        setCameraAutoOrbitEnabled(false);
        cameraViewMode = "free";
        earthViewFocusDistance = null;
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

      function viewFromEarth() {
        const earthPoint = new THREE.Vector3(0.0, 0.0, 20.8);
        let targetPoint = earthViewTargetPoint();
        if (!targetPoint || !Number.isFinite(targetPoint.x) || !Number.isFinite(targetPoint.y) || !Number.isFinite(targetPoint.z)) {
          targetPoint = new THREE.Vector3(8122.0, 0.0, 0.0);
        }
        const direction = new THREE.Vector3().subVectors(targetPoint, earthPoint);
        if (direction.lengthSq() <= 1e-12) {
          direction.set(1.0, 0.0, 0.0);
          targetPoint = earthPoint.clone().add(direction);
        }
        earthViewFocusDistance = Math.max(direction.length(), 1e-6);
        direction.normalize();
        const orbitRadius = Math.max(1e-3, Math.min(0.05, earthViewFocusDistance * 1e-6));
        cameraViewMode = "earth";
        controls.target.copy(earthPoint);
        camera.position.copy(earthPoint.clone().sub(direction.clone().multiplyScalar(orbitRadius)));
        camera.up.set(sceneUp.x ?? 0.0, sceneUp.y ?? 0.0, sceneUp.z ?? 1.0);
        camera.fov = 90.0;
        applyGlobalControlState();
        applyCameraViewMode();
        controls.update();
        renderSceneControls();
        updateScaleBar();
      }

      function renderSceneControls() {
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
          cameraFovLabelEl.textContent = `Camera FOV (${Math.round(Number(camera.fov))} deg)`;
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
          applyCameraViewMode();
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
        if (galacticSimpleModeEnabled && galacticSimpleTracksOrbitTargetToSun) {
          enableGalacticSimpleOrbitTargetTracking();
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
