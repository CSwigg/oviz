from __future__ import annotations


THREEJS_WIDGET_RUNTIME_JS = """
      function widgetPanelForKey(widgetKey) {
        if (widgetKey === "sky") {
          return skyPanelEl;
        }
        if (widgetKey === "box_metrics") {
          return boxMetricsPanelEl;
        }
        if (widgetKey === "age_kde") {
          return ageKdePanelEl;
        }
        if (widgetKey === "cluster_filter") {
          return clusterFilterPanelEl;
        }
        if (widgetKey === "dendrogram") {
          return dendrogramPanelEl;
        }
        return null;
      }

      function widgetEnabled(widgetKey) {
        if (widgetKey === "sky") {
          return Boolean(skySpec.enabled);
        }
        if (widgetKey === "box_metrics") {
          return Boolean(selectionBoxSpec.enabled);
        }
        if (widgetKey === "age_kde") {
          return Boolean(ageKdeSpec.enabled);
        }
        if (widgetKey === "cluster_filter") {
          return Boolean(clusterFilterSpec.enabled);
        }
        if (widgetKey === "dendrogram") {
          return Boolean(dendrogramSpec.enabled);
        }
        return false;
      }

      function widgetModeForKey(widgetKey) {
        if (widgetKey === "sky") {
          return skyPanelMode;
        }
        if (widgetKey === "box_metrics") {
          return boxMetricsPanelMode;
        }
        if (widgetKey === "age_kde") {
          return ageKdePanelMode;
        }
        if (widgetKey === "cluster_filter") {
          return clusterFilterPanelMode;
        }
        if (widgetKey === "dendrogram") {
          return dendrogramPanelMode;
        }
        return "hidden";
      }

      function setWidgetModeValue(widgetKey, mode) {
        if (widgetKey === "sky") {
          skyPanelMode = mode;
        } else if (widgetKey === "box_metrics") {
          boxMetricsPanelMode = mode;
        } else if (widgetKey === "age_kde") {
          ageKdePanelMode = mode;
        } else if (widgetKey === "cluster_filter") {
          clusterFilterPanelMode = mode;
        } else if (widgetKey === "dendrogram") {
          dendrogramPanelMode = mode;
        }
      }

      function widgetDefaultRect(widgetKey) {
        if (widgetKey === "sky") {
          return {
            left: Math.max(12, window.innerWidth - Math.min(window.innerWidth * 0.38, 560) - 12),
            top: 72,
            width: Math.min(window.innerWidth * 0.38, 560),
            height: Math.min(window.innerHeight * 0.56, 560),
          };
        }
        if (widgetKey === "box_metrics") {
          return {
            left: Math.max(12, window.innerWidth - Math.min(window.innerWidth * 0.40, 580) - 18),
            top: Math.max(72, window.innerHeight - Math.min(window.innerHeight * 0.46, 420) - 92),
            width: Math.min(window.innerWidth * 0.40, 580),
            height: Math.min(window.innerHeight * 0.46, 420),
          };
        }
        if (widgetKey === "age_kde") {
          return {
            left: Math.max(12, window.innerWidth - Math.min(window.innerWidth * 0.36, 540) - 44),
            top: 96,
            width: Math.min(window.innerWidth * 0.36, 540),
            height: Math.min(window.innerHeight * 0.42, 360),
          };
        }
        if (widgetKey === "cluster_filter") {
          return {
            left: Math.max(12, window.innerWidth - Math.min(window.innerWidth * 0.34, 480) - 72),
            top: 108,
            width: Math.min(window.innerWidth * 0.34, 480),
            height: Math.min(window.innerHeight * 0.40, 340),
          };
        }
        if (widgetKey === "dendrogram") {
          return {
            left: Math.max(12, window.innerWidth - Math.min(window.innerWidth * 0.40, 560) - 88),
            top: 104,
            width: Math.min(window.innerWidth * 0.40, 560),
            height: Math.min(window.innerHeight * 0.52, 440),
          };
        }
        return { left: 12, top: 72, width: 360, height: 260 };
      }

      function raiseWidget(widgetKey) {
        const panelEl = widgetPanelForKey(widgetKey);
        if (!panelEl) {
          return;
        }
        widgetZIndexCounter += 1;
        panelEl.style.zIndex = String(widgetZIndexCounter);
      }

      function clampWidgetPosition(left, top, width, height) {
        const margin = 6;
        return {
          left: Math.min(Math.max(margin, left), Math.max(margin, window.innerWidth - width - margin)),
          top: Math.min(Math.max(margin, top), Math.max(margin, window.innerHeight - height - margin)),
        };
      }

      function resizeWidgetRect(state, clientX, clientY) {
        const margin = 6;
        const minWidth = Math.min(220, Math.max(80, window.innerWidth - 2 * margin));
        const minHeight = Math.min(220, Math.max(80, window.innerHeight - 2 * margin));
        let left = state.startLeft;
        let right = state.startRight;
        let top = state.startTop;
        let bottom = state.startBottom;
        const dx = clientX - state.startX;
        const dy = clientY - state.startY;
        const dir = state.dir || "se";

        if (dir.includes("w")) {
          left = Math.max(margin, Math.min(state.startLeft + dx, state.startRight - minWidth));
        } else if (dir.includes("e")) {
          right = Math.min(window.innerWidth - margin, Math.max(state.startRight + dx, state.startLeft + minWidth));
        }
        if (dir.includes("n")) {
          top = Math.max(margin, Math.min(state.startTop + dy, state.startBottom - minHeight));
        } else if (dir.includes("s")) {
          bottom = Math.min(window.innerHeight - margin, Math.max(state.startBottom + dy, state.startTop + minHeight));
        }

        return {
          left,
          top,
          width: Math.max(minWidth, right - left),
          height: Math.max(minHeight, bottom - top),
        };
      }

      function storeWidgetRect(widgetKey) {
        const panelEl = widgetPanelForKey(widgetKey);
        if (!panelEl) {
          return;
        }
        const rect = panelEl.getBoundingClientRect();
        panelEl.dataset.normalLeft = String(rect.left);
        panelEl.dataset.normalTop = String(rect.top);
        panelEl.dataset.normalWidth = String(rect.width);
        panelEl.dataset.normalHeight = String(rect.height);
      }

      function restoreWidgetRect(widgetKey) {
        const panelEl = widgetPanelForKey(widgetKey);
        if (!panelEl) {
          return;
        }
        const left = Number(panelEl.dataset.normalLeft);
        const top = Number(panelEl.dataset.normalTop);
        const width = Number(panelEl.dataset.normalWidth);
        const height = Number(panelEl.dataset.normalHeight);
        const defaults = widgetDefaultRect(widgetKey);
        const next = [left, top, width, height].every(Number.isFinite)
          ? clampWidgetPosition(left, top, width, height)
          : clampWidgetPosition(defaults.left, defaults.top, defaults.width, defaults.height);
        const nextWidth = Number.isFinite(width) ? width : defaults.width;
        const nextHeight = Number.isFinite(height) ? height : defaults.height;
        panelEl.style.left = `${next.left}px`;
        panelEl.style.top = `${next.top}px`;
        panelEl.style.right = "auto";
        panelEl.style.bottom = "auto";
        panelEl.style.width = `${nextWidth}px`;
        panelEl.style.height = `${nextHeight}px`;
      }

      function applyBoxMetricsPanelMode() {
        applyWidgetMode("box_metrics");
        if (boxMetricsFullButtonEl) {
          boxMetricsFullButtonEl.setAttribute(
            "title",
            boxMetricsPanelMode === "fullscreen" ? "Restore box metrics panel" : "Maximize box metrics panel"
          );
          boxMetricsFullButtonEl.setAttribute(
            "aria-label",
            boxMetricsPanelMode === "fullscreen" ? "Restore box metrics panel" : "Maximize box metrics panel"
          );
        }
      }

      function applySkyPanelMode() {
        applyWidgetMode("sky");
        skyFullButtonEl.setAttribute("title", skyPanelMode === "fullscreen" ? "Restore sky panel" : "Maximize sky panel");
        skyFullButtonEl.setAttribute("aria-label", skyPanelMode === "fullscreen" ? "Restore sky panel" : "Maximize sky panel");
      }

      function applyAgeKdePanelMode() {
        applyWidgetMode("age_kde");
        ageKdeFullButtonEl.setAttribute("title", ageKdePanelMode === "fullscreen" ? "Restore age KDE panel" : "Maximize age KDE panel");
        ageKdeFullButtonEl.setAttribute("aria-label", ageKdePanelMode === "fullscreen" ? "Restore age KDE panel" : "Maximize age KDE panel");
      }

      function applyClusterFilterPanelMode() {
        applyWidgetMode("cluster_filter");
        clusterFilterFullButtonEl.setAttribute("title", clusterFilterPanelMode === "fullscreen" ? "Restore filter panel" : "Maximize filter panel");
        clusterFilterFullButtonEl.setAttribute("aria-label", clusterFilterPanelMode === "fullscreen" ? "Restore filter panel" : "Maximize filter panel");
      }

      function applyDendrogramPanelMode() {
        applyWidgetMode("dendrogram");
        dendrogramFullButtonEl.setAttribute("title", dendrogramPanelMode === "fullscreen" ? "Restore dendrogram panel" : "Maximize dendrogram panel");
        dendrogramFullButtonEl.setAttribute("aria-label", dendrogramPanelMode === "fullscreen" ? "Restore dendrogram panel" : "Maximize dendrogram panel");
        if (widgetModeForKey("dendrogram") === "hidden") {
          clearDendrogramSelectionState();
        }
      }

      function applyWidgetMode(widgetKey) {
        const panelEl = widgetPanelForKey(widgetKey);
        if (!panelEl) {
          return;
        }
        const mode = widgetEnabled(widgetKey) ? widgetModeForKey(widgetKey) : "hidden";
        panelEl.dataset.mode = mode;
      }

      function setWidgetMode(widgetKey, mode) {
        if (minimalModeEnabled) {
          setWidgetModeValue(widgetKey, "hidden");
          applyWidgetMode(widgetKey);
          return;
        }
        if (!widgetEnabled(widgetKey)) {
          return;
        }
        const panelEl = widgetPanelForKey(widgetKey);
        if (!panelEl) {
          return;
        }
        const currentMode = widgetModeForKey(widgetKey);
        const nextMode = ["normal", "fullscreen", "hidden"].includes(mode) ? mode : "normal";
        if (nextMode === "fullscreen" && currentMode === "normal") {
          storeWidgetRect(widgetKey);
        }
        setWidgetModeValue(widgetKey, nextMode);
        if (widgetKey === "sky") {
          applySkyPanelMode();
        } else if (widgetKey === "box_metrics") {
          applyBoxMetricsPanelMode();
        } else if (widgetKey === "age_kde") {
          applyAgeKdePanelMode();
        } else if (widgetKey === "cluster_filter") {
          applyClusterFilterPanelMode();
        } else if (widgetKey === "dendrogram") {
          applyDendrogramPanelMode();
        } else {
          applyWidgetMode(widgetKey);
        }
        if (nextMode === "normal") {
          restoreWidgetRect(widgetKey);
          raiseWidget(widgetKey);
        } else if (nextMode === "fullscreen") {
          raiseWidget(widgetKey);
          panelEl.style.left = "0px";
          panelEl.style.top = "0px";
          panelEl.style.right = "0px";
          panelEl.style.bottom = "0px";
          panelEl.style.width = "auto";
          panelEl.style.height = "auto";
        }
        resize();
        renderBoxMetricsWidget();
        if (widgetKey === "sky" && nextMode !== "hidden") {
          updateSkyPanel();
        }
        renderBoxMetricsWidget();
        renderAgeKdeWidget();
        renderClusterFilterWidget();
        renderDendrogramWidget();
        if (widgetKey === "box_metrics" || widgetKey === "dendrogram") {
          renderFrame(currentFrameIndex);
        }
      }

      function renderWidgetMenu() {
        if (!widgetSelectEl) {
          return;
        }
        if (minimalModeEnabled) {
          widgetSelectEl.innerHTML = "";
          widgetSelectEl.style.display = "none";
          widgetSelectEl.value = "";
          return;
        }
        widgetSelectEl.innerHTML = "";
        const placeholder = document.createElement("option");
        placeholder.value = "";
        placeholder.textContent = "Widgets";
        placeholder.selected = true;
        widgetSelectEl.appendChild(placeholder);

        const items = [];
        if (skySpec.enabled) {
          items.push({ key: "sky", label: "Sky View" });
        }
        if (selectionBoxSpec.enabled) {
          items.push({ key: "box_metrics", label: "Box Metrics" });
        }
        if (ageKdeSpec.enabled) {
          items.push({ key: "age_kde", label: "Age KDE" });
        }
        if (clusterFilterSpec.enabled) {
          items.push({ key: "cluster_filter", label: "Cluster Filter" });
        }
        if (dendrogramSpec.enabled) {
          items.push({ key: "dendrogram", label: "Dendrogram" });
        }

        items.forEach((item) => {
          const option = document.createElement("option");
          option.value = item.key;
          option.textContent = item.label;
          widgetSelectEl.appendChild(option);
        });
        widgetSelectEl.style.display = items.length ? "block" : "none";
        widgetSelectEl.value = "";
      }

      function onWidgetPointerStart(event) {
        const panelEl = event.target.closest(".oviz-three-widget-panel");
        if (!panelEl) {
          return;
        }
        const widgetKey = String(panelEl.dataset.widgetKey || "");
        if (widgetModeForKey(widgetKey) !== "normal") {
          return;
        }
        const resizeHandle = event.target.closest(".oviz-three-widget-resize");
        const dragHandle = event.target.closest(".oviz-three-widget-drag");
        if (!resizeHandle && !dragHandle) {
          return;
        }
        if (event.target.closest(".oviz-three-widget-window-controls")) {
          return;
        }
        const rect = panelEl.getBoundingClientRect();
        widgetPointerState = {
          widgetKey,
          panelEl,
          mode: resizeHandle ? "resize" : "drag",
          dir: resizeHandle ? (resizeHandle.dataset.dir || "se").toLowerCase() : null,
          startX: event.clientX,
          startY: event.clientY,
          startLeft: rect.left,
          startTop: rect.top,
          startWidth: rect.width,
          startHeight: rect.height,
          startRight: rect.right,
          startBottom: rect.bottom,
          handle: resizeHandle || dragHandle,
        };
        panelEl.style.left = `${rect.left}px`;
        panelEl.style.top = `${rect.top}px`;
        panelEl.style.right = "auto";
        panelEl.style.bottom = "auto";
        panelEl.style.width = `${rect.width}px`;
        panelEl.style.height = `${rect.height}px`;
        raiseWidget(widgetKey);
        controls.enabled = false;
        if (widgetPointerState.handle && typeof widgetPointerState.handle.setPointerCapture === "function" && event.pointerId !== undefined) {
          try {
            widgetPointerState.handle.setPointerCapture(event.pointerId);
          } catch (_err) {
          }
        }
        document.body.style.userSelect = "none";
        event.preventDefault();
        event.stopPropagation();
      }

      function onWidgetPointerMove(event) {
        if (!widgetPointerState) {
          return;
        }
        if (widgetPointerState.mode === "drag") {
          const left = widgetPointerState.startLeft + (event.clientX - widgetPointerState.startX);
          const top = widgetPointerState.startTop + (event.clientY - widgetPointerState.startY);
          const next = clampWidgetPosition(left, top, widgetPointerState.startWidth, widgetPointerState.startHeight);
          widgetPointerState.panelEl.style.left = `${next.left}px`;
          widgetPointerState.panelEl.style.top = `${next.top}px`;
        } else {
          const next = resizeWidgetRect(widgetPointerState, event.clientX, event.clientY);
          widgetPointerState.panelEl.style.left = `${next.left}px`;
          widgetPointerState.panelEl.style.top = `${next.top}px`;
          widgetPointerState.panelEl.style.width = `${next.width}px`;
          widgetPointerState.panelEl.style.height = `${next.height}px`;
        }
        resize();
        renderBoxMetricsWidget();
        renderAgeKdeWidget();
        event.preventDefault();
        event.stopPropagation();
      }

      function onWidgetPointerEnd(event) {
        if (!widgetPointerState) {
          return;
        }
        if (widgetPointerState.handle && typeof widgetPointerState.handle.releasePointerCapture === "function" && event.pointerId !== undefined) {
          try {
            widgetPointerState.handle.releasePointerCapture(event.pointerId);
          } catch (_err) {
          }
        }
        controls.enabled = true;
        document.body.style.userSelect = "";
        storeWidgetRect(widgetPointerState.widgetKey);
        widgetPointerState = null;
      }

      function initSkyPanel() {
        if (!skySpec.enabled) {
          applySkyPanelMode();
          return;
        }
        skyFrameEl.addEventListener("load", () => {
          lastSentSkyHoverClusterKey = null;
          postParentHoverToSkyFrame();
        });
        skyFrameEl.srcdoc = buildEmptySkySrcdoc();
        skyHideButtonEl.addEventListener("click", () => setWidgetMode("sky", "hidden"));
        skyFullButtonEl.addEventListener("click", () => {
          setWidgetMode("sky", skyPanelMode === "fullscreen" ? "normal" : "fullscreen");
        });
        applySkyPanelMode();
      }

      function initBoxMetricsPanel() {
        if (!selectionBoxSpec.enabled) {
          applyBoxMetricsPanelMode();
          return;
        }
        boxMetricsHideButtonEl.addEventListener("click", () => setWidgetMode("box_metrics", "hidden"));
        boxMetricsFullButtonEl.addEventListener("click", () => {
          setWidgetMode("box_metrics", boxMetricsPanelMode === "fullscreen" ? "normal" : "fullscreen");
        });
        boxMetricsResetButtonEl.addEventListener("click", () => {
          selectionBoxState = buildDefaultSelectionBoxState();
          syncSelectionBoxVisibilityInput(true);
          renderFrame(currentFrameIndex);
          scheduleSelectionBoxMetricsRecompute(true);
        });
        if (boxMetricsVisibleEl) {
          boxMetricsVisibleEl.addEventListener("change", () => {
            selectionBoxState.visible = Boolean(boxMetricsVisibleEl.checked);
            renderFrame(currentFrameIndex);
            renderBoxMetricsWidget();
          });
        }
        syncSelectionBoxVisibilityInput(true);
        applyBoxMetricsPanelMode();
        scheduleSelectionBoxMetricsRecompute(true);
      }

      function initAgeKdePanel() {
        if (!ageKdeSpec.enabled) {
          applyAgeKdePanelMode();
          return;
        }
        ageKdeHideButtonEl.addEventListener("click", () => setWidgetMode("age_kde", "hidden"));
        ageKdeFullButtonEl.addEventListener("click", () => {
          setWidgetMode("age_kde", ageKdePanelMode === "fullscreen" ? "normal" : "fullscreen");
        });
        if (ageKdeFilterRangeMinEl && ageKdeFilterRangeMaxEl) {
          ageKdeFilterRangeMinEl.addEventListener("input", () => {
            const parameter = ageKdeFilterParameterSpec();
            if (!parameter) {
              return;
            }
            const currentAxisRange = ageKdeAxisFilterRange();
            const minAxisValue = ageKdeSliderValueToAxisValue(ageKdeFilterRangeMinEl.value);
            const maxAxisValue = currentAxisRange ? Number(currentAxisRange.max) : ageKdeSliderValueToAxisValue(ageKdeFilterRangeMaxEl.value);
            setClusterAgeFilterFromKdeAxisRange(Math.min(minAxisValue, maxAxisValue), maxAxisValue);
          });
          ageKdeFilterRangeMaxEl.addEventListener("input", () => {
            const parameter = ageKdeFilterParameterSpec();
            if (!parameter) {
              return;
            }
            const currentAxisRange = ageKdeAxisFilterRange();
            const minAxisValue = currentAxisRange ? Number(currentAxisRange.min) : ageKdeSliderValueToAxisValue(ageKdeFilterRangeMinEl.value);
            const maxAxisValue = ageKdeSliderValueToAxisValue(ageKdeFilterRangeMaxEl.value);
            setClusterAgeFilterFromKdeAxisRange(minAxisValue, Math.max(minAxisValue, maxAxisValue));
          });
        }
        renderAgeKdeWidget();
        applyAgeKdePanelMode();
      }

      function initClusterFilterPanel() {
        if (!clusterFilterSpec.enabled) {
          applyClusterFilterPanelMode();
          return;
        }
        clusterFilterParameterEl.innerHTML = "";
        clusterFilterParameters.forEach((parameter) => {
          const option = document.createElement("option");
          option.value = String(parameter.key || "");
          option.textContent = String(parameter.label || parameter.key || "");
          clusterFilterParameterEl.appendChild(option);
        });
        if (!clusterFilterParameterKey && clusterFilterParameters.length) {
          clusterFilterParameterKey = String(clusterFilterParameters[0].key || "");
        }
        clusterFilterHideButtonEl.addEventListener("click", () => setWidgetMode("cluster_filter", "hidden"));
        clusterFilterFullButtonEl.addEventListener("click", () => {
          setWidgetMode("cluster_filter", clusterFilterPanelMode === "fullscreen" ? "normal" : "fullscreen");
        });
        clusterFilterParameterEl.addEventListener("change", () => {
          clusterFilterParameterKey = String(clusterFilterParameterEl.value || "");
          clampClusterFilterRangeForParameter(activeClusterFilterParameterSpec());
          applyClusterFilterState();
        });
        clusterFilterRangeMinEl.addEventListener("input", () => {
          const parameter = activeClusterFilterParameterSpec();
          if (!parameter) {
            return;
          }
          const rangeState = clampClusterFilterRangeForParameter(parameter);
          rangeState.min = Math.min(clusterFilterSliderValueToActual(clusterFilterRangeMinEl.value, parameter), Number(rangeState.max));
          clusterFilterRangeStateByKey[String(parameter.key)] = rangeState;
          applyClusterFilterState();
        });
        clusterFilterRangeMaxEl.addEventListener("input", () => {
          const parameter = activeClusterFilterParameterSpec();
          if (!parameter) {
            return;
          }
          const rangeState = clampClusterFilterRangeForParameter(parameter);
          rangeState.max = Math.max(clusterFilterSliderValueToActual(clusterFilterRangeMaxEl.value, parameter), Number(rangeState.min));
          clusterFilterRangeStateByKey[String(parameter.key)] = rangeState;
          applyClusterFilterState();
        });
        renderClusterFilterWidget();
        applyClusterFilterPanelMode();
      }

      function dendrogramHitRegionAtCanvasPoint(x, y) {
        let bestNode = null;
        let bestNodeDistanceSq = Infinity;
        let bestBranch = null;
        let bestBranchDistanceSq = Infinity;
        for (let index = dendrogramHitRegions.length - 1; index >= 0; index -= 1) {
          const region = dendrogramHitRegions[index];
          if (!region) {
            continue;
          }
          if (region.type === "node") {
            const dx = x - Number(region.centerX);
            const dy = y - Number(region.centerY);
            const distanceSq = dx * dx + dy * dy;
            if (distanceSq <= Math.pow(Number(region.radius) || 0.0, 2.0) && distanceSq < bestNodeDistanceSq) {
              bestNode = region;
              bestNodeDistanceSq = distanceSq;
            }
            continue;
          }
          if (!Array.isArray(region.segments)) {
            continue;
          }
          let minDistanceSq = Infinity;
          region.segments.forEach((segment) => {
            if (!Array.isArray(segment) || segment.length < 4) {
              return;
            }
            minDistanceSq = Math.min(minDistanceSq, pointSegmentDistanceSq(x, y, segment[0], segment[1], segment[2], segment[3]));
          });
          const hitRadius = Math.max(Number(region.hitRadius) || 0.0, 0.0);
          if (minDistanceSq <= hitRadius * hitRadius && minDistanceSq < bestBranchDistanceSq) {
            bestBranch = region;
            bestBranchDistanceSq = minDistanceSq;
          }
        }
        if (bestNode && bestNodeDistanceSq <= 18.0) {
          return bestNode;
        }
        return bestBranch || bestNode || null;
      }

      function onDendrogramPointerMove(event) {
        if (!dendrogramSpec.enabled || widgetModeForKey("dendrogram") === "hidden" || !dendrogramCanvasEl) {
          return;
        }
        const rect = dendrogramCanvasEl.getBoundingClientRect();
        const x = event.clientX - rect.left;
        const y = event.clientY - rect.top;
        const hitRegion = dendrogramHitRegionAtCanvasPoint(x, y);
        if (hitRegion) {
          setDendrogramHoveredSelectionKeys(hitRegion.selectionKeys, hitRegion.label, hitRegion.count, hitRegion.key);
        } else {
          clearDendrogramHoverState();
        }
      }

      function onDendrogramClick(event) {
        if (!dendrogramSpec.enabled || widgetModeForKey("dendrogram") === "hidden" || !dendrogramCanvasEl) {
          return;
        }
        const rect = dendrogramCanvasEl.getBoundingClientRect();
        const x = event.clientX - rect.left;
        const y = event.clientY - rect.top;
        const hitRegion = dendrogramHitRegionAtCanvasPoint(x, y);
        if (!hitRegion) {
          clearDendrogramPinnedState();
          return;
        }
        if (String(dendrogramPinnedRegionKey || "") === String(hitRegion.key || "")) {
          clearDendrogramPinnedState();
          return;
        }
        setDendrogramPinnedSelectionKeys(hitRegion.selectionKeys, hitRegion.label, hitRegion.count, hitRegion.key);
      }

      function initDendrogramPanel() {
        if (!dendrogramSpec.enabled) {
          applyDendrogramPanelMode();
          return;
        }
        dendrogramHideButtonEl.addEventListener("click", () => setWidgetMode("dendrogram", "hidden"));
        dendrogramFullButtonEl.addEventListener("click", () => {
          setWidgetMode("dendrogram", dendrogramPanelMode === "fullscreen" ? "normal" : "fullscreen");
        });
        dendrogramTraceEl.addEventListener("change", () => {
          dendrogramTraceKey = String(dendrogramTraceEl.value || "");
          clearDendrogramSelectionState();
          renderFrame(currentFrameIndex);
        });
        dendrogramConnectionEl.addEventListener("change", () => {
          dendrogramConnectionMode = String(dendrogramConnectionEl.value || "birth_to_older_track");
          clearDendrogramSelectionState();
          renderDendrogramWidget();
        });
        dendrogramModeEl.addEventListener("change", () => {
          dendrogramThresholdMode = String(dendrogramModeEl.value || "distance_pc");
          clearDendrogramSelectionState();
          renderDendrogramWidget();
        });
        dendrogramThresholdEl.addEventListener("change", () => {
          if (currentDendrogramThresholdMode() === "birth_age_myr") {
            dendrogramThresholdAgeMyr = Math.max(Number(dendrogramThresholdEl.value) || 0.0, 0.0);
          } else {
            dendrogramThresholdPc = Math.max(Number(dendrogramThresholdEl.value) || 0.0, 0.0);
          }
          clearDendrogramSelectionState();
          renderDendrogramWidget();
        });
        dendrogramCanvasEl.addEventListener("pointermove", onDendrogramPointerMove);
        dendrogramCanvasEl.addEventListener("pointerleave", clearDendrogramHoverState);
        dendrogramCanvasEl.addEventListener("click", onDendrogramClick);
        renderDendrogramWidget();
        applyDendrogramPanelMode();
      }
""".strip()
