from __future__ import annotations


THREEJS_WIDGET_CONTENT_RUNTIME_JS = """
      function clusterFilterParameterSpecForKey(parameterKey) {
        return clusterFilterParameters.find((parameter) => String(parameter.key) === String(parameterKey)) || null;
      }

      function activeClusterFilterParameterSpec() {
        return clusterFilterParameterSpecForKey(clusterFilterParameterKey)
          || (clusterFilterParameters.length ? clusterFilterParameters[0] : null);
      }

      function clampClusterFilterRangeForParameter(parameter) {
        if (!parameter) {
          return null;
        }
        const key = String(parameter.key || "");
        const parameterMin = Number(parameter.min);
        const parameterMax = Number(parameter.max);
        const rangeState = clusterFilterRangeStateByKey[key] || {
          min: parameterMin,
          max: parameterMax,
        };
        let minValue = Number(rangeState.min);
        let maxValue = Number(rangeState.max);
        if (!Number.isFinite(minValue)) {
          minValue = parameterMin;
        }
        if (!Number.isFinite(maxValue)) {
          maxValue = parameterMax;
        }
        minValue = clampRange(minValue, parameterMin, parameterMax);
        maxValue = clampRange(maxValue, parameterMin, parameterMax);
        if (minValue > maxValue) {
          const middle = 0.5 * (minValue + maxValue);
          minValue = middle;
          maxValue = middle;
        }
        rangeState.min = minValue;
        rangeState.max = maxValue;
        clusterFilterRangeStateByKey[key] = rangeState;
        return rangeState;
      }

      function formatClusterFilterValue(value, parameter) {
        const numericValue = Number(value);
        if (!Number.isFinite(numericValue)) {
          return "";
        }
        const unit = String((parameter && parameter.unit) || "").trim();
        const valueText = numericValue >= 1000.0
          ? formatCompactNumber(numericValue)
          : formatCompactNumber(numericValue);
        return unit ? `${valueText} ${unit}` : valueText;
      }

      function clusterFilterEntryValue(entry, parameterKey) {
        if (!entry || typeof entry !== "object") {
          return NaN;
        }
        return Number(entry[String(parameterKey || "")]);
      }

      function clusterFilterSelectionKeyForPoint(point) {
        if (point && point.motion && point.motion.key) {
          return normalizeMemberKey(point.motion.key);
        }
        if (point && point.selection) {
          return normalizedSelectionKeyFor(point.selection);
        }
        return "";
      }

      function clusterFilterPassesSelectionKey(selectionKey) {
        if (!clusterFilterSpec.enabled) {
          return true;
        }
        const key = normalizeMemberKey(selectionKey);
        if (!key) {
          return true;
        }
        const parameter = activeClusterFilterParameterSpec();
        if (!parameter) {
          return true;
        }
        const rangeState = clampClusterFilterRangeForParameter(parameter);
        const entry = clusterFilterEntryByKey.get(key);
        if (!entry || !rangeState) {
          return true;
        }
        const value = clusterFilterEntryValue(entry, parameter.key);
        if (!Number.isFinite(value)) {
          return true;
        }
        return value >= Number(rangeState.min) && value <= Number(rangeState.max);
      }

      function clusterFilterPassesSelectionKeyExcludingParameter(selectionKey, excludedParameterKey) {
        if (!clusterFilterSpec.enabled) {
          return true;
        }
        const key = normalizeMemberKey(selectionKey);
        if (!key) {
          return true;
        }
        const parameter = activeClusterFilterParameterSpec();
        if (!parameter) {
          return true;
        }
        if (excludedParameterKey && String(parameter.key) === String(excludedParameterKey)) {
          return true;
        }
        const rangeState = clampClusterFilterRangeForParameter(parameter);
        const entry = clusterFilterEntryByKey.get(key);
        if (!entry || !rangeState) {
          return true;
        }
        const value = clusterFilterEntryValue(entry, parameter.key);
        if (!Number.isFinite(value)) {
          return true;
        }
        return value >= Number(rangeState.min) && value <= Number(rangeState.max);
      }

      function clusterFilterPassesSelection(selection) {
        return clusterFilterPassesSelectionKey(normalizedSelectionKeyFor(selection));
      }

      function clusterFilterPassesPoint(point) {
        return clusterFilterPassesSelectionKey(clusterFilterSelectionKeyForPoint(point));
      }

      function pruneSelectionsToActiveClusterFilter() {
        currentSelections = uniqueSelections(currentSelections.filter((selection) => clusterFilterPassesSelection(selection)));
        if (currentSelection && !clusterFilterPassesSelection(currentSelection)) {
          currentSelection = null;
        }
        const nextSelectedKeys = new Set();
        currentSelections.forEach((selection) => {
          const key = normalizedSelectionKeyFor(selection);
          if (key) {
            nextSelectedKeys.add(key);
          }
        });
        if (currentSelection) {
          const focusKey = normalizedSelectionKeyFor(currentSelection);
          if (focusKey) {
            nextSelectedKeys.add(focusKey);
          }
        }
        selectedClusterKeys = nextSelectedKeys;
        currentSelectionMode = currentSelection ? "click" : ((currentSelections.length || hasActiveLassoSelectionMask()) ? "lasso" : "none");
        if (!clusterFilterPassesSelectionKey(localHoveredClusterKey)) {
          setLocalHoveredClusterKey("");
        }
        if (!clusterFilterPassesSelectionKey(skyHoveredClusterKey)) {
          skyHoveredClusterKey = "";
          lastSentSkyHoverClusterKey = null;
        }
      }

      function ageKdeGridValues() {
        const traces = Array.isArray(ageKdeSpec.traces) ? ageKdeSpec.traces : [];
        if (traces.length && Array.isArray(traces[0].x) && traces[0].x.length >= 2) {
          return traces[0].x.map((value) => Number(value)).filter(Number.isFinite);
        }
        const xRange = Array.isArray(ageKdeSpec.x_range) ? ageKdeSpec.x_range : [-1.0, 0.0];
        const xMin = Number(xRange[0]);
        const xMax = Number(xRange[1]);
        const count = 300;
        return Array.from({ length: count }, (_, index) => {
          const frac = count <= 1 ? 0.0 : index / (count - 1);
          return xMin + frac * (xMax - xMin);
        });
      }

      function gaussianKde(values, xGrid, bandwidth) {
        const finiteValues = (Array.isArray(values) ? values : []).map((value) => Number(value)).filter(Number.isFinite);
        if (!finiteValues.length) {
          return xGrid.map(() => 0.0);
        }
        const bw = Math.max(Number(bandwidth) || 1.0, 1e-6);
        const norm = finiteValues.length * bw * Math.sqrt(2.0 * Math.PI);
        return xGrid.map((xValue) => {
          let sum = 0.0;
          for (const value of finiteValues) {
            const u = (Number(xValue) - value) / bw;
            sum += Math.exp(-0.5 * u * u);
          }
          return sum / norm;
        });
      }

      function ageKdeFilterParameterSpec() {
        return clusterFilterParameterSpecForKey("age_now_myr");
      }

      function ageKdeAxisRange() {
        const xRange = Array.isArray(ageKdeSpec.x_range) ? ageKdeSpec.x_range : [-1.0, 0.0];
        let minValue = Number(xRange[0]);
        let maxValue = Number(xRange[1]);
        if (!Number.isFinite(minValue)) {
          minValue = -1.0;
        }
        if (!Number.isFinite(maxValue)) {
          maxValue = 0.0;
        }
        if (minValue > maxValue) {
          const swap = minValue;
          minValue = maxValue;
          maxValue = swap;
        }
        return { min: minValue, max: maxValue };
      }

      function ageNowToKdeAxisValue(ageNowMyr) {
        const value = Number(ageNowMyr);
        if (!Number.isFinite(value)) {
          return NaN;
        }
        return -Math.abs(value);
      }

      function kdeAxisValueToAgeNow(axisValue) {
        const value = Number(axisValue);
        if (!Number.isFinite(value)) {
          return NaN;
        }
        return Math.abs(value);
      }

      function ageKdeAxisValueToSlider(axisValue) {
        const axisRange = ageKdeAxisRange();
        const denom = Math.max(axisRange.max - axisRange.min, 1e-9);
        return Math.round(1000.0 * clampRange((Number(axisValue) - axisRange.min) / denom, 0.0, 1.0));
      }

      function ageKdeSliderValueToAxisValue(sliderValue) {
        const axisRange = ageKdeAxisRange();
        const normalized = clampRange(Number(sliderValue) / 1000.0, 0.0, 1.0);
        return axisRange.min + normalized * (axisRange.max - axisRange.min);
      }

      function ageKdeAxisFilterRange() {
        const parameter = ageKdeFilterParameterSpec();
        if (!parameter) {
          return null;
        }
        const rangeState = clampClusterFilterRangeForParameter(parameter);
        const axisRange = ageKdeAxisRange();
        let minAxis = clampRange(ageNowToKdeAxisValue(rangeState.max), axisRange.min, axisRange.max);
        let maxAxis = clampRange(ageNowToKdeAxisValue(rangeState.min), axisRange.min, axisRange.max);
        if (minAxis > maxAxis) {
          const swap = minAxis;
          minAxis = maxAxis;
          maxAxis = swap;
        }
        return { min: minAxis, max: maxAxis };
      }

      function formatAgeKdeAxisValue(axisValue) {
        const value = Number(axisValue);
        if (!Number.isFinite(value)) {
          return "";
        }
        return `${formatCompactNumber(value)} Myr`;
      }

      function setClusterAgeFilterFromKdeAxisRange(minAxisValue, maxAxisValue) {
        const parameter = ageKdeFilterParameterSpec();
        if (!parameter) {
          return;
        }
        const axisRange = ageKdeAxisRange();
        const clampedMinAxis = clampRange(Math.min(Number(minAxisValue), Number(maxAxisValue)), axisRange.min, axisRange.max);
        const clampedMaxAxis = clampRange(Math.max(Number(minAxisValue), Number(maxAxisValue)), axisRange.min, axisRange.max);
        const youngerAge = clampRange(
          kdeAxisValueToAgeNow(clampedMaxAxis),
          Number(parameter.min),
          Number(parameter.max),
        );
        const olderAge = clampRange(
          kdeAxisValueToAgeNow(clampedMinAxis),
          Number(parameter.min),
          Number(parameter.max),
        );
        clusterFilterParameterKey = String(parameter.key || "");
        clusterFilterRangeStateByKey[String(parameter.key)] = {
          min: Math.min(youngerAge, olderAge),
          max: Math.max(youngerAge, olderAge),
        };
        applyClusterFilterState();
      }

      function ageKdeBaseEntries() {
        const ageParameter = ageKdeFilterParameterSpec();
        if (!ageParameter) {
          return [];
        }
        const selectedKeys = selectedClusterKeys.size ? selectedClusterKeys : null;
        return clusterFilterEntries.filter((entry) => {
          if (!entry || typeof entry !== "object") {
            return false;
          }
          const selectionKey = normalizeMemberKey(entry.selection_key || "");
          if (!selectionKey) {
            return false;
          }
          if (!clusterFilterEntryVisibleInScene(entry)) {
            return false;
          }
          if (selectedKeys && selectedKeys.size && !selectedKeys.has(selectionKey)) {
            return false;
          }
          if (!clusterFilterPassesSelectionKeyExcludingParameter(selectionKey, ageParameter.key)) {
            return false;
          }
          return Number.isFinite(clusterFilterEntryValue(entry, ageParameter.key));
        });
      }

      function filteredAgeKdeSeries() {
        const defaults = groupDefaults(currentGroup);
        const selectedKeys = selectedClusterKeys.size ? selectedClusterKeys : null;
        const traceMetaByKey = new Map();
        const traceMetaByName = new Map();
        (ageKdeSpec.traces || []).forEach((traceSpec) => {
          if (traceSpec.trace_key) {
            traceMetaByKey.set(String(traceSpec.trace_key), traceSpec);
          }
          if (traceSpec.trace_name) {
            traceMetaByName.set(String(traceSpec.trace_name), traceSpec);
          }
        });

        const groupedAges = new Map();
        (ageKdeSpec.cluster_points || []).forEach((point) => {
          const traceKey = point && point.trace_key ? String(point.trace_key) : null;
          const traceName = point && point.trace_name ? String(point.trace_name) : "";
          const traceLookupKey = traceKey || traceName;
          if (!traceLookupKey) {
            return;
          }
          if (traceKey) {
            const mode = defaults[traceKey];
            if (mode !== true || legendState[traceKey] === false) {
              return;
            }
          }
          if (selectedKeys && selectedKeys.size) {
            const selectionLike = {
              cluster_name: point.cluster_name,
              trace_name: traceName,
            };
            const pointKey = normalizedSelectionKeyFor(selectionLike);
            if (!pointKey || !selectedKeys.has(pointKey)) {
              return;
            }
          }
          const pointSelectionKey = normalizedSelectionKeyFor({
            cluster_name: point.cluster_name,
            trace_name: traceName,
          });
          if (!clusterFilterPassesSelectionKey(pointSelectionKey)) {
            return;
          }
          const ageNow = Number(point.age_now_myr);
          if (!Number.isFinite(ageNow)) {
            return;
          }
          if (!groupedAges.has(traceLookupKey)) {
            groupedAges.set(traceLookupKey, []);
          }
          groupedAges.get(traceLookupKey).push(-Math.abs(ageNow));
        });

        const xGrid = ageKdeGridValues();
        const series = [];
        groupedAges.forEach((lookbackValues, traceLookupKey) => {
          if (!lookbackValues.length) {
            return;
          }
          const meta = traceMetaByKey.get(String(traceLookupKey)) || traceMetaByName.get(String(traceLookupKey)) || {};
          const traceKey = meta.trace_key ? String(meta.trace_key) : (String(traceLookupKey).startsWith("trace-") ? String(traceLookupKey) : null);
          const styleState = traceKey ? traceStyleStateForKey(traceKey) : null;
          const densityRaw = gaussianKde(lookbackValues, xGrid, ageKdeSpec.bandwidth_myr);
          const densityMax = Math.max(...densityRaw, 0.0);
          const density = densityMax > 0.0 ? densityRaw.map((value) => value / densityMax) : densityRaw;
          series.push({
            traceName: String(meta.trace_name || traceLookupKey),
            traceKey,
            color: styleState && styleState.color ? styleState.color : (meta.color || axisSpec.x?.linecolor || theme.axis_color || "#808080"),
            opacity: clamp01((styleState && Number.isFinite(styleState.opacity) ? styleState.opacity : 1.0) * Number(meta.opacity ?? 1.0)),
            x: xGrid,
            y: density,
            count: lookbackValues.length,
          });
        });

        series.sort((a, b) => String(a.traceName).localeCompare(String(b.traceName)));
        return series;
      }

      function renderAgeKdeWidget() {
        if (!ageKdeSpec.enabled || !ageKdeCanvasEl || widgetModeForKey("age_kde") === "hidden") {
          return;
        }

        const rect = ageKdeCanvasEl.getBoundingClientRect();
        const cssWidth = Math.max(1, Math.round(rect.width));
        const cssHeight = Math.max(1, Math.round(rect.height));
        const dpr = Math.max(window.devicePixelRatio || 1, 1);
        if (ageKdeCanvasEl.width !== Math.round(cssWidth * dpr) || ageKdeCanvasEl.height !== Math.round(cssHeight * dpr)) {
          ageKdeCanvasEl.width = Math.round(cssWidth * dpr);
          ageKdeCanvasEl.height = Math.round(cssHeight * dpr);
        }

        const ctx = ageKdeCanvasEl.getContext("2d");
        if (!ctx) {
          return;
        }
        ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
        ctx.clearRect(0, 0, cssWidth, cssHeight);
        ctx.fillStyle = theme.scene_bgcolor || theme.paper_bgcolor || "#000000";
        ctx.fillRect(0, 0, cssWidth, cssHeight);

        const margin = { left: 44, right: 16, top: 18, bottom: 28 };
        const plotWidth = Math.max(40, cssWidth - margin.left - margin.right);
        const plotHeight = Math.max(40, cssHeight - margin.top - margin.bottom);
        const xRange = Array.isArray(ageKdeSpec.x_range) ? ageKdeSpec.x_range : [-1.0, 0.0];
        const xMin = Number(xRange[0]);
        const xMax = Number(xRange[1]);
        const axisColor = String(ageKdeSpec.axis_color || theme.axis_color || "#808080");
        const filteredSeries = filteredAgeKdeSeries();
        const ageFilterParameter = ageKdeFilterParameterSpec();
        const ageRangeState = ageFilterParameter ? clampClusterFilterRangeForParameter(ageFilterParameter) : null;
        const ageAxisRangeState = ageFilterParameter ? ageKdeAxisFilterRange() : null;
        const ageEntries = ageKdeBaseEntries();
        const yMax = Math.max(
          ...filteredSeries.map((series) => Math.max(...series.y, 0.0)),
          Number((Array.isArray(ageKdeSpec.y_range) ? ageKdeSpec.y_range[1] : 1.0)) || 1.0,
          1e-6,
        );

        function xToPx(value) {
          const denom = Math.max(xMax - xMin, 1e-6);
          return margin.left + ((Number(value) - xMin) / denom) * plotWidth;
        }

        function yToPx(value) {
          return margin.top + plotHeight - (Math.max(Number(value), 0.0) / yMax) * plotHeight;
        }

        ctx.strokeStyle = axisColor;
        ctx.lineWidth = 1.5;
        ctx.beginPath();
        ctx.moveTo(margin.left, margin.top);
        ctx.lineTo(margin.left, margin.top + plotHeight);
        ctx.lineTo(margin.left + plotWidth, margin.top + plotHeight);
        ctx.stroke();

        ctx.fillStyle = axisColor;
        ctx.font = "11px Helvetica, Arial, sans-serif";
        ctx.textBaseline = "middle";
        ctx.textAlign = "right";
        [0.0, 0.5, 1.0].forEach((frac) => {
          const yValue = frac * yMax;
          const yPx = yToPx(yValue);
          ctx.globalAlpha = 0.16;
          ctx.beginPath();
          ctx.moveTo(margin.left, yPx);
          ctx.lineTo(margin.left + plotWidth, yPx);
          ctx.stroke();
          ctx.globalAlpha = 1.0;
          ctx.fillText(formatCompactNumber(yValue), margin.left - 8, yPx);
        });

        ctx.textAlign = "center";
        ctx.textBaseline = "top";
        [xMin, 0.5 * (xMin + xMax), xMax].forEach((xValue) => {
          const xPx = xToPx(xValue);
          ctx.fillText(formatCompactNumber(xValue), xPx, margin.top + plotHeight + 8);
        });

        const visibleTraceNames = [];
        filteredSeries.forEach((traceSpec) => {
          const lineColor = traceSpec.color || axisColor;
          const lineOpacity = clamp01(traceSpec.opacity);
          const xValues = Array.isArray(traceSpec.x) ? traceSpec.x : [];
          const yValues = Array.isArray(traceSpec.y) ? traceSpec.y : [];
          if (xValues.length < 2 || yValues.length !== xValues.length) {
            return;
          }
          visibleTraceNames.push(String(traceSpec.traceName || traceSpec.traceKey || "trace"));

          ctx.beginPath();
          ctx.moveTo(xToPx(xValues[0]), margin.top + plotHeight);
          for (let i = 0; i < xValues.length; i += 1) {
            ctx.lineTo(xToPx(xValues[i]), yToPx(yValues[i]));
          }
          ctx.lineTo(xToPx(xValues[xValues.length - 1]), margin.top + plotHeight);
          ctx.closePath();
          ctx.fillStyle = cssColorWithAlpha(lineColor, Math.min(0.28, 0.14 + 0.18 * lineOpacity), axisColor);
          ctx.fill();

          ctx.beginPath();
          for (let i = 0; i < xValues.length; i += 1) {
            const xPx = xToPx(xValues[i]);
            const yPx = yToPx(yValues[i]);
            if (i === 0) {
              ctx.moveTo(xPx, yPx);
            } else {
              ctx.lineTo(xPx, yPx);
            }
          }
          ctx.strokeStyle = cssColorWithAlpha(lineColor, lineOpacity, axisColor);
          ctx.lineWidth = 2.0;
          ctx.stroke();
        });

        const frame = currentFrame();
        const timeValue = frame ? Number(frame.time) : 0.0;
        const markerX = xToPx(Math.min(Math.max(timeValue, xMin), xMax));
        ctx.save();
        ctx.setLineDash([6, 6]);
        ctx.strokeStyle = axisColor;
        ctx.lineWidth = 1.5;
        ctx.beginPath();
        ctx.moveTo(markerX, margin.top);
        ctx.lineTo(markerX, margin.top + plotHeight);
        ctx.stroke();
        ctx.restore();

        if (ageFilterParameter && ageRangeState && ageAxisRangeState) {
          if (ageKdeFilterRangeMinEl) {
            ageKdeFilterRangeMinEl.value = String(ageKdeAxisValueToSlider(ageAxisRangeState.min));
          }
          if (ageKdeFilterRangeMaxEl) {
            ageKdeFilterRangeMaxEl.value = String(ageKdeAxisValueToSlider(ageAxisRangeState.max));
          }
          if (ageKdeFilterRangeReadoutMinEl) {
            ageKdeFilterRangeReadoutMinEl.textContent = formatAgeKdeAxisValue(ageAxisRangeState.min);
          }
          if (ageKdeFilterRangeReadoutMaxEl) {
            ageKdeFilterRangeReadoutMaxEl.textContent = formatAgeKdeAxisValue(ageAxisRangeState.max);
          }
        }

      }

      function renderBoxMetricsWidget() {
        syncSelectionBoxVisibilityInput();
        if (boxMetricsSummaryEl) {
          boxMetricsSummaryEl.textContent = selectionBoxSummaryText(selectionBoxMetricsCache, selectionBoxMetricsPending);
        }
        if (!selectionBoxSpec.enabled || !boxMetricsCanvasEl || widgetModeForKey("box_metrics") === "hidden") {
          return;
        }

        const rect = boxMetricsCanvasEl.getBoundingClientRect();
        const cssWidth = Math.max(1, Math.round(rect.width));
        const cssHeight = Math.max(1, Math.round(rect.height));
        const dpr = Math.max(window.devicePixelRatio || 1, 1);
        if (boxMetricsCanvasEl.width !== Math.round(cssWidth * dpr) || boxMetricsCanvasEl.height !== Math.round(cssHeight * dpr)) {
          boxMetricsCanvasEl.width = Math.round(cssWidth * dpr);
          boxMetricsCanvasEl.height = Math.round(cssHeight * dpr);
        }
        const ctx = boxMetricsCanvasEl.getContext("2d");
        if (!ctx) {
          return;
        }
        ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
        ctx.clearRect(0, 0, cssWidth, cssHeight);
        ctx.fillStyle = theme.scene_bgcolor || theme.paper_bgcolor || "#000000";
        ctx.fillRect(0, 0, cssWidth, cssHeight);

        const axisColor = String(theme.axis_color || "#808080");
        const topLabelColor = String(theme.text_color || "#d0d0d0");
        const frame = currentFrame();
        const currentLookback = Math.abs(Number(frame && frame.time) || 0.0);
        const margin = { left: 48, right: 14, top: 18, bottom: 28 };
        const gap = 28;
        const plotWidth = Math.max(80, cssWidth - margin.left - margin.right);
        const totalPlotHeight = Math.max(80, cssHeight - margin.top - margin.bottom);
        const panelHeight = Math.max(56, (totalPlotHeight - gap) / 2.0);
        const rateRect = { x: margin.left, y: margin.top, width: plotWidth, height: panelHeight };
        const corrRect = { x: margin.left, y: margin.top + panelHeight + gap, width: plotWidth, height: panelHeight };
        const lookbackMax = Math.max(Number(selectionBoxMetricsSpec.lookback_max_myr) || 0.0, 1.0);

        function xToPx(value) {
          return rateRect.x + (clampRange(Number(value), 0.0, lookbackMax) / lookbackMax) * rateRect.width;
        }

        function yToPx(rectDef, value, yMin, yMax) {
          const denom = Math.max(yMax - yMin, 1e-6);
          return rectDef.y + rectDef.height - ((Number(value) - yMin) / denom) * rectDef.height;
        }

        function drawPanelAxes(rectDef, yMin, yMax, title, showXLabels, includeZeroLine = false) {
          ctx.strokeStyle = axisColor;
          ctx.lineWidth = 1.3;
          ctx.beginPath();
          ctx.moveTo(rectDef.x, rectDef.y);
          ctx.lineTo(rectDef.x, rectDef.y + rectDef.height);
          ctx.lineTo(rectDef.x + rectDef.width, rectDef.y + rectDef.height);
          ctx.stroke();
          if (includeZeroLine && yMin < 0.0 && yMax > 0.0) {
            const zeroY = yToPx(rectDef, 0.0, yMin, yMax);
            ctx.save();
            ctx.setLineDash([4, 4]);
            ctx.globalAlpha = 0.45;
            ctx.beginPath();
            ctx.moveTo(rectDef.x, zeroY);
            ctx.lineTo(rectDef.x + rectDef.width, zeroY);
            ctx.stroke();
            ctx.restore();
          }
          ctx.fillStyle = topLabelColor;
          ctx.font = "11px Helvetica, Arial, sans-serif";
          ctx.textAlign = "left";
          ctx.textBaseline = "bottom";
          ctx.fillText(title, rectDef.x, rectDef.y - 4);

          ctx.textAlign = "right";
          ctx.textBaseline = "middle";
          [0.0, 0.5, 1.0].forEach((fraction) => {
            const value = yMin + fraction * (yMax - yMin);
            const yPx = yToPx(rectDef, value, yMin, yMax);
            ctx.globalAlpha = 0.14;
            ctx.beginPath();
            ctx.moveTo(rectDef.x, yPx);
            ctx.lineTo(rectDef.x + rectDef.width, yPx);
            ctx.stroke();
            ctx.globalAlpha = 1.0;
            ctx.fillText(formatCompactNumber(value), rectDef.x - 8, yPx);
          });

          if (showXLabels) {
            ctx.textAlign = "center";
            ctx.textBaseline = "top";
            [0.0, 0.5 * lookbackMax, lookbackMax].forEach((value) => {
              const xPx = xToPx(value);
              ctx.fillText(formatCompactNumber(value), xPx, rectDef.y + rectDef.height + 8);
            });
            ctx.textAlign = "right";
            ctx.fillText("Lookback Time (Myr)", rectDef.x + rectDef.width, rectDef.y + rectDef.height + 20);
          }

          const markerX = xToPx(currentLookback);
          ctx.save();
          ctx.setLineDash([6, 6]);
          ctx.strokeStyle = axisColor;
          ctx.lineWidth = 1.3;
          ctx.beginPath();
          ctx.moveTo(markerX, rectDef.y);
          ctx.lineTo(markerX, rectDef.y + rectDef.height);
          ctx.stroke();
          ctx.restore();
        }

        if (!selectionBoxMetricsCache) {
          ctx.fillStyle = topLabelColor;
          ctx.font = "12px Menlo, Monaco, Consolas, monospace";
          ctx.textAlign = "center";
          ctx.textBaseline = "middle";
          ctx.fillText(selectionBoxMetricsPending ? "Updating box metrics..." : "No box metrics available.", cssWidth / 2, cssHeight / 2);
          return;
        }

        const rate = selectionBoxMetricsCache.rate || { centers: [], mean: [], lo: [], hi: [] };
        const clustering = selectionBoxMetricsCache.clustering || { centers: [], median: [], lo: [], hi: [] };
        const finiteRate = rate.mean.filter((value) => Number.isFinite(value));
        const finiteRateHi = rate.hi.filter((value) => Number.isFinite(value));
        const rateMax = Math.max(0.5, ...finiteRate, ...finiteRateHi, 0.0);
        drawPanelAxes(rateRect, 0.0, rateMax * 1.08, "Mean ccSN Volume Density in Box (Myr^-1 kpc^-3)", false, false);

        if (Array.isArray(rate.centers) && rate.centers.length && Array.isArray(rate.mean) && rate.mean.length === rate.centers.length) {
          const fillStarted = [];
          ctx.beginPath();
          rate.centers.forEach((centerValue, index) => {
            const loValue = Number(rate.lo[index]);
            const hiValue = Number(rate.hi[index]);
            if (!Number.isFinite(loValue) || !Number.isFinite(hiValue)) {
              return;
            }
            const xPx = xToPx(centerValue);
            const yPx = yToPx(rateRect, hiValue, 0.0, rateMax * 1.08);
            if (!fillStarted.length) {
              ctx.moveTo(xPx, yPx);
              fillStarted.push(index);
            } else {
              ctx.lineTo(xPx, yPx);
            }
          });
          for (let index = rate.centers.length - 1; index >= 0; index -= 1) {
            const loValue = Number(rate.lo[index]);
            const hiValue = Number(rate.hi[index]);
            if (!Number.isFinite(loValue) || !Number.isFinite(hiValue)) {
              continue;
            }
            const xPx = xToPx(rate.centers[index]);
            const yPx = yToPx(rateRect, loValue, 0.0, rateMax * 1.08);
            ctx.lineTo(xPx, yPx);
          }
          if (fillStarted.length) {
            ctx.closePath();
            ctx.fillStyle = cssColorWithAlpha("#73a3ff", 0.20, "#73a3ff");
            ctx.fill();
          }
          ctx.beginPath();
          rate.centers.forEach((centerValue, index) => {
            const meanValue = Number(rate.mean[index]);
            if (!Number.isFinite(meanValue)) {
              return;
            }
            const xPx = xToPx(centerValue);
            const yPx = yToPx(rateRect, meanValue, 0.0, rateMax * 1.08);
            if (index === 0) {
              ctx.moveTo(xPx, yPx);
            } else {
              ctx.lineTo(xPx, yPx);
            }
          });
          ctx.strokeStyle = "#73a3ff";
          ctx.lineWidth = 2.0;
          ctx.stroke();
        }

        const finiteClustering = clustering.median.filter((value) => Number.isFinite(value));
        const finiteClusteringHi = clustering.hi.filter((value) => Number.isFinite(value));
        const clusteringMax = Math.max(1.25, ...finiteClustering, ...finiteClusteringHi, 1.0);
        drawPanelAxes(corrRect, 0.0, clusteringMax * 1.08, "NN Clustering Enhancement in Box", false, true);

        const unityY = yToPx(corrRect, 1.0, 0.0, clusteringMax * 1.08);
        ctx.save();
        ctx.setLineDash([5, 5]);
        ctx.strokeStyle = cssColorWithAlpha(axisColor, 0.55, axisColor);
        ctx.lineWidth = 1.2;
        ctx.beginPath();
        ctx.moveTo(corrRect.x, unityY);
        ctx.lineTo(corrRect.x + corrRect.width, unityY);
        ctx.stroke();
        ctx.restore();

        if (
          Array.isArray(clustering.centers)
          && clustering.centers.length
          && Array.isArray(clustering.lo)
          && Array.isArray(clustering.hi)
          && clustering.lo.length === clustering.centers.length
          && clustering.hi.length === clustering.centers.length
        ) {
          let started = false;
          ctx.beginPath();
          clustering.centers.forEach((centerValue, index) => {
            const loValue = Number(clustering.lo[index]);
            const hiValue = Number(clustering.hi[index]);
            if (!Number.isFinite(loValue) || !Number.isFinite(hiValue)) {
              return;
            }
            const xPx = xToPx(centerValue);
            const yPx = yToPx(corrRect, hiValue, 0.0, clusteringMax * 1.08);
            if (!started) {
              ctx.moveTo(xPx, yPx);
              started = true;
            } else {
              ctx.lineTo(xPx, yPx);
            }
          });
          for (let index = clustering.centers.length - 1; index >= 0; index -= 1) {
            const loValue = Number(clustering.lo[index]);
            const hiValue = Number(clustering.hi[index]);
            if (!Number.isFinite(loValue) || !Number.isFinite(hiValue)) {
              continue;
            }
            const xPx = xToPx(clustering.centers[index]);
            const yPx = yToPx(corrRect, loValue, 0.0, clusteringMax * 1.08);
            ctx.lineTo(xPx, yPx);
          }
          if (started) {
            ctx.closePath();
            ctx.fillStyle = cssColorWithAlpha("#9FB3C8", 0.24, "#9FB3C8");
            ctx.fill();
          }
        }

        if (
          Array.isArray(clustering.centers)
          && clustering.centers.length
          && Array.isArray(clustering.median)
          && clustering.median.length === clustering.centers.length
        ) {
          ctx.beginPath();
          let moved = false;
          clustering.centers.forEach((centerValue, index) => {
            const metricValue = Number(clustering.median[index]);
            if (!Number.isFinite(metricValue)) {
              return;
            }
            const xPx = xToPx(centerValue);
            const yPx = yToPx(corrRect, metricValue, 0.0, clusteringMax * 1.08);
            if (!moved) {
              ctx.moveTo(xPx, yPx);
              moved = true;
            } else {
              ctx.lineTo(xPx, yPx);
            }
          });
          if (moved) {
            ctx.strokeStyle = "#1F2A44";
            ctx.lineWidth = 2.0;
            ctx.stroke();
          }
        }
      }

      function clusterFilterSliderValueToActual(sliderValue, parameter) {
        const spec = parameter || activeClusterFilterParameterSpec();
        if (!spec) {
          return 0.0;
        }
        const minValue = Number(spec.min);
        const maxValue = Number(spec.max);
        if (!Number.isFinite(minValue) || !Number.isFinite(maxValue) || Math.abs(maxValue - minValue) <= 1e-12) {
          return minValue;
        }
        const normalized = clampRange(Number(sliderValue) / 1000.0, 0.0, 1.0);
        return minValue + normalized * (maxValue - minValue);
      }

      function clusterFilterActualValueToSlider(actualValue, parameter) {
        const spec = parameter || activeClusterFilterParameterSpec();
        if (!spec) {
          return 0;
        }
        const minValue = Number(spec.min);
        const maxValue = Number(spec.max);
        if (!Number.isFinite(minValue) || !Number.isFinite(maxValue) || Math.abs(maxValue - minValue) <= 1e-12) {
          return 0;
        }
        return Math.round(1000.0 * clampRange((Number(actualValue) - minValue) / (maxValue - minValue), 0.0, 1.0));
      }

      function clusterFilterEntryVisibleInScene(entry) {
        if (!entry || typeof entry !== "object") {
          return false;
        }
        const traceKey = entry.trace_key ? String(entry.trace_key) : "";
        if (!traceKey) {
          return true;
        }
        const defaults = groupDefaults(currentGroup);
        const mode = defaults[traceKey];
        if (mode === false || mode === undefined) {
          return false;
        }
        return legendState[traceKey] !== false;
      }

      function applyClusterFilterState() {
        clampClusterFilterRangeForParameter(activeClusterFilterParameterSpec());
        pruneSelectionsToActiveClusterFilter();
        updateSelectionUI();
        updateSkyPanel();
        renderFrame(currentFrameIndex);
      }

      function renderClusterFilterWidget() {
        if (!clusterFilterSpec.enabled || !clusterFilterCanvasEl || widgetModeForKey("cluster_filter") === "hidden") {
          return;
        }

        const parameter = activeClusterFilterParameterSpec();
        if (!parameter) {
          return;
        }
        const rangeState = clampClusterFilterRangeForParameter(parameter);
        clusterFilterParameterKey = String(parameter.key);
        if (clusterFilterParameterEl) {
          clusterFilterParameterEl.value = clusterFilterParameterKey;
        }
        if (clusterFilterRangeMinEl) {
          clusterFilterRangeMinEl.value = String(clusterFilterActualValueToSlider(rangeState.min, parameter));
        }
        if (clusterFilterRangeMaxEl) {
          clusterFilterRangeMaxEl.value = String(clusterFilterActualValueToSlider(rangeState.max, parameter));
        }
        if (clusterFilterRangeReadoutMinEl) {
          clusterFilterRangeReadoutMinEl.textContent = formatClusterFilterValue(rangeState.min, parameter);
        }
        if (clusterFilterRangeReadoutMaxEl) {
          clusterFilterRangeReadoutMaxEl.textContent = formatClusterFilterValue(rangeState.max, parameter);
        }

        const finiteEntries = clusterFilterEntries.filter((entry) => Number.isFinite(clusterFilterEntryValue(entry, parameter.key)));
        const rect = clusterFilterCanvasEl.getBoundingClientRect();
        const cssWidth = Math.max(1, Math.round(rect.width));
        const cssHeight = Math.max(1, Math.round(rect.height));
        const dpr = Math.max(window.devicePixelRatio || 1, 1);
        if (clusterFilterCanvasEl.width !== Math.round(cssWidth * dpr) || clusterFilterCanvasEl.height !== Math.round(cssHeight * dpr)) {
          clusterFilterCanvasEl.width = Math.round(cssWidth * dpr);
          clusterFilterCanvasEl.height = Math.round(cssHeight * dpr);
        }

        const ctx = clusterFilterCanvasEl.getContext("2d");
        if (!ctx) {
          return;
        }
        ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
        ctx.clearRect(0, 0, cssWidth, cssHeight);
        ctx.fillStyle = theme.scene_bgcolor || theme.paper_bgcolor || "#000000";
        ctx.fillRect(0, 0, cssWidth, cssHeight);

        const margin = { left: 34, right: 12, top: 10, bottom: 22 };
        const plotWidth = Math.max(40, cssWidth - margin.left - margin.right);
        const plotHeight = Math.max(40, cssHeight - margin.top - margin.bottom);
        const xMin = Number(parameter.min);
        const xMax = Number(parameter.max);
        const axisColor = String(theme.axis_color || "#808080");
        const bins = Math.max(12, Math.min(28, Math.round(plotWidth / 18)));
        const counts = new Array(bins).fill(0);
        const highlightedCounts = new Array(bins).fill(0);
        const denom = Math.max(xMax - xMin, 1e-9);

        finiteEntries.forEach((entry) => {
          const value = clusterFilterEntryValue(entry, parameter.key);
          const frac = clampRange((value - xMin) / denom, 0.0, 0.999999);
          const binIndex = Math.max(0, Math.min(bins - 1, Math.floor(frac * bins)));
          counts[binIndex] += 1;
          if (value >= rangeState.min && value <= rangeState.max) {
            highlightedCounts[binIndex] += 1;
          }
        });

        const yMax = Math.max(...counts, 1);
        const binWidth = plotWidth / bins;
        ctx.strokeStyle = axisColor;
        ctx.lineWidth = 1.0;
        ctx.beginPath();
        ctx.moveTo(margin.left, margin.top);
        ctx.lineTo(margin.left, margin.top + plotHeight);
        ctx.lineTo(margin.left + plotWidth, margin.top + plotHeight);
        ctx.stroke();

        for (let i = 0; i < bins; i += 1) {
          const totalCount = counts[i];
          if (!totalCount) {
            continue;
          }
          const x = margin.left + i * binWidth + 1;
          const fullHeight = (totalCount / yMax) * plotHeight;
          const selectedHeight = (highlightedCounts[i] / yMax) * plotHeight;
          ctx.fillStyle = cssColorWithAlpha(theme.axis_color || "#6f7f8f", 0.22, axisColor);
          ctx.fillRect(x, margin.top + plotHeight - fullHeight, Math.max(1, binWidth - 2), fullHeight);
          if (selectedHeight > 0) {
            ctx.fillStyle = cssColorWithAlpha(theme.text_color || "#ffffff", 0.80, "#ffffff");
            ctx.fillRect(x, margin.top + plotHeight - selectedHeight, Math.max(1, binWidth - 2), selectedHeight);
          }
        }

        ctx.fillStyle = axisColor;
        ctx.font = "10px Menlo, Monaco, Consolas, monospace";
        ctx.textAlign = "left";
        ctx.textBaseline = "top";
        ctx.fillText(formatCompactNumber(xMin), margin.left, margin.top + plotHeight + 6);
        ctx.textAlign = "right";
        ctx.fillText(formatCompactNumber(xMax), margin.left + plotWidth, margin.top + plotHeight + 6);

      }

      function dendrogramTraceOptionsForCurrentGroup() {
        const defaults = groupDefaults(currentGroup);
        return (Array.isArray(dendrogramSpec.traces) ? dendrogramSpec.traces : [])
          .filter((traceOption) => {
            const traceKey = String(traceOption && traceOption.trace_key ? traceOption.trace_key : "");
            if (!traceKey) {
              return false;
            }
            return defaults[traceKey] === true;
          })
          .sort((a, b) => String(a.trace_name || a.trace_key || "").localeCompare(String(b.trace_name || b.trace_key || "")));
      }

      function activeDendrogramEntries() {
        const traceKey = dendrogramFocusTraceKey();
        if (!traceKey) {
          return [];
        }
        return (Array.isArray(dendrogramSpec.entries) ? dendrogramSpec.entries : []).filter((entry) => {
          if (!entry || typeof entry !== "object" || String(entry.trace_key || "") !== traceKey) {
            return false;
          }
          return clusterFilterPassesSelectionKey(entry.selection_key);
        });
      }

      function currentDendrogramThresholdMode() {
        return String(dendrogramThresholdMode || "distance_pc") === "birth_age_myr"
          ? "birth_age_myr"
          : "distance_pc";
      }

      function currentDendrogramConnectionMode() {
        return String(dendrogramConnectionMode || "birth_to_older_track") === "birth_to_birth"
          ? "birth_to_birth"
          : "birth_to_older_track";
      }

      function currentDendrogramThresholdValue() {
        return currentDendrogramThresholdMode() === "birth_age_myr"
          ? Math.max(Number(dendrogramThresholdAgeMyr) || 0.0, 0.0)
          : Math.max(Number(dendrogramThresholdPc) || 0.0, 0.0);
      }

      function syncDendrogramThresholdControls() {
        const mode = currentDendrogramThresholdMode();
        if (dendrogramModeEl) {
          dendrogramModeEl.value = mode;
        }
        if (dendrogramThresholdLabelEl) {
          dendrogramThresholdLabelEl.textContent = mode === "birth_age_myr" ? "Threshold (Myr)" : "Threshold (pc)";
        }
        if (!dendrogramThresholdEl) {
          return;
        }
        if (mode === "birth_age_myr") {
          const minValue = Math.max(Number(dendrogramSpec.threshold_min_age_myr) || 0.0, 0.0);
          const maxValue = Math.max(Number(dendrogramSpec.threshold_max_age_myr) || 10.0, minValue + 1.0);
          dendrogramThresholdEl.min = String(minValue);
          dendrogramThresholdEl.max = String(maxValue);
          dendrogramThresholdEl.step = "0.5";
          dendrogramThresholdEl.value = Number(dendrogramThresholdAgeMyr).toFixed(1);
          return;
        }
        const minValue = Math.max(Number(dendrogramSpec.threshold_min_pc) || 0.0, 0.0);
        const maxValue = Math.max(Number(dendrogramSpec.threshold_max_pc) || 1000.0, minValue + 1.0);
        dendrogramThresholdEl.min = String(minValue);
        dendrogramThresholdEl.max = String(maxValue);
        dendrogramThresholdEl.step = "1";
        dendrogramThresholdEl.value = Number(dendrogramThresholdPc).toFixed(0);
      }

      function pointSegmentDistanceSq(px, py, x1, y1, x2, y2) {
        const dx = x2 - x1;
        const dy = y2 - y1;
        if ((dx * dx + dy * dy) <= 1e-12) {
          const ddx = px - x1;
          const ddy = py - y1;
          return ddx * ddx + ddy * ddy;
        }
        const t = clampRange(((px - x1) * dx + (py - y1) * dy) / (dx * dx + dy * dy), 0.0, 1.0);
        const cx = x1 + t * dx;
        const cy = y1 + t * dy;
        const ddx = px - cx;
        const ddy = py - cy;
        return ddx * ddx + ddy * ddy;
      }

      function buildDendrogramModel(entries, thresholdValue, thresholdMode, connectionMode) {
        const mode = String(thresholdMode || "distance_pc") === "birth_age_myr"
          ? "birth_age_myr"
          : "distance_pc";
        const linksMode = String(connectionMode || "birth_to_older_track") === "birth_to_birth"
          ? "birth_to_birth"
          : "birth_to_older_track";
        const threshold = Math.max(Number(thresholdValue) || 0.0, 0.0);
        const sortedEntries = (Array.isArray(entries) ? entries : [])
          .filter((entry) => (
            entry
            && typeof entry === "object"
            && Number.isFinite(Number(entry.age_now_myr))
            && Number.isFinite(Number(entry.birth_time_myr))
            && Number.isFinite(Number(entry.x_birth))
            && Number.isFinite(Number(entry.y_birth))
            && Number.isFinite(Number(entry.z_birth))
            && normalizeMemberKey(entry.selection_key)
          ))
          .map((entry) => ({
            key: normalizeMemberKey(entry.selection_key),
            selection_key: normalizeMemberKey(entry.selection_key),
            cluster_name: String(entry.cluster_name || entry.selection_key || ""),
            trace_name: String(entry.trace_name || ""),
            trace_key: String(entry.trace_key || ""),
            color: String(entry.color || "#ffffff"),
            age_now_myr: Number(entry.age_now_myr),
            birth_time_myr: Number(entry.birth_time_myr),
            x_birth: Number(entry.x_birth),
            y_birth: Number(entry.y_birth),
            z_birth: Number(entry.z_birth),
            time_samples: Array.isArray(entry.time_samples) ? entry.time_samples.map((value) => Number(value)) : [],
            x_samples: Array.isArray(entry.x_samples) ? entry.x_samples.map((value) => Number(value)) : [],
            y_samples: Array.isArray(entry.y_samples) ? entry.y_samples.map((value) => Number(value)) : [],
            z_samples: Array.isArray(entry.z_samples) ? entry.z_samples.map((value) => Number(value)) : [],
            parent_key: "",
            parent_distance_pc: NaN,
            parent_age_gap_myr: NaN,
            axis_value: 0.0,
            children: [],
            plot_order: 0,
          }))
          .sort((a, b) => (
            Number(b.age_now_myr) - Number(a.age_now_myr)
            || String(a.cluster_name).localeCompare(String(b.cluster_name))
          ));

        const nodeByKey = new Map(sortedEntries.map((entry) => [entry.key, entry]));
        for (let index = 0; index < sortedEntries.length; index += 1) {
          const node = sortedEntries[index];
          let bestParent = null;
          let bestDistance = Infinity;
          let bestAgeGap = Infinity;
          for (let olderIndex = 0; olderIndex < index; olderIndex += 1) {
            const candidate = sortedEntries[olderIndex];
            const ageGap = Math.max(0.0, Number(candidate.age_now_myr) - Number(node.age_now_myr));
            const candidatePosition = linksMode === "birth_to_birth"
              ? {
                x: Number(candidate.x_birth),
                y: Number(candidate.y_birth),
                z: Number(candidate.z_birth),
              }
              : interpolateDendrogramPosition(candidate, node.birth_time_myr);
            if (!candidatePosition) {
              continue;
            }
            const dx = node.x_birth - candidatePosition.x;
            const dy = node.y_birth - candidatePosition.y;
            const dz = node.z_birth - candidatePosition.z;
            const distance = Math.sqrt(dx * dx + dy * dy + dz * dz);
            if (mode === "birth_age_myr") {
              if (!(ageGap <= threshold) || distance >= bestDistance) {
                continue;
              }
            } else if (!(distance <= threshold) || distance >= bestDistance) {
              continue;
            }
            bestDistance = distance;
            bestAgeGap = ageGap;
            bestParent = candidate;
          }
          if (bestParent) {
            node.parent_key = bestParent.key;
            node.parent_distance_pc = bestDistance;
            node.parent_age_gap_myr = bestAgeGap;
            bestParent.children.push(node.key);
          }
          node.axis_value = mode === "birth_age_myr"
            ? 0.0
            : Math.max(Number(node.age_now_myr) || 0.0, 0.0);
        }

        if (mode === "birth_age_myr") {
          const cumulativeAxisByKey = new Map();
          function cumulativeBirthDistanceFor(node) {
            if (!node) {
              return 0.0;
            }
            if (cumulativeAxisByKey.has(node.key)) {
              return cumulativeAxisByKey.get(node.key);
            }
            const localDistance = Number.isFinite(node.parent_distance_pc)
              ? Math.max(Number(node.parent_distance_pc), 0.0)
              : 0.0;
            const parentDistance = node.parent_key
              ? cumulativeBirthDistanceFor(nodeByKey.get(node.parent_key))
              : 0.0;
            const cumulativeDistance = parentDistance + localDistance;
            cumulativeAxisByKey.set(node.key, cumulativeDistance);
            return cumulativeDistance;
          }

          sortedEntries.forEach((node) => {
            node.axis_value = cumulativeBirthDistanceFor(node);
          });
        }

        const roots = sortedEntries
          .filter((entry) => !entry.parent_key)
          .sort((a, b) => (
            Number(b.age_now_myr) - Number(a.age_now_myr)
            || String(a.cluster_name).localeCompare(String(b.cluster_name))
          ));

        let nextLeafOrder = 0;
        function assignPlotOrder(node) {
          const children = node.children
            .map((childKey) => nodeByKey.get(childKey))
            .filter(Boolean)
            .sort((a, b) => (
              Number(b.age_now_myr) - Number(a.age_now_myr)
              || String(a.cluster_name).localeCompare(String(b.cluster_name))
            ));
          if (!children.length) {
            node.plot_order = nextLeafOrder;
            nextLeafOrder += 1;
            return node.plot_order;
          }
          const childOrders = children.map(assignPlotOrder);
          node.plot_order = childOrders.reduce((total, value) => total + value, 0.0) / childOrders.length;
          return node.plot_order;
        }
        roots.forEach(assignPlotOrder);

        const descendantKeysByNode = new Map();
        function descendantKeysFor(node) {
          if (!node) {
            return [];
          }
          if (descendantKeysByNode.has(node.key)) {
            return descendantKeysByNode.get(node.key);
          }
          const keys = [node.key];
          node.children.forEach((childKey) => {
            const child = nodeByKey.get(childKey);
            descendantKeysFor(child).forEach((selectionKey) => keys.push(selectionKey));
          });
          const uniqueKeys = Array.from(new Set(keys.map((value) => normalizeMemberKey(value)).filter(Boolean)));
          descendantKeysByNode.set(node.key, uniqueKeys);
          return uniqueKeys;
        }

        const branches = [];
        sortedEntries.forEach((node) => {
          if (!node.parent_key) {
            return;
          }
          const parent = nodeByKey.get(node.parent_key);
          if (!parent) {
            return;
          }
          branches.push({
            key: `${node.key}->${parent.key}`,
            child: node,
            parent,
            label: `${node.cluster_name} from ${parent.cluster_name}`,
            selectionKeys: descendantKeysFor(node),
            count: descendantKeysFor(node).length,
            distance_pc: Number(node.parent_distance_pc),
            age_gap_myr: Number(node.parent_age_gap_myr),
          });
        });

        const axisValues = sortedEntries
          .map((entry) => Number(entry.axis_value))
          .filter((value) => Number.isFinite(value) && value >= 0.0);

        return {
          nodes: sortedEntries,
          roots,
          branches,
          leafCount: Math.max(nextLeafOrder, 1),
          maxAgeMyr: Math.max(...sortedEntries.map((entry) => Number(entry.age_now_myr)), 0.0),
          maxAxisValue: axisValues.length ? Math.max(...axisValues, 0.0) : 0.0,
          connectionMode: linksMode,
          thresholdMode: mode,
        };
      }

      function renderDendrogramWidget() {
        if (!dendrogramSpec.enabled || !dendrogramCanvasEl || widgetModeForKey("dendrogram") === "hidden") {
          dendrogramHitRegions = [];
          return;
        }

        const traceOptions = dendrogramTraceOptionsForCurrentGroup();
        dendrogramTraceEl.innerHTML = "";
        traceOptions.forEach((traceOption) => {
          const option = document.createElement("option");
          option.value = String(traceOption.trace_key || "");
          option.textContent = String(traceOption.trace_name || traceOption.trace_key || "");
          dendrogramTraceEl.appendChild(option);
        });
        if (!traceOptions.some((traceOption) => String(traceOption.trace_key || "") === String(dendrogramTraceKey || ""))) {
          dendrogramTraceKey = traceOptions.length ? String(traceOptions[0].trace_key || "") : "";
          clearDendrogramHoverState();
        }
        dendrogramTraceEl.value = dendrogramTraceKey;
        dendrogramConnectionMode = currentDendrogramConnectionMode();
        if (dendrogramConnectionEl) {
          dendrogramConnectionEl.value = dendrogramConnectionMode;
        }
        dendrogramThresholdPc = Math.max(Number(dendrogramThresholdPc) || 0.0, 0.0);
        dendrogramThresholdAgeMyr = Math.max(Number(dendrogramThresholdAgeMyr) || 0.0, 0.0);
        syncDendrogramThresholdControls();

        const entries = activeDendrogramEntries();
        const rect = dendrogramCanvasEl.getBoundingClientRect();
        const cssWidth = Math.max(1, Math.round(rect.width));
        const cssHeight = Math.max(1, Math.round(rect.height));
        const dpr = Math.max(window.devicePixelRatio || 1, 1);
        if (dendrogramCanvasEl.width !== Math.round(cssWidth * dpr) || dendrogramCanvasEl.height !== Math.round(cssHeight * dpr)) {
          dendrogramCanvasEl.width = Math.round(cssWidth * dpr);
          dendrogramCanvasEl.height = Math.round(cssHeight * dpr);
        }
        const ctx = dendrogramCanvasEl.getContext("2d");
        if (!ctx) {
          dendrogramHitRegions = [];
          return;
        }
        ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
        ctx.clearRect(0, 0, cssWidth, cssHeight);
        ctx.fillStyle = theme.scene_bgcolor || theme.paper_bgcolor || "#000000";
        ctx.fillRect(0, 0, cssWidth, cssHeight);

        if (!traceOptions.length || !entries.length) {
          dendrogramHitRegions = [];
          clearDendrogramSelectionState();
          return;
        }

        const thresholdMode = currentDendrogramThresholdMode();
        const thresholdValue = currentDendrogramThresholdValue();
        const connectionMode = currentDendrogramConnectionMode();
        const model = buildDendrogramModel(entries, thresholdValue, thresholdMode, connectionMode);
        const margin = { left: 52, right: 16, top: 18, bottom: 28 };
        const plotWidth = Math.max(40, cssWidth - margin.left - margin.right);
        const plotHeight = Math.max(40, cssHeight - margin.top - margin.bottom);
        const currentTraceMeta = traceOptions.find((traceOption) => String(traceOption.trace_key || "") === String(dendrogramTraceKey || "")) || traceOptions[0];
        const traceMaxAge = Math.max(Number(currentTraceMeta && currentTraceMeta.max_age_myr) || 0.0, Number(model.maxAgeMyr) || 0.0, 0.0);
        const axisMode = model.thresholdMode || thresholdMode;
        const axisLabel = axisMode === "birth_age_myr" ? "Birth distance" : "Birth age";
        const axisUnit = axisMode === "birth_age_myr" ? "pc" : "Myr";
        const axisMax = axisMode === "birth_age_myr"
          ? Math.max((Number(model.maxAxisValue) || 0.0) * 1.08, 25.0)
          : Math.max(traceMaxAge + 5.0, 1.0);
        const traceState = traceStyleStateForKey(dendrogramTraceKey);
        const traceColor = (traceState && traceState.color) || String(currentTraceMeta.color || theme.text_color || "#ffffff");
        const axisColor = String(theme.axis_color || "#808080");
        const activeKeys = activeDendrogramSelectionKeys();
        const hoveredClusterKeys = activeHoveredClusterKeys();
        dendrogramHitRegions = [];

        function xToPx(orderValue) {
          if (model.leafCount <= 1) {
            return margin.left + plotWidth * 0.5;
          }
          return margin.left + (Number(orderValue) / Math.max(model.leafCount - 1, 1)) * plotWidth;
        }

        function yToPx(axisValue) {
          return margin.top + (1.0 - clampRange(Number(axisValue) / axisMax, 0.0, 1.0)) * plotHeight;
        }

        ctx.strokeStyle = axisColor;
        ctx.lineWidth = 1.2;
        ctx.beginPath();
        ctx.moveTo(margin.left, margin.top);
        ctx.lineTo(margin.left, margin.top + plotHeight);
        ctx.lineTo(margin.left + plotWidth, margin.top + plotHeight);
        ctx.stroke();

        ctx.fillStyle = axisColor;
        ctx.font = "10px Menlo, Monaco, Consolas, monospace";
        ctx.textAlign = "right";
        ctx.textBaseline = "middle";
        [0.0, 0.33, 0.66, 1.0].forEach((fraction) => {
          const axisValue = fraction * axisMax;
          const yPx = yToPx(axisValue);
          ctx.globalAlpha = 0.14;
          ctx.beginPath();
          ctx.moveTo(margin.left, yPx);
          ctx.lineTo(margin.left + plotWidth, yPx);
          ctx.stroke();
          ctx.globalAlpha = 1.0;
          ctx.fillText(`${formatCompactNumber(axisValue)} ${axisUnit}`, margin.left - 8, yPx);
        });

        ctx.save();
        ctx.translate(14, margin.top + plotHeight * 0.5);
        ctx.rotate(-Math.PI * 0.5);
        ctx.textAlign = "center";
        ctx.textBaseline = "middle";
        ctx.fillText(axisLabel, 0, 0);
        ctx.restore();

        const formationAge = Math.max(0.0, -(currentFrame() ? Number(currentFrame().time) : 0.0));
        if (axisMode === "distance_pc") {
          const markerY = yToPx(Math.min(formationAge, axisMax));
          ctx.save();
          ctx.setLineDash([6, 6]);
          ctx.strokeStyle = axisColor;
          ctx.lineWidth = 1.4;
          ctx.beginPath();
          ctx.moveTo(margin.left, markerY);
          ctx.lineTo(margin.left + plotWidth, markerY);
          ctx.stroke();
          ctx.restore();
        } else {
          ctx.save();
          ctx.textAlign = "right";
          ctx.textBaseline = "top";
          ctx.fillStyle = cssColorWithAlpha(axisColor, 0.9, axisColor);
          ctx.fillText(`Formation age: ${formatCompactNumber(formationAge)} Myr`, margin.left + plotWidth, margin.top + 2);
          ctx.restore();
        }

        model.branches.forEach((branch) => {
          const childX = xToPx(branch.child.plot_order);
          const childY = yToPx(branch.child.axis_value);
          const parentX = xToPx(branch.parent.plot_order);
          const parentY = yToPx(branch.parent.axis_value);
          const isPinned = dendrogramPinnedRegionKey && dendrogramPinnedRegionKey === branch.key;
          const isHovered = !isPinned && dendrogramHoveredRegionKey === branch.key;
          const isActive = activeKeys.size && branch.selectionKeys.some((selectionKey) => activeKeys.has(selectionKey));
          ctx.strokeStyle = cssColorWithAlpha(traceColor, isPinned ? 1.0 : (isHovered ? 0.96 : (isActive ? 0.86 : 0.62)), traceColor);
          ctx.lineWidth = isPinned ? 3.6 : (isHovered ? 3.0 : (isActive ? 2.4 : 1.7));
          ctx.beginPath();
          ctx.moveTo(childX, childY);
          ctx.lineTo(childX, parentY);
          ctx.lineTo(parentX, parentY);
          ctx.stroke();
          dendrogramHitRegions.push({
            type: "branch",
            key: branch.key,
            label: branch.label,
            count: branch.count,
            selectionKeys: branch.selectionKeys,
            hitRadius: isPinned ? 10.0 : (isHovered ? 9.0 : 8.0),
            segments: [
              [childX, childY, childX, parentY],
              [childX, parentY, parentX, parentY],
            ],
          });
        });

        model.nodes.forEach((node) => {
          const nodeX = xToPx(node.plot_order);
          const nodeY = yToPx(node.axis_value);
          const nodeRegionKey = `node:${node.key}`;
          const isPinned = dendrogramPinnedRegionKey && dendrogramPinnedRegionKey === nodeRegionKey;
          const isHovered = !isPinned && dendrogramHoveredRegionKey === nodeRegionKey;
          const isActive = activeKeys.size && activeKeys.has(node.selection_key);
          const isSceneHovered = hoveredClusterKeys.has(node.selection_key);
          const nodeRadius = isPinned ? 5.2 : (isHovered ? 4.8 : (isSceneHovered ? 4.6 : (isActive ? 4.1 : 3.4)));
          ctx.beginPath();
          ctx.fillStyle = cssColorWithAlpha(traceColor, isPinned ? 1.0 : (isHovered ? 0.96 : (isSceneHovered ? 0.94 : (isActive ? 0.88 : 0.82))), traceColor);
          ctx.arc(nodeX, nodeY, nodeRadius, 0, Math.PI * 2.0);
          ctx.fill();
          if (isSceneHovered && !isPinned) {
            ctx.save();
            ctx.strokeStyle = cssColorWithAlpha(traceColor, 0.95, traceColor);
            ctx.lineWidth = 1.4;
            ctx.beginPath();
            ctx.arc(nodeX, nodeY, nodeRadius + 2.6, 0, Math.PI * 2.0);
            ctx.stroke();
            ctx.restore();
          }
          dendrogramHitRegions.push({
            type: "node",
            key: nodeRegionKey,
            label: node.cluster_name,
            count: 1,
            selectionKeys: [node.selection_key],
            centerX: nodeX,
            centerY: nodeY,
            radius: Math.max(7.5, nodeRadius + 2.5),
          });
        });

      }
""".strip()
