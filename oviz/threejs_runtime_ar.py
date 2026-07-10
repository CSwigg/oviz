from .threejs_runtime_usdz import THREEJS_AR_USDZ_RUNTIME_JS


THREEJS_AR_RUNTIME_JS = THREEJS_AR_USDZ_RUNTIME_JS + "\n\n" + r"""
      const OVIZ_AR_MAX_POINTS = 240;
      const OVIZ_AR_MAX_VOLUME_SAMPLES = 900;
      const OVIZ_AR_VOLUME_SCAN_TARGET = 120000;
      const OVIZ_AR_VOLUME_COLOR_BINS = 8;
      const OVIZ_AR_VOLUME_MIN_STRETCHED_VALUE = 0.02;
      const OVIZ_AR_VOLUME_BUCKET_SAMPLES = 2;
      const OVIZ_AR_MIN_USDZ_BYTES = 1024;
      const OVIZ_AR_QUICKLOOK_CACHE = "oviz-ar-quicklook-v2";
      const OVIZ_AR_QUICKLOOK_WORKER_VERSION = "20260710_soft_volume_v4";
      const OVIZ_AR_SKY_TEXTURE_HIGH = { width: 4096, height: 2048 };
      const OVIZ_AR_SKY_TEXTURE_LOW = { width: 2048, height: 1024 };
      let ovizArDialogEl = null;
      let ovizArObjectUrl = "";
      let ovizArQuickLookUrl = "";
      let ovizArExportBusy = false;

      function ovizArFiniteNumber(value, fallback = NaN) {
        const number = Number(value);
        return Number.isFinite(number) ? number : fallback;
      }

      function ovizArInitialFrameIndexFallback() {
        return Math.max(
          0,
          Math.min(
            Math.round(Number(sceneSpec.initial_frame_index) || 0),
            Math.max((Array.isArray(frameSpecs) ? frameSpecs.length : 0) - 1, 0)
          )
        );
      }

      function ovizArPresentFrameIndex() {
        const frames = Array.isArray(frameSpecs) ? frameSpecs : [];
        for (let index = 0; index < frames.length; index += 1) {
          const frameTime = Number(frames[index] && frames[index].time);
          if (Number.isFinite(frameTime) && Math.abs(frameTime) <= 1e-9) {
            return index;
          }
        }
        return ovizArInitialFrameIndexFallback();
      }

      function ovizArPresentFrame() {
        const frames = Array.isArray(frameSpecs) ? frameSpecs : [];
        return frames[ovizArPresentFrameIndex()] || frames[0] || null;
      }

      function ovizArSelectedClusterKeys() {
        const keys = new Set();
        if (selectedClusterKeys && typeof selectedClusterKeys.forEach === "function") {
          selectedClusterKeys.forEach((key) => {
            const normalized = normalizeMemberKey(key);
            if (normalized) {
              keys.add(normalized);
            }
          });
        }
        const addSelection = (selection) => {
          const key = normalizedSelectionKeyFor(selection);
          if (key) {
            keys.add(key);
          }
        };
        if (currentSelection) {
          addSelection(currentSelection);
        }
        (Array.isArray(currentSelections) ? currentSelections : []).forEach(addSelection);
        return keys;
      }

      function ovizArSelectionObjectsByKey() {
        const selections = new Map();
        const addSelection = (selection) => {
          const key = normalizedSelectionKeyFor(selection);
          if (key && !selections.has(key)) {
            selections.set(key, selection);
          }
        };
        if (currentSelection) {
          addSelection(currentSelection);
        }
        (Array.isArray(currentSelections) ? currentSelections : []).forEach(addSelection);
        ovizArSelectedClusterKeys().forEach((key) => {
          if (!selections.has(key)) {
            const metadata = selectionMetadataForKey(key);
            if (metadata) {
              selections.set(key, metadata);
            }
          }
        });
        return selections;
      }

      function ovizArCanExportSelection() {
        return !minimalModeEnabled;
      }

      function renderArSnapshotButtonState() {
        if (!mobileArButtonEl) {
          return;
        }
        const selectedCount = ovizArSelectedClusterKeys().size;
        const enabled = !minimalModeEnabled && !ovizArExportBusy;
        mobileArButtonEl.disabled = !enabled;
        mobileArButtonEl.dataset.hasSelection = selectedCount > 0 ? "true" : "false";
        mobileArButtonEl.dataset.active = ovizArDialogEl && ovizArDialogEl.dataset.open === "true" ? "true" : "false";
        mobileArButtonEl.setAttribute("aria-expanded", mobileArButtonEl.dataset.active === "true" ? "true" : "false");
        mobileArButtonEl.title = selectedCount
          ? `Export ${selectedCount} selected cluster${selectedCount === 1 ? "" : "s"} to AR Quick Look`
          : "Open AR export options for the present-day scene.";
      }

      function ovizArSelectionLabel(selection, fallback = "Cluster") {
        return String(
          (selection && (selection.cluster_name || selection.trace_name))
          || fallback
          || "Cluster"
        ).trim() || "Cluster";
      }

      function ovizArPointSelectionKey(point, trace) {
        const selection = selectionForPoint(point, trace);
        return selection ? normalizedSelectionKeyFor(selection) : "";
      }

      function ovizArPointCoordinates(point, selection) {
        const x = ovizArFiniteNumber(point && point.x, ovizArFiniteNumber(selection && selection.x0));
        const y = ovizArFiniteNumber(point && point.y, ovizArFiniteNumber(selection && selection.y0));
        const z = ovizArFiniteNumber(point && point.z, ovizArFiniteNumber(selection && selection.z0));
        if (![x, y, z].every(Number.isFinite)) {
          return null;
        }
        return { x, y, z };
      }

      function ovizArPlotGroupOffset() {
        const offset = (typeof plotGroup !== "undefined" && plotGroup && plotGroup.position)
          ? plotGroup.position
          : { x: 0.0, y: 0.0, z: 0.0 };
        return {
          x: ovizArFiniteNumber(offset.x, 0.0),
          y: ovizArFiniteNumber(offset.y, 0.0),
          z: ovizArFiniteNumber(offset.z, 0.0),
        };
      }

      function ovizArActiveSelectionMask() {
        if (
          typeof hasActiveLassoSelectionMask === "function"
          && hasActiveLassoSelectionMask()
          && typeof currentLassoSelectionMask !== "undefined"
        ) {
          return currentLassoSelectionMask;
        }
        return null;
      }

      function ovizArPointInsideProjectedMask(point, mask) {
        if (!mask || !point || typeof pointInsideProjectedLassoMask !== "function") {
          return true;
        }
        const offset = ovizArPlotGroupOffset();
        return pointInsideProjectedLassoMask(
          ovizArFiniteNumber(point.x, 0.0) + offset.x,
          ovizArFiniteNumber(point.y, 0.0) + offset.y,
          ovizArFiniteNumber(point.z, 0.0) + offset.z,
          mask
        );
      }

      function ovizArSelectionSkyCoordinates(selection, point = null) {
        const source = selection && typeof selection === "object" ? selection : {};
        const fallback = point && typeof point === "object" ? point : {};
        const ra = ovizArFiniteNumber(source.ra_deg, ovizArFiniteNumber(fallback.ra_deg, ovizArFiniteNumber(fallback.ra)));
        const dec = ovizArFiniteNumber(source.dec_deg, ovizArFiniteNumber(fallback.dec_deg, ovizArFiniteNumber(fallback.dec)));
        const l = ovizArFiniteNumber(source.l_deg, ovizArFiniteNumber(fallback.l_deg, ovizArFiniteNumber(fallback.l)));
        const b = ovizArFiniteNumber(source.b_deg, ovizArFiniteNumber(fallback.b_deg, ovizArFiniteNumber(fallback.b)));
        return { ra, dec, l, b };
      }

      function ovizArPointColor(point, trace, selection) {
        return String(
          (selection && selection.cluster_color)
          || pointBaseColorForTrace(point, trace)
          || "#ffffff"
        );
      }

      function ovizArNormalizePointRecord(point, trace, key, selection, frame, frameIndex) {
        const coords = ovizArPointCoordinates(point, selection);
        if (!coords) {
          return null;
        }
        const skyCoords = ovizArSelectionSkyCoordinates(selection, point);
        const traceName = String((trace && trace.name) || "");
        const label = ovizArSelectionLabel(selection, traceName || key);
        const nStars = ovizArFiniteNumber(
          selection && selection.n_stars,
          ovizArFiniteNumber(point && point.n_stars, NaN)
        );
        return {
          key,
          label,
          traceName,
          color: ovizArPointColor(point, trace, selection),
          x: coords.x,
          y: coords.y,
          z: coords.z,
          ra: skyCoords.ra,
          dec: skyCoords.dec,
          l: skyCoords.l,
          b: skyCoords.b,
          nStars,
          frameIndex,
          timeMyr: ovizArFiniteNumber(frame && frame.time, 0.0),
        };
      }

      function ovizArPointInsideVolumeBounds(point, layer) {
        const bounds = layer && layer.bounds ? layer.bounds : {};
        return ["x", "y", "z"].every((axis) => {
          const axisBounds = Array.isArray(bounds[axis]) ? bounds[axis] : null;
          const value = ovizArFiniteNumber(point && point[axis]);
          if (!axisBounds || axisBounds.length < 2 || !Number.isFinite(value)) {
            return false;
          }
          const low = Math.min(ovizArFiniteNumber(axisBounds[0], 0.0), ovizArFiniteNumber(axisBounds[1], 0.0));
          const high = Math.max(ovizArFiniteNumber(axisBounds[0], 0.0), ovizArFiniteNumber(axisBounds[1], 0.0));
          return value >= low && value <= high;
        });
      }

      function ovizArStringHash(value) {
        const text = String(value || "");
        let hash = 2166136261;
        for (let index = 0; index < text.length; index += 1) {
          hash ^= text.charCodeAt(index);
          hash = Math.imul(hash, 16777619);
        }
        return hash >>> 0;
      }

      function ovizArRepresentativePoints(points, maxPoints) {
        const source = Array.isArray(points) ? points : [];
        const cap = Math.max(1, Math.round(Number(maxPoints) || 1));
        if (source.length <= cap) {
          return source.slice();
        }
        const bounds = {
          x: [Infinity, -Infinity],
          y: [Infinity, -Infinity],
          z: [Infinity, -Infinity],
        };
        source.forEach((point) => {
          ["x", "y", "z"].forEach((axis) => {
            const value = ovizArFiniteNumber(point && point[axis], 0.0);
            bounds[axis][0] = Math.min(bounds[axis][0], value);
            bounds[axis][1] = Math.max(bounds[axis][1], value);
          });
        });
        const cellsPerAxis = Math.max(2, Math.ceil(Math.cbrt(cap * 1.5)));
        const cellWinners = new Map();
        source.forEach((point) => {
          const indices = ["x", "y", "z"].map((axis) => {
            const span = Math.max(bounds[axis][1] - bounds[axis][0], 1e-9);
            const normalized = (ovizArFiniteNumber(point && point[axis], bounds[axis][0]) - bounds[axis][0]) / span;
            return Math.max(0, Math.min(cellsPerAxis - 1, Math.floor(normalized * cellsPerAxis)));
          });
          const cellKey = indices.join(":");
          const incumbent = cellWinners.get(cellKey);
          const pointStars = ovizArFiniteNumber(point && point.nStars, -Infinity);
          const incumbentStars = ovizArFiniteNumber(incumbent && incumbent.nStars, -Infinity);
          if (
            !incumbent
            || pointStars > incumbentStars
            || (pointStars === incumbentStars && ovizArStringHash(point.key) < ovizArStringHash(incumbent.key))
          ) {
            cellWinners.set(cellKey, point);
          }
        });
        const selected = Array.from(cellWinners.values())
          .sort((left, right) => ovizArStringHash(left.key) - ovizArStringHash(right.key));
        if (selected.length >= cap) {
          return selected.slice(0, cap);
        }
        const selectedKeys = new Set(selected.map((point) => point.key));
        source
          .filter((point) => !selectedKeys.has(point.key))
          .sort((left, right) => ovizArStringHash(left.key) - ovizArStringHash(right.key))
          .slice(0, cap - selected.length)
          .forEach((point) => selected.push(point));
        return selected;
      }

      function collectOvizArSnapshot(mode = "3d") {
        const selectedKeys = ovizArSelectedClusterKeys();
        const selectionMask = ovizArActiveSelectionMask();
        const maskOnlySelection = selectedKeys.size === 0 && Boolean(selectionMask);
        const selectionObjects = ovizArSelectionObjectsByKey();
        const frameIndex = ovizArPresentFrameIndex();
        const frame = ovizArPresentFrame();
        const points = [];
        const seen = new Set();
        const visibleVolumeLayers = (
          String(mode || "3d") !== "sky"
          && selectedKeys.size === 0
          && !selectionMask
        ) ? ovizArVolumeLayersForSnapshot(frame) : [];
        if (frame) {
          const traces = Array.isArray(frame.traces) ? frame.traces : [];
          traces.forEach((trace) => {
            const tracePoints = trace && Array.isArray(trace.points) ? trace.points : [];
            tracePoints.forEach((point) => {
              const selection = selectionForPoint(point, trace);
              const key = selection ? normalizedSelectionKeyFor(selection) : "";
              if (!key || (selectedKeys.size > 0 && !selectedKeys.has(key)) || seen.has(key)) {
                return;
              }
              const mergedSelection = Object.assign(
                {},
                selectionObjects.get(key) || {},
                selection || {}
              );
              const record = ovizArNormalizePointRecord(point, trace, key, mergedSelection, frame, frameIndex);
              if (!record) {
                return;
              }
              if (maskOnlySelection && !ovizArPointInsideProjectedMask(record, selectionMask)) {
                return;
              }
              if (
                visibleVolumeLayers.length
                && !visibleVolumeLayers.some((layer) => ovizArPointInsideVolumeBounds(record, layer))
              ) {
                return;
              }
              seen.add(key);
              points.push(record);
            });
          });
        }
        const cappedPoints = ovizArRepresentativePoints(points, OVIZ_AR_MAX_POINTS);
        return {
          mode: String(mode || "3d"),
          selectionMode: selectedKeys.size > 0 ? "selection" : (maskOnlySelection ? "volume-lasso" : "present-day-scene"),
          selectedCount: selectedKeys.size,
          hasLassoMask: Boolean(selectionMask),
          presentFrameIndex: frameIndex,
          presentTimeMyr: ovizArFiniteNumber(frame && frame.time, 0.0),
          points: cappedPoints,
          trails: [],
          pointLimit: OVIZ_AR_MAX_POINTS,
          truncated: points.length > cappedPoints.length,
        };
      }

      function ovizArVolumeBoundsPoints(layers) {
        const points = [];
        (Array.isArray(layers) ? layers : []).forEach((layer) => {
          const bounds = layer && layer.bounds ? layer.bounds : {};
          const x = Array.isArray(bounds.x) ? bounds.x.map(Number) : [];
          const y = Array.isArray(bounds.y) ? bounds.y.map(Number) : [];
          const z = Array.isArray(bounds.z) ? bounds.z.map(Number) : [];
          if (x.length < 2 || y.length < 2 || z.length < 2 || !x.concat(y, z).every(Number.isFinite)) {
            return;
          }
          [x[0], x[1]].forEach((xValue) => {
            [y[0], y[1]].forEach((yValue) => {
              [z[0], z[1]].forEach((zValue) => points.push({ x: xValue, y: yValue, z: zValue }));
            });
          });
        });
        return points;
      }

      function ovizArSceneTransform(points, targetRadiusMeters = 0.75, volumeLayers = []) {
        const source = (Array.isArray(points) ? points : []).concat(ovizArVolumeBoundsPoints(volumeLayers));
        const bounds = {
          x: [Infinity, -Infinity],
          y: [Infinity, -Infinity],
          z: [Infinity, -Infinity],
        };
        source.forEach((point) => {
          ["x", "y", "z"].forEach((axis) => {
            const value = ovizArFiniteNumber(point && point[axis]);
            if (Number.isFinite(value)) {
              bounds[axis][0] = Math.min(bounds[axis][0], value);
              bounds[axis][1] = Math.max(bounds[axis][1], value);
            }
          });
        });
        ["x", "y", "z"].forEach((axis) => {
          if (!bounds[axis].every(Number.isFinite)) {
            bounds[axis] = [-0.5, 0.5];
          }
        });
        const center = {
          x: 0.5 * (bounds.x[0] + bounds.x[1]),
          y: 0.5 * (bounds.y[0] + bounds.y[1]),
          z: 0.5 * (bounds.z[0] + bounds.z[1]),
        };
        const halfX = 0.5 * Math.max(bounds.x[1] - bounds.x[0], 1e-6);
        const halfY = 0.5 * Math.max(bounds.y[1] - bounds.y[0], 1e-6);
        const halfZ = 0.5 * Math.max(bounds.z[1] - bounds.z[0], 1e-6);
        const maxDistance = Math.max(Math.sqrt((halfX * halfX) + (halfY * halfY) + (halfZ * halfZ)), 1e-6);
        const scale = targetRadiusMeters / maxDistance;
        return {
          bounds,
          center,
          scale,
          yOffset: 0.045 - ((bounds.z[0] - center.z) * scale),
        };
      }

      function ovizArVectorFromPoint(point, transform) {
        const center = transform && transform.center ? transform.center : { x: 0, y: 0, z: 0 };
        const scale = Number(transform && transform.scale) || 1.0;
        const x = (ovizArFiniteNumber(point && point.x, 0.0) - center.x) * scale;
        const y = (ovizArFiniteNumber(point && point.y, 0.0) - center.y) * scale;
        const z = (ovizArFiniteNumber(point && point.z, 0.0) - center.z) * scale;
        return new THREE.Vector3(x, z + ovizArFiniteNumber(transform && transform.yOffset, 0.0), -y);
      }

      function ovizArColor(value, fallback = "#ffffff") {
        try {
          return new THREE.Color(value || fallback);
        } catch (_err) {
          return new THREE.Color(fallback);
        }
      }

      function ovizArVolumeStateForLayer(layer) {
        const stateKey = typeof volumeStateKeyForLayer === "function"
          ? volumeStateKeyForLayer(layer)
          : String((layer && (layer.state_key || layer.key)) || "");
        if (!stateKey || typeof volumeStateByKey === "undefined") {
          return null;
        }
        return volumeStateByKey[stateKey] || null;
      }

      function ovizArVolumeLayerVisible(layer, frame) {
        if (!layer) {
          return false;
        }
        const state = ovizArVolumeStateForLayer(layer);
        if (!state || state.visible === false) {
          return false;
        }
        const stateKey = typeof volumeStateKeyForLayer === "function"
          ? volumeStateKeyForLayer(layer)
          : String((layer && (layer.state_key || layer.key)) || "");
        if (
          typeof legendState !== "undefined"
          && stateKey
          && (legendState[stateKey] === false || legendState[String(layer.key || "")] === false)
        ) {
          return false;
        }
        return typeof volumeVisibleForFrame === "function"
          ? volumeVisibleForFrame(layer, state, frame)
          : true;
      }

      function ovizArVolumeLayersForSnapshot(frame = ovizArPresentFrame()) {
        const layers = [];
        const seen = new Set();
        const addLayer = (layer) => {
          if (!layer) {
            return;
          }
          const stateKey = typeof volumeStateKeyForLayer === "function"
            ? volumeStateKeyForLayer(layer)
            : String(layer.state_key || layer.key || "");
          const key = stateKey || String(layer.key || "");
          if (!key || seen.has(key) || !ovizArVolumeLayerVisible(layer, frame)) {
            return;
          }
          seen.add(key);
          layers.push(layer);
        };
        if (typeof activeVolumeKey !== "undefined" && activeVolumeKey) {
          if (typeof frameVolumeLayerForStateKey === "function") {
            addLayer(frameVolumeLayerForStateKey(activeVolumeKey, frame));
          }
          if (typeof volumeLayerForKey === "function") {
            addLayer(volumeLayerForKey(activeVolumeKey));
          }
        }
        if (typeof frameVolumeLayers === "function") {
          frameVolumeLayers(frame).forEach(addLayer);
        }
        return layers;
      }

      function ovizArHasVolumeSnapshot(frame = ovizArPresentFrame()) {
        return ovizArVolumeLayersForSnapshot(frame).length > 0;
      }

      async function ovizArVolumeScalarArrayFor(layer) {
        if (!layer || typeof volumeScalarArrayFor !== "function") {
          return null;
        }
        const layerKey = String(layer.key || "");
        const data = volumeScalarArrayFor(layer);
        if (String(layer.data_encoding || "uint8") === "png_atlas_uint8") {
          const pending = (
            typeof volumeScalarDataPendingCache !== "undefined"
            && volumeScalarDataPendingCache
            && typeof volumeScalarDataPendingCache.get === "function"
          ) ? volumeScalarDataPendingCache.get(layerKey) : null;
          if (pending && typeof pending.then === "function") {
            const decoded = await pending;
            if (decoded && decoded.length) {
              return decoded;
            }
          }
          if (
            typeof volumeScalarDataCache !== "undefined"
            && volumeScalarDataCache
            && typeof volumeScalarDataCache.get === "function"
          ) {
            const cached = volumeScalarDataCache.get(layerKey);
            if (cached && cached.length) {
              return cached;
            }
          }
        }
        return data && data.length ? data : null;
      }

      function ovizArVolumeWindowFor(layer, state) {
        if (typeof normalizedVolumeWindowFor === "function") {
          return normalizedVolumeWindowFor(layer, state);
        }
        return { low: 0.0, high: 1.0 };
      }

      function ovizArVolumeStretchValue(value, stretchName) {
        const clamped = Math.min(Math.max(Number(value) || 0.0, 0.0), 1.0);
        const stretch = String(stretchName || "linear").toLowerCase();
        if (stretch === "log10") {
          const strength = 999.0;
          return Math.log(1.0 + strength * clamped) / Math.log(1.0 + strength);
        }
        if (stretch === "asinh") {
          const strength = 10.0;
          return Math.log((strength * clamped) + Math.sqrt((strength * clamped * strength * clamped) + 1.0))
            / Math.log(strength + Math.sqrt((strength * strength) + 1.0));
        }
        return clamped;
      }

      function ovizArVolumeColorBytesForLayer(layer, state) {
        const option = typeof volumeColormapOptionFor === "function"
          ? volumeColormapOptionFor(layer, state && state.colormap)
          : (((layer && layer.colormap_options) || [])[0] || null);
        if (!option) {
          return null;
        }
        if (typeof volumeColorBytesForOption === "function") {
          return volumeColorBytesForOption(option);
        }
        if (typeof base64ToUint8Array === "function") {
          return base64ToUint8Array(option.lut_b64 || "");
        }
        return null;
      }

      function ovizArHexColorFromBytes(bytes, offset) {
        const r = Math.max(0, Math.min(255, Math.round(Number(bytes && bytes[offset]) || 0)));
        const g = Math.max(0, Math.min(255, Math.round(Number(bytes && bytes[offset + 1]) || 0)));
        const b = Math.max(0, Math.min(255, Math.round(Number(bytes && bytes[offset + 2]) || 0)));
        return `#${r.toString(16).padStart(2, "0")}${g.toString(16).padStart(2, "0")}${b.toString(16).padStart(2, "0")}`;
      }

      function ovizArVolumeSamplePosition(layer, ix, iy, iz, nx, ny, nz) {
        const bounds = layer.bounds || {};
        const xBounds = Array.isArray(bounds.x) ? bounds.x : [-0.5, 0.5];
        const yBounds = Array.isArray(bounds.y) ? bounds.y : [-0.5, 0.5];
        const zBounds = Array.isArray(bounds.z) ? bounds.z : [-0.5, 0.5];
        const xSpan = ovizArFiniteNumber(xBounds[1], 0.5) - ovizArFiniteNumber(xBounds[0], -0.5);
        const ySpan = ovizArFiniteNumber(yBounds[1], 0.5) - ovizArFiniteNumber(yBounds[0], -0.5);
        const zSpan = ovizArFiniteNumber(zBounds[1], 0.5) - ovizArFiniteNumber(zBounds[0], -0.5);
        const cellSizeXPc = Math.max(1e-6, Math.abs(xSpan / Math.max(nx, 1)));
        const cellSizeYPc = Math.max(1e-6, Math.abs(ySpan / Math.max(ny, 1)));
        const cellSizeZPc = Math.max(1e-6, Math.abs(zSpan / Math.max(nz, 1)));
        return {
          x: ovizArFiniteNumber(xBounds[0], -0.5) + ((ix + 0.5) / Math.max(nx, 1)) * xSpan,
          y: ovizArFiniteNumber(yBounds[0], -0.5) + ((iy + 0.5) / Math.max(ny, 1)) * ySpan,
          z: ovizArFiniteNumber(zBounds[0], -0.5) + ((iz + 0.5) / Math.max(nz, 1)) * zSpan,
          cellSizeXPc,
          cellSizeYPc,
          cellSizeZPc,
        };
      }

      function ovizArVolumeSampleSource(layer, scalarData) {
        const proxy = layer && layer.ar_proxy && typeof layer.ar_proxy === "object" ? layer.ar_proxy : null;
        const proxyShape = proxy && proxy.shape ? proxy.shape : {};
        const proxyNx = Math.max(0, Math.round(Number(proxyShape.x) || 0));
        const proxyNy = Math.max(0, Math.round(Number(proxyShape.y) || 0));
        const proxyNz = Math.max(0, Math.round(Number(proxyShape.z) || 0));
        if (
          proxy
          && String(proxy.data_encoding || "uint8") === "uint8"
          && proxy.data_b64
          && proxyNx * proxyNy * proxyNz > 0
          && typeof base64ToUint8Array === "function"
        ) {
          const proxyData = base64ToUint8Array(proxy.data_b64);
          if (proxyData && proxyData.length >= proxyNx * proxyNy * proxyNz) {
            return {
              data: proxyData,
              nx: proxyNx,
              ny: proxyNy,
              nz: proxyNz,
              stride: 1,
              kind: "block-max proxy",
            };
          }
        }
        const shape = layer && layer.shape ? layer.shape : {};
        const nx = Math.max(1, Math.round(Number(shape.x) || 0));
        const ny = Math.max(1, Math.round(Number(shape.y) || 0));
        const nz = Math.max(1, Math.round(Number(shape.z) || 0));
        return {
          data: scalarData,
          nx,
          ny,
          nz,
          stride: Math.max(1, Math.ceil(Math.cbrt((nx * ny * nz) / OVIZ_AR_VOLUME_SCAN_TARGET))),
          kind: "volume fallback",
        };
      }

      async function ovizArCollectVolumeSamples(snapshot, options = {}) {
        const frame = ovizArPresentFrame();
        const layers = ovizArVolumeLayersForSnapshot(frame);
        const maxSamples = Math.max(1, Math.round(Number(options.maxSamples) || OVIZ_AR_MAX_VOLUME_SAMPLES));
        const scanTarget = Math.max(1, Math.round(Number(options.scanTarget) || OVIZ_AR_VOLUME_SCAN_TARGET));
        const mask = (
          typeof activeVolumeLassoSelectionMask === "function"
          && activeVolumeLassoSelectionMask()
        ) || (snapshot && snapshot.hasLassoMask ? ovizArActiveSelectionMask() : null);
        const candidateBuckets = new Map();
        const layerSummaries = [];
        const bucketAxisCount = Math.max(2, Math.ceil(Math.cbrt(maxSamples * 0.85)));

        for (let layerIndex = 0; layerIndex < layers.length; layerIndex += 1) {
          const layer = layers[layerIndex];
          const state = ovizArVolumeStateForLayer(layer);
          if (!state) {
            continue;
          }
          let source = ovizArVolumeSampleSource(layer, null);
          if (!source.data || !source.data.length) {
            const scalarData = await ovizArVolumeScalarArrayFor(layer);
            source = ovizArVolumeSampleSource(layer, scalarData);
          }
          if (!source.data || !source.data.length) {
            continue;
          }
          const nx = source.nx;
          const ny = source.ny;
          const nz = source.nz;
          const totalVoxels = nx * ny * nz;
          if (!(totalVoxels > 0)) {
            continue;
          }
          const colorBytes = ovizArVolumeColorBytesForLayer(layer, state);
          const lutSamples = Math.max(1, Math.floor((colorBytes && colorBytes.length ? colorBytes.length : 0) / 4));
          if (!colorBytes || lutSamples < 1) {
            continue;
          }
          const windowState = ovizArVolumeWindowFor(layer, state);
          const low = Math.min(Math.max(ovizArFiniteNumber(windowState.low, 0.0), 0.0), 1.0);
          const high = Math.min(Math.max(ovizArFiniteNumber(windowState.high, 1.0), 0.0), 1.0);
          const span = Math.max(high - low, 1e-6);
          const stride = source.kind === "block-max proxy"
            ? 1
            : Math.max(source.stride, Math.ceil(Math.cbrt(totalVoxels / scanTarget)));
          const layerName = typeof volumeBaseNameForLayer === "function"
            ? volumeBaseNameForLayer(layer)
            : String(layer.name || layer.key || "Volume");
          let layerSampleCount = 0;

          for (let iz = 0; iz < nz; iz += stride) {
            for (let iy = 0; iy < ny; iy += stride) {
              for (let ix = 0; ix < nx; ix += stride) {
                const dataIndex = (iz * ny * nx) + (iy * nx) + ix;
                const normalizedValue = (Number(source.data[dataIndex]) || 0) / 255.0;
                if (!(normalizedValue > low)) {
                  continue;
                }
                const scaledValue = Math.min(Math.max((normalizedValue - low) / span, 0.0), 1.0);
                const stretchedValue = ovizArVolumeStretchValue(scaledValue, state.stretch);
                if (!(stretchedValue >= OVIZ_AR_VOLUME_MIN_STRETCHED_VALUE)) {
                  continue;
                }
                const position = ovizArVolumeSamplePosition(layer, ix, iy, iz, nx, ny, nz);
                if (mask && !ovizArPointInsideProjectedMask(position, mask)) {
                  continue;
                }
                const colorBin = Math.max(
                  0,
                  Math.min(OVIZ_AR_VOLUME_COLOR_BINS - 1, Math.round(stretchedValue * (OVIZ_AR_VOLUME_COLOR_BINS - 1)))
                );
                const colorIndex = Math.max(
                  0,
                  Math.min(lutSamples - 1, Math.round((colorBin / Math.max(OVIZ_AR_VOLUME_COLOR_BINS - 1, 1)) * (lutSamples - 1)))
                ) * 4;
                const sample = {
                  key: String(layer.key || layerName),
                  label: layerName || "Volume",
                  x: position.x,
                  y: position.y,
                  z: position.z,
                  sizeXPc: position.cellSizeXPc * Math.max(stride * 1.5, (nx / bucketAxisCount) * 1.35),
                  sizeYPc: position.cellSizeYPc * Math.max(stride * 1.5, (ny / bucketAxisCount) * 1.35),
                  sizeZPc: position.cellSizeZPc * Math.max(stride * 1.5, (nz / bucketAxisCount) * 1.35),
                  value: normalizedValue,
                  weight: stretchedValue * Math.max(ovizArFiniteNumber(state.opacity, 0.25), 0.08),
                  color: ovizArHexColorFromBytes(colorBytes, colorIndex),
                  colorBin,
                };
                const bucketX = Math.min(bucketAxisCount - 1, Math.floor((ix / Math.max(nx, 1)) * bucketAxisCount));
                const bucketY = Math.min(bucketAxisCount - 1, Math.floor((iy / Math.max(ny, 1)) * bucketAxisCount));
                const bucketZ = Math.min(bucketAxisCount - 1, Math.floor((iz / Math.max(nz, 1)) * bucketAxisCount));
                const bucketKey = `${layerIndex}:${bucketX}:${bucketY}:${bucketZ}`;
                const bucket = candidateBuckets.get(bucketKey) || [];
                bucket.push(sample);
                bucket.sort((left, right) => ovizArFiniteNumber(right.weight, 0.0) - ovizArFiniteNumber(left.weight, 0.0));
                if (bucket.length > OVIZ_AR_VOLUME_BUCKET_SAMPLES) {
                  bucket.length = OVIZ_AR_VOLUME_BUCKET_SAMPLES;
                }
                candidateBuckets.set(bucketKey, bucket);
                layerSampleCount += 1;
              }
            }
          }

          if (layerSampleCount > 0) {
            layerSummaries.push({
              key: String(layer.key || ""),
              stateKey: typeof volumeStateKeyForLayer === "function" ? volumeStateKeyForLayer(layer) : String(layer.state_key || layer.key || ""),
              label: layerName || "Volume",
              samples: layerSampleCount,
              source: source.kind,
            });
          }
        }

        const samples = [];
        for (let rank = 0; rank < OVIZ_AR_VOLUME_BUCKET_SAMPLES; rank += 1) {
          const rankedSamples = [];
          candidateBuckets.forEach((bucket) => {
            if (bucket[rank]) {
              rankedSamples.push(bucket[rank]);
            }
          });
          rankedSamples.sort((left, right) => ovizArFiniteNumber(right.weight, 0.0) - ovizArFiniteNumber(left.weight, 0.0));
          rankedSamples.forEach((sample) => samples.push(sample));
        }
        const cappedSamples = samples.slice(0, maxSamples);
        return {
          samples: cappedSamples,
          layers: layerSummaries,
          truncated: samples.length > cappedSamples.length,
          totalCandidates: samples.length,
        };
      }

      function ovizArCreateVolumeCloudTexture() {
        const canvasEl = document.createElement("canvas");
        canvasEl.width = 64;
        canvasEl.height = 64;
        const context = canvasEl.getContext("2d", { alpha: true });
        if (!context) {
          throw new Error("AR volume texture canvas is unavailable.");
        }
        context.clearRect(0, 0, canvasEl.width, canvasEl.height);
        const gradient = context.createRadialGradient(32, 32, 0, 32, 32, 31.5);
        gradient.addColorStop(0.0, "rgba(255, 255, 255, 1.0)");
        gradient.addColorStop(0.34, "rgba(255, 255, 255, 0.82)");
        gradient.addColorStop(0.68, "rgba(255, 255, 255, 0.26)");
        gradient.addColorStop(1.0, "rgba(255, 255, 255, 0.0)");
        context.fillStyle = gradient;
        context.fillRect(0, 0, canvasEl.width, canvasEl.height);
        const texture = new THREE.CanvasTexture(canvasEl);
        texture.minFilter = THREE.LinearFilter;
        texture.magFilter = THREE.LinearFilter;
        texture.generateMipmaps = false;
        texture.needsUpdate = true;
        return texture;
      }

      function ovizArAppendCloudQuad(vertices, normals, uvs, indices, corners, normal) {
        const baseIndex = vertices.length / 3;
        corners.forEach((corner) => {
          vertices.push(corner[0], corner[1], corner[2]);
          normals.push(normal[0], normal[1], normal[2]);
        });
        uvs.push(0, 0, 1, 0, 1, 1, 0, 1);
        indices.push(
          baseIndex, baseIndex + 1, baseIndex + 2,
          baseIndex, baseIndex + 2, baseIndex + 3
        );
      }

      function ovizArAppendCloudletGeometry(vertices, normals, uvs, indices, center, size) {
        const cx = ovizArFiniteNumber(center && center.x, 0.0);
        const cy = ovizArFiniteNumber(center && center.y, 0.0);
        const cz = ovizArFiniteNumber(center && center.z, 0.0);
        const halfX = Math.max(ovizArFiniteNumber(size && size.x, 0.01), 1e-5) * 0.5;
        const halfY = Math.max(ovizArFiniteNumber(size && size.y, 0.01), 1e-5) * 0.5;
        const halfZ = Math.max(ovizArFiniteNumber(size && size.z, 0.01), 1e-5) * 0.5;
        ovizArAppendCloudQuad(vertices, normals, uvs, indices, [
          [cx - halfX, cy - halfY, cz],
          [cx + halfX, cy - halfY, cz],
          [cx + halfX, cy + halfY, cz],
          [cx - halfX, cy + halfY, cz],
        ], [0, 0, 1]);
        ovizArAppendCloudQuad(vertices, normals, uvs, indices, [
          [cx - halfX, cy, cz - halfZ],
          [cx + halfX, cy, cz - halfZ],
          [cx + halfX, cy, cz + halfZ],
          [cx - halfX, cy, cz + halfZ],
        ], [0, 1, 0]);
        ovizArAppendCloudQuad(vertices, normals, uvs, indices, [
          [cx, cy - halfY, cz - halfZ],
          [cx, cy - halfY, cz + halfZ],
          [cx, cy + halfY, cz + halfZ],
          [cx, cy + halfY, cz - halfZ],
        ], [1, 0, 0]);
      }

      function ovizArAddVolumeProxy(group, volumeSamples, transform) {
        const samples = Array.isArray(volumeSamples) ? volumeSamples : [];
        if (!samples.length) {
          return;
        }
        const cloudTexture = ovizArCreateVolumeCloudTexture();
        const samplesByColor = new Map();
        samples.forEach((sample) => {
          const colorBin = Math.max(0, Math.min(OVIZ_AR_VOLUME_COLOR_BINS - 1, Math.round(Number(sample.colorBin) || 0)));
          if (!samplesByColor.has(colorBin)) {
            samplesByColor.set(colorBin, []);
          }
          samplesByColor.get(colorBin).push(sample);
        });
        samplesByColor.forEach((colorSamples) => {
          const vertices = [];
          const normals = [];
          const uvs = [];
          const indices = [];
          colorSamples.forEach((sample) => {
            const center = ovizArVectorFromPoint(sample, transform);
            const scale = ovizArFiniteNumber(transform && transform.scale, 1.0);
            const size = {
              x: Math.max(0.008, Math.min(0.18, ovizArFiniteNumber(sample.sizeXPc, 1.0) * scale)),
              y: Math.max(0.006, Math.min(0.12, ovizArFiniteNumber(sample.sizeZPc, 1.0) * scale)),
              z: Math.max(0.008, Math.min(0.18, ovizArFiniteNumber(sample.sizeYPc, 1.0) * scale)),
            };
            ovizArAppendCloudletGeometry(vertices, normals, uvs, indices, center, size);
          });
          if (!vertices.length || !indices.length) {
            return;
          }
          const geometry = new THREE.BufferGeometry();
          geometry.setAttribute("position", new THREE.Float32BufferAttribute(vertices, 3));
          geometry.setAttribute("normal", new THREE.Float32BufferAttribute(normals, 3));
          geometry.setAttribute("uv", new THREE.Float32BufferAttribute(uvs, 2));
          geometry.setIndex(indices);
          geometry.userData = geometry.userData || {};
          geometry.userData.ovizArDoubleSided = true;
          const representative = colorSamples.reduce((best, sample) => (
            !best || ovizArFiniteNumber(sample.weight, 0.0) > ovizArFiniteNumber(best.weight, 0.0) ? sample : best
          ), null);
          const color = String((representative && representative.color) || "#ffffff");
          const meanWeight = colorSamples.reduce((sum, sample) => sum + ovizArFiniteNumber(sample.weight, 0.0), 0.0)
            / Math.max(colorSamples.length, 1);
          const material = new THREE.MeshStandardMaterial({
            map: cloudTexture,
            color: ovizArColor(color, "#ffffff"),
            emissive: ovizArColor(color, "#ffffff"),
            emissiveIntensity: 0.20,
            roughness: 0.90,
            metalness: 0.0,
            transparent: true,
            depthWrite: false,
            side: THREE.DoubleSide,
            opacity: Math.max(0.06, Math.min(0.24, 0.045 + (0.18 * Math.sqrt(Math.max(meanWeight, 0.0))))),
          });
          const mesh = new THREE.Mesh(geometry, material);
          mesh.userData.ovizArVolumeProxy = true;
          group.add(mesh);
        });
      }

      async function buildOvizAr3DScene(snapshot) {
        const volumeResult = await ovizArCollectVolumeSamples(snapshot);
        const volumeLayers = ovizArVolumeLayersForSnapshot(ovizArPresentFrame());
        const sceneAr = new THREE.Scene();
        sceneAr.add(new THREE.AmbientLight(0xffffff, 0.9));
        const light = new THREE.DirectionalLight(0xffffff, 0.7);
        light.position.set(1.2, 2.0, 1.4);
        sceneAr.add(light);
        const group = new THREE.Group();
        sceneAr.add(group);
        const transform = ovizArSceneTransform(snapshot.points || [], 0.72, volumeLayers);
        const sphereGeometry = new THREE.SphereGeometry(0.0085, 16, 10);
        const markerMaterials = new Map();
        ovizArAddVolumeProxy(group, volumeResult.samples, transform);
        (snapshot.points || []).forEach((point) => {
          const color = ovizArColor(point.color, "#ffffff");
          const materialKey = String(point.color || "#ffffff").toLowerCase();
          if (!markerMaterials.has(materialKey)) {
            markerMaterials.set(materialKey, new THREE.MeshStandardMaterial({
              color,
              emissive: color,
              emissiveIntensity: 0.35,
              roughness: 0.45,
              metalness: 0.0,
            }));
          }
          const material = markerMaterials.get(materialKey);
          const marker = new THREE.Mesh(sphereGeometry, material);
          const position = ovizArVectorFromPoint(point, transform);
          const nStars = ovizArFiniteNumber(point.nStars, NaN);
          const sizeScale = Number.isFinite(nStars)
            ? Math.max(0.72, Math.min(1.35, Math.sqrt(Math.max(nStars, 1.0)) / 12.0))
            : 1.0;
          marker.scale.setScalar(sizeScale);
          marker.position.copy(position);
          marker.userData.ovizArKey = point.key;
          group.add(marker);
        });
        sceneAr.updateMatrixWorld(true);
        sceneAr.userData.ovizArSummary = {
          mode: "3d",
          pointCount: snapshot.points.length,
          volumeSampleCount: volumeResult.samples.length,
          volumeLayerCount: volumeResult.layers.length,
          volumeTruncated: volumeResult.truncated,
          presentTimeMyr: snapshot.presentTimeMyr,
          scalePcToMeters: transform.scale,
          volumeRepresentation: "soft-gaussian-cloudlets",
          volumeSources: volumeResult.layers.map((layer) => layer.source).filter(Boolean),
        };
        return sceneAr;
      }

      function ovizArSkyDirectionForLonLatDeg(lonDeg, latDeg, radius = 1.0) {
        const lon = normalizeSkyLongitude(lonDeg) * Math.PI / 180.0;
        const lat = Math.max(-90.0, Math.min(90.0, Number(latDeg))) * Math.PI / 180.0;
        if (!Number.isFinite(lon) || !Number.isFinite(lat)) {
          return null;
        }
        const theta = (Math.PI / 2.0) - lat;
        const sinTheta = Math.sin(theta);
        return {
          x: -radius * Math.cos(lon) * sinTheta,
          y: radius * Math.cos(theta),
          z: radius * Math.sin(lon) * sinTheta,
        };
      }

      function ovizArSkyDirectionForPoint(point, radius = 1.0) {
        const coordsys = typeof skyDomeHips2FitsCoordsys === "function" ? skyDomeHips2FitsCoordsys() : "galactic";
        if (coordsys === "icrs") {
          let ra = ovizArFiniteNumber(point && point.ra);
          let dec = ovizArFiniteNumber(point && point.dec);
          if (!Number.isFinite(ra) || !Number.isFinite(dec)) {
            const icrs = icrsDegFromGalacticDeg(point && point.l, point && point.b);
            ra = ovizArFiniteNumber(icrs && icrs.ra);
            dec = ovizArFiniteNumber(icrs && icrs.dec);
          }
          return ovizArSkyDirectionForLonLatDeg(ra, dec, radius);
        }
        let l = ovizArFiniteNumber(point && point.l);
        let b = ovizArFiniteNumber(point && point.b);
        if (!Number.isFinite(l) || !Number.isFinite(b)) {
          const galactic = galacticDegFromIcrsDeg(point && point.ra, point && point.dec);
          l = ovizArFiniteNumber(galactic && galactic.l);
          b = ovizArFiniteNumber(galactic && galactic.b);
        }
        return ovizArSkyDirectionForLonLatDeg(l, b, radius);
      }

      function ovizArTextureLoaderLoad(url) {
        return new Promise((resolve, reject) => {
          const loader = new THREE.TextureLoader();
          loader.setCrossOrigin("anonymous");
          loader.load(
            url,
            (texture) => {
              texture.minFilter = THREE.LinearFilter;
              texture.magFilter = THREE.LinearFilter;
              texture.generateMipmaps = false;
              if (THREE.sRGBEncoding != null) {
                texture.encoding = THREE.sRGBEncoding;
              }
              resolve(texture);
            },
            undefined,
            reject
          );
        });
      }

      function ovizArSkyLayerForExport() {
        if (typeof activeSkyLayer === "function") {
          const layer = activeSkyLayer();
          if (layer && layer.visible !== false && layer.survey) {
            return layer;
          }
        }
        if (typeof visibleSkyLayers === "function") {
          const layers = visibleSkyLayers();
          if (Array.isArray(layers) && layers.length && layers[0].survey) {
            return layers[0];
          }
        }
        return null;
      }

      function ovizArSkyTextureUrl(width, height) {
        const baseUrl = skyDomeHips2FitsUrl(width, height);
        const layer = ovizArSkyLayerForExport();
        if (!layer || !layer.survey) {
          return baseUrl;
        }
        const url = new URL(baseUrl);
        url.searchParams.set("hips", String(layer.survey));
        if (layer.stretch) {
          url.searchParams.set("stretch", String(layer.stretch));
        }
        if (layer.cutMin !== "" && layer.cutMin != null) {
          url.searchParams.set("min_cut", String(layer.cutMin));
        } else {
          url.searchParams.delete("min_cut");
        }
        if (layer.cutMax !== "" && layer.cutMax != null) {
          url.searchParams.set("max_cut", String(layer.cutMax));
        } else {
          url.searchParams.delete("max_cut");
        }
        return url.href;
      }

      async function ovizArLoadSkyTexture(target) {
        const width = Math.round(Number(target && target.width) || OVIZ_AR_SKY_TEXTURE_LOW.width);
        const height = Math.round(Number(target && target.height) || OVIZ_AR_SKY_TEXTURE_LOW.height);
        const url = ovizArSkyTextureUrl(width, height);
        const texture = await ovizArTextureLoaderLoad(url);
        const layer = ovizArSkyLayerForExport();
        return { texture, width, height, url, survey: layer && layer.survey ? String(layer.survey) : "" };
      }

      function ovizArInvertGeometryFaces(geometry) {
        if (!geometry) {
          return geometry;
        }
        const index = geometry.index;
        if (index && index.array) {
          for (let offset = 0; offset + 2 < index.array.length; offset += 3) {
            const swap = index.array[offset + 1];
            index.array[offset + 1] = index.array[offset + 2];
            index.array[offset + 2] = swap;
          }
          index.needsUpdate = true;
        }
        const normal = geometry.attributes && geometry.attributes.normal;
        if (normal) {
          for (let offset = 0; offset < normal.count; offset += 1) {
            normal.setXYZ(offset, -normal.getX(offset), -normal.getY(offset), -normal.getZ(offset));
          }
          normal.needsUpdate = true;
        }
        return geometry;
      }

      async function buildOvizArSkyDomeScene(snapshot, target = OVIZ_AR_SKY_TEXTURE_HIGH) {
        const loaded = await ovizArLoadSkyTexture(target);
        const radius = 1.45;
        const centerY = radius + 0.04;
        const sceneAr = new THREE.Scene();
        sceneAr.add(new THREE.AmbientLight(0xffffff, 1.0));
        loaded.texture.userData = loaded.texture.userData || {};
        loaded.texture.userData.ovizArTextureFormat = "jpeg";
        const domeGeometry = ovizArInvertGeometryFaces(new THREE.SphereGeometry(radius, 96, 48));
        const domeMaterial = new THREE.MeshStandardMaterial({
          map: loaded.texture,
          roughness: 1.0,
          metalness: 0.0,
        });
        const dome = new THREE.Mesh(domeGeometry, domeMaterial);
        dome.position.y = centerY;
        sceneAr.add(dome);
        const markerGeometry = new THREE.SphereGeometry(0.010, 16, 10);
        const markerMaterials = new Map();
        (snapshot.points || []).forEach((point) => {
          const direction = ovizArSkyDirectionForPoint(point, radius * 0.965);
          if (!direction) {
            return;
          }
          const color = ovizArColor(point.color, "#ffffff");
          const materialKey = String(point.color || "#ffffff").toLowerCase();
          if (!markerMaterials.has(materialKey)) {
            markerMaterials.set(materialKey, new THREE.MeshStandardMaterial({
              color,
              emissive: color,
              emissiveIntensity: 0.65,
              roughness: 0.5,
              metalness: 0.0,
            }));
          }
          const material = markerMaterials.get(materialKey);
          const marker = new THREE.Mesh(markerGeometry, material);
          marker.position.set(direction.x, direction.y + centerY, direction.z);
          marker.userData.ovizArKey = point.key;
          sceneAr.add(marker);
        });
        sceneAr.updateMatrixWorld(true);
        sceneAr.userData.ovizArSummary = {
          mode: "sky",
          pointCount: snapshot.points.length,
          textureWidth: loaded.width,
          textureHeight: loaded.height,
          textureUrl: loaded.url,
          survey: loaded.survey,
          coordsys: typeof skyDomeHips2FitsCoordsys === "function" ? skyDomeHips2FitsCoordsys() : "galactic",
        };
        return sceneAr;
      }

      function ovizArFormatByteSize(byteCount) {
        const bytes = Math.max(0, Math.round(Number(byteCount) || 0));
        if (bytes < 1024) {
          return `${bytes} B`;
        }
        if (bytes < 1024 * 1024) {
          return `${(bytes / 1024).toFixed(1)} KB`;
        }
        return `${(bytes / (1024 * 1024)).toFixed(1)} MB`;
      }

      async function ovizArArrayBufferFromExporterResult(result) {
        if (result instanceof ArrayBuffer) {
          return result;
        }
        if (typeof ArrayBuffer !== "undefined" && ArrayBuffer.isView && ArrayBuffer.isView(result)) {
          return result.buffer.slice(result.byteOffset, result.byteOffset + result.byteLength);
        }
        if (result && result.buffer instanceof ArrayBuffer) {
          return result.buffer;
        }
        if (typeof Blob !== "undefined" && result instanceof Blob) {
          return await result.arrayBuffer();
        }
        throw new Error("USDZ exporter returned an unsupported payload.");
      }

      function ovizArValidateUsdZArrayBuffer(arrayBuffer) {
        const byteLength = Number(arrayBuffer && arrayBuffer.byteLength) || 0;
        if (!(arrayBuffer instanceof ArrayBuffer) || byteLength < OVIZ_AR_MIN_USDZ_BYTES) {
          throw new Error(`USDZ exporter produced an empty or incomplete package (${ovizArFormatByteSize(byteLength)}).`);
        }
        return arrayBuffer;
      }

      async function ovizArParseUsdZ(sceneAr) {
        const exporter = new OvizUSDZExporter();
        if (sceneAr && typeof sceneAr.updateMatrixWorld === "function") {
          sceneAr.updateMatrixWorld(true);
        }
        const result = await exporter.parse(sceneAr, {
          maxTextureSize: OVIZ_AR_SKY_TEXTURE_HIGH.width,
          jpegQuality: 0.88,
        });
        return ovizArValidateUsdZArrayBuffer(await ovizArArrayBufferFromExporterResult(result));
      }

      function ovizArClearObjectUrl() {
        if (ovizArObjectUrl) {
          URL.revokeObjectURL(ovizArObjectUrl);
          ovizArObjectUrl = "";
        }
        ovizArQuickLookUrl = "";
      }

      function ovizArQuickLookAvailable() {
        const nav = window.navigator || {};
        const userAgent = String(nav.userAgent || "");
        const platform = String(nav.platform || "");
        return /iPad|iPhone|iPod/i.test(userAgent)
          || (platform === "MacIntel" && Number(nav.maxTouchPoints) > 1);
      }

      function ovizArHostedQuickLookAvailable() {
        return Boolean(
          typeof window !== "undefined"
          && window.isSecureContext
          && window.location
          && /^https?:$/i.test(String(window.location.protocol || ""))
          && typeof navigator !== "undefined"
          && navigator.serviceWorker
          && window.caches
        );
      }

      function ovizArQuickLookWorkerUrl() {
        const url = new URL("oviz-ar-quicklook-sw.js", window.location.href);
        url.searchParams.set("v", OVIZ_AR_QUICKLOOK_WORKER_VERSION);
        return url.href;
      }

      function ovizArQuickLookScopeUrl() {
        return new URL("./", window.location.href).href;
      }

      function ovizArQuickLookAssetUrl(filename) {
        const safeFilename = String(filename || "oviz_ar_snapshot.usdz")
          .replace(/[^a-z0-9_.-]+/gi, "_")
          .replace(/^_+/, "")
          || "oviz_ar_snapshot.usdz";
        const cryptoObj = typeof window !== "undefined" ? window.crypto : null;
        const token = cryptoObj && typeof cryptoObj.randomUUID === "function"
          ? cryptoObj.randomUUID()
          : `${Date.now()}_${Math.random().toString(36).slice(2)}`;
        return new URL(`oviz-ar-quicklook/${token}/${safeFilename}`, window.location.href).href;
      }

      function ovizArWaitForServiceWorkerControl(timeoutMs = 2500) {
        if (!navigator.serviceWorker || navigator.serviceWorker.controller) {
          return Promise.resolve(Boolean(navigator.serviceWorker && navigator.serviceWorker.controller));
        }
        return new Promise((resolve) => {
          let settled = false;
          const finish = (value) => {
            if (settled) {
              return;
            }
            settled = true;
            window.clearTimeout(timer);
            navigator.serviceWorker.removeEventListener("controllerchange", onControllerChange);
            resolve(value);
          };
          const onControllerChange = () => finish(true);
          const timer = window.setTimeout(() => finish(Boolean(navigator.serviceWorker.controller)), timeoutMs);
          navigator.serviceWorker.addEventListener("controllerchange", onControllerChange);
        });
      }

      function ovizArWithTimeout(promise, timeoutMs, message) {
        return new Promise((resolve, reject) => {
          let settled = false;
          const finish = (callback, value) => {
            if (settled) {
              return;
            }
            settled = true;
            window.clearTimeout(timer);
            callback(value);
          };
          const timer = window.setTimeout(
            () => finish(reject, new Error(message || "Timed out.")),
            Math.max(1, Number(timeoutMs) || 1)
          );
          Promise.resolve(promise).then(
            (value) => finish(resolve, value),
            (err) => finish(reject, err)
          );
        });
      }

      async function ovizArPrepareQuickLookAsset(blob, filename) {
        if (!ovizArHostedQuickLookAvailable()) {
          throw new Error("Hosted Quick Look handoff is unavailable.");
        }
        const byteLength = Number(blob && blob.size) || 0;
        if (byteLength < OVIZ_AR_MIN_USDZ_BYTES) {
          throw new Error(`USDZ package is too small for Quick Look (${ovizArFormatByteSize(byteLength)}).`);
        }
        const registration = await ovizArWithTimeout(
          navigator.serviceWorker.register(
            ovizArQuickLookWorkerUrl(),
            { scope: ovizArQuickLookScopeUrl() }
          ),
          2500,
          "Quick Look handoff worker did not register in time."
        );
        if (registration && typeof registration.update === "function") {
          try {
            await registration.update();
          } catch (_err) {
            // A stale but active worker is still usable for this short-lived cache handoff.
          }
        }
        try {
          await ovizArWithTimeout(
            navigator.serviceWorker.ready,
            2500,
            "Quick Look handoff worker did not become ready in time."
          );
        } catch (err) {
          if (!registration || !registration.active) {
            throw err;
          }
        }
        const assetUrl = ovizArQuickLookAssetUrl(filename);
        const cache = await ovizArWithTimeout(
          window.caches.open(OVIZ_AR_QUICKLOOK_CACHE),
          2000,
          "Quick Look cache did not open in time."
        );
        const staleRequests = await cache.keys();
        await Promise.all(staleRequests.map((request) => cache.delete(request)));
        const safeFilename = String(filename || "oviz_ar_snapshot.usdz").replace(/"/g, "");
        await ovizArWithTimeout(
          cache.put(assetUrl, new Response(blob, {
            headers: {
              "Content-Type": "model/vnd.usdz+zip",
              "Content-Disposition": `inline; filename="${safeFilename}"`,
              "Cache-Control": "no-store",
              "Accept-Ranges": "bytes",
              "X-Oviz-AR-Bytes": String(byteLength),
            },
          })),
          3000,
          "Quick Look package did not stage in time."
        );
        const controlled = await ovizArWaitForServiceWorkerControl(2500);
        if (!controlled && !navigator.serviceWorker.controller) {
          throw new Error("Quick Look handoff is not ready yet. Close and reopen the AR chooser once.");
        }
        return assetUrl;
      }

      function ovizArUsdZFilename(mode) {
        const cleanMode = String(mode || "scene").replace(/[^a-z0-9_-]+/gi, "_").toLowerCase();
        const count = collectOvizArSnapshot(mode).points.length;
        return `oviz_${cleanMode}_${count}_clusters.usdz`;
      }

      function ovizArSetStatus(message, tone = "neutral") {
        const dialog = ensureOvizArDialog();
        const statusEl = dialog.querySelector(".oviz-three-ar-status");
        if (statusEl) {
          statusEl.textContent = String(message || "");
          statusEl.dataset.tone = tone;
        }
      }

      function ovizArSetBusy(isBusy) {
        ovizArExportBusy = Boolean(isBusy);
        const dialog = ensureOvizArDialog();
        dialog.querySelectorAll(".oviz-three-ar-mode").forEach((button) => {
          button.disabled = ovizArExportBusy;
        });
        renderArSnapshotButtonState();
      }

      async function ovizArSetReady(blob, filename, message) {
        const byteLength = Number(blob && blob.size) || 0;
        if (byteLength < OVIZ_AR_MIN_USDZ_BYTES) {
          throw new Error(`USDZ exporter produced an empty or incomplete package (${ovizArFormatByteSize(byteLength)}).`);
        }
        ovizArClearObjectUrl();
        ovizArObjectUrl = URL.createObjectURL(blob);
        let quickLookError = "";
        if (ovizArQuickLookAvailable()) {
          try {
            ovizArQuickLookUrl = await ovizArPrepareQuickLookAsset(blob, filename);
          } catch (err) {
            quickLookError = err && err.message ? String(err.message) : "direct Quick Look handoff failed";
            ovizArQuickLookUrl = "";
          }
        }
        const dialog = ensureOvizArDialog();
        const hiddenLink = dialog.querySelector(".oviz-three-ar-hidden-link");
        const openButton = dialog.querySelector(".oviz-three-ar-open-link");
        const downloadLink = dialog.querySelector(".oviz-three-ar-download");
        const quickLookHref = ovizArQuickLookUrl || ovizArObjectUrl;
        if (hiddenLink) {
          hiddenLink.href = quickLookHref;
          hiddenLink.removeAttribute("download");
        }
        if (openButton) {
          openButton.hidden = false;
          openButton.disabled = false;
        }
        if (downloadLink) {
          downloadLink.hidden = false;
          downloadLink.href = ovizArObjectUrl;
          downloadLink.download = filename;
        }
        const sizeText = ovizArFormatByteSize(byteLength);
        const baseMessage = message || (ovizArQuickLookAvailable() ? "AR snapshot is ready." : "USDZ snapshot is ready.");
        const canOpenQuickLook = Boolean(ovizArQuickLookAvailable() && ovizArQuickLookUrl);
        let statusMessage = `${baseMessage} (${sizeText}). Tap Open AR Snapshot.`;
        let tone = canOpenQuickLook ? "ready" : "warn";
        if (ovizArQuickLookAvailable() && !ovizArQuickLookUrl) {
          statusMessage = `${baseMessage} (${sizeText}). Tap Open AR Snapshot to use the direct AR fallback.`;
          if (quickLookError) {
            statusMessage += ` ${quickLookError}`;
          }
        } else if (!ovizArQuickLookAvailable()) {
          statusMessage = `${baseMessage} (${sizeText}). Download the USDZ on this browser.`;
        }
        ovizArSetStatus(
          statusMessage,
          tone
        );
      }

      async function ovizArBuildAndExport(mode) {
        const snapshot = collectOvizArSnapshot(mode);
        if (mode === "sky" && !snapshot.points.length) {
          ovizArSetStatus("No AR-ready clusters were found in the t=0 Myr frame.", "warn");
          return;
        }
        if (mode !== "sky" && !snapshot.points.length && !ovizArHasVolumeSnapshot()) {
          ovizArSetStatus("No AR-ready clusters or visible volume were found at t=0 Myr.", "warn");
          return;
        }
        ovizArSetBusy(true);
        ovizArSetStatus(mode === "sky" ? "Generating AR snapshot..." : "Sampling AR scene...", "neutral");
        const filename = ovizArUsdZFilename(mode);
        try {
          let sceneAr = null;
          let arrayBuffer = null;
          if (mode === "sky") {
            const primarySkyTarget = ovizArQuickLookAvailable()
              ? OVIZ_AR_SKY_TEXTURE_LOW
              : OVIZ_AR_SKY_TEXTURE_HIGH;
            const fallbackSkyTarget = ovizArQuickLookAvailable()
              ? { width: 1024, height: 512 }
              : OVIZ_AR_SKY_TEXTURE_LOW;
            try {
              sceneAr = await buildOvizArSkyDomeScene(snapshot, primarySkyTarget);
              arrayBuffer = await ovizArParseUsdZ(sceneAr);
            } catch (_err) {
              sceneAr = await buildOvizArSkyDomeScene(snapshot, fallbackSkyTarget);
              arrayBuffer = await ovizArParseUsdZ(sceneAr);
            }
          } else {
            sceneAr = await buildOvizAr3DScene(snapshot);
            if (!snapshot.points.length && !(Number((sceneAr.userData.ovizArSummary || {}).volumeSampleCount) > 0)) {
              throw new Error("No volume samples survived AR downsampling.");
            }
            arrayBuffer = await ovizArParseUsdZ(sceneAr);
          }
          const blob = new Blob([arrayBuffer], { type: "model/vnd.usdz+zip" });
          const pointText = `${snapshot.points.length} cluster${snapshot.points.length === 1 ? "" : "s"}`;
          const summary = sceneAr && sceneAr.userData ? sceneAr.userData.ovizArSummary || {} : {};
          const volumeCount = Math.round(Number(summary.volumeSampleCount) || 0);
          const volumeText = volumeCount ? ` plus ${volumeCount} volume sample${volumeCount === 1 ? "" : "s"}` : "";
          const selectionText = snapshot.selectionMode === "selection" ? "selected " : "";
          const suffix = (snapshot.truncated || summary.volumeTruncated) ? " AR limits applied." : "";
          await ovizArSetReady(blob, filename, `${selectionText}${pointText}${volumeText} exported at t=0 Myr.${suffix}`);
        } catch (err) {
          const message = err && err.message ? String(err.message) : "AR export failed.";
          ovizArSetStatus(`AR export failed: ${message}`, "error");
        } finally {
          ovizArSetBusy(false);
        }
      }

      function ovizArResetDialogActions() {
        ovizArQuickLookUrl = "";
        const dialog = ensureOvizArDialog();
        const openButton = dialog.querySelector(".oviz-three-ar-open-link");
        const downloadLink = dialog.querySelector(".oviz-three-ar-download");
        const hiddenLink = dialog.querySelector(".oviz-three-ar-hidden-link");
        if (openButton) {
          openButton.hidden = true;
          openButton.disabled = true;
        }
        if (downloadLink) {
          downloadLink.hidden = true;
          downloadLink.removeAttribute("href");
          downloadLink.removeAttribute("download");
        }
        if (hiddenLink) {
          hiddenLink.removeAttribute("href");
          hiddenLink.removeAttribute("download");
        }
      }

      function ensureOvizArDialog() {
        if (ovizArDialogEl) {
          return ovizArDialogEl;
        }
        const dialog = document.createElement("div");
        dialog.className = "oviz-three-ar-dialog";
        dialog.dataset.open = "false";
        dialog.setAttribute("aria-hidden", "true");
        dialog.innerHTML = `
          <div class="oviz-three-ar-backdrop"></div>
          <div class="oviz-three-ar-card" role="dialog" aria-modal="true" aria-label="AR Snapshot">
            <div class="oviz-three-ar-heading">AR Snapshot</div>
            <div class="oviz-three-ar-actions">
              <button class="oviz-three-ar-mode" type="button" data-mode="3d">3D Scene</button>
              <button class="oviz-three-ar-mode" type="button" data-mode="sky">Sky Dome</button>
            </div>
            <div class="oviz-three-ar-status" data-tone="neutral"></div>
            <div class="oviz-three-ar-ready-actions">
              <button class="oviz-three-ar-open-link" type="button" hidden disabled>Open AR Snapshot</button>
              <a class="oviz-three-ar-download" hidden>Download USDZ</a>
              <a class="oviz-three-ar-hidden-link" rel="ar" hidden><img alt="" /></a>
            </div>
            <button class="oviz-three-ar-close" type="button">Close</button>
          </div>
        `;
        root.appendChild(dialog);
        dialog.querySelector(".oviz-three-ar-backdrop").addEventListener("click", closeOvizArSnapshotChooser);
        dialog.querySelector(".oviz-three-ar-close").addEventListener("click", closeOvizArSnapshotChooser);
        dialog.querySelectorAll(".oviz-three-ar-mode").forEach((button) => {
          button.addEventListener("click", () => {
            ovizArResetDialogActions();
            ovizArBuildAndExport(button.dataset.mode === "sky" ? "sky" : "3d");
          });
        });
        const openButton = dialog.querySelector(".oviz-three-ar-open-link");
        const hiddenLink = dialog.querySelector(".oviz-three-ar-hidden-link");
        if (openButton && hiddenLink) {
          openButton.addEventListener("click", () => hiddenLink.click());
        }
        ovizArDialogEl = dialog;
        return ovizArDialogEl;
      }

      function openOvizArSnapshotChooser() {
        const dialog = ensureOvizArDialog();
        ovizArResetDialogActions();
        if (!ovizArCanExportSelection()) {
          ovizArSetStatus("AR export is unavailable for this view.", "warn");
        } else {
          const selectedCount = ovizArSelectedClusterKeys().size;
          ovizArSetStatus(
            selectedCount
              ? `${selectedCount} selected cluster${selectedCount === 1 ? "" : "s"}.`
              : `No selection; exporting a capped t=0 Myr scene snapshot.`,
            selectedCount ? "neutral" : "warn"
          );
        }
        dialog.dataset.open = "true";
        dialog.setAttribute("aria-hidden", "false");
        if (mobileArButtonEl) {
          mobileArButtonEl.dataset.active = "true";
          mobileArButtonEl.setAttribute("aria-expanded", "true");
        }
      }

      function closeOvizArSnapshotChooser() {
        if (!ovizArDialogEl) {
          return;
        }
        ovizArDialogEl.dataset.open = "false";
        ovizArDialogEl.setAttribute("aria-hidden", "true");
        renderArSnapshotButtonState();
        focusViewer();
      }

      function ovizArWarmQuickLookWorker() {
        if (!ovizArQuickLookAvailable() || !ovizArHostedQuickLookAvailable()) {
          return;
        }
        navigator.serviceWorker.register(
          ovizArQuickLookWorkerUrl(),
          { scope: ovizArQuickLookScopeUrl() }
        ).catch(() => null);
      }

      if (typeof window !== "undefined") {
        window.setTimeout(ovizArWarmQuickLookWorker, 0);
      }
""".strip()


THREEJS_AR_QUICKLOOK_SERVICE_WORKER_JS = r"""
const OVIZ_AR_QUICKLOOK_CACHE = "oviz-ar-quicklook-v2";

async function ovizArRangeResponse(request, response) {
  const range = request.headers.get("range");
  if (!range) {
    return response;
  }
  const match = /^bytes=(\d*)-(\d*)$/i.exec(range);
  if (!match) {
    return response;
  }
  const buffer = await response.arrayBuffer();
  const size = buffer.byteLength;
  let start = match[1] ? Number.parseInt(match[1], 10) : 0;
  let end = match[2] ? Number.parseInt(match[2], 10) : size - 1;
  if (!match[1] && match[2]) {
    const suffixLength = Math.max(0, Number.parseInt(match[2], 10) || 0);
    start = Math.max(size - suffixLength, 0);
    end = size - 1;
  }
  start = Math.max(0, Math.min(start, size - 1));
  end = Math.max(start, Math.min(end, size - 1));
  const headers = new Headers(response.headers);
  headers.set("Content-Range", `bytes ${start}-${end}/${size}`);
  headers.set("Accept-Ranges", "bytes");
  headers.set("Content-Length", String(end - start + 1));
  return new Response(buffer.slice(start, end + 1), {
    status: 206,
    statusText: "Partial Content",
    headers,
  });
}

self.addEventListener("install", (event) => {
  event.waitUntil(self.skipWaiting());
});

self.addEventListener("activate", (event) => {
  event.waitUntil(self.clients.claim());
});

self.addEventListener("fetch", (event) => {
  const url = new URL(event.request.url);
  if (!url.pathname.includes("/oviz-ar-quicklook/")) {
    return;
  }
  event.respondWith((async () => {
    const cache = await caches.open(OVIZ_AR_QUICKLOOK_CACHE);
    const cached = await cache.match(event.request, { ignoreSearch: false });
    if (cached) {
      return ovizArRangeResponse(event.request, cached);
    }
    return new Response("Missing Oviz AR snapshot.", {
      status: 404,
      headers: {
        "Content-Type": "text/plain; charset=utf-8",
        "Cache-Control": "no-store",
      },
    });
  })());
});
""".strip()
