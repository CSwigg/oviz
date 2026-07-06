THREEJS_AR_RUNTIME_JS = r"""
      const OVIZ_AR_MAX_POINTS = 320;
      const OVIZ_AR_MAX_LABELS = 36;
      const OVIZ_AR_MAX_TRAIL_CLUSTERS = 48;
      const OVIZ_AR_MAX_TRAIL_POINTS_PER_CLUSTER = 32;
      const OVIZ_AR_MAX_VOLUME_SAMPLES = 1200;
      const OVIZ_AR_VOLUME_SCAN_TARGET = 180000;
      const OVIZ_AR_VOLUME_COLOR_BINS = 8;
      const OVIZ_AR_VOLUME_MIN_STRETCHED_VALUE = 0.055;
      const OVIZ_AR_SKY_TEXTURE_HIGH = { width: 4096, height: 2048 };
      const OVIZ_AR_SKY_TEXTURE_LOW = { width: 2048, height: 1024 };
      let ovizArDialogEl = null;
      let ovizArObjectUrl = "";
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
          : "Open AR export options. Select clusters before exporting.";
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

      function ovizArDownsampleTrail(points, maxPoints) {
        const source = Array.isArray(points) ? points : [];
        const safeMax = Math.max(Math.round(Number(maxPoints) || 0), 2);
        if (source.length <= safeMax) {
          return source.slice();
        }
        const selected = [];
        const used = new Set();
        for (let index = 0; index < safeMax; index += 1) {
          const sourceIndex = Math.round(index * (source.length - 1) / (safeMax - 1));
          if (!used.has(sourceIndex)) {
            used.add(sourceIndex);
            selected.push(source[sourceIndex]);
          }
        }
        return selected;
      }

      function ovizArCollectOrbitTrails(selectedKeys) {
        const allowedKeys = new Set(Array.from(selectedKeys || []).slice(0, OVIZ_AR_MAX_TRAIL_CLUSTERS));
        const byKey = new Map();
        const seenTimeByKey = new Map();
        (Array.isArray(frameSpecs) ? frameSpecs : []).forEach((frame, frameIndex) => {
          const traces = frame && Array.isArray(frame.traces) ? frame.traces : [];
          traces.forEach((trace) => {
            const points = trace && Array.isArray(trace.points) ? trace.points : [];
            points.forEach((point) => {
              const selection = selectionForPoint(point, trace);
              const key = selection ? normalizedSelectionKeyFor(selection) : "";
              if (!key || !allowedKeys.has(key)) {
                return;
              }
              const coords = ovizArPointCoordinates(point, selection);
              if (!coords) {
                return;
              }
              if (!byKey.has(key)) {
                byKey.set(key, []);
                seenTimeByKey.set(key, new Set());
              }
              const timeMyr = ovizArFiniteNumber(frame && frame.time, frameIndex);
              const timeKey = Number.isFinite(timeMyr) ? timeMyr.toFixed(9) : String(frameIndex);
              if (seenTimeByKey.get(key).has(timeKey)) {
                return;
              }
              seenTimeByKey.get(key).add(timeKey);
              byKey.get(key).push({
                x: coords.x,
                y: coords.y,
                z: coords.z,
                timeMyr,
                color: ovizArPointColor(point, trace, selection),
              });
            });
          });
        });
        const trails = [];
        byKey.forEach((points, key) => {
          const ordered = points
            .slice()
            .sort((left, right) => ovizArFiniteNumber(left.timeMyr, 0.0) - ovizArFiniteNumber(right.timeMyr, 0.0));
          if (ordered.length >= 2) {
            trails.push({
              key,
              points: ovizArDownsampleTrail(ordered, OVIZ_AR_MAX_TRAIL_POINTS_PER_CLUSTER),
            });
          }
        });
        return trails;
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
              seen.add(key);
              points.push(record);
            });
          });
        }
        const cappedPoints = points.slice(0, OVIZ_AR_MAX_POINTS);
        const includedKeys = new Set(cappedPoints.map((point) => point.key));
        return {
          mode: String(mode || "3d"),
          selectionMode: selectedKeys.size > 0 ? "selection" : (maskOnlySelection ? "volume-lasso" : "present-day-scene"),
          selectedCount: selectedKeys.size,
          hasLassoMask: Boolean(selectionMask),
          presentFrameIndex: frameIndex,
          presentTimeMyr: ovizArFiniteNumber(frame && frame.time, 0.0),
          points: cappedPoints,
          trails: ovizArCollectOrbitTrails(includedKeys),
          pointLimit: OVIZ_AR_MAX_POINTS,
          labelLimit: OVIZ_AR_MAX_LABELS,
          truncated: points.length > cappedPoints.length,
        };
      }

      function ovizArSceneCenter(points) {
        const source = Array.isArray(points) ? points : [];
        if (!source.length) {
          return { x: 0.0, y: 0.0, z: 0.0 };
        }
        const sum = source.reduce((acc, point) => ({
          x: acc.x + ovizArFiniteNumber(point.x, 0.0),
          y: acc.y + ovizArFiniteNumber(point.y, 0.0),
          z: acc.z + ovizArFiniteNumber(point.z, 0.0),
        }), { x: 0.0, y: 0.0, z: 0.0 });
        return {
          x: sum.x / source.length,
          y: sum.y / source.length,
          z: sum.z / source.length,
        };
      }

      function ovizArSceneTransform(points, targetRadiusMeters = 0.75) {
        const center = ovizArSceneCenter(points);
        let maxDistance = 1.0;
        (Array.isArray(points) ? points : []).forEach((point) => {
          const dx = ovizArFiniteNumber(point.x, center.x) - center.x;
          const dy = ovizArFiniteNumber(point.y, center.y) - center.y;
          const dz = ovizArFiniteNumber(point.z, center.z) - center.z;
          const distance = Math.sqrt((dx * dx) + (dy * dy) + (dz * dz));
          if (Number.isFinite(distance)) {
            maxDistance = Math.max(maxDistance, distance);
          }
        });
        return {
          center,
          scale: targetRadiusMeters / maxDistance,
        };
      }

      function ovizArVectorFromPoint(point, transform) {
        const center = transform && transform.center ? transform.center : { x: 0, y: 0, z: 0 };
        const scale = Number(transform && transform.scale) || 1.0;
        const x = (ovizArFiniteNumber(point && point.x, 0.0) - center.x) * scale;
        const y = (ovizArFiniteNumber(point && point.y, 0.0) - center.y) * scale;
        const z = (ovizArFiniteNumber(point && point.z, 0.0) - center.z) * scale;
        return new THREE.Vector3(x, z, y);
      }

      function ovizArColor(value, fallback = "#ffffff") {
        try {
          return new THREE.Color(value || fallback);
        } catch (_err) {
          return new THREE.Color(fallback);
        }
      }

      function ovizArCreateLabelTexture(label, color) {
        const canvasEl = document.createElement("canvas");
        canvasEl.width = 512;
        canvasEl.height = 128;
        const ctx = canvasEl.getContext("2d");
        ctx.clearRect(0, 0, canvasEl.width, canvasEl.height);
        ctx.font = "700 42px Helvetica, Arial, sans-serif";
        ctx.textBaseline = "middle";
        ctx.fillStyle = "rgba(0, 0, 0, 0.58)";
        ctx.fillRect(0, 18, canvasEl.width, 92);
        ctx.strokeStyle = String(color || "#ffffff");
        ctx.lineWidth = 5;
        ctx.strokeRect(2.5, 20.5, canvasEl.width - 5, 87);
        ctx.fillStyle = "#ffffff";
        const text = String(label || "Cluster").slice(0, 32);
        ctx.fillText(text, 28, 64, canvasEl.width - 56);
        const texture = new THREE.CanvasTexture(canvasEl);
        texture.minFilter = THREE.LinearFilter;
        texture.magFilter = THREE.LinearFilter;
        texture.generateMipmaps = false;
        texture.needsUpdate = true;
        return texture;
      }

      function ovizArAddLabelPlane(group, label, position, color, index) {
        const texture = ovizArCreateLabelTexture(label, color);
        const material = new THREE.MeshBasicMaterial({
          map: texture,
          transparent: true,
          side: THREE.DoubleSide,
          depthWrite: false,
        });
        const plane = new THREE.Mesh(new THREE.PlaneGeometry(0.34, 0.085), material);
        plane.position.copy(position).add(new THREE.Vector3(0.05, 0.06 + ((index % 3) * 0.025), 0.035));
        plane.rotation.set(-0.35, 0.15, 0.0);
        plane.userData.ovizArLabel = true;
        group.add(plane);
      }

      function ovizArAddCylinderBetween(group, start, end, radius, material) {
        const delta = new THREE.Vector3().subVectors(end, start);
        const length = delta.length();
        if (!(length > 1e-6)) {
          return;
        }
        const geometry = new THREE.CylinderGeometry(radius, radius, 1.0, 8, 1, false);
        const mesh = new THREE.Mesh(geometry, material);
        mesh.position.copy(start).add(end).multiplyScalar(0.5);
        mesh.quaternion.setFromUnitVectors(new THREE.Vector3(0, 1, 0), delta.clone().normalize());
        mesh.scale.set(1.0, length, 1.0);
        group.add(mesh);
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
        return {
          x: ovizArFiniteNumber(xBounds[0], -0.5) + ((ix + 0.5) / Math.max(nx, 1)) * xSpan,
          y: ovizArFiniteNumber(yBounds[0], -0.5) + ((iy + 0.5) / Math.max(ny, 1)) * ySpan,
          z: ovizArFiniteNumber(zBounds[0], -0.5) + ((iz + 0.5) / Math.max(nz, 1)) * zSpan,
          cellSizePc: Math.max(
            1e-6,
            Math.min(
              Math.abs(xSpan / Math.max(nx, 1)),
              Math.abs(ySpan / Math.max(ny, 1)),
              Math.abs(zSpan / Math.max(nz, 1))
            )
          ),
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
        const samples = [];
        const layerSummaries = [];

        for (const layer of layers) {
          const state = ovizArVolumeStateForLayer(layer);
          const scalarData = await ovizArVolumeScalarArrayFor(layer);
          if (!state || !scalarData || !scalarData.length) {
            continue;
          }
          const shape = layer.shape || {};
          const nx = Math.max(1, Math.round(Number(shape.x) || 0));
          const ny = Math.max(1, Math.round(Number(shape.y) || 0));
          const nz = Math.max(1, Math.round(Number(shape.z) || 0));
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
          const stride = Math.max(1, Math.ceil(Math.cbrt(totalVoxels / scanTarget)));
          const layerName = typeof volumeBaseNameForLayer === "function"
            ? volumeBaseNameForLayer(layer)
            : String(layer.name || layer.key || "Volume");
          let layerSampleCount = 0;

          for (let iz = 0; iz < nz; iz += stride) {
            for (let iy = 0; iy < ny; iy += stride) {
              for (let ix = 0; ix < nx; ix += stride) {
                const dataIndex = (iz * ny * nx) + (iy * nx) + ix;
                const normalizedValue = (Number(scalarData[dataIndex]) || 0) / 255.0;
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
                  Math.min(lutSamples - 1, Math.round(stretchedValue * (lutSamples - 1)))
                ) * 4;
                const sample = {
                  key: String(layer.key || layerName),
                  label: layerName || "Volume",
                  x: position.x,
                  y: position.y,
                  z: position.z,
                  sizePc: position.cellSizePc * stride * 1.12,
                  value: normalizedValue,
                  weight: stretchedValue * Math.max(ovizArFiniteNumber(state.opacity, 0.25), 0.08),
                  color: ovizArHexColorFromBytes(colorBytes, colorIndex),
                  colorBin,
                };
                samples.push(sample);
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
            });
          }
        }

        samples.sort((left, right) => ovizArFiniteNumber(right.weight, 0.0) - ovizArFiniteNumber(left.weight, 0.0));
        const cappedSamples = samples.slice(0, maxSamples);
        return {
          samples: cappedSamples,
          layers: layerSummaries,
          truncated: samples.length > cappedSamples.length,
          totalCandidates: samples.length,
        };
      }

      function ovizArAppendBoxGeometry(vertices, indices, center, size) {
        const cx = ovizArFiniteNumber(center && center.x, 0.0);
        const cy = ovizArFiniteNumber(center && center.y, 0.0);
        const cz = ovizArFiniteNumber(center && center.z, 0.0);
        const half = Math.max(ovizArFiniteNumber(size, 0.01), 1e-5) * 0.5;
        const baseIndex = vertices.length / 3;
        [
          [cx - half, cy - half, cz - half],
          [cx + half, cy - half, cz - half],
          [cx + half, cy + half, cz - half],
          [cx - half, cy + half, cz - half],
          [cx - half, cy - half, cz + half],
          [cx + half, cy - half, cz + half],
          [cx + half, cy + half, cz + half],
          [cx - half, cy + half, cz + half],
        ].forEach((vertex) => vertices.push(vertex[0], vertex[1], vertex[2]));
        [
          0, 1, 2, 0, 2, 3,
          4, 6, 5, 4, 7, 6,
          0, 4, 5, 0, 5, 1,
          1, 5, 6, 1, 6, 2,
          2, 6, 7, 2, 7, 3,
          3, 7, 4, 3, 4, 0,
        ].forEach((index) => indices.push(baseIndex + index));
      }

      function ovizArAddVolumeProxy(group, volumeSamples, transform) {
        const samples = Array.isArray(volumeSamples) ? volumeSamples : [];
        if (!samples.length) {
          return;
        }
        const samplesByColor = new Map();
        samples.forEach((sample) => {
          const color = String(sample.color || "#ffffff");
          if (!samplesByColor.has(color)) {
            samplesByColor.set(color, []);
          }
          samplesByColor.get(color).push(sample);
        });
        samplesByColor.forEach((colorSamples, color) => {
          const vertices = [];
          const indices = [];
          colorSamples.forEach((sample) => {
            const center = ovizArVectorFromPoint(sample, transform);
            const size = Math.max(0.006, Math.min(0.045, ovizArFiniteNumber(sample.sizePc, 1.0) * ovizArFiniteNumber(transform && transform.scale, 1.0)));
            ovizArAppendBoxGeometry(vertices, indices, center, size);
          });
          if (!vertices.length || !indices.length) {
            return;
          }
          const geometry = new THREE.BufferGeometry();
          geometry.setAttribute("position", new THREE.Float32BufferAttribute(vertices, 3));
          geometry.setIndex(indices);
          geometry.computeVertexNormals();
          const material = new THREE.MeshStandardMaterial({
            color: ovizArColor(color, "#ffffff"),
            emissive: ovizArColor(color, "#ffffff"),
            emissiveIntensity: 0.10,
            roughness: 0.82,
            metalness: 0.0,
            transparent: true,
            opacity: 0.38,
          });
          const mesh = new THREE.Mesh(geometry, material);
          mesh.userData.ovizArVolumeProxy = true;
          group.add(mesh);
        });
      }

      async function buildOvizAr3DScene(snapshot, options = {}) {
        const includeLabels = options.includeLabels !== false;
        const volumeResult = await ovizArCollectVolumeSamples(snapshot);
        const sceneAr = new THREE.Scene();
        sceneAr.add(new THREE.AmbientLight(0xffffff, 0.9));
        const light = new THREE.DirectionalLight(0xffffff, 0.7);
        light.position.set(1.2, 2.0, 1.4);
        sceneAr.add(light);
        const group = new THREE.Group();
        sceneAr.add(group);
        const transformPoints = (snapshot.points || []).concat(volumeResult.samples || []);
        const transform = ovizArSceneTransform(transformPoints, 0.72);
        const sphereGeometry = new THREE.SphereGeometry(0.026, 24, 16);
        ovizArAddVolumeProxy(group, volumeResult.samples, transform);
        (snapshot.trails || []).forEach((trail) => {
          const trailPoints = Array.isArray(trail.points) ? trail.points : [];
          if (trailPoints.length < 2) {
            return;
          }
          const trailColor = ovizArColor((trailPoints[trailPoints.length - 1] || {}).color, "#ffffff");
          const trailMaterial = new THREE.MeshStandardMaterial({
            color: trailColor,
            roughness: 0.72,
            metalness: 0.0,
            transparent: true,
            opacity: 0.58,
          });
          for (let index = 1; index < trailPoints.length; index += 1) {
            const start = ovizArVectorFromPoint(trailPoints[index - 1], transform);
            const end = ovizArVectorFromPoint(trailPoints[index], transform);
            ovizArAddCylinderBetween(group, start, end, 0.0045, trailMaterial);
          }
        });
        (snapshot.points || []).forEach((point, index) => {
          const color = ovizArColor(point.color, "#ffffff");
          const material = new THREE.MeshStandardMaterial({
            color,
            emissive: color,
            emissiveIntensity: 0.35,
            roughness: 0.45,
            metalness: 0.0,
          });
          const marker = new THREE.Mesh(sphereGeometry, material);
          const position = ovizArVectorFromPoint(point, transform);
          const nStars = ovizArFiniteNumber(point.nStars, NaN);
          const sizeScale = Number.isFinite(nStars)
            ? Math.max(0.85, Math.min(1.85, Math.sqrt(Math.max(nStars, 1.0)) / 8.0))
            : 1.0;
          marker.scale.setScalar(sizeScale);
          marker.position.copy(position);
          marker.userData.ovizArKey = point.key;
          group.add(marker);
          if (includeLabels && index < OVIZ_AR_MAX_LABELS) {
            ovizArAddLabelPlane(group, point.label, position, point.color, index);
          }
        });
        const base = new THREE.Mesh(
          new THREE.CircleGeometry(0.86, 96),
          new THREE.MeshBasicMaterial({
            color: 0x27313d,
            transparent: true,
            opacity: 0.18,
            side: THREE.DoubleSide,
          })
        );
        base.rotation.x = -Math.PI / 2;
        base.position.y = -0.04;
        group.add(base);
        sceneAr.userData.ovizArSummary = {
          mode: "3d",
          pointCount: snapshot.points.length,
          volumeSampleCount: volumeResult.samples.length,
          volumeLayerCount: volumeResult.layers.length,
          volumeTruncated: volumeResult.truncated,
          presentTimeMyr: snapshot.presentTimeMyr,
          scalePcToMeters: transform.scale,
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

      async function ovizArLoadSkyTexture(target) {
        const width = Math.round(Number(target && target.width) || OVIZ_AR_SKY_TEXTURE_LOW.width);
        const height = Math.round(Number(target && target.height) || OVIZ_AR_SKY_TEXTURE_LOW.height);
        const url = skyDomeHips2FitsUrl(width, height);
        const texture = await ovizArTextureLoaderLoad(url);
        return { texture, width, height, url };
      }

      async function buildOvizArSkyDomeScene(snapshot, target = OVIZ_AR_SKY_TEXTURE_HIGH, options = {}) {
        const includeLabels = options.includeLabels === true;
        const loaded = await ovizArLoadSkyTexture(target);
        const radius = 1.45;
        const sceneAr = new THREE.Scene();
        sceneAr.add(new THREE.AmbientLight(0xffffff, 1.0));
        const domeGeometry = new THREE.SphereGeometry(radius, 96, 48);
        const domeMaterial = new THREE.MeshBasicMaterial({
          map: loaded.texture,
          side: THREE.BackSide,
        });
        const dome = new THREE.Mesh(domeGeometry, domeMaterial);
        sceneAr.add(dome);
        const markerGeometry = new THREE.SphereGeometry(0.018, 18, 12);
        (snapshot.points || []).forEach((point, index) => {
          const direction = ovizArSkyDirectionForPoint(point, radius * 0.985);
          if (!direction) {
            return;
          }
          const color = ovizArColor(point.color, "#ffffff");
          const material = new THREE.MeshBasicMaterial({ color });
          const marker = new THREE.Mesh(markerGeometry, material);
          marker.position.set(direction.x, direction.y, direction.z);
          marker.userData.ovizArKey = point.key;
          sceneAr.add(marker);
          if (includeLabels && index < OVIZ_AR_MAX_LABELS) {
            ovizArAddLabelPlane(sceneAr, point.label, marker.position, point.color, index);
          }
        });
        sceneAr.userData.ovizArSummary = {
          mode: "sky",
          pointCount: snapshot.points.length,
          textureWidth: loaded.width,
          textureHeight: loaded.height,
          textureUrl: loaded.url,
          coordsys: typeof skyDomeHips2FitsCoordsys === "function" ? skyDomeHips2FitsCoordsys() : "galactic",
        };
        return sceneAr;
      }

      function ovizArLoadScript(url) {
        return new Promise((resolve, reject) => {
          const script = document.createElement("script");
          script.src = url;
          script.async = true;
          script.onload = resolve;
          script.onerror = reject;
          document.head.appendChild(script);
        });
      }

      async function loadOvizFflateGlobal() {
        if (typeof window !== "undefined" && window.fflate && window.fflate.zipSync && window.fflate.strToU8) {
          return window.fflate;
        }
        const urls = [
          "https://cdn.jsdelivr.net/npm/three@0.128.0/examples/js/libs/fflate.min.js",
          "https://unpkg.com/three@0.128.0/examples/js/libs/fflate.min.js",
        ];
        let lastError = null;
        for (const url of urls) {
          try {
            await ovizArLoadScript(url);
            if (typeof window !== "undefined" && window.fflate && window.fflate.zipSync && window.fflate.strToU8) {
              return window.fflate;
            }
          } catch (err) {
            lastError = err;
          }
        }
        throw lastError || new Error("fflate is unavailable for USDZ export.");
      }

      async function loadOvizUSDZExporter() {
        if (THREE.USDZExporter) {
          await loadOvizFflateGlobal();
          return THREE.USDZExporter;
        }
        const legacyUrls = [
          "https://cdn.jsdelivr.net/npm/three@0.128.0/examples/js/exporters/USDZExporter.js",
          "https://unpkg.com/three@0.128.0/examples/js/exporters/USDZExporter.js",
        ];
        let legacyDependencyLoaded = false;
        try {
          await loadOvizFflateGlobal();
          legacyDependencyLoaded = true;
        } catch (_err) {
          legacyDependencyLoaded = false;
        }
        for (const url of legacyUrls) {
          try {
            if (!legacyDependencyLoaded) {
              break;
            }
            await ovizArLoadScript(url);
            if (THREE.USDZExporter) {
              return THREE.USDZExporter;
            }
          } catch (_err) {
            // Try the module exporter below.
          }
        }
        const moduleUrls = [
          "https://cdn.jsdelivr.net/npm/three@0.128.0/examples/jsm/exporters/USDZExporter.js",
          "https://unpkg.com/three@0.128.0/examples/jsm/exporters/USDZExporter.js",
        ];
        let lastError = null;
        for (const url of moduleUrls) {
          try {
            const module = await import(url);
            if (module && module.USDZExporter) {
              return module.USDZExporter;
            }
          } catch (err) {
            lastError = err;
          }
        }
        throw lastError || new Error("USDZExporter is unavailable.");
      }

      async function ovizArParseUsdZ(sceneAr) {
        const USDZExporter = await loadOvizUSDZExporter();
        const exporter = new USDZExporter();
        normalizeOvizArSceneForUSDZ(sceneAr);
        const result = await Promise.resolve(exporter.parse(sceneAr, { quickLookCompatible: true }));
        if (result instanceof ArrayBuffer) {
          return result;
        }
        if (result && result.buffer instanceof ArrayBuffer) {
          return result.buffer;
        }
        if (result instanceof Blob) {
          return await result.arrayBuffer();
        }
        throw new Error("USDZ exporter returned an unsupported payload.");
      }

      function normalizeOvizArMaterialForUSDZ(material) {
        if (!material || typeof material !== "object") {
          return;
        }
        const textureSlots = ["map", "normalMap", "aoMap", "roughnessMap", "metalnessMap", "emissiveMap"];
        textureSlots.forEach((slot) => {
          if (material[slot] === undefined) {
            material[slot] = null;
          }
        });
        if (!material.color) {
          material.color = new THREE.Color(0xffffff);
        }
        if (!material.emissive) {
          material.emissive = new THREE.Color(0x000000);
        }
        if (material.roughness === undefined) {
          material.roughness = 0.7;
        }
        if (material.metalness === undefined) {
          material.metalness = 0.0;
        }
      }

      function normalizeOvizArSceneForUSDZ(sceneAr) {
        if (!sceneAr || typeof sceneAr.traverse !== "function") {
          return;
        }
        sceneAr.traverse((object) => {
          if (!object || !object.isMesh) {
            return;
          }
          const materials = Array.isArray(object.material) ? object.material : [object.material];
          materials.forEach(normalizeOvizArMaterialForUSDZ);
        });
      }

      function ovizArClearObjectUrl() {
        if (ovizArObjectUrl) {
          URL.revokeObjectURL(ovizArObjectUrl);
          ovizArObjectUrl = "";
        }
      }

      function ovizArQuickLookAvailable() {
        const nav = window.navigator || {};
        const userAgent = String(nav.userAgent || "");
        const platform = String(nav.platform || "");
        return /iPad|iPhone|iPod/i.test(userAgent)
          || (platform === "MacIntel" && Number(nav.maxTouchPoints) > 1);
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

      function ovizArSetReady(blob, filename, message) {
        ovizArClearObjectUrl();
        ovizArObjectUrl = URL.createObjectURL(blob);
        const dialog = ensureOvizArDialog();
        const hiddenLink = dialog.querySelector(".oviz-three-ar-hidden-link");
        const openButton = dialog.querySelector(".oviz-three-ar-open-link");
        const downloadLink = dialog.querySelector(".oviz-three-ar-download");
        if (hiddenLink) {
          hiddenLink.href = ovizArObjectUrl;
          hiddenLink.download = filename;
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
        ovizArSetStatus(
          message || (ovizArQuickLookAvailable() ? "AR snapshot is ready." : "USDZ snapshot is ready."),
          ovizArQuickLookAvailable() ? "ready" : "warn"
        );
        if (ovizArQuickLookAvailable() && hiddenLink) {
          window.setTimeout(() => {
            try {
              hiddenLink.click();
            } catch (_err) {
              // The explicit Open button remains available.
            }
          }, 80);
        }
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
            try {
              sceneAr = await buildOvizArSkyDomeScene(snapshot, OVIZ_AR_SKY_TEXTURE_HIGH, { includeLabels: false });
              arrayBuffer = await ovizArParseUsdZ(sceneAr);
            } catch (_err) {
              sceneAr = await buildOvizArSkyDomeScene(snapshot, OVIZ_AR_SKY_TEXTURE_LOW, { includeLabels: false });
              arrayBuffer = await ovizArParseUsdZ(sceneAr);
            }
          } else {
            try {
              sceneAr = await buildOvizAr3DScene(snapshot, { includeLabels: true });
              if (!snapshot.points.length && !(Number((sceneAr.userData.ovizArSummary || {}).volumeSampleCount) > 0)) {
                throw new Error("No volume samples survived AR downsampling.");
              }
              arrayBuffer = await ovizArParseUsdZ(sceneAr);
            } catch (_err) {
              sceneAr = await buildOvizAr3DScene(snapshot, { includeLabels: false });
              if (!snapshot.points.length && !(Number((sceneAr.userData.ovizArSummary || {}).volumeSampleCount) > 0)) {
                throw new Error("No volume samples survived AR downsampling.");
              }
              arrayBuffer = await ovizArParseUsdZ(sceneAr);
            }
          }
          const blob = new Blob([arrayBuffer], { type: "model/vnd.usdz+zip" });
          const pointText = `${snapshot.points.length} cluster${snapshot.points.length === 1 ? "" : "s"}`;
          const summary = sceneAr && sceneAr.userData ? sceneAr.userData.ovizArSummary || {} : {};
          const volumeCount = Math.round(Number(summary.volumeSampleCount) || 0);
          const volumeText = volumeCount ? ` plus ${volumeCount} volume sample${volumeCount === 1 ? "" : "s"}` : "";
          const selectionText = snapshot.selectionMode === "selection" ? "selected " : "";
          const suffix = (snapshot.truncated || summary.volumeTruncated) ? " AR limits applied." : "";
          ovizArSetReady(blob, filename, `${selectionText}${pointText}${volumeText} exported at t=0 Myr.${suffix}`);
        } catch (err) {
          const message = err && err.message ? String(err.message) : "AR export failed.";
          ovizArSetStatus(`AR export failed: ${message}`, "error");
        } finally {
          ovizArSetBusy(false);
        }
      }

      function ovizArResetDialogActions() {
        const dialog = ensureOvizArDialog();
        const openButton = dialog.querySelector(".oviz-three-ar-open-link");
        const downloadLink = dialog.querySelector(".oviz-three-ar-download");
        if (openButton) {
          openButton.hidden = true;
          openButton.disabled = true;
        }
        if (downloadLink) {
          downloadLink.hidden = true;
          downloadLink.removeAttribute("href");
          downloadLink.removeAttribute("download");
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
""".strip()
