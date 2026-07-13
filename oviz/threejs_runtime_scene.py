from __future__ import annotations


THREEJS_SCENE_RUNTIME_JS = """
      function selectedVolumeLayer() {
        if (!activeVolumeKey) {
          return null;
        }
        return frameVolumeLayerForStateKey(activeVolumeKey);
      }

      function selectedVolumeState() {
        if (!activeVolumeKey) {
          return null;
        }
        return volumeStateByKey[String(activeVolumeKey)] || null;
      }

      function clampVolumeStateForLayer(layer, state) {
        if (!layer || !state) {
          return;
        }
        const dataRange = Array.isArray(layer.data_range) ? layer.data_range : [0.0, 1.0];
        const dataMin = Number(dataRange[0]);
        const dataMax = Number(dataRange[1]);
        if (!Number.isFinite(state.vmin)) {
          state.vmin = Number((layer.default_controls || {}).vmin);
        }
        if (!Number.isFinite(state.vmax)) {
          state.vmax = Number((layer.default_controls || {}).vmax);
        }
        if (Number.isFinite(dataMin) && Number.isFinite(dataMax)) {
          state.vmin = Math.min(Math.max(state.vmin, dataMin), dataMax);
          state.vmax = Math.min(Math.max(state.vmax, dataMin), dataMax);
        }
        if (!(state.vmax > state.vmin)) {
          if (Number.isFinite(dataMax) && dataMax > state.vmin) {
            state.vmax = dataMax;
          } else {
            state.vmax = state.vmin + 1e-6;
          }
        }
        state.opacity = Math.min(Math.max(Number(state.opacity), 0.0), 1.0);
        state.steps = Math.round(Math.min(Math.max(Number(state.steps), 24.0), 768.0));
        if (!Number.isFinite(state.steps) || state.steps < 24) {
          state.steps = Number((layer.default_controls || {}).steps || 100);
        }
        state.alphaCoef = Math.min(Math.max(Number(state.alphaCoef), 1.0), 200.0);
        if (!Number.isFinite(state.alphaCoef)) {
          state.alphaCoef = Number((layer.default_controls || {}).alpha_coef || 50.0);
        }
        state.gradientStep = Math.min(Math.max(Number(state.gradientStep), 1e-4), 0.05);
        if (!Number.isFinite(state.gradientStep)) {
          state.gradientStep = Number((layer.default_controls || {}).gradient_step || 0.005);
        }
        state.stretch = normalizeVolumeStretch(
          state.stretch !== undefined ? state.stretch : (layer.default_controls || {}).stretch
        );
      }

      function volumeColormapOptionFor(layer, colormapName) {
        const options = (layer && layer.colormap_options) || [];
        const requested = String(colormapName || "").trim().toLowerCase();
        for (const option of options) {
          if (String(option.name || "").trim().toLowerCase() === requested) {
            return option;
          }
        }
        return options.length ? options[0] : null;
      }

      function volumeStretchOptions() {
        return [
          { value: "linear", label: "Linear" },
          { value: "log10", label: "log10" },
          { value: "asinh", label: "asinh" },
        ];
      }

      function normalizeVolumeStretch(stretchName) {
        const requested = String(stretchName || "").trim().toLowerCase();
        const option = volumeStretchOptions().find((item) => item.value === requested);
        return option ? option.value : "linear";
      }

      function volumeStretchModeValue(stretchName) {
        const stretch = normalizeVolumeStretch(stretchName);
        if (stretch === "log10") {
          return 1.0;
        }
        if (stretch === "asinh") {
          return 2.0;
        }
        return 0.0;
      }

      function startPngAtlasVolumeDecode(layer, cacheKey, onDecoded) {
        if (volumeScalarDataPendingCache.has(cacheKey)) {
          return;
        }
        const shape = layer.shape || {};
        const nx = Math.max(1, Math.round(Number(shape.x) || 0));
        const ny = Math.max(1, Math.round(Number(shape.y) || 0));
        const nz = Math.max(1, Math.round(Number(shape.z) || 0));
        const tiles = layer.data_atlas_tiles || {};
        const tileCols = Math.max(1, Math.round(Number(tiles.x) || Math.ceil(Math.sqrt(nz))));
        const tileRows = Math.max(1, Math.round(Number(tiles.y) || Math.ceil(nz / tileCols)));
        const pending = new Promise((resolve) => {
          const image = new Image();
          image.onload = () => {
            const atlasWidth = Math.max(1, Number(image.naturalWidth || image.width || 0));
            const atlasHeight = Math.max(1, Number(image.naturalHeight || image.height || 0));
            const canvasEl = document.createElement("canvas");
            canvasEl.width = atlasWidth;
            canvasEl.height = atlasHeight;
            const ctx = canvasEl.getContext("2d");
            if (!ctx) {
              volumeScalarDataPendingCache.delete(cacheKey);
              resolve(null);
              return;
            }
            ctx.drawImage(image, 0, 0, atlasWidth, atlasHeight);
            const imageData = ctx.getImageData(0, 0, atlasWidth, atlasHeight);
            const rgba = imageData.data || [];
            const values = new Uint8Array(nx * ny * nz);
            for (let zIndex = 0; zIndex < nz; zIndex += 1) {
              const tileRow = Math.floor(zIndex / tileCols);
              const tileCol = zIndex % tileCols;
              if (tileRow >= tileRows) {
                break;
              }
              const xOffset = tileCol * nx;
              const yOffset = tileRow * ny;
              for (let yIndex = 0; yIndex < ny; yIndex += 1) {
                for (let xIndex = 0; xIndex < nx; xIndex += 1) {
                  const atlasIndex = (((yOffset + yIndex) * atlasWidth) + (xOffset + xIndex)) * 4;
                  const voxelIndex = (zIndex * ny * nx) + (yIndex * nx) + xIndex;
                  values[voxelIndex] = Number(rgba[atlasIndex] || 0);
                }
              }
            }
            volumeScalarDataCache.set(cacheKey, values);
            volumeScalarDataPendingCache.delete(cacheKey);
            if (typeof onDecoded === "function") {
              onDecoded(values);
            }
            resolve(values);
          };
          image.onerror = () => {
            volumeScalarDataPendingCache.delete(cacheKey);
            resolve(null);
          };
          image.src = `data:image/png;base64,${String(layer.data_b64 || "")}`;
        });
        volumeScalarDataPendingCache.set(cacheKey, pending);
      }

      function volumeScalarArrayFor(layer) {
        const layerKey = String(layer.key);
        if (volumeScalarDataCache.has(layerKey)) {
          return volumeScalarDataCache.get(layerKey);
        }
        const encoding = String(layer.data_encoding || "uint16_le");
        if (encoding === "png_atlas_uint8") {
          const shape = layer.shape || {};
          const nx = Math.max(1, Math.round(Number(shape.x) || 0));
          const ny = Math.max(1, Math.round(Number(shape.y) || 0));
          const nz = Math.max(1, Math.round(Number(shape.z) || 0));
          const placeholder = new Uint8Array(nx * ny * nz);
          volumeScalarDataCache.set(layerKey, placeholder);
          startPngAtlasVolumeDecode(layer, layerKey, (values) => {
            const volumeTexture = volumeTextureCache.get(layerKey);
            if (volumeTexture && volumeTexture.texture && values) {
              const image = volumeTexture.texture.image || {};
              image.data = values;
              image.width = volumeTexture.nx;
              image.height = volumeTexture.ny;
              image.depth = volumeTexture.nz;
              volumeTexture.texture.image = image;
              volumeTexture.texture.needsUpdate = true;
            }
          });
          return placeholder;
        }
        let data = null;
        if (encoding === "uint16_le") {
          data = base64ToUint16Array(layer.data_b64 || "");
        } else {
          data = base64ToUint8Array(layer.data_b64 || "");
        }
        volumeScalarDataCache.set(layerKey, data);
        return data;
      }

      function volumeSkyScalarArrayFor(layer) {
        const layerKey = `${String(layer.key)}::sky-overlay`;
        if (volumeScalarDataCache.has(layerKey)) {
          return volumeScalarDataCache.get(layerKey);
        }
        if (!layer || !layer.sky_overlay_data_b64) {
          const fallback = volumeScalarArrayFor(layer);
          volumeScalarDataCache.set(layerKey, fallback);
          return fallback;
        }
        const encoding = String(layer.sky_overlay_data_encoding || layer.data_encoding || "uint8");
        if (encoding === "png_atlas_uint8") {
          if (!volumeScalarDataPendingCache.has(layerKey)) {
            const shape = layer.sky_overlay_shape || layer.shape || {};
            const nx = Math.max(1, Math.round(Number(shape.x) || 0));
            const ny = Math.max(1, Math.round(Number(shape.y) || 0));
            const nz = Math.max(1, Math.round(Number(shape.z) || 0));
            const tiles = layer.sky_overlay_atlas_tiles || {};
            const tileCols = Math.max(1, Math.round(Number(tiles.x) || Math.ceil(Math.sqrt(nz))));
            const tileRows = Math.max(1, Math.round(Number(tiles.y) || Math.ceil(nz / tileCols)));
            const pending = new Promise((resolve) => {
              const image = new Image();
              image.onload = () => {
                const atlasWidth = Math.max(1, Number(image.naturalWidth || image.width || 0));
                const atlasHeight = Math.max(1, Number(image.naturalHeight || image.height || 0));
                const canvasEl = document.createElement("canvas");
                canvasEl.width = atlasWidth;
                canvasEl.height = atlasHeight;
                const ctx = canvasEl.getContext("2d");
                if (!ctx) {
                  volumeScalarDataPendingCache.delete(layerKey);
                  resolve(null);
                  return;
                }
                ctx.drawImage(image, 0, 0, atlasWidth, atlasHeight);
                const imageData = ctx.getImageData(0, 0, atlasWidth, atlasHeight);
                const rgba = imageData.data || [];
                const values = new Uint8Array(nx * ny * nz);
                for (let zIndex = 0; zIndex < nz; zIndex += 1) {
                  const tileRow = Math.floor(zIndex / tileCols);
                  const tileCol = zIndex % tileCols;
                  if (tileRow >= tileRows) {
                    break;
                  }
                  const xOffset = tileCol * nx;
                  const yOffset = tileRow * ny;
                  for (let yIndex = 0; yIndex < ny; yIndex += 1) {
                    for (let xIndex = 0; xIndex < nx; xIndex += 1) {
                      const atlasIndex = (((yOffset + yIndex) * atlasWidth) + (xOffset + xIndex)) * 4;
                      const voxelIndex = (zIndex * ny * nx) + (yIndex * nx) + xIndex;
                      values[voxelIndex] = Number(rgba[atlasIndex] || 0);
                    }
                  }
                }
                volumeScalarDataCache.set(layerKey, values);
                volumeScalarDataPendingCache.delete(layerKey);
                if (skySpec.enabled) {
                  updateSkyPanel();
                }
                resolve(values);
              };
              image.onerror = () => {
                volumeScalarDataPendingCache.delete(layerKey);
                resolve(null);
              };
              image.src = `data:image/png;base64,${String(layer.sky_overlay_data_b64 || "")}`;
            });
            volumeScalarDataPendingCache.set(layerKey, pending);
          }
          return null;
        }
        let data = null;
        if (encoding === "uint16_le") {
          data = base64ToUint16Array(layer.sky_overlay_data_b64 || "");
        } else {
          data = base64ToUint8Array(layer.sky_overlay_data_b64 || "");
        }
        volumeScalarDataCache.set(layerKey, data);
        return data;
      }

      function volumeColorTextureFor(option) {
        const optionKey = String(option.name || "volume-colormap");
        if (volumeColorTextureCache.has(optionKey)) {
          return volumeColorTextureCache.get(optionKey);
        }
        const bytes = base64ToUint8Array(option.lut_b64 || "");
        const width = Math.max(1, Math.floor(bytes.length / 4));
        const texture = new THREE.DataTexture(bytes, width, 1, THREE.RGBAFormat);
        texture.minFilter = THREE.NearestFilter;
        texture.magFilter = THREE.NearestFilter;
        texture.wrapS = THREE.ClampToEdgeWrapping;
        texture.wrapT = THREE.ClampToEdgeWrapping;
        texture.unpackAlignment = 1;
        texture.generateMipmaps = false;
        texture.needsUpdate = true;
        if ("colorSpace" in texture && THREE.SRGBColorSpace) {
          texture.colorSpace = THREE.SRGBColorSpace;
        } else if ("encoding" in texture && THREE.sRGBEncoding) {
          texture.encoding = THREE.sRGBEncoding;
        }
        volumeColorTextureCache.set(optionKey, texture);
        return texture;
      }

      function volumeColorBytesForOption(option) {
        const optionKey = String((option && option.name) || "volume-colormap");
        if (volumeColorBytesCache.has(optionKey)) {
          return volumeColorBytesCache.get(optionKey);
        }
        const bytes = base64ToUint8Array((option && option.lut_b64) || "");
        volumeColorBytesCache.set(optionKey, bytes);
        return bytes;
      }

      function volumeTextureFor(layer) {
        const layerKey = String(layer.key);
        if (volumeTextureCache.has(layerKey)) {
          return volumeTextureCache.get(layerKey);
        }
        const data = volumeScalarArrayFor(layer);
        const shape = layer.shape || {};
        const nx = Math.max(1, Number(shape.x || 1));
        const ny = Math.max(1, Number(shape.y || 1));
        const nz = Math.max(1, Number(shape.z || 1));
        const VolumeTextureCtor = THREE.Data3DTexture || THREE.DataTexture3D;
        if (!VolumeTextureCtor) {
          throw new Error("Three.js volume textures are unavailable in this browser build.");
        }
        const texture = new VolumeTextureCtor(data, nx, ny, nz);
        texture.format = THREE.RedFormat;
        texture.type = THREE.UnsignedByteType;
        texture.minFilter = layer.interpolation === false ? THREE.NearestFilter : THREE.LinearFilter;
        texture.magFilter = layer.interpolation === false ? THREE.NearestFilter : THREE.LinearFilter;
        texture.wrapS = THREE.ClampToEdgeWrapping;
        texture.wrapT = THREE.ClampToEdgeWrapping;
        texture.wrapR = THREE.ClampToEdgeWrapping;
        texture.unpackAlignment = 1;
        texture.generateMipmaps = false;
        texture.needsUpdate = true;
        const volumeTexture = { texture, nx, ny, nz };
        volumeTextureCache.set(layerKey, volumeTexture);
        return volumeTexture;
      }

      function volumeJitterTextureFor() {
        if (volumeJitterTexture) {
          return volumeJitterTexture;
        }
        const bytes = new Uint8Array(64 * 64);
        for (let i = 0; i < bytes.length; i += 1) {
          bytes[i] = Math.floor(Math.random() * 256.0);
        }
        volumeJitterTexture = new THREE.DataTexture(bytes, 64, 64, THREE.RedFormat, THREE.UnsignedByteType);
        volumeJitterTexture.minFilter = THREE.LinearFilter;
        volumeJitterTexture.magFilter = THREE.LinearFilter;
        volumeJitterTexture.wrapS = THREE.MirroredRepeatWrapping;
        volumeJitterTexture.wrapT = THREE.MirroredRepeatWrapping;
        volumeJitterTexture.generateMipmaps = false;
        volumeJitterTexture.unpackAlignment = 1;
        volumeJitterTexture.needsUpdate = true;
        return volumeJitterTexture;
      }

      function normalizedVolumeWindowFor(layer, state) {
        const dataRange = Array.isArray(layer.data_range) ? layer.data_range : [0.0, 1.0];
        const dataMin = Number(dataRange[0]);
        const dataMax = Number(dataRange[1]);
        if (!Number.isFinite(dataMin) || !Number.isFinite(dataMax) || !(dataMax > dataMin)) {
          return { low: 0.0, high: 1.0 };
        }
        const span = dataMax - dataMin;
        const low = Math.min(Math.max((Number(state.vmin) - dataMin) / span, 0.0), 1.0);
        let high = Math.min(Math.max((Number(state.vmax) - dataMin) / span, 0.0), 1.0);
        if (!(high > low)) {
          high = Math.min(1.0, low + 1e-6);
        }
        return { low, high };
      }

      function volumeLayerTimeMyr(layer) {
        if (!layer) {
          return null;
        }
        const rawTime = layer.time_myr;
        if (rawTime === null || rawTime === undefined || rawTime === "" || rawTime === false) {
          return null;
        }
        const timeValue = Number(rawTime);
        return Number.isFinite(timeValue) ? timeValue : null;
      }

      function volumeSupportsShowAllTimes(layer) {
        return Boolean(
          layer
          && layer.supports_show_all_times
          && volumeLayerTimeMyr(layer) === null
        );
      }

      function volumeVisibleForFrame(layer, state, frame = currentFrame()) {
        if (!layer || !state || state.visible === false) {
          return false;
        }
        const frameTime = frame ? Number(frame.time) : 0.0;
        const layerTime = volumeLayerTimeMyr(layer);
        if (layerTime !== null) {
          return approximatelyZero(frameTime - layerTime);
        }
        if (state.showAllTimes && volumeSupportsShowAllTimes(layer)) {
          return true;
        }
        if (layer.only_at_t0 === false) {
          return true;
        }
        return approximatelyZero(frameTime);
      }

      function volumeRotationAngleForFrame(layer, state, frame = currentFrame()) {
        if (
          !layer
          || !state
          || !state.showAllTimes
          || !layer.co_rotate_with_frame
          || !volumeSupportsShowAllTimes(layer)
          || !Number.isFinite(volumeCoRotationRateRadPerMyr)
        ) {
          return 0.0;
        }
        const frameTime = frame ? Number(frame.time) : 0.0;
        if (!Number.isFinite(frameTime)) {
          return 0.0;
        }
        const referenceTime = Number.isFinite(Number(layer.reference_time_myr))
          ? Number(layer.reference_time_myr)
          : 0.0;
        return volumeCoRotationRateRadPerMyr * (frameTime - referenceTime);
      }

      function volumeQuaternionForZRotation(angleRad) {
        const halfAngle = 0.5 * (Number.isFinite(angleRad) ? angleRad : 0.0);
        return new THREE.Vector4(0.0, 0.0, Math.sin(halfAngle), Math.cos(halfAngle));
      }

      const VOLUME_VERTEX_SHADER = `
        uniform vec4 rotation;
        uniform vec4 translation;

        varying vec3 localPosition;
        varying vec3 transformedCameraPosition;
        varying vec3 transformedWorldPosition;

        vec3 rotate_vertex_position(vec3 pos, vec3 t, vec4 q) {
          vec3 p = pos.xyz - t.xyz;
          return p.xyz + 2.0 * cross(cross(p.xyz, q.xyz) + q.w * p.xyz, q.xyz) + t.xyz;
        }

        void main() {
          vec3 transformed = position;
          localPosition = position;
          vec4 worldPosition = modelMatrix * vec4(transformed, 1.0);
          transformedCameraPosition = rotate_vertex_position(cameraPosition.xyz, translation.xyz, rotation);
          transformedWorldPosition = rotate_vertex_position(worldPosition.xyz, translation.xyz, rotation);
          gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0);
        }
      `;

      const VOLUME_FRAGMENT_SHADER = `
        #include <common>
        #include <lights_pars_begin>

        precision highp sampler3D;

        uniform sampler3D volumeTexture;
        uniform sampler2D colormap;
        uniform sampler2D jitterTexture;
        uniform sampler2D selectionMaskTexture;
        uniform sampler2D selectionSourceSecondaryMaskTexture;
        uniform sampler2D selectionSourceTertiaryMaskTexture;
        uniform sampler2D selectionSourceQuaternaryMaskTexture;
        uniform sampler2D selectionTargetMaskTexture;
        uniform float low;
        uniform float high;
        uniform float opacity;
        uniform float samples;
        uniform float alpha_coef;
        uniform float gradient_step;
        uniform float stretch_mode;
        uniform vec4 scale;
        uniform vec4 translation;
        uniform vec4 rotation;
        uniform bool useSelectionPolygon;
        uniform mat4 selectionViewProjectionMatrix;
        uniform float selectionDimOutside;
        uniform bool useSelectionSourceSecondaryPolygon;
        uniform mat4 selectionSourceSecondaryViewProjectionMatrix;
        uniform bool useSelectionSourceTertiaryPolygon;
        uniform mat4 selectionSourceTertiaryViewProjectionMatrix;
        uniform bool useSelectionSourceQuaternaryPolygon;
        uniform mat4 selectionSourceQuaternaryViewProjectionMatrix;
        uniform float selectionSourceBlend;
        uniform vec4 selectionSourceWeights;
        uniform bool useSelectionTargetPolygon;
        uniform mat4 selectionTargetViewProjectionMatrix;
        uniform bool selectionTransitionActive;
        uniform float selectionTransitionProgress;

        varying vec3 localPosition;
        varying vec3 transformedCameraPosition;
        varying vec3 transformedWorldPosition;

        float inv_range;

        struct Ray {
          vec3 origin;
          vec3 direction;
          vec3 inv_direction;
          int sign[3];
        };

        vec3 aabb[2] = vec3[2](
          vec3(-0.5, -0.5, -0.5),
          vec3(0.5, 0.5, 0.5)
        );

        Ray makeRay(vec3 origin, vec3 direction) {
          vec3 inv_direction = vec3(1.0) / direction;
          return Ray(
            origin,
            direction,
            inv_direction,
            int[3](
              ((inv_direction.x < 0.0) ? 1 : 0),
              ((inv_direction.y < 0.0) ? 1 : 0),
              ((inv_direction.z < 0.0) ? 1 : 0)
            )
          );
        }

        void intersect(in Ray ray, in vec3 bounds[2], out float tmin, out float tmax) {
          float tymin;
          float tymax;
          float tzmin;
          float tzmax;
          tmin = (bounds[ray.sign[0]].x - ray.origin.x) * ray.inv_direction.x;
          tmax = (bounds[1 - ray.sign[0]].x - ray.origin.x) * ray.inv_direction.x;
          tymin = (bounds[ray.sign[1]].y - ray.origin.y) * ray.inv_direction.y;
          tymax = (bounds[1 - ray.sign[1]].y - ray.origin.y) * ray.inv_direction.y;
          tzmin = (bounds[ray.sign[2]].z - ray.origin.z) * ray.inv_direction.z;
          tzmax = (bounds[1 - ray.sign[2]].z - ray.origin.z) * ray.inv_direction.z;
          tmin = max(max(tmin, tymin), tzmin);
          tmax = min(min(tmax, tymax), tzmax);
        }

        float sampleVolume(vec3 pos) {
          return texture(volumeTexture, clamp(pos, vec3(0.0), vec3(1.0))).x;
        }

        vec3 rotate_vertex_position(vec3 pos, vec3 t, vec4 q) {
          vec3 p = pos.xyz - t.xyz;
          return p.xyz + 2.0 * cross(cross(p.xyz, q.xyz) + q.w * p.xyz, q.xyz) + t.xyz;
        }

        vec3 inverse_rotate_vertex_position(vec3 pos, vec3 t, vec4 q) {
          return rotate_vertex_position(pos, t, vec4(-q.xyz, q.w));
        }

        float asinhStretch(float value) {
          return log(value + sqrt(value * value + 1.0));
        }

        float applyStretch(float value) {
          float clampedValue = clamp(value, 0.0, 1.0);
          if (stretch_mode < 0.5) {
            return clampedValue;
          }
          if (stretch_mode < 1.5) {
            float strength = 999.0;
            return log(1.0 + strength * clampedValue) / log(1.0 + strength);
          }
          float asinhStrength = 10.0;
          return asinhStretch(asinhStrength * clampedValue) / asinhStretch(asinhStrength);
        }

        bool sampleInsideSelectionPolygon(vec3 worldPos) {
          if (!useSelectionPolygon) {
            return true;
          }
          vec4 clip = selectionViewProjectionMatrix * vec4(worldPos, 1.0);
          if (clip.w <= 0.0) {
            return false;
          }
          vec3 ndc = clip.xyz / clip.w;
          if (ndc.z < -1.0 || ndc.z > 1.0) {
            return false;
          }
          vec2 uv = vec2(ndc.x * 0.5 + 0.5, 0.5 - ndc.y * 0.5);
          if (uv.x <= 0.0 || uv.x >= 1.0 || uv.y <= 0.0 || uv.y >= 1.0) {
            return false;
          }
          return texture2D(selectionMaskTexture, uv).r > 0.5;
        }

        bool sampleInsideTargetSelectionPolygon(vec3 worldPos) {
          if (!useSelectionTargetPolygon) {
            return true;
          }
          vec4 clip = selectionTargetViewProjectionMatrix * vec4(worldPos, 1.0);
          if (clip.w <= 0.0) {
            return false;
          }
          vec3 ndc = clip.xyz / clip.w;
          if (ndc.z < -1.0 || ndc.z > 1.0) {
            return false;
          }
          vec2 uv = vec2(ndc.x * 0.5 + 0.5, 0.5 - ndc.y * 0.5);
          if (uv.x <= 0.0 || uv.x >= 1.0 || uv.y <= 0.0 || uv.y >= 1.0) {
            return false;
          }
          return texture2D(selectionTargetMaskTexture, uv).r > 0.5;
        }

        bool sampleInsideSourceSecondarySelectionPolygon(vec3 worldPos) {
          if (!useSelectionSourceSecondaryPolygon) {
            return true;
          }
          vec4 clip = selectionSourceSecondaryViewProjectionMatrix * vec4(worldPos, 1.0);
          if (clip.w <= 0.0) {
            return false;
          }
          vec3 ndc = clip.xyz / clip.w;
          if (ndc.z < -1.0 || ndc.z > 1.0) {
            return false;
          }
          vec2 uv = vec2(ndc.x * 0.5 + 0.5, 0.5 - ndc.y * 0.5);
          if (uv.x <= 0.0 || uv.x >= 1.0 || uv.y <= 0.0 || uv.y >= 1.0) {
            return false;
          }
          return texture2D(selectionSourceSecondaryMaskTexture, uv).r > 0.5;
        }

        bool sampleInsideSourceTertiarySelectionPolygon(vec3 worldPos) {
          if (!useSelectionSourceTertiaryPolygon) {
            return true;
          }
          vec4 clip = selectionSourceTertiaryViewProjectionMatrix * vec4(worldPos, 1.0);
          if (clip.w <= 0.0) {
            return false;
          }
          vec3 ndc = clip.xyz / clip.w;
          if (ndc.z < -1.0 || ndc.z > 1.0) {
            return false;
          }
          vec2 uv = vec2(ndc.x * 0.5 + 0.5, 0.5 - ndc.y * 0.5);
          if (uv.x <= 0.0 || uv.x >= 1.0 || uv.y <= 0.0 || uv.y >= 1.0) {
            return false;
          }
          return texture2D(selectionSourceTertiaryMaskTexture, uv).r > 0.5;
        }

        bool sampleInsideSourceQuaternarySelectionPolygon(vec3 worldPos) {
          if (!useSelectionSourceQuaternaryPolygon) {
            return true;
          }
          vec4 clip = selectionSourceQuaternaryViewProjectionMatrix * vec4(worldPos, 1.0);
          if (clip.w <= 0.0) {
            return false;
          }
          vec3 ndc = clip.xyz / clip.w;
          if (ndc.z < -1.0 || ndc.z > 1.0) {
            return false;
          }
          vec2 uv = vec2(ndc.x * 0.5 + 0.5, 0.5 - ndc.y * 0.5);
          if (uv.x <= 0.0 || uv.x >= 1.0 || uv.y <= 0.0 || uv.y >= 1.0) {
            return false;
          }
          return texture2D(selectionSourceQuaternaryMaskTexture, uv).r > 0.5;
        }

        vec3 worldGetNormal(in float px, in vec3 pos) {
          return normalize(vec3(
            px - sampleVolume(pos + vec3(gradient_step, 0.0, 0.0)),
            px - sampleVolume(pos + vec3(0.0, gradient_step, 0.0)),
            px - sampleVolume(pos + vec3(0.0, 0.0, gradient_step))
          ));
        }

        void main() {
          float jitter = texture2D(jitterTexture, gl_FragCoord.xy / 64.0).r;
          float tmin = 0.0;
          float tmax = 0.0;
          float px = 0.0;
          vec4 pxColor = vec4(0.0);
          vec4 value = vec4(0.0);
          vec3 direction = normalize(transformedWorldPosition - transformedCameraPosition);

          inv_range = 1.0 / max(high - low, 1e-6);
          aabb[0] = aabb[0] * scale.xyz + translation.xyz;
          aabb[1] = aabb[1] * scale.xyz + translation.xyz;
          intersect(makeRay(transformedCameraPosition, direction), aabb, tmin, tmax);

          if (tmax <= max(0.0, tmin)) {
            discard;
          }

          vec3 textcoord_end = localPosition + vec3(0.5);
          vec3 textcoord_start = textcoord_end - (tmax - max(0.0, tmin)) * direction / scale.xyz;
          vec3 textcoord_delta = textcoord_end - textcoord_start;
          int sampleCount = min(int(length(textcoord_delta) * samples), int(samples * 1.8));
          if (sampleCount <= 0) {
            discard;
          }

          textcoord_delta = textcoord_delta / float(sampleCount);
          textcoord_start = textcoord_start - textcoord_delta * (0.01 + 0.98 * jitter);
          vec3 textcoord = textcoord_start - textcoord_delta;
          float step = length(textcoord_delta);

          for (int count = 0; count < 2048; count++) {
            if (count >= sampleCount) {
              break;
            }

            textcoord += textcoord_delta;
            px = texture(volumeTexture, textcoord).x;
            float scaled_px = (px - low) * inv_range;

            if (scaled_px > 0.0) {
              scaled_px = applyStretch(min(scaled_px, 0.999));
              pxColor = texture(colormap, vec2(scaled_px, 0.5));
              vec3 worldPos = inverse_rotate_vertex_position(
                (textcoord - vec3(0.5)) * scale.xyz + translation.xyz,
                translation.xyz,
                rotation
              );
              bool insideSelection = sampleInsideSelectionPolygon(worldPos);
              float selectionWeight = insideSelection ? 1.0 : selectionDimOutside;
              if (selectionTransitionActive) {
                bool insideSourceSecondarySelection = sampleInsideSourceSecondarySelectionPolygon(worldPos);
                bool insideSourceTertiarySelection = sampleInsideSourceTertiarySelectionPolygon(worldPos);
                bool insideSourceQuaternarySelection = sampleInsideSourceQuaternarySelectionPolygon(worldPos);
                bool insideTargetSelection = sampleInsideTargetSelectionPolygon(worldPos);
                float sourceSelectionWeight = mix(
                  insideSelection ? 1.0 : 0.0,
                  insideSourceSecondarySelection ? 1.0 : 0.0,
                  clamp(selectionSourceBlend, 0.0, 1.0)
                );
                sourceSelectionWeight = dot(
                  selectionSourceWeights,
                  vec4(
                    insideSelection ? 1.0 : 0.0,
                    insideSourceSecondarySelection ? 1.0 : 0.0,
                    insideSourceTertiarySelection ? 1.0 : 0.0,
                    insideSourceQuaternarySelection ? 1.0 : 0.0
                  )
                );
                selectionWeight = mix(
                  sourceSelectionWeight,
                  insideTargetSelection ? 1.0 : 0.0,
                  clamp(selectionTransitionProgress, 0.0, 1.0)
                );
              }
              if (selectionWeight <= 0.001) {
                continue;
              }
              pxColor.a = 1.0 - pow(1.0 - clamp(pxColor.a * opacity, 0.0, 0.999), step * alpha_coef);
              pxColor.a *= selectionWeight;
              pxColor.a *= (1.0 - value.a);
              pxColor.rgb *= pxColor.a;

              #if NUM_DIR_LIGHTS > 0
              if (pxColor.a > 0.0) {
                vec4 addedLights = vec4(ambientLightColor / PI, 1.0);
                vec3 specularColor = vec3(0.0);
                vec3 normal = worldGetNormal(px, textcoord);
                vec3 lightDirection;
                float lightingIntensity;
                vec3 lightReflect;
                float specularFactor;

                #pragma unroll_loop_start
                for (int i = 0; i < NUM_DIR_LIGHTS; i++) {
                  lightDirection = directionalLights[i].direction;
                  lightingIntensity = clamp(dot(lightDirection, normal), 0.0, 1.0);
                  addedLights.rgb += directionalLights[i].color / PI * (0.2 + 0.8 * lightingIntensity);
                  lightReflect = normalize(reflect(lightDirection, normal));
                  specularFactor = dot(direction, lightReflect);
                  if (specularFactor > 0.0) {
                    specularColor += 0.002 * scaled_px * (1.0 / max(step, 1e-6))
                      * directionalLights[i].color / PI * pow(specularFactor, 250.0) * pxColor.a;
                  }
                }
                #pragma unroll_loop_end

                pxColor.rgb = pxColor.rgb * addedLights.xyz + specularColor;
              }
              #endif

              value += pxColor;
              if (value.a >= 0.99) {
                value.a = 1.0;
                break;
              }
            }
          }

          if (value.a <= 0.001) {
            discard;
          }
          gl_FragColor = value;
        }
      `;

      function createVolumeRuntime(layer) {
        if (!volumeSupported) {
          return null;
        }
        const state = volumeStateByKey[volumeStateKeyForLayer(layer)];
        if (!state) {
          return null;
        }
        clampVolumeStateForLayer(layer, state);
        const bounds = layer.bounds || {};
        const xBounds = Array.isArray(bounds.x) ? bounds.x : [-0.5, 0.5];
        const yBounds = Array.isArray(bounds.y) ? bounds.y : [-0.5, 0.5];
        const zBounds = Array.isArray(bounds.z) ? bounds.z : [-0.5, 0.5];
        const sizeX = Math.max(1e-6, Number(xBounds[1]) - Number(xBounds[0]));
        const sizeY = Math.max(1e-6, Number(yBounds[1]) - Number(yBounds[0]));
        const sizeZ = Math.max(1e-6, Number(zBounds[1]) - Number(zBounds[0]));
        const centerX = 0.5 * (Number(xBounds[0]) + Number(xBounds[1]));
        const centerY = 0.5 * (Number(yBounds[0]) + Number(yBounds[1]));
        const centerZ = 0.5 * (Number(zBounds[0]) + Number(zBounds[1]));
        const option = volumeColormapOptionFor(layer, state.colormap);
        if (!option) {
          return null;
        }
        const volumeTexture = volumeTextureFor(layer);
        const windowState = normalizedVolumeWindowFor(layer, state);
        const uniforms = THREE.UniformsUtils.merge([
          THREE.UniformsLib.lights,
          {
            volumeTexture: { value: volumeTexture.texture },
            colormap: { value: volumeColorTextureFor(option) },
            jitterTexture: { value: volumeJitterTextureFor() },
            low: { value: Number(windowState.low) },
            high: { value: Number(windowState.high) },
            opacity: { value: Number(state.opacity) },
            samples: { value: Number(state.steps) },
            alpha_coef: { value: Number(state.alphaCoef) },
            gradient_step: { value: Number(state.gradientStep) },
            stretch_mode: { value: volumeStretchModeValue(state.stretch) },
            scale: { value: new THREE.Vector4(sizeX, sizeY, sizeZ, 1.0) },
            translation: { value: new THREE.Vector4(centerX, centerY, centerZ, 1.0) },
            rotation: { value: new THREE.Vector4(0.0, 0.0, 0.0, 1.0) },
            useSelectionPolygon: { value: false },
            selectionViewProjectionMatrix: { value: new THREE.Matrix4() },
            selectionMaskTexture: { value: null },
            selectionDimOutside: { value: 1.0 },
            useSelectionSourceSecondaryPolygon: { value: false },
            selectionSourceSecondaryViewProjectionMatrix: { value: new THREE.Matrix4() },
            selectionSourceSecondaryMaskTexture: { value: null },
            useSelectionSourceTertiaryPolygon: { value: false },
            selectionSourceTertiaryViewProjectionMatrix: { value: new THREE.Matrix4() },
            selectionSourceTertiaryMaskTexture: { value: null },
            useSelectionSourceQuaternaryPolygon: { value: false },
            selectionSourceQuaternaryViewProjectionMatrix: { value: new THREE.Matrix4() },
            selectionSourceQuaternaryMaskTexture: { value: null },
            selectionSourceBlend: { value: 0.0 },
            selectionSourceWeights: { value: new THREE.Vector4(1.0, 0.0, 0.0, 0.0) },
            useSelectionTargetPolygon: { value: false },
            selectionTargetViewProjectionMatrix: { value: new THREE.Matrix4() },
            selectionTargetMaskTexture: { value: null },
            selectionTransitionActive: { value: false },
            selectionTransitionProgress: { value: 0.0 },
          },
        ]);
        applyLassoSelectionTransitionUniforms(uniforms);

        const material = new THREE.ShaderMaterial({
          uniforms,
          vertexShader: VOLUME_VERTEX_SHADER,
          fragmentShader: VOLUME_FRAGMENT_SHADER,
          transparent: true,
          side: THREE.BackSide,
          depthTest: true,
          depthWrite: false,
          lights: true,
        });

        const geometry = new THREE.BoxBufferGeometry(1, 1, 1);
        const mesh = new THREE.Mesh(geometry, material);
        mesh.position.set(centerX, centerY, centerZ);
        mesh.scale.set(sizeX, sizeY, sizeZ);
        mesh.renderOrder = -30;
        mesh.frustumCulled = false;
        const runtime = { mesh, material, geometry, layer };
        mesh.userData.ovizVolumeRuntime = runtime;
        applyVolumeStateToRuntime(layer, runtime);
        return runtime;
      }

      function applyVolumeStateToRuntime(layer, runtime, frame = currentFrame()) {
        if (!runtime || !runtime.material || !runtime.mesh) {
          return;
        }
        const state = volumeStateByKey[volumeStateKeyForLayer(layer)] || {};
        const option = volumeColormapOptionFor(layer, state.colormap);
        const windowState = normalizedVolumeWindowFor(layer, state);
        runtime.mesh.visible = volumeVisibleForFrame(layer, state, frame);
        runtime.material.uniforms.low.value = Number(windowState.low);
        runtime.material.uniforms.high.value = Number(windowState.high);
        runtime.material.uniforms.opacity.value = Number(state.opacity);
        runtime.material.uniforms.samples.value = Number(state.steps);
        runtime.material.uniforms.alpha_coef.value = Number(state.alphaCoef);
        runtime.material.uniforms.gradient_step.value = Number(state.gradientStep);
        runtime.material.uniforms.stretch_mode.value = volumeStretchModeValue(state.stretch);
        runtime.material.uniforms.rotation.value.copy(
          volumeQuaternionForZRotation(volumeRotationAngleForFrame(layer, state, frame))
        );
        applyLassoSelectionTransitionUniforms(runtime.material.uniforms);
        if (option) {
          runtime.material.uniforms.colormap.value = volumeColorTextureFor(option);
        }
      }

      function renderVolumeControls() {
        if (
          !volumePanelEl
          || !volumeSelectEl
          || !volumeVisibleEl
          || !volumeColormapEl
          || !volumeStretchEl
          || !volumeVMinEl
          || !volumeVMaxEl
          || !volumeOpacityEl
          || !volumeAlphaEl
          || !volumeStepsEl
          || !volumeOpacityLabelEl
          || !volumeAlphaLabelEl
          || !volumeStepsLabelEl
          || !volumeSummaryEl
        ) {
          return;
        }
        const enabled = volumeStateKeys.length > 0;
        volumePanelEl.dataset.enabled = enabled ? "true" : "false";
        if (!enabled) {
          return;
        }

        if (!activeVolumeKey || !Object.prototype.hasOwnProperty.call(volumeStateByKey, String(activeVolumeKey))) {
          activeVolumeKey = String(volumeStateKeys[0]);
        }

        const layer = selectedVolumeLayer();
        const state = selectedVolumeState();
        if (!state) {
          return;
        }
        if (layer) {
          clampVolumeStateForLayer(layer, state);
        }

        volumeSelectEl.innerHTML = "";
        const controlOptions = volumeControlOptions();
        const selectedControlKey = activeVolumeControlKey();
        controlOptions.forEach((option) => {
          const optionEl = document.createElement("option");
          optionEl.value = String(option.controlKey);
          optionEl.textContent = String(option.label);
          if (String(option.controlKey) === String(selectedControlKey)) {
            optionEl.selected = true;
          }
          volumeSelectEl.appendChild(optionEl);
        });
        volumeSelectEl.disabled = !volumeSupported || controlOptions.length <= 1;

        const controlLayer = layer || volumeLayerForKey(activeVolumeKey);
        if (!controlLayer) {
          return;
        }
        const controlVariantGroup = volumeVariantGroupForLayer(controlLayer);
        const smoothingLayers = controlVariantGroup
          ? volumeVariantLayersForGroup(controlVariantGroup)
          : [];
        const showSmoothingControl = smoothingLayers.length > 1;
        if (volumeSmoothingFieldEl) {
          volumeSmoothingFieldEl.style.display = showSmoothingControl ? "" : "none";
        }
        if (volumeSmoothingEl) {
          volumeSmoothingEl.innerHTML = "";
          smoothingLayers.forEach((variantLayer) => {
            const optionEl = document.createElement("option");
            optionEl.value = volumeStateKeyForLayer(variantLayer);
            optionEl.textContent = volumeVariantLabelForLayer(variantLayer);
            if (String(optionEl.value) === String(activeVolumeKey)) {
              optionEl.selected = true;
            }
            volumeSmoothingEl.appendChild(optionEl);
          });
          volumeSmoothingEl.disabled = !volumeSupported || !showSmoothingControl;
        }
        volumeVisibleEl.checked = state.visible !== false;
        volumeColormapEl.innerHTML = "";
        ((controlLayer.colormap_options || [])).forEach((option) => {
          const optionEl = document.createElement("option");
          optionEl.value = String(option.name);
          optionEl.textContent = String(option.label || option.name);
          if (String(option.name) === String(state.colormap)) {
            optionEl.selected = true;
          }
          volumeColormapEl.appendChild(optionEl);
        });
        volumeColormapEl.value = String(state.colormap);
        volumeStretchEl.innerHTML = "";
        volumeStretchOptions().forEach((option) => {
          const optionEl = document.createElement("option");
          optionEl.value = String(option.value);
          optionEl.textContent = String(option.label);
          if (String(option.value) === String(state.stretch)) {
            optionEl.selected = true;
          }
          volumeStretchEl.appendChild(optionEl);
        });
        volumeStretchEl.value = String(state.stretch);

        syncVolumeWindowInput(volumeVMinEl, state.vmin, controlLayer);
        syncVolumeWindowInput(volumeVMaxEl, state.vmax, controlLayer);
        volumeOpacityEl.value = String(state.opacity);
        volumeAlphaEl.value = String(state.alphaCoef);
        volumeStepsEl.value = String(state.steps);
        volumeOpacityLabelEl.textContent = `Opacity (${Number(state.opacity).toFixed(2)})`;
        volumeAlphaLabelEl.textContent = `Alpha coef (${Math.round(Number(state.alphaCoef))})`;
        volumeStepsLabelEl.textContent = `Samples (${Math.round(Number(state.steps))})`;
        volumeVisibleEl.disabled = !volumeSupported;
        volumeColormapEl.disabled = !volumeSupported;
        volumeStretchEl.disabled = !volumeSupported;
        volumeVMinEl.disabled = !volumeSupported;
        volumeVMaxEl.disabled = !volumeSupported;
        volumeOpacityEl.disabled = !volumeSupported;
        volumeAlphaEl.disabled = !volumeSupported;
        volumeStepsEl.disabled = !volumeSupported;

        volumeSummaryEl.textContent = volumeSummaryTextFor(controlLayer, state);
      }

      function interpolateNumber(fromValue, toValue, alpha, fallbackValue = 0.0) {
        const fromNumber = Number(fromValue);
        const toNumber = Number(toValue);
        if (Number.isFinite(fromNumber) && Number.isFinite(toNumber)) {
          return fromNumber + (toNumber - fromNumber) * alpha;
        }
        if (Number.isFinite(fromNumber)) {
          return fromNumber;
        }
        if (Number.isFinite(toNumber)) {
          return toNumber;
        }
        return Number(fallbackValue) || 0.0;
      }

      function cloneTracePoint(point) {
        if (!point || typeof point !== "object") {
          return point;
        }
        const clone = Object.assign({}, point);
        if (point.motion && typeof point.motion === "object") {
          clone.motion = Object.assign({}, point.motion);
        }
        if (point.selection && typeof point.selection === "object") {
          clone.selection = Object.assign({}, point.selection);
        }
        return clone;
      }

      function cloneTraceLabel(label) {
        return label && typeof label === "object" ? Object.assign({}, label) : label;
      }

      function ovizTransitionStablePointKey(point, trace = null) {
        const motionKey = normalizeMemberKey(motionKeyForPoint(point));
        if (motionKey) {
          return `motion:${motionKey}`;
        }
        const selection = point && point.selection && typeof point.selection === "object"
          ? point.selection
          : selectionForPoint(point, trace);
        const selectionKey = normalizedSelectionKeyFor(selection);
        return selectionKey ? `selection:${selectionKey}` : "";
      }

      function ovizTransitionIdentityRecords(items, stableKeyForItem, allowIndex, side) {
        const occurrences = new Map();
        return (Array.isArray(items) ? items : []).map((item, index) => {
          const stableKey = stableKeyForItem(item, index);
          if (stableKey) {
            const occurrence = Number(occurrences.get(stableKey) || 0);
            occurrences.set(stableKey, occurrence + 1);
            return { item, index, key: `${stableKey}#${occurrence}`, stable: true };
          }
          return {
            item,
            index,
            key: allowIndex ? `index:${index}` : `${side}-only:${index}`,
            stable: false,
          };
        });
      }

      function cloneTracePointWithPresence(point, presence) {
        const clone = cloneTracePoint(point);
        if (clone && typeof clone === "object") {
          clone.oviz_presence_opacity = clampRange(Number(presence), 0.0, 1.0);
        }
        return clone;
      }

      function cloneTraceLabelWithPresence(label, presence) {
        const clone = cloneTraceLabel(label);
        if (clone && typeof clone === "object") {
          clone.oviz_presence_opacity = clampRange(Number(presence), 0.0, 1.0);
        }
        return clone;
      }

      function cloneTraceWithPresence(trace, presence) {
        if (!trace || typeof trace !== "object") {
          return null;
        }
        const clone = Object.assign({}, trace);
        clone.oviz_presence_opacity = clampRange(Number(presence), 0.0, 1.0);
        if (Array.isArray(trace.points)) clone.points = trace.points.map(cloneTracePoint);
        if (Array.isArray(trace.labels)) clone.labels = trace.labels.map(cloneTraceLabel);
        if (Array.isArray(trace.segments)) clone.segments = trace.segments.map((segment) => (
          Array.isArray(segment) ? segment.slice() : segment
        ));
        return clone;
      }

      function interpolateTracePoint(pointA, pointB, alpha, timeValue) {
        const fallbackPoint = alpha < 0.5 ? pointA : pointB;
        if (!pointA || !pointB) {
          return cloneTracePoint(fallbackPoint);
        }
        const keyA = clusterFilterSelectionKeyForPoint(pointA);
        const keyB = clusterFilterSelectionKeyForPoint(pointB);
        if (keyA && keyB && keyA !== keyB) {
          return cloneTracePoint(fallbackPoint);
        }
        const blended = cloneTracePoint(fallbackPoint);
        blended.x = interpolateNumber(pointA.x, pointB.x, alpha, fallbackPoint && fallbackPoint.x);
        blended.y = interpolateNumber(pointA.y, pointB.y, alpha, fallbackPoint && fallbackPoint.y);
        blended.z = interpolateNumber(pointA.z, pointB.z, alpha, fallbackPoint && fallbackPoint.z);
        blended.size = interpolateNumber(pointA.size, pointB.size, alpha, fallbackPoint && fallbackPoint.size);
        blended.opacity = interpolateNumber(pointA.opacity, pointB.opacity, alpha, fallbackPoint && fallbackPoint.opacity);
        blended.oviz_presence_opacity = 1.0;
        if (blended.motion && typeof blended.motion === "object") {
          blended.motion.time_myr = Number.isFinite(timeValue) ? timeValue : Number(blended.motion.time_myr) || 0.0;
          const ageNowMyr = Number(blended.motion.age_now_myr);
          if (Number.isFinite(ageNowMyr) && Number.isFinite(blended.motion.time_myr)) {
            blended.motion.age_at_t_myr = ageNowMyr + blended.motion.time_myr;
          }
        }
        return blended;
      }

      function interpolateTraceLabel(labelA, labelB, alpha) {
        const fallbackLabel = alpha < 0.5 ? labelA : labelB;
        if (!labelA || !labelB) {
          return cloneTraceLabel(fallbackLabel);
        }
        if (String(labelA.text || "") !== String(labelB.text || "")) {
          return cloneTraceLabel(fallbackLabel);
        }
        const blended = cloneTraceLabel(fallbackLabel);
        blended.x = interpolateNumber(labelA.x, labelB.x, alpha, fallbackLabel && fallbackLabel.x);
        blended.y = interpolateNumber(labelA.y, labelB.y, alpha, fallbackLabel && fallbackLabel.y);
        blended.z = interpolateNumber(labelA.z, labelB.z, alpha, fallbackLabel && fallbackLabel.z);
        blended.size = interpolateNumber(labelA.size, labelB.size, alpha, fallbackLabel && fallbackLabel.size);
        blended.oviz_presence_opacity = 1.0;
        return blended;
      }

      function interpolateTraceSpec(traceA, traceB, alpha, timeValue) {
        const fallbackTrace = alpha < 0.5 ? traceA : traceB;
        if (!traceA || !traceB) {
          return fallbackTrace || null;
        }
        const blended = Object.assign({}, fallbackTrace);
        blended.oviz_presence_opacity = 1.0;
        blended.opacity = interpolateNumber(traceA.opacity, traceB.opacity, alpha, fallbackTrace && fallbackTrace.opacity);
        blended.default_opacity = interpolateNumber(
          traceA.default_opacity,
          traceB.default_opacity,
          alpha,
          fallbackTrace && fallbackTrace.default_opacity
        );
        blended.default_point_size = interpolateNumber(
          traceA.default_point_size,
          traceB.default_point_size,
          alpha,
          fallbackTrace && fallbackTrace.default_point_size
        );

        if (Array.isArray(traceA.segments) && Array.isArray(traceB.segments) && traceA.segments.length === traceB.segments.length) {
          blended.segments = traceA.segments.map((segmentA, index) => {
            const segmentB = traceB.segments[index];
            if (!Array.isArray(segmentA) || !Array.isArray(segmentB) || segmentA.length !== 6 || segmentB.length !== 6) {
              return alpha < 0.5 ? segmentA : segmentB;
            }
            return segmentA.map((value, coordinateIndex) => (
              interpolateNumber(value, segmentB[coordinateIndex], alpha, value)
            ));
          });
        }

        if (Array.isArray(traceA.points) || Array.isArray(traceB.points)) {
          const pointsA = Array.isArray(traceA.points) ? traceA.points : [];
          const pointsB = Array.isArray(traceB.points) ? traceB.points : [];
          const allowIndex = pointsA.length === pointsB.length;
          const recordsA = ovizTransitionIdentityRecords(
            pointsA,
            (point) => ovizTransitionStablePointKey(point, traceA),
            allowIndex,
            "source",
          );
          const recordsB = ovizTransitionIdentityRecords(
            pointsB,
            (point) => ovizTransitionStablePointKey(point, traceB),
            allowIndex,
            "target",
          );
          const targetByKey = new Map(recordsB.map((record) => [record.key, record]));
          const usedTargetKeys = new Set();
          blended.points = [];
          recordsA.forEach((recordA) => {
            const recordB = targetByKey.get(recordA.key);
            if (recordB) {
              usedTargetKeys.add(recordA.key);
              blended.points.push(interpolateTracePoint(recordA.item, recordB.item, alpha, timeValue));
            } else {
              blended.points.push(cloneTracePointWithPresence(recordA.item, 1.0 - alpha));
            }
          });
          recordsB.forEach((recordB) => {
            if (!usedTargetKeys.has(recordB.key)) {
              blended.points.push(cloneTracePointWithPresence(recordB.item, alpha));
            }
          });
        }

        if (Array.isArray(traceA.labels) || Array.isArray(traceB.labels)) {
          const labelsA = Array.isArray(traceA.labels) ? traceA.labels : [];
          const labelsB = Array.isArray(traceB.labels) ? traceB.labels : [];
          const allowIndex = labelsA.length === labelsB.length;
          const labelKey = (label) => String(label && label.text || "").trim();
          const recordsA = ovizTransitionIdentityRecords(labelsA, labelKey, allowIndex, "source");
          const recordsB = ovizTransitionIdentityRecords(labelsB, labelKey, allowIndex, "target");
          const targetByKey = new Map(recordsB.map((record) => [record.key, record]));
          const usedTargetKeys = new Set();
          blended.labels = [];
          recordsA.forEach((recordA) => {
            const recordB = targetByKey.get(recordA.key);
            if (recordB) {
              usedTargetKeys.add(recordA.key);
              blended.labels.push(interpolateTraceLabel(recordA.item, recordB.item, alpha));
            } else {
              blended.labels.push(cloneTraceLabelWithPresence(recordA.item, 1.0 - alpha));
            }
          });
          recordsB.forEach((recordB) => {
            if (!usedTargetKeys.has(recordB.key)) {
              blended.labels.push(cloneTraceLabelWithPresence(recordB.item, alpha));
            }
          });
        }

        return blended;
      }

      function interpolatedFrameSpecForValue(frameValue, timeValue) {
        const state = frameValueState(frameValue);
        const lowerFrame = frameSpecs[state.lowerIndex] || null;
        const upperFrame = frameSpecs[state.upperIndex] || lowerFrame;
        const easedAlpha = smoothstep01(state.alpha);
        if (!lowerFrame) {
          return null;
        }
        if (!upperFrame || state.upperIndex === state.lowerIndex || easedAlpha <= 1e-6) {
          return lowerFrame;
        }
        if (easedAlpha >= 1.0 - 1e-6) {
          return upperFrame;
        }

        const lowerTraces = Array.isArray(lowerFrame.traces) ? lowerFrame.traces : [];
        const upperTraces = Array.isArray(upperFrame.traces) ? upperFrame.traces : [];
        const upperTraceByKey = new Map();
        upperTraces.forEach((trace) => {
          const traceKey = trace && trace.key !== undefined && trace.key !== null
            ? String(trace.key)
            : "";
          if (traceKey) upperTraceByKey.set(traceKey, trace);
        });
        const usedUpperKeys = new Set();
        const blendedTraces = [];

        lowerTraces.forEach((lowerTrace) => {
          const traceKey = lowerTrace && lowerTrace.key !== undefined && lowerTrace.key !== null
            ? String(lowerTrace.key)
            : "";
          const upperTrace = traceKey ? upperTraceByKey.get(traceKey) : null;
          if (traceKey && upperTrace) {
            usedUpperKeys.add(traceKey);
          }
          const blendedTrace = upperTrace
            ? interpolateTraceSpec(lowerTrace, upperTrace, easedAlpha, timeValue)
            : cloneTraceWithPresence(lowerTrace, 1.0 - easedAlpha);
          if (blendedTrace) {
            blendedTraces.push(blendedTrace);
          }
        });

        upperTraces.forEach((upperTrace) => {
          const traceKey = upperTrace && upperTrace.key !== undefined && upperTrace.key !== null
            ? String(upperTrace.key)
            : "";
          if (traceKey && usedUpperKeys.has(traceKey)) {
            return;
          }
          if (upperTrace) {
            blendedTraces.push(cloneTraceWithPresence(upperTrace, easedAlpha));
          }
        });

        const decorationFrame = easedAlpha < 0.5 ? lowerFrame : upperFrame;
        const lowerDecorations = Array.isArray(lowerFrame.decorations) ? lowerFrame.decorations : [];
        const upperDecorations = Array.isArray(upperFrame.decorations) ? upperFrame.decorations : [];
        const decorationIdentity = (decoration, index) => {
          const key = String((decoration && decoration.key) || "");
          const kind = String((decoration && decoration.kind) || "");
          return key ? `${kind}:${key}` : `${kind}:#${index}`;
        };
        const lowerDecorationByIdentity = new Map(
          lowerDecorations.map((decoration, index) => [decorationIdentity(decoration, index), decoration])
        );
        const upperDecorationByIdentity = new Map(
          upperDecorations.map((decoration, index) => [decorationIdentity(decoration, index), decoration])
        );
        const blendedDecorations = (
          decorationFrame && Array.isArray(decorationFrame.decorations)
            ? decorationFrame.decorations
            : []
        ).map((decoration, index) => {
          const identity = decorationIdentity(decoration, index);
          const lowerDecoration = lowerDecorationByIdentity.get(identity);
          const upperDecoration = upperDecorationByIdentity.get(identity);
          if (!lowerDecoration || !upperDecoration) {
            return decoration;
          }
          const blended = Object.assign({}, decoration);
          ["opacity", "opacity_scale"].forEach((fieldName) => {
            if (
              Number.isFinite(Number(lowerDecoration[fieldName]))
              && Number.isFinite(Number(upperDecoration[fieldName]))
            ) {
              blended[fieldName] = interpolateNumber(
                lowerDecoration[fieldName],
                upperDecoration[fieldName],
                easedAlpha,
                decoration[fieldName]
              );
            }
          });
          return blended;
        });
        return {
          name: formatTick(timeValue),
          time: timeValue,
          traces: blendedTraces,
          decorations: blendedDecorations,
        };
      }

      function makeVectorObject(trace, materialBucket) {
        const options = arguments.length > 2 && arguments[2] ? arguments[2] : {};
        const vectors = Array.isArray(trace && trace.vectors) ? trace.vectors : [];
        if (!vectors.length) {
          return null;
        }
        const traceState = traceStyleStateForKey(trace.key);
        const vectorConfig = trace.vector && typeof trace.vector === "object" ? trace.vector : {};
        const color = (traceState && traceState.color)
          || trace.default_color
          || trace.color
          || ((trace.line || {}).color)
          || "#ffffff";
        let opacity = traceState ? clamp01(traceState.opacity) : clamp01(trace.opacity ?? trace.default_opacity ?? 1.0);
        opacity *= traceVisibilityOpacityMultiplier(trace);
        opacity *= clamp01(Number(trace.oviz_presence_opacity ?? 1.0));
        const focusedTraceKey = dendrogramFocusTraceKey();
        if (focusedTraceKey && String(trace.key) !== focusedTraceKey) {
          opacity *= 0.16;
        }
        if (opacity <= 0.001 && options.forceResident !== true) {
          return null;
        }
        const sizeScale = traceState ? Math.max(Number(traceState.sizeScale), 0.05) : 1.0;
        const pcPerKms = Math.max(Number(vectorConfig.pc_per_kms ?? trace.pc_per_kms ?? 1.0), 0.0) * sizeScale;
        if (!(pcPerKms > 0.0)) {
          return null;
        }
        const referenceFrames = (trace.reference_frames && typeof trace.reference_frames === "object")
          ? trace.reference_frames
          : {};
        const requestedReferenceKey = normalizeMemberKey(focusSelectionKey);
        const referenceFrame = requestedReferenceKey ? referenceFrames[requestedReferenceKey] : null;
        const subtractReferenceVelocity = Boolean(referenceFrame && vectorConfig.reference_velocity !== false);
        const referenceVx = subtractReferenceVelocity ? Number(referenceFrame.vx ?? 0.0) : 0.0;
        const referenceVy = subtractReferenceVelocity ? Number(referenceFrame.vy ?? 0.0) : 0.0;
        const referenceVz = subtractReferenceVelocity ? Number(referenceFrame.vz ?? 0.0) : 0.0;
        const headFraction = clampRange(Number(vectorConfig.head_fraction ?? 0.22), 0.02, 0.80);
        const headMinPc = Math.max(Number(vectorConfig.head_min_pc ?? 4.0), 0.0);
        const headMaxPc = Math.max(Number(vectorConfig.head_max_pc ?? 30.0), headMinPc);
        const headWidthFraction = clampRange(Number(vectorConfig.head_width_fraction ?? 0.48), 0.05, 1.5);
        const maxLengthBasePc = Math.max(Number(vectorConfig.max_length_pc ?? 0.0), 0.0);
        const maxLengthPc = maxLengthBasePc > 0.0 ? maxLengthBasePc * Math.max(sizeScale, 1.0) : 0.0;
        const lineWidth = Math.max(Number((trace.line || {}).width ?? vectorConfig.line_width ?? 1.35), 0.25);
        const positions = [];
        const fallbackUp = sceneUpVector instanceof THREE.Vector3
          ? sceneUpVector.clone()
          : new THREE.Vector3(0.0, 0.0, 1.0);
        const fallbackSideAxis = new THREE.Vector3(1.0, 0.0, 0.0);

        vectors.forEach((item) => {
          const x = Number(item && item.x);
          const y = Number(item && item.y);
          const z = Number(item && item.z);
          const vx = Number(item && item.vx);
          const vy = Number(item && item.vy);
          const vz = Number(item && item.vz);
          if (
            !Number.isFinite(x)
            || !Number.isFinite(y)
            || !Number.isFinite(z)
            || !Number.isFinite(vx)
            || !Number.isFinite(vy)
            || !Number.isFinite(vz)
          ) {
            return;
          }
          const velocity = new THREE.Vector3(vx - referenceVx, vy - referenceVy, vz - referenceVz);
          const speed = velocity.length();
          if (!(speed > 1e-8)) {
            return;
          }
          const direction = velocity.multiplyScalar(1.0 / speed);
          let lengthPc = speed * pcPerKms;
          if (maxLengthPc > 0.0) {
            lengthPc = Math.min(lengthPc, maxLengthPc);
          }
          if (!(lengthPc > 1e-6)) {
            return;
          }
          const start = new THREE.Vector3(x, y, z);
          const end = start.clone().add(direction.clone().multiplyScalar(lengthPc));
          positions.push(start.x, start.y, start.z, end.x, end.y, end.z);

          let side = new THREE.Vector3().crossVectors(direction, fallbackUp);
          if (side.lengthSq() <= 1e-10) {
            side = new THREE.Vector3().crossVectors(direction, fallbackSideAxis);
          }
          if (side.lengthSq() <= 1e-10) {
            return;
          }
          side.normalize();
          const up = new THREE.Vector3().crossVectors(side, direction).normalize();
          const headLength = Math.min(
            Math.max(lengthPc * headFraction, headMinPc * Math.max(sizeScale, 0.35)),
            Math.min(headMaxPc * Math.max(sizeScale, 0.35), lengthPc * 0.85)
          );
          const headWidth = headLength * headWidthFraction;
          const headBase = end.clone().add(direction.clone().multiplyScalar(-headLength));
          [side, side.clone().multiplyScalar(-1.0), up, up.clone().multiplyScalar(-1.0)].forEach((axis) => {
            const headPoint = headBase.clone().add(axis.multiplyScalar(headWidth));
            positions.push(end.x, end.y, end.z, headPoint.x, headPoint.y, headPoint.z);
          });
        });

        if (!positions.length) {
          return null;
        }
        const geometry = new LineSegmentsGeometry();
        geometry.setPositions(positions);
        const material = new LineMaterial({
          color,
          linewidth: lineWidth,
          dashed: false,
          transparent: opacity < 1.0,
          opacity,
          worldUnits: false,
        });
        material.resolution.set(root.clientWidth, root.clientHeight);
        const line = new LineSegments2(geometry, material);
        line.computeLineDistances();
        line.userData.ovizRetainedTrace = trace;
        line.userData.ovizRetainedKind = "vectors";
        materialBucket.push(material);
        return line;
      }

      function renderFrameScene(frame, displayedTimeMyr, options = {}) {
        const updateWidgets = options.updateWidgets !== false;
        const preserveCamera = options.preserveCamera === true;
        const renderRoot = options.targetGroup instanceof THREE.Group ? options.targetGroup : plotGroup;
        const resetRegistries = options.resetRegistries !== false;
        const clearScene = options.clearScene !== false;
        const forceResident = options.forceResident === true;
        const includeOverlays = options.includeOverlays !== false;
        if (!frame) {
          return;
        }
        if (!preserveCamera && zoomAnchorTracksFrame !== false) {
          currentZoomAnchorPoint = trackedZoomAnchorPointForFrame(frame);
        }
        if (!preserveCamera && galacticSimpleOrbitTargetTrackingActive && currentZoomAnchorPoint instanceof THREE.Vector3) {
          const orbitTargetDelta = currentZoomAnchorPoint.clone().sub(controls.target);
          if (orbitTargetDelta.lengthSq() > 1e-12) {
            controls.target.add(orbitTargetDelta);
            camera.position.add(orbitTargetDelta);
          }
        }
        tooltipEl.style.display = "none";
        if (resetRegistries) {
          hoverTargets.length = 0;
          cameraResponsivePointEntries.length = 0;
          cameraResponsiveImagePlaneEntries.length = 0;
          galacticReferenceOpacityGroups.length = 0;
          selectionSpriteEntriesByKey.clear();
          screenStableTextSprites.length = 0;
          frameLineMaterials.length = 0;
          volumeRuntimeByKey.clear();
        }
        if (clearScene) {
          clearGroup(renderRoot);
        }
        renderRoot.position.set(0.0, 0.0, 0.0);

        (Array.isArray(frame.traces) ? frame.traces : []).forEach((trace) => {
          if (!forceResident && !traceVisible(trace)) {
            return;
          }
          const galacticReferenceGroup = isGalacticReferenceTrace(trace) ? new THREE.Group() : null;
          const traceParent = galacticReferenceGroup || renderRoot;

          if (trace.vectors && trace.vectors.length) {
            const vectorTrace = makeVectorObject(trace, frameLineMaterials, { forceResident });
            if (vectorTrace) {
              traceParent.add(vectorTrace);
            }
          }
          if (trace.segments && trace.segments.length) {
            const line = makeLineObject(trace, frameLineMaterials);
            if (line) {
              traceParent.add(line);
            }
          }
          if (trace.points && trace.points.length) {
            addMarkerTrace(traceParent, trace, { forceResident });
          }
          if (trace.labels && trace.labels.length) {
            addTextTrace(traceParent, trace, { forceResident });
          }
          if (galacticReferenceGroup) {
            galacticReferenceGroup.traverse((object) => {
              const material = object && object.material;
              if (!material) return;
              const materials = Array.isArray(material) ? material : [material];
              materials.forEach((item) => {
                if (!item || !item.userData) return;
                item.userData.ovizTimelineBaseOpacity = Number(item.opacity ?? 1.0);
              });
            });
            galacticReferenceOpacityGroups.push(galacticReferenceGroup);
            renderRoot.add(galacticReferenceGroup);
          }
        });

        const selectionTransition = (
          typeof ovizStateTransition !== "undefined"
          && ovizStateTransition
          && ovizStateTransition.phasePlan
          && ovizStateTransition.phasePlan.changed.appearance
        ) ? ovizStateTransition : null;
        if (includeOverlays && selectionTransition) {
          const selectionProgress = clampRange(
            Number(selectionTransition.currentAppearanceProgress) || 0.0,
            0.0,
            1.0,
          );
          const fromSelection = selectionTransition.fromSnapshot
            && selectionTransition.fromSnapshot.current_selection;
          const toSelection = selectionTransition.targetSnapshot
            && selectionTransition.targetSnapshot.current_selection;
          const fromFootprint = buildSelectionFootprint(
            fromSelection,
            frameLineMaterials,
            1.0 - selectionProgress,
          );
          const toFootprint = buildSelectionFootprint(
            toSelection,
            frameLineMaterials,
            selectionProgress,
          );
          if (fromFootprint) renderRoot.add(fromFootprint);
          if (toFootprint) renderRoot.add(toFootprint);
        } else if (includeOverlays && currentSelectionMode === "click" && currentSelection) {
          const footprint = buildSelectionFootprint(currentSelection, frameLineMaterials);
          if (footprint) renderRoot.add(footprint);
        }

        (frame.decorations || []).forEach((decoration) => {
          addDecoration(renderRoot, decoration, { forceResident });
        });

        const selectionBoxGroup = includeOverlays ? buildSelectionBoxGroupForFrame(displayedTimeMyr) : null;
        if (selectionBoxGroup) {
          renderRoot.add(selectionBoxGroup);
        }

        if (includeOverlays) {
          addManualLabels(renderRoot);
        }
        updateTimelineMotionOpacity();

        const focusOffset = focusTrackingOffsetForFrame(frame);
        if (focusOffset) {
          renderRoot.position.copy(focusOffset.multiplyScalar(-1.0));
        }

        applySceneHoverState();
        updateClusterInfoTooltipPosition();
        if (updateWidgets) {
          resize();
          renderVolumeControls();
          renderBoxMetricsWidget();
          renderAgeKdeWidget();
          renderClusterFilterWidget();
          renderDendrogramWidget();
        }
      }

      let ovizRetainedTransitionScene = null;
      let ovizRetainedSceneBuildSerial = 0;
      let ovizRetainedSceneUpdateSerial = 0;

      function ovizRetainedTransitionOwner(options = {}) {
        if (options.preserveCamera !== true) {
          return "";
        }
        if (options.transitionOwnerToken !== undefined && options.transitionOwnerToken !== null) {
          return String(options.transitionOwnerToken || "");
        }
        if (typeof ovizStateTransition !== "undefined" && ovizStateTransition) {
          return String(
            ovizStateTransition.transitionId
            || ovizStateTransition.targetId
            || "state-transition"
          );
        }
        if (typeof activeActionRun !== "undefined" && activeActionRun) {
          return `action:${String(activeActionRun.id || "")}`;
        }
        if (typeof timeActionTrack !== "undefined" && timeActionTrack) {
          return [
            "action-time",
            Number(timeActionTrack.stepIndex) || 0,
            Number(timeActionTrack.startTime) || 0,
          ].join(":");
        }
        return "";
      }

      function disposeRetainedTransitionScene(ownerToken = "", options = {}) {
        const runtime = ovizRetainedTransitionScene;
        if (!runtime) return false;
        const requestedOwner = String(ownerToken || "");
        if (requestedOwner && requestedOwner !== runtime.ownerToken) return false;
        ovizRetainedTransitionScene = null;
        if (options.clear === true) {
          clearGroup(plotGroup);
        }
        return true;
      }

      function ovizRetainedMaterials(object) {
        if (!object || !object.material) {
          return [];
        }
        return Array.isArray(object.material) ? object.material : [object.material];
      }

      function ovizPrepareRetainedEndpoint(rootGroup, frame, frameIndex) {
        const pointEntries = [];
        let materialCount = 0;
        rootGroup.traverse((object) => {
          if (!object || !object.material) {
            return;
          }
          if (Array.isArray(object.material)) {
            object.material = object.material.map((material) => {
              if (!material || !material.userData || !material.userData.cached) {
                return material;
              }
              return ovizRetainedSpriteMaterial(material);
            });
          } else if (object.material.userData && object.material.userData.cached) {
            object.material = ovizRetainedSpriteMaterial(object.material);
          }
          ovizRetainedMaterials(object).forEach((material) => {
            if (!material) return;
            materialCount += 1;
            material.userData = Object.assign({}, material.userData || {});
            material.userData.ovizRetainedBaseOpacity = Number(material.opacity ?? 1.0);
          });
          const pointMetadata = object.userData && object.userData.ovizRetainedPoint;
          if (pointMetadata) {
            pointEntries.push({ object, metadata: pointMetadata });
          }
        });
        return {
          root: rootGroup,
          focusPosition: rootGroup.position.clone(),
          frame,
          frameIndex,
          pointEntries,
          materialCount,
        };
      }

      function ovizPointComponentIdentity(metadata, allowIndex) {
        if (!metadata) return "";
        const pointIdentity = metadata.stablePointKey
          ? `stable:${metadata.stablePointKey}#${Number(metadata.stableOccurrence) || 0}`
          : (allowIndex ? `index:${Number(metadata.pointIndex) || 0}` : "");
        if (!pointIdentity) return "";
        return `${String(metadata.traceKey || "")}\n${pointIdentity}\n${String(metadata.component || "")}`;
      }

      function ovizFramePointTopologyByTrace(frame) {
        const topology = new Map();
        (frame && Array.isArray(frame.traces) ? frame.traces : []).forEach((trace) => {
          const points = Array.isArray(trace && trace.points) ? trace.points : [];
          topology.set(String(trace && trace.key || ""), {
            count: points.length,
            positional: points.every((point) => !ovizTransitionStablePointKey(point)),
          });
        });
        return topology;
      }

      function ovizPairRetainedPoints(fromEndpoint, toEndpoint) {
        const fromTopology = ovizFramePointTopologyByTrace(fromEndpoint.frame);
        const toTopology = ovizFramePointTopologyByTrace(toEndpoint.frame);
        const canUseIndex = (traceKey) => {
          const from = fromTopology.get(traceKey);
          const to = toTopology.get(traceKey);
          return Boolean(
            from
            && to
            && from.positional
            && to.positional
            && Number(from.count) === Number(to.count)
          );
        };
        const targetByIdentity = new Map();
        toEndpoint.pointEntries.forEach((entry) => {
          const traceKey = String(entry.metadata.traceKey || "");
          const identity = ovizPointComponentIdentity(entry.metadata, canUseIndex(traceKey));
          if (identity) targetByIdentity.set(identity, entry);
        });
        const usedTargets = new Set();
        const pairs = [];
        fromEndpoint.pointEntries.forEach((fromEntry) => {
          const traceKey = String(fromEntry.metadata.traceKey || "");
          const identity = ovizPointComponentIdentity(fromEntry.metadata, canUseIndex(traceKey));
          const toEntry = identity ? targetByIdentity.get(identity) : null;
          if (toEntry) usedTargets.add(toEntry);
          pairs.push({
            from: fromEntry,
            to: toEntry || null,
            livePoint: Object.assign({}, fromEntry.metadata.point || {}),
          });
        });
        toEndpoint.pointEntries.forEach((toEntry) => {
          if (!usedTargets.has(toEntry)) {
            pairs.push({
              from: null,
              to: toEntry,
              livePoint: Object.assign({}, toEntry.metadata.point || {}),
            });
          }
        });
        return pairs;
      }

      function ovizRetainedCloneLivePoint(pair, alpha, displayedTimeMyr) {
        const fromPoint = pair.from && pair.from.metadata.point;
        const toPoint = pair.to && pair.to.metadata.point;
        const fallback = fromPoint || toPoint || {};
        const livePoint = pair.livePoint;
        Object.assign(livePoint, fallback);
        if (fromPoint && toPoint) {
          livePoint.x = interpolateNumber(fromPoint.x, toPoint.x, alpha, fallback.x);
          livePoint.y = interpolateNumber(fromPoint.y, toPoint.y, alpha, fallback.y);
          livePoint.z = interpolateNumber(fromPoint.z, toPoint.z, alpha, fallback.z);
          livePoint.size = interpolateNumber(fromPoint.size, toPoint.size, alpha, fallback.size);
          livePoint.opacity = interpolateNumber(fromPoint.opacity, toPoint.opacity, alpha, fallback.opacity);
        }
        if (fallback.motion && typeof fallback.motion === "object") {
          livePoint.motion = Object.assign({}, fallback.motion, { time_myr: displayedTimeMyr });
        }
        return livePoint;
      }

      function ovizPointBirthVisibility(pointState, point, trace) {
        // In the normal birth-time mode animatedPointState() expresses a
        // not-yet-born point by shrinking its size to zero.  Retained points
        // must carry that same presence into opacity and into the responsive
        // size floor; otherwise forceResident keeps them visible until the
        // exact final frame removes them in a single-frame pop.
        if (fadeOpacityByBirthTimeEnabled) {
          return 1.0;
        }
        const baseSize = Math.max(Number(pointSizeForTrace(point, trace)) || 0.0, 0.0);
        const animatedSize = Math.max(Number(pointState && pointState.size) || 0.0, 0.0);
        if (baseSize <= 1e-12) {
          return animatedSize > 1e-12 ? 1.0 : 0.0;
        }
        return clamp01(animatedSize / baseSize);
      }

      function ovizRetainedPointVisual(entry, livePoint, displayedTimeMyr) {
        if (!entry || !entry.object || !entry.metadata) {
          return { opacity: 0.0, pointScale: 0.0, color: "#ffffff" };
        }
        const metadata = entry.metadata;
        const trace = metadata.trace || {};
        const point = metadata.point || livePoint || {};
        const traceState = traceStyleStateForKey(trace.key);
        const pointState = animatedPointState(point, trace, displayedTimeMyr);
        const birthVisibility = ovizPointBirthVisibility(pointState, point, trace);
        const traceOpacityMultiplier = traceState
          ? clamp01(traceState.opacity) / Math.max(clamp01(Number(trace.default_opacity ?? 1.0)), 1e-6)
          : 1.0;
        let opacityMultiplier = clusterFilterPassesPoint(livePoint) ? 1.0 : 0.0;
        const pointKey = clusterFilterSelectionKeyForPoint(livePoint)
          || normalizedSelectionKeyFor(livePoint && livePoint.selection)
          || clusterFilterSelectionKeyForPoint(point)
          || normalizedSelectionKeyFor(point && point.selection);
        const focusedTraceKey = dendrogramFocusTraceKey();
        const dendrogramActiveKeys = activeDendrogramSelectionKeys();
        if (focusedTraceKey) {
          if (String(trace.key) !== focusedTraceKey) {
            opacityMultiplier *= 0.14;
          } else if (dendrogramActiveKeys.size && pointKey) {
            opacityMultiplier *= dendrogramActiveKeys.has(pointKey) ? 1.0 : 0.24;
          }
        }
        if (typeof ovizSelectionMembershipOpacity === "function") {
          opacityMultiplier *= ovizSelectionMembershipOpacity(pointKey, livePoint);
        } else if (lassoSelectionFilterActive()) {
          opacityMultiplier *= pointKey && selectedClusterKeys.has(pointKey) ? 1.0 : 0.0;
        }
        const baseOpacity = Number(pointState.opacity ?? pointOpacityForTrace(point, trace));
        const presenceOpacity = clamp01(Number(point.oviz_presence_opacity ?? 1.0))
          * clamp01(Number(trace.oviz_presence_opacity ?? 1.0));
        const effectiveOpacity = clamp01(
          baseOpacity
          * traceOpacityMultiplier
          * (traceVisible(trace) ? traceVisibilityOpacityMultiplier(trace) : 0.0)
          * opacityMultiplier
          * globalPointOpacityScale
          * presenceOpacity
          * birthVisibility
        );
        const sizeScaleFactor = traceState ? Math.max(Number(traceState.sizeScale), 0.05) : 1.0;
        const starsFactor = sizeByStarsFactorForPoint(point, trace, traceState);
        const scaleFloor = pointScale * 0.5 * Math.max(globalPointSizeScale, 0.05) * birthVisibility;
        const baseScale = Math.max(
          Math.max(Number(pointState.size) || 0.0, 0.0)
            * sizeScaleFactor
            * starsFactor
            * globalPointSizeScale
            * pointScale,
          scaleFloor,
        );
        const glowStrength = Math.max(Number(globalPointGlowStrength) || 0.0, 0.0);
        const glowMix = clampRange(glowStrength / 0.15, 0.0, 1.0);
        let componentOpacity = effectiveOpacity;
        let responsivePointScale = nonGlowMarkerScaleForPoint(baseScale, entry.object.position);
        if (metadata.component === "glow") {
          componentOpacity = clampRange(
            effectiveOpacity * (0.34 + 0.18 * glowStrength) * glowMix,
            0.0,
            0.78,
          );
          responsivePointScale = baseScale;
        } else if (metadata.component === "core") {
          componentOpacity = clampRange(
            effectiveOpacity * (1.00 + 0.24 * glowStrength) * glowMix,
            0.0,
            1.0,
          );
          responsivePointScale = baseScale;
        } else {
          componentOpacity *= 1.0 - glowMix;
        }
        return {
          opacity: componentOpacity,
          pointScale: responsivePointScale,
          color: metadata.component === "core"
            ? "#ffffff"
            : pointColorForTrace(point, trace, traceState),
        };
      }

      function ovizApplyRetainedPointEntry(entry, visual, endpointWeight, livePoint) {
        if (!entry || !entry.object) return;
        const object = entry.object;
        object.position.set(
          Number(livePoint && livePoint.x) || 0.0,
          Number(livePoint && livePoint.y) || 0.0,
          Number(livePoint && livePoint.z) || 0.0,
        );
        const material = object.material;
        if (material) {
          material.transparent = true;
          material.opacity = clamp01(Number(visual.opacity) || 0.0) * clamp01(endpointWeight);
          if (material.color && typeof material.color.set === "function") {
            material.color.set(String(visual.color || "#ffffff"));
          }
        }
        const responsive = object.userData && object.userData.ovizResponsivePointEntry;
        if (responsive) {
          responsive.position.copy(object.position);
          responsive.pointScale = Math.max(Number(visual.pointScale) || 0.0, 0.0);
        }
        // Match addMarkerTrace(), which omits exact-frame sprites at or below
        // 0.001 opacity.  Keeping a fainter retained sprite visible would
        // create a tiny endpoint pop and a false final-fidelity failure.
        object.visible = Boolean(material && material.opacity > 0.001);
      }

      function ovizApplyRetainedEndpointWeight(endpoint, weight, displayedTimeMyr) {
        const endpointWeight = clamp01(weight);
        endpoint.root.traverse((object) => {
          if (!object || !object.material || (object.userData && object.userData.ovizRetainedPoint)) {
            return;
          }
          const responsiveImage = object.userData && object.userData.ovizResponsiveImagePlaneEntry;
          if (responsiveImage) {
            responsiveImage.transitionOpacityScale = endpointWeight;
            return;
          }
          const volumeRuntime = object.userData && object.userData.ovizVolumeRuntime;
          if (volumeRuntime) {
            applyVolumeStateToRuntime(volumeRuntime.layer, volumeRuntime, endpoint.frame);
            const opacityUniform = volumeRuntime.material
              && volumeRuntime.material.uniforms
              && volumeRuntime.material.uniforms.opacity;
            if (opacityUniform) {
              opacityUniform.value = Number(opacityUniform.value || 0.0) * endpointWeight;
              object.visible = endpointWeight > 0.0001 && object.visible !== false;
            }
            return;
          }
          const trace = object.userData && object.userData.ovizRetainedTrace;
          const labelMetadata = object.userData && object.userData.ovizRetainedLabel;
          ovizRetainedMaterials(object).forEach((material) => {
            if (!material) return;
            let baseOpacity = Number(material.userData && material.userData.ovizRetainedBaseOpacity);
            if (!Number.isFinite(baseOpacity)) baseOpacity = Number(material.opacity ?? 1.0);
            if (trace) {
              const traceState = traceStyleStateForKey(trace.key);
              baseOpacity = (traceState
                ? clamp01(traceState.opacity)
                : clamp01(trace.opacity ?? trace.default_opacity ?? 1.0))
                * (traceVisible(trace) ? traceVisibilityOpacityMultiplier(trace) : 0.0)
                * clamp01(Number(trace.oviz_presence_opacity ?? 1.0));
              const color = (traceState && traceState.color)
                || trace.default_color
                || trace.color
                || ((trace.line || {}).color)
                || "#ffffff";
              if (material.color && typeof material.color.set === "function") material.color.set(color);
            } else if (labelMetadata) {
              const labelTrace = labelMetadata.trace || {};
              const label = labelMetadata.label || {};
              const traceState = traceStyleStateForKey(labelTrace.key);
              baseOpacity = (traceState ? clamp01(traceState.opacity) : 1.0)
                * (traceVisible(labelTrace) ? traceVisibilityOpacityMultiplier(labelTrace) : 0.0)
                * clamp01(Number(labelTrace.oviz_presence_opacity ?? 1.0))
                * clamp01(Number(label.oviz_presence_opacity ?? 1.0));
              const color = (traceState && traceState.color) || label.color || theme.axis_color;
              if (material.color && typeof material.color.set === "function") material.color.set(color);
            }
            material.transparent = true;
            material.opacity = clamp01(baseOpacity) * endpointWeight;
          });
        });
      }

      function ovizPrepareRetainedSelectionOverlay(runtime) {
        if (!runtime || !runtime.overlayRoot) return null;
        const materialBucket = [];
        const registerEndpoint = (group) => {
          if (!group) return null;
          runtime.overlayRoot.add(group);
          const materials = [];
          group.traverse((object) => {
            ovizRetainedMaterials(object).forEach((material) => {
              if (!material) return;
              material.userData = Object.assign({}, material.userData || {});
              material.userData.ovizRetainedBaseOpacity = Number(material.opacity ?? 1.0);
              materials.push(material);
            });
          });
          return { group, materials };
        };
        const createSelectionEndpoint = (selection) => registerEndpoint(
          buildSelectionFootprint(selection, materialBucket, 1.0)
        );
        const createManualLabelEndpoint = (labels) => {
          const group = new THREE.Group();
          addManualLabels(group, { labels, interactive: false });
          return group.children.length ? registerEndpoint(group) : null;
        };
        const createSelectionBoxEndpoint = (state, timeMyr) => registerEndpoint(
          buildSelectionBoxGroupForFrame(timeMyr, state, { interactive: false })
        );
        const stateAppearanceTransition = (
          typeof ovizStateTransition !== "undefined"
          && ovizStateTransition
          && ovizStateTransition.phasePlan
          && ovizStateTransition.phasePlan.changed.appearance
        ) ? ovizStateTransition : null;
        const selectionTransition = stateAppearanceTransition || (
          typeof actionHeldAppearanceRollback !== "undefined"
          && actionHeldAppearanceRollback
          && Number(actionHeldAppearanceRollback.frozenAppearanceProgress) > 0.0
            ? actionHeldAppearanceRollback
            : null
        );
        const sourceComposite = stateAppearanceTransition
          && stateAppearanceTransition.sourceAppearanceComposite
          ? stateAppearanceTransition.sourceAppearanceComposite
          : null;
        const compositeWeight = sourceComposite
          ? clampRange(Number(sourceComposite.targetWeight) || 0.0, 0.0, 1.0)
          : 0.0;
        let from = null;
        let composite = null;
        let to = null;
        if (selectionTransition) {
          const sourceSnapshot = sourceComposite
            ? sourceComposite.fromSnapshot
            : selectionTransition.fromSnapshot;
          const fromSelection = sourceSnapshot && sourceSnapshot.current_selection;
          const compositeSelection = sourceComposite
            && sourceComposite.targetSnapshot
            && sourceComposite.targetSnapshot.current_selection;
          const toSelection = selectionTransition.targetSnapshot
            && selectionTransition.targetSnapshot.current_selection;
          from = createSelectionEndpoint(fromSelection);
          composite = createSelectionEndpoint(compositeSelection);
          to = createSelectionEndpoint(toSelection);
        } else if (currentSelectionMode === "click" && currentSelection) {
          from = createSelectionEndpoint(currentSelection);
        }
        const fromSnapshot = sourceComposite
          ? sourceComposite.fromSnapshot
          : (selectionTransition ? selectionTransition.fromSnapshot : null);
        const compositeSnapshot = sourceComposite ? sourceComposite.targetSnapshot : null;
        const toSnapshot = selectionTransition ? selectionTransition.targetSnapshot : null;
        const fromLabels = fromSnapshot && Array.isArray(fromSnapshot.manual_labels)
          ? fromSnapshot.manual_labels
          : manualLabels;
        const toLabels = toSnapshot && Array.isArray(toSnapshot.manual_labels)
          ? toSnapshot.manual_labels
          : fromLabels;
        const compositeLabels = compositeSnapshot && Array.isArray(compositeSnapshot.manual_labels)
          ? compositeSnapshot.manual_labels
          : fromLabels;
        const manualCompositeSameFrom = JSON.stringify(fromLabels) === JSON.stringify(compositeLabels);
        const manualToSameFrom = JSON.stringify(fromLabels) === JSON.stringify(toLabels);
        const manualToSameComposite = JSON.stringify(compositeLabels) === JSON.stringify(toLabels);
        const manualFrom = createManualLabelEndpoint(fromLabels);
        const manualComposite = sourceComposite
          && !manualCompositeSameFrom
          ? createManualLabelEndpoint(compositeLabels)
          : null;
        const manualTo = selectionTransition
          && !manualToSameFrom
          && !(sourceComposite && manualToSameComposite)
          ? createManualLabelEndpoint(toLabels)
          : null;
        const fromBoxState = fromSnapshot && fromSnapshot.selection_box_state
          ? fromSnapshot.selection_box_state
          : {
              center_local_pc: {
                x: Number(selectionBoxState.center && selectionBoxState.center.x) || 0.0,
                y: Number(selectionBoxState.center && selectionBoxState.center.y) || 0.0,
                z: Number(selectionBoxState.center && selectionBoxState.center.z) || 0.0,
              },
              visible: selectionBoxState.visible !== false,
              half_width_pc: clampSelectionBoxHalfWidth(selectionBoxState.halfWidthPc),
            };
        const toBoxState = toSnapshot && toSnapshot.selection_box_state
          ? toSnapshot.selection_box_state
          : fromBoxState;
        const compositeBoxState = compositeSnapshot && compositeSnapshot.selection_box_state
          ? compositeSnapshot.selection_box_state
          : fromBoxState;
        const compositeBoxSameFrom = JSON.stringify(fromBoxState) === JSON.stringify(compositeBoxState);
        const toBoxSameFrom = JSON.stringify(fromBoxState) === JSON.stringify(toBoxState);
        const toBoxSameComposite = JSON.stringify(compositeBoxState) === JSON.stringify(toBoxState);
        const lowerTime = Number(runtime.fromEndpoint.frame && runtime.fromEndpoint.frame.time) || 0.0;
        const upperTime = Number(runtime.toEndpoint.frame && runtime.toEndpoint.frame.time) || lowerTime;
        const boxes = {
          fromLower: createSelectionBoxEndpoint(fromBoxState, lowerTime),
          fromUpper: createSelectionBoxEndpoint(fromBoxState, upperTime),
          compositeLower: null,
          compositeUpper: null,
          toLower: null,
          toUpper: null,
        };
        if (sourceComposite && !compositeBoxSameFrom) {
          boxes.compositeLower = createSelectionBoxEndpoint(compositeBoxState, lowerTime);
          boxes.compositeUpper = createSelectionBoxEndpoint(compositeBoxState, upperTime);
        }
        if (selectionTransition && !toBoxSameFrom && !(sourceComposite && toBoxSameComposite)) {
          boxes.toLower = createSelectionBoxEndpoint(toBoxState, lowerTime);
          boxes.toUpper = createSelectionBoxEndpoint(toBoxState, upperTime);
        }
        runtime.overlayLineMaterials = materialBucket;
        return {
          from,
          composite,
          to,
          manualFrom,
          manualComposite,
          manualTo,
          boxes,
          transition: selectionTransition,
          compositeWeight,
          manualCompositeSameFrom,
          manualToSameFrom,
          manualToSameComposite,
          compositeBoxSameFrom,
          toBoxSameFrom,
          toBoxSameComposite,
        };
      }

      function ovizApplyRetainedSelectionOverlay(runtime, frameAlpha = 0.0) {
        const overlay = runtime && runtime.selectionOverlay;
        if (!overlay) return;
        const progress = overlay.transition
          ? clampRange(Number(overlay.transition.currentAppearanceProgress) || 0.0, 0.0, 1.0)
          : 0.0;
        const applyEndpoint = (endpoint, weight) => {
          if (!endpoint) return;
          const clampedWeight = clamp01(weight);
          endpoint.group.visible = clampedWeight > 0.0001;
          endpoint.materials.forEach((material) => {
            const baseOpacity = Number(material.userData && material.userData.ovizRetainedBaseOpacity);
            material.transparent = true;
            material.opacity = (Number.isFinite(baseOpacity) ? baseOpacity : 1.0) * clampedWeight;
          });
        };
        const sourceWeight = overlay.transition ? 1.0 - progress : 1.0;
        const compositeWeight = clamp01(Number(overlay.compositeWeight) || 0.0);
        applyEndpoint(overlay.from, sourceWeight * (1.0 - compositeWeight));
        applyEndpoint(overlay.composite, sourceWeight * compositeWeight);
        applyEndpoint(overlay.to, progress);
        let manualFromWeight = sourceWeight * (1.0 - compositeWeight);
        let manualCompositeWeight = sourceWeight * compositeWeight;
        let manualToWeight = progress;
        if (overlay.manualCompositeSameFrom) {
          manualFromWeight += manualCompositeWeight;
          manualCompositeWeight = 0.0;
        }
        if (overlay.manualToSameFrom) {
          manualFromWeight += manualToWeight;
          manualToWeight = 0.0;
        } else if (overlay.manualToSameComposite) {
          manualCompositeWeight += manualToWeight;
          manualToWeight = 0.0;
        }
        applyEndpoint(overlay.manualFrom, manualFromWeight);
        applyEndpoint(overlay.manualComposite, manualCompositeWeight);
        applyEndpoint(overlay.manualTo, manualToWeight);
        const lowerWeight = 1.0 - clamp01(frameAlpha);
        const upperWeight = clamp01(frameAlpha);
        let boxFromWeight = sourceWeight * (1.0 - compositeWeight);
        let boxCompositeWeight = sourceWeight * compositeWeight;
        let boxToWeight = progress;
        if (overlay.compositeBoxSameFrom) {
          boxFromWeight += boxCompositeWeight;
          boxCompositeWeight = 0.0;
        }
        if (overlay.toBoxSameFrom) {
          boxFromWeight += boxToWeight;
          boxToWeight = 0.0;
        } else if (overlay.toBoxSameComposite) {
          boxCompositeWeight += boxToWeight;
          boxToWeight = 0.0;
        }
        applyEndpoint(overlay.boxes.fromLower, boxFromWeight * lowerWeight);
        applyEndpoint(overlay.boxes.fromUpper, boxFromWeight * upperWeight);
        applyEndpoint(overlay.boxes.compositeLower, boxCompositeWeight * lowerWeight);
        applyEndpoint(overlay.boxes.compositeUpper, boxCompositeWeight * upperWeight);
        applyEndpoint(overlay.boxes.toLower, boxToWeight * lowerWeight);
        applyEndpoint(overlay.boxes.toUpper, boxToWeight * upperWeight);
      }

      function ovizPrepareRetainedTransitionScene(ownerToken, frameState) {
        clearGroup(plotGroup);
        plotGroup.position.set(0.0, 0.0, 0.0);
        const fromRoot = new THREE.Group();
        const toRoot = new THREE.Group();
        const overlayRoot = new THREE.Group();
        plotGroup.add(fromRoot);
        plotGroup.add(toRoot);
        plotGroup.add(overlayRoot);
        const fromFrame = frameSpecs[frameState.lowerIndex] || null;
        const toFrame = frameSpecs[frameState.upperIndex] || fromFrame;
        renderFrameScene(fromFrame, Number(fromFrame && fromFrame.time) || 0.0, {
          updateWidgets: false,
          preserveCamera: true,
          targetGroup: fromRoot,
          resetRegistries: true,
          clearScene: false,
          forceResident: true,
          includeOverlays: false,
        });
        renderFrameScene(toFrame, Number(toFrame && toFrame.time) || 0.0, {
          updateWidgets: false,
          preserveCamera: true,
          targetGroup: toRoot,
          resetRegistries: false,
          clearScene: false,
          forceResident: true,
          includeOverlays: false,
        });
        const fromEndpoint = ovizPrepareRetainedEndpoint(fromRoot, fromFrame, frameState.lowerIndex);
        const toEndpoint = ovizPrepareRetainedEndpoint(toRoot, toFrame, frameState.upperIndex);
        ovizRetainedSceneBuildSerial += 1;
        const runtime = {
          ownerToken,
          intervalKey: `${frameState.lowerIndex}:${frameState.upperIndex}`,
          fromEndpoint,
          toEndpoint,
          pointPairs: ovizPairRetainedPoints(fromEndpoint, toEndpoint),
          overlayRoot,
          overlayLineMaterials: [],
          selectionOverlay: null,
          buildSerial: ovizRetainedSceneBuildSerial,
          updateCount: 0,
          materialCount: fromEndpoint.materialCount + toEndpoint.materialCount,
        };
        runtime.selectionOverlay = ovizPrepareRetainedSelectionOverlay(runtime);
        ovizRetainedTransitionScene = runtime;
        return runtime;
      }

      function ovizRetainedDebugSnapshot(runtime) {
        const effectiveTraceOpacities = {};
        let renderedObjectCount = 0;
        let hiddenObjectCount = 0;
        [runtime.fromEndpoint.root, runtime.toEndpoint.root, runtime.overlayRoot].forEach((endpointRoot) => {
          endpointRoot.traverse((object) => {
            if (!object || !object.material) return;
            const materials = ovizRetainedMaterials(object);
            const opacity = materials.reduce(
              (maximum, material) => Math.max(maximum, Number(material && material.opacity) || 0.0),
              0.0,
            );
            if (object.visible !== false && opacity > 0.0001) renderedObjectCount += 1;
            else hiddenObjectCount += 1;
            const pointMetadata = object.userData && object.userData.ovizRetainedPoint;
            const traceMetadata = object.userData && object.userData.ovizRetainedTrace;
            const labelMetadata = object.userData && object.userData.ovizRetainedLabel;
            const traceKey = String(
              (pointMetadata && pointMetadata.traceKey)
              || (traceMetadata && traceMetadata.key)
              || (labelMetadata && labelMetadata.trace && labelMetadata.trace.key)
              || ""
            );
            if (traceKey) {
              effectiveTraceOpacities[traceKey] = Math.max(
                Number(effectiveTraceOpacities[traceKey]) || 0.0,
                opacity,
              );
            }
          });
        });
        const selection = (
          typeof ovizStateSelectionTransition !== "undefined" && ovizStateSelectionTransition
        ) || (
          typeof ovizHeldSelectionTransition !== "undefined" && ovizHeldSelectionTransition
        ) || null;
        return {
          effectiveTraceOpacities,
          renderedObjectCount,
          hiddenObjectCount,
          lassoWeights: selection
            ? Object.fromEntries(Array.from(selection.fromWeightByKey || []).slice(0, 256))
            : {},
          lassoProgress: selection ? Number(selection.progress) || 0.0 : null,
          canvasOpacity: Number(renderer && renderer.domElement && renderer.domElement.style.opacity || 1.0),
        };
      }

      function ovizExpectedTraceOpacity(trace, displayedTimeMyr) {
        if (!trace || !trace.key || !traceVisible(trace)) return 0.0;
        const traceState = traceStyleStateForKey(trace.key);
        const visibility = traceVisibilityOpacityMultiplier(trace);
        const presence = clamp01(Number(trace.oviz_presence_opacity ?? 1.0));
        let maximum = 0.0;
        if ((trace.segments && trace.segments.length) || (trace.vectors && trace.vectors.length)) {
          maximum = Math.max(
            maximum,
            (traceState ? clamp01(traceState.opacity) : clamp01(trace.opacity ?? 1.0))
              * visibility
              * presence,
          );
        }
        if (trace.labels && trace.labels.some((label) => label && label.text)) {
          maximum = Math.max(
            maximum,
            (traceState ? clamp01(traceState.opacity) : 1.0) * visibility * presence,
          );
        }
        (Array.isArray(trace.points) ? trace.points : []).forEach((point) => {
          maximum = Math.max(maximum, ovizExpectedPointOpacity(trace, point, displayedTimeMyr));
        });
        return maximum;
      }

      function ovizExpectedPointOpacity(trace, point, displayedTimeMyr) {
          if (!trace || !point || !traceVisible(trace) || !clusterFilterPassesPoint(point)) return 0.0;
          const traceState = traceStyleStateForKey(trace.key);
          const visibility = traceVisibilityOpacityMultiplier(trace);
          const presence = clamp01(Number(trace.oviz_presence_opacity ?? 1.0));
          const pointState = animatedPointState(point, trace, displayedTimeMyr);
          const birthVisibility = ovizPointBirthVisibility(pointState, point, trace);
          const defaultOpacity = Math.max(clamp01(Number(trace.default_opacity ?? 1.0)), 1e-6);
          const traceOpacityMultiplier = traceState
            ? clamp01(traceState.opacity) / defaultOpacity
            : 1.0;
          let selectionOpacity = 1.0;
          const pointKey = clusterFilterSelectionKeyForPoint(point)
            || normalizedSelectionKeyFor(point && point.selection);
          const focusedTraceKey = dendrogramFocusTraceKey();
          const dendrogramActiveKeys = activeDendrogramSelectionKeys();
          if (focusedTraceKey) {
            if (String(trace.key) !== focusedTraceKey) {
              selectionOpacity *= 0.14;
            } else if (dendrogramActiveKeys.size && pointKey) {
              selectionOpacity *= dendrogramActiveKeys.has(pointKey) ? 1.0 : 0.24;
            }
          }
          if (typeof ovizSelectionMembershipOpacity === "function") {
            selectionOpacity *= ovizSelectionMembershipOpacity(pointKey, point);
          } else if (lassoSelectionFilterActive()) {
            selectionOpacity *= pointKey && selectedClusterKeys.has(pointKey) ? 1.0 : 0.0;
          }
          return clamp01(
            Number(pointState.opacity ?? pointOpacityForTrace(point, trace))
            * traceOpacityMultiplier
            * visibility
            * selectionOpacity
            * globalPointOpacityScale
            * presence
            * clamp01(Number(point.oviz_presence_opacity ?? 1.0))
            * birthVisibility
          );
      }

      function ovizRenderedSceneFidelityDifferences(snapshot = {}) {
        const actualOpacityByTrace = new Map();
        const actualOpacityByPoint = new Map();
        let renderedManualLabelCount = 0;
        plotGroup.traverse((object) => {
          if (!object || !object.material || object.visible === false) return;
          const pointMetadata = object.userData && object.userData.ovizRetainedPoint;
          const traceMetadata = object.userData && object.userData.ovizRetainedTrace;
          const labelMetadata = object.userData && object.userData.ovizRetainedLabel;
          if (object.userData && object.userData.ovizRetainedManualLabel) {
            renderedManualLabelCount += 1;
          }
          const traceKey = String(
            (pointMetadata && pointMetadata.traceKey)
            || (traceMetadata && traceMetadata.key)
            || (labelMetadata && labelMetadata.trace && labelMetadata.trace.key)
            || ""
          );
          if (!traceKey) return;
          const opacity = ovizRetainedMaterials(object).reduce(
            (maximum, material) => Math.max(maximum, Number(material && material.opacity) || 0.0),
            0.0,
          );
          actualOpacityByTrace.set(
            traceKey,
            Math.max(Number(actualOpacityByTrace.get(traceKey)) || 0.0, opacity),
          );
          if (pointMetadata) {
            const pointKey = `${traceKey}:${Number(pointMetadata.pointIndex) || 0}`;
            actualOpacityByPoint.set(
              pointKey,
              Math.max(Number(actualOpacityByPoint.get(pointKey)) || 0.0, opacity),
            );
          }
        });
        const displayedTimeMyr = frameTimeForValue(displayedFrameValue);
        const renderedFrame = interpolatedFrameSpecForValue(displayedFrameValue, displayedTimeMyr);
        const expectedOpacityByTrace = {};
        const expectedOpacityByPoint = {};
        const differences = [];
        (renderedFrame && Array.isArray(renderedFrame.traces) ? renderedFrame.traces : []).forEach((trace) => {
          const expectedOpacity = ovizExpectedTraceOpacity(trace, displayedTimeMyr);
          expectedOpacityByTrace[String(trace.key || "")] = expectedOpacity;
          const requiredTraceOpacity = Math.min(0.0001, expectedOpacity * 0.01);
          if (
            expectedOpacity > 0.0001
            && (Number(actualOpacityByTrace.get(String(trace.key || ""))) || 0.0) <= requiredTraceOpacity
          ) {
            differences.push(`rendered_trace:${String(trace.key || trace.name || "unknown")}`);
          } else if (
            expectedOpacity > 0.0001
            && (Number(actualOpacityByTrace.get(String(trace.key || ""))) || 0.0) < expectedOpacity * 0.02
          ) {
            differences.push(`rendered_trace_opacity:${String(trace.key || trace.name || "unknown")}`);
          }
          (Array.isArray(trace.points) ? trace.points : []).forEach((point, pointIndex) => {
            const pointKey = `${String(trace.key || "")}:${pointIndex}`;
            const expectedPointOpacity = ovizExpectedPointOpacity(trace, point, displayedTimeMyr);
            const actualPointOpacity = Number(actualOpacityByPoint.get(pointKey)) || 0.0;
            expectedOpacityByPoint[pointKey] = expectedPointOpacity;
            const requiredPointOpacity = Math.min(0.001, expectedPointOpacity * 0.01);
            if (expectedPointOpacity > 0.001 && actualPointOpacity <= requiredPointOpacity) {
              differences.push(`rendered_point:${pointKey}`);
            } else if (expectedPointOpacity <= 0.001 && actualPointOpacity > 0.001) {
              differences.push(`rendered_selection_extra:${pointKey}`);
            }
          });
        });
        const expectedManualLabelCount = Array.isArray(snapshot.manual_labels)
          ? snapshot.manual_labels.filter((label) => label && label.text).length
          : manualLabels.filter((label) => label && label.text).length;
        if (renderedManualLabelCount !== expectedManualLabelCount) {
          differences.push("rendered_manual_labels");
        }
        volumeRuntimeByKey.forEach((runtime) => {
          const uniforms = runtime && runtime.material && runtime.material.uniforms;
          if (
            uniforms
            && uniforms.selectionTransitionActive
            && uniforms.selectionTransitionActive.value
          ) {
            differences.push(`volume_mask_transition:${String(runtime.layer && runtime.layer.key || "unknown")}`);
          }
        });
        // Exact State restoration must not inherit any Action-owned opacity
        // multiplier.  Such a leftover can make the retained transition look
        // correct, then hide destination points as soon as the final frame is
        // rebuilt.  Report it explicitly instead of letting traceVisible()
        // treat the stale multiplier as part of the expected result.
        if (actionHeldTraceOpacityByKey) {
          differences.push("transient_action_trace_opacity");
        }
        if (Math.abs(Number(renderer.domElement.style.opacity || 1.0) - 1.0) > 1e-6) {
          differences.push("canvas_opacity");
        }
        if (root && root.dataset) {
          root.dataset.renderedTraceFidelity = JSON.stringify({
            exact: differences.length === 0,
            differences,
            expectedOpacityByTrace,
            actualOpacityByTrace: Object.fromEntries(actualOpacityByTrace),
            expectedVisiblePointCount: Object.values(expectedOpacityByPoint).filter((value) => value > 0.001).length,
            actualVisiblePointCount: Array.from(actualOpacityByPoint.values()).filter((value) => value > 0.001).length,
            renderedManualLabelCount,
            expectedManualLabelCount,
          });
        }
        return differences;
      }

      function ovizUpdateRetainedTransitionScene(ownerToken, frameValue, displayedTimeMyr) {
        const frameState = frameValueState(frameValue);
        const intervalKey = `${frameState.lowerIndex}:${frameState.upperIndex}`;
        let runtime = ovizRetainedTransitionScene;
        if (!runtime || runtime.ownerToken !== ownerToken || runtime.intervalKey !== intervalKey) {
          runtime = ovizPrepareRetainedTransitionScene(ownerToken, frameState);
        }
        const alpha = smoothstep01(frameState.alpha);
        const fromWeight = 1.0 - alpha;
        const toWeight = alpha;
        const focusPosition = runtime.fromEndpoint.focusPosition.clone().lerp(
          runtime.toEndpoint.focusPosition,
          alpha,
        );
        runtime.fromEndpoint.root.position.copy(focusPosition);
        runtime.toEndpoint.root.position.copy(focusPosition);
        ovizApplyRetainedSelectionOverlay(runtime, alpha);
        updateTimelineMotionOpacity();
        // Timeline opacity and image-plane camera rules run independently;
        // reapply the endpoint weight after they have updated their base state.
        ovizApplyRetainedEndpointWeight(runtime.fromEndpoint, fromWeight, displayedTimeMyr);
        ovizApplyRetainedEndpointWeight(runtime.toEndpoint, toWeight, displayedTimeMyr);
        runtime.pointPairs.forEach((pair) => {
          const livePoint = ovizRetainedCloneLivePoint(pair, alpha, displayedTimeMyr);
          if (pair.from) {
            const visual = ovizRetainedPointVisual(pair.from, livePoint, displayedTimeMyr);
            ovizApplyRetainedPointEntry(pair.from, visual, fromWeight, livePoint);
          }
          if (pair.to) {
            const visual = ovizRetainedPointVisual(pair.to, livePoint, displayedTimeMyr);
            ovizApplyRetainedPointEntry(pair.to, visual, toWeight, livePoint);
          }
        });
        updateCameraResponsiveImagePlanes();
        runtime.updateCount += 1;
        ovizRetainedSceneUpdateSerial += 1;
        if (root && root.dataset) {
          root.dataset.retainedSceneMetrics = JSON.stringify({
            owner: ownerToken,
            interval: intervalKey,
            builds: ovizRetainedSceneBuildSerial,
            updates: ovizRetainedSceneUpdateSerial,
            intervalUpdates: runtime.updateCount,
            materials: runtime.materialCount,
            pointEntries: runtime.pointPairs.length,
          });
          if (typeof ovizTransitionDebugEnabled === "function" && ovizTransitionDebugEnabled()) {
            root.dataset.retainedSceneDebug = JSON.stringify(ovizRetainedDebugSnapshot(runtime));
          }
        }
        return runtime;
      }

      function renderInterpolatedFrameValue(frameValue, options = {}) {
        displayedFrameValue = clampFrameValue(frameValue);
        currentFrameIndex = clampFrameIndex(displayedFrameValue);
        const displayedTimeMyr = frameTimeForValue(displayedFrameValue);
        updateTimelineUi(displayedFrameValue, displayedTimeMyr);
        const retainedOwner = ovizRetainedTransitionOwner(options);
        if (retainedOwner) {
          ovizUpdateRetainedTransitionScene(retainedOwner, displayedFrameValue, displayedTimeMyr);
          return;
        }
        ovizRetainedTransitionScene = null;
        const frame = interpolatedFrameSpecForValue(displayedFrameValue, displayedTimeMyr);
        renderFrameScene(frame, displayedTimeMyr, options);
      }

      function renderFrame(index) {
        ovizRetainedTransitionScene = null;
        currentFrameIndex = clampFrameIndex(index);
        displayedFrameValue = currentFrameIndex;
        frameTransitionState = null;
        lastTransitionRenderTimestamp = null;
        pendingSliderFrameIndex = null;
        if (sliderScrubRenderHandle !== null) {
          window.cancelAnimationFrame(sliderScrubRenderHandle);
          sliderScrubRenderHandle = null;
        }
        const frame = frameSpecs[currentFrameIndex] || null;
        const displayedTimeMyr = Number(frame && frame.time) || 0.0;
        updateTimelineUi(displayedFrameValue, displayedTimeMyr);
        renderFrameScene(frame, displayedTimeMyr, { updateWidgets: true });
      }

      function updateTopbarDensity() {
        if (!root) {
          return;
        }
        if (!widgetMenuEl) {
          root.dataset.topbarDensity = "regular";
          return;
        }

        const width = Math.max(Number(root.clientWidth) || 0, 1);
        const height = Math.max(Number(root.clientHeight) || 0, 1);
        const aspect = width / Math.max(height, 1);
        const menuWidth = Math.max(
          Number(widgetMenuEl.scrollWidth) || 0,
          Number(widgetMenuEl.getBoundingClientRect().width) || 0
        );
        const titleWidth = titleEl && titleEl.offsetParent !== null
          ? Math.max(Number(titleEl.scrollWidth) || 0, Number(titleEl.getBoundingClientRect().width) || 0)
          : 0;
        const horizontalBudget = Math.max(width - titleWidth - 72.0, 180.0);

        let density = "regular";
        if (width < 780 || height < 560 || aspect < 1.18) {
          density = "stacked";
        } else if (
          width < 1180
          || height < 760
          || aspect < 1.50
          || menuWidth > horizontalBudget
        ) {
          density = "compact";
        }

        if (menuWidth > Math.max(320.0, horizontalBudget * 1.08)) {
          density = "stacked";
        }

        root.dataset.topbarDensity = density;
      }

      function resize() {
        const width = root.clientWidth;
        const height = root.clientHeight;
        updateTopbarDensity();
        renderer.setSize(width, height, false);
        camera.aspect = width / Math.max(height, 1);
        if (typeof applyActionCameraViewOffset === "function") {
          applyActionCameraViewOffset(
            typeof currentActionCameraViewOffset === "object" ? currentActionCameraViewOffset : null
          );
        } else {
          camera.updateProjectionMatrix();
        }
        const retainedOverlayLineMaterials = (
          ovizRetainedTransitionScene
          && Array.isArray(ovizRetainedTransitionScene.overlayLineMaterials)
        ) ? ovizRetainedTransitionScene.overlayLineMaterials : [];
        [...axisLineMaterials, ...frameLineMaterials, ...retainedOverlayLineMaterials].forEach((material) => {
          material.resolution.set(width, height);
        });
        if (legendPanelEl) {
          applyLegendPanelRect(legendPanelRectState || defaultLegendPanelRect());
        }
        if (activeLegendEditorKey && legendEditButtonByKey.has(activeLegendEditorKey)) {
          positionLegendPopover(legendEditButtonByKey.get(activeLegendEditorKey));
        }
        renderTimeSliderTicks();
        updateScaleBar();
        applyScaleBarPosition();
        renderBoxMetricsWidget();
        renderAgeKdeWidget();
        renderClusterFilterWidget();
        renderDendrogramWidget();
        if (typeof updateSkyPanel === "function") {
          updateSkyPanel();
        }
      }

      function animateToFrame(targetIndex, options = {}) {
        const clampedTarget = clampFrameIndex(targetIndex);
        renderFrame(clampedTarget);
      }

      function updateActiveVolumeRuntime(options = {}) {
        const syncControls = options.syncControls !== false;
        const layer = selectedVolumeLayer();
        const state = selectedVolumeState();
        if (!state) {
          if (syncControls) {
            renderVolumeControls();
          }
          return;
        }
        if (layer) {
          clampVolumeStateForLayer(layer, state);
        }
        if (syncControls) {
          renderVolumeControls();
        }
        if (!layer) {
          return;
        }
        const runtime = volumeRuntimeByKey.get(String(layer.key));
        if (runtime) {
          applyVolumeStateToRuntime(layer, runtime);
          if (skySpec.enabled && !currentSelection) {
            updateSkyPanel();
          }
          return;
        }
        const frame = currentFrame();
        if (frame && frameVolumeLayerForStateKey(activeVolumeKey, frame)) {
          renderFrame(currentFrameIndex);
          if (skySpec.enabled && !currentSelection) {
            updateSkyPanel();
          }
        }
      }

      function initVolumeControls() {
        if (
          !volumePanelEl
          || !volumeSelectEl
          || !volumeVisibleEl
          || !volumeColormapEl
          || !volumeStretchEl
          || !volumeVMinEl
          || !volumeVMaxEl
          || !volumeOpacityEl
          || !volumeAlphaEl
          || !volumeStepsEl
        ) {
          return;
        }
        renderVolumeControls();
        if (!volumeStateKeys.length) {
          return;
        }

        volumeSelectEl.addEventListener("change", () => {
          const selectedControlKey = String(volumeSelectEl.value || "");
          if (selectedControlKey.startsWith("variant:")) {
            const variantGroup = selectedControlKey.slice("variant:".length);
            const activeLayer = volumeLayerForKey(activeVolumeKey);
            if (volumeVariantGroupForLayer(activeLayer) !== variantGroup) {
              const nextStateKey = firstVolumeVariantStateKey(variantGroup);
              if (nextStateKey) {
                activeVolumeKey = nextStateKey;
              }
            }
          } else if (selectedControlKey.startsWith("state:")) {
            activeVolumeKey = selectedControlKey.slice("state:".length);
          } else {
            activeVolumeKey = selectedControlKey;
          }
          setExclusiveVolumeVariantSelection(activeVolumeKey);
          updateActiveVolumeRuntime();
        });
        if (volumeSmoothingEl) {
          volumeSmoothingEl.addEventListener("change", () => {
            activeVolumeKey = String(volumeSmoothingEl.value || activeVolumeKey || "");
            setExclusiveVolumeVariantSelection(activeVolumeKey);
            updateActiveVolumeRuntime();
          });
        }
        volumeVisibleEl.addEventListener("change", () => {
          const state = selectedVolumeState();
          if (!state) {
            return;
          }
          state.visible = Boolean(volumeVisibleEl.checked);
          updateActiveVolumeRuntime();
        });
        volumeColormapEl.addEventListener("change", () => {
          const state = selectedVolumeState();
          if (!state) {
            return;
          }
          state.colormap = String(volumeColormapEl.value);
          updateActiveVolumeRuntime();
        });
        volumeStretchEl.addEventListener("change", () => {
          const state = selectedVolumeState();
          if (!state) {
            return;
          }
          state.stretch = normalizeVolumeStretch(volumeStretchEl.value);
          updateActiveVolumeRuntime();
        });
        volumeOpacityEl.addEventListener("input", () => {
          const state = selectedVolumeState();
          if (!state) {
            return;
          }
          state.opacity = Number(volumeOpacityEl.value);
          updateActiveVolumeRuntime();
        });
        volumeAlphaEl.addEventListener("input", () => {
          const state = selectedVolumeState();
          if (!state) {
            return;
          }
          state.alphaCoef = Number(volumeAlphaEl.value);
          updateActiveVolumeRuntime();
        });
        volumeStepsEl.addEventListener("input", () => {
          const state = selectedVolumeState();
          if (!state) {
            return;
          }
          state.steps = Number(volumeStepsEl.value);
          updateActiveVolumeRuntime();
        });
        function updateVolumeWindowFromInput(inputEl, key, options = {}) {
          const state = selectedVolumeState();
          if (!state) {
            return;
          }
          const value = finiteNumberInputValue(inputEl);
          if (value === null) {
            return;
          }
          state[key] = value;
          updateActiveVolumeRuntime(options);
        }
        volumeVMinEl.addEventListener("input", () => {
          updateVolumeWindowFromInput(volumeVMinEl, "vmin", { syncControls: false });
        });
        volumeVMinEl.addEventListener("change", () => {
          updateVolumeWindowFromInput(volumeVMinEl, "vmin");
        });
        volumeVMaxEl.addEventListener("input", () => {
          updateVolumeWindowFromInput(volumeVMaxEl, "vmax", { syncControls: false });
        });
        volumeVMaxEl.addEventListener("change", () => {
          updateVolumeWindowFromInput(volumeVMaxEl, "vmax");
        });
      }
""".strip()
