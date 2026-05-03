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

      function volumeScalarArrayFor(layer) {
        const layerKey = String(layer.key);
        if (volumeScalarDataCache.has(layerKey)) {
          return volumeScalarDataCache.get(layerKey);
        }
        const encoding = String(layer.data_encoding || "uint16_le");
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
              if (!insideSelection && selectionDimOutside <= 0.001) {
                continue;
              }
              pxColor.a = 1.0 - pow(1.0 - clamp(pxColor.a * opacity, 0.0, 0.999), step * alpha_coef);
              if (!insideSelection) {
                pxColor.a *= selectionDimOutside;
              }
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
          },
        ]);
        applyLassoSelectionMaskUniforms(uniforms, activeVolumeLassoSelectionMask());

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
        applyLassoSelectionMaskUniforms(runtime.material.uniforms, activeVolumeLassoSelectionMask());
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
        return blended;
      }

      function interpolateTraceSpec(traceA, traceB, alpha, timeValue) {
        const fallbackTrace = alpha < 0.5 ? traceA : traceB;
        if (!traceA || !traceB) {
          return fallbackTrace || null;
        }
        const blended = Object.assign({}, fallbackTrace);
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

        if (Array.isArray(traceA.points) && Array.isArray(traceB.points) && traceA.points.length === traceB.points.length) {
          blended.points = traceA.points.map((pointA, index) => (
            interpolateTracePoint(pointA, traceB.points[index], alpha, timeValue)
          ));
        }

        if (Array.isArray(traceA.labels) && Array.isArray(traceB.labels) && traceA.labels.length === traceB.labels.length) {
          blended.labels = traceA.labels.map((labelA, index) => (
            interpolateTraceLabel(labelA, traceB.labels[index], alpha)
          ));
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
        const upperTraceByKey = new Map(upperTraces.map((trace) => [String(trace && trace.key), trace]));
        const usedUpperKeys = new Set();
        const blendedTraces = [];

        lowerTraces.forEach((lowerTrace, index) => {
          const traceKey = String(lowerTrace && lowerTrace.key);
          let upperTrace = traceKey ? upperTraceByKey.get(traceKey) : null;
          if (!upperTrace && index < upperTraces.length) {
            upperTrace = upperTraces[index];
          }
          if (traceKey && upperTrace) {
            usedUpperKeys.add(traceKey);
          }
          const blendedTrace = interpolateTraceSpec(lowerTrace, upperTrace, easedAlpha, timeValue);
          if (blendedTrace) {
            blendedTraces.push(blendedTrace);
          }
        });

        upperTraces.forEach((upperTrace, index) => {
          const traceKey = String(upperTrace && upperTrace.key);
          if ((traceKey && usedUpperKeys.has(traceKey)) || index < lowerTraces.length) {
            return;
          }
          if (upperTrace) {
            blendedTraces.push(upperTrace);
          }
        });

        const decorationFrame = easedAlpha < 0.5 ? lowerFrame : upperFrame;
        return {
          name: formatTick(timeValue),
          time: timeValue,
          traces: blendedTraces,
          decorations: decorationFrame && Array.isArray(decorationFrame.decorations)
            ? decorationFrame.decorations
            : [],
        };
      }

      function renderFrameScene(frame, displayedTimeMyr, options = {}) {
        const updateWidgets = options.updateWidgets !== false;
        if (!frame) {
          return;
        }
        currentZoomAnchorPoint = trackedZoomAnchorPointForFrame(frame);
        if (galacticSimpleOrbitTargetTrackingActive && currentZoomAnchorPoint instanceof THREE.Vector3) {
          const orbitTargetDelta = currentZoomAnchorPoint.clone().sub(controls.target);
          if (orbitTargetDelta.lengthSq() > 1e-12) {
            controls.target.add(orbitTargetDelta);
            camera.position.add(orbitTargetDelta);
          }
        }
        tooltipEl.style.display = "none";
        hoverTargets.length = 0;
        cameraResponsivePointEntries.length = 0;
        cameraResponsiveImagePlaneEntries.length = 0;
        selectionSpriteEntriesByKey.clear();
        screenStableTextSprites.length = 0;
        clearGroup(plotGroup);
        plotGroup.position.set(0.0, 0.0, 0.0);
        frameLineMaterials.length = 0;
        volumeRuntimeByKey.clear();

        frame.traces.forEach((trace) => {
          if (!traceVisible(trace)) {
            return;
          }

          if (trace.segments && trace.segments.length) {
            const line = makeLineObject(trace, frameLineMaterials);
            if (line) {
              plotGroup.add(line);
            }
          }
          if (trace.points && trace.points.length) {
            addMarkerTrace(plotGroup, trace);
          }
          if (trace.labels && trace.labels.length) {
            addTextTrace(plotGroup, trace);
          }
        });

        if (currentSelectionMode === "click" && currentSelection && approximatelyZero(displayedTimeMyr)) {
          const footprint = buildSelectionFootprint(currentSelection, frameLineMaterials);
          if (footprint) {
            plotGroup.add(footprint);
          }
        }

        (frame.decorations || []).forEach((decoration) => {
          addDecoration(plotGroup, decoration);
        });

        const selectionBoxGroup = buildSelectionBoxGroupForFrame(displayedTimeMyr);
        if (selectionBoxGroup) {
          plotGroup.add(selectionBoxGroup);
        }

        addManualLabels(plotGroup);

        const focusOffset = focusTrackingOffsetForFrame(frame);
        if (focusOffset) {
          plotGroup.position.copy(focusOffset.multiplyScalar(-1.0));
        }

        applySceneHoverState();
        if (updateWidgets) {
          resize();
          renderVolumeControls();
          renderBoxMetricsWidget();
          renderAgeKdeWidget();
          renderClusterFilterWidget();
          renderDendrogramWidget();
        }
      }

      function renderInterpolatedFrameValue(frameValue, options = {}) {
        displayedFrameValue = clampFrameValue(frameValue);
        currentFrameIndex = clampFrameIndex(displayedFrameValue);
        const displayedTimeMyr = frameTimeForValue(displayedFrameValue);
        const frame = interpolatedFrameSpecForValue(displayedFrameValue, displayedTimeMyr);
        updateTimelineUi(displayedFrameValue, displayedTimeMyr);
        renderFrameScene(frame, displayedTimeMyr, options);
      }

      function renderFrame(index) {
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
        [...axisLineMaterials, ...frameLineMaterials].forEach((material) => {
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
