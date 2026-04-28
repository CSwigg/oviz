from __future__ import annotations


THREEJS_SKY_RUNTIME_JS = """
      function normalizeMemberKey(value) {
        return String(value || "")
          .trim()
          .toLowerCase()
          .replace(/[_\s]+/g, " ");
      }

      const GALACTIC_TO_ICRS_MATRIX = [
        -0.0548755604162154, 0.4941094278755837, -0.8676661490190047,
        -0.8734370902348850, -0.4448296299600112, -0.1980763734312015,
        -0.4838350155487132, 0.7469822444972189, 0.4559837761750669,
      ];

      function normalizeSkyLongitude(value) {
        let lon = Number(value);
        if (!Number.isFinite(lon)) {
          return NaN;
        }
        lon %= 360.0;
        if (lon < 0.0) {
          lon += 360.0;
        }
        return lon;
      }

      function wrapLongitudeDeltaDeg(value) {
        let delta = Number(value);
        if (!Number.isFinite(delta)) {
          return NaN;
        }
        while (delta <= -180.0) {
          delta += 360.0;
        }
        while (delta > 180.0) {
          delta -= 360.0;
        }
        return delta;
      }

      function galacticLonLatDegFromCartesian(x, y, z) {
        const xx = Number(x);
        const yy = Number(y);
        const zz = Number(z);
        if (!Number.isFinite(xx) || !Number.isFinite(yy) || !Number.isFinite(zz)) {
          return null;
        }
        const distance = Math.sqrt(xx * xx + yy * yy + zz * zz);
        if (!(distance > 1e-9)) {
          return null;
        }
        return {
          l: normalizeSkyLongitude(Math.atan2(yy, xx) * 180.0 / Math.PI),
          b: Math.asin(Math.min(1.0, Math.max(-1.0, zz / distance))) * 180.0 / Math.PI,
          distance,
        };
      }

      function applyVolumeSkyAxisTransform(x, y, z, transform) {
        const coords = [Number(x), Number(y), Number(z)];
        if (!coords.every(Number.isFinite)) {
          return null;
        }
        const xyPermutation = Array.isArray(transform && transform.xyPermutation) && transform.xyPermutation.length === 2
          ? transform.xyPermutation
          : [0, 1];
        const xySigns = Array.isArray(transform && transform.xySigns) && transform.xySigns.length === 2
          ? transform.xySigns
          : [1, 1];
        const zSign = Number(transform && transform.zSign);
        const tx = (Number(xySigns[0]) || 1) * coords[Number(xyPermutation[0]) || 0];
        const ty = (Number(xySigns[1]) || 1) * coords[Number(xyPermutation[1]) || 1];
        const tz = (Number.isFinite(zSign) && Math.abs(zSign) > 0.5 ? zSign : 1) * coords[2];
        return { x: tx, y: ty, z: tz };
      }

      function deriveVolumeSkyAxisTransform() {
        const identity = { xyPermutation: [0, 1], xySigns: [1, 1], zSign: 1 };
        const referenceFrame = frameSpecs.find((frame) => approximatelyZero(Number(frame && frame.time))) || frameSpecs[0] || null;
        if (!referenceFrame || !Array.isArray(referenceFrame.traces)) {
          return identity;
        }

        const samples = [];
        referenceFrame.traces.forEach((trace) => {
          if (!trace || !Array.isArray(trace.points)) {
            return;
          }
          trace.points.forEach((point) => {
            const selection = point && typeof point === "object" ? point.selection : null;
            if (!selection || typeof selection !== "object") {
              return;
            }
            const x0 = Number(selection.x0);
            const y0 = Number(selection.y0);
            const z0 = Number(selection.z0);
            const lDeg = Number(selection.l_deg);
            const bDeg = Number(selection.b_deg);
            if (
              !Number.isFinite(x0)
              || !Number.isFinite(y0)
              || !Number.isFinite(z0)
              || !Number.isFinite(lDeg)
              || !Number.isFinite(bDeg)
            ) {
              return;
            }
            const distance = Math.sqrt(x0 * x0 + y0 * y0 + z0 * z0);
            if (!(distance > 1e-6)) {
              return;
            }
            samples.push({ x: x0, y: y0, z: z0, l: lDeg, b: bDeg });
          });
        });

        if (samples.length < 3) {
          return identity;
        }

        const cappedSamples = samples.length > 256
          ? samples.filter((_, idx) => (idx % Math.ceil(samples.length / 256)) === 0)
          : samples;
        const xyCandidates = [
          { xyPermutation: [0, 1], xySigns: [1, 1] },
          { xyPermutation: [0, 1], xySigns: [1, -1] },
          { xyPermutation: [0, 1], xySigns: [-1, 1] },
          { xyPermutation: [0, 1], xySigns: [-1, -1] },
          { xyPermutation: [1, 0], xySigns: [1, 1] },
          { xyPermutation: [1, 0], xySigns: [1, -1] },
          { xyPermutation: [1, 0], xySigns: [-1, 1] },
          { xyPermutation: [1, 0], xySigns: [-1, -1] },
        ];
        const zCandidates = [1, -1];
        let bestTransform = identity;
        let bestScore = Infinity;

        xyCandidates.forEach((xyCandidate) => {
          zCandidates.forEach((zSign) => {
            const transform = {
              xyPermutation: xyCandidate.xyPermutation,
              xySigns: xyCandidate.xySigns,
              zSign,
            };
            let score = 0.0;
            let matchedCount = 0;
            cappedSamples.forEach((sample) => {
              const transformed = applyVolumeSkyAxisTransform(sample.x, sample.y, sample.z, transform);
              const projected = transformed
                ? galacticLonLatDegFromCartesian(transformed.x, transformed.y, transformed.z)
                : null;
              if (!projected) {
                return;
              }
              const lonError = wrapLongitudeDeltaDeg(Number(projected.l) - Number(sample.l));
              const latError = Number(projected.b) - Number(sample.b);
              if (!Number.isFinite(lonError) || !Number.isFinite(latError)) {
                return;
              }
              score += (lonError * lonError) + (latError * latError);
              matchedCount += 1;
            });
            if (!matchedCount) {
              return;
            }
            const normalizedScore = score / matchedCount;
            if (normalizedScore < bestScore) {
              bestScore = normalizedScore;
              bestTransform = transform;
            }
          });
        });

        return bestTransform;
      }

      const volumeSkyAxisTransform = deriveVolumeSkyAxisTransform();

      function volumeSkyGalacticLonLatDegFromCartesian(x, y, z) {
        const transformed = applyVolumeSkyAxisTransform(x, y, z, volumeSkyAxisTransform);
        if (!transformed) {
          return null;
        }
        return galacticLonLatDegFromCartesian(transformed.x, transformed.y, transformed.z);
      }

      function icrsDegFromGalacticDeg(lDeg, bDeg) {
        const lon = Number(lDeg) * Math.PI / 180.0;
        const lat = Number(bDeg) * Math.PI / 180.0;
        if (!Number.isFinite(lon) || !Number.isFinite(lat)) {
          return null;
        }
        const cosLat = Math.cos(lat);
        const xGal = cosLat * Math.cos(lon);
        const yGal = cosLat * Math.sin(lon);
        const zGal = Math.sin(lat);
        const xEq = (
          GALACTIC_TO_ICRS_MATRIX[0] * xGal
          + GALACTIC_TO_ICRS_MATRIX[1] * yGal
          + GALACTIC_TO_ICRS_MATRIX[2] * zGal
        );
        const yEq = (
          GALACTIC_TO_ICRS_MATRIX[3] * xGal
          + GALACTIC_TO_ICRS_MATRIX[4] * yGal
          + GALACTIC_TO_ICRS_MATRIX[5] * zGal
        );
        const zEq = (
          GALACTIC_TO_ICRS_MATRIX[6] * xGal
          + GALACTIC_TO_ICRS_MATRIX[7] * yGal
          + GALACTIC_TO_ICRS_MATRIX[8] * zGal
        );
        const dec = Math.asin(Math.min(1.0, Math.max(-1.0, zEq))) * 180.0 / Math.PI;
        const ra = normalizeSkyLongitude(Math.atan2(yEq, xEq) * 180.0 / Math.PI);
        if (!Number.isFinite(ra) || !Number.isFinite(dec)) {
          return null;
        }
        return { ra, dec };
      }

      function galacticDegFromIcrsDeg(raDeg, decDeg) {
        const ra = Number(raDeg) * Math.PI / 180.0;
        const dec = Number(decDeg) * Math.PI / 180.0;
        if (!Number.isFinite(ra) || !Number.isFinite(dec)) {
          return null;
        }
        const cosDec = Math.cos(dec);
        const xEq = cosDec * Math.cos(ra);
        const yEq = cosDec * Math.sin(ra);
        const zEq = Math.sin(dec);
        const xGal = (
          GALACTIC_TO_ICRS_MATRIX[0] * xEq
          + GALACTIC_TO_ICRS_MATRIX[3] * yEq
          + GALACTIC_TO_ICRS_MATRIX[6] * zEq
        );
        const yGal = (
          GALACTIC_TO_ICRS_MATRIX[1] * xEq
          + GALACTIC_TO_ICRS_MATRIX[4] * yEq
          + GALACTIC_TO_ICRS_MATRIX[7] * zEq
        );
        const zGal = (
          GALACTIC_TO_ICRS_MATRIX[2] * xEq
          + GALACTIC_TO_ICRS_MATRIX[5] * yEq
          + GALACTIC_TO_ICRS_MATRIX[8] * zEq
        );
        const b = Math.asin(Math.min(1.0, Math.max(-1.0, zGal))) * 180.0 / Math.PI;
        const l = normalizeSkyLongitude(Math.atan2(yGal, xGal) * 180.0 / Math.PI);
        if (!Number.isFinite(l) || !Number.isFinite(b)) {
          return null;
        }
        return { l, b };
      }

      function volumeOverlayStretchValue(value, stretchName) {
        const clamped = Math.min(Math.max(Number(value), 0.0), 1.0);
        const stretch = normalizeVolumeStretch(stretchName);
        if (stretch === "log10") {
          const strength = 999.0;
          return Math.log(1.0 + strength * clamped) / Math.log(1.0 + strength);
        }
        if (stretch === "asinh") {
          const strength = 10.0;
          const numer = Math.log(strength * clamped + Math.sqrt((strength * clamped) * (strength * clamped) + 1.0));
          const denom = Math.log(strength + Math.sqrt(strength * strength + 1.0));
          return denom > 0.0 ? numer / denom : clamped;
        }
        return clamped;
      }

      function galacticDirectionVectorFromLonLatDeg(lDeg, bDeg) {
        const lonRad = Number(lDeg) * Math.PI / 180.0;
        const latRad = Number(bDeg) * Math.PI / 180.0;
        if (!Number.isFinite(lonRad) || !Number.isFinite(latRad)) {
          return null;
        }
        const cosLat = Math.cos(latRad);
        return {
          x: cosLat * Math.cos(lonRad),
          y: cosLat * Math.sin(lonRad),
          z: Math.sin(latRad),
        };
      }

      function intersectRayWithBounds(direction, xBounds, yBounds, zBounds) {
        if (!direction) {
          return null;
        }
        const dir = [
          Number(direction.x) || 0.0,
          Number(direction.y) || 0.0,
          Number(direction.z) || 0.0,
        ];
        const bounds = [
          Array.isArray(xBounds) ? xBounds : [-0.5, 0.5],
          Array.isArray(yBounds) ? yBounds : [-0.5, 0.5],
          Array.isArray(zBounds) ? zBounds : [-0.5, 0.5],
        ];
        let tMin = -Infinity;
        let tMax = Infinity;
        for (let axis = 0; axis < 3; axis += 1) {
          const d = dir[axis];
          const low = Number(bounds[axis][0]);
          const high = Number(bounds[axis][1]);
          if (!Number.isFinite(low) || !Number.isFinite(high)) {
            return null;
          }
          if (Math.abs(d) < 1e-8) {
            if (0.0 < Math.min(low, high) || 0.0 > Math.max(low, high)) {
              return null;
            }
            continue;
          }
          let t0 = low / d;
          let t1 = high / d;
          if (t0 > t1) {
            const swap = t0;
            t0 = t1;
            t1 = swap;
          }
          tMin = Math.max(tMin, t0);
          tMax = Math.min(tMax, t1);
          if (!(tMax > tMin)) {
            return null;
          }
        }
        const near = Math.max(tMin, 0.0);
        if (!(tMax > near)) {
          return null;
        }
        return { tMin: near, tMax };
      }

      function sampleVolumeScalarTrilinear(scalarData, nx, ny, nz, xNorm, yNorm, zNorm) {
        if (!scalarData || !scalarData.length) {
          return 0.0;
        }
        const x = Math.min(Math.max(Number(xNorm), 0.0), 1.0) * Math.max(nx - 1, 0);
        const y = Math.min(Math.max(Number(yNorm), 0.0), 1.0) * Math.max(ny - 1, 0);
        const z = Math.min(Math.max(Number(zNorm), 0.0), 1.0) * Math.max(nz - 1, 0);
        const x0 = Math.floor(x);
        const y0 = Math.floor(y);
        const z0 = Math.floor(z);
        const x1 = Math.min(x0 + 1, nx - 1);
        const y1 = Math.min(y0 + 1, ny - 1);
        const z1 = Math.min(z0 + 1, nz - 1);
        const tx = x - x0;
        const ty = y - y0;
        const tz = z - z0;
        const strideY = nx;
        const strideZ = nx * ny;
        const index000 = z0 * strideZ + y0 * strideY + x0;
        const index100 = z0 * strideZ + y0 * strideY + x1;
        const index010 = z0 * strideZ + y1 * strideY + x0;
        const index110 = z0 * strideZ + y1 * strideY + x1;
        const index001 = z1 * strideZ + y0 * strideY + x0;
        const index101 = z1 * strideZ + y0 * strideY + x1;
        const index011 = z1 * strideZ + y1 * strideY + x0;
        const index111 = z1 * strideZ + y1 * strideY + x1;
        const c000 = Number(scalarData[index000] || 0) / 255.0;
        const c100 = Number(scalarData[index100] || 0) / 255.0;
        const c010 = Number(scalarData[index010] || 0) / 255.0;
        const c110 = Number(scalarData[index110] || 0) / 255.0;
        const c001 = Number(scalarData[index001] || 0) / 255.0;
        const c101 = Number(scalarData[index101] || 0) / 255.0;
        const c011 = Number(scalarData[index011] || 0) / 255.0;
        const c111 = Number(scalarData[index111] || 0) / 255.0;
        const c00 = c000 * (1.0 - tx) + c100 * tx;
        const c10 = c010 * (1.0 - tx) + c110 * tx;
        const c01 = c001 * (1.0 - tx) + c101 * tx;
        const c11 = c011 * (1.0 - tx) + c111 * tx;
        const c0 = c00 * (1.0 - ty) + c10 * ty;
        const c1 = c01 * (1.0 - ty) + c11 * ty;
        return c0 * (1.0 - tz) + c1 * tz;
      }

      function pointInsideProjectedLassoMask(worldX, worldY, worldZ, mask) {
        if (
          !mask
          || !mask.viewProjectionMatrix
          || !Array.isArray(mask.polygonNdc)
          || mask.polygonNdc.length < 3
        ) {
          return false;
        }
        const e = mask.viewProjectionMatrix.elements || [];
        if (e.length !== 16) {
          return false;
        }
        const clipX = e[0] * worldX + e[4] * worldY + e[8] * worldZ + e[12];
        const clipY = e[1] * worldX + e[5] * worldY + e[9] * worldZ + e[13];
        const clipZ = e[2] * worldX + e[6] * worldY + e[10] * worldZ + e[14];
        const clipW = e[3] * worldX + e[7] * worldY + e[11] * worldZ + e[15];
        if (!(clipW > 0.0)) {
          return false;
        }
        const ndcPoint = {
          x: clipX / clipW,
          y: clipY / clipW,
          z: clipZ / clipW,
        };
        if (!Number.isFinite(ndcPoint.x) || !Number.isFinite(ndcPoint.y) || !Number.isFinite(ndcPoint.z)) {
          return false;
        }
        if (ndcPoint.z < -1.0 || ndcPoint.z > 1.0) {
          return false;
        }
        if (ndcPoint.x < -1.0 || ndcPoint.x > 1.0 || ndcPoint.y < -1.0 || ndcPoint.y > 1.0) {
          return false;
        }
        const maskSize = Math.max(0, Math.round(Number(mask.maskSize) || 0));
        const maskAlphaData = mask.maskAlphaData;
        if (maskSize > 0 && maskAlphaData && maskAlphaData.length >= maskSize * maskSize * 4) {
          const maskX = Math.max(0, Math.min(maskSize - 1, Math.round((ndcPoint.x * 0.5 + 0.5) * (maskSize - 1))));
          const maskY = Math.max(0, Math.min(maskSize - 1, Math.round((0.5 - ndcPoint.y * 0.5) * (maskSize - 1))));
          const alphaIndex = ((maskY * maskSize) + maskX) * 4 + 3;
          return Number(maskAlphaData[alphaIndex] || 0) > 0;
        }
        return pointInPolygon(ndcPoint, mask.polygonNdc);
      }

      function buildVolumeSkyImageOverlaySpec(mode = "overview") {
        if (mode === "click") {
          return null;
        }
        const activeMask = activeVolumeLassoSelectionMask();
        if (
          !activeMask
          || !activeMask.viewProjectionMatrix
          || !Array.isArray(activeMask.polygonNdc)
          || activeMask.polygonNdc.length < 3
        ) {
          return null;
        }

        const frame = currentFrame();
        const frameTime = frame ? Number(frame.time) : NaN;
        const displayOffsetX = Number(plotGroup.position.x) || 0.0;
        const displayOffsetY = Number(plotGroup.position.y) || 0.0;
        const displayOffsetZ = Number(plotGroup.position.z) || 0.0;
        const collectedSamples = [];
        const usedLayerNames = [];
        let usedLayerCount = 0;

        frameVolumeLayers(frame).forEach((layer) => {
          const stateKey = volumeStateKeyForLayer(layer);
          const state = volumeStateByKey[stateKey];
          if (!state || state.visible === false || legendState[stateKey] === false) {
            return;
          }

          const shape = layer.sky_overlay_shape || layer.shape || {};
          const nx = Math.max(1, Math.round(Number(shape.x) || 0));
          const ny = Math.max(1, Math.round(Number(shape.y) || 0));
          const nz = Math.max(1, Math.round(Number(shape.z) || 0));
          const totalVoxels = nx * ny * nz;
          if (!(totalVoxels > 0)) {
            return;
          }

          const scalarData = volumeSkyScalarArrayFor(layer);
          if (!scalarData || !scalarData.length) {
            return;
          }

          const option = volumeColormapOptionFor(layer, state.colormap);
          const colorBytes = option ? volumeColorBytesForOption(option) : null;
          if (!colorBytes || colorBytes.length < 4) {
            return;
          }

          const bounds = layer.bounds || {};
          const xBounds = Array.isArray(bounds.x) ? bounds.x : [-0.5, 0.5];
          const yBounds = Array.isArray(bounds.y) ? bounds.y : [-0.5, 0.5];
          const zBounds = Array.isArray(bounds.z) ? bounds.z : [-0.5, 0.5];
          const xSpan = Number(xBounds[1]) - Number(xBounds[0]);
          const ySpan = Number(yBounds[1]) - Number(yBounds[0]);
          const zSpan = Number(zBounds[1]) - Number(zBounds[0]);
          if (!(Number.isFinite(xSpan) && Number.isFinite(ySpan) && Number.isFinite(zSpan))) {
            return;
          }

          const windowState = normalizedVolumeWindowFor(layer, state);
          const low = Number(windowState.low);
          const high = Number(windowState.high);
          const span = Math.max(high - low, 1e-6);
          const desiredScanCount = 2500000;
          const stride = Math.max(1, Math.ceil(Math.cbrt(totalVoxels / desiredScanCount)));
          const minScaledThreshold = 0.02;
          const layerOpacityWeight = Math.max(0.05, Number(state.opacity) || 0.0);
          const lutSamples = Math.max(1, Math.floor(colorBytes.length / 4));
          const cellSizeX = Math.abs(xSpan / Math.max(nx, 1));
          const cellSizeY = Math.abs(ySpan / Math.max(ny, 1));
          const cellSizeZ = Math.abs(zSpan / Math.max(nz, 1));
          const cellSizeMin = Math.max(1e-6, Math.min(cellSizeX, cellSizeY, cellSizeZ));
          let layerUsed = false;

          for (let iz = 0; iz < nz; iz += stride) {
            const localZ = Number(zBounds[0]) + ((iz + 0.5) / nz) * zSpan;
            const displayedZ = localZ + displayOffsetZ;
            for (let iy = 0; iy < ny; iy += stride) {
              const localY = Number(yBounds[0]) + ((iy + 0.5) / ny) * ySpan;
              const displayedY = localY + displayOffsetY;
              for (let ix = 0; ix < nx; ix += stride) {
                const dataIndex = iz * nx * ny + iy * nx + ix;
                const normalizedValue = Number(scalarData[dataIndex] || 0) / 255.0;
                if (!(normalizedValue > low)) {
                  continue;
                }
                const scaledValue = Math.min(Math.max((normalizedValue - low) / span, 0.0), 1.0);
                const stretchedValue = volumeOverlayStretchValue(scaledValue, state.stretch);
                if (!(stretchedValue > minScaledThreshold)) {
                  continue;
                }

                const localX = Number(xBounds[0]) + ((ix + 0.5) / nx) * xSpan;
                const displayedX = localX + displayOffsetX;
                if (!pointInsideProjectedLassoMask(displayedX, displayedY, displayedZ, activeMask)) {
                  continue;
                }

                const galactic = volumeSkyGalacticLonLatDegFromCartesian(localX, localY, localZ);
                if (!galactic) {
                  continue;
                }

                const lDeg = normalizeSkyLongitude(galactic.l);
                const bDeg = Number(galactic.b);
                if (!Number.isFinite(lDeg) || !Number.isFinite(bDeg)) {
                  continue;
                }

                const weight = stretchedValue * layerOpacityWeight;
                const icrs = icrsDegFromGalacticDeg(lDeg, bDeg);
                if (!icrs) {
                  continue;
                }
                const colorIndex = Math.max(
                  0,
                  Math.min(
                    lutSamples - 1,
                    Math.round(stretchedValue * (lutSamples - 1))
                  )
                ) * 4;
                const angularFootprintDeg = Math.atan2(cellSizeMin * 0.9, Math.max(Number(galactic.distance) || 0.0, 1e-6)) * 180.0 / Math.PI;
                collectedSamples.push({
                  l: lDeg,
                  b: bDeg,
                  ra: Number(icrs.ra),
                  dec: Number(icrs.dec),
                  weight,
                  r: Number(colorBytes[colorIndex]) || 0,
                  g: Number(colorBytes[colorIndex + 1]) || 0,
                  bColor: Number(colorBytes[colorIndex + 2]) || 0,
                  sigmaDeg: Math.min(Math.max(angularFootprintDeg * 0.85, 0.02), 6.0),
                });
                layerUsed = true;
              }
            }
          }

          if (layerUsed) {
            usedLayerCount += 1;
            const baseName = volumeBaseNameForLayer(layer) || volumeStateNameForLayer(layer);
            if (baseName && !usedLayerNames.includes(baseName)) {
              usedLayerNames.push(baseName);
            }
          }
        });

        if (!collectedSamples.length || !usedLayerCount) {
          return null;
        }

        let sumLonSin = 0.0;
        let sumLonCos = 0.0;
        let sumLat = 0.0;
        let sumWeight = 0.0;
        collectedSamples.forEach((sample) => {
          const weight = Math.max(Number(sample.weight) || 0.0, 1e-6);
          const lonRad = Number(sample.l) * Math.PI / 180.0;
          sumLonSin += Math.sin(lonRad) * weight;
          sumLonCos += Math.cos(lonRad) * weight;
          sumLat += Number(sample.b) * weight;
          sumWeight += weight;
        });

        const baseCenterLon = normalizeSkyLongitude(Math.atan2(sumLonSin, sumLonCos) * 180.0 / Math.PI);
        const baseCenterLat = sumWeight > 0.0 ? (sumLat / sumWeight) : 0.0;
        let minRelLon = Infinity;
        let maxRelLon = -Infinity;
        let minLat = Infinity;
        let maxLat = -Infinity;
        collectedSamples.forEach((sample) => {
          const relLon = wrapLongitudeDeltaDeg(Number(sample.l) - baseCenterLon);
          sample.relLonBase = relLon;
          minRelLon = Math.min(minRelLon, relLon);
          maxRelLon = Math.max(maxRelLon, relLon);
          minLat = Math.min(minLat, Number(sample.b));
          maxLat = Math.max(maxLat, Number(sample.b));
        });

        const lonSpanRaw = Math.max(maxRelLon - minRelLon, 0.25);
        const latSpanRaw = Math.max(maxLat - minLat, 0.25);
        const lonSpan = Math.min(Math.max(lonSpanRaw * 1.12, 2.0), 180.0);
        const latSpan = Math.min(Math.max(latSpanRaw * 1.12, 2.0), 180.0);
        const lonMidOffset = 0.5 * (minRelLon + maxRelLon);
        const overlayCenterLon = normalizeSkyLongitude(baseCenterLon + lonMidOffset);
        const overlayCenterLat = Math.min(Math.max(baseCenterLat + 0.5 * ((minLat + maxLat) - (2.0 * baseCenterLat)), -89.0), 89.0);
        const aspect = Math.max(latSpan / lonSpan, 0.35);
        const width = lonSpan <= 10.0 ? 2048 : (lonSpan <= 28.0 ? 1536 : (lonSpan <= 70.0 ? 1024 : 768));
        const height = Math.max(384, Math.min(1536, Math.round(width * aspect)));
        const pixelCount = width * height;
        const weightGrid = new Float32Array(pixelCount);
        const redGrid = new Float32Array(pixelCount);
        const greenGrid = new Float32Array(pixelCount);
        const blueGrid = new Float32Array(pixelCount);
        const halfLonSpan = lonSpan * 0.5;
        const halfLatSpan = latSpan * 0.5;
        const pixelScaleLon = lonSpan / Math.max(width, 1);
        const pixelScaleLat = latSpan / Math.max(height, 1);
        const pixelScale = Math.max(Math.min(pixelScaleLon, pixelScaleLat), 1e-6);
        const maxSamplesForOverlay = 120000;
        const overlaySamples = collectedSamples.length > maxSamplesForOverlay
          ? collectedSamples.filter((_, idx) => (idx % Math.ceil(collectedSamples.length / maxSamplesForOverlay)) === 0)
          : collectedSamples;

        overlaySamples.forEach((sample) => {
          const relLon = wrapLongitudeDeltaDeg(Number(sample.l) - overlayCenterLon);
          const relLat = Number(sample.b) - overlayCenterLat;
          if (Math.abs(relLon) > halfLonSpan || Math.abs(relLat) > halfLatSpan) {
            return;
          }
          const xCenter = ((halfLonSpan - relLon) / lonSpan) * width - 0.5;
          const yCenter = ((halfLatSpan - relLat) / latSpan) * height - 0.5;
          const sigmaPx = Math.max(1.1, Number(sample.sigmaDeg || pixelScale) / pixelScale);
          const radiusPx = Math.max(2, Math.min(18, Math.ceil(sigmaPx * 2.8)));
          const xMin = Math.max(0, Math.floor(xCenter - radiusPx));
          const xMax = Math.min(width - 1, Math.ceil(xCenter + radiusPx));
          const yMin = Math.max(0, Math.floor(yCenter - radiusPx));
          const yMax = Math.min(height - 1, Math.ceil(yCenter + radiusPx));
          const baseWeight = Math.max(Number(sample.weight) || 0.0, 0.0);
          if (!(baseWeight > 0.0)) {
            return;
          }
          const sigmaSq = sigmaPx * sigmaPx;
          for (let yIndex = yMin; yIndex <= yMax; yIndex += 1) {
            const dy = yIndex - yCenter;
            for (let xIndex = xMin; xIndex <= xMax; xIndex += 1) {
              const dx = xIndex - xCenter;
              const distanceSq = dx * dx + dy * dy;
              const kernelWeight = Math.exp(-0.5 * distanceSq / sigmaSq);
              if (!(kernelWeight > 1e-4)) {
                continue;
              }
              const contribution = baseWeight * kernelWeight;
              const pixelIndex = yIndex * width + xIndex;
              weightGrid[pixelIndex] += contribution;
              redGrid[pixelIndex] += Number(sample.r || 0) * contribution;
              greenGrid[pixelIndex] += Number(sample.g || 0) * contribution;
              blueGrid[pixelIndex] += Number(sample.bColor || 0) * contribution;
            }
          }
        });

        let maxWeight = 0.0;
        let nonZeroPixels = 0;
        for (let index = 0; index < pixelCount; index += 1) {
          const weight = Number(weightGrid[index]);
          if (weight > 0.0) {
            nonZeroPixels += 1;
            if (weight > maxWeight) {
              maxWeight = weight;
            }
          }
        }

        if (!(maxWeight > 0.0) || !nonZeroPixels) {
          return null;
        }

        const overlayCanvas = document.createElement("canvas");
        overlayCanvas.width = width;
        overlayCanvas.height = height;
        const overlayCtx = overlayCanvas.getContext("2d");
        if (!overlayCtx) {
          return null;
        }
        const imageData = overlayCtx.createImageData(width, height);
        const rgba = imageData.data;
        const normalizer = maxWeight > 0.0 ? maxWeight : 1.0;

        for (let index = 0; index < pixelCount; index += 1) {
          const weight = Number(weightGrid[index]);
          if (!(weight > 0.0)) {
            continue;
          }
          const normalized = Math.min(weight / normalizer, 1.0);
          const intensity = Math.pow(normalized, 0.72);
          const outIndex = index * 4;
          const invWeight = 1.0 / weight;
          rgba[outIndex] = Math.max(
            0,
            Math.min(255, Math.round((Number(redGrid[index]) * invWeight) * (0.30 + 0.70 * intensity)))
          );
          rgba[outIndex + 1] = Math.max(
            0,
            Math.min(255, Math.round((Number(greenGrid[index]) * invWeight) * (0.30 + 0.70 * intensity)))
          );
          rgba[outIndex + 2] = Math.max(
            0,
            Math.min(255, Math.round((Number(blueGrid[index]) * invWeight) * (0.30 + 0.70 * intensity)))
          );
          rgba[outIndex + 3] = Math.max(0, Math.min(255, Math.round(255.0 * Math.pow(normalized, 0.78))));
        }

        overlayCtx.putImageData(imageData, 0, 0);
        const overlayName = usedLayerNames.length === 1
          ? `Selected ${usedLayerNames[0]}`
          : (usedLayerCount > 1 ? "Selected Volume Intensity" : "Selected Volume");
        return {
          kind: "volume_image",
          name: overlayName,
          data_url: overlayCanvas.toDataURL("image/png"),
          width,
          height,
          sample_count: collectedSamples.length,
          non_zero_pixels: nonZeroPixels,
          wcs: {
            NAXIS: 2,
            NAXIS1: width,
            NAXIS2: height,
            CTYPE1: "GLON-CAR",
            CTYPE2: "GLAT-CAR",
            CUNIT1: "deg",
            CUNIT2: "deg",
            CRPIX1: (width / 2.0) + 0.5,
            CRPIX2: (height / 2.0) + 0.5,
            CRVAL1: overlayCenterLon,
            CRVAL2: overlayCenterLat,
            CDELT1: -lonSpan / width,
            // The PNG is written in browser image coordinates with a top-left origin,
            // so latitude should decrease as rows move downward.
            CDELT2: latSpan / height,
            CROTA2: 0.0,
            LONPOLE: 180.0,
            LATPOLE: 90.0,
          },
        };
      }

      function resolveMemberPoints(selection) {
        const byCluster = skySpec.members_by_cluster || {};
        const traceName = selection && selection.trace_name ? String(selection.trace_name) : "";
        const clusterName = selection && selection.cluster_name ? String(selection.cluster_name) : "";
        const candidates = [clusterName, traceName]
          .filter(Boolean)
          .flatMap((name) => [name, name.replace(/_/g, " "), name.replace(/\s+/g, "_")]);

        for (const key of candidates) {
          if (Array.isArray(byCluster[key]) && byCluster[key].length) {
            return { key, points: byCluster[key] };
          }
        }

        const normalizedCandidates = new Set(candidates.map((name) => normalizeMemberKey(name)));
        for (const key of Object.keys(byCluster)) {
          if (normalizedCandidates.has(normalizeMemberKey(key)) && Array.isArray(byCluster[key]) && byCluster[key].length) {
            return { key, points: byCluster[key] };
          }
        }

        return { key: traceName || clusterName || "Selection", points: null };
      }

      function buildAladinCatalogPayload(selections, mode = "overview") {
        const activeSelections = uniqueSelections(selections);
        const payload = activeSelections.map((selection) => {
          const resolvedMembers = resolveMemberPoints(selection);
          const lookupName = resolvedMembers.key || selectionKeyFor(selection) || "Selection";
          const traceName = selection && selection.trace_name ? String(selection.trace_name) : lookupName;
          const clusterKey = normalizedSelectionKeyFor(selection);
          const clusterColor = selection.cluster_color ? String(selection.cluster_color) : "#ffffff";
          const members = Array.isArray(resolvedMembers.points) ? resolvedMembers.points : null;
          let points = [];
          if (members && members.length) {
            points = members.map((pt) => ({
              l: Number(pt.l),
              b: Number(pt.b),
              ra: Number(pt.ra),
              dec: Number(pt.dec),
              label: pt.label || lookupName,
              clusterKey,
            }));
          } else if (Number.isFinite(Number(selection.ra_deg)) && Number.isFinite(Number(selection.dec_deg))) {
            points = [{
              l: Number(selection.l_deg),
              b: Number(selection.b_deg),
              ra: Number(selection.ra_deg),
              dec: Number(selection.dec_deg),
              label: lookupName,
              clusterKey,
            }];
          }
          return {
            name: lookupName,
            traceName,
            color: clusterColor,
            opacity: 1.0,
            sourceSize: points.length > 1 ? 4 : 7,
            points,
          };
        }).filter((catalog) => (catalog.points || []).length);

        if (mode === "click" || payload.length <= 1) {
          return payload;
        }

        const grouped = new Map();
        payload.forEach((catalog) => {
          const groupKey = String(catalog.traceName || catalog.name || "Selected trace");
          if (!grouped.has(groupKey)) {
            grouped.set(groupKey, {
              name: groupKey,
              color: catalog.color,
              opacity: 1.0,
              sourceSize: 4,
              points: [],
              seen: new Set(),
            });
          }
          const group = grouped.get(groupKey);
          (catalog.points || []).forEach((point) => {
            const ra = Number(point.ra);
            const dec = Number(point.dec);
            const label = String(point.label || catalog.name || "Selection");
            const clusterKey = String(point.clusterKey || "");
            const key = `${ra.toFixed(8)}|${dec.toFixed(8)}|${label}|${clusterKey}`;
            if (!Number.isFinite(ra) || !Number.isFinite(dec) || group.seen.has(key)) {
              return;
            }
            group.seen.add(key);
            group.points.push({
              l: Number(point.l),
              b: Number(point.b),
              ra,
              dec,
              label,
              clusterKey,
            });
          });
        });

        return Array.from(grouped.values())
          .map((group) => ({
            name: group.name,
            color: group.color,
            opacity: group.opacity,
            sourceSize: group.sourceSize,
            points: group.points,
          }))
          .filter((group) => group.points.length);
      }

      function angularSeparationDeg(ra1Deg, dec1Deg, ra2Deg, dec2Deg) {
        const rad = Math.PI / 180.0;
        const sin1 = Math.sin(dec1Deg * rad);
        const sin2 = Math.sin(dec2Deg * rad);
        const cos1 = Math.cos(dec1Deg * rad);
        const cos2 = Math.cos(dec2Deg * rad);
        const deltaRa = (ra1Deg - ra2Deg) * rad;
        const cosSep = Math.min(1.0, Math.max(-1.0, sin1 * sin2 + cos1 * cos2 * Math.cos(deltaRa)));
        return Math.acos(cosSep) * 180.0 / Math.PI;
      }

      function skyFocusFromPayload(selections, catalogPayload) {
        const focusPoints = [];
        (catalogPayload || []).forEach((catalog) => {
          (catalog.points || []).forEach((point) => {
            const ra = Number(point.ra);
            const dec = Number(point.dec);
            if (Number.isFinite(ra) && Number.isFinite(dec)) {
              focusPoints.push({ ra, dec });
            }
          });
        });
        if (!focusPoints.length) {
          uniqueSelections(selections).forEach((selection) => {
            const ra = Number(selection.ra_deg);
            const dec = Number(selection.dec_deg);
            if (Number.isFinite(ra) && Number.isFinite(dec)) {
              focusPoints.push({ ra, dec });
            }
          });
        }
        if (!focusPoints.length) {
          return null;
        }

        const rad = Math.PI / 180.0;
        let sx = 0.0;
        let sy = 0.0;
        let sz = 0.0;
        focusPoints.forEach((point) => {
          const ra = point.ra * rad;
          const dec = point.dec * rad;
          const cosDec = Math.cos(dec);
          sx += cosDec * Math.cos(ra);
          sy += cosDec * Math.sin(ra);
          sz += Math.sin(dec);
        });

        const norm = Math.sqrt(sx * sx + sy * sy + sz * sz);
        let centerRa = focusPoints[0].ra;
        let centerDec = focusPoints[0].dec;
        if (norm > 1e-9) {
          centerRa = Math.atan2(sy, sx) * 180.0 / Math.PI;
          if (centerRa < 0.0) {
            centerRa += 360.0;
          }
          centerDec = Math.asin(Math.min(1.0, Math.max(-1.0, sz / norm))) * 180.0 / Math.PI;
        }

        let maxSep = 0.0;
        focusPoints.forEach((point) => {
          maxSep = Math.max(maxSep, angularSeparationDeg(centerRa, centerDec, point.ra, point.dec));
        });

        const radiusDeg = Number(skySpec.radius_deg || 1.0);
        return {
          ra: centerRa,
          dec: centerDec,
          fovDeg: Math.min(Math.max(radiusDeg * 2.4, maxSep * 2.8, 1.2), 180.0),
        };
      }
""".strip()
