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

      function skyDomeIsEnabled() {
        return Boolean(skyDomeSpec && skyDomeSpec.enabled);
      }

      function skyDomeImageDataUrl() {
        const dataUrl = skyDomeSpec && skyDomeSpec.image_data_url;
        return typeof dataUrl === "string" && dataUrl.startsWith("data:image/") ? dataUrl : "";
      }

      function normalizeSkyDomeProjectionName(projectionName) {
        const normalized = String(projectionName || "").trim().toUpperCase();
        if (normalized === "MOL" || normalized === "MOLLWEIDE") {
          return "MOL";
        }
        return "CAR";
      }

      function skyDomeProjectionMode(projectionName) {
        return normalizeSkyDomeProjectionName(projectionName) === "MOL" ? 0 : 1;
      }

      function skyDomeTextureCoordinateFrameName() {
        const frame = String(
          (skyDomeSpec && (
            skyDomeSpec.texture_coordinate_frame
            || skyDomeSpec.texture_frame
            || skyDomeSpec.hips2fits_coordsys
            || skyDomeSpec.hips2fits_coordinate_system
            || skyDomeSpec.coordinate_frame
            || (skyDomeSpec.projection_metadata && skyDomeSpec.projection_metadata.coordinate_frame)
          ))
          || "galactic"
        ).trim().toLowerCase();
        return frame === "icrs" || frame === "equatorial" ? "icrs" : "galactic";
      }

      function skyDomeTextureFrameMode() {
        return skyDomeTextureCoordinateFrameName() === "icrs" ? 1 : 0;
      }

      function skyDomeGalacticToIcrsMatrix() {
        const matrix = new THREE.Matrix3();
        matrix.set(
          GALACTIC_TO_ICRS_MATRIX[0], GALACTIC_TO_ICRS_MATRIX[1], GALACTIC_TO_ICRS_MATRIX[2],
          GALACTIC_TO_ICRS_MATRIX[3], GALACTIC_TO_ICRS_MATRIX[4], GALACTIC_TO_ICRS_MATRIX[5],
          GALACTIC_TO_ICRS_MATRIX[6], GALACTIC_TO_ICRS_MATRIX[7], GALACTIC_TO_ICRS_MATRIX[8]
        );
        return matrix;
      }

      function skyDomeAxisTransformMatrix() {
        const transform = volumeSkyAxisTransform || {};
        const xyPermutation = Array.isArray(transform.xyPermutation) && transform.xyPermutation.length === 2
          ? transform.xyPermutation
          : [0, 1];
        const xySigns = Array.isArray(transform.xySigns) && transform.xySigns.length === 2
          ? transform.xySigns
          : [1, 1];
        const zSign = Number.isFinite(Number(transform.zSign)) && Math.abs(Number(transform.zSign)) > 0.5
          ? Number(transform.zSign)
          : 1;
        const rows = [
          [0.0, 0.0, 0.0],
          [0.0, 0.0, 0.0],
          [0.0, 0.0, 0.0],
        ];
        rows[0][Math.min(Math.max(Number(xyPermutation[0]) || 0, 0), 2)] = Number(xySigns[0]) || 1;
        rows[1][Math.min(Math.max(Number(xyPermutation[1]) || 1, 0), 2)] = Number(xySigns[1]) || 1;
        rows[2][2] = zSign;
        const matrix = new THREE.Matrix3();
        matrix.set(
          rows[0][0], rows[0][1], rows[0][2],
          rows[1][0], rows[1][1], rows[1][2],
          rows[2][0], rows[2][1], rows[2][2]
        );
        return matrix;
      }

      function skyDomeLocalSourceName(fallback = "") {
        const candidates = [
          fallback,
          skyDomeSpec && skyDomeSpec.name,
          skyDomeSpec && skyDomeSpec.label,
          skyDomeSpec && skyDomeSpec.source_name,
          skyDomeSpec && skyDomeSpec.survey,
          skyDomeSpec && typeof skyDomeSpec.source === "string" ? skyDomeSpec.source : "",
        ];
        for (const candidate of candidates) {
          const label = String(candidate || "").trim();
          if (label) {
            return label;
          }
        }
        return "local";
      }

      function skyDomeBooleanOption(name, fallback = false) {
        if (!skyDomeSpec || !Object.prototype.hasOwnProperty.call(skyDomeSpec, name)) {
          return Boolean(fallback);
        }
        const value = skyDomeSpec[name];
        if (typeof value === "string") {
          return value.trim().toLowerCase() === "true";
        }
        return Boolean(value);
      }

      function skyDomeLongitudeOffsetRad() {
        const offsetDeg = Number(
          (skyDomeSpec && (
            skyDomeSpec.longitude_offset_deg
            ?? skyDomeSpec.lon_offset_deg
            ?? skyDomeSpec.rotation_deg
          ))
          ?? 0.0
        );
        return Number.isFinite(offsetDeg) ? offsetDeg * Math.PI / 180.0 : 0.0;
      }

      function setSkyDomeSnapshotStatus(status, message) {
        skyDomeSnapshotStatus = String(status || "idle");
        if (root && root.dataset) {
          root.dataset.skyDomeSnapshotStatus = skyDomeSnapshotStatus;
          root.dataset.skyDomeProjection = skyDomeUsesNativeHips()
            ? "HEALPIX"
            : normalizeSkyDomeProjectionName(skyDomeProjection || skyDomeSpec.projection || "CAR");
          root.dataset.skyDomeSurvey = String(skyDomeSurvey || skyDomeLocalSourceName(""));
          if (message) {
            root.dataset.skyDomeSnapshotMessage = String(message);
          } else {
            delete root.dataset.skyDomeSnapshotMessage;
          }
        }
        if (typeof refreshSkyDomeControlStatus === "function") {
          refreshSkyDomeControlStatus();
        }
      }

      function skyDomeCenterForCurrentFrame() {
        const frame = currentFrame();
        if (frame && Array.isArray(frame.traces)) {
          for (const trace of frame.traces) {
            const traceName = String((trace && trace.name) || "").trim().toLowerCase();
            if (traceName !== "sun" && traceName !== "solar system") {
              continue;
            }
            const points = Array.isArray(trace.points) ? trace.points : [];
            for (const point of points) {
              const x = Number(point && point.x);
              const y = Number(point && point.y);
              const z = Number(point && point.z);
              if (Number.isFinite(x) && Number.isFinite(y) && Number.isFinite(z)) {
                return new THREE.Vector3(x, y, z).add(plotGroup.position);
              }
            }
          }
        }
        if (currentZoomAnchorPoint instanceof THREE.Vector3) {
          return currentZoomAnchorPoint.clone().add(plotGroup.position);
        }
        return controls.target.clone();
      }

      function skyDomeWorldCenterForCurrentView() {
        if (skyDomeUsesNativeHips() && camera && camera.position) {
          return camera.position.clone();
        }
        return skyDomeCenterForCurrentFrame();
      }

      function skyDomeOpacityForCurrentView() {
        if (!skyDomeIsEnabled()) {
          return 0.0;
        }
        if (cameraViewMode !== "earth") {
          return 0.0;
        }
        const maxOpacity = Math.min(Math.max(Number(skyDomeSpec.opacity ?? 0.55), 0.0), 1.0);
        return maxOpacity * Math.min(Math.max(Number(skyDomeViewOpacityScale) || 0.0, 0.0), 1.0);
      }

      function skyDomeBackgroundViewForCamera() {
        if (!camera || !canvas) {
          return null;
        }
        const direction = new THREE.Vector3();
        camera.getWorldDirection(direction);
        const galactic = volumeSkyGalacticLonLatDegFromCartesian(direction.x, direction.y, direction.z);
        if (!galactic) {
          return null;
        }
        const icrs = icrsDegFromGalacticDeg(galactic.l, galactic.b);
        if (!icrs) {
          return null;
        }
        const cameraFovDeg = Math.min(Math.max(Number(camera.fov) || 60.0, 0.05), 120.0);
        const width = Math.max(Number(canvas && canvas.clientWidth) || 1.0, 1.0);
        const height = Math.max(Number(canvas && canvas.clientHeight) || 1.0, 1.0);
        const aspect = width / height;
        const horizontalFovDeg = 2.0 * Math.atan(
          Math.tan(THREE.MathUtils.degToRad(cameraFovDeg * 0.5)) * aspect
        ) * 180.0 / Math.PI;
        const fovDeg = Math.min(Math.max(horizontalFovDeg, 0.05), 179.0);
        return {
          l: Number(galactic.l),
          b: Number(galactic.b),
          ra: Number(icrs.ra),
          dec: Number(icrs.dec),
          fovDeg,
          cameraFovDeg,
        };
      }

      function holdSkyDomeBackgroundForCameraMotion(holdMs = 500.0) {
        if (!skyDomeUsesAladinBackground()) {
          return;
        }
        skyDomeBackgroundMotionHoldUntil = 0.0;
      }

      function setSkyDomeBackgroundCameraActive(active) {
        if (!skyDomeUsesAladinBackground()) {
          return;
        }
        skyDomeBackgroundUserCameraActive = Boolean(active);
        holdSkyDomeBackgroundForCameraMotion(0.0);
      }

      function setSkyDomeBackgroundDebugState(reason) {
        if (!root || !root.dataset) {
          return;
        }
        root.dataset.skyDomeBackgroundReady = skyDomeBackgroundFrameReady ? "true" : "false";
        root.dataset.skyDomeBackgroundHasContentWindow = (
          skyDomeFrameEl && skyDomeFrameEl.contentWindow
        ) ? "true" : "false";
        root.dataset.skyDomeBackgroundBlocker = String(reason || "");
      }

      function skyDomeBackgroundCameraDirection() {
        if (!camera) {
          return null;
        }
        const direction = new THREE.Vector3();
        camera.getWorldDirection(direction);
        if (!Number.isFinite(direction.x) || !Number.isFinite(direction.y) || !Number.isFinite(direction.z)) {
          return null;
        }
        return direction.normalize();
      }

      function clearSkyDomeBackgroundPredictiveTransform() {
        if (!skyDomeFrameEl) {
          return;
        }
        skyDomeFrameEl.style.transform = "";
      }

      function updateSkyDomeBackgroundPredictiveTransform(currentView = null) {
        if (!skyDomeFrameEl || !skyDomeUsesAladinBackground()) {
          return;
        }
        if (!skyDomeBackgroundAlignedView) {
          clearSkyDomeBackgroundPredictiveTransform();
          return;
        }
        if (!skyDomeBackgroundUserCameraActive && skyDomeBackgroundSentViews.size === 0) {
          clearSkyDomeBackgroundPredictiveTransform();
          return;
        }
        const alignedDirection = skyDomeBackgroundAlignedView.direction;
        if (!(alignedDirection instanceof THREE.Vector3)) {
          clearSkyDomeBackgroundPredictiveTransform();
          return;
        }
        camera.updateMatrixWorld(true);
        camera.updateProjectionMatrix();
        const distance = Math.max(Number(sceneSpec.max_span) || 1.0, 1000.0);
        const projected = alignedDirection.clone()
          .multiplyScalar(distance)
          .add(camera.position)
          .project(camera);
        if (
          !Number.isFinite(projected.x)
          || !Number.isFinite(projected.y)
          || Math.abs(projected.x) > 3.0
          || Math.abs(projected.y) > 3.0
        ) {
          clearSkyDomeBackgroundPredictiveTransform();
          return;
        }
        const width = Math.max(Number(canvas && canvas.clientWidth) || 1.0, 1.0);
        const height = Math.max(Number(canvas && canvas.clientHeight) || 1.0, 1.0);
        const maxShiftX = width * 0.22;
        const maxShiftY = height * 0.22;
        const dx = Math.min(Math.max(projected.x * width * 0.5, -maxShiftX), maxShiftX);
        const dy = Math.min(Math.max(-projected.y * height * 0.5, -maxShiftY), maxShiftY);
        const alignedFov = Math.min(Math.max(Number(skyDomeBackgroundAlignedView.fovDeg) || 1.0, 0.05), 179.0);
        const currentFov = Math.min(Math.max(Number(currentView && currentView.fovDeg) || alignedFov, 0.05), 179.0);
        const alignedTan = Math.max(Math.tan(THREE.MathUtils.degToRad(alignedFov * 0.5)), 1e-6);
        const currentTan = Math.max(Math.tan(THREE.MathUtils.degToRad(currentFov * 0.5)), 1e-6);
        const scale = Math.min(Math.max(alignedTan / currentTan, 0.55), 1.9);
        if (Math.abs(dx) < 0.35 && Math.abs(dy) < 0.35 && Math.abs(scale - 1.0) < 0.002) {
          clearSkyDomeBackgroundPredictiveTransform();
          return;
        }
        skyDomeFrameEl.style.transform = `translate3d(${dx.toFixed(2)}px, ${dy.toFixed(2)}px, 0) scale(${scale.toFixed(5)})`;
      }

      function recordSkyDomeBackgroundSentView(seq, view) {
        const direction = skyDomeBackgroundCameraDirection();
        if (!direction || !view) {
          return;
        }
        skyDomeBackgroundSentViews.set(seq, {
          direction,
          fovDeg: Number(view.fovDeg) || 0.0,
          cameraFovDeg: Number(view.cameraFovDeg) || Number(view.fovDeg) || 0.0,
        });
        if (!skyDomeBackgroundAlignedView) {
          skyDomeBackgroundAlignedView = skyDomeBackgroundSentViews.get(seq);
        }
        if (skyDomeBackgroundSentViews.size > 8) {
          const keys = Array.from(skyDomeBackgroundSentViews.keys()).sort((a, b) => a - b);
          keys.slice(0, Math.max(0, keys.length - 8)).forEach((key) => skyDomeBackgroundSentViews.delete(key));
        }
      }

      function markSkyDomeBackgroundViewApplied(seq) {
        const safeSeq = Math.round(Number(seq) || 0);
        const appliedView = skyDomeBackgroundSentViews.get(safeSeq);
        if (appliedView) {
          skyDomeBackgroundAlignedView = appliedView;
          Array.from(skyDomeBackgroundSentViews.keys()).forEach((key) => {
            if (key <= safeSeq) {
              skyDomeBackgroundSentViews.delete(key);
            }
          });
        }
        updateSkyDomeBackgroundPredictiveTransform();
      }

      function updateSkyDomeBackgroundFrame(timestampMs = 0.0, options = {}) {
        if (!skyDomeFrameEl || !skyDomeUsesAladinBackground()) {
          setSkyDomeBackgroundDebugState("disabled");
          return;
        }
        const baseOpacity = skyDomeOpacityForCurrentView();
        if (root && root.dataset) {
          root.dataset.skyDomeMode = "aladin-background";
        }
        applySkyDomeFrameVisualSettings();
        if (baseOpacity <= 0.002 || !skyDomeBackgroundFrameReady || !skyDomeFrameEl.contentWindow) {
          skyDomeFrameEl.style.opacity = "0";
          setSkyDomeBackgroundDebugState(
            baseOpacity <= 0.002
              ? "transparent"
              : (!skyDomeBackgroundFrameReady ? "waiting-for-aladin" : "missing-content-window")
          );
          return;
        }
        const view = skyDomeBackgroundViewForCamera();
        if (!view) {
          skyDomeFrameEl.style.opacity = "0";
          setSkyDomeBackgroundDebugState("missing-camera-view");
          return;
        }
        const signature = [
          view.ra.toFixed(3),
          view.dec.toFixed(3),
          view.fovDeg.toFixed(2),
        ].join("|");
        const now = Number(timestampMs) || 0.0;
        if (signature !== skyDomeBackgroundLatestViewSignature) {
          skyDomeBackgroundLatestViewSignature = signature;
        }
        skyDomeFrameEl.style.opacity = String(baseOpacity);
        setSkyDomeBackgroundDebugState("visible");
        updateSkyDomeBackgroundPredictiveTransform(view);
        if (root && root.dataset) {
          root.dataset.skyDomeMotion = skyDomeBackgroundUserCameraActive ? "camera-moving" : "settled";
        }
        const forceUpdate = Boolean(options && options.force);
        if (!forceUpdate && signature === skyDomeBackgroundViewSignature && (now - skyDomeBackgroundLastSentAt) < 500.0) {
          return;
        }
        const minUpdateIntervalMs = skyDomeBackgroundUserCameraActive ? 16.0 : 50.0;
        if (!forceUpdate && (now - skyDomeBackgroundLastSentAt) < minUpdateIntervalMs) {
          return;
        }
        skyDomeBackgroundViewSignature = signature;
        skyDomeBackgroundLastSentAt = now;
        skyDomeBackgroundViewSequence += 1;
        const viewSeq = skyDomeBackgroundViewSequence;
        recordSkyDomeBackgroundSentView(viewSeq, view);
        skyDomeFrameEl.contentWindow.postMessage({
          type: "oviz-sky-background-view",
          seq: viewSeq,
          ra: view.ra,
          dec: view.dec,
          l: view.l,
          b: view.b,
          fovDeg: view.fovDeg,
          cameraFovDeg: view.cameraFovDeg,
        }, "*");
      }

      function skyDomeHipsSurveyName() {
        return String(
          (skyDomeSpec && (skyDomeSpec.hips_survey || skyDomeSpec.survey))
          || (skySpec && skySpec.survey)
          || "P/DSS2/color"
        );
      }

      function skyDomeHipsBaseUrl() {
        const explicitUrl = String(
          (skyDomeSpec && (
            skyDomeSpec.hips_base_url
            || skyDomeSpec.hipsUrl
            || skyDomeSpec.hips_url
            || skyDomeSpec.base_url
          ))
          || ""
        ).trim();
        if (explicitUrl) {
          return explicitUrl.replace(/\/+$/, "");
        }
        const survey = skyDomeHipsSurveyName().toLowerCase();
        if (survey === "p/dss2/color" || survey.includes("dss2/color")) {
          return "https://alasky.cds.unistra.fr/DSS/DSSColor";
        }
        return "";
      }

      function skyDomeHips2FitsServiceUrl() {
        return String(
          (skyDomeSpec && (
            skyDomeSpec.hips2fits_service_url
            || skyDomeSpec.hips_2_fits_service_url
            || skyDomeSpec.hips2fits_url
          ))
          || "https://alasky.cds.unistra.fr/hips-image-services/hips2fits"
        ).trim().replace(/\?+$/, "");
      }

      function skyDomeHips2FitsMaxTextureSize() {
        const rendererLimit = renderer && renderer.capabilities
          ? Number(renderer.capabilities.maxTextureSize)
          : NaN;
        return Math.max(
          1024,
          Math.min(12000, Number.isFinite(rendererLimit) ? Math.floor(rendererLimit) : 12000)
        );
      }

      function skyDomeHips2FitsRequestedWidth() {
        return Math.round(Number(
          skyDomeSpec && (skyDomeSpec.hips2fits_width || skyDomeSpec.hips2fits_width_px)
        ) || 8192);
      }

      function skyDomeHips2FitsWidth() {
        return Math.max(
          1024,
          Math.min(skyDomeHips2FitsMaxTextureSize(), skyDomeHips2FitsRequestedWidth())
        );
      }

      function skyDomeHips2FitsHeight() {
        const requestedWidth = Math.max(1, skyDomeHips2FitsRequestedWidth());
        const safeWidth = skyDomeHips2FitsWidth();
        const requestedHeight = Math.round(Number(
          skyDomeSpec && (skyDomeSpec.hips2fits_height || skyDomeSpec.hips2fits_height_px)
        ) || Math.round(requestedWidth / 2));
        const scaledHeight = Math.round(requestedHeight * Math.min(1.0, safeWidth / requestedWidth));
        return Math.max(
          512,
          Math.min(6000, skyDomeHips2FitsMaxTextureSize(), scaledHeight)
        );
      }

      function skyDomeHips2FitsProjection() {
        return normalizeSkyDomeProjectionName(
          (skyDomeSpec && (
            skyDomeSpec.hips2fits_projection
            || skyDomeSpec.projection
          ))
          || "CAR"
        );
      }

      function skyDomeHips2FitsCoordsys() {
        const coordsys = String(
          (skyDomeSpec && (
            skyDomeSpec.hips2fits_coordsys
            || skyDomeSpec.hips2fits_coordinate_system
            || skyDomeSpec.coordinate_frame
            || skyDomeSpec.hips_frame
          ))
          || "galactic"
        ).trim().toLowerCase();
        return coordsys === "icrs" || coordsys === "equatorial" ? "icrs" : "galactic";
      }

      function skyDomeHips2FitsPreviewWidth() {
        const finalWidth = skyDomeHips2FitsWidth();
        const configuredWidth = Number(
          skyDomeSpec && (
            skyDomeSpec.hips2fits_preview_width
            || skyDomeSpec.hips2fits_preview_width_px
          )
        );
        const previewWidth = Number.isFinite(configuredWidth) && configuredWidth > 0.0
          ? configuredWidth
          : Math.min(1024, finalWidth);
        return Math.max(512, Math.min(finalWidth, Math.round(previewWidth)));
      }

      function skyDomeHips2FitsPreviewHeight() {
        const finalHeight = skyDomeHips2FitsHeight();
        const configuredHeight = Number(
          skyDomeSpec && (
            skyDomeSpec.hips2fits_preview_height
            || skyDomeSpec.hips2fits_preview_height_px
          )
        );
        const previewHeight = Number.isFinite(configuredHeight) && configuredHeight > 0.0
          ? configuredHeight
          : Math.round(skyDomeHips2FitsPreviewWidth() / 2);
        return Math.max(256, Math.min(finalHeight, Math.round(previewHeight)));
      }

      function skyDomeHips2FitsMediumWidth() {
        const finalWidth = skyDomeHips2FitsWidth();
        const previewWidth = skyDomeHips2FitsPreviewWidth();
        const configuredWidth = Number(
          skyDomeSpec && (
            skyDomeSpec.hips2fits_medium_width
            || skyDomeSpec.hips2fits_medium_width_px
          )
        );
        const mediumWidth = Number.isFinite(configuredWidth) && configuredWidth > 0.0
          ? configuredWidth
          : Math.min(4096, finalWidth);
        return Math.max(previewWidth, Math.min(finalWidth, Math.round(mediumWidth)));
      }

      function skyDomeHips2FitsMediumHeight() {
        const finalHeight = skyDomeHips2FitsHeight();
        const previewHeight = skyDomeHips2FitsPreviewHeight();
        const configuredHeight = Number(
          skyDomeSpec && (
            skyDomeSpec.hips2fits_medium_height
            || skyDomeSpec.hips2fits_medium_height_px
          )
        );
        const mediumHeight = Number.isFinite(configuredHeight) && configuredHeight > 0.0
          ? configuredHeight
          : Math.round(skyDomeHips2FitsMediumWidth() / 2);
        return Math.max(previewHeight, Math.min(finalHeight, Math.round(mediumHeight)));
      }

      function skyDomeHips2FitsFormat() {
        const format = String(
          (skyDomeSpec && (
            skyDomeSpec.hips2fits_format
            || skyDomeSpec.hips2fits_image_format
          ))
          || "jpg"
        ).trim().replace(/^\./, "").toLowerCase();
        return format === "png" ? "png" : "jpg";
      }

      function skyDomeHips2FitsCenterFrame() {
        const frame = String(
          (skyDomeSpec && (
            skyDomeSpec.hips2fits_center_frame
            || skyDomeSpec.hips2fits_center_coordsys
          ))
          || (skyDomeHips2FitsCoordsys() === "galactic" ? "galactic" : "icrs")
        ).trim().toLowerCase();
        return frame === "galactic" || frame === "gal" ? "galactic" : "icrs";
      }

      function skyDomeHips2FitsCenterIcrsDeg() {
        if (skyDomeHips2FitsCenterFrame() === "galactic") {
          const lDeg = Number(
            skyDomeSpec && (
              skyDomeSpec.hips2fits_l_deg
              ?? skyDomeSpec.hips2fits_lon_deg
              ?? skyDomeSpec.hips2fits_glon_deg
            )
          );
          const bDeg = Number(
            skyDomeSpec && (
              skyDomeSpec.hips2fits_b_deg
              ?? skyDomeSpec.hips2fits_lat_deg
              ?? skyDomeSpec.hips2fits_glat_deg
            )
          );
          const icrs = icrsDegFromGalacticDeg(
            Number.isFinite(lDeg) ? lDeg : 0.0,
            Number.isFinite(bDeg) ? bDeg : 0.0
          );
          if (icrs) {
            return icrs;
          }
        }
        return {
          ra: Number(skyDomeSpec && skyDomeSpec.hips2fits_ra_deg) || 0.0,
          dec: Number(skyDomeSpec && skyDomeSpec.hips2fits_dec_deg) || 0.0,
        };
      }

      function skyDomeHips2FitsUrl(width = null, height = null) {
        const serviceUrl = skyDomeHips2FitsServiceUrl();
        const params = new URLSearchParams();
        const maxTextureSize = skyDomeHips2FitsMaxTextureSize();
        const safeWidth = Math.max(512, Math.min(maxTextureSize, Math.round(Number(width) || skyDomeHips2FitsWidth())));
        const safeHeight = Math.max(256, Math.min(6000, maxTextureSize, Math.round(Number(height) || skyDomeHips2FitsHeight())));
        const centerIcrs = skyDomeHips2FitsCenterIcrsDeg();
        params.set("hips", skyDomeHipsSurveyName());
        params.set("width", String(safeWidth));
        params.set("height", String(safeHeight));
        params.set("projection", skyDomeHips2FitsProjection());
        params.set("coordsys", skyDomeHips2FitsCoordsys());
        params.set("ra", String(Number(centerIcrs.ra) || 0.0));
        params.set("dec", String(Number(centerIcrs.dec) || 0.0));
        params.set("fov", String(Number(skyDomeSpec && skyDomeSpec.hips2fits_fov_deg) || 360.0));
        params.set("format", skyDomeHips2FitsFormat());
        if (skyDomeSpec && skyDomeSpec.hips2fits_min_cut != null) {
          params.set("min_cut", String(skyDomeSpec.hips2fits_min_cut));
        }
        if (skyDomeSpec && skyDomeSpec.hips2fits_max_cut != null) {
          params.set("max_cut", String(skyDomeSpec.hips2fits_max_cut));
        }
        if (skyDomeSpec && skyDomeSpec.hips2fits_stretch != null) {
          params.set("stretch", String(skyDomeSpec.hips2fits_stretch));
        }
        return `${serviceUrl}?${params.toString()}`;
      }

      function skyDomeHipsTileFormat() {
        const format = String(
          (skyDomeSpec && (skyDomeSpec.hips_tile_format || skyDomeSpec.tile_format))
          || "jpg"
        ).trim().replace(/^\./, "").toLowerCase();
        return format === "png" ? "png" : "jpg";
      }

      function skyDomeHipsFrameName() {
        const frame = String(
          (skyDomeSpec && (skyDomeSpec.hips_frame || skyDomeSpec.frame || skyDomeSpec.coordinate_frame))
          || "icrs"
        ).trim().toLowerCase();
        return frame === "gal" || frame === "galactic" ? "galactic" : "icrs";
      }

      function skyDomeHipsAllskyOrder() {
        return Math.max(0, Math.min(6, Math.round(Number(skyDomeSpec && skyDomeSpec.hips_allsky_order) || 3)));
      }

      function skyDomeHipsTileOrder() {
        const allskyOrder = skyDomeHipsAllskyOrder();
        return Math.max(
          allskyOrder,
          Math.min(9, Math.round(Number(skyDomeSpec && skyDomeSpec.hips_tile_order) || Math.max(allskyOrder + 1, 4)))
        );
      }

      function skyDomeHipsTileOrderForFov(fovDeg) {
        const maxOrder = skyDomeHipsTileOrder();
        const allskyOrder = skyDomeHipsAllskyOrder();
        const fov = Math.max(Number(fovDeg) || 60.0, 0.01);
        let order = allskyOrder;
        if (fov <= 28.0) {
          order = 4;
        }
        if (fov <= 12.0) {
          order = 5;
        }
        if (fov <= 5.0) {
          order = 6;
        }
        if (fov <= 2.2) {
          order = 7;
        }
        if (fov <= 1.0) {
          order = 8;
        }
        if (fov <= 0.45) {
          order = 9;
        }
        return Math.max(allskyOrder, Math.min(maxOrder, order));
      }

      function skyDomeHipsTileSubdivisions() {
        return Math.max(2, Math.min(64, Math.round(Number(skyDomeSpec && skyDomeSpec.hips_tile_subdivisions) || 16)));
      }

      function skyDomeHipsAllskyTileSubdivisions() {
        return Math.max(
          3,
          Math.min(64, Math.round(Number(skyDomeSpec && skyDomeSpec.hips_allsky_tile_subdivisions) || 16))
        );
      }

      function skyDomeHipsMaxActiveTiles() {
        return Math.max(12, Math.min(512, Math.round(Number(skyDomeSpec && skyDomeSpec.hips_max_active_tiles) || 160)));
      }

      function skyDomeHipsMaxConcurrentTileLoads() {
        return Math.max(1, Math.min(32, Math.round(Number(skyDomeSpec && skyDomeSpec.hips_max_concurrent_tile_loads) || 8)));
      }

      function skyDomeHipsStartupPreloadTiles() {
        const configured = Number(skyDomeSpec && skyDomeSpec.hips_startup_preload_tiles);
        const fallback = Math.min(skyDomeHipsMaxActiveTiles(), 96);
        return Math.max(0, Math.min(256, Math.round(Number.isFinite(configured) ? configured : fallback)));
      }

      function skyDomeHipsStartupWaitMs() {
        const configured = Number(skyDomeSpec && skyDomeSpec.hips_startup_wait_ms);
        return Math.max(0.0, Math.min(3000.0, Number.isFinite(configured) ? configured : 900.0));
      }

      function skyDomeHipsBrightness() {
        return Math.max(0.1, Math.min(8.0, Number(skyDomeSpec && skyDomeSpec.hips_brightness) || 2.4));
      }

      function skyDomeHipsContrast() {
        return Math.max(0.1, Math.min(4.0, Number(skyDomeSpec && skyDomeSpec.hips_contrast) || 1.25));
      }

      function skyDomeHipsGamma() {
        return Math.max(0.2, Math.min(4.0, Number(skyDomeSpec && skyDomeSpec.hips_gamma) || 1.35));
      }

      function skyDomeFrameVisualFilterCss() {
        const defaultBrightness = Math.max(Number(skyDomeDefaultHipsBrightness) || 1.0, 0.1);
        const defaultContrast = Math.max(Number(skyDomeDefaultHipsContrast) || 1.0, 0.1);
        const defaultGamma = Math.max(Number(skyDomeDefaultHipsGamma) || 1.0, 0.2);
        const brightness = Math.min(Math.max(skyDomeHipsBrightness() / defaultBrightness, 0.08), 8.0);
        const contrast = Math.min(Math.max(skyDomeHipsContrast() / defaultContrast, 0.1), 4.0);
        const gammaTone = Math.pow(Math.min(Math.max(skyDomeHipsGamma() / defaultGamma, 0.2), 4.0), -0.25);
        return `brightness(${(brightness * gammaTone).toFixed(3)}) contrast(${contrast.toFixed(3)})`;
      }

      function applySkyDomeFrameVisualSettings() {
        if (!skyDomeFrameEl) {
          return;
        }
        skyDomeFrameEl.style.filter = skyDomeFrameVisualFilterCss();
      }

      function skyDomeHipsUpdateIntervalMs() {
        return Math.max(50.0, Math.min(2000.0, Number(skyDomeSpec && skyDomeSpec.hips_update_interval_ms) || 250.0));
      }

      function skyDomeHipsTileUrl(order, pix) {
        const baseUrl = skyDomeHipsBaseUrl();
        const safeOrder = Math.max(0, Math.round(Number(order) || 0));
        const safePix = Math.max(0, Math.round(Number(pix) || 0));
        const dirIndex = Math.floor(safePix / 10000) * 10000;
        return `${baseUrl}/Norder${safeOrder}/Dir${dirIndex}/Npix${safePix}.${skyDomeHipsTileFormat()}`;
      }

      function skyDomeHipsAllskyUrl() {
        return `${skyDomeHipsBaseUrl()}/Norder${skyDomeHipsAllskyOrder()}/Allsky.${skyDomeHipsTileFormat()}`;
      }

      const HEALPIX_JRLL = [2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4];
      const HEALPIX_JPLL = [1, 3, 5, 7, 0, 2, 4, 6, 1, 3, 5, 7];

      function healpixNestedPixelCount(order) {
        const safeOrder = Math.max(0, Math.round(Number(order) || 0));
        return 12 * Math.pow(4, safeOrder);
      }

      function healpixNestedDeinterleaveXY(value, order) {
        const safeOrder = Math.max(0, Math.round(Number(order) || 0));
        const safeValue = Math.max(0, Math.round(Number(value) || 0));
        let ix = 0;
        let iy = 0;
        for (let bit = 0; bit < safeOrder; bit += 1) {
          ix += (((safeValue >> (2 * bit)) & 1) << bit);
          iy += (((safeValue >> ((2 * bit) + 1)) & 1) << bit);
        }
        return { ix, iy };
      }

      function healpixNestedInterleaveXY(ix, iy, order) {
        const safeOrder = Math.max(0, Math.round(Number(order) || 0));
        const safeIx = Math.max(0, Math.round(Number(ix) || 0));
        const safeIy = Math.max(0, Math.round(Number(iy) || 0));
        let value = 0;
        for (let bit = 0; bit < safeOrder; bit += 1) {
          value += (((safeIx >> bit) & 1) << (2 * bit));
          value += (((safeIy >> bit) & 1) << ((2 * bit) + 1));
        }
        return value;
      }

      function healpixNestedXyf(order, pix) {
        const safeOrder = Math.max(0, Math.round(Number(order) || 0));
        const nside = Math.pow(2, safeOrder);
        const npface = nside * nside;
        const safePix = Math.max(0, Math.min(healpixNestedPixelCount(safeOrder) - 1, Math.round(Number(pix) || 0)));
        const face = Math.floor(safePix / npface);
        const xy = healpixNestedDeinterleaveXY(safePix % npface, safeOrder);
        return { face, ix: xy.ix, iy: xy.iy };
      }

      function healpixNestedPixelFromVector(order, vector) {
        if (!vector) {
          return null;
        }
        const safeOrder = Math.max(0, Math.round(Number(order) || 0));
        const nside = Math.pow(2, safeOrder);
        const npface = nside * nside;
        const x = Number(vector.x) || 0.0;
        const y = Number(vector.y) || 0.0;
        const zRaw = Number(vector.z) || 0.0;
        const norm = Math.sqrt((x * x) + (y * y) + (zRaw * zRaw));
        if (!(norm > 1e-12)) {
          return null;
        }
        const z = Math.min(Math.max(zRaw / norm, -1.0), 1.0);
        const za = Math.abs(z);
        let phi = Math.atan2(y, x);
        if (phi < 0.0) {
          phi += Math.PI * 2.0;
        }
        let tt = phi / (Math.PI * 0.5);
        while (tt >= 4.0) {
          tt -= 4.0;
        }
        let face;
        let ix;
        let iy;
        if (za <= (2.0 / 3.0)) {
          const temp1 = nside * (0.5 + tt);
          const temp2 = nside * z * 0.75;
          const jp = Math.floor(temp1 - temp2);
          const jm = Math.floor(temp1 + temp2);
          const ifp = Math.floor(jp / nside);
          const ifm = Math.floor(jm / nside);
          if (ifp === ifm) {
            face = (ifp % 4) + 4;
          } else if (ifp < ifm) {
            face = ifp % 4;
          } else {
            face = (ifm % 4) + 8;
          }
          ix = jm % nside;
          iy = nside - (jp % nside) - 1;
        } else {
          let ntt = Math.floor(tt);
          if (ntt >= 4) {
            ntt = 3;
          }
          const tp = tt - ntt;
          const tmp = nside * Math.sqrt(Math.max(0.0, 3.0 * (1.0 - za)));
          const jp = Math.min(nside - 1, Math.floor(tp * tmp));
          const jm = Math.min(nside - 1, Math.floor((1.0 - tp) * tmp));
          if (z >= 0.0) {
            face = ntt;
            ix = nside - jm - 1;
            iy = nside - jp - 1;
          } else {
            face = ntt + 8;
            ix = jp;
            iy = jm;
          }
        }
        const safeFace = Math.max(0, Math.min(11, Math.round(Number(face) || 0)));
        const safeIx = Math.max(0, Math.min(nside - 1, Math.round(Number(ix) || 0)));
        const safeIy = Math.max(0, Math.min(nside - 1, Math.round(Number(iy) || 0)));
        return (safeFace * npface) + healpixNestedInterleaveXY(safeIx, safeIy, safeOrder);
      }

      function healpixNestedVectorFromXyf(order, face, ix, iy) {
        const safeOrder = Math.max(0, Math.round(Number(order) || 0));
        const nside = Math.pow(2, safeOrder);
        const safeFace = Math.max(0, Math.min(11, Math.round(Number(face) || 0)));
        const safeIx = Math.max(0, Math.min(nside - 1, Math.round(Number(ix) || 0)));
        const safeIy = Math.max(0, Math.min(nside - 1, Math.round(Number(iy) || 0)));
        const jr = (HEALPIX_JRLL[safeFace] * nside) - safeIx - safeIy - 1;
        const nl4 = 4 * nside;
        const fact1 = 1.0 / (3.0 * nside * nside);
        const fact2 = 2.0 / (3.0 * nside);
        let nr;
        let z;
        let kshift;
        if (jr < nside) {
          nr = Math.max(jr, 1);
          z = 1.0 - (nr * nr * fact1);
          kshift = 0;
        } else if (jr > 3 * nside) {
          nr = Math.max(nl4 - jr, 1);
          z = -1.0 + (nr * nr * fact1);
          kshift = 0;
        } else {
          nr = nside;
          z = (2 * nside - jr) * fact2;
          kshift = (jr - nside) & 1;
        }
        let jp = Math.floor(((HEALPIX_JPLL[safeFace] * nr) + safeIx - safeIy + 1 + kshift) / 2);
        if (jp > nl4) {
          jp -= nl4;
        }
        if (jp < 1) {
          jp += nl4;
        }
        const phi = (jp - ((kshift + 1) * 0.5)) * (Math.PI / (2.0 * nr));
        const zz = Math.min(Math.max(z, -1.0), 1.0);
        const sinTheta = Math.sqrt(Math.max(0.0, 1.0 - (zz * zz)));
        return {
          x: sinTheta * Math.cos(phi),
          y: sinTheta * Math.sin(phi),
          z: zz,
        };
      }

      function healpixNestedVectorFromFaceXY(order, face, xFace, yFace) {
        const safeOrder = Math.max(0, Math.round(Number(order) || 0));
        const nside = Math.pow(2, safeOrder);
        const safeFace = Math.max(0, Math.min(11, Math.round(Number(face) || 0)));
        const x = Math.min(Math.max(Number(xFace), 0.0), 1.0);
        const y = Math.min(Math.max(Number(yFace), 0.0), 1.0);
        const jr = HEALPIX_JRLL[safeFace] - x - y;
        let nr = 1.0;
        let z = 0.0;
        let sinTheta = 1.0;
        if (jr < 1.0) {
          nr = Math.max(jr, 0.0);
          const tmp = (nr * nr) / 3.0;
          z = 1.0 - tmp;
          sinTheta = Math.sqrt(Math.max(0.0, tmp * (2.0 - tmp)));
        } else if (jr > 3.0) {
          nr = Math.max(4.0 - jr, 0.0);
          const tmp = (nr * nr) / 3.0;
          z = tmp - 1.0;
          sinTheta = Math.sqrt(Math.max(0.0, tmp * (2.0 - tmp)));
        } else {
          nr = 1.0;
          z = (2.0 - jr) * (2.0 / 3.0);
          sinTheta = Math.sqrt(Math.max(0.0, 1.0 - (z * z)));
        }
        let tmpPhi = (HEALPIX_JPLL[safeFace] * nr) + x - y;
        while (tmpPhi < 0.0) {
          tmpPhi += 8.0;
        }
        while (tmpPhi >= 8.0) {
          tmpPhi -= 8.0;
        }
        const phi = nr < 1e-15 ? 0.0 : (Math.PI * 0.25 * tmpPhi) / nr;
        return {
          x: sinTheta * Math.cos(phi),
          y: sinTheta * Math.sin(phi),
          z: Math.min(Math.max(z, -1.0), 1.0),
        };
      }

      function healpixNestedPixelCenterVector(order, pix) {
        const xyf = healpixNestedXyf(order, pix);
        const nside = Math.pow(2, Math.max(0, Math.round(Number(order) || 0)));
        return healpixNestedVectorFromFaceXY(
          order,
          xyf.face,
          (xyf.ix + 0.5) / nside,
          (xyf.iy + 0.5) / nside
        );
      }

      function lonLatDegFromUnitVector(vector) {
        if (!vector) {
          return null;
        }
        const x = Number(vector.x);
        const y = Number(vector.y);
        const z = Number(vector.z);
        const norm = Math.sqrt((x * x) + (y * y) + (z * z));
        if (!(norm > 1e-9)) {
          return null;
        }
        return {
          lon: normalizeSkyLongitude(Math.atan2(y, x) * 180.0 / Math.PI),
          lat: Math.asin(Math.min(1.0, Math.max(-1.0, z / norm))) * 180.0 / Math.PI,
        };
      }

      function unitVectorFromLonLatDeg(lonDeg, latDeg) {
        const lon = Number(lonDeg) * Math.PI / 180.0;
        const lat = Number(latDeg) * Math.PI / 180.0;
        if (!Number.isFinite(lon) || !Number.isFinite(lat)) {
          return null;
        }
        const cosLat = Math.cos(lat);
        return {
          x: cosLat * Math.cos(lon),
          y: cosLat * Math.sin(lon),
          z: Math.sin(lat),
        };
      }

      function hipsFrameVectorToLocalVector(vector) {
        let galacticVector = vector;
        if (skyDomeHipsFrameName() !== "galactic") {
          const lonLat = lonLatDegFromUnitVector(vector);
          const galactic = lonLat ? galacticDegFromIcrsDeg(lonLat.lon, lonLat.lat) : null;
          galacticVector = galactic ? unitVectorFromLonLatDeg(galactic.l, galactic.b) : null;
        }
        if (!galacticVector) {
          return new THREE.Vector3(1, 0, 0);
        }
        const axisMatrix = skyDomeAxisTransformMatrix().clone().transpose();
        return new THREE.Vector3(
          Number(galacticVector.x) || 0.0,
          Number(galacticVector.y) || 0.0,
          Number(galacticVector.z) || 0.0
        ).applyMatrix3(axisMatrix).normalize();
      }

      function cameraHipsFrameDirection() {
        if (!camera) {
          return null;
        }
        const direction = new THREE.Vector3();
        camera.getWorldDirection(direction);
        const galactic = volumeSkyGalacticLonLatDegFromCartesian(direction.x, direction.y, direction.z);
        if (!galactic) {
          return null;
        }
        if (skyDomeHipsFrameName() === "galactic") {
          return unitVectorFromLonLatDeg(galactic.l, galactic.b);
        }
        const icrs = icrsDegFromGalacticDeg(galactic.l, galactic.b);
        return icrs ? unitVectorFromLonLatDeg(icrs.ra, icrs.dec) : null;
      }

      function screenHipsFrameDirection(ndcX, ndcY) {
        if (!camera) {
          return null;
        }
        const direction = new THREE.Vector3(Number(ndcX) || 0.0, Number(ndcY) || 0.0, 0.5)
          .unproject(camera)
          .sub(camera.position)
          .normalize();
        const galactic = volumeSkyGalacticLonLatDegFromCartesian(direction.x, direction.y, direction.z);
        if (!galactic) {
          return null;
        }
        if (skyDomeHipsFrameName() === "galactic") {
          return unitVectorFromLonLatDeg(galactic.l, galactic.b);
        }
        const icrs = icrsDegFromGalacticDeg(galactic.l, galactic.b);
        return icrs ? unitVectorFromLonLatDeg(icrs.ra, icrs.dec) : null;
      }

      function skyDomeHipsCameraFovDeg() {
        if (!camera || !canvas) {
          return 60.0;
        }
        const verticalFovRad = (Number(camera.fov) || 60.0) * Math.PI / 180.0;
        const aspect = Math.max(Number(canvas.clientWidth) || 1.0, 1.0) / Math.max(Number(canvas.clientHeight) || 1.0, 1.0);
        const horizontalFovDeg = 2.0 * Math.atan(Math.tan(verticalFovRad * 0.5) * aspect) * 180.0 / Math.PI;
        return Math.min(Math.max(Math.max(Number(camera.fov) || 60.0, horizontalFovDeg) * 1.04, 0.05), 140.0);
      }

      function nativeHipsTileRadiusDeg(order) {
        const count = Math.max(1, healpixNestedPixelCount(order));
        return Math.sqrt(41252.96124941927 / count) * 0.85;
      }

      function nativeHipsAtlasInfo(order) {
        const npix = healpixNestedPixelCount(order);
        const cols = Math.max(1, Math.floor(Math.sqrt(npix)));
        return {
          order,
          npix,
          cols,
          rows: Math.max(1, Math.ceil(npix / cols)),
        };
      }

      function nativeHipsLocalTextureUv(u, v) {
        const uu = skyDomeBooleanOption("hips_flip_x") ? 1.0 - u : u;
        const vv = skyDomeBooleanOption("hips_flip_y") ? 1.0 - v : v;
        return { u: uu, v: vv };
      }

      function createNativeHipsTileGeometry(order, pix, radiusPc, divisions, atlasInfo = null) {
        const geometry = new THREE.BufferGeometry();
        const safeDivisions = Math.max(2, Math.min(64, Math.round(Number(divisions) || 16)));
        const vertexCount = (safeDivisions + 1) * (safeDivisions + 1);
        const positions = new Float32Array(vertexCount * 3);
        const uvs = new Float32Array(vertexCount * 2);
        const indices = [];
        const xyf = healpixNestedXyf(order, pix);
        const nside = Math.pow(2, Math.max(0, Math.round(Number(order) || 0)));
        const atlas = atlasInfo || null;
        const atlasCol = atlas ? (pix % atlas.cols) : 0;
        const atlasRow = atlas ? Math.floor(pix / atlas.cols) : 0;
        let vertexIndex = 0;
        for (let y = 0; y <= safeDivisions; y += 1) {
          const v = y / safeDivisions;
          for (let x = 0; x <= safeDivisions; x += 1) {
            const u = x / safeDivisions;
            const vector = healpixNestedVectorFromFaceXY(
              order,
              xyf.face,
              (xyf.ix + u) / nside,
              (xyf.iy + v) / nside
            );
            const localVector = hipsFrameVectorToLocalVector(vector);
            positions[(vertexIndex * 3)] = localVector.x * radiusPc;
            positions[(vertexIndex * 3) + 1] = localVector.y * radiusPc;
            positions[(vertexIndex * 3) + 2] = localVector.z * radiusPc;
            const localUv = nativeHipsLocalTextureUv(u, v);
            if (atlas) {
              uvs[(vertexIndex * 2)] = (atlasCol + localUv.u) / atlas.cols;
              uvs[(vertexIndex * 2) + 1] = 1.0 - ((atlasRow + (1.0 - localUv.v)) / atlas.rows);
            } else {
              uvs[(vertexIndex * 2)] = localUv.u;
              uvs[(vertexIndex * 2) + 1] = localUv.v;
            }
            vertexIndex += 1;
          }
        }
        const rowStride = safeDivisions + 1;
        for (let y = 0; y < safeDivisions; y += 1) {
          for (let x = 0; x < safeDivisions; x += 1) {
            const a = (y * rowStride) + x;
            const b = a + 1;
            const c = a + rowStride;
            const d = c + 1;
            indices.push(a, c, b, b, c, d);
          }
        }
        geometry.setAttribute("position", new THREE.BufferAttribute(positions, 3));
        geometry.setAttribute("uv", new THREE.BufferAttribute(uvs, 2));
        geometry.setIndex(indices);
        geometry.computeBoundingSphere();
        return geometry;
      }

      function createNativeHipsAllskyGeometry(order, radiusPc) {
        const atlas = nativeHipsAtlasInfo(order);
        const divisions = skyDomeHipsAllskyTileSubdivisions();
        const positions = [];
        const uvs = [];
        const indices = [];
        let vertexOffset = 0;
        for (let pix = 0; pix < atlas.npix; pix += 1) {
          const tileGeometry = createNativeHipsTileGeometry(order, pix, radiusPc, divisions, atlas);
          const positionArray = tileGeometry.getAttribute("position").array;
          const uvArray = tileGeometry.getAttribute("uv").array;
          const indexArray = tileGeometry.index.array;
          for (let i = 0; i < positionArray.length; i += 1) {
            positions.push(positionArray[i]);
          }
          for (let i = 0; i < uvArray.length; i += 1) {
            uvs.push(uvArray[i]);
          }
          for (let i = 0; i < indexArray.length; i += 1) {
            indices.push(indexArray[i] + vertexOffset);
          }
          vertexOffset += positionArray.length / 3;
          tileGeometry.dispose();
        }
        const geometry = new THREE.BufferGeometry();
        geometry.setAttribute("position", new THREE.Float32BufferAttribute(positions, 3));
        geometry.setAttribute("uv", new THREE.Float32BufferAttribute(uvs, 2));
        geometry.setIndex(indices);
        geometry.computeBoundingSphere();
        return geometry;
      }

      function configureNativeHipsTexture(texture) {
        if (!texture) {
          return;
        }
        texture.minFilter = THREE.LinearMipmapLinearFilter || THREE.LinearFilter;
        texture.magFilter = THREE.LinearFilter;
        texture.generateMipmaps = true;
        texture.wrapS = THREE.ClampToEdgeWrapping;
        texture.wrapT = THREE.ClampToEdgeWrapping;
        if (renderer && renderer.capabilities && typeof renderer.capabilities.getMaxAnisotropy === "function") {
          texture.anisotropy = Math.max(1, renderer.capabilities.getMaxAnisotropy());
        }
        if ("colorSpace" in texture && THREE.SRGBColorSpace) {
          texture.colorSpace = THREE.SRGBColorSpace;
        } else if ("encoding" in texture && THREE.sRGBEncoding) {
          texture.encoding = THREE.sRGBEncoding;
        }
      }

      function createNativeHipsMaterial(texture, opacity = 0.0) {
        const material = new THREE.ShaderMaterial({
          uniforms: {
            skyTexture: { value: texture || null },
            opacity: { value: opacity },
            brightness: { value: skyDomeHipsBrightness() },
            contrast: { value: skyDomeHipsContrast() },
            gamma: { value: skyDomeHipsGamma() },
          },
          vertexShader: `
            varying vec2 vUv;
            void main() {
              vUv = uv;
              gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0);
            }
          `,
          fragmentShader: `
            uniform sampler2D skyTexture;
            uniform float opacity;
            uniform float brightness;
            uniform float contrast;
            uniform float gamma;
            varying vec2 vUv;
            void main() {
              vec4 sampleColor = texture2D(skyTexture, vUv);
              vec3 color = sampleColor.rgb;
              color = pow(max(color, vec3(0.0)), vec3(1.0 / max(gamma, 0.001)));
              color = ((color - vec3(0.5)) * contrast) + vec3(0.5);
              color *= brightness;
              color = clamp(color, 0.0, 1.0);
              gl_FragColor = vec4(color, sampleColor.a * opacity);
            }
          `,
          transparent: true,
          side: THREE.DoubleSide,
          depthWrite: false,
          depthTest: false,
          toneMapped: false,
        });
        material.name = "oviz-native-hips-material";
        return material;
      }

      function nativeHipsTileCenters(order) {
        if (!skyDomeHipsState) {
          return [];
        }
        if (!skyDomeHipsState.centerCacheByOrder.has(order)) {
          const centers = [];
          const count = healpixNestedPixelCount(order);
          for (let pix = 0; pix < count; pix += 1) {
            centers.push({
              pix,
              vector: healpixNestedPixelCenterVector(order, pix),
            });
          }
          skyDomeHipsState.centerCacheByOrder.set(order, centers);
        }
        return skyDomeHipsState.centerCacheByOrder.get(order) || [];
      }

      function refreshNativeHipsTileDataset() {
        if (!root || !root.dataset || !skyDomeHipsState) {
          return;
        }
        const entries = Array.from(skyDomeHipsState.tileCache.values());
        root.dataset.skyDomeHipsLoadedTiles = String(entries.filter((entry) => entry.loaded).length);
        root.dataset.skyDomeHipsLoadingTiles = String(entries.filter((entry) => entry.loading).length);
        root.dataset.skyDomeHipsPendingTiles = String(skyDomeHipsState.tileLoadQueue.length);
        root.dataset.skyDomeHipsCachedTiles = String(entries.length);
        if (typeof refreshSkyDomeControlStatus === "function") {
          refreshSkyDomeControlStatus();
        }
      }

      function enqueueNativeHipsTile(entry) {
        if (
          !skyDomeHipsState
          || !entry
          || entry.loaded
          || entry.loading
          || entry.failed
          || skyDomeHipsState.tileLoadQueueKeys.has(entry.key)
        ) {
          return;
        }
        skyDomeHipsState.tileLoadQueue.push(entry);
        skyDomeHipsState.tileLoadQueueKeys.add(entry.key);
      }

      function applyNativeHipsMaterialSettings(material, texture = null, opacity = null) {
        if (!material || !material.uniforms) {
          return;
        }
        if (texture && material.uniforms.skyTexture) {
          material.uniforms.skyTexture.value = texture;
        }
        if (opacity !== null && material.uniforms.opacity) {
          material.uniforms.opacity.value = opacity;
        }
        if (material.uniforms.brightness) {
          material.uniforms.brightness.value = skyDomeHipsBrightness();
        }
        if (material.uniforms.contrast) {
          material.uniforms.contrast.value = skyDomeHipsContrast();
        }
        if (material.uniforms.gamma) {
          material.uniforms.gamma.value = skyDomeHipsGamma();
        }
      }

      function pumpNativeHipsTileLoads() {
        if (!skyDomeHipsState) {
          return;
        }
        const maxConcurrent = skyDomeHipsMaxConcurrentTileLoads();
        while (
          skyDomeHipsState.activeTileLoads < maxConcurrent
          && skyDomeHipsState.tileLoadQueue.length
        ) {
          const entry = skyDomeHipsState.tileLoadQueue.shift();
          skyDomeHipsState.tileLoadQueueKeys.delete(entry.key);
          if (!entry || entry.loaded || entry.loading || entry.failed || !skyDomeHipsState.tileCache.has(entry.key)) {
            continue;
          }
          entry.loading = true;
          skyDomeHipsState.activeTileLoads += 1;
          const texture = skyDomeHipsState.loader.load(
            skyDomeHipsTileUrl(entry.order, entry.pix),
            () => {
              texture.needsUpdate = true;
              entry.loaded = true;
              entry.loading = false;
              entry.failed = false;
              skyDomeHipsState.activeTileLoads = Math.max(0, skyDomeHipsState.activeTileLoads - 1);
              if (entry.mesh) {
                entry.mesh.visible = skyDomeHipsState.selectedTileKeys.has(entry.key);
                applyNativeHipsMaterialSettings(entry.mesh.material, texture, skyDomeOpacityForCurrentView());
                entry.mesh.material.needsUpdate = true;
              }
              setSkyDomeSnapshotStatus("loaded", "");
              refreshNativeHipsTileDataset();
              pumpNativeHipsTileLoads();
            },
            undefined,
            (err) => {
              entry.loading = false;
              entry.failed = true;
              skyDomeHipsState.activeTileLoads = Math.max(0, skyDomeHipsState.activeTileLoads - 1);
              if (entry.mesh) {
                entry.mesh.visible = false;
              }
              if (!skyDomeHipsState.reportedTileError) {
                skyDomeHipsState.reportedTileError = true;
                setSkyDomeSnapshotStatus(
                  "hips-tile-error",
                  err && err.message ? err.message : String(err || "Unable to load one or more HiPS tiles.")
                );
              }
              refreshNativeHipsTileDataset();
              pumpNativeHipsTileLoads();
            }
          );
          configureNativeHipsTexture(texture);
          entry.texture = texture;
          if (entry.mesh) {
            applyNativeHipsMaterialSettings(entry.mesh.material, texture, null);
            entry.mesh.material.needsUpdate = true;
          }
        }
        refreshNativeHipsTileDataset();
      }

      function ensureNativeHipsTile(order, pix, opacity) {
        if (!skyDomeHipsState) {
          return null;
        }
        const key = `${order}/${pix}`;
        const now = (typeof performance !== "undefined" && performance.now) ? performance.now() : Date.now();
        let entry = skyDomeHipsState.tileCache.get(key);
        if (entry) {
          entry.lastUsedAt = now;
          if (entry.mesh && entry.loaded) {
            entry.mesh.visible = true;
            applyNativeHipsMaterialSettings(entry.mesh.material, entry.texture, opacity);
          } else {
            enqueueNativeHipsTile(entry);
            pumpNativeHipsTileLoads();
          }
          return entry;
        }
        const geometry = createNativeHipsTileGeometry(
          order,
          pix,
          skyDomeHipsState.radiusPc,
          skyDomeHipsTileSubdivisions()
        );
        const material = createNativeHipsMaterial(null, 0.0);
        const mesh = new THREE.Mesh(geometry, material);
        mesh.name = `oviz-native-hips-tile-${order}-${pix}`;
        mesh.renderOrder = -998;
        mesh.frustumCulled = false;
        mesh.visible = false;
        skyDomeHipsState.detailGroup.add(mesh);
        entry = {
          key,
          order,
          pix,
          mesh,
          texture: null,
          loaded: false,
          loading: false,
          failed: false,
          lastUsedAt: now,
        };
        skyDomeHipsState.tileCache.set(key, entry);
        enqueueNativeHipsTile(entry);
        pumpNativeHipsTileLoads();
        return entry;
      }

      function nativeHipsTilesForCurrentView(order, maxActiveTiles) {
        const centerDirection = cameraHipsFrameDirection();
        if (!centerDirection) {
          return [];
        }
        const fovDeg = skyDomeHipsCameraFovDeg();
        const samplesX = Math.max(9, Math.min(33, Math.ceil(Math.sqrt(maxActiveTiles) * 1.55)));
        const samplesY = Math.max(7, Math.min(25, Math.ceil(samplesX * 0.62)));
        const margin = fovDeg > 50.0 ? 1.12 : 1.06;
        const selectedByPix = new Map();
        for (let yIndex = 0; yIndex < samplesY; yIndex += 1) {
          const yT = samplesY <= 1 ? 0.5 : yIndex / (samplesY - 1);
          const ndcY = (1.0 - (2.0 * yT)) * margin;
          for (let xIndex = 0; xIndex < samplesX; xIndex += 1) {
            const xT = samplesX <= 1 ? 0.5 : xIndex / (samplesX - 1);
            const ndcX = ((2.0 * xT) - 1.0) * margin;
            const direction = screenHipsFrameDirection(ndcX, ndcY);
            const pix = healpixNestedPixelFromVector(order, direction);
            if (pix === null) {
              continue;
            }
            selectedByPix.set(pix, true);
          }
        }
        const items = Array.from(selectedByPix.keys()).map((pix) => {
          const center = healpixNestedPixelCenterVector(order, pix);
          const dot = (
            (Number(centerDirection.x) || 0.0) * (Number(center.x) || 0.0)
            + (Number(centerDirection.y) || 0.0) * (Number(center.y) || 0.0)
            + (Number(centerDirection.z) || 0.0) * (Number(center.z) || 0.0)
          );
          return { pix, dot };
        });
        items.sort((a, b) => b.dot - a.dot);
        return items.slice(0, maxActiveTiles);
      }

      function preloadNativeHipsStartupTiles() {
        if (
          !skyDomeIsEnabled()
          || typeof skyDomeUsesNativeHips !== "function"
          || !skyDomeUsesNativeHips()
          || !camera
          || !canvas
        ) {
          return;
        }
        const preloadCount = skyDomeHipsStartupPreloadTiles();
        const urls = [];
        const allskyUrl = skyDomeHipsAllskyUrl();
        if (allskyUrl) {
          urls.push(allskyUrl);
        }
        const fovDeg = skyDomeHipsCameraFovDeg();
        const tileOrder = skyDomeHipsTileOrderForFov(fovDeg);
        if (tileOrder > skyDomeHipsAllskyOrder()) {
          nativeHipsTilesForCurrentView(tileOrder, preloadCount).forEach((item) => {
            urls.push(skyDomeHipsTileUrl(tileOrder, item.pix));
          });
        }
        const seenUrls = new Set();
        skyDomeHipsWarmupImages = urls.filter((url) => {
          if (!url || seenUrls.has(url)) {
            return false;
          }
          seenUrls.add(url);
          return true;
        }).map((url) => {
          const image = new Image();
          image.crossOrigin = "anonymous";
          image.decoding = "async";
          image.src = url;
          return image;
        });
        if (root && root.dataset) {
          root.dataset.skyDomeHipsPreloadOrder = String(tileOrder);
          root.dataset.skyDomeHipsPreloadedTiles = String(Math.max(0, skyDomeHipsWarmupImages.length - 1));
        }
      }

      async function waitForNativeHipsStartup() {
        if (!skyDomeHipsState || !skyDomeUsesNativeHips()) {
          return;
        }
        const maxWaitMs = skyDomeHipsStartupWaitMs();
        if (!(maxWaitMs > 0.0)) {
          return;
        }
        const startedAt = (typeof performance !== "undefined" && performance.now) ? performance.now() : Date.now();
        while (true) {
          const entries = Array.from(skyDomeHipsState.tileCache.values());
          const selectedKeys = skyDomeHipsState.selectedTileKeys || new Set();
          const activeCount = selectedKeys.size;
          const loadedActiveCount = entries.filter((entry) => entry.loaded && selectedKeys.has(entry.key)).length;
          const allskyLoaded = Boolean(skyDomeHipsState.allskyTexture && skyDomeHipsState.allskyTexture.image);
          if (
            (activeCount <= 0 && allskyLoaded)
            ||
            (activeCount > 0 && loadedActiveCount >= activeCount)
            || (loadedActiveCount > 0 && allskyLoaded)
          ) {
            return;
          }
          const now = (typeof performance !== "undefined" && performance.now) ? performance.now() : Date.now();
          if ((now - startedAt) >= maxWaitMs) {
            return;
          }
          await new Promise((resolve) => setTimeout(resolve, 40));
        }
      }

      function pruneNativeHipsTileCache() {
        if (!skyDomeHipsState) {
          return;
        }
        const maxCacheEntries = Math.max(skyDomeHipsMaxActiveTiles() * 3, skyDomeHipsMaxActiveTiles() + 24);
        if (skyDomeHipsState.tileCache.size <= maxCacheEntries) {
          return;
        }
        const entries = Array.from(skyDomeHipsState.tileCache.values())
          .filter((entry) => !entry.loading && !skyDomeHipsState.selectedTileKeys.has(entry.key))
          .sort((a, b) => (Number(a.lastUsedAt) || 0) - (Number(b.lastUsedAt) || 0));
        while (skyDomeHipsState.tileCache.size > maxCacheEntries && entries.length) {
          const entry = entries.shift();
          if (!entry) {
            break;
          }
          if (entry.mesh) {
            skyDomeHipsState.detailGroup.remove(entry.mesh);
            if (entry.mesh.geometry && typeof entry.mesh.geometry.dispose === "function") {
              entry.mesh.geometry.dispose();
            }
            if (entry.mesh.material) {
              if (entry.texture && typeof entry.texture.dispose === "function") {
                entry.texture.dispose();
              }
              entry.mesh.material.dispose();
            }
          }
          skyDomeHipsState.tileCache.delete(entry.key);
        }
      }

      function updateNativeHipsTiles(timestampMs = 0.0, force = false) {
        if (!skyDomeHipsState || !skyDomeUsesNativeHips()) {
          return;
        }
        const now = Number(timestampMs) || ((typeof performance !== "undefined" && performance.now) ? performance.now() : Date.now());
        if (!force && (now - skyDomeHipsState.lastTileUpdateAt) < skyDomeHipsUpdateIntervalMs()) {
          return;
        }
        skyDomeHipsState.lastTileUpdateAt = now;
        const direction = cameraHipsFrameDirection();
        if (!direction) {
          return;
        }
        const fovDeg = skyDomeHipsCameraFovDeg();
        const tileOrder = skyDomeHipsTileOrderForFov(fovDeg);
        const maxActiveTiles = skyDomeHipsMaxActiveTiles();
        const active = tileOrder > skyDomeHipsAllskyOrder()
          ? nativeHipsTilesForCurrentView(tileOrder, maxActiveTiles)
          : [];
        const activeKeys = new Set(active.map((item) => `${tileOrder}/${item.pix}`));
        skyDomeHipsState.selectedTileKeys = activeKeys;
        skyDomeHipsState.tileCache.forEach((entry) => {
          if (entry.mesh) {
            entry.mesh.visible = activeKeys.has(entry.key) && entry.loaded;
          }
        });
        const opacity = skyDomeOpacityForCurrentView();
        active.forEach((item) => {
          ensureNativeHipsTile(tileOrder, item.pix, opacity);
        });
        if (root && root.dataset) {
          root.dataset.skyDomeMode = "native-hips";
          root.dataset.skyDomeHipsOrder = String(tileOrder);
          root.dataset.skyDomeHipsActiveTiles = String(active.length);
        }
        pumpNativeHipsTileLoads();
        refreshNativeHipsTileDataset();
        pruneNativeHipsTileCache();
      }

      function initializeNativeHipsSkyDome() {
        if (!skyDomeIsEnabled()) {
          setSkyDomeSnapshotStatus("disabled", "");
          return false;
        }
        const baseUrl = skyDomeHipsBaseUrl();
        if (!baseUrl) {
          setSkyDomeSnapshotStatus(
            "missing-hips-url",
            "Native HiPS sky mode needs sky_dome_hips_base_url for this survey."
          );
          return false;
        }
        if (skyDomeMesh && skyDomeMesh.userData && skyDomeMesh.userData.kind === "native-hips") {
          updateSkyDome();
          return true;
        }
        disposeSkyDomeObject();
        disposeSkyDomeTextures(skyDomeTexture);
        skyDomeTexture = null;

        const radiusPc = Math.max(Number(skyDomeSpec.radius_pc) || 40000.0, 100.0);
        const group = new THREE.Group();
        group.name = "oviz-native-hips-sky-dome";
        group.userData.kind = "native-hips";
        group.renderOrder = -1000;
        group.frustumCulled = false;
        group.visible = false;
        const detailGroup = new THREE.Group();
        detailGroup.name = "oviz-native-hips-detail-tiles";
        group.add(detailGroup);
        skyDomeMesh = group;
        skyDomeHipsState = {
          group,
          detailGroup,
          radiusPc,
          loader: new THREE.TextureLoader(),
          tileCache: new Map(),
          selectedTileKeys: new Set(),
          centerCacheByOrder: new Map(),
          tileLoadQueue: [],
          tileLoadQueueKeys: new Set(),
          activeTileLoads: 0,
          lastTileUpdateAt: 0.0,
          reportedTileError: false,
        };
        if (typeof skyDomeHipsState.loader.setCrossOrigin === "function") {
          skyDomeHipsState.loader.setCrossOrigin("anonymous");
        }

        skyDomeSurvey = skyDomeHipsSurveyName();
        skyDomeProjection = "HEALPIX";
        if (root && root.dataset) {
          root.dataset.skyDomeMode = "native-hips";
          root.dataset.skyDomeSurvey = skyDomeSurvey;
          root.dataset.skyDomeProjection = "HEALPIX";
          root.dataset.skyDomeHipsBaseUrl = baseUrl;
        }
        setSkyDomeSnapshotStatus("loading", "Loading native HiPS sky tiles.");

        const allskyOrder = skyDomeHipsAllskyOrder();
        const allskyTexture = skyDomeHipsState.loader.load(
          skyDomeHipsAllskyUrl(),
          () => {
            allskyTexture.needsUpdate = true;
            setSkyDomeSnapshotStatus("loaded", "");
            updateSkyDome();
          },
          undefined,
          (err) => {
            setSkyDomeSnapshotStatus(
              "hips-allsky-error",
              err && err.message ? err.message : String(err || "Unable to load the HiPS all-sky preview.")
            );
          }
        );
        configureNativeHipsTexture(allskyTexture);
        skyDomeHipsState.allskyTexture = allskyTexture;
        const allskyGeometry = createNativeHipsAllskyGeometry(allskyOrder, radiusPc);
        const allskyMaterial = createNativeHipsMaterial(allskyTexture, 0.0);
        const allskyMesh = new THREE.Mesh(allskyGeometry, allskyMaterial);
        allskyMesh.name = `oviz-native-hips-allsky-order-${allskyOrder}`;
        allskyMesh.renderOrder = -999;
        allskyMesh.frustumCulled = false;
        allskyMesh.visible = true;
        skyDomeHipsState.allskyMesh = allskyMesh;
        group.add(allskyMesh);
        scene.add(group);
        updateNativeHipsTiles(0.0, true);
        updateSkyDome();
        return true;
      }

      function applySkyDomeMaterialUniforms(material, texture = null) {
        if (!material || !material.uniforms) {
          return;
        }
        if (texture && material.uniforms.skyTexture) {
          material.uniforms.skyTexture.value = texture;
        }
        if (material.uniforms.projectionMode) {
          material.uniforms.projectionMode.value = skyDomeProjectionMode(skyDomeProjection);
        }
        if (material.uniforms.longitudeOffsetRad) {
          material.uniforms.longitudeOffsetRad.value = skyDomeLongitudeOffsetRad();
        }
        if (material.uniforms.galacticAxisMatrix) {
          material.uniforms.galacticAxisMatrix.value.copy(skyDomeAxisTransformMatrix());
        }
        if (material.uniforms.galacticToIcrsMatrix) {
          material.uniforms.galacticToIcrsMatrix.value.copy(skyDomeGalacticToIcrsMatrix());
        }
        if (material.uniforms.textureFrameMode) {
          material.uniforms.textureFrameMode.value = skyDomeTextureFrameMode();
        }
        if (material.uniforms.flipX) {
          material.uniforms.flipX.value = skyDomeBooleanOption("flip_x") ? 1.0 : 0.0;
        }
        if (material.uniforms.flipY) {
          material.uniforms.flipY.value = skyDomeBooleanOption("flip_y") ? 1.0 : 0.0;
        }
      }

      function createSkyDomeMaterial(texture) {
        const material = new THREE.ShaderMaterial({
          uniforms: {
            skyTexture: { value: texture },
            opacity: { value: 0.0 },
            projectionMode: { value: skyDomeProjectionMode(skyDomeProjection) },
            longitudeOffsetRad: { value: skyDomeLongitudeOffsetRad() },
            galacticAxisMatrix: { value: skyDomeAxisTransformMatrix() },
            galacticToIcrsMatrix: { value: skyDomeGalacticToIcrsMatrix() },
            textureFrameMode: { value: skyDomeTextureFrameMode() },
            flipX: { value: skyDomeBooleanOption("flip_x") ? 1.0 : 0.0 },
            flipY: { value: skyDomeBooleanOption("flip_y") ? 1.0 : 0.0 },
          },
          vertexShader: `
            varying vec3 vLocalDirection;
            void main() {
              vLocalDirection = normalize(position);
              gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0);
            }
          `,
          fragmentShader: `
            uniform sampler2D skyTexture;
            uniform float opacity;
            uniform int projectionMode;
            uniform float longitudeOffsetRad;
            uniform mat3 galacticAxisMatrix;
            uniform mat3 galacticToIcrsMatrix;
            uniform int textureFrameMode;
            uniform float flipX;
            uniform float flipY;
            varying vec3 vLocalDirection;
            const float PI = 3.1415926535897932384626433832795;
            const float ROOT_TWO = 1.4142135623730951;

            float wrap01(float value) {
              return value - floor(value);
            }

            float solveMollweideTheta(float lat) {
              if (abs(abs(lat) - (0.5 * PI)) < 0.0001) {
                return sign(lat) * 0.5 * PI;
              }
              float theta = lat;
              for (int i = 0; i < 6; i += 1) {
                float f = (2.0 * theta) + sin(2.0 * theta) - (PI * sin(lat));
                float fp = max(2.0 + (2.0 * cos(2.0 * theta)), 0.0001);
                theta -= f / fp;
              }
              return theta;
            }

            void main() {
              vec3 galacticDir = normalize(galacticAxisMatrix * normalize(vLocalDirection));
              vec3 dir = textureFrameMode == 1
                ? normalize(galacticToIcrsMatrix * galacticDir)
                : galacticDir;
              float lon = atan(dir.y, dir.x) + longitudeOffsetRad;
              float lat = asin(clamp(dir.z, -1.0, 1.0));
              vec2 uv;
              if (projectionMode == 1) {
                uv = vec2(
                  wrap01(0.5 - (lon / (2.0 * PI))),
                  clamp(0.5 - (lat / PI), 0.0, 1.0)
                );
              } else {
                float theta = solveMollweideTheta(lat);
                float x = (2.0 * ROOT_TWO / PI) * lon * cos(theta);
                float y = ROOT_TWO * sin(theta);
                uv = vec2(
                  0.5 - (x / (4.0 * ROOT_TWO)),
                  0.5 - (y / (2.0 * ROOT_TWO))
                );
                if (uv.x < 0.0 || uv.x > 1.0 || uv.y < 0.0 || uv.y > 1.0) {
                  discard;
                }
              }
              if (flipX > 0.5) {
                uv.x = 1.0 - uv.x;
              }
              if (flipY > 0.5) {
                uv.y = 1.0 - uv.y;
              }
              vec4 sampleColor = texture2D(skyTexture, uv);
              gl_FragColor = vec4(sampleColor.rgb, sampleColor.a * opacity);
            }
          `,
          side: THREE.BackSide,
          transparent: true,
          depthWrite: false,
          depthTest: false,
          toneMapped: false,
        });
        applySkyDomeMaterialUniforms(material, texture);
        return material;
      }

      function configureSkyDomeTexture(texture) {
        if (!texture) {
          return;
        }
        texture.minFilter = THREE.LinearFilter;
        texture.magFilter = THREE.LinearFilter;
        texture.generateMipmaps = true;
        texture.wrapS = normalizeSkyDomeProjectionName(skyDomeProjection) === "CAR"
          ? THREE.RepeatWrapping
          : THREE.ClampToEdgeWrapping;
        texture.wrapT = THREE.ClampToEdgeWrapping;
        if (renderer && renderer.capabilities && typeof renderer.capabilities.getMaxAnisotropy === "function") {
          texture.anisotropy = Math.max(1, renderer.capabilities.getMaxAnisotropy());
        }
        if ("colorSpace" in texture && THREE.SRGBColorSpace) {
          texture.colorSpace = THREE.SRGBColorSpace;
        } else if ("encoding" in texture && THREE.sRGBEncoding) {
          texture.encoding = THREE.sRGBEncoding;
        }
      }

      function disposeSkyDomeTextures(value) {
        if (!value) {
          return;
        }
        if (value && typeof value.dispose === "function") {
          value.dispose();
        }
      }

      function disposeSkyDomeObject() {
        if (!skyDomeMesh) {
          return;
        }
        scene.remove(skyDomeMesh);
        skyDomeMesh.traverse((child) => {
          if (child.geometry && typeof child.geometry.dispose === "function") {
            child.geometry.dispose();
          }
          const materials = Array.isArray(child.material) ? child.material : [child.material];
          materials.forEach((material) => {
            if (material && material.map && typeof material.map.dispose === "function") {
              material.map.dispose();
            }
            if (material && typeof material.dispose === "function") {
              material.dispose();
            }
          });
        });
        skyDomeMesh = null;
        skyDomeHipsState = null;
        skyDomeTexturePriority = -1;
      }

      function setSkyDomeObjectOpacity(opacity) {
        if (!skyDomeMesh) {
          return;
        }
        skyDomeMesh.traverse((child) => {
          const materials = Array.isArray(child.material) ? child.material : [child.material];
          materials.forEach((material) => {
            if (!material) {
              return;
            }
            if (material.uniforms && material.uniforms.opacity) {
              material.uniforms.opacity.value = opacity;
            } else if (Object.prototype.hasOwnProperty.call(material, "opacity")) {
              material.opacity = opacity;
              material.transparent = true;
            }
          });
        });
      }

      function installSkyDomeImageTexture(texture, imageUrl, sourceName, projectionName, modeName = "remote-image", options = {}) {
        const priority = Number.isFinite(Number(options.priority)) ? Number(options.priority) : 0;
        if (skyDomeTexture && priority < skyDomeTexturePriority) {
          disposeSkyDomeTextures(texture);
          return false;
        }
        const normalizedProjection = normalizeSkyDomeProjectionName(projectionName || skyDomeSpec.projection || "CAR");
        skyDomeProjection = normalizedProjection;
        skyDomeSurvey = skyDomeLocalSourceName(sourceName || (skySpec && skySpec.survey) || "remote");
        if (root && root.dataset) {
          root.dataset.skyDomeMode = String(modeName || "remote-image");
          root.dataset.skyDomeImageUrl = imageUrl;
          root.dataset.skyDomeTexturePriority = String(priority);
        }

        const previousTexture = skyDomeTexture;
        skyDomeTexture = texture;
        skyDomeTexturePriority = priority;

        if (!skyDomeMesh || !(skyDomeMesh.userData && skyDomeMesh.userData.kind === "image")) {
          disposeSkyDomeObject();
          const radiusPc = Math.max(Number(skyDomeSpec.radius_pc) || 40000.0, 100.0);
          const geometry = new THREE.SphereGeometry(radiusPc, 160, 80);
          const material = createSkyDomeMaterial(texture);
          skyDomeMesh = new THREE.Mesh(geometry, material);
          skyDomeMesh.userData.kind = "image";
          skyDomeMesh.name = "oviz-sky-dome";
          skyDomeMesh.renderOrder = -1000;
          skyDomeMesh.frustumCulled = false;
          skyDomeMesh.visible = false;
          scene.add(skyDomeMesh);
        } else {
          applySkyDomeMaterialUniforms(skyDomeMesh.material, texture);
        }
        if (previousTexture !== texture) {
          disposeSkyDomeTextures(previousTexture);
        }
        updateSkyDome();
        return true;
      }

      function setSkyDomeTextureFromImageUrl(imageUrl, sourceName, projectionName, modeName = "remote-image", options = {}) {
        if (!skyDomeIsEnabled()) {
          return false;
        }
        if (typeof imageUrl !== "string" || !imageUrl.trim()) {
          setSkyDomeSnapshotStatus(
            "missing-image-url",
            "The sky dome needs a valid remote image URL."
          );
          return false;
        }

        const normalizedProjection = normalizeSkyDomeProjectionName(projectionName || skyDomeSpec.projection || "CAR");
        skyDomeProjection = normalizedProjection;
        skyDomeSurvey = skyDomeLocalSourceName(sourceName || (skySpec && skySpec.survey) || "remote");
        if (root && root.dataset) {
          root.dataset.skyDomePendingImageUrl = imageUrl;
        }
        if (
          skyDomeTexture
          && skyDomeTexture.userData
          && skyDomeTexture.userData.sourceImageUrl === imageUrl
          && skyDomeProjection === normalizedProjection
        ) {
          updateSkyDome();
          return true;
        }

        setSkyDomeSnapshotStatus(
          options.status || "loading",
          options.message || `Loading ${skyDomeSurvey} sky texture.`
        );
        const loader = new THREE.TextureLoader();
        if (typeof loader.setCrossOrigin === "function") {
          loader.setCrossOrigin("anonymous");
        }
        const texture = loader.load(
          imageUrl,
          () => {
            texture.needsUpdate = true;
            let installed = !options.deferInstall
              ? skyDomeTexture === texture
              : false;
            if (options.deferInstall) {
              installed = installSkyDomeImageTexture(
                texture,
                imageUrl,
                sourceName,
                projectionName,
                modeName,
                { priority: options.priority }
              );
            }
            if (installed || skyDomeTexture === texture) {
              setSkyDomeSnapshotStatus(options.loadedStatus || "loaded", options.loadedMessage || "");
              updateSkyDome();
            }
            if ((installed || skyDomeTexture === texture) && typeof options.onLoaded === "function") {
              options.onLoaded(texture);
            }
          },
          undefined,
          (err) => {
            const fallbackMessage = options.errorMessage || "Unable to load the sky texture.";
            setSkyDomeSnapshotStatus(
              options.errorStatus || "texture-error",
              err && typeof err.message === "string" && err.message.trim()
                ? err.message
                : fallbackMessage
            );
          }
        );
        texture.userData = Object.assign({}, texture.userData || {}, { sourceImageUrl: imageUrl });
        configureSkyDomeTexture(texture);

        if (!options.deferInstall) {
          installSkyDomeImageTexture(
            texture,
            imageUrl,
            sourceName,
            projectionName,
            modeName,
            { priority: options.priority }
          );
        }
        return true;
      }

      function setSkyDomeTextureFromDataUrl(dataUrl, sourceName, projectionName) {
        if (!skyDomeIsEnabled()) {
          return false;
        }
        if (typeof dataUrl !== "string" || !dataUrl.startsWith("data:image/")) {
          setSkyDomeSnapshotStatus(
            "missing-image",
            "The Aladin sky capture did not return a data:image URL."
          );
          return false;
        }

        const normalizedProjection = normalizeSkyDomeProjectionName(projectionName || skyDomeSpec.projection || "CAR");
        skyDomeProjection = normalizedProjection;
        skyDomeSurvey = skyDomeLocalSourceName(sourceName || (skySpec && skySpec.survey) || "local");
        if (root && root.dataset) {
          root.dataset.skyDomeMode = skyDomeUsesAladinBackground() ? "aladin-background" : "snapshot";
        }
        if (
          skyDomeTexture
          && skyDomeTexture.userData
          && skyDomeTexture.userData.sourceDataUrl === dataUrl
          && skyDomeProjection === normalizedProjection
        ) {
          updateSkyDome();
          return true;
        }

        setSkyDomeSnapshotStatus("loading", "");
        const texture = new THREE.TextureLoader().load(
          dataUrl,
          () => {
            texture.needsUpdate = true;
            setSkyDomeSnapshotStatus("loaded", "");
            updateSkyDome();
          },
          undefined,
          (err) => {
            setSkyDomeSnapshotStatus(
              "texture-error",
              err && err.message ? err.message : String(err || "")
            );
          }
        );
        texture.userData = Object.assign({}, texture.userData || {}, { sourceDataUrl: dataUrl });
        configureSkyDomeTexture(texture);

        const previousTexture = skyDomeTexture;
        skyDomeTexture = texture;

        if (!skyDomeMesh || !(skyDomeMesh.userData && skyDomeMesh.userData.kind === "image")) {
          disposeSkyDomeObject();
          const radiusPc = Math.max(Number(skyDomeSpec.radius_pc) || 40000.0, 100.0);
          const geometry = new THREE.SphereGeometry(radiusPc, 160, 80);
          const material = createSkyDomeMaterial(texture);
          skyDomeMesh = new THREE.Mesh(geometry, material);
          skyDomeMesh.userData.kind = "image";
          skyDomeMesh.name = "oviz-sky-dome";
          skyDomeMesh.renderOrder = -1000;
          skyDomeMesh.frustumCulled = false;
          skyDomeMesh.visible = false;
          scene.add(skyDomeMesh);
        } else {
          applySkyDomeMaterialUniforms(skyDomeMesh.material, texture);
        }
        disposeSkyDomeTextures(previousTexture);
        updateSkyDome();
        return true;
      }

      function initializeHips2FitsSkyDome() {
        if (!skyDomeIsEnabled()) {
          setSkyDomeSnapshotStatus("disabled", "");
          return false;
        }
        const projectionName = skyDomeHips2FitsProjection();
        skyDomeProjection = projectionName;
        if (skyDomeSpec) {
          skyDomeSpec.projection = projectionName;
        }
        const fullWidth = skyDomeHips2FitsWidth();
        const fullHeight = skyDomeHips2FitsHeight();
        const previewWidth = skyDomeHips2FitsPreviewWidth();
        const previewHeight = skyDomeHips2FitsPreviewHeight();
        const mediumWidth = skyDomeHips2FitsMediumWidth();
        const mediumHeight = skyDomeHips2FitsMediumHeight();
        const fullUrl = skyDomeHips2FitsUrl(fullWidth, fullHeight);
        const previewUrl = skyDomeHips2FitsUrl(previewWidth, previewHeight);
        const mediumUrl = skyDomeHips2FitsUrl(mediumWidth, mediumHeight);
        const sourceName = `${skyDomeHipsSurveyName()} HiPS`;

        if (previewWidth < fullWidth || previewHeight < fullHeight) {
          const previewStarted = setSkyDomeTextureFromImageUrl(
            previewUrl,
            sourceName,
            projectionName,
            "hips2fits-preview",
            {
              priority: 0,
              message: `Loading ${previewWidth}x${previewHeight} sky preview.`,
              loadedStatus: "loaded-preview",
              loadedMessage: `Loading ${fullWidth}x${fullHeight} high-resolution sky texture.`,
            }
          );
          if (mediumWidth < fullWidth || mediumHeight < fullHeight) {
            setSkyDomeTextureFromImageUrl(
              mediumUrl,
              sourceName,
              projectionName,
              "hips2fits-medium",
              {
                deferInstall: true,
                priority: 1,
                status: "loading-medium",
                message: `Loading ${mediumWidth}x${mediumHeight} sky texture.`,
                loadedStatus: "loaded-medium",
                loadedMessage: `Loading ${fullWidth}x${fullHeight} high-resolution sky texture.`,
                errorStatus: "medium-texture-error",
                errorMessage: `Unable to load the ${mediumWidth}x${mediumHeight} sky texture; keeping the preview sky.`,
              }
            );
          }
          setSkyDomeTextureFromImageUrl(
            fullUrl,
            sourceName,
            projectionName,
            "hips2fits",
            {
              deferInstall: true,
              priority: 2,
              status: "loading-hires",
              message: `Loading ${fullWidth}x${fullHeight} high-resolution sky texture.`,
              errorStatus: "hires-texture-error",
              errorMessage: `Unable to load the ${fullWidth}x${fullHeight} high-resolution sky texture; keeping the best available sky.`,
            }
          );
          return previewStarted;
        }
        return setSkyDomeTextureFromImageUrl(
          fullUrl,
          sourceName,
          projectionName,
          "hips2fits",
          { priority: 2 }
        );
      }

      function initializeSkyDomeFromSceneSpec() {
        if (!skyDomeIsEnabled()) {
          setSkyDomeSnapshotStatus("disabled", "");
          return false;
        }
        if (typeof skyDomeUsesHips2Fits === "function" && skyDomeUsesHips2Fits()) {
          return initializeHips2FitsSkyDome();
        }
        if (skyDomeUsesNativeHips()) {
          return initializeNativeHipsSkyDome();
        }
        const dataUrl = skyDomeImageDataUrl();
        if (!dataUrl) {
          if (skyDomeUsesAladinBackground()) {
            if (skyDomeBackgroundFrameReady) {
              setSkyDomeSnapshotStatus("loaded", "");
              return true;
            }
            setSkyDomeSnapshotStatus(
              "waiting-for-aladin-background",
              "Waiting for Aladin to initialize the live sky background."
            );
            return false;
          }
          setSkyDomeSnapshotStatus(
            "waiting-for-aladin",
            "Waiting for Aladin to capture the sky background."
          );
          return false;
        }
        return setSkyDomeTextureFromDataUrl(
          dataUrl,
          skyDomeLocalSourceName((skySpec && skySpec.survey) || "local"),
          skyDomeSpec.projection || "CAR"
        );
      }

      function updateSkyDome(timestampMs = 0.0) {
        if (!skyDomeMesh) {
          return;
        }
        if (!skyDomeIsEnabled()) {
          skyDomeMesh.visible = false;
          if (root && root.dataset) {
            root.dataset.skyDomeOpacity = "0.000";
            root.dataset.skyDomeForceVisible = skyDomeForceVisible ? "true" : "false";
          }
          if (typeof refreshSkyDomeControlStatus === "function") {
            refreshSkyDomeControlStatus();
          }
          return;
        }
        const opacity = skyDomeOpacityForCurrentView();
        skyDomeMesh.visible = opacity > 0.002;
        skyDomeMesh.position.copy(skyDomeWorldCenterForCurrentView());
        if (root && root.dataset) {
          root.dataset.skyDomeOpacity = opacity.toFixed(3);
          root.dataset.skyDomeForceVisible = skyDomeForceVisible ? "true" : "false";
        }
        if (skyDomeMesh.material) {
          applySkyDomeMaterialUniforms(skyDomeMesh.material);
        }
        setSkyDomeObjectOpacity(opacity);
        if (skyDomeUsesNativeHips() && opacity > 0.002) {
          updateNativeHipsTiles(timestampMs);
        }
      }
""".strip()
