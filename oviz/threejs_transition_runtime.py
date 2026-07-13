from __future__ import annotations


THREEJS_TRANSITION_RUNTIME_JS = r"""
      const OvizTransitionRuntime = (() => {
        let sessionSerial = 0;

        function clamp01(value) {
          return Math.min(Math.max(Number(value) || 0.0, 0.0), 1.0);
        }

        function easingValue(name, value) {
          const t = clamp01(value);
          const easing = String(name || "easeInOutCubic");
          if (easing === "linear") return t;
          if (easing === "easeInCubic") return t * t * t;
          if (easing === "easeOutCubic") return 1.0 - Math.pow(1.0 - t, 3.0);
          return t < 0.5
            ? 4.0 * t * t * t
            : 1.0 - Math.pow(-2.0 * t + 2.0, 3.0) / 2.0;
        }

        function own(object, key) {
          return Boolean(
            object
            && typeof object === "object"
            && Object.prototype.hasOwnProperty.call(object, key)
          );
        }

        function resolveTraceVisibility(trace, defaults, overrides) {
          if (!trace || !trace.key) return false;
          const key = String(trace.key);
          const mode = defaults && typeof defaults === "object" ? defaults[key] : undefined;
          if (mode === false || mode === undefined) return false;
          if (trace.showlegend) {
            return own(overrides, key) ? Boolean(overrides[key]) : mode === true;
          }
          return mode === true;
        }

        function traceMap(traces) {
          const result = new Map();
          (Array.isArray(traces) ? traces : []).forEach((trace) => {
            const key = String(trace && trace.key || "").trim();
            if (!key || result.has(key)) return;
            result.set(key, trace);
          });
          return result;
        }

        function buildTraceUnion(fromTraces, toTraces) {
          const fromByKey = traceMap(fromTraces);
          const toByKey = traceMap(toTraces);
          const orderedKeys = [];
          fromByKey.forEach((_trace, key) => orderedKeys.push(key));
          toByKey.forEach((_trace, key) => {
            if (!fromByKey.has(key)) orderedKeys.push(key);
          });
          return orderedKeys.map((key) => ({
            key,
            from: fromByKey.get(key) || null,
            to: toByKey.get(key) || null,
          }));
        }

        function normalizedIdentity(value) {
          return String(value === undefined || value === null ? "" : value)
            .trim()
            .toLowerCase()
            .replace(/\s+/g, "_");
        }

        function stablePointIdentity(point) {
          if (!point || typeof point !== "object") return "";
          const motionKey = normalizedIdentity(point.motion && point.motion.key);
          if (motionKey) return `motion:${motionKey}`;
          const selection = point.selection && typeof point.selection === "object"
            ? point.selection
            : {};
          const selectionKey = normalizedIdentity(
            selection.selection_key
            ?? selection.key
            ?? selection.cluster_key
            ?? selection.cluster_name
            ?? point.selection_key
          );
          return selectionKey ? `selection:${selectionKey}` : "";
        }

        function pointEntries(points) {
          const occurrences = new Map();
          return (Array.isArray(points) ? points : []).map((point, index) => {
            const base = stablePointIdentity(point);
            if (!base) return { point, index, base: "", key: "" };
            const occurrence = occurrences.get(base) || 0;
            occurrences.set(base, occurrence + 1);
            return {
              point,
              index,
              base,
              key: occurrence ? `${base}#${occurrence}` : base,
            };
          });
        }

        function matchTracePoints(traceKey, fromPoints, toPoints) {
          const fromEntries = pointEntries(fromPoints);
          const toEntries = pointEntries(toPoints);
          const fromAllStable = fromEntries.every((entry) => Boolean(entry.key));
          const toAllStable = toEntries.every((entry) => Boolean(entry.key));
          const compatibleIndexTopology = (
            fromEntries.length === toEntries.length
            && !fromEntries.some((entry) => entry.key)
            && !toEntries.some((entry) => entry.key)
          );
          if (compatibleIndexTopology) {
            return fromEntries.map((entry, index) => ({
              key: `${String(traceKey)}:index:${index}`,
              from: entry.point,
              to: toEntries[index].point,
              fromIndex: index,
              toIndex: index,
              matchedBy: "index",
            }));
          }

          const fromByKey = new Map();
          const toByKey = new Map();
          fromEntries.forEach((entry) => {
            const key = entry.key || `${String(traceKey)}:from:${entry.index}`;
            fromByKey.set(key, entry);
          });
          toEntries.forEach((entry) => {
            const key = entry.key || `${String(traceKey)}:to:${entry.index}`;
            toByKey.set(key, entry);
          });
          const orderedKeys = [];
          fromByKey.forEach((_entry, key) => orderedKeys.push(key));
          toByKey.forEach((_entry, key) => {
            if (!fromByKey.has(key)) orderedKeys.push(key);
          });
          return orderedKeys.map((key) => {
            const from = fromByKey.get(key) || null;
            const to = toByKey.get(key) || null;
            return {
              key,
              from: from && from.point,
              to: to && to.point,
              fromIndex: from ? from.index : -1,
              toIndex: to ? to.index : -1,
              matchedBy: from && to && from.key && to.key && fromAllStable && toAllStable
                ? "stable"
                : (from && to && from.key && to.key ? "stable" : "crossfade"),
            };
          });
        }

        function membershipContains(membership, key, point) {
          if (membership === null || membership === undefined || membership === "all") return true;
          if (typeof membership === "function") return Boolean(membership(key, point));
          if (membership instanceof Set || membership instanceof Map) return membership.has(key);
          if (Array.isArray(membership)) return membership.includes(key);
          return own(membership, key) ? Boolean(membership[key]) : false;
        }

        function membershipWeight(key, point, fromMembership, toMembership, progress) {
          const from = membershipContains(fromMembership, key, point) ? 1.0 : 0.0;
          const to = membershipContains(toMembership, key, point) ? 1.0 : 0.0;
          const t = clamp01(progress);
          return from + (to - from) * t;
        }

        function normalizePhases(phases) {
          const result = [];
          let cursor = 0.0;
          (Array.isArray(phases) ? phases : []).forEach((phase, index) => {
            const durationMs = Math.max(Number(phase && phase.durationMs) || 0.0, 0.0);
            result.push({
              name: String(phase && phase.name || `phase-${index + 1}`),
              durationMs,
              startMs: cursor,
              endMs: cursor + durationMs,
            });
            cursor += durationMs;
          });
          if (!result.length) {
            result.push({ name: "settle", durationMs: 0.0, startMs: 0.0, endMs: 0.0 });
          }
          return { phases: result, effectiveDurationMs: cursor };
        }

        function phaseAt(plan, elapsedMs) {
          const elapsed = Math.min(Math.max(Number(elapsedMs) || 0.0, 0.0), plan.effectiveDurationMs);
          let phaseIndex = plan.phases.findIndex((item) => elapsed <= item.endMs + 1e-9);
          if (phaseIndex < 0) phaseIndex = plan.phases.length - 1;
          const phase = plan.phases[phaseIndex];
          const duration = Math.max(phase.durationMs, 1e-9);
          return {
            ...phase,
            index: phaseIndex,
            progress: phase.durationMs <= 0.0 ? 1.0 : clamp01((elapsed - phase.startMs) / duration),
          };
        }

        function abortControllerForRun() {
          if (typeof AbortController === "function") return new AbortController();
          const signal = { aborted: false, reason: undefined };
          return {
            signal,
            abort(reason) {
              signal.aborted = true;
              signal.reason = reason;
            },
          };
        }

        function createSession(adapter = {}) {
          const now = typeof adapter.now === "function" ? adapter.now : (() => Date.now());
          const requestFrame = typeof adapter.requestFrame === "function" ? adapter.requestFrame : null;
          const cancelFrame = typeof adapter.cancelFrame === "function" ? adapter.cancelFrame : null;
          let active = null;
          let scheduledFrame = null;

          function emit(type, detail) {
            if (typeof adapter.emit === "function") adapter.emit(type, detail);
          }

          function callPhase(name, progress, context) {
            const callback = name === "camera"
              ? adapter.applyCamera
              : (name === "appearance"
                ? adapter.applyAppearance
                : (name === "time" ? adapter.applyTime : null));
            if (typeof callback === "function") callback(progress, context);
          }

          function terminalDetail(run, extra = {}) {
            return {
              owner: run.owner,
              runId: run.runId,
              phase: run.phase,
              phaseProgress: run.phaseProgress,
              effectiveDurationMs: run.plan.effectiveDurationMs,
              ...extra,
            };
          }

          function runIsCurrent(run) {
            return Boolean(
              run
              && run === active
              && !run.cancelled
              && !(run.signal && run.signal.aborted)
            );
          }

          function phaseContext(run, phase, rawProgress, easedProgress) {
            run.phase = phase.name;
            run.phaseProgress = phase.progress;
            return terminalDetail(run, {
              target: run.target,
              rawProgress,
              easedProgress,
              source: run.source,
            });
          }

          function applyCrossedPhaseTerminals(run, phase, rawProgress) {
            for (let index = 0; index < phase.index; index += 1) {
              if (run.terminalPhaseIndexes.has(index)) continue;
              const crossed = run.plan.phases[index];
              const terminal = { ...crossed, index, progress: 1.0 };
              const context = phaseContext(run, terminal, rawProgress, 1.0);
              callPhase(crossed.name, 1.0, context);
              run.terminalPhaseIndexes.add(index);
            }
          }

          function settle(run) {
            if (!runIsCurrent(run) || run.settling) return;
            run.settling = true;
            Promise.resolve(typeof adapter.finalize === "function" ? adapter.finalize(run.target, run) : undefined)
              .then(() => {
                if (!runIsCurrent(run)) return undefined;
                return Promise.resolve(
                  typeof adapter.validate === "function"
                    ? adapter.validate(run.target, run)
                    : undefined
                );
              })
              .then((validation) => {
                if (!runIsCurrent(run)) return;
                active = null;
                emit("transition-complete", terminalDetail(run, { validation }));
                run.resolve({ cancelled: false, owner: run.owner, runId: run.runId, validation });
              })
              .catch((error) => {
                if (!runIsCurrent(run)) return;
                active = null;
                emit("transition-error", terminalDetail(run, { message: String(error && error.message || error) }));
                run.reject(error);
              });
          }

          function tick(timestamp) {
            scheduledFrame = null;
            const run = active;
            if (!run || run.cancelled || run.settling) return false;
            try {
              const currentNow = Number.isFinite(Number(timestamp)) ? Number(timestamp) : Number(now());
              const elapsedMs = Math.max(0.0, currentNow - run.startedAt);
              const phase = phaseAt(run.plan, elapsedMs);
              const rawProgress = run.plan.effectiveDurationMs <= 0.0
                ? 1.0
                : clamp01(elapsedMs / run.plan.effectiveDurationMs);
              applyCrossedPhaseTerminals(run, phase, rawProgress);
              const eased = easingValue(run.easing, phase.progress);
              const context = phaseContext(run, phase, rawProgress, eased);
              callPhase(phase.name, eased, context);
              if (phase.progress >= 1.0) run.terminalPhaseIndexes.add(phase.index);
              emit("transition-progress", context);
              if (rawProgress >= 1.0) {
                settle(run);
              } else if (requestFrame) {
                scheduledFrame = requestFrame(tick);
              }
              return true;
            } catch (error) {
              if (run === active) active = null;
              emit("transition-error", terminalDetail(run, { message: String(error && error.message || error) }));
              run.reject(error);
              return false;
            }
          }

          function cancel(reason = "cancelled") {
            const run = active;
            if (!run) return false;
            active = null;
            run.cancelled = true;
            if (run.abortController) run.abortController.abort(reason);
            if (scheduledFrame !== null && cancelFrame) cancelFrame(scheduledFrame);
            scheduledFrame = null;
            let cleanupError = null;
            if (typeof adapter.cancel === "function") {
              try {
                adapter.cancel(reason, run);
              } catch (error) {
                cleanupError = error;
                emit("transition-error", terminalDetail(run, {
                  stage: "cancel",
                  message: String(error && error.message || error),
                }));
              }
            }
            emit("transition-cancel", terminalDetail(run, {
              reason,
              cleanupError: cleanupError ? String(cleanupError && cleanupError.message || cleanupError) : "",
            }));
            run.resolve({
              cancelled: true,
              owner: run.owner,
              runId: run.runId,
              reason,
              cleanupError: cleanupError ? String(cleanupError && cleanupError.message || cleanupError) : "",
            });
            return true;
          }

          function start(spec = {}) {
            if (active) cancel("retargeted");
            const plan = normalizePhases(spec.phases);
            const startedAt = Number.isFinite(Number(spec.startedAt)) ? Number(spec.startedAt) : Number(now());
            let resolve;
            let reject;
            const promise = new Promise((onResolve, onReject) => {
              resolve = onResolve;
              reject = onReject;
            });
            const abortController = abortControllerForRun();
            const run = {
              owner: String(spec.owner || "transition"),
              runId: String(spec.runId || `transition-${++sessionSerial}`),
              target: spec.target,
              source: spec.source !== undefined
                ? spec.source
                : (typeof adapter.capture === "function" ? adapter.capture() : undefined),
              easing: String(spec.easing || "easeInOutCubic"),
              plan,
              startedAt,
              phase: plan.phases[0].name,
              phaseProgress: 0.0,
              cancelled: false,
              settling: false,
              terminalPhaseIndexes: new Set(),
              abortController,
              signal: abortController.signal,
              resolve,
              reject,
              promise,
            };
            run.isCurrent = () => runIsCurrent(run);
            run.commit = (callback) => {
              if (!runIsCurrent(run) || typeof callback !== "function") return false;
              callback();
              return true;
            };
            active = run;
            try {
              if (typeof adapter.prepare === "function") adapter.prepare(run.target, run);
              emit("transition-start", terminalDetail(run, { target: run.target }));
              if (requestFrame) scheduledFrame = requestFrame(tick);
            } catch (error) {
              if (run === active) active = null;
              run.cancelled = true;
              abortController.abort("prepare-error");
              emit("transition-error", terminalDetail(run, {
                stage: "prepare",
                message: String(error && error.message || error),
              }));
              run.reject(error);
            }
            return promise;
          }

          function retarget(spec = {}) {
            const source = typeof adapter.capture === "function" ? adapter.capture() : undefined;
            cancel("retargeted");
            return start({ ...spec, source });
          }

          return {
            start,
            tick,
            cancel,
            retarget,
            current: () => active,
          };
        }

        function createFrameCoordinator(adapter = {}) {
          let lastFrameToken = null;
          let lastResult = false;
          let activeOwner = "";
          let activeRunId = "";
          let updateCount = 0;

          function update(frameToken, descriptor = {}, callback = null) {
            const token = frameToken === undefined || frameToken === null
              ? `frame-${updateCount + 1}`
              : String(frameToken);
            if (token === lastFrameToken) return lastResult;
            lastFrameToken = token;
            activeOwner = String(descriptor.owner || "idle");
            activeRunId = String(descriptor.runId || "");
            updateCount += 1;
            try {
              lastResult = typeof callback === "function"
                ? Boolean(callback({
                  owner: activeOwner,
                  runId: activeRunId,
                  frameToken: token,
                  updateCount,
                }))
                : false;
              return lastResult;
            } catch (error) {
              lastResult = false;
              if (typeof adapter.onError === "function") {
                try {
                  adapter.onError(error, {
                    owner: activeOwner,
                    runId: activeRunId,
                    frameToken: token,
                    updateCount,
                  });
                } catch (_handlerError) {
                }
              }
              return false;
            }
          }

          function snapshot() {
            return {
              owner: activeOwner,
              runId: activeRunId,
              frameToken: lastFrameToken,
              updateCount,
            };
          }

          return { update, snapshot };
        }

        return {
          clamp01,
          easingValue,
          resolveTraceVisibility,
          buildTraceUnion,
          stablePointIdentity,
          matchTracePoints,
          membershipWeight,
          createSession,
          createFrameCoordinator,
        };
      })();
      globalThis.OvizTransitionRuntime = OvizTransitionRuntime;
""".strip()
