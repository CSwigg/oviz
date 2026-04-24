THREEJS_LEGEND_RUNTIME_JS = """
      let legendGroupChooserGroups = [];
      let legendGroupDropdownOpen = false;

      function closeLegendPopover() {
        activeLegendEditorKey = "";
        if (legendPopoverEl) {
          legendPopoverEl.dataset.open = "false";
          legendPopoverEl.innerHTML = "";
        }
      }

      function setLegendGroupDropdownOpen(isOpen, options = {}) {
        legendGroupDropdownOpen = Boolean(isOpen);
        if (groupDropdownEl) {
          groupDropdownEl.dataset.open = legendGroupDropdownOpen ? "true" : "false";
        }
        if (groupDropdownTriggerEl) {
          groupDropdownTriggerEl.setAttribute("aria-expanded", legendGroupDropdownOpen ? "true" : "false");
        }
        if (legendGroupDropdownOpen && options.focusCurrent) {
          focusLegendGroupOption(String(currentGroup || defaultGroup || ""));
        }
      }

      function syncLegendGroupChooser(options = {}) {
        const activeGroup = String(currentGroup || defaultGroup || "");
        if (groupSelectEl) {
          groupSelectEl.value = activeGroup;
        }
        if (groupDropdownLabelEl) {
          groupDropdownLabelEl.textContent = activeGroup || "Group";
        }
        if (groupDropdownTriggerEl) {
          groupDropdownTriggerEl.disabled = legendGroupChooserGroups.length <= 1;
          groupDropdownTriggerEl.setAttribute("aria-expanded", legendGroupDropdownOpen ? "true" : "false");
        }
        if (groupDropdownEl) {
          groupDropdownEl.dataset.open = legendGroupDropdownOpen ? "true" : "false";
        }
        if (groupDropdownListEl) {
          groupDropdownListEl.querySelectorAll(".oviz-three-group-option").forEach((option) => {
            const isActive = String(option.dataset.group || "") === activeGroup;
            option.dataset.active = isActive ? "true" : "false";
            option.setAttribute("aria-selected", isActive ? "true" : "false");
            option.setAttribute("tabindex", legendGroupDropdownOpen ? "0" : "-1");
          });
        }
      }

      function normalizedLegendGroupOrder(groups) {
        const sourceGroups = Array.isArray(groups) && groups.length ? groups : [defaultGroup];
        const seenGroups = new Set();
        const normalizedGroups = sourceGroups
          .map((groupName) => String(groupName || ""))
          .filter((groupName) => {
            if (!groupName || seenGroups.has(groupName)) {
              return false;
            }
            seenGroups.add(groupName);
            return true;
          });
        const allGroupIndex = normalizedGroups.findIndex((groupName) => groupName.toLowerCase() === "all");
        if (allGroupIndex > 0) {
          const allGroup = normalizedGroups[allGroupIndex];
          normalizedGroups.splice(allGroupIndex, 1);
          normalizedGroups.unshift(allGroup);
        }
        return normalizedGroups;
      }

      function wrappedLegendGroupIndex(index, count) {
        if (count <= 0) {
          return 0;
        }
        return ((index % count) + count) % count;
      }


      function setLegendGroup(nextGroup, options = {}) {
        const requestedGroup = String(nextGroup || "");
        if (!requestedGroup || !Object.prototype.hasOwnProperty.call(groupVisibility, requestedGroup)) {
          syncLegendGroupChooser(options);
          return;
        }
        const interrupt = options.interrupt !== false;
        if (requestedGroup === String(currentGroup || "")) {
          setLegendGroupDropdownOpen(false);
          syncLegendGroupChooser(options);
          return;
        }
        if (interrupt && !actionInterruptsMuted()) {
          interruptActionRun("legend", { disableOrbit: false });
        }
        currentGroup = requestedGroup;
        resetLegendState(currentGroup);
        setLegendGroupDropdownOpen(false);
        syncLegendGroupChooser({ centerActive: true, smooth: true });
        renderLegend();
        renderFrame(currentFrameIndex);
      }

      function shiftLegendGroupByOffset(offset) {
        const count = legendGroupChooserGroups.length;
        if (!count) {
          return;
        }
        const activeGroup = String(currentGroup || "");
        const activeIndex = legendGroupChooserGroups.findIndex((groupName) => String(groupName || "") === activeGroup);
        const nextIndex = wrappedLegendGroupIndex((activeIndex >= 0 ? activeIndex : 0) + offset, count);
        const nextGroup = legendGroupChooserGroups[nextIndex];
        if (nextGroup) {
          setLegendGroup(String(nextGroup || ""), { centerActive: true, smooth: true });
        }
      }

      function focusLegendGroupOption(groupName) {
        if (!groupDropdownListEl) {
          return;
        }
        window.requestAnimationFrame(() => {
          const cleanGroupName = String(groupName || "");
          const option = Array.from(groupDropdownListEl.querySelectorAll(".oviz-three-group-option"))
            .find((candidate) => String(candidate.dataset.group || "") === cleanGroupName)
            || groupDropdownListEl.querySelector(".oviz-three-group-option");
          if (option) {
            option.focus();
          }
        });
      }

      function renderLegendGroupDropdown(groups) {
        if (!legendGroupFieldEl) {
          return;
        }
        legendGroupChooserGroups = normalizedLegendGroupOrder(groups);
        const showChooser = legendGroupChooserGroups.length > 1;
        legendGroupFieldEl.style.display = showChooser ? "" : "none";
        if (groupDropdownEl) {
          groupDropdownEl.style.display = showChooser ? "block" : "none";
        }
        if (groupSelectEl) {
          groupSelectEl.style.display = "none";
        }
        if (groupDropdownListEl) {
          groupDropdownListEl.innerHTML = "";
          legendGroupChooserGroups.forEach((groupName, optionIndex) => {
            const option = document.createElement("button");
            option.type = "button";
            option.className = "oviz-three-group-option";
            option.dataset.group = groupName;
            option.style.setProperty("--group-option-delay", `${Math.min(optionIndex * 26, 112)}ms`);
            option.setAttribute("role", "option");
            option.textContent = groupName;
            option.title = `Show ${groupName}`;
            option.addEventListener("click", (event) => {
              event.preventDefault();
              event.stopPropagation();
              setLegendGroup(groupName, { centerActive: true, smooth: true });
            });
            option.addEventListener("keydown", (event) => {
              if (event.key === "Escape") {
                event.preventDefault();
                setLegendGroupDropdownOpen(false);
                if (groupDropdownTriggerEl) {
                  groupDropdownTriggerEl.focus();
                }
                return;
              }
              if (event.key === "ArrowDown" || event.key === "ArrowUp") {
                event.preventDefault();
                const direction = event.key === "ArrowDown" ? 1 : -1;
                const currentIndex = legendGroupChooserGroups.findIndex((candidate) => String(candidate || "") === String(groupName || ""));
                const nextGroup = legendGroupChooserGroups[wrappedLegendGroupIndex(currentIndex + direction, legendGroupChooserGroups.length)];
                focusLegendGroupOption(nextGroup);
                return;
              }
              if (event.key !== "Enter" && event.key !== " " && event.key !== "Spacebar") {
                return;
              }
              event.preventDefault();
              event.stopPropagation();
              setLegendGroup(groupName, { centerActive: true, smooth: true });
            });
            groupDropdownListEl.appendChild(option);
          });
        }
        if (groupDropdownTriggerEl && !groupDropdownTriggerEl.dataset.bound) {
          groupDropdownTriggerEl.dataset.bound = "true";
          groupDropdownTriggerEl.addEventListener("click", (event) => {
            event.preventDefault();
            event.stopPropagation();
            setLegendGroupDropdownOpen(!legendGroupDropdownOpen, { focusCurrent: !legendGroupDropdownOpen });
            syncLegendGroupChooser();
          });
          groupDropdownTriggerEl.addEventListener("keydown", (event) => {
            if (event.key === "ArrowDown" || event.key === "Enter" || event.key === " " || event.key === "Spacebar") {
              event.preventDefault();
              setLegendGroupDropdownOpen(true, { focusCurrent: true });
              syncLegendGroupChooser();
              return;
            }
            if (event.key === "ArrowUp") {
              event.preventDefault();
              setLegendGroupDropdownOpen(true);
              syncLegendGroupChooser();
              const activeIndex = legendGroupChooserGroups.findIndex((candidate) => String(candidate || "") === String(currentGroup || ""));
              const previousGroup = legendGroupChooserGroups[wrappedLegendGroupIndex((activeIndex >= 0 ? activeIndex : 0) - 1, legendGroupChooserGroups.length)];
              focusLegendGroupOption(previousGroup);
            }
          });
        }
        if (groupSelectEl && !groupSelectEl.dataset.bound) {
          groupSelectEl.dataset.bound = "true";
          groupSelectEl.addEventListener("change", () => {
            setLegendGroup(groupSelectEl.value, { centerActive: true, smooth: true });
          });
        }
        if (groupDropdownEl && !groupDropdownEl.dataset.closeBound) {
          groupDropdownEl.dataset.closeBound = "true";
          window.addEventListener("pointerdown", (event) => {
            if (!legendGroupDropdownOpen || (groupDropdownEl && groupDropdownEl.contains(event.target))) {
              return;
            }
            setLegendGroupDropdownOpen(false);
            syncLegendGroupChooser();
          });
        }
        syncLegendGroupChooser({ centerActive: true, smooth: false });
      }

      function setLegendPanelOpen(isOpen) {
        legendPanelOpen = minimalModeEnabled ? true : Boolean(isOpen);
        if (legendPanelEl) {
          legendPanelEl.dataset.open = legendPanelOpen ? "true" : "false";
        }
        if (legendPanelBodyEl) {
          legendPanelBodyEl.style.display = legendPanelOpen ? "flex" : "none";
        }
        if (legendPanelToggleEl) {
          legendPanelToggleEl.textContent = legendPanelOpen ? "▾" : "▸";
          legendPanelToggleEl.setAttribute(
            "title",
            legendPanelOpen ? "Collapse the legend" : "Expand the legend"
          );
          legendPanelToggleEl.setAttribute("aria-expanded", legendPanelOpen ? "true" : "false");
        }
        if (legendPanelEl) {
          const rect = legendPanelRectState || defaultLegendPanelRect();
          applyLegendPanelRect(rect);
        }
        if (!legendPanelOpen) {
          closeLegendPopover();
        }
        legendResizeEls.forEach((handle) => {
          handle.style.display = legendPanelOpen ? "" : "none";
        });
      }

      function visibleLegendItemsForCurrentGroup() {
        const defaults = groupDefaults(currentGroup);
        return legendItems.filter((item) => {
          const mode = defaults[item.key];
          return !(mode === false || mode === undefined);
        });
      }

      function collapseLegendEntryControls(entry, onCollapsed) {
        const controls = entry ? entry.querySelector(".oviz-three-legend-controls") : null;
        if (!controls) {
          onCollapsed();
          return;
        }
        controls.dataset.visible = "false";
        entry.dataset.editorOpen = "closing";
        window.setTimeout(onCollapsed, 280);
      }

      function renderLegend() {
        if (!legendTraceListEl || !legendVolumeListEl) {
          return;
        }
        syncLegendGroupChooser();
        const items = visibleLegendItemsForCurrentGroup();
        const traceItems = minimalModeEnabled
          ? items
          : items.filter((item) => !volumeLayerForKey(item.key));
        const volumeItems = minimalModeEnabled
          ? []
          : items.filter((item) => volumeLayerForKey(item.key));
        legendTraceListEl.innerHTML = "";
        legendVolumeListEl.innerHTML = "";
        legendEditButtonByKey.clear();

        if (legendTraceSectionEl) {
          legendTraceSectionEl.dataset.empty = traceItems.length ? "false" : "true";
        }
        if (legendVolumeSectionEl) {
          legendVolumeSectionEl.dataset.empty = volumeItems.length ? "false" : "true";
        }
        applyLegendSectionState("traces");
        applyLegendSectionState("volumes");

        if (!traceItems.length && !volumeItems.length) {
          const emptyState = document.createElement("div");
          emptyState.className = "oviz-three-legend-title";
          emptyState.textContent = "No visible traces for this group.";
          legendTraceListEl.appendChild(emptyState);
          closeLegendPopover();
          return;
        }

        let activeEditorStillVisible = false;

        function appendLegendEntries(targetEl, sectionItems) {
          sectionItems.forEach((item) => {
            const itemKey = String(item.key);
            const entry = document.createElement("div");
            const itemColor = legendColorForItem(item);
            const editorOpen = activeLegendEditorKey === itemKey;
            entry.className = "oviz-three-legend-entry";
            entry.dataset.active = legendState[itemKey] ? "true" : "false";
            entry.dataset.editorOpen = editorOpen ? "true" : "false";

            const row = document.createElement("div");
            row.className = "oviz-three-legend-row";

            const toggleButton = document.createElement("button");
            toggleButton.type = "button";
            toggleButton.className = "oviz-three-legend-item";
            toggleButton.dataset.active = legendState[itemKey] ? "true" : "false";
            toggleButton.title = minimalModeEnabled
              ? String(item.name || "")
              : `${String(item.name || "")}: click to show or hide`;
            toggleButton.style.color = itemColor;

            const swatch = document.createElement("span");
            swatch.className = "oviz-three-legend-swatch";
            swatch.style.background = itemColor;
            toggleButton.appendChild(swatch);

            const meta = document.createElement("span");
            meta.className = "oviz-three-legend-meta";

            const name = document.createElement("span");
            name.className = "oviz-three-legend-name";
            name.textContent = String(item.name || "");
            meta.appendChild(name);
            toggleButton.appendChild(meta);

            toggleButton.addEventListener("click", () => {
              legendState[itemKey] = !legendState[itemKey];
              renderLegend();
              renderFrame(currentFrameIndex);
            });
            if (!minimalModeEnabled && !volumeLayerForKey(item.key)) {
              toggleButton.addEventListener("dblclick", () => {
                const traceItems = items.filter((candidate) => !volumeLayerForKey(candidate.key));
                const onlyThisVisible = traceItems.every((candidate) => {
                  const candidateKey = String(candidate.key);
                  return legendState[candidateKey] === (candidateKey === itemKey);
                });
                traceItems.forEach((candidate) => {
                  const candidateKey = String(candidate.key);
                  legendState[candidateKey] = onlyThisVisible ? true : candidateKey === itemKey;
                });
                renderLegend();
                renderFrame(currentFrameIndex);
              });
            }
            if (!minimalModeEnabled) {
              const editButton = document.createElement("button");
              editButton.type = "button";
              editButton.className = "oviz-three-legend-edit";
              editButton.dataset.open = editorOpen ? "true" : "false";
              editButton.textContent = editorOpen ? "⌄" : "›";
              editButton.style.color = itemColor;
              editButton.style.borderColor = itemColor;
              editButton.title = editorOpen
                ? "Collapse legend controls"
                : `Expand ${String(item.name || "item")} controls`;
              editButton.addEventListener("click", (event) => {
                event.preventDefault();
                event.stopPropagation();
                if (activeLegendEditorKey === itemKey) {
                  activeLegendEditorKey = "";
                  editButton.dataset.open = "false";
                  editButton.textContent = "›";
                  collapseLegendEntryControls(entry, () => {
                    closeLegendPopover();
                    renderLegend();
                  });
                  return;
                }
                activeLegendEditorKey = itemKey;
                renderLegend();
              });
              row.appendChild(editButton);
              if (editorOpen) {
                const controls = volumeLayerForKey(item.key)
                  ? buildVolumeLegendControls(item, toggleButton)
                  : buildTraceLegendControls(item, toggleButton);
                if (controls) {
                  activeEditorStillVisible = true;
                  controls.dataset.visible = "false";
                  entry.appendChild(controls);
                  window.requestAnimationFrame(() => {
                    if (entry.isConnected && entry.dataset.editorOpen === "true") {
                      controls.dataset.visible = "true";
                    }
                  });
                }
              }
            }
            row.appendChild(toggleButton);
            entry.insertBefore(row, entry.firstChild);
            targetEl.appendChild(entry);
          });
        }

        appendLegendEntries(legendTraceListEl, traceItems);
        appendLegendEntries(legendVolumeListEl, volumeItems);

        if (!activeEditorStillVisible) {
          closeLegendPopover();
        }
        if (legendPanelOpen && legendPanelEl) {
          applyLegendPanelRect(legendPanelRectState || defaultLegendPanelRect());
        }
      }
""".strip()
