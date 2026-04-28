THREEJS_LEGEND_RUNTIME_JS = """
      function closeLegendPopover() {
        activeLegendEditorKey = "";
        if (legendPopoverEl) {
          legendPopoverEl.dataset.open = "false";
          legendPopoverEl.innerHTML = "";
        }
      }

      function positionLegendPopover(anchorEl) {
        if (!legendPopoverEl || !anchorEl) {
          return;
        }
        const rootRect = root.getBoundingClientRect();
        const anchorRect = anchorEl.getBoundingClientRect();
        const legendRect = legendPanelEl ? legendPanelEl.getBoundingClientRect() : anchorRect;
        const popoverWidth = legendPopoverEl.offsetWidth || 244;
        const popoverHeight = legendPopoverEl.offsetHeight || 320;
        let left = Math.max(anchorRect.right, legendRect.right) - rootRect.left + 12;
        if (left + popoverWidth > root.clientWidth - 12) {
          left = anchorRect.left - rootRect.left - popoverWidth - 12;
        }
        left = clampRange(left, 12, Math.max(12, root.clientWidth - popoverWidth - 12));
        const top = clampRange(
          anchorRect.top - rootRect.top - 6,
          12,
          Math.max(12, root.clientHeight - popoverHeight - 12)
        );
        legendPopoverEl.style.left = `${left}px`;
        legendPopoverEl.style.top = `${top}px`;
      }

      function renderLegendPopover(item, toggleButton) {
        if (minimalModeEnabled) {
          closeLegendPopover();
          return;
        }
        if (!legendPopoverEl || !item || !toggleButton) {
          closeLegendPopover();
          return;
        }
        legendPopoverEl.innerHTML = "";
        const head = document.createElement("div");
        head.className = "oviz-three-legend-popover-head";

        const title = document.createElement("div");
        title.className = "oviz-three-legend-popover-title";
        const eyebrow = document.createElement("div");
        eyebrow.className = "oviz-three-legend-popover-eyebrow";
        eyebrow.textContent = volumeLayerForKey(item.key) ? "Volume" : "Trace";
        const name = document.createElement("div");
        name.className = "oviz-three-legend-popover-name";
        name.textContent = String(item.name || "");
        name.style.color = legendColorForItem(item);
        title.appendChild(eyebrow);
        title.appendChild(name);
        head.appendChild(title);

        const closeButton = document.createElement("button");
        closeButton.type = "button";
        closeButton.className = "oviz-three-legend-popover-close";
        closeButton.textContent = "×";
        closeButton.title = "Close legend editor";
        closeButton.addEventListener("click", () => {
          closeLegendPopover();
          renderLegend();
        });
        head.appendChild(closeButton);
        legendPopoverEl.appendChild(head);

        const controls = volumeLayerForKey(item.key)
          ? buildVolumeLegendControls(item, toggleButton)
          : buildTraceLegendControls(item, toggleButton);
        if (controls) {
          setLegendPanelOpen(true);
          legendPopoverEl.appendChild(controls);
          legendPopoverEl.dataset.open = "true";
          window.requestAnimationFrame(() => positionLegendPopover(toggleButton));
        } else {
          closeLegendPopover();
        }
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

      function renderLegend() {
        if (!legendTraceListEl || !legendVolumeListEl) {
          return;
        }
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

        function appendLegendEntries(targetEl, sectionItems) {
          sectionItems.forEach((item) => {
          const itemKey = String(item.key);
          const entry = document.createElement("div");
          const itemColor = legendColorForItem(item);
          entry.className = "oviz-three-legend-entry";
          entry.dataset.active = legendState[itemKey] ? "true" : "false";
          entry.dataset.editorOpen = activeLegendEditorKey === itemKey ? "true" : "false";
          entry.style.borderColor = legendState[itemKey] ? itemColor : "rgba(255, 255, 255, 0.06)";

          const toggleButton = document.createElement("button");
          toggleButton.type = "button";
          toggleButton.className = "oviz-three-legend-item";
          toggleButton.dataset.active = legendState[itemKey] ? "true" : "false";
          toggleButton.title = minimalModeEnabled
            ? String(item.name || "")
            : `${String(item.name || "")}: click to show or hide`;

          const swatch = document.createElement("span");
          swatch.className = "oviz-three-legend-swatch";
          swatch.style.background = itemColor;
          toggleButton.appendChild(swatch);

          const meta = document.createElement("span");
          meta.className = "oviz-three-legend-meta";

          const name = document.createElement("span");
          name.className = "oviz-three-legend-name";
          name.textContent = String(item.name || "");
          name.style.color = itemColor;
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
          entry.appendChild(toggleButton);

          if (!minimalModeEnabled) {
            const editButton = document.createElement("button");
            editButton.type = "button";
            editButton.className = "oviz-three-legend-edit";
            editButton.textContent = activeLegendEditorKey === itemKey ? "×" : "›";
            editButton.title = activeLegendEditorKey === itemKey
              ? "Close legend editor"
              : `Edit ${String(item.name || "item")}`;
            editButton.addEventListener("click", (event) => {
              event.preventDefault();
              event.stopPropagation();
              if (activeLegendEditorKey === itemKey) {
                closeLegendPopover();
                renderLegend();
                return;
              }
              activeLegendEditorKey = itemKey;
              renderLegend();
            });
            entry.appendChild(editButton);
            legendEditButtonByKey.set(itemKey, editButton);
          }
            targetEl.appendChild(entry);
          });
        }

        appendLegendEntries(legendTraceListEl, traceItems);
        appendLegendEntries(legendVolumeListEl, volumeItems);

        if (activeLegendEditorKey && legendEditButtonByKey.has(activeLegendEditorKey)) {
          const activeItem = items.find((item) => String(item.key) === activeLegendEditorKey);
          renderLegendPopover(activeItem, legendEditButtonByKey.get(activeLegendEditorKey));
        } else {
          closeLegendPopover();
        }
        if (legendPanelOpen && legendPanelEl) {
          applyLegendPanelRect(legendPanelRectState || defaultLegendPanelRect());
        }
      }
""".strip()
