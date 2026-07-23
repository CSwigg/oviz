"""Reveal.js deck authoring and presentation runtime for Oviz."""

THREEJS_DECK_RUNTIME_JS = r"""
      const OVIZ_DECK_VERSION = 2;
      const OVIZ_DECK_WIDTH = 1600;
      const OVIZ_DECK_HEIGHT = 900;
      const OVIZ_DECK_REVEAL_VERSION = "5.2.1";
      let ovizDeckProject = null;
      let ovizDeckEditorOpen = false;
      let ovizDeckPresenting = false;
      let ovizDeckActiveIndex = 0;
      let ovizDeckSelectedBlockId = "";
      let ovizDeckSelectedObjectIds = new Set();
      let ovizDeckObjectClipboard = [];
      let ovizDeckAddMenuOpen = false;
      let ovizDeckObjectsPanelOpen = false;
      let ovizDeckMarqueeState = null;
      let ovizDeckGuideEls = [];
      let ovizDeckDirty = false;
      let ovizDeckReveal = null;
      let ovizDeckRevealPromise = null;
      let ovizDeckRenderGeneration = 0;
      let ovizDeckDragState = null;
      let ovizDeckButtonEl = null;
      let ovizDeckEditorEl = null;
      let ovizDeckAuthoringLayerEl = null;
      let ovizDeckRevealRootEl = null;
      let ovizDeckStatusEl = null;
      let ovizDeckUndoStack = [];
      let ovizDeckRedoStack = [];
      let ovizDeckHistoryPresent = null;
      let ovizDeckLastHistoryKey = "";
      let ovizDeckLastHistoryAt = 0;

      function ovizDeckClone(value, fallback = null) {
        try {
          return JSON.parse(JSON.stringify(value));
        } catch (_err) {
          return fallback;
        }
      }

      function ovizDeckUuid(prefix = "deck") {
        if (window.crypto && typeof window.crypto.randomUUID === "function") {
          return `${prefix}-${window.crypto.randomUUID()}`;
        }
        return `${prefix}-${Date.now().toString(36)}-${Math.random().toString(36).slice(2, 10)}`;
      }

      function ovizDeckFinite(value, fallback) {
        const number = Number(value);
        return Number.isFinite(number) ? number : Number(fallback);
      }

      function ovizDeckClamp(value, minimum, maximum, fallback) {
        return Math.min(Math.max(ovizDeckFinite(value, fallback), minimum), maximum);
      }

      function ovizDeckNormalizeColor(value, fallback = "#ffffff") {
        const text = String(value || "").trim();
        return /^#[0-9a-f]{6}$/i.test(text) ? text.toLowerCase() : fallback;
      }

      function ovizDeckNormalizeBlock(value, index = 0, seenIds = new Set(), legacyPercent = false) {
        const source = value && typeof value === "object" ? value : {};
        let id = String(source.id || "").trim() || ovizDeckUuid("object");
        while (seenIds.has(id)) id = ovizDeckUuid("object");
        seenIds.add(id);
        const kind = ["text", "shape"].includes(String(source.kind || "").toLowerCase())
          ? String(source.kind).toLowerCase()
          : "text";
        const type = ["title", "subtitle", "text"].includes(String(source.type || "").toLowerCase())
          ? String(source.type).toLowerCase()
          : "text";
        const shapeTypes = ["rectangle", "rounded_rectangle", "ellipse", "triangle", "diamond", "line", "arrow"];
        const shapeType = shapeTypes.includes(String(source.shape_type || source.shape || "").toLowerCase())
          ? String(source.shape_type || source.shape).toLowerCase()
          : "rectangle";
        const defaultSize = type === "title" ? 64 : (type === "subtitle" ? 34 : 26);
        const defaultWeight = type === "title" ? 700 : (type === "subtitle" ? 500 : 400);
        const align = ["left", "center", "right"].includes(String(source.align || "").toLowerCase())
          ? String(source.align).toLowerCase()
          : "left";
        const fontFamily = String(source.font_family || "sans").trim().slice(0, 100) || "sans";
        const fontStyle = String(source.font_style || "").toLowerCase() === "italic" ? "italic" : "normal";
        const listStyle = ["none", "bullet", "number"].includes(String(source.list_style || "").toLowerCase())
          ? String(source.list_style).toLowerCase()
          : "none";
        const legacyX = ovizDeckFinite(source.x, 7);
        const legacyY = ovizDeckFinite(source.y, 14 + index * 12);
        const legacyWidth = ovizDeckFinite(source.width, kind === "shape" ? 18 : 64);
        const defaultHeight = kind === "shape" ? 180 : (type === "title" ? 120 : (type === "subtitle" ? 76 : 110));
        const x = legacyPercent ? legacyX * 16 : ovizDeckFinite(source.x, 112);
        const y = legacyPercent ? legacyY * 9 : ovizDeckFinite(source.y, 126 + index * 90);
        const width = legacyPercent ? legacyWidth * 16 : ovizDeckFinite(source.width, kind === "shape" ? 280 : 1024);
        const height = legacyPercent && source.height !== undefined
          ? ovizDeckFinite(source.height, defaultHeight / 9) * 9
          : ovizDeckFinite(source.height, defaultHeight);
        const defaultName = kind === "text" && String(source.text || "").trim()
          ? String(source.text).trim().split(/\r?\n/)[0].slice(0, 36)
          : shapeType.replace(/_/g, " ").replace(/\b\w/g, (letter) => letter.toUpperCase());
        const borderStyle = ["solid", "dashed", "dotted"].includes(String(source.border_style || "").toLowerCase())
          ? String(source.border_style).toLowerCase()
          : "solid";
        return {
          id,
          name: String(source.name || defaultName || `Object ${index + 1}`).trim().slice(0, 80) || `Object ${index + 1}`,
          kind,
          type,
          shape_type: shapeType,
          text: String(source.text || ""),
          x: ovizDeckClamp(x, 0, OVIZ_DECK_WIDTH - 20, 112),
          y: ovizDeckClamp(y, 0, OVIZ_DECK_HEIGHT - 20, 126 + index * 90),
          width: ovizDeckClamp(width, 20, OVIZ_DECK_WIDTH, kind === "shape" ? 280 : 1024),
          height: ovizDeckClamp(height, 20, OVIZ_DECK_HEIGHT, defaultHeight),
          rotation: ovizDeckClamp(source.rotation, -360, 360, 0),
          opacity: ovizDeckClamp(source.opacity, 0, 1, 1),
          locked: Boolean(source.locked),
          group_id: String(source.group_id || "").trim() || null,
          font_size: ovizDeckClamp(source.font_size, 10, 240, defaultSize),
          font_weight: Math.round(ovizDeckClamp(source.font_weight, 100, 900, defaultWeight)),
          font_family: fontFamily,
          font_style: fontStyle,
          underline: Boolean(source.underline),
          strikethrough: Boolean(source.strikethrough),
          line_height: ovizDeckClamp(source.line_height, 0.8, 2.4, 1.08),
          character_spacing: ovizDeckClamp(source.character_spacing, -5, 30, 0),
          list_style: listStyle,
          color: ovizDeckNormalizeColor(source.color),
          align,
          text_outline: Boolean(source.text_outline),
          text_outline_color: ovizDeckNormalizeColor(source.text_outline_color, "#000000"),
          text_outline_width: ovizDeckClamp(source.text_outline_width, 0, 8, 1),
          text_shadow: Boolean(source.text_shadow),
          fill_color: ovizDeckNormalizeColor(source.fill_color, "#2b6cb0"),
          fill_opacity: ovizDeckClamp(source.fill_opacity, 0, 1, 0.7),
          border_color: ovizDeckNormalizeColor(source.border_color, "#ffffff"),
          border_width: ovizDeckClamp(source.border_width, 0, 24, 2),
          border_style: borderStyle,
          corner_radius: ovizDeckClamp(source.corner_radius, 0, 200, 24),
          shadow: Boolean(source.shadow),
        };
      }

      function ovizDeckNormalizeProject(value) {
        const source = value && typeof value === "object" ? value : {};
        const rawSlides = Array.isArray(source.slides) ? source.slides : [];
        const seenSlides = new Set();
        const seenBlocks = new Set();
        const reveal = source.reveal && typeof source.reveal === "object" ? source.reveal : {};
        const guides = source.guides && typeof source.guides === "object" ? source.guides : {};
        const sourceVersion = Math.max(1, Math.floor(Number(source.schema_version) || 1));
        return {
          schema_version: OVIZ_DECK_VERSION,
          available: source.available === undefined ? true : Boolean(source.available),
          enabled: Boolean(source.enabled || rawSlides.length),
          embedded: Boolean(source.embedded),
          revision: Math.max(0, Math.floor(Number(source.revision) || 0)),
          aspect_ratio: "16:9",
          guides: {
            smart: guides.smart === undefined ? true : Boolean(guides.smart),
            grid: Boolean(guides.grid),
            grid_size: ovizDeckClamp(guides.grid_size, 5, 200, 20),
          },
          reveal: {
            version: String(reveal.version || OVIZ_DECK_REVEAL_VERSION),
            transition: String(reveal.transition || "fade"),
            background_transition: String(reveal.background_transition || "fade"),
          },
          slides: rawSlides.map((rawSlide, slideIndex) => {
            const slide = rawSlide && typeof rawSlide === "object" ? rawSlide : {};
            let id = String(slide.id || "").trim() || ovizDeckUuid("slide");
            while (seenSlides.has(id)) id = ovizDeckUuid("slide");
            seenSlides.add(id);
            const stateId = slide.state_id;
            const usesLegacyBlocks = !Array.isArray(slide.objects);
            const rawObjects = usesLegacyBlocks
              ? (Array.isArray(slide.blocks) ? slide.blocks : [])
              : slide.objects;
            return {
              id,
              name: String(slide.name || `Slide ${slideIndex + 1}`).trim() || `Slide ${slideIndex + 1}`,
              state_id: stateId === null || stateId === undefined || stateId === "" || stateId === "original"
                ? null
                : String(stateId),
              objects: rawObjects.map(
                (block, blockIndex) => ovizDeckNormalizeBlock(
                  block,
                  blockIndex,
                  seenBlocks,
                  usesLegacyBlocks || sourceVersion < 2,
                )
              ),
              notes: String(slide.notes || ""),
            };
          }),
        };
      }

      function ovizDeckExportSpec(options = {}) {
        const output = ovizDeckClone(ovizDeckProject || ovizDeckNormalizeProject(null), {});
        output.enabled = Boolean(output.available && output.slides && output.slides.length);
        output.embedded = options.embedded === undefined ? true : Boolean(options.embedded);
        return output;
      }

      function ovizDeckEvent(name, detail = {}) {
        const payload = Object.assign({
          rootId: root.id,
          revision: ovizDeckProject ? ovizDeckProject.revision : 0,
        }, detail || {});
        root.dispatchEvent(new CustomEvent(name, { detail: payload }));
        try {
          if (window.parent && window.parent !== window) {
            window.parent.postMessage({ source: "oviz", type: name, detail: payload }, "*");
          }
        } catch (_err) {
          // Parent messaging is optional; root events remain authoritative.
        }
      }

      function ovizDeckHistorySnapshot() {
        return {
          project: ovizDeckClone(ovizDeckProject, {}),
          active_index: ovizDeckActiveIndex,
          selected_block_id: ovizDeckSelectedBlockId,
          selected_object_ids: Array.from(ovizDeckSelectedObjectIds),
        };
      }

      function ovizDeckHistoryKey(reason) {
        const continuous = ["edit-text", "resize-text", "color-text", "line-height-text"].includes(reason);
        return continuous ? `${reason}:${ovizDeckSelectedBlockId}` : "";
      }

      function ovizDeckSyncHistoryUi() {
        if (!ovizDeckEditorEl) return;
        const undo = ovizDeckEditorEl.querySelector(".oviz-deck-undo");
        const redo = ovizDeckEditorEl.querySelector(".oviz-deck-redo");
        if (undo) undo.disabled = ovizDeckUndoStack.length === 0;
        if (redo) redo.disabled = ovizDeckRedoStack.length === 0;
        root.dataset.deckCanUndo = ovizDeckUndoStack.length ? "true" : "false";
        root.dataset.deckCanRedo = ovizDeckRedoStack.length ? "true" : "false";
      }

      function ovizDeckRecordHistory(reason) {
        const current = ovizDeckHistorySnapshot();
        if (!ovizDeckHistoryPresent) {
          ovizDeckHistoryPresent = current;
          ovizDeckSyncHistoryUi();
          return;
        }
        if (JSON.stringify(current) === JSON.stringify(ovizDeckHistoryPresent)) return;
        const now = performance.now();
        const historyKey = ovizDeckHistoryKey(reason);
        const coalesce = Boolean(
          historyKey
          && historyKey === ovizDeckLastHistoryKey
          && now - ovizDeckLastHistoryAt < 700
        );
        if (!coalesce) {
          ovizDeckUndoStack.push(ovizDeckHistoryPresent);
          if (ovizDeckUndoStack.length > 100) ovizDeckUndoStack.shift();
        }
        ovizDeckHistoryPresent = current;
        ovizDeckRedoStack = [];
        ovizDeckLastHistoryKey = historyKey;
        ovizDeckLastHistoryAt = now;
        ovizDeckSyncHistoryUi();
      }

      function ovizDeckApplyHistorySnapshot(snapshot, reason) {
        if (!snapshot || !snapshot.project) return false;
        const previousRevision = ovizDeckProject ? Number(ovizDeckProject.revision) || 0 : 0;
        ovizDeckProject = ovizDeckNormalizeProject(snapshot.project);
        ovizDeckProject.revision = Math.max(previousRevision, Number(ovizDeckProject.revision) || 0) + 1;
        ovizDeckActiveIndex = Math.min(
          Math.max(Number(snapshot.active_index) || 0, 0),
          Math.max(ovizDeckProject.slides.length - 1, 0),
        );
        const active = ovizDeckActiveSlide();
        ovizDeckSelectedBlockId = active && active.objects.some((block) => block.id === snapshot.selected_block_id)
          ? snapshot.selected_block_id
          : "";
        const restoredSelection = Array.isArray(snapshot.selected_object_ids)
          ? snapshot.selected_object_ids.filter((id) => active && active.objects.some((object) => object.id === id))
          : (ovizDeckSelectedBlockId ? [ovizDeckSelectedBlockId] : []);
        ovizDeckSelectedObjectIds = new Set(restoredSelection);
        ovizDeckDirty = true;
        sceneSpec.deck = ovizDeckExportSpec({ embedded: false });
        ovizDeckHistoryPresent = ovizDeckHistorySnapshot();
        ovizDeckLastHistoryKey = "";
        ovizDeckLastHistoryAt = 0;
        ovizDeckRenderEditor();
        ovizDeckRenderAuthoringSlide();
        ovizDeckRenderRevealSlides();
        ovizDeckEvent("deck-changed", { reason, slides: ovizDeckList() });
        ovizDeckSyncHistoryUi();
        return true;
      }

      function ovizDeckUndo() {
        if (!ovizDeckUndoStack.length) return false;
        const target = ovizDeckUndoStack.pop();
        if (ovizDeckHistoryPresent) ovizDeckRedoStack.push(ovizDeckHistoryPresent);
        return ovizDeckApplyHistorySnapshot(target, "undo");
      }

      function ovizDeckRedo() {
        if (!ovizDeckRedoStack.length) return false;
        const target = ovizDeckRedoStack.pop();
        if (ovizDeckHistoryPresent) ovizDeckUndoStack.push(ovizDeckHistoryPresent);
        return ovizDeckApplyHistorySnapshot(target, "redo");
      }

      function ovizDeckHandleHistoryShortcut(direction) {
        if (!ovizDeckEditorOpen || ovizDeckPresenting) return false;
        return direction === "redo" ? ovizDeckRedo() : ovizDeckUndo();
      }

      function ovizDeckChanged(reason, options = {}) {
        if (!ovizDeckProject) return;
        ovizDeckProject.enabled = ovizDeckProject.slides.length > 0;
        ovizDeckProject.revision += 1;
        ovizDeckDirty = true;
        sceneSpec.deck = ovizDeckExportSpec({ embedded: false });
        ovizDeckRecordHistory(reason);
        if (options.render !== false) {
          ovizDeckRenderEditor();
          ovizDeckRenderAuthoringSlide();
        } else {
          ovizDeckSyncStatus();
        }
        ovizDeckEvent("deck-changed", {
          reason,
          slides: ovizDeckList(),
        });
        if (/(block|object|text|shape|geometry|align|distribute|arrange|group|lock|font|list)/.test(String(reason))) {
          ovizDeckEvent("deck-object-changed", {
            reason,
            slide_id: ovizDeckActiveSlide() ? ovizDeckActiveSlide().id : null,
            object_ids: Array.from(ovizDeckSelectedObjectIds),
            objects: ovizDeckSelectedObjects().map((object) => ovizDeckClone(object, {})),
          });
        }
      }

      function ovizDeckList() {
        if (!ovizDeckProject) return [];
        return ovizDeckProject.slides.map((slide, index) => ({
          id: slide.id,
          index,
          number: index + 1,
          name: slide.name,
          state_id: slide.state_id,
          block_count: slide.objects.length,
          object_count: slide.objects.length,
          active: index === ovizDeckActiveIndex,
        }));
      }

      function ovizDeckActiveSlide() {
        if (!ovizDeckProject || !ovizDeckProject.slides.length) return null;
        ovizDeckActiveIndex = Math.min(Math.max(ovizDeckActiveIndex, 0), ovizDeckProject.slides.length - 1);
        return ovizDeckProject.slides[ovizDeckActiveIndex];
      }

      function ovizDeckSlideIndex(idOrIndex) {
        if (!ovizDeckProject) return -1;
        if (typeof idOrIndex === "number" && Number.isFinite(idOrIndex)) {
          const index = Math.floor(idOrIndex);
          return index >= 0 && index < ovizDeckProject.slides.length ? index : -1;
        }
        const text = String(idOrIndex || "").trim();
        if (/^\d+$/.test(text)) {
          const index = Number(text) - 1;
          return index >= 0 && index < ovizDeckProject.slides.length ? index : -1;
        }
        return ovizDeckProject.slides.findIndex((slide) => slide.id === text);
      }

      function ovizDeckStateExists(stateId) {
        if (!stateId) return true;
        return Boolean(ovizStatesProject && ovizStatesProject.items.some((state) => state.id === stateId));
      }

      function ovizDeckNavigateToState(slide) {
        if (!slide) return Promise.resolve({ boundary: true });
        const target = slide.state_id || "original";
        if (slide.state_id && !ovizDeckStateExists(slide.state_id)) {
          const error = new Error(`Slide state is missing: ${slide.state_id}`);
          ovizDeckSetStatus(error.message, true);
          return Promise.reject(error);
        }
        return Promise.resolve(ovizGoToState(target)).catch((error) => {
          ovizDeckSetStatus(`State transition failed: ${String(error && error.message || error)}`, true);
          throw error;
        });
      }

      async function ovizDeckGoTo(idOrIndex, options = {}) {
        const index = ovizDeckSlideIndex(idOrIndex);
        if (index < 0) throw new Error("Unknown slide.");
        ovizDeckActiveIndex = index;
        ovizDeckSelectedBlockId = "";
        ovizDeckSelectedObjectIds.clear();
        ovizDeckRenderEditor();
        ovizDeckRenderAuthoringSlide();
        if (ovizDeckPresenting && ovizDeckReveal && options.fromReveal !== true) {
          ovizDeckReveal.slide(index);
        }
        ovizDeckEvent("deck-slide-changed", {
          index,
          slide: ovizDeckClone(ovizDeckProject.slides[index], {}),
        });
        return ovizDeckNavigateToState(ovizDeckProject.slides[index]);
      }

      function ovizDeckEnsureStateEditMode() {
        if (ovizStatesMode !== "edit") ovizSetStatesMode("edit");
      }

      function ovizDeckRequireAvailable() {
        if (!ovizDeckProject || !ovizDeckProject.available) {
          throw new Error("Slides are not available in this figure.");
        }
      }

      function ovizDeckDefaultBlocks(slideNumber) {
        return [
          ovizDeckNormalizeBlock({
            type: "title",
            text: `Slide ${slideNumber}`,
            x: 112,
            y: 558,
            width: 1088,
            height: 120,
            font_size: 64,
            font_weight: 700,
          }, 0, new Set()),
          ovizDeckNormalizeBlock({
            type: "subtitle",
            text: "Double-click to edit this text",
            x: 112,
            y: 684,
            width: 992,
            height: 76,
            font_size: 30,
            font_weight: 400,
            color: "#e6eef8",
          }, 1, new Set()),
        ];
      }

      function ovizDeckAdd(options = {}) {
        ovizDeckRequireAvailable();
        ovizDeckEnsureStateEditMode();
        const number = ovizDeckProject.slides.length + 1;
        const state = options.state_id !== undefined
          ? null
          : ovizAddState({ name: options.state_name || `Slide ${number} view` });
        const slide = {
          id: String(options.id || "").trim() || ovizDeckUuid("slide"),
          name: String(options.name || `Slide ${number}`).trim() || `Slide ${number}`,
          state_id: options.state_id === null || options.state_id === "original"
            ? null
            : String(options.state_id || (state && state.id) || "") || null,
          objects: Array.isArray(options.objects || options.blocks)
            ? (options.objects || options.blocks).map((block, index) => ovizDeckNormalizeBlock(
              block,
              index,
              new Set(),
              !Array.isArray(options.objects),
            ))
            : ovizDeckDefaultBlocks(number),
          notes: String(options.notes || ""),
        };
        ovizDeckProject.slides.push(slide);
        ovizDeckActiveIndex = ovizDeckProject.slides.length - 1;
        ovizDeckSelectedBlockId = slide.objects.length ? slide.objects[0].id : "";
        ovizDeckSelectedObjectIds = new Set(ovizDeckSelectedBlockId ? [ovizDeckSelectedBlockId] : []);
        ovizDeckChanged("add-slide");
        return ovizDeckClone(slide, {});
      }

      function ovizDeckUpdateState(idOrIndex = ovizDeckActiveIndex) {
        ovizDeckEnsureStateEditMode();
        const index = ovizDeckSlideIndex(idOrIndex);
        if (index < 0) throw new Error("Unknown slide.");
        const slide = ovizDeckProject.slides[index];
        if (slide.state_id && ovizDeckStateExists(slide.state_id)) {
          ovizUpdateState(slide.state_id, { name: `${slide.name} view` });
        } else {
          const state = ovizAddState({ name: `${slide.name} view` });
          slide.state_id = state.id;
        }
        ovizDeckChanged("update-slide-state");
        return ovizDeckClone(slide, {});
      }

      function ovizDeckRename(idOrIndex, name) {
        const index = ovizDeckSlideIndex(idOrIndex);
        if (index < 0) throw new Error("Unknown slide.");
        const slide = ovizDeckProject.slides[index];
        slide.name = String(name || "").trim() || slide.name;
        ovizDeckChanged("rename-slide");
        return ovizDeckClone(slide, {});
      }

      function ovizDeckDuplicate(idOrIndex = ovizDeckActiveIndex) {
        const index = ovizDeckSlideIndex(idOrIndex);
        if (index < 0) throw new Error("Unknown slide.");
        const source = ovizDeckProject.slides[index];
        const copy = ovizDeckClone(source, {});
        copy.id = ovizDeckUuid("slide");
        copy.name = `${source.name} copy`;
        const groupMap = new Map();
        copy.objects = copy.objects.map((block) => {
          const groupId = block.group_id
            ? (groupMap.get(block.group_id) || (() => {
              const next = ovizDeckUuid("group");
              groupMap.set(block.group_id, next);
              return next;
            })())
            : null;
          return Object.assign({}, block, { id: ovizDeckUuid("object"), group_id: groupId });
        });
        ovizDeckProject.slides.splice(index + 1, 0, copy);
        ovizDeckActiveIndex = index + 1;
        ovizDeckSelectedBlockId = "";
        ovizDeckSelectedObjectIds.clear();
        ovizDeckChanged("duplicate-slide");
        return ovizDeckClone(copy, {});
      }

      function ovizDeckMove(idOrIndex, destinationIndex) {
        const index = ovizDeckSlideIndex(idOrIndex);
        if (index < 0) throw new Error("Unknown slide.");
        const destination = Math.min(
          Math.max(Math.floor(Number(destinationIndex)), 0),
          Math.max(ovizDeckProject.slides.length - 1, 0),
        );
        const [slide] = ovizDeckProject.slides.splice(index, 1);
        ovizDeckProject.slides.splice(destination, 0, slide);
        ovizDeckActiveIndex = destination;
        ovizDeckChanged("move-slide");
        return ovizDeckList();
      }

      function ovizDeckRemove(idOrIndex = ovizDeckActiveIndex) {
        const index = ovizDeckSlideIndex(idOrIndex);
        if (index < 0) throw new Error("Unknown slide.");
        const [slide] = ovizDeckProject.slides.splice(index, 1);
        ovizDeckActiveIndex = Math.min(index, Math.max(ovizDeckProject.slides.length - 1, 0));
        ovizDeckSelectedBlockId = "";
        ovizDeckSelectedObjectIds.clear();
        ovizDeckChanged("remove-slide");
        return ovizDeckClone(slide, {});
      }

      function ovizDeckAddBlock(type = "text") {
        ovizDeckRequireAvailable();
        const slide = ovizDeckActiveSlide();
        if (!slide) throw new Error("Add a slide before adding text.");
        const seenIds = new Set(slide.objects.map((block) => block.id));
        const block = ovizDeckNormalizeBlock({
          kind: "text",
          type,
          text: type === "title" ? "Title" : (type === "subtitle" ? "Subtitle" : "Text"),
          x: 192,
          y: 162 + slide.objects.length * 72,
          width: 928,
          height: type === "title" ? 120 : 100,
        }, slide.objects.length, seenIds);
        slide.objects.push(block);
        ovizDeckSelectedBlockId = block.id;
        ovizDeckSelectedObjectIds = new Set([block.id]);
        ovizDeckChanged("add-block");
        return ovizDeckClone(block, {});
      }

      function ovizDeckRemoveBlock(blockId = ovizDeckSelectedBlockId) {
        const slide = ovizDeckActiveSlide();
        if (!slide) throw new Error("No active slide.");
        const index = slide.objects.findIndex((block) => block.id === String(blockId));
        if (index < 0) throw new Error("Unknown text block.");
        const [block] = slide.objects.splice(index, 1);
        ovizDeckSelectedBlockId = "";
        ovizDeckSelectedObjectIds.clear();
        ovizDeckChanged("remove-block");
        return ovizDeckClone(block, {});
      }

      function ovizDeckObjectById(objectId, slide = ovizDeckActiveSlide()) {
        if (!slide) return null;
        return slide.objects.find((object) => object.id === String(objectId)) || null;
      }

      function ovizDeckSelectedObjects(options = {}) {
        const slide = ovizDeckActiveSlide();
        if (!slide) return [];
        return slide.objects.filter((object) => (
          ovizDeckSelectedObjectIds.has(object.id)
          && (!options.unlockedOnly || !object.locked)
        ));
      }

      function ovizDeckSetSelection(ids, options = {}) {
        const slide = ovizDeckActiveSlide();
        const valid = new Set();
        (Array.isArray(ids) ? ids : [ids]).forEach((id) => {
          const object = ovizDeckObjectById(id, slide);
          if (!object) return;
          if (options.expandGroup !== false && object.group_id) {
            slide.objects.forEach((candidate) => {
              if (candidate.group_id === object.group_id) valid.add(candidate.id);
            });
          } else {
            valid.add(object.id);
          }
        });
        ovizDeckSelectedObjectIds = valid;
        const requestedPrimary = String(options.primary || "");
        ovizDeckSelectedBlockId = valid.has(requestedPrimary)
          ? requestedPrimary
          : (valid.values().next().value || "");
        if (ovizDeckHistoryPresent) {
          ovizDeckHistoryPresent.active_index = ovizDeckActiveIndex;
          ovizDeckHistoryPresent.selected_block_id = ovizDeckSelectedBlockId;
          ovizDeckHistoryPresent.selected_object_ids = Array.from(valid);
        }
        if (options.render !== false) {
          ovizDeckRenderEditor();
          ovizDeckRenderAuthoringSlide();
        }
        ovizDeckEvent("deck-selection-changed", {
          slide_id: slide ? slide.id : null,
          object_ids: Array.from(valid),
          primary_id: ovizDeckSelectedBlockId || null,
        });
        return Array.from(valid);
      }

      function ovizDeckSelectObjectFromPointer(object, event, options = {}) {
        if (!object) return [];
        const toggle = Boolean(event && (event.shiftKey || event.metaKey || event.ctrlKey));
        if (!toggle) return ovizDeckSetSelection([object.id], { primary: object.id, render: options.render, expandGroup: true });
        const members = object.group_id && ovizDeckActiveSlide()
          ? ovizDeckActiveSlide().objects.filter((candidate) => candidate.group_id === object.group_id).map((candidate) => candidate.id)
          : [object.id];
        const ids = new Set(ovizDeckSelectedObjectIds);
        const remove = members.every((id) => ids.has(id));
        members.forEach((id) => { if (remove) ids.delete(id); else ids.add(id); });
        return ovizDeckSetSelection(Array.from(ids), { primary: object.id, render: options.render, expandGroup: false });
      }

      function ovizDeckObjectList() {
        const slide = ovizDeckActiveSlide();
        return slide ? ovizDeckClone(slide.objects, []) : [];
      }

      function ovizDeckObjectGet(objectId) {
        return ovizDeckClone(ovizDeckObjectById(objectId), null);
      }

      function ovizDeckObjectUpdate(objectId, patch = {}) {
        const object = ovizDeckObjectById(objectId);
        if (!object) throw new Error("Unknown slide object.");
        const index = ovizDeckActiveSlide().objects.indexOf(object);
        const normalized = ovizDeckNormalizeBlock(Object.assign({}, object, patch, { id: object.id }), index, new Set());
        Object.assign(object, normalized, { id: object.id });
        ovizDeckSelectedBlockId = object.id;
        ovizDeckSelectedObjectIds.add(object.id);
        ovizDeckChanged("update-object");
        return ovizDeckClone(object, {});
      }

      function ovizDeckAddShape(shapeType = "rectangle", options = {}) {
        ovizDeckRequireAvailable();
        const slide = ovizDeckActiveSlide();
        if (!slide) throw new Error("Add a slide before adding a shape.");
        const object = ovizDeckNormalizeBlock(Object.assign({
          kind: "shape",
          shape_type: shapeType,
          name: String(shapeType).replace(/_/g, " ").replace(/\b\w/g, (letter) => letter.toUpperCase()),
          x: 660,
          y: 310,
          width: ["line", "arrow"].includes(shapeType) ? 360 : 280,
          height: ["line", "arrow"].includes(shapeType) ? 60 : 180,
          text: "",
        }, options), slide.objects.length, new Set(slide.objects.map((item) => item.id)));
        slide.objects.push(object);
        ovizDeckSetSelection([object.id], { primary: object.id, render: false });
        ovizDeckAddMenuOpen = false;
        ovizDeckChanged("add-shape");
        return ovizDeckClone(object, {});
      }

      function ovizDeckRemoveObjects(ids = Array.from(ovizDeckSelectedObjectIds)) {
        const slide = ovizDeckActiveSlide();
        if (!slide) throw new Error("No active slide.");
        const wanted = new Set(Array.isArray(ids) ? ids.map(String) : [String(ids)]);
        const removed = slide.objects.filter((object) => wanted.has(object.id) && !object.locked);
        if (!removed.length) return [];
        slide.objects = slide.objects.filter((object) => !removed.includes(object));
        ovizDeckSelectedObjectIds.clear();
        ovizDeckSelectedBlockId = "";
        ovizDeckChanged("remove-object");
        return ovizDeckClone(removed, []);
      }

      function ovizDeckDuplicateObjects(ids = Array.from(ovizDeckSelectedObjectIds), offset = { x: 20, y: 20 }) {
        const slide = ovizDeckActiveSlide();
        if (!slide) throw new Error("No active slide.");
        const wanted = new Set(Array.isArray(ids) ? ids.map(String) : [String(ids)]);
        const source = slide.objects.filter((object) => wanted.has(object.id));
        if (!source.length) return [];
        const groupMap = new Map();
        const copies = source.map((object) => {
          const copy = ovizDeckClone(object, {});
          copy.id = ovizDeckUuid("object");
          copy.name = `${object.name || object.type || "Object"} copy`;
          copy.x = ovizDeckClamp(object.x + ovizDeckFinite(offset.x, 20), 0, OVIZ_DECK_WIDTH - object.width, object.x);
          copy.y = ovizDeckClamp(object.y + ovizDeckFinite(offset.y, 20), 0, OVIZ_DECK_HEIGHT - object.height, object.y);
          copy.group_id = object.group_id
            ? (groupMap.get(object.group_id) || (() => {
              const groupId = ovizDeckUuid("group");
              groupMap.set(object.group_id, groupId);
              return groupId;
            })())
            : null;
          return copy;
        });
        slide.objects.push(...copies);
        ovizDeckSetSelection(copies.map((object) => object.id), { primary: copies[0].id, render: false, expandGroup: false });
        ovizDeckChanged("duplicate-object");
        return ovizDeckClone(copies, []);
      }

      function ovizDeckGroupObjects(ids = Array.from(ovizDeckSelectedObjectIds)) {
        const slide = ovizDeckActiveSlide();
        const wanted = new Set(Array.isArray(ids) ? ids.map(String) : [String(ids)]);
        const objects = slide ? slide.objects.filter((object) => wanted.has(object.id)) : [];
        if (objects.length < 2) return null;
        const groupId = ovizDeckUuid("group");
        objects.forEach((object) => { object.group_id = groupId; });
        ovizDeckChanged("group-objects");
        return groupId;
      }

      function ovizDeckUngroupObjects(ids = Array.from(ovizDeckSelectedObjectIds)) {
        const slide = ovizDeckActiveSlide();
        if (!slide) return [];
        const selectedGroups = new Set();
        (Array.isArray(ids) ? ids : [ids]).forEach((id) => {
          const object = ovizDeckObjectById(id, slide);
          if (object && object.group_id) selectedGroups.add(object.group_id);
        });
        const changed = slide.objects.filter((object) => selectedGroups.has(object.group_id));
        changed.forEach((object) => { object.group_id = null; });
        if (changed.length) ovizDeckChanged("ungroup-objects");
        return changed.map((object) => object.id);
      }

      function ovizDeckLockObjects(ids = Array.from(ovizDeckSelectedObjectIds), locked = true) {
        const wanted = new Set(Array.isArray(ids) ? ids.map(String) : [String(ids)]);
        const changed = (ovizDeckActiveSlide() ? ovizDeckActiveSlide().objects : [])
          .filter((object) => wanted.has(object.id));
        changed.forEach((object) => { object.locked = Boolean(locked); });
        if (changed.length) ovizDeckChanged(locked ? "lock-objects" : "unlock-objects");
        return changed.map((object) => object.id);
      }

      function ovizDeckArrangeObjects(action, ids = Array.from(ovizDeckSelectedObjectIds)) {
        const slide = ovizDeckActiveSlide();
        if (!slide) return [];
        const wanted = new Set(Array.isArray(ids) ? ids.map(String) : [String(ids)]);
        const selected = slide.objects.filter((object) => wanted.has(object.id) && !object.locked);
        if (!selected.length) return [];
        if (action === "front") {
          slide.objects = slide.objects.filter((object) => !wanted.has(object.id)).concat(selected);
        } else if (action === "back") {
          slide.objects = selected.concat(slide.objects.filter((object) => !wanted.has(object.id)));
        } else {
          const direction = action === "forward" ? 1 : (action === "backward" ? -1 : 0);
          if (!direction) throw new Error("Unknown arrange action.");
          const ordered = direction > 0 ? [...selected].reverse() : selected;
          ordered.forEach((object) => {
            const index = slide.objects.indexOf(object);
            const destination = Math.min(Math.max(index + direction, 0), slide.objects.length - 1);
            if (destination === index) return;
            slide.objects.splice(index, 1);
            slide.objects.splice(destination, 0, object);
          });
        }
        ovizDeckChanged(`arrange-${action}`);
        return ovizDeckObjectList();
      }

      function ovizDeckAlignObjects(alignment, ids = Array.from(ovizDeckSelectedObjectIds)) {
        const wanted = new Set(Array.isArray(ids) ? ids.map(String) : [String(ids)]);
        const objects = (ovizDeckActiveSlide() ? ovizDeckActiveSlide().objects : [])
          .filter((object) => wanted.has(object.id) && !object.locked);
        if (!objects.length) return [];
        const left = Math.min(...objects.map((object) => object.x));
        const top = Math.min(...objects.map((object) => object.y));
        const right = Math.max(...objects.map((object) => object.x + object.width));
        const bottom = Math.max(...objects.map((object) => object.y + object.height));
        const centerX = (left + right) / 2;
        const centerY = (top + bottom) / 2;
        objects.forEach((object) => {
          if (alignment === "left") object.x = left;
          else if (alignment === "center") object.x = centerX - object.width / 2;
          else if (alignment === "right") object.x = right - object.width;
          else if (alignment === "top") object.y = top;
          else if (alignment === "middle") object.y = centerY - object.height / 2;
          else if (alignment === "bottom") object.y = bottom - object.height;
          else throw new Error("Unknown alignment.");
        });
        ovizDeckChanged(`align-${alignment}`);
        return objects.map((object) => object.id);
      }

      function ovizDeckDistributeObjects(axis, ids = Array.from(ovizDeckSelectedObjectIds)) {
        const wanted = new Set(Array.isArray(ids) ? ids.map(String) : [String(ids)]);
        const objects = (ovizDeckActiveSlide() ? ovizDeckActiveSlide().objects : [])
          .filter((object) => wanted.has(object.id) && !object.locked)
          .sort((a, b) => axis === "horizontal" ? a.x - b.x : a.y - b.y);
        if (objects.length < 3) return [];
        const start = axis === "horizontal" ? objects[0].x : objects[0].y;
        const endObject = objects[objects.length - 1];
        const end = axis === "horizontal" ? endObject.x + endObject.width : endObject.y + endObject.height;
        const totalSize = objects.reduce((sum, object) => sum + (axis === "horizontal" ? object.width : object.height), 0);
        const gap = (end - start - totalSize) / (objects.length - 1);
        let cursor = start;
        objects.forEach((object) => {
          if (axis === "horizontal") { object.x = cursor; cursor += object.width + gap; }
          else { object.y = cursor; cursor += object.height + gap; }
        });
        ovizDeckChanged(`distribute-${axis}`);
        return objects.map((object) => object.id);
      }

      function ovizDeckGuidesGet() {
        return ovizDeckClone(ovizDeckProject && ovizDeckProject.guides, { smart: true, grid: false, grid_size: 20 });
      }

      function ovizDeckGuidesSet(patch = {}) {
        ovizDeckProject.guides = Object.assign({}, ovizDeckProject.guides, patch);
        ovizDeckProject.guides.smart = Boolean(ovizDeckProject.guides.smart);
        ovizDeckProject.guides.grid = Boolean(ovizDeckProject.guides.grid);
        ovizDeckProject.guides.grid_size = ovizDeckClamp(ovizDeckProject.guides.grid_size, 5, 200, 20);
        ovizDeckChanged("change-guides");
        return ovizDeckGuidesGet();
      }

      function ovizDeckSetStatus(text, isError = false) {
        if (!ovizDeckStatusEl) return;
        ovizDeckStatusEl.textContent = String(text || "");
        ovizDeckStatusEl.dataset.error = isError ? "true" : "false";
      }

      function ovizDeckSyncStatus() {
        if (!ovizDeckStatusEl || !ovizDeckProject) return;
        const slide = ovizDeckActiveSlide();
        const missing = slide && slide.state_id && !ovizDeckStateExists(slide.state_id);
        if (missing) {
          ovizDeckSetStatus("This slide's State was removed. Capture or link a new view.", true);
        } else {
          ovizDeckSetStatus(ovizDeckDirty ? "Unexported changes" : "", false);
        }
      }

      function ovizDeckMakeButton(label, callback, className = "") {
        const button = document.createElement("button");
        button.type = "button";
        button.textContent = label;
        if (className) button.className = className;
        button.addEventListener("click", (event) => {
          event.preventDefault();
          event.stopPropagation();
          Promise.resolve().then(callback).catch((error) => {
            ovizDeckSetStatus(String(error && error.message || error), true);
          });
        });
        return button;
      }

      function ovizDeckRenderEditor() {
        if (!ovizDeckEditorEl || !ovizDeckProject) return;
        ovizDeckEditorEl.dataset.open = ovizDeckEditorOpen ? "true" : "false";
        ovizDeckEditorEl.setAttribute("aria-hidden", ovizDeckEditorOpen ? "false" : "true");
        if (!ovizDeckEditorOpen) return;
        ovizDeckEditorEl.innerHTML = "";

        const head = document.createElement("div");
        head.className = "oviz-deck-editor-head";
        const title = document.createElement("strong");
        title.textContent = "Slides";
        const count = document.createElement("small");
        count.textContent = `${ovizDeckProject.slides.length} slide${ovizDeckProject.slides.length === 1 ? "" : "s"}`;
        const undo = ovizDeckMakeButton("↶", () => ovizDeckUndo(), "oviz-deck-history-button oviz-deck-undo");
        undo.title = "Undo slide change (⌘Z)";
        undo.setAttribute("aria-label", "Undo slide change");
        const redo = ovizDeckMakeButton("↷", () => ovizDeckRedo(), "oviz-deck-history-button oviz-deck-redo");
        redo.title = "Redo slide change (⌘⇧Z)";
        redo.setAttribute("aria-label", "Redo slide change");
        const close = ovizDeckMakeButton("Close", () => ovizDeckSetEditorOpen(false));
        head.append(title, count, undo, redo, close);
        ovizDeckEditorEl.append(head);
        ovizDeckSyncHistoryUi();

        const toolbar = document.createElement("div");
        toolbar.className = "oviz-deck-toolbar";
        toolbar.append(
          ovizDeckMakeButton("+ Slide", () => ovizDeckAdd()),
          ovizDeckMakeButton("+ Object", () => {
            ovizDeckAddMenuOpen = !ovizDeckAddMenuOpen;
            ovizDeckRenderEditor();
          }),
          ovizDeckMakeButton("Present", () => setPresentationMode(true)),
          ovizDeckMakeButton("Export", () => ovizPromptExportStatesHtml()),
        );
        ovizDeckEditorEl.append(toolbar);
        if (ovizDeckAddMenuOpen) {
          const addMenu = document.createElement("div");
          addMenu.className = "oviz-deck-add-menu";
          const addTextButton = ovizDeckMakeButton("Text", () => {
            ovizDeckAddMenuOpen = false;
            return ovizDeckAddBlock("text");
          });
          addMenu.append(addTextButton);
          [
            ["rectangle", "Rectangle"],
            ["rounded_rectangle", "Rounded rectangle"],
            ["ellipse", "Ellipse"],
            ["triangle", "Triangle"],
            ["diamond", "Diamond"],
            ["line", "Line"],
            ["arrow", "Arrow"],
          ].forEach(([value, label]) => addMenu.append(ovizDeckMakeButton(label, () => ovizDeckAddShape(value))));
          ovizDeckEditorEl.append(addMenu);
        }
        const hint = document.createElement("div");
        hint.className = "oviz-deck-hint";
        hint.textContent = "Drag to move; Shift constrains; ⌘ disables snapping. Double-click shapes to type.";
        ovizDeckEditorEl.append(hint);

        const rows = document.createElement("div");
        rows.className = "oviz-deck-slide-list";
        ovizDeckProject.slides.forEach((slide, index) => {
          const row = document.createElement("div");
          row.className = "oviz-deck-slide-row";
          row.dataset.active = index === ovizDeckActiveIndex ? "true" : "false";
          const number = document.createElement("span");
          number.className = "oviz-deck-slide-number";
          number.textContent = String(index + 1);
          const name = ovizDeckMakeButton(slide.name, () => ovizDeckGoTo(index), "oviz-deck-slide-name");
          name.title = slide.name;
          const actions = document.createElement("details");
          actions.className = "oviz-deck-slide-actions";
          const actionsToggle = document.createElement("summary");
          actionsToggle.textContent = "•••";
          actionsToggle.title = `Options for ${slide.name}`;
          actionsToggle.setAttribute("aria-label", `Options for ${slide.name}`);
          actionsToggle.addEventListener("click", (event) => event.stopPropagation());
          const actionsMenu = document.createElement("div");
          actionsMenu.className = "oviz-deck-slide-action-menu";
          const menuAction = (label, callback, danger = false) => {
            const button = ovizDeckMakeButton(label, () => {
              actions.open = false;
              return callback();
            });
            if (danger) button.dataset.danger = "true";
            return button;
          };
          actionsMenu.append(
            menuAction("Move up", () => ovizDeckMove(index, Math.max(index - 1, 0))),
            menuAction("Move down", () => ovizDeckMove(index, Math.min(index + 1, ovizDeckProject.slides.length - 1))),
            menuAction("Rename", () => {
              const next = window.prompt("Slide name", slide.name);
              if (next !== null) ovizDeckRename(index, next);
            }),
            menuAction("Duplicate", () => ovizDeckDuplicate(index)),
            menuAction("Delete", () => ovizDeckRemove(index), true),
          );
          actions.append(actionsToggle, actionsMenu);
          row.append(number, name, actions);
          rows.append(row);
        });
        ovizDeckEditorEl.append(rows);

        const activeSlide = ovizDeckActiveSlide();
        if (activeSlide) {
          const sceneField = document.createElement("label");
          sceneField.className = "oviz-deck-field";
          const sceneLabel = document.createElement("span");
          sceneLabel.textContent = "Linked view";
          const sceneSelect = document.createElement("select");
          const originalOption = document.createElement("option");
          originalOption.value = "";
          originalOption.textContent = "Original scene";
          sceneSelect.append(originalOption);
          (ovizStatesProject ? ovizStatesProject.items : []).forEach((state, index) => {
            const option = document.createElement("option");
            option.value = state.id;
            option.textContent = `${index + 1}. ${state.name}`;
            sceneSelect.append(option);
          });
          sceneSelect.value = activeSlide.state_id || "";
          sceneSelect.addEventListener("change", () => {
            activeSlide.state_id = sceneSelect.value || null;
            ovizDeckChanged("link-state");
            ovizDeckNavigateToState(activeSlide).catch(() => {});
          });
          sceneField.append(sceneLabel, sceneSelect);
          const updateState = ovizDeckMakeButton("Update linked view", () => ovizDeckUpdateState(ovizDeckActiveIndex));
          updateState.className = "oviz-deck-capture-view";
          ovizDeckEditorEl.append(sceneField, updateState);

          const selectedBlock = activeSlide.objects.find((block) => block.id === ovizDeckSelectedBlockId);
          if (selectedBlock && ovizDeckSelectedObjectIds.size === 1) ovizDeckRenderBlockInspector(selectedBlock);
          if (ovizDeckSelectedObjectIds.size > 1) ovizDeckRenderMultiInspector();
          ovizDeckRenderObjectsList(activeSlide);
        } else {
          const empty = document.createElement("div");
          empty.className = "oviz-deck-empty";
          empty.textContent = "Add a slide to capture the current view and place text over the figure.";
          ovizDeckEditorEl.append(empty);
        }

        ovizDeckStatusEl = document.createElement("div");
        ovizDeckStatusEl.className = "oviz-deck-status";
        ovizDeckEditorEl.append(ovizDeckStatusEl);
        ovizDeckSyncStatus();
      }

      function ovizDeckRenderBlockInspector(block) {
        const inspector = document.createElement("div");
        inspector.className = "oviz-deck-block-inspector";
        const heading = document.createElement("strong");
        heading.textContent = block.kind === "shape" ? "Shape format" : "Text format";

        const nameField = document.createElement("label");
        nameField.className = "oviz-deck-field";
        const nameLabel = document.createElement("span");
        nameLabel.textContent = "Object name";
        const nameInput = document.createElement("input");
        nameInput.type = "text";
        nameInput.maxLength = 80;
        nameInput.value = block.name || "";
        nameInput.addEventListener("change", () => {
          block.name = String(nameInput.value || "").trim().slice(0, 80) || block.name;
          ovizDeckChanged("rename-object");
        });
        nameField.append(nameLabel, nameInput);

        const typeField = document.createElement("label");
        typeField.className = "oviz-deck-field";
        const typeLabel = document.createElement("span");
        typeLabel.textContent = "Style";
        const typeSelect = document.createElement("select");
        [["title", "Title"], ["subtitle", "Subtitle"], ["text", "Body"]].forEach(([value, label]) => {
          const option = document.createElement("option");
          option.value = value;
          option.textContent = label;
          typeSelect.append(option);
        });
        typeSelect.value = block.type;
        typeSelect.addEventListener("change", () => {
          block.type = typeSelect.value;
          if (block.type === "title") { block.font_size = 64; block.font_weight = 700; }
          else if (block.type === "subtitle") { block.font_size = 34; block.font_weight = 500; }
          else { block.font_size = 26; block.font_weight = 400; }
          ovizDeckChanged("change-block-style");
        });
        typeField.append(typeLabel, typeSelect);

        const fontField = document.createElement("label");
        fontField.className = "oviz-deck-field";
        const fontLabel = document.createElement("span");
        fontLabel.textContent = "Font";
        const fontSelect = document.createElement("select");
        const fontCatalog = [
          ["sans", "System Sans"], ["serif", "System Serif"], ["mono", "System Mono"],
          ["SF Pro Display", "SF Pro"], ["Avenir Next", "Avenir Next"],
          ["Helvetica Neue", "Helvetica Neue"], ["Futura", "Futura"],
          ["Gill Sans", "Gill Sans"], ["Baskerville", "Baskerville"],
          ["Georgia", "Georgia"], ["Palatino", "Palatino"],
          ["Times New Roman", "Times New Roman"], ["Menlo", "Menlo"],
          ["Monaco", "Monaco"], ["Courier New", "Courier New"],
        ];
        if (!fontCatalog.some(([value]) => value === block.font_family)) {
          fontCatalog.push([block.font_family, `${block.font_family} (fallback if unavailable)`]);
        }
        fontCatalog.forEach(([value, label]) => {
          const option = document.createElement("option");
          option.value = value;
          const available = ["sans", "serif", "mono"].includes(value)
            || !document.fonts
            || document.fonts.check(`12px "${String(value).replace(/["\\]/g, "")}"`);
          option.textContent = available ? label : `${label} (unavailable)`;
          option.dataset.available = available ? "true" : "false";
          fontSelect.append(option);
        });
        const customFontOption = document.createElement("option");
        customFontOption.value = "__custom__";
        customFontOption.textContent = "Custom font…";
        fontSelect.append(customFontOption);
        fontSelect.value = block.font_family;
        fontSelect.addEventListener("change", () => {
          if (fontSelect.value === "__custom__") {
            const candidate = window.prompt("Font family installed on this device", block.font_family === "sans" ? "" : block.font_family);
            if (candidate === null) {
              fontSelect.value = block.font_family;
              return;
            }
            const cleaned = String(candidate).trim().replace(/[;{}]/g, "").slice(0, 100);
            if (!cleaned) {
              fontSelect.value = block.font_family;
              return;
            }
            block.font_family = cleaned;
            const available = !document.fonts || document.fonts.check(`12px "${cleaned.replace(/["\\]/g, "")}"`);
            if (!available) ovizDeckSetStatus(`${cleaned} is not installed here; a system fallback will be used.`, false);
          } else {
            block.font_family = fontSelect.value;
          }
          ovizDeckChanged("font-family-text");
        });
        fontField.append(fontLabel, fontSelect);

        const sizeField = document.createElement("label");
        sizeField.className = "oviz-deck-field";
        const sizeLabel = document.createElement("span");
        sizeLabel.textContent = `Size ${Math.round(block.font_size)} px`;
        const sizeInput = document.createElement("input");
        sizeInput.type = "range";
        sizeInput.min = "10";
        sizeInput.max = "160";
        sizeInput.step = "1";
        sizeInput.value = String(block.font_size);
        sizeInput.addEventListener("input", () => {
          block.font_size = Number(sizeInput.value);
          sizeLabel.textContent = `Size ${Math.round(block.font_size)} px`;
          ovizDeckChanged("resize-text", { render: false });
          ovizDeckRenderAuthoringSlide();
        });
        sizeField.append(sizeLabel, sizeInput);

        const lineHeightField = document.createElement("label");
        lineHeightField.className = "oviz-deck-field";
        const lineHeightLabel = document.createElement("span");
        lineHeightLabel.textContent = `Line spacing ${Number(block.line_height).toFixed(2)}`;
        const lineHeightInput = document.createElement("input");
        lineHeightInput.type = "range";
        lineHeightInput.min = "0.8";
        lineHeightInput.max = "2.4";
        lineHeightInput.step = "0.05";
        lineHeightInput.value = String(block.line_height);
        lineHeightInput.addEventListener("input", () => {
          block.line_height = Number(lineHeightInput.value);
          lineHeightLabel.textContent = `Line spacing ${block.line_height.toFixed(2)}`;
          ovizDeckChanged("line-height-text", { render: false });
          ovizDeckRenderAuthoringSlide();
        });
        lineHeightField.append(lineHeightLabel, lineHeightInput);

        const styleRow = document.createElement("div");
        styleRow.className = "oviz-deck-inspector-row";
        const alignSelect = document.createElement("select");
        ["left", "center", "right"].forEach((value) => {
          const option = document.createElement("option");
          option.value = value;
          option.textContent = value[0].toUpperCase() + value.slice(1);
          alignSelect.append(option);
        });
        alignSelect.value = block.align;
        alignSelect.addEventListener("change", () => {
          block.align = alignSelect.value;
          ovizDeckChanged("align-text");
        });
        const colorInput = document.createElement("input");
        colorInput.type = "color";
        colorInput.value = block.color;
        colorInput.title = "Text color";
        colorInput.addEventListener("input", () => {
          block.color = colorInput.value;
          ovizDeckChanged("color-text", { render: false });
          ovizDeckRenderAuthoringSlide();
        });
        const styleToggle = (label, title, active, callback) => {
          const button = ovizDeckMakeButton(label, callback, "oviz-deck-text-toggle");
          button.title = title;
          button.setAttribute("aria-label", title);
          button.setAttribute("aria-pressed", active ? "true" : "false");
          button.dataset.active = active ? "true" : "false";
          return button;
        };
        const weightButton = styleToggle("B", "Bold", block.font_weight >= 600, () => {
          block.font_weight = block.font_weight >= 600 ? 400 : 700;
          ovizDeckChanged("weight-text");
        });
        const italicButton = styleToggle("I", "Italic", block.font_style === "italic", () => {
          block.font_style = block.font_style === "italic" ? "normal" : "italic";
          ovizDeckChanged("italic-text");
        });
        const underlineButton = styleToggle("U", "Underline", Boolean(block.underline), () => {
          block.underline = !block.underline;
          ovizDeckChanged("underline-text");
        });
        const strikeButton = styleToggle("S", "Strikethrough", Boolean(block.strikethrough), () => {
          block.strikethrough = !block.strikethrough;
          ovizDeckChanged("strikethrough-text");
        });
        styleRow.append(alignSelect, colorInput, weightButton, italicButton, underlineButton, strikeButton);

        const listField = document.createElement("label");
        listField.className = "oviz-deck-field";
        const listLabel = document.createElement("span");
        listLabel.textContent = "List";
        const listSelect = document.createElement("select");
        [["none", "None"], ["bullet", "Bullets"], ["number", "Numbered"]].forEach(([value, label]) => {
          const option = document.createElement("option");
          option.value = value;
          option.textContent = label;
          listSelect.append(option);
        });
        listSelect.value = block.list_style;
        listSelect.addEventListener("change", () => {
          block.list_style = listSelect.value;
          ovizDeckChanged("list-style-text");
        });
        listField.append(listLabel, listSelect);

        const typographyAdvanced = document.createElement("details");
        typographyAdvanced.className = "oviz-deck-inspector-section";
        const typographySummary = document.createElement("summary");
        typographySummary.textContent = "Advanced text";
        const spacingField = ovizDeckNumericField("Character spacing", block.character_spacing, -5, 30, 0.25, (value) => {
          block.character_spacing = value;
          ovizDeckChanged("character-spacing-text");
        });
        const textFxRow = document.createElement("div");
        textFxRow.className = "oviz-deck-inspector-row";
        textFxRow.append(
          styleToggle("Outline", "Text outline", Boolean(block.text_outline), () => {
            block.text_outline = !block.text_outline;
            ovizDeckChanged("outline-text");
          }),
          styleToggle("Shadow", "Text shadow", Boolean(block.text_shadow), () => {
            block.text_shadow = !block.text_shadow;
            ovizDeckChanged("shadow-text");
          }),
        );
        typographyAdvanced.append(typographySummary, spacingField, textFxRow);

        const geometry = ovizDeckRenderGeometryControls(block);
        const shapeControls = block.kind === "shape" ? ovizDeckRenderShapeControls(block) : null;
        const arrange = ovizDeckRenderArrangeControls();

        const remove = ovizDeckMakeButton(`Delete ${block.kind === "shape" ? "shape" : "text box"}`, () => ovizDeckRemoveObjects([block.id]));
        remove.className = "oviz-deck-delete-text";
        remove.title = "Delete selected object (Backspace)";

        inspector.append(heading, nameField);
        if (shapeControls) inspector.append(shapeControls);
        inspector.append(typeField, fontField, sizeField, lineHeightField, styleRow, listField, typographyAdvanced, geometry, arrange, remove);
        ovizDeckEditorEl.append(inspector);
      }

      function ovizDeckNumericField(label, value, minimum, maximum, step, callback) {
        const field = document.createElement("label");
        field.className = "oviz-deck-field oviz-deck-number-field";
        const title = document.createElement("span");
        title.textContent = label;
        const input = document.createElement("input");
        input.type = "number";
        input.min = String(minimum);
        input.max = String(maximum);
        input.step = String(step);
        input.value = String(Math.round(Number(value) * 100) / 100);
        input.addEventListener("change", () => callback(ovizDeckClamp(input.value, minimum, maximum, value)));
        field.append(title, input);
        return field;
      }

      function ovizDeckRenderGeometryControls(block) {
        const details = document.createElement("details");
        details.className = "oviz-deck-inspector-section";
        const summary = document.createElement("summary");
        summary.textContent = "Position & size";
        const grid = document.createElement("div");
        grid.className = "oviz-deck-geometry-grid";
        [
          ["X", "x", 0, OVIZ_DECK_WIDTH], ["Y", "y", 0, OVIZ_DECK_HEIGHT],
          ["W", "width", 20, OVIZ_DECK_WIDTH], ["H", "height", 20, OVIZ_DECK_HEIGHT],
          ["Rotate", "rotation", -360, 360], ["Opacity", "opacity", 0, 1, 0.05],
        ].forEach(([label, key, minimum, maximum, step = 1]) => {
          grid.append(ovizDeckNumericField(label, block[key], minimum, maximum, step, (value) => {
            block[key] = value;
            ovizDeckChanged(`geometry-${key}`);
          }));
        });
        details.append(summary, grid);
        return details;
      }

      function ovizDeckRenderShapeControls(block) {
        const details = document.createElement("details");
        details.className = "oviz-deck-inspector-section";
        details.open = true;
        const summary = document.createElement("summary");
        summary.textContent = "Fill & border";
        const colors = document.createElement("div");
        colors.className = "oviz-deck-inspector-row";
        const colorControl = (label, key) => {
          const field = document.createElement("label");
          field.className = "oviz-deck-color-field";
          const title = document.createElement("span");
          title.textContent = label;
          const input = document.createElement("input");
          input.type = "color";
          input.value = block[key];
          input.addEventListener("input", () => {
            block[key] = input.value;
            ovizDeckChanged(`shape-${key}`, { render: false });
            ovizDeckRenderAuthoringSlide();
          });
          field.append(title, input);
          return field;
        };
        colors.append(colorControl("Fill", "fill_color"), colorControl("Border", "border_color"));
        const fillOpacity = ovizDeckNumericField("Fill opacity", block.fill_opacity, 0, 1, 0.05, (value) => {
          block.fill_opacity = value;
          ovizDeckChanged("shape-fill-opacity");
        });
        const borderWidth = ovizDeckNumericField("Border width", block.border_width, 0, 24, 1, (value) => {
          block.border_width = value;
          ovizDeckChanged("shape-border-width");
        });
        const styleField = document.createElement("label");
        styleField.className = "oviz-deck-field";
        const styleLabel = document.createElement("span");
        styleLabel.textContent = "Border style";
        const styleSelect = document.createElement("select");
        ["solid", "dashed", "dotted"].forEach((value) => {
          const option = document.createElement("option");
          option.value = value;
          option.textContent = value[0].toUpperCase() + value.slice(1);
          styleSelect.append(option);
        });
        styleSelect.value = block.border_style;
        styleSelect.addEventListener("change", () => {
          block.border_style = styleSelect.value;
          ovizDeckChanged("shape-border-style");
        });
        styleField.append(styleLabel, styleSelect);
        const shapeFx = document.createElement("div");
        shapeFx.className = "oviz-deck-inspector-row";
        const shadow = ovizDeckMakeButton("Shadow", () => {
          block.shadow = !block.shadow;
          ovizDeckChanged("shape-shadow");
        }, "oviz-deck-text-toggle");
        shadow.dataset.active = block.shadow ? "true" : "false";
        shapeFx.append(shadow);
        if (block.shape_type === "rounded_rectangle") {
          shapeFx.append(ovizDeckNumericField("Radius", block.corner_radius, 0, 200, 1, (value) => {
            block.corner_radius = value;
            ovizDeckChanged("shape-corner-radius");
          }));
        }
        details.append(summary, colors, fillOpacity, borderWidth, styleField, shapeFx);
        return details;
      }

      function ovizDeckRenderArrangeControls() {
        const details = document.createElement("details");
        details.className = "oviz-deck-inspector-section";
        const summary = document.createElement("summary");
        summary.textContent = "Arrange";
        const buttons = document.createElement("div");
        buttons.className = "oviz-deck-arrange-grid";
        [
          ["Front", () => ovizDeckArrangeObjects("front")],
          ["Back", () => ovizDeckArrangeObjects("back")],
          ["Forward", () => ovizDeckArrangeObjects("forward")],
          ["Backward", () => ovizDeckArrangeObjects("backward")],
          ["Align left", () => ovizDeckAlignObjects("left")],
          ["Align center", () => ovizDeckAlignObjects("center")],
          ["Align right", () => ovizDeckAlignObjects("right")],
          ["Align top", () => ovizDeckAlignObjects("top")],
          ["Align middle", () => ovizDeckAlignObjects("middle")],
          ["Align bottom", () => ovizDeckAlignObjects("bottom")],
          ["Distribute H", () => ovizDeckDistributeObjects("horizontal")],
          ["Distribute V", () => ovizDeckDistributeObjects("vertical")],
        ].forEach(([label, callback]) => buttons.append(ovizDeckMakeButton(label, callback)));
        const groupRow = document.createElement("div");
        groupRow.className = "oviz-deck-inspector-row";
        groupRow.append(
          ovizDeckMakeButton("Group", () => ovizDeckGroupObjects()),
          ovizDeckMakeButton("Ungroup", () => ovizDeckUngroupObjects()),
          ovizDeckMakeButton("Lock", () => ovizDeckLockObjects()),
          ovizDeckMakeButton("Unlock", () => ovizDeckLockObjects(Array.from(ovizDeckSelectedObjectIds), false)),
        );
        details.append(summary, buttons, groupRow);
        return details;
      }

      function ovizDeckRenderMultiInspector() {
        const inspector = document.createElement("div");
        inspector.className = "oviz-deck-block-inspector";
        const heading = document.createElement("strong");
        heading.textContent = `${ovizDeckSelectedObjectIds.size} objects selected`;
        const format = document.createElement("details");
        format.className = "oviz-deck-inspector-section";
        format.open = true;
        const summary = document.createElement("summary");
        summary.textContent = "Format";
        const row = document.createElement("div");
        row.className = "oviz-deck-inspector-row";
        const colorField = document.createElement("label");
        colorField.className = "oviz-deck-color-field";
        const colorLabel = document.createElement("span");
        colorLabel.textContent = "Primary color";
        const color = document.createElement("input");
        color.type = "color";
        const first = ovizDeckSelectedObjects()[0];
        color.value = first && first.kind === "shape" ? first.fill_color : (first ? first.color : "#ffffff");
        color.addEventListener("input", () => {
          ovizDeckSelectedObjects({ unlockedOnly: true }).forEach((object) => {
            if (object.kind === "shape") object.fill_color = color.value;
            else object.color = color.value;
          });
          ovizDeckChanged("multi-color", { render: false });
          ovizDeckRenderAuthoringSlide();
        });
        colorField.append(colorLabel, color);
        row.append(colorField);
        const opacity = ovizDeckNumericField("Opacity", first ? first.opacity : 1, 0, 1, 0.05, (value) => {
          ovizDeckSelectedObjects({ unlockedOnly: true }).forEach((object) => { object.opacity = value; });
          ovizDeckChanged("multi-opacity");
        });
        format.append(summary, row, opacity);
        inspector.append(heading, format, ovizDeckRenderArrangeControls());
        ovizDeckEditorEl.append(inspector);
      }

      function ovizDeckRenderObjectsList(slide) {
        const details = document.createElement("details");
        details.className = "oviz-deck-objects-panel";
        details.open = ovizDeckObjectsPanelOpen;
        details.addEventListener("toggle", () => { ovizDeckObjectsPanelOpen = details.open; });
        const summary = document.createElement("summary");
        summary.textContent = `Objects (${slide.objects.length})`;
        const list = document.createElement("div");
        list.className = "oviz-deck-object-list";
        [...slide.objects].reverse().forEach((object) => {
          const row = document.createElement("div");
          row.className = "oviz-deck-object-row";
          row.draggable = !object.locked;
          row.dataset.objectId = object.id;
          row.dataset.selected = ovizDeckSelectedObjectIds.has(object.id) ? "true" : "false";
          const grip = document.createElement("span");
          grip.textContent = "⋮⋮";
          grip.className = "oviz-deck-object-grip";
          const name = document.createElement("button");
          name.type = "button";
          name.textContent = object.name || object.type;
          name.addEventListener("click", (event) => {
            event.preventDefault();
            event.stopPropagation();
            ovizDeckSelectObjectFromPointer(object, event);
          });
          name.title = object.name || object.type;
          const meta = document.createElement("small");
          meta.textContent = `${object.kind === "shape" ? object.shape_type.replace(/_/g, " ") : object.type}${object.group_id ? " · group" : ""}`;
          const lock = ovizDeckMakeButton(object.locked ? "🔒" : "○", () => ovizDeckLockObjects([object.id], !object.locked));
          lock.title = object.locked ? "Unlock object" : "Lock object";
          row.addEventListener("dragstart", (event) => {
            event.dataTransfer.setData("text/oviz-object", object.id);
            event.dataTransfer.effectAllowed = "move";
          });
          row.addEventListener("dragover", (event) => {
            event.preventDefault();
            event.dataTransfer.dropEffect = "move";
          });
          row.addEventListener("drop", (event) => {
            event.preventDefault();
            const sourceId = event.dataTransfer.getData("text/oviz-object");
            const source = ovizDeckObjectById(sourceId, slide);
            if (!source || source === object) return;
            const sourceIndex = slide.objects.indexOf(source);
            slide.objects.splice(sourceIndex, 1);
            const targetIndex = slide.objects.indexOf(object);
            slide.objects.splice(targetIndex + 1, 0, source);
            ovizDeckChanged("reorder-object");
          });
          const label = document.createElement("span");
          label.className = "oviz-deck-object-label";
          label.append(name, meta);
          row.append(grip, label, lock);
          list.append(row);
        });
        const guides = document.createElement("div");
        guides.className = "oviz-deck-guides-controls";
        const toggle = (label, key) => {
          const field = document.createElement("label");
          const input = document.createElement("input");
          input.type = "checkbox";
          input.checked = Boolean(ovizDeckProject.guides[key]);
          input.addEventListener("change", () => ovizDeckGuidesSet({ [key]: input.checked }));
          field.append(input, document.createTextNode(label));
          return field;
        };
        guides.append(toggle("Smart guides", "smart"), toggle("20 px grid", "grid"));
        details.append(summary, list, guides);
        ovizDeckEditorEl.append(details);
      }

      function ovizDeckBlockLines(block) {
        const lines = String(block && block.text || "").replace(/\r/g, "").split("\n");
        return lines.length ? lines : [""];
      }

      function ovizDeckPopulateBlockContent(element, block) {
        element.innerHTML = "";
        if (!block || block.list_style === "none") {
          element.textContent = block ? block.text : "";
          return;
        }
        const list = document.createElement(block.list_style === "number" ? "ol" : "ul");
        list.className = "oviz-deck-text-list";
        ovizDeckBlockLines(block).forEach((line) => {
          const item = document.createElement("li");
          item.textContent = line;
          list.append(item);
        });
        element.append(list);
      }

      function ovizDeckReadBlockContent(element, block) {
        if (block && block.list_style !== "none") {
          const items = Array.from(element.querySelectorAll(".oviz-deck-text-list > li"));
          if (items.length) return items.map((item) => item.textContent || "").join("\n");
        }
        return String(element.innerText !== undefined ? element.innerText : element.textContent || "").replace(/\r/g, "");
      }

      function ovizDeckApplyBlockStyle(element, block, scale = 1) {
        const scaleX = typeof scale === "object" ? ovizDeckFinite(scale.scaleX, 1) : scale;
        const scaleY = typeof scale === "object" ? ovizDeckFinite(scale.scaleY, 1) : scale;
        const fontScale = typeof scale === "object" ? ovizDeckFinite(scale.scale, Math.min(scaleX, scaleY)) : scale;
        element.style.left = `${block.x * scaleX}px`;
        element.style.top = `${block.y * scaleY}px`;
        element.style.width = `${block.width * scaleX}px`;
        element.style.height = `${block.height * scaleY}px`;
        element.style.transform = `rotate(${block.rotation || 0}deg)`;
        element.style.opacity = String(block.opacity === undefined ? 1 : block.opacity);
        element.style.color = block.color;
        element.style.fontSize = `${block.font_size * fontScale}px`;
        element.style.fontWeight = String(block.font_weight);
        const fontFallback = block.font_family === "serif"
          ? 'Georgia, "Times New Roman", serif'
          : (block.font_family === "mono"
            ? 'Menlo, Monaco, "Courier New", monospace'
            : '-apple-system, BlinkMacSystemFont, "SF Pro Display", "Helvetica Neue", sans-serif');
        element.style.fontFamily = ["sans", "serif", "mono"].includes(block.font_family)
          ? fontFallback
          : `"${String(block.font_family).replace(/["\\]/g, "")}", ${fontFallback}`;
        element.style.fontStyle = block.font_style;
        element.style.textDecoration = [block.underline ? "underline" : "", block.strikethrough ? "line-through" : ""]
          .filter(Boolean).join(" ") || "none";
        element.style.lineHeight = String(block.line_height);
        element.style.letterSpacing = `${ovizDeckFinite(block.character_spacing, 0) * fontScale}px`;
        element.style.textAlign = block.align;
        element.style.webkitTextStroke = block.text_outline
          ? `${ovizDeckFinite(block.text_outline_width, 1) * fontScale}px ${block.text_outline_color}`
          : "0 transparent";
        element.style.textShadow = block.text_shadow ? `0 ${2 * fontScale}px ${8 * fontScale}px rgba(0,0,0,.65)` : "none";
      }

      function ovizDeckAuthoringMetrics() {
        const rect = ovizDeckAuthoringLayerEl.getBoundingClientRect();
        return {
          rect,
          scaleX: Math.max(rect.width, 1) / OVIZ_DECK_WIDTH,
          scaleY: Math.max(rect.height, 1) / OVIZ_DECK_HEIGHT,
          scale: Math.min(Math.max(rect.width, 1) / OVIZ_DECK_WIDTH, Math.max(rect.height, 1) / OVIZ_DECK_HEIGHT),
        };
      }

      function ovizDeckBounds(objects) {
        if (!objects.length) return { x: 0, y: 0, width: 0, height: 0 };
        const left = Math.min(...objects.map((object) => object.x));
        const top = Math.min(...objects.map((object) => object.y));
        const right = Math.max(...objects.map((object) => object.x + object.width));
        const bottom = Math.max(...objects.map((object) => object.y + object.height));
        return { x: left, y: top, width: right - left, height: bottom - top };
      }

      function ovizDeckClearGuides() {
        ovizDeckGuideEls.forEach((element) => element.remove());
        ovizDeckGuideEls = [];
      }

      function ovizDeckShowGuide(axis, value, label = "") {
        if (!ovizDeckAuthoringLayerEl) return;
        const metrics = ovizDeckAuthoringMetrics();
        const line = document.createElement("div");
        line.className = `oviz-deck-guide-line oviz-deck-guide-${axis}`;
        if (axis === "x") line.style.left = `${value * metrics.scaleX}px`;
        else line.style.top = `${value * metrics.scaleY}px`;
        if (label) line.dataset.label = label;
        ovizDeckAuthoringLayerEl.append(line);
        ovizDeckGuideEls.push(line);
      }

      function ovizDeckSnapMove(bounds, x, y, excludedIds, event, gesture) {
        if ((event.metaKey || event.ctrlKey) || !ovizDeckProject.guides) return { x, y };
        const metrics = ovizDeckAuthoringMetrics();
        const thresholdX = 6 / metrics.scaleX;
        const thresholdY = 6 / metrics.scaleY;
        const hysteresisX = 8 / metrics.scaleX;
        const hysteresisY = 8 / metrics.scaleY;
        const movingX = [x, x + bounds.width / 2, x + bounds.width];
        const movingY = [y, y + bounds.height / 2, y + bounds.height];
        const xCandidates = [0, OVIZ_DECK_WIDTH / 2, OVIZ_DECK_WIDTH];
        const yCandidates = [0, OVIZ_DECK_HEIGHT / 2, OVIZ_DECK_HEIGHT];
        const others = (ovizDeckActiveSlide() ? ovizDeckActiveSlide().objects : [])
          .filter((object) => !excludedIds.has(object.id));
        others.forEach((object) => {
          xCandidates.push(object.x, object.x + object.width / 2, object.x + object.width);
          yCandidates.push(object.y, object.y + object.height / 2, object.y + object.height);
        });
        const sortedX = [...others].sort((a, b) => a.x - b.x);
        for (let index = 0; index < sortedX.length - 1; index += 1) {
          const left = sortedX[index];
          const right = sortedX[index + 1];
          if (right.x >= left.x + left.width) xCandidates.push((left.x + left.width + right.x - bounds.width) / 2);
        }
        const sortedY = [...others].sort((a, b) => a.y - b.y);
        for (let index = 0; index < sortedY.length - 1; index += 1) {
          const top = sortedY[index];
          const bottom = sortedY[index + 1];
          if (bottom.y >= top.y + top.height) yCandidates.push((top.y + top.height + bottom.y - bounds.height) / 2);
        }
        const closest = (moving, candidates, threshold) => {
          let match = null;
          moving.forEach((anchor, anchorIndex) => candidates.forEach((candidate) => {
            const delta = candidate - anchor;
            if (Math.abs(delta) <= threshold && (!match || Math.abs(delta) < Math.abs(match.delta))) {
              match = { delta, candidate, anchorIndex };
            }
          }));
          return match;
        };
        let snapX = null;
        let snapY = null;
        if (ovizDeckProject.guides.smart) {
          snapX = closest(movingX, xCandidates, thresholdX);
          snapY = closest(movingY, yCandidates, thresholdY);
          if (gesture.snapX && Math.abs(gesture.snapX.candidate - movingX[gesture.snapX.anchorIndex]) <= hysteresisX) snapX = gesture.snapX;
          if (gesture.snapY && Math.abs(gesture.snapY.candidate - movingY[gesture.snapY.anchorIndex]) <= hysteresisY) snapY = gesture.snapY;
        }
        const gridSize = ovizDeckFinite(ovizDeckProject.guides.grid_size, 20);
        if (!snapX && ovizDeckProject.guides.grid) {
          const grid = Math.round(x / gridSize) * gridSize;
          if (Math.abs(grid - x) <= thresholdX) snapX = { delta: grid - x, candidate: grid, anchorIndex: 0, grid: true };
        }
        if (!snapY && ovizDeckProject.guides.grid) {
          const grid = Math.round(y / gridSize) * gridSize;
          if (Math.abs(grid - y) <= thresholdY) snapY = { delta: grid - y, candidate: grid, anchorIndex: 0, grid: true };
        }
        gesture.snapX = snapX;
        gesture.snapY = snapY;
        ovizDeckClearGuides();
        if (snapX && !snapX.grid) ovizDeckShowGuide("x", snapX.candidate);
        if (snapY && !snapY.grid) ovizDeckShowGuide("y", snapY.candidate);
        return { x: x + (snapX ? snapX.delta : 0), y: y + (snapY ? snapY.delta : 0) };
      }

      function ovizDeckStartBlockGesture(event, block, element, mode, handle = "se") {
        if (!ovizDeckAuthoringLayerEl || block.locked) return;
        event.preventDefault();
        event.stopPropagation();
        if (!ovizDeckSelectedObjectIds.has(block.id)) {
          ovizDeckSetSelection([block.id], { primary: block.id, render: false });
        }
        if (mode === "move" && event.altKey) {
          const copies = ovizDeckDuplicateObjects(Array.from(ovizDeckSelectedObjectIds), { x: 0, y: 0 });
          block = ovizDeckObjectById(copies[0] && copies[0].id) || block;
        }
        const objects = ovizDeckSelectedObjects({ unlockedOnly: true });
        const bounds = ovizDeckBounds(objects);
        const metrics = ovizDeckAuthoringMetrics();
        ovizDeckDragState = {
          block,
          element,
          mode,
          handle,
          pointerId: event.pointerId,
          startX: event.clientX,
          startY: event.clientY,
          metrics,
          bounds,
          objects: objects.map((object) => ({ object, start: ovizDeckClone(object, {}) })),
          startAngle: Math.atan2(
            event.clientY - (metrics.rect.top + (bounds.y + bounds.height / 2) * metrics.scaleY),
            event.clientX - (metrics.rect.left + (bounds.x + bounds.width / 2) * metrics.scaleX),
          ),
          snapX: null,
          snapY: null,
          moved: false,
        };
        if (event.currentTarget && typeof event.currentTarget.setPointerCapture === "function") {
          try { event.currentTarget.setPointerCapture(event.pointerId); } catch (_err) {}
        }
      }

      function ovizDeckUpdateGestureElements(gesture) {
        gesture.objects.forEach(({ object }) => {
          const element = ovizDeckAuthoringLayerEl.querySelector(`[data-oviz-deck-object-id="${CSS.escape(object.id)}"]`);
          if (element) ovizDeckApplyBlockStyle(element, object, gesture.metrics);
        });
        ovizDeckRenderSelectionFrame();
      }

      function ovizDeckOnPointerMove(event) {
        if (ovizDeckMarqueeState && event.pointerId === ovizDeckMarqueeState.pointerId) {
          const state = ovizDeckMarqueeState;
          state.currentX = event.clientX;
          state.currentY = event.clientY;
          ovizDeckRenderMarquee();
          event.preventDefault();
          return;
        }
        const gesture = ovizDeckDragState;
        if (!gesture || event.pointerId !== gesture.pointerId) return;
        let dx = (event.clientX - gesture.startX) / gesture.metrics.scaleX;
        let dy = (event.clientY - gesture.startY) / gesture.metrics.scaleY;
        if (Math.abs(event.clientX - gesture.startX) > 1 || Math.abs(event.clientY - gesture.startY) > 1) {
          gesture.moved = true;
        }
        if (gesture.mode === "move") {
          if (event.shiftKey) {
            if (Math.abs(dx) >= Math.abs(dy)) dy = 0;
            else dx = 0;
          }
          let targetX = ovizDeckClamp(gesture.bounds.x + dx, 0, OVIZ_DECK_WIDTH - gesture.bounds.width, gesture.bounds.x);
          let targetY = ovizDeckClamp(gesture.bounds.y + dy, 0, OVIZ_DECK_HEIGHT - gesture.bounds.height, gesture.bounds.y);
          const snapped = ovizDeckSnapMove(
            gesture.bounds,
            targetX,
            targetY,
            new Set(gesture.objects.map(({ object }) => object.id)),
            event,
            gesture,
          );
          dx = snapped.x - gesture.bounds.x;
          dy = snapped.y - gesture.bounds.y;
          gesture.objects.forEach(({ object, start }) => {
            object.x = start.x + dx;
            object.y = start.y + dy;
          });
        } else if (gesture.mode === "rotate") {
          const centerX = gesture.metrics.rect.left + (gesture.bounds.x + gesture.bounds.width / 2) * gesture.metrics.scaleX;
          const centerY = gesture.metrics.rect.top + (gesture.bounds.y + gesture.bounds.height / 2) * gesture.metrics.scaleY;
          const angle = Math.atan2(event.clientY - centerY, event.clientX - centerX);
          let deltaDegrees = (angle - gesture.startAngle) * 180 / Math.PI;
          if (event.shiftKey) deltaDegrees = Math.round(deltaDegrees / 15) * 15;
          const radians = deltaDegrees * Math.PI / 180;
          const centerDesignX = gesture.bounds.x + gesture.bounds.width / 2;
          const centerDesignY = gesture.bounds.y + gesture.bounds.height / 2;
          gesture.objects.forEach(({ object, start }) => {
            const objectCenterX = start.x + start.width / 2;
            const objectCenterY = start.y + start.height / 2;
            const localX = objectCenterX - centerDesignX;
            const localY = objectCenterY - centerDesignY;
            object.x = centerDesignX + localX * Math.cos(radians) - localY * Math.sin(radians) - start.width / 2;
            object.y = centerDesignY + localX * Math.sin(radians) + localY * Math.cos(radians) - start.height / 2;
            object.rotation = start.rotation + deltaDegrees;
          });
        } else {
          const west = gesture.handle.includes("w");
          const east = gesture.handle.includes("e");
          const north = gesture.handle.includes("n");
          const south = gesture.handle.includes("s");
          let left = gesture.bounds.x + (west ? dx : 0);
          let right = gesture.bounds.x + gesture.bounds.width + (east ? dx : 0);
          let top = gesture.bounds.y + (north ? dy : 0);
          let bottom = gesture.bounds.y + gesture.bounds.height + (south ? dy : 0);
          if (event.altKey) {
            if (west) right -= dx;
            if (east) left -= dx;
            if (north) bottom -= dy;
            if (south) top -= dy;
          }
          let width = Math.max(20, right - left);
          let height = Math.max(20, bottom - top);
          if (event.shiftKey) {
            const aspect = gesture.bounds.width / Math.max(gesture.bounds.height, 1);
            if (Math.abs(width - gesture.bounds.width) >= Math.abs(height - gesture.bounds.height)) height = width / aspect;
            else width = height * aspect;
            if (west && !east) left = right - width;
            else right = left + width;
            if (north && !south) top = bottom - height;
            else bottom = top + height;
          }
          if (!(event.metaKey || event.ctrlKey) && ovizDeckProject.guides) {
            const thresholdX = 6 / gesture.metrics.scaleX;
            const thresholdY = 6 / gesture.metrics.scaleY;
            const others = (ovizDeckActiveSlide() ? ovizDeckActiveSlide().objects : [])
              .filter((object) => !gesture.objects.some((entry) => entry.object.id === object.id));
            let widthMatch = null;
            let heightMatch = null;
            if (ovizDeckProject.guides.smart) {
              others.forEach((object) => {
                if (Math.abs(object.width - width) <= thresholdX && (!widthMatch || Math.abs(object.width - width) < Math.abs(widthMatch - width))) widthMatch = object.width;
                if (Math.abs(object.height - height) <= thresholdY && (!heightMatch || Math.abs(object.height - height) < Math.abs(heightMatch - height))) heightMatch = object.height;
              });
            }
            const gridSize = ovizDeckFinite(ovizDeckProject.guides.grid_size, 20);
            if (widthMatch === null && ovizDeckProject.guides.grid) {
              const gridWidth = Math.round(width / gridSize) * gridSize;
              if (Math.abs(gridWidth - width) <= thresholdX) widthMatch = gridWidth;
            }
            if (heightMatch === null && ovizDeckProject.guides.grid) {
              const gridHeight = Math.round(height / gridSize) * gridSize;
              if (Math.abs(gridHeight - height) <= thresholdY) heightMatch = gridHeight;
            }
            if (widthMatch !== null) {
              width = Math.max(20, widthMatch);
              if (west && !east) left = right - width;
              else right = left + width;
            }
            if (heightMatch !== null) {
              height = Math.max(20, heightMatch);
              if (north && !south) top = bottom - height;
              else bottom = top + height;
            }
          }
          left = ovizDeckClamp(left, 0, OVIZ_DECK_WIDTH - 20, gesture.bounds.x);
          top = ovizDeckClamp(top, 0, OVIZ_DECK_HEIGHT - 20, gesture.bounds.y);
          width = Math.min(width, OVIZ_DECK_WIDTH - left);
          height = Math.min(height, OVIZ_DECK_HEIGHT - top);
          const scaleX = width / Math.max(gesture.bounds.width, 1);
          const scaleY = height / Math.max(gesture.bounds.height, 1);
          gesture.objects.forEach(({ object, start }) => {
            object.x = left + (start.x - gesture.bounds.x) * scaleX;
            object.y = top + (start.y - gesture.bounds.y) * scaleY;
            object.width = Math.max(20, start.width * scaleX);
            object.height = Math.max(20, start.height * scaleY);
          });
          ovizDeckClearGuides();
        }
        ovizDeckUpdateGestureElements(gesture);
        event.preventDefault();
      }

      function ovizDeckEndPointerGesture(event) {
        if (ovizDeckMarqueeState && (event.pointerId === undefined || event.pointerId === ovizDeckMarqueeState.pointerId)) {
          const state = ovizDeckMarqueeState;
          const metrics = ovizDeckAuthoringMetrics();
          const x1 = (Math.min(state.startX, state.currentX) - metrics.rect.left) / metrics.scaleX;
          const y1 = (Math.min(state.startY, state.currentY) - metrics.rect.top) / metrics.scaleY;
          const x2 = (Math.max(state.startX, state.currentX) - metrics.rect.left) / metrics.scaleX;
          const y2 = (Math.max(state.startY, state.currentY) - metrics.rect.top) / metrics.scaleY;
          const ids = (ovizDeckActiveSlide() ? ovizDeckActiveSlide().objects : []).filter((object) => (
            object.x < x2 && object.x + object.width > x1 && object.y < y2 && object.y + object.height > y1
          )).map((object) => object.id);
          ovizDeckMarqueeState = null;
          ovizDeckSetSelection(ids, { expandGroup: true });
          return;
        }
        if (!ovizDeckDragState || (event.pointerId !== undefined && event.pointerId !== ovizDeckDragState.pointerId)) return;
        if (!ovizDeckDragState.moved) {
          ovizDeckDragState = null;
          ovizDeckClearGuides();
          ovizDeckRenderEditor();
          ovizDeckRenderSelectionFrame();
          return;
        }
        const reason = ovizDeckDragState.mode === "move"
          ? "move-object"
          : (ovizDeckDragState.mode === "rotate" ? "rotate-object" : "resize-object");
        ovizDeckDragState = null;
        ovizDeckClearGuides();
        ovizDeckChanged(reason);
      }

      function ovizDeckRenderMarquee() {
        const state = ovizDeckMarqueeState;
        if (!state || !ovizDeckAuthoringLayerEl) return;
        let marquee = ovizDeckAuthoringLayerEl.querySelector(".oviz-deck-marquee");
        if (!marquee) {
          marquee = document.createElement("div");
          marquee.className = "oviz-deck-marquee";
          ovizDeckAuthoringLayerEl.append(marquee);
        }
        const rect = ovizDeckAuthoringLayerEl.getBoundingClientRect();
        marquee.style.left = `${Math.min(state.startX, state.currentX) - rect.left}px`;
        marquee.style.top = `${Math.min(state.startY, state.currentY) - rect.top}px`;
        marquee.style.width = `${Math.abs(state.currentX - state.startX)}px`;
        marquee.style.height = `${Math.abs(state.currentY - state.startY)}px`;
      }

      function ovizDeckMakeShapeVisual(block) {
        const svg = document.createElementNS("http://www.w3.org/2000/svg", "svg");
        svg.classList.add("oviz-deck-shape-visual");
        svg.setAttribute("viewBox", "0 0 100 100");
        svg.setAttribute("preserveAspectRatio", "none");
        const make = (tag, attributes) => {
          const node = document.createElementNS("http://www.w3.org/2000/svg", tag);
          Object.entries(attributes).forEach(([key, value]) => node.setAttribute(key, String(value)));
          return node;
        };
        let shape;
        if (block.shape_type === "ellipse") shape = make("ellipse", { cx: 50, cy: 50, rx: 48, ry: 48 });
        else if (block.shape_type === "triangle") shape = make("polygon", { points: "50,2 98,98 2,98" });
        else if (block.shape_type === "diamond") shape = make("polygon", { points: "50,2 98,50 50,98 2,50" });
        else if (["line", "arrow"].includes(block.shape_type)) {
          if (block.shape_type === "arrow") {
            const defs = make("defs", {});
            const marker = make("marker", { id: `arrow-${block.id}`, markerWidth: 10, markerHeight: 8, refX: 9, refY: 4, orient: "auto", markerUnits: "strokeWidth" });
            const head = make("path", { d: "M0,0 L10,4 L0,8 Z", fill: block.border_color });
            marker.append(head);
            defs.append(marker);
            svg.append(defs);
          }
          shape = make("line", { x1: 2, y1: 50, x2: 96, y2: 50 });
          if (block.shape_type === "arrow") shape.setAttribute("marker-end", `url(#arrow-${block.id})`);
        } else {
          shape = make("rect", {
            x: 2, y: 2, width: 96, height: 96,
            rx: block.shape_type === "rounded_rectangle" ? Math.min(48, block.corner_radius / Math.max(block.width, 1) * 100) : 0,
          });
        }
        shape.setAttribute("fill", ["line", "arrow"].includes(block.shape_type) ? "none" : block.fill_color);
        shape.setAttribute("fill-opacity", String(block.fill_opacity));
        shape.setAttribute("stroke", block.border_color);
        shape.setAttribute("stroke-width", String(Math.max(0.5, block.border_width)));
        shape.setAttribute("vector-effect", "non-scaling-stroke");
        if (block.border_style === "dashed") shape.setAttribute("stroke-dasharray", "9 6");
        if (block.border_style === "dotted") shape.setAttribute("stroke-dasharray", "2 5");
        if (block.shadow) svg.style.filter = "drop-shadow(0 7px 12px rgba(0,0,0,.36))";
        svg.append(shape);
        return svg;
      }

      function ovizDeckRenderSelectionFrame() {
        if (!ovizDeckAuthoringLayerEl) return;
        ovizDeckAuthoringLayerEl.querySelectorAll(".oviz-deck-selection-frame").forEach((element) => element.remove());
        const objects = ovizDeckSelectedObjects();
        if (!objects.length) return;
        const bounds = ovizDeckBounds(objects);
        const metrics = ovizDeckAuthoringMetrics();
        const frame = document.createElement("div");
        frame.className = "oviz-deck-selection-frame";
        frame.style.left = `${bounds.x * metrics.scaleX}px`;
        frame.style.top = `${bounds.y * metrics.scaleY}px`;
        frame.style.width = `${bounds.width * metrics.scaleX}px`;
        frame.style.height = `${bounds.height * metrics.scaleY}px`;
        if (objects.some((object) => object.locked)) frame.dataset.locked = "true";
        if (!objects.some((object) => object.locked)) {
          ["nw", "n", "ne", "e", "se", "s", "sw", "w"].forEach((handle) => {
            const button = document.createElement("button");
            button.type = "button";
            button.className = `oviz-deck-resize-handle oviz-deck-resize-${handle}`;
            button.setAttribute("aria-label", `Resize ${handle}`);
            button.addEventListener("pointerdown", (event) => ovizDeckStartBlockGesture(event, objects[0], frame, "resize", handle));
            frame.append(button);
          });
          const rotate = document.createElement("button");
          rotate.type = "button";
          rotate.className = "oviz-deck-rotate-handle";
          rotate.title = "Rotate (hold Shift for 15° increments)";
          rotate.addEventListener("pointerdown", (event) => ovizDeckStartBlockGesture(event, objects[0], frame, "rotate"));
          frame.append(rotate);
        }
        ovizDeckAuthoringLayerEl.append(frame);
      }

      function ovizDeckRenderAuthoringSlide() {
        if (!ovizDeckAuthoringLayerEl || !ovizDeckProject) return;
        const visible = ovizDeckEditorOpen && !ovizDeckPresenting && Boolean(ovizDeckActiveSlide());
        ovizDeckAuthoringLayerEl.dataset.visible = visible ? "true" : "false";
        ovizDeckAuthoringLayerEl.setAttribute("aria-hidden", visible ? "false" : "true");
        ovizDeckAuthoringLayerEl.innerHTML = "";
        ovizDeckGuideEls = [];
        if (!visible) return;
        const slide = ovizDeckActiveSlide();
        const metrics = ovizDeckAuthoringMetrics();
        ovizDeckAuthoringLayerEl.dataset.grid = ovizDeckProject.guides.grid ? "true" : "false";
        ovizDeckAuthoringLayerEl.style.setProperty("--oviz-deck-grid-x", `${ovizDeckProject.guides.grid_size * metrics.scaleX}px`);
        ovizDeckAuthoringLayerEl.style.setProperty("--oviz-deck-grid-y", `${ovizDeckProject.guides.grid_size * metrics.scaleY}px`);
        slide.objects.forEach((block) => {
          const wrapper = document.createElement("div");
          wrapper.className = "oviz-deck-author-block";
          wrapper.tabIndex = -1;
          wrapper.dataset.ovizDeckObjectId = block.id;
          wrapper.dataset.kind = block.kind;
          wrapper.dataset.selected = ovizDeckSelectedObjectIds.has(block.id) ? "true" : "false";
          wrapper.dataset.locked = block.locked ? "true" : "false";
          ovizDeckApplyBlockStyle(wrapper, block, metrics);
          wrapper.addEventListener("pointerdown", (event) => {
            event.stopPropagation();
            if (event.target && event.target.closest && event.target.closest(".oviz-deck-author-content[contenteditable='true']")) return;
            if ((event.shiftKey || event.metaKey || event.ctrlKey) || !ovizDeckSelectedObjectIds.has(block.id)) {
              ovizDeckSelectObjectFromPointer(block, event, { render: false });
            }
            if (block.kind === "shape" && !block.locked) ovizDeckStartBlockGesture(event, block, wrapper, "move");
          });
          wrapper.addEventListener("click", (event) => {
            event.stopPropagation();
            ovizDeckSetSelection(Array.from(ovizDeckSelectedObjectIds), { primary: block.id, render: false, expandGroup: false });
            ovizDeckRenderEditor();
            ovizDeckRenderSelectionFrame();
          });

          if (block.kind === "shape") wrapper.append(ovizDeckMakeShapeVisual(block));
          const content = document.createElement("div");
          content.className = "oviz-deck-author-content";
          content.contentEditable = block.kind === "text" ? "true" : "false";
          content.spellcheck = true;
          ovizDeckPopulateBlockContent(content, block);
          const selectContentBlock = () => {
            if (!ovizDeckSelectedObjectIds.has(block.id)) {
              ovizDeckSetSelection([block.id], { primary: block.id, render: false });
              ovizDeckRenderEditor();
              ovizDeckRenderSelectionFrame();
            }
          };
          content.addEventListener("pointerdown", (event) => {
            event.stopPropagation();
            if (block.kind === "shape" && content.contentEditable !== "true") {
              if ((event.shiftKey || event.metaKey || event.ctrlKey) || !ovizDeckSelectedObjectIds.has(block.id)) {
                ovizDeckSelectObjectFromPointer(block, event, { render: false });
              }
              if (!block.locked && event.detail < 2) ovizDeckStartBlockGesture(event, block, wrapper, "move");
            } else {
              selectContentBlock();
            }
          });
          content.addEventListener("click", (event) => { event.stopPropagation(); selectContentBlock(); });
          content.addEventListener("dblclick", (event) => {
            if (block.kind !== "shape") return;
            event.stopPropagation();
            content.contentEditable = "true";
            content.focus();
          });
          content.addEventListener("focus", selectContentBlock);
          content.addEventListener("input", () => {
            block.text = ovizDeckReadBlockContent(content, block);
            if (block.kind === "text" && block.text.trim()) block.name = block.text.trim().split(/\r?\n/)[0].slice(0, 36);
            ovizDeckChanged("edit-text", { render: false });
          });
          content.addEventListener("keydown", (event) => {
            if (event.key === "Escape") {
              event.preventDefault();
              event.stopPropagation();
              if (block.kind === "shape") content.contentEditable = "false";
              try { wrapper.focus({ preventScroll: true }); } catch (_err) { wrapper.focus(); }
              return;
            }
            event.stopPropagation();
          });
          wrapper.addEventListener("keydown", (event) => {
            if (event.key === "Enter" && block.kind === "shape") {
              event.preventDefault();
              content.contentEditable = "true";
              content.focus();
            }
          });

          if (block.kind === "text" && !block.locked) {
            const move = document.createElement("button");
            move.type = "button";
            move.className = "oviz-deck-block-move";
            move.textContent = "Drag";
            move.title = "Drag text box (Option-drag to duplicate)";
            move.addEventListener("pointerdown", (event) => ovizDeckStartBlockGesture(event, block, wrapper, "move"));
            wrapper.append(content, move);
          } else {
            wrapper.append(content);
          }
          ovizDeckAuthoringLayerEl.append(wrapper);
        });
        ovizDeckRenderSelectionFrame();
      }

      function ovizDeckRevealBlock(block) {
        const element = document.createElement("div");
        element.className = `oviz-deck-present-block oviz-deck-present-${block.type} oviz-deck-present-${block.kind}`;
        ovizDeckPopulateBlockContent(element, block);
        if (block.kind === "shape") element.prepend(ovizDeckMakeShapeVisual(block));
        ovizDeckApplyBlockStyle(element, block);
        return element;
      }

      function ovizDeckRenderRevealSlides() {
        if (!ovizDeckRevealRootEl || !ovizDeckProject) return;
        const slidesEl = ovizDeckRevealRootEl.querySelector(".slides");
        if (!slidesEl) return;
        slidesEl.innerHTML = "";
        ovizDeckProject.slides.forEach((slide, index) => {
          const section = document.createElement("section");
          section.dataset.ovizDeckSlideId = slide.id;
          section.setAttribute("aria-label", slide.name);
          slide.objects.forEach((block) => section.append(ovizDeckRevealBlock(block)));
          if (slide.notes) {
            const notes = document.createElement("aside");
            notes.className = "notes";
            notes.textContent = slide.notes;
            section.append(notes);
          }
          slidesEl.append(section);
        });
        ovizDeckRenderGeneration += 1;
        if (ovizDeckReveal && typeof ovizDeckReveal.sync === "function") ovizDeckReveal.sync();
      }

      function ovizDeckLoadRevealAssets() {
        if (window.Reveal) return Promise.resolve(window.Reveal);
        if (ovizDeckRevealPromise) return ovizDeckRevealPromise;
        const version = String(
          ovizDeckProject && ovizDeckProject.reveal && ovizDeckProject.reveal.version
          || OVIZ_DECK_REVEAL_VERSION
        );
        const cssUrl = `https://cdn.jsdelivr.net/npm/reveal.js@${version}/dist/reveal.css`;
        const jsUrl = `https://cdn.jsdelivr.net/npm/reveal.js@${version}/dist/reveal.js`;
        if (!document.querySelector('link[data-oviz-deck-reveal="true"]')) {
          const link = document.createElement("link");
          link.rel = "stylesheet";
          link.href = cssUrl;
          link.dataset.ovizDeckReveal = "true";
          document.head.append(link);
        }
        ovizDeckRevealPromise = new Promise((resolve, reject) => {
          const existing = document.querySelector('script[data-oviz-deck-reveal="true"]');
          if (existing) {
            existing.addEventListener("load", () => window.Reveal ? resolve(window.Reveal) : reject(new Error("Reveal.js did not initialize.")), { once: true });
            existing.addEventListener("error", () => reject(new Error("Reveal.js could not be loaded.")), { once: true });
            return;
          }
          const script = document.createElement("script");
          script.src = jsUrl;
          script.dataset.ovizDeckReveal = "true";
          script.onload = () => window.Reveal ? resolve(window.Reveal) : reject(new Error("Reveal.js did not initialize."));
          script.onerror = () => reject(new Error("Reveal.js could not be loaded."));
          document.head.append(script);
        });
        return ovizDeckRevealPromise;
      }

      async function ovizDeckEnsureReveal() {
        if (ovizDeckReveal) return ovizDeckReveal;
        const RevealConstructor = await ovizDeckLoadRevealAssets();
        ovizDeckReveal = new RevealConstructor(ovizDeckRevealRootEl, {
          embedded: true,
          controls: false,
          controlsTutorial: false,
          progress: false,
          hash: false,
          keyboard: false,
          touch: true,
          overview: false,
          center: false,
          transition: String(ovizDeckProject.reveal.transition || "fade"),
          backgroundTransition: String(ovizDeckProject.reveal.background_transition || "fade"),
          width: 1600,
          height: 900,
          margin: 0,
        });
        ovizDeckReveal.on("slidechanged", (event) => {
          if (!ovizDeckPresenting) return;
          const index = Number(event && event.indexh);
          if (!Number.isFinite(index) || index === ovizDeckActiveIndex) return;
          ovizDeckGoTo(index, { fromReveal: true }).catch(() => {});
        });
        await ovizDeckReveal.initialize();
        return ovizDeckReveal;
      }

      function ovizDeckHasSlides() {
        return Boolean(ovizDeckProject && ovizDeckProject.available && ovizDeckProject.enabled && ovizDeckProject.slides.length);
      }

      function ovizDeckIsPresenting() {
        return Boolean(ovizDeckPresenting && ovizDeckHasSlides());
      }

      async function ovizDeckSetPresenting(enabled) {
        const next = Boolean(enabled);
        if (!next) {
          ovizDeckPresenting = false;
          root.dataset.deckPresenting = "false";
          if (ovizDeckRevealRootEl) {
            ovizDeckRevealRootEl.dataset.visible = "false";
            ovizDeckRevealRootEl.setAttribute("aria-hidden", "true");
          }
          ovizDeckRenderAuthoringSlide();
          ovizDeckEvent("deck-present-complete", { presenting: false, index: ovizDeckActiveIndex });
          return false;
        }
        if (!ovizDeckHasSlides()) return false;
        ovizDeckEditorOpen = false;
        ovizDeckPresenting = true;
        root.dataset.deckPresenting = "true";
        ovizDeckRenderEditor();
        ovizDeckRenderAuthoringSlide();
        ovizDeckRenderRevealSlides();
        if (ovizDeckRevealRootEl) {
          ovizDeckRevealRootEl.dataset.visible = "true";
          ovizDeckRevealRootEl.setAttribute("aria-hidden", "false");
        }
        ovizDeckEvent("deck-present-start", { presenting: true, index: ovizDeckActiveIndex });
        let reveal;
        try {
          reveal = await ovizDeckEnsureReveal();
          reveal.slide(ovizDeckActiveIndex);
        } catch (error) {
          ovizDeckPresenting = false;
          root.dataset.deckPresenting = "false";
          root.dataset.deckError = String(error && error.message || error);
          if (ovizDeckRevealRootEl) ovizDeckRevealRootEl.dataset.visible = "false";
          ovizDeckSetStatus(String(error && error.message || error), true);
          ovizDeckEvent("deck-error", { message: String(error && error.message || error) });
          console.error("Oviz Reveal presentation failed.", error);
          if (presentationModeEnabled) setPresentationMode(false);
          return false;
        }
        root.dataset.deckError = "";
        try {
          await ovizDeckGoTo(ovizDeckActiveIndex, { fromReveal: true });
        } catch (error) {
          const message = String(error && error.message || error);
          root.dataset.deckStateError = message;
          ovizDeckSetStatus(message, true);
          ovizDeckEvent("deck-error", { message, domain: "state" });
        }
        return true;
      }

      function ovizDeckNext() {
        if (!ovizDeckIsPresenting()) return Promise.resolve({ boundary: true });
        const next = Math.min(ovizDeckActiveIndex + 1, ovizDeckProject.slides.length - 1);
        if (next === ovizDeckActiveIndex) return Promise.resolve({ boundary: true, index: next });
        return ovizDeckGoTo(next);
      }

      function ovizDeckPrevious() {
        if (!ovizDeckIsPresenting()) return Promise.resolve({ boundary: true });
        const previous = Math.max(ovizDeckActiveIndex - 1, 0);
        if (previous === ovizDeckActiveIndex) return Promise.resolve({ boundary: true, index: previous });
        return ovizDeckGoTo(previous);
      }

      function ovizDeckSetEditorOpen(enabled) {
        ovizDeckEditorOpen = Boolean(enabled)
          && Boolean(ovizDeckProject && ovizDeckProject.available)
          && !ovizDeckPresenting
          && !mobileModeEnabled;
        if (ovizDeckEditorOpen) {
          ovizDeckEnsureStateEditMode();
          if (typeof setManualLabelMenuOpen === "function") setManualLabelMenuOpen(false);
          if (ovizStatesShellEl) {
            ovizStatesShellEl.dataset.open = "false";
            const statesToggle = ovizStatesShellEl.querySelector(".oviz-states-toggle");
            if (statesToggle) statesToggle.textContent = "States ▸";
          }
          if (typeof setControlsDrawerOpen === "function") setControlsDrawerOpen(false);
        }
        root.dataset.deckEditing = ovizDeckEditorOpen ? "true" : "false";
        if (ovizDeckButtonEl) {
          ovizDeckButtonEl.dataset.active = ovizDeckEditorOpen ? "true" : "false";
          ovizDeckButtonEl.setAttribute("aria-pressed", ovizDeckEditorOpen ? "true" : "false");
          ovizDeckButtonEl.setAttribute("aria-expanded", ovizDeckEditorOpen ? "true" : "false");
          ovizDeckButtonEl.textContent = ovizDeckEditorOpen ? "Slides ▾" : "Slides ▸";
        }
        if (ovizDeckEditorEl) {
          ovizDeckEditorEl.addEventListener("pointerdown", (event) => event.stopPropagation());
        }
        ovizDeckRenderEditor();
        ovizDeckRenderAuthoringSlide();
        return ovizDeckEditorOpen;
      }

      function ovizDeckHandleObjectShortcut(event) {
        if (!ovizDeckEditorOpen || ovizDeckPresenting) return false;
        const key = String(event.key || "");
        const lower = key.toLowerCase();
        const target = event.target;
        if (target && (
          target.isContentEditable
          || (typeof target.closest === "function" && target.closest("[contenteditable='true'],input,textarea,select"))
        )) return false;
        const command = event.metaKey || event.ctrlKey;
        const selected = ovizDeckSelectedObjects({ unlockedOnly: true });
        if ((key === "Backspace" || key === "Delete") && selected.length) {
          ovizDeckRemoveObjects(selected.map((object) => object.id));
        } else if (["ArrowLeft", "ArrowRight", "ArrowUp", "ArrowDown"].includes(key) && selected.length && !command && !event.altKey) {
          const amount = event.shiftKey ? 10 : 1;
          selected.forEach((object) => {
            if (key === "ArrowLeft") object.x = Math.max(0, object.x - amount);
            if (key === "ArrowRight") object.x = Math.min(OVIZ_DECK_WIDTH - object.width, object.x + amount);
            if (key === "ArrowUp") object.y = Math.max(0, object.y - amount);
            if (key === "ArrowDown") object.y = Math.min(OVIZ_DECK_HEIGHT - object.height, object.y + amount);
          });
          ovizDeckChanged("nudge-object");
        } else if (command && lower === "d" && selected.length) {
          ovizDeckDuplicateObjects(selected.map((object) => object.id));
        } else if (command && lower === "c" && selected.length) {
          ovizDeckObjectClipboard = ovizDeckClone(selected, []);
        } else if (command && lower === "x" && selected.length) {
          ovizDeckObjectClipboard = ovizDeckClone(selected, []);
          ovizDeckRemoveObjects(selected.map((object) => object.id));
        } else if (command && lower === "v" && ovizDeckObjectClipboard.length) {
          const slide = ovizDeckActiveSlide();
          const copies = ovizDeckObjectClipboard.map((source) => {
            const copy = Object.assign({}, source, {
              id: ovizDeckUuid("object"),
              group_id: null,
              x: ovizDeckClamp(source.x + 20, 0, OVIZ_DECK_WIDTH - source.width, source.x),
              y: ovizDeckClamp(source.y + 20, 0, OVIZ_DECK_HEIGHT - source.height, source.y),
            });
            return copy;
          });
          slide.objects.push(...copies);
          ovizDeckSetSelection(copies.map((object) => object.id), { primary: copies[0].id, render: false, expandGroup: false });
          ovizDeckChanged("paste-object");
        } else if (command && lower === "a") {
          ovizDeckSetSelection((ovizDeckActiveSlide() ? ovizDeckActiveSlide().objects : []).map((object) => object.id), { expandGroup: false });
        } else if (command && event.altKey && lower === "g" && event.shiftKey) {
          ovizDeckUngroupObjects();
        } else if (command && event.altKey && lower === "g") {
          ovizDeckGroupObjects();
        } else if (command && lower === "l" && selected.length) {
          ovizDeckLockObjects(selected.map((object) => object.id), !event.altKey);
        } else if (command && event.shiftKey && lower === "f") {
          ovizDeckArrangeObjects(event.altKey ? "forward" : "front");
        } else if (command && event.shiftKey && lower === "b") {
          ovizDeckArrangeObjects(event.altKey ? "backward" : "back");
        } else if (key === "Tab") {
          const objects = ovizDeckActiveSlide() ? ovizDeckActiveSlide().objects : [];
          if (!objects.length) return false;
          const current = objects.findIndex((object) => object.id === ovizDeckSelectedBlockId);
          const next = (current + (event.shiftKey ? -1 : 1) + objects.length) % objects.length;
          ovizDeckSetSelection([objects[next].id], { primary: objects[next].id });
        } else {
          return false;
        }
        event.preventDefault();
        event.stopPropagation();
        return true;
      }

      function ovizDeckOnKeyDown(event) {
        if (event.defaultPrevented) return;
        ovizDeckHandleObjectShortcut(event);
      }

      function ovizDeckInstallPublicApi() {
        const viewer = window.Oviz && window.Oviz.get && window.Oviz.get(root.id);
        if (!viewer) return;
        viewer.deck = {
          list: ovizDeckList,
          current: () => ovizDeckClone(ovizDeckActiveSlide(), null),
          goTo: ovizDeckGoTo,
          add: ovizDeckAdd,
          updateState: ovizDeckUpdateState,
          rename: ovizDeckRename,
          duplicate: ovizDeckDuplicate,
          move: ovizDeckMove,
          remove: ovizDeckRemove,
          addText: ovizDeckAddBlock,
          removeText: ovizDeckRemoveBlock,
          undo: ovizDeckUndo,
          redo: ovizDeckRedo,
          setEditorOpen: ovizDeckSetEditorOpen,
          present: () => setPresentationMode(true),
          exit: () => setPresentationMode(false),
          exportHtml: (options = {}) => ovizExportStatesHtml(options),
          objects: {
            list: ovizDeckObjectList,
            get: ovizDeckObjectGet,
            select: (ids) => ovizDeckSetSelection(ids),
            update: ovizDeckObjectUpdate,
            remove: ovizDeckRemoveObjects,
            duplicate: ovizDeckDuplicateObjects,
            addShape: ovizDeckAddShape,
            group: ovizDeckGroupObjects,
            ungroup: ovizDeckUngroupObjects,
            lock: (ids) => ovizDeckLockObjects(ids, true),
            unlock: (ids) => ovizDeckLockObjects(ids, false),
            arrange: ovizDeckArrangeObjects,
            align: ovizDeckAlignObjects,
            distribute: ovizDeckDistributeObjects,
          },
          guides: {
            get: ovizDeckGuidesGet,
            set: ovizDeckGuidesSet,
          },
        };
      }

      async function initializeOvizDeck() {
        ovizDeckProject = ovizDeckNormalizeProject(sceneSpec.deck);
        sceneSpec.deck = ovizDeckExportSpec({ embedded: Boolean(sceneSpec.deck && sceneSpec.deck.embedded) });
        ovizDeckUndoStack = [];
        ovizDeckRedoStack = [];
        ovizDeckHistoryPresent = ovizDeckHistorySnapshot();
        ovizDeckButtonEl = root.querySelector(".oviz-three-deck-editor");
        ovizDeckEditorEl = root.querySelector(".oviz-deck-editor");
        ovizDeckAuthoringLayerEl = root.querySelector(".oviz-deck-authoring-layer");
        ovizDeckRevealRootEl = root.querySelector(".oviz-deck-reveal");
        root.dataset.deckPresenting = "false";
        root.dataset.deckEditing = "false";
        root.dataset.deckAvailable = ovizDeckProject.available ? "true" : "false";
        if (ovizDeckButtonEl) {
          ovizDeckButtonEl.addEventListener("click", () => {
            ovizDeckSetEditorOpen(!ovizDeckEditorOpen);
          });
        }
        if (ovizDeckAuthoringLayerEl) {
          ovizDeckAuthoringLayerEl.addEventListener("pointerdown", (event) => {
            event.stopPropagation();
            if (event.target !== ovizDeckAuthoringLayerEl) return;
            ovizDeckSetSelection([], { render: false });
            ovizDeckMarqueeState = {
              pointerId: event.pointerId,
              startX: event.clientX,
              startY: event.clientY,
              currentX: event.clientX,
              currentY: event.clientY,
            };
            try { ovizDeckAuthoringLayerEl.setPointerCapture(event.pointerId); } catch (_err) {}
            ovizDeckRenderEditor();
            ovizDeckRenderSelectionFrame();
            event.preventDefault();
          });
        }
        window.addEventListener("pointermove", ovizDeckOnPointerMove);
        window.addEventListener("pointerup", ovizDeckEndPointerGesture);
        window.addEventListener("pointercancel", ovizDeckEndPointerGesture);
        window.addEventListener("keydown", ovizDeckOnKeyDown);
        window.addEventListener("resize", () => {
          if (ovizDeckEditorOpen && !ovizDeckPresenting) ovizDeckRenderAuthoringSlide();
        });
        ovizDeckInstallPublicApi();
        ovizDeckRenderEditor();
        ovizDeckRenderAuthoringSlide();
        root.dataset.deckReady = "true";
        ovizDeckEvent("deck-ready", { slides: ovizDeckList() });
      }
""".strip()
