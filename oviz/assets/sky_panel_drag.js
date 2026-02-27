(function () {
  "use strict";

  var DRAG_HANDLE_ID = "oviz-sky-drag-handle";
  var PANEL_ID = "oviz-sky-container";
  var RESIZE_HANDLE_CLASS = "oviz-sky-resize-handle";
  var MARGIN_PX = 6;
  var MIN_WIDTH_PX = 220;
  var MIN_HEIGHT_PX = 220;

  var pointerState = null;
  var disabledPlots = [];
  var panelObserver = null;
  var panelWatchTimer = null;

  function getPanel() {
    return document.getElementById(PANEL_ID);
  }

  function isDraggablePanel(panel) {
    if (!panel) {
      return false;
    }
    var cs = window.getComputedStyle(panel);
    if (cs.display === "none" || cs.visibility === "hidden") {
      return false;
    }
    var rect = panel.getBoundingClientRect();
    if (rect.width <= 0 || rect.height <= 0) {
      return false;
    }
    var nearFullscreen =
      Math.abs(rect.width - window.innerWidth) <= 2 &&
      Math.abs(rect.height - window.innerHeight) <= 2;
    return !nearFullscreen;
  }

  function clampPosition(left, top, width, height) {
    var maxLeft = Math.max(MARGIN_PX, window.innerWidth - width - MARGIN_PX);
    var maxTop = Math.max(MARGIN_PX, window.innerHeight - height - MARGIN_PX);
    return {
      left: Math.min(Math.max(MARGIN_PX, left), maxLeft),
      top: Math.min(Math.max(MARGIN_PX, top), maxTop),
    };
  }

  function disablePlotInteractions() {
    disabledPlots = Array.from(document.querySelectorAll(".js-plotly-plot"));
    disabledPlots.forEach(function (el) {
      el.setAttribute("data-oviz-prev-pointer-events", el.style.pointerEvents || "");
      el.style.pointerEvents = "none";
    });
  }

  function restorePlotInteractions() {
    disabledPlots.forEach(function (el) {
      if (!el) {
        return;
      }
      var prev = el.getAttribute("data-oviz-prev-pointer-events");
      if (prev === null) {
        el.style.pointerEvents = "";
      } else {
        el.style.pointerEvents = prev;
      }
      el.removeAttribute("data-oviz-prev-pointer-events");
    });
    disabledPlots = [];
  }

  function setResizeHandlesVisibility() {
    var panel = getPanel();
    var handles = document.querySelectorAll("." + RESIZE_HANDLE_CLASS);
    var visible = isDraggablePanel(panel);
    handles.forEach(function (h) {
      h.style.display = visible ? "block" : "none";
    });
  }

  function startPanelWatch() {
    if (panelObserver || panelWatchTimer) {
      return;
    }
    panelWatchTimer = window.setInterval(function () {
      var panel = getPanel();
      if (!panel) {
        return;
      }
      setResizeHandlesVisibility();
      panelObserver = new MutationObserver(function () {
        setResizeHandlesVisibility();
      });
      panelObserver.observe(panel, { attributes: true, attributeFilter: ["style", "class"] });
      window.clearInterval(panelWatchTimer);
      panelWatchTimer = null;
    }, 200);
  }

  function startPointer(ev) {
    if (typeof ev.button === "number" && ev.button !== 0) {
      return;
    }

    var resizeHandle = ev.target.closest("." + RESIZE_HANDLE_CLASS);
    var dragHandle = ev.target.closest("#" + DRAG_HANDLE_ID);
    if (!resizeHandle && !dragHandle) {
      return;
    }

    var panel = getPanel();
    if (!isDraggablePanel(panel)) {
      return;
    }

    var rect = panel.getBoundingClientRect();
    pointerState = {
      mode: resizeHandle ? "resize" : "drag",
      dir: resizeHandle ? (resizeHandle.dataset.dir || "se").toLowerCase() : null,
      startX: ev.clientX,
      startY: ev.clientY,
      startLeft: rect.left,
      startTop: rect.top,
      startWidth: rect.width,
      startHeight: rect.height,
      startRight: rect.right,
      startBottom: rect.bottom,
      pointerId: ev.pointerId,
      handle: resizeHandle || dragHandle,
    };

    panel.style.left = rect.left + "px";
    panel.style.top = rect.top + "px";
    panel.style.right = "auto";
    panel.style.bottom = "auto";
    panel.style.transform = "none";
    panel.style.width = rect.width + "px";
    panel.style.height = rect.height + "px";

    if (pointerState.mode === "resize") {
      panel.style.aspectRatio = "auto";
    }

    if (pointerState.handle && typeof pointerState.handle.setPointerCapture === "function" && ev.pointerId !== undefined) {
      try {
        pointerState.handle.setPointerCapture(ev.pointerId);
      } catch (_err) {
      }
    }

    disablePlotInteractions();
    document.body.style.userSelect = "none";
    ev.preventDefault();
    ev.stopPropagation();
  }

  function resizeRectForPointer(state, clientX, clientY) {
    var dx = clientX - state.startX;
    var dy = clientY - state.startY;
    var dir = state.dir || "se";
    var availableW = Math.max(80, window.innerWidth - 2 * MARGIN_PX);
    var availableH = Math.max(80, window.innerHeight - 2 * MARGIN_PX);
    var minW = Math.min(MIN_WIDTH_PX, availableW);
    var minH = Math.min(MIN_HEIGHT_PX, availableH);

    var minX = MARGIN_PX;
    var maxX = window.innerWidth - MARGIN_PX;
    var minY = MARGIN_PX;
    var maxY = window.innerHeight - MARGIN_PX;

    var left = state.startLeft;
    var right = state.startRight;
    var top = state.startTop;
    var bottom = state.startBottom;

    if (dir.indexOf("w") !== -1) {
      left = state.startLeft + dx;
      left = Math.min(left, state.startRight - minW);
      left = Math.max(minX, left);
      right = state.startRight;
    } else if (dir.indexOf("e") !== -1) {
      left = state.startLeft;
      right = state.startRight + dx;
      right = Math.max(right, state.startLeft + minW);
      right = Math.min(maxX, right);
    }

    if (dir.indexOf("n") !== -1) {
      top = state.startTop + dy;
      top = Math.min(top, state.startBottom - minH);
      top = Math.max(minY, top);
      bottom = state.startBottom;
    } else if (dir.indexOf("s") !== -1) {
      top = state.startTop;
      bottom = state.startBottom + dy;
      bottom = Math.max(bottom, state.startTop + minH);
      bottom = Math.min(maxY, bottom);
    }

    return {
      left: left,
      top: top,
      width: Math.max(minW, right - left),
      height: Math.max(minH, bottom - top),
    };
  }

  function movePointer(ev) {
    if (!pointerState) {
      return;
    }
    var panel = getPanel();
    if (!panel) {
      return;
    }

    if (pointerState.mode === "drag") {
      var left = pointerState.startLeft + (ev.clientX - pointerState.startX);
      var top = pointerState.startTop + (ev.clientY - pointerState.startY);
      var clamped = clampPosition(left, top, pointerState.startWidth, pointerState.startHeight);
      panel.style.left = clamped.left + "px";
      panel.style.top = clamped.top + "px";
    } else {
      var resized = resizeRectForPointer(pointerState, ev.clientX, ev.clientY);
      panel.style.left = resized.left + "px";
      panel.style.top = resized.top + "px";
      panel.style.width = resized.width + "px";
      panel.style.height = resized.height + "px";
      panel.style.right = "auto";
      panel.style.bottom = "auto";
      panel.style.transform = "none";
      pointerState.startWidth = resized.width;
      pointerState.startHeight = resized.height;
    }

    ev.preventDefault();
    ev.stopPropagation();
  }

  function endPointer(ev) {
    if (!pointerState) {
      return;
    }

    var handle = pointerState.handle;
    var pointerId = pointerState.pointerId;
    if (handle && typeof handle.releasePointerCapture === "function" && pointerId !== undefined) {
      try {
        handle.releasePointerCapture(pointerId);
      } catch (_err) {
      }
    }

    pointerState = null;
    document.body.style.userSelect = "";
    restorePlotInteractions();
    if (ev) {
      ev.preventDefault();
      ev.stopPropagation();
    }
  }

  function onResize() {
    var panel = getPanel();
    setResizeHandlesVisibility();
    if (!isDraggablePanel(panel)) {
      return;
    }

    var rect = panel.getBoundingClientRect();
    var maxWidth = Math.max(80, window.innerWidth - 2 * MARGIN_PX);
    var maxHeight = Math.max(80, window.innerHeight - 2 * MARGIN_PX);
    var width = Math.min(rect.width, maxWidth);
    var height = Math.min(rect.height, maxHeight);
    var left = rect.left;
    var top = rect.top;

    if (left + width > window.innerWidth - MARGIN_PX) {
      left = window.innerWidth - MARGIN_PX - width;
    }
    if (top + height > window.innerHeight - MARGIN_PX) {
      top = window.innerHeight - MARGIN_PX - height;
    }
    left = Math.max(MARGIN_PX, left);
    top = Math.max(MARGIN_PX, top);

    if (
      Math.abs(left - rect.left) > 0.5 ||
      Math.abs(top - rect.top) > 0.5 ||
      Math.abs(width - rect.width) > 0.5 ||
      Math.abs(height - rect.height) > 0.5
    ) {
      panel.style.left = left + "px";
      panel.style.top = top + "px";
      panel.style.width = width + "px";
      panel.style.height = height + "px";
      panel.style.right = "auto";
      panel.style.bottom = "auto";
      panel.style.transform = "none";
      panel.style.aspectRatio = "auto";
    }
  }

  document.addEventListener("pointerdown", startPointer, true);
  document.addEventListener("pointermove", movePointer, true);
  document.addEventListener("pointerup", endPointer, true);
  document.addEventListener("pointercancel", endPointer, true);
  window.addEventListener("blur", endPointer);
  window.addEventListener("resize", onResize);
  startPanelWatch();
  setResizeHandlesVisibility();
})();
