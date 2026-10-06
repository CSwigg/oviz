The Oviz viewer
===============

Every figure is one standalone HTML file. Open it in a modern browser (Chrome,
Edge, Firefox, Safari 16.4+) from disk or any static host. The 3D scene renders
offline; Sky backgrounds stream from HiPS servers through Aladin Lite.

``make_plot`` writes this viewer by default. ``viewer="classic"`` writes the
previous Three.js runtime byte for byte; Slides and Paper exist only there.


Controls
--------

=====================  ==========================================================
Navigate               Drag to orbit, right- or Ctrl-drag to pan, scroll or
                       middle-drag to zoom. Hold **W A S D** to orbit and tilt,
                       **Shift + W A S D** to fly, **Q E** to zoom, **R F** to move
                       up and down. **O** auto-orbits, **Home** resets.
Time                   **Space** plays, **← →** step (Shift: five frames),
                       **< >** change speed, **0** returns to today, **J** draws
                       motion trails.
Inspect                Click an object for its details; Shift-click a second one to
                       measure. Double-click an object to follow it, elsewhere to
                       return to the LSR.
Select                 Shift-drag (or **L**) to lasso; **C** toggles the filter,
                       **⌘Z** undoes. A lasso also clips the dust outside it.
Search                 **⌘K** or **/** searches objects (with catalogue aliases),
                       layers, views and actions; *Find … on the sky* asks SIMBAD.
Layers                 **1–9** toggle, Shift + number or double-click solos, **T**
                       hides or restores all, **[ ]** resize points, **Shift + L**
                       opens the layers panel.
Sky                    **V** flies to the Sun and looks out; **Shift + B** picks the
                       background survey.
Views                  **N** saves a view, **Y** opens the story, **P** presents.
Capture and share      **I** saves a PNG (up to 4×); the capture menu records MP4 or
                       WebM; *Copy link* shares the current view.
Display                **G** grid, **B** Sky background, **Z** hides the interface,
                       **M** fullscreen, **U** switches Focus/Detailed mode,
                       **?** lists every shortcut.
=====================  ==========================================================

On touch screens panels become bottom sheets, key rows take a tap (toggle), a
long press (solo) or *Hide all*, a lasso is one finger stroke, and a button
beside search hides the interface.


Concepts
--------

**Time.** Positions are interpolated between stored frames on the GPU, so
playback is smooth at any speed. The time slider shows when the visible objects
were born and marks today and every saved view.

**Camera anchor.** Oviz scenes are centred on the Local Standard of Rest: the
frame follows a circular reference orbit, so the LSR stays at the origin. The
camera orbits and zooms about its anchor and moves with it through time, while
the Sun and the clusters move around it. The anchor is the LSR by default; the
Sun, any followed object, or none (*Free*) are alternatives under *Display
settings* or ``make_plot(camera_anchor=...)``. Pans keep the anchor; **Home**
re-anchors.

**Sky view.** The 3D scene is registered exactly to Aladin Lite (TAN
projection, matching centre, field and Galactic-north-up orientation). Cluster
markers crossfade into their member stars, which inherit their cluster's orbit,
birth time, visibility and colour, and can show motion arrows. The background
picker orders whole-sky surveys by wavelength, crossfades between neighbours,
and searches every CDS survey; the sky lens shows another survey inside a
circle; right-click asks SIMBAD what is there. The Sun is not drawn in Sky view.
Today's surveys fade out as time moves away from the present.

**Views and presenting.** A view (a *State*) captures the whole viewer:
camera and anchor, time, mode, layer and volume styles, Sky layers, filters,
lasso selection, widgets and display settings. Transitions animate every
continuous property and land exactly; lasso selections crossfade. While
presenting, only the layers in view are listed in a quiet key. ⌘S saves the
figure with its views; *Export presentation* writes a present-only copy.
Unsaved edits autosave in the browser.

**Links.** *Copy link* stores the current view as a compressed difference from
how the figure opens (typically a few hundred characters) in the URL hash.
When views are saved, links carry them too (``w=``), and *Copy presentation
link* opens presenting them (``p=1``).

**Modes.** *Focus* (the default) shows the figure, its key and one quiet bar
that fades while the pointer rests. *Detailed* keeps the layers panel, the
inspector and the full transport in view. Choose with
``make_plot(viewer_mode=...)`` or ``?mode=detailed``.

**Widgets and actions.** The widgets button opens the *Birth tree*, *Relative
SFH* and *Notes* panels; they are saved with views. Figures made with
``make_plot(actions=...)`` show an action bar of scripted camera and time
moves.

**Trace options for authors.** Set them with the trace's ``meta``:

- ``{"oviz_role": "guide"}``: an annotation (labels, a length bar, a model
  curve) with its own switch under *Guides*, not a row of the key.
- ``{"oviz_pickable": False}``: drawn, but never hovered, clicked, lassoed or
  searched.
- ``{"oviz_hide_in_sky": True}``: left out of Sky view and its key, like the
  Sun and spiral-arm traces.

``make_plot(show_sun=True)`` (the default) adds the Sun when the data have no
"Sun" trace. ``make_plot(spiral_arm_models=...)`` adds published spiral arms
(see :mod:`oviz.spiral_models`), hidden until switched on.


Upgrading existing figures
--------------------------

Figures written by earlier Oviz releases carry their full scene, so they
convert without re-running the analysis; saved States migrate too.

.. code-block:: bash

   python -m oviz.viewer.upgrade old_figure.html new_figure.html [--mode detailed] [--camera-anchor sun]

.. code-block:: python

   import oviz
   oviz.upgrade_html("old_figure.html", "new_figure.html")


AR on iPhone and iPad
---------------------

**View in AR** opens the figure in Apple's AR Quick Look as a tabletop
time-lapse: visible objects move along their tracks, visible dust volumes
become see-through image slices that appear at the present day, and a base
plate carries distance rings and a time readout. *View this moment in AR*
gives a still of the current time. The Share menu saves the same USDZ.

To ship a prebuilt model next to a published page, pass
``OvizFigure(..., ar_model="figure.usdz")`` and save the file from the Share
menu (or ``window.Oviz.viewer.arModel()``). Serve ``.usdz`` as
``model/vnd.usdz+zip``; check it with ``usdchecker --arkit``.


Browser API
-----------

``window.Oviz.viewer`` (or ``window.Oviz.get(rootId)``) exposes:

``setTime(t)``, ``play()``, ``pause()``, ``time``
    Timeline control.
``setViewMode("3d" | "sky")``
    Fly between Galactic 3D and Sky.
``cameraAnchor``, ``setCameraAnchor(anchor)``
    ``"lsr"``, ``"sun"``, ``"free"``, or ``{trace, index}`` for an object.
``screenshot({scale})``
    A PNG ``Blob`` of the WebGL layer.
``arModel({moment})``, ``viewInAr({moment})``
    The AR model as ``{blob, summary}``, and the Quick Look hand-off.
``viewLink({views, present})``, ``presentationLink()``
    Promises of shareable links.
``getState()``, ``applyState(state, {instant, duration_ms, easing, keepCamera})``
    Capture the whole viewer, and restore a State exactly (animated unless
    ``instant``; the promise resolves on arrival).
``states``
    ``list()``, ``goTo(indexOrId)``, ``next()``, ``previous()``, ``add()``,
    ``capture()``, ``present(on)``, ``exportHtml({presentOnly, currentOnly,
    download})`` (``download: false`` returns the HTML text).
``stats()``
    Frame and volume timings.

Embedding pages can drive a figure in an iframe with ``postMessage``
(``{type: "oviz", command, args}``; commands ``goTo``, ``next``,
``previous``, ``present``, ``setTime``, ``play``, ``pause``,
``setViewMode``) and listen for ``oviz:transition-start`` and
``oviz:transition-end`` events.


How it works (for maintainers)
------------------------------

**Bundle.** :func:`oviz.viewer.compile_scene_spec` turns the scene spec into a
JSON manifest plus binary blobs. Positions over time are ``(frames, objects,
3)`` float32, collapsed when static or rigidly moving; per-object attributes
are stored once; volumes keep uint8 voxels plus a block-max occupancy grid
(event-density volumes keep their events and are rebuilt in the browser);
member stars become position and tangential-velocity columns. Blobs are
byte-shuffled, gzipped and base64-encoded, and decode in priority order so
traces draw before volumes arrive. The July 25 figure is 24 MB (48 MB
classic) and draws its first frame in about 70 ms.

**Engine.** Plain WebGL2: one instanced draw per trace; frame positions read
from a float texture and spline-interpolated in the vertex shader; the star
PSF as one premultiplied halo + core pass; ray-marched volumes at half
resolution while moving and cached when still; GPU picking on each point's
core; demand-driven rendering; full recovery after a lost WebGL context.

**Source layout.** ``oviz/viewer/web/src`` holds ES modules that also run in
Node for tests (``node --test tests/viewer_js/*.test.mjs``).
:mod:`oviz.viewer.build` bundles them into one script; each module keeps its
own scope and exported names must be unique. ``styles/`` holds the design
system (tokens in ``00-tokens.css``) and ``template.html`` the page skeleton,
mirrored by ``src/app/export.js``.

- ``core/``: maths, colour, bundle loader, events.
- ``engine/``: GL helpers, camera, anchors, controls, renderer, frame textures.
- ``layers/``: points, lines, labels, image planes, volumes, KDE volumes.
- ``sky/``: Aladin, survey catalogue, layer stack, picker, lens, SIMBAD,
  member stars.
- ``app/``: viewer, timeline, state model, States, links, drafts, export, API.
- ``ui/``: shell, layers panel, dock, inspector, palette, story, filter,
  lasso, recorder, widget host and widgets, actions, axes.

Per-object GPU state is a bitfield: 1 dimmed and 2 hidden by the filter, 4
replaced by member stars in Sky view, 8 dimmed and 16 hidden by a lasso (for
CPU readers; the shaders draw the lasso from a separate per-object pair so
selections can crossfade), 32 dimmed by the Birth tree.
