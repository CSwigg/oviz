The Oviz viewer
===============

Every figure is one standalone HTML file. Open it in a modern browser (Chrome,
Edge, Firefox, Safari 16.4+) from disk or any static host. The 3D scene
renders fully offline. Sky backgrounds stream from HiPS servers through
Aladin Lite.

``Animate3D.make_plot`` writes this viewer by default. Pass
``viewer="classic"`` to get the previous Three.js runtime byte-for-byte.
Slides, Paper, AR and the dendrogram widget remain classic-only.


Using a figure
--------------

**Navigate.** Drag to orbit, Shift-drag or right-drag to pan, and scroll to
zoom toward the cursor. The classic keyboard controls are kept: hold
**W A S D** to orbit and tilt, **Shift + W A S D** to fly, **Q E** to zoom out
and in, and **R F** to move up and down. In Sky view the same keys look around
and **Q E** change the field of view. Double-click an object to fly to it, or
empty space to reset the view (also **Home**). **O** auto-orbits.

**Time.** Time is continuous: positions are interpolated between the stored
frames on the GPU, so playback is smooth at any speed. Use **Space** to play,
**←/→** to step a frame (Shift for five), **< >** to change speed and **0** to
return to the present day. The timeline marks t = 0 and each saved view.

**Inspect.** Click an object for its details (age today and at time t, member
count, distance, Galactic and ICRS coordinates, position). From the inspector
you can fly to it, *follow* it through time, draw its full *orbit trail*, or
copy its coordinates. Shift-click a second object to measure the 3D separation
and the angular separation seen from the Sun.

**Select.** Press **L** (or **X**) to lasso objects; the rest dims. **C**
turns the dimming off and on while keeping the selection, and **⌘Z** undoes
the last selection change. Frame, isolate or export the selection as CSV
(names, ages, members, coordinates and positions at time t).

**Search.** Press **⌘K** (Ctrl K), or **/**, to search every object, including
catalogue aliases, as well as layers, saved views and actions.

**Layers.** Toggle layers (**1–9**), solo them with a double-click (or
Shift + 1–9), show or hide everything with **T**, and edit their colour,
colour-by-value colormap, opacity and size. **[ ]** shrink and grow every
point; **Shift + L** opens the layers panel. Volumes expose
colormap, stretch, data window, opacity, density gain and samples. The
*Distribution* section shows a live histogram of age, members, distance or
height; brush a range to dim or hide everything else.

**Sky.** Press **V** or use the 3D/Sky switch to fly to the Sun and look out.
The WebGL scene is registered to Aladin Lite exactly (TAN projection, matching
centre, horizontal field and Galactic-north-up orientation). Cluster markers
crossfade into their member stars, which inherit the parent's orbit, birth time,
visibility and colour. Survey layers can be shown, faded and reordered.

**Views & story (Y).** Save the complete viewer as a view with **N**. That
includes camera and 3D/Sky mode, time, layer styles, volumes, Sky layers,
filters and display settings. Transitions animate every continuous property
together and then assign the saved values exactly. Views get thumbnails and
captions, can be reordered by dragging, and play as a presentation
(**P**, arrows to navigate, **Esc** or **P** to exit).

Your unsaved edits autosave in the browser. **⌘S** saves the figure as a new
HTML file. *Export presentation* writes a present-only file that opens
straight into view 1.

**Capture.** **I** saves a PNG, at 1×, 2× or 4× resolution, or copies it to
the clipboard. The capture menu records video to MP4 or WebM: live, a
time-lapse of the whole timeline, or a tour of your saved views. The Sky
background is composited into both.

**Share.** *Copy link to this view* encodes time, mode, group and camera in the
URL hash. Opening the link restores that view.

**Display.** **G** toggles the Galactic grid, **B** the Sky background,
**Z** hides the interface and **M** goes fullscreen. **?** lists every
shortcut.

**Themes and devices.** A light theme restyles the interface; the data stage
stays a night sky because starlight and emission volumes are calibrated against
it. On phones, panels become bottom sheets and touch targets grow. The desktop
layout is unchanged.


Upgrading existing figures
--------------------------

Any HTML figure written by an earlier Oviz release can be converted without
re-running its pipeline. Saved States are migrated too.

.. code-block:: bash

   python -m oviz.viewer.upgrade old_figure.html new_figure.html

.. code-block:: python

   import oviz
   oviz.upgrade_html("old_figure.html", "new_figure.html")


Browser API
-----------

``window.Oviz.viewer`` (or ``window.Oviz.get(rootId)``) exposes:

``setTime(t)``, ``play()``, ``pause()``, ``time``
    Timeline control.
``setViewMode("3d" | "sky")``
    Fly between Galactic 3D and Sky.
``screenshot({scale})``
    A PNG ``Blob`` of the WebGL layer.
``states``
    ``list()``, ``goTo(indexOrId)``, ``next()``, ``previous()``, ``add()``,
    ``capture()``, ``present(on)``, and ``exportHtml({presentOnly,
    currentOnly, download})``. With ``download: false`` it returns the HTML
    text.
``stats()``
    Frame and volume timings, and the boot timeline.
``viewer`` / ``ui``
    The underlying objects, for debugging.


How it works (for maintainers)
------------------------------

**Bundle.** :func:`oviz.viewer.compile_scene_spec` converts the scene spec into
a JSON manifest plus binary blobs:

- *Positions over time* are stored as ``(frames, objects, 3)`` float32. They
  are collapsed to one frame when the data are static, or to one frame plus a
  per-frame offset when every vertex moves rigidly (e.g. re-centred reference
  circles).
- *Per-object attributes* (colour, age, member count, sky coordinates, names)
  are stored once rather than per frame.
- *Volumes* keep their raw uint8 voxels plus a precomputed block-max occupancy
  grid.
- *Member stars* become Cartesian position and tangential-velocity columns.

Each blob is gzipped, byte-shuffled for floats, and base64-encoded once. The
browser decodes blobs with native primitives in priority order, so traces
render before volumes arrive. The July 25 figure shrinks from 48 MB to 24 MB,
its first frame appears in about 70 ms, and all volumes are resident in about
170 ms. The classic viewer took about 1.8 s to start.

**Engine.** The engine is plain WebGL2 with no three.js:

- one instanced draw per trace;
- a vertex shader that fetches frame positions from a float texture and
  interpolates them with Catmull-Rom splines;
- the legacy star PSF, drawn as a single premultiplied halo + core pass;
- analytic ray-box volume marching with the legacy optical model, drawn at
  half resolution while moving and cached when settled;
- GPU picking;
- demand-driven rendering, so idle figures cost nothing;
- rebuilding of every GPU resource after a lost WebGL context.

**Source layout.** ``oviz/viewer/web/src`` holds ordinary ES modules, runnable
by ``node --test tests/viewer_js/*.test.mjs``. :mod:`oviz.viewer.build`
bundles them into one classic script: each module gets its own scope, and
exported names must be unique. ``oviz/viewer/web/styles`` holds the design
system; its tokens live in ``00-tokens.css``. ``oviz/viewer/web/template.html``
is the page skeleton, mirrored by ``src/app/export.js`` for self-export.

- ``core/``: maths, colour, bundle loader, events.
- ``engine/``: GL helpers, camera, controls, renderer, frame textures.
- ``layers/``: points, lines, labels, image planes, volumes.
- ``sky/``: Aladin integration and member stars.
- ``app/``: the viewer, timeline, state model, States, export, the public API.
- ``ui/``: the shell, layers panel, dock, inspector, palette, story, filter,
  recorder.

Per-object GPU state is a bitfield: 1 = dimmed, 2 = hidden, 4 = replaced by
member stars in Sky view.
