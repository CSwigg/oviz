The Oviz viewer
===============

Every figure is one standalone HTML file. Open it in a modern browser (Chrome,
Edge, Firefox, Safari 16.4+) from disk or any static host. The 3D scene
renders fully offline. Sky backgrounds stream from HiPS servers through
Aladin Lite.

``Animate3D.make_plot`` writes this viewer by default. Pass
``viewer="classic"`` to get the previous Three.js runtime byte-for-byte.
Slides and Paper remain classic-only; the classic Dendrogram and Age KDE
widgets are the *Birth tree* and *Relative SFH* widgets here. AR has moved to
the new viewer (see *AR on iPhone and iPad* below).


Using a figure
--------------

**Navigate.** Drag to orbit, right-drag (or Ctrl-drag) to pan, and scroll or
middle-drag to zoom. The classic keyboard controls are kept: hold
**W A S D** to orbit and tilt, **Shift + W A S D** to fly, **Q E** to zoom out
and in, and **R F** to move up and down. In Sky view the same keys look around
and **Q E** change the field of view. Double-click a cluster to follow it
(see *Camera anchor*), or anywhere else to return to the LSR; **Home** resets
the view. **O** auto-orbits.
*Display settings* set the scroll and key speed and the auto-orbit speed and
direction.

**Camera anchor.** The camera is anchored to the Local Standard of Rest by
default. Oviz scenes are LSR-centred: positions are relative to a reference
orbit that starts at the Sun with the LSR's velocity, so the LSR is the origin
at every time. The camera orbits the LSR, zooming keeps it centred, and the
camera stays with it through time while the Sun and the clusters move around
it. The scale readout says what the camera is anchored to (for example
"Eye 6 kpc from the LSR"). You can change the anchor under *Display settings*:

- **Sun** rides along with the Sun.
- **Follow** in the inspector rides along with any object.
- **Free** lets the wheel zoom toward the pointer.

The anchor is the camera's frame. Panning and flying keep it: the orbit
centre can sit off the anchor, and the camera keeps that offset as the anchor
moves through time. Double-click a cluster to make it the anchor (the camera
centres on it, keeping the zoom, and moves with it through time), and
double-click anywhere else to return to the LSR, centred on the dust or line
under the pointer, or on the LSR itself. **Home** re-anchors and resets the
view. Views and shared links remember the anchor. Figures choose
theirs with ``make_plot(..., camera_anchor="lsr" | "sun" | "free")`` or
``OvizFigure(..., camera_anchor=...)``. Scenes centred on a focus group's orbit
have no LSR at their origin, so they offer no LSR anchor.

**Time.** Time is continuous: positions are interpolated between the stored
frames on the GPU, so playback is smooth at any speed. Use **Space** to play,
**←/→** to step a frame (Shift for five), **< >** to change speed and **0** to
return to the present day. The timeline marks t = 0 and each saved view.
**J** turns on motion trails: every object draws a fading streak along its
own orbit for the last few Myr, pointing back along the playback direction
(set the length under *Display*). Trails are saved with views.

**Inspect.** Hovering an object only highlights it; its name tag appears
when you click it. Click an object for its details (age today and at time t, member
count, distance, Galactic and ICRS coordinates, position). A reticle locks on
to the object and a leader line ties it to the inspector. The inspector charts
the object's distance from the Sun across the whole timeline; drag along the
chart to scrub time. From the inspector you can fly to it, *follow* it through
time, draw its full *orbit trail*, or copy its coordinates. Shift-click a second object to measure the 3D separation
and the angular separation seen from the Sun.

**Select.** Shift-drag to lasso objects, or press **L** (or **X**) to keep
lassoing until **L** or **Esc**. On touch screens choose *Lasso select* in
the more menu and draw one stroke with a finger. The lasso selects the
volumes too: everything else hides (or dims, from the selection bar),
including the dust outside the lasso, and a lasso around dust alone selects
just that dust. *Lasso selects volumes too* in *Display settings* turns that
off. **C** turns the filter off and on while keeping the selection, and **⌘Z** undoes the last selection change. Frame the
selection or export it as CSV (names, ages, members, coordinates and
positions at time t).

**Search.** Press **⌘K** (Ctrl K), or **/**, to search every object, including
catalogue aliases, as well as layers, saved views and actions. In figures
with a Sky view, *Find … on the sky* looks any name up in SIMBAD and turns
the Sky view to it.

**Layers.** Volumes are layers like the data traces: they share one list
and the key, are numbered in the same order, and follow the group presets.
Toggle layers (**1–9**), solo them with a double-click (or Shift + 1–9),
show or hide everything with **T** (press it again to bring back the mix
that was on), and edit their colour,
colour-by-value colormap, opacity and size. **[ ]** shrink and grow every
point; **Shift + L** opens the layers panel. Volumes expose
colormap, stretch, data window, opacity, density gain and samples; an
event-density volume (such as the ccSN KDE) also has its trailing *Time
window*, and volumes with a colour scale show it beside the key. The
*Distribution* section shows a live histogram of age, members, distance or
height; brush a range to dim or hide everything else.

**Guides.** Annotations that are not data (structure names, a length bar, a
model curve, axis markers) can be marked as guides with
``meta={"oviz_role": "guide"}`` on their static trace. Guides stay out of
the key and of **T** / **1–9**; each gets its own switch in the layers
panel's *Guides* section, next to the Galactic grid, and saved views can
show or hide them like any layer.

Two more trace options suit points that sample a surface or a model rather
than objects. ``meta={"oviz_pickable": False}`` draws the points but never
picks them: no hover, no click, and the lasso, search and distribution filter
leave them out. ``meta={"oviz_hide_in_sky": True}`` keeps a layer out of Sky
view and its key, the way the Sun is (a shell around the Sun, seen from
inside, would cover the sky); it comes back in 3D as it was.

**Sky.** Press **V** or use the 3D/Sky switch to fly to the Sun and look out.
The corner readout gives the Galactic l, b at the centre of the view and the
field of view. Since Sky view looks out from the Sun, the Sun is not drawn
there and leaves the key and the layers panel. It comes back, as you left
it, in 3D.

A figure can open in Sky: give the opening state the classic Sun view,
``global_controls={"camera_view_mode": "earth", "camera_fov": ...,
"earth_view_return_camera_state": {...}}`` with a camera just behind the Sun
looking along the line of sight. It opens with that gaze and field of view,
and **V** (or Home in 3D) leads to the return view. Sky-start figures made by
earlier releases open this way once upgraded.

Every figure ``make_plot`` writes has the Sun: when the data bring no trace
named "Sun", Oviz adds one at the heliocentric origin, integrated like the
clusters (so it moves with the Sun's peculiar motion in the LSR frame) and
shown in every layer group. Pass ``show_sun=False`` to leave it out.
The WebGL scene is registered to Aladin Lite exactly (TAN projection, matching
centre, horizontal field and Galactic-north-up orientation). Cluster markers
crossfade into their member stars, which inherit the parent's orbit, birth time,
visibility and colour, and can show their motion as arrows. The survey is
today's sky, so it fades out as time moves away from the present.

The background button in Sky view shows the survey behind the data (its
thumbnail and name). It, **Shift + B** from either view, or *Sky
background…* in the more menu opens the picker, which chooses the background
the way a maps app picks its map type:
whole-sky thumbnails in wavelength order, a wavelength slider that crossfades
between neighbouring surveys, a blend strength for each survey in view, and a
search over every survey the CDS publishes (a HiPS ID or URL works too). The
layers panel's *Sky background* section lists what is in view, with opacity,
stretch, colours, cuts and order for each, and *Remove* for surveys you
added. The *sky lens* shows another wavelength inside a circle: drag its
ring, scroll over it to resize, step through bands with its arrows, and
double-click it to make that survey the background. Right-click the sky to
ask SIMBAD what is there.

**Widgets.** The widgets button opens floating panels: the *Birth tree*
(which clusters were born near older ones; hover a branch to highlight it in
3D, click to pin it), the *Relative SFH* (a density of birth times per layer;
brush it to filter by age) and *Notes* (text labels in the scene). Panels
drag by their title, resize from the corner, become bottom sheets on phones,
and are saved with views.

**Actions.** Figures made with ``make_plot(actions=...)`` show an action bar.
Each button runs its steps (go to a view, crossfade to a legend group, move
the camera, play time); click it again to go back to where it started. Any
camera input interrupts it.

**Views & story (Y).** Save the complete viewer as a view with **N**. That
includes camera and 3D/Sky mode, time, layer styles, volumes, Sky layers,
filters and display settings. Transitions animate every continuous property
together and then assign the saved values exactly. When two views have
different lasso selections, the objects, member stars and dust that leave
the selection fade out over the first three quarters of the flight. Those
that join fade in over the last three quarters, so the view never blinks
empty. Moving on mid-flight (the next arrow, or grabbing the camera) carries
the fade on from wherever it was, without a jump. Views get thumbnails and
captions, can be reordered by dragging, and play as a presentation
(**P**, arrows to navigate, **Esc** or **P** to exit). While presenting,
progress segments along the top fill as the camera flies to each view, and a
3D reel of view thumbnails rises above the caption whenever the pointer
moves; click a segment or a card to jump. The key stays in the corner as a
quiet figure key: only the layers in view (present-day dust maps drop out
away from t = 0), in small type, and not clickable. On phones the key and
scale bar sit just above the caption card.

Your unsaved edits autosave in the browser. **⌘S** saves the figure as a new
HTML file. *Export presentation* writes a present-only file that opens
straight into view 1. *Export current view only* writes a file that opens at
this view, and its Home returns there.

**Capture.** **I** saves a PNG, at 1×, 2× or 4× resolution, or copies it to
the clipboard. The capture menu records video to MP4 or WebM: live, a
time-lapse of the whole timeline, or a tour of your saved views. The Sky
background is composited into both.

**Share.** *Copy link to this view* puts the whole view in the URL hash, the
same snapshot a saved view keeps: time, mode, camera and anchor, which layers
and volumes are shown and how they are styled, the Sky background (which
surveys, their blend, stretch, colours and cuts, and the lens), filters, a
lasso selection with the outline that clips the dust, notes, open widgets,
the selected object and the display settings. The link stores only what
differs from how the figure opens, compressed, so a typical one is a few
hundred to about 1,500 characters. Opening it restores the view exactly. The
camera and time are also kept in readable form (``t=``, ``v=``, ``c=``), and
links made by older figures still open.

When you have saved views, the link carries them too (``w=``, compressed,
without thumbnails, which are redrawn as each view is visited). The receiver
gets the same story to explore or present, kept apart from views their own
browser has saved for the figure. *Copy presentation link* (Share menu, the
story panel's *Export* menu, or search) makes a link that opens presenting
those views from the first one (``p=1``). A typical link with a handful of
views is under 2,000 characters. For a file, use *Export presentation* or
**⌘S**.

**Modes.** Figures open in *Focus* mode: just the figure, its key and one
quiet, see-through bar with play, time, layers, views, 3D/Sky and a "more"
menu. The bar fades while the pointer rests and returns the moment it moves.
Clicking an object shows a small label beside it; *Details* opens the full
inspector. Press **U** (or the panel button in the bar) for *Detailed* mode:
the same design with everything in view, including the layers panel, a
details panel for the selection, the toolbar, the full transport (frame
steps, speed, loop, labelled ticks) and the 3D/Sky switch. On phones both
modes use bottom sheets that swipe down to dismiss.

The time slider shows a faint silhouette of when the visible objects were
born (brighter where time has already played), the present day, and each
saved view as a dot you can click; a bubble reads out the time under the
pointer. Choose the mode a figure opens in with
``make_plot(..., viewer_mode="detailed")`` (or ``OvizFigure(..., mode=...)``);
``?mode=detailed`` in the URL overrides it.

**Motion.** Figures open with a short camera dolly into the home view;
panels fade and slide a few pixels. Everything animated respects the system
*reduce motion* setting.

**Display.** **G** toggles the Galactic grid (the axes are in *Display
settings*). Its longitude lines and labels are drawn from today's Sun, so
they fade out within 1 Myr of the present day. **B** toggles the Sky background,
**Z** (or *Hide interface* in the ⋯ menu) hides the interface and **M** goes
fullscreen. On touch screens a button beside search hides the interface in one
tap, and a *Show interface* pill brings it back. **?** lists every shortcut.

**Themes and devices.** A light theme restyles the interface; the data stage
stays pure black because starlight and emission volumes are calibrated against
it. On phones, panels become bottom sheets and touch targets grow. On touch
screens the key's rows are finger-sized: tap a layer to show or hide it,
press and hold to show it alone, and use *Hide all* at the top of the key to
hide or show everything. The desktop layout is unchanged.


Upgrading existing figures
--------------------------

Any HTML figure written by an earlier Oviz release can be converted without
re-running its pipeline. Saved States are migrated too, with their lasso
selections (matched by cluster name), filters, open widgets, sky layers and
apertures, and notes.

.. code-block:: bash

   python -m oviz.viewer.upgrade old_figure.html new_figure.html [--mode detailed] [--camera-anchor sun]

.. code-block:: python

   import oviz
   oviz.upgrade_html("old_figure.html", "new_figure.html")


AR on iPhone and iPad
---------------------

On an iPhone or iPad, **View in AR** (in the ⋯ and Share menus, and in
search) opens the figure in Apple's AR Quick Look as a tabletop model. On any
device, the Share menu saves the same model as a USDZ file.

The model is what is on screen, as a **time-lapse** through the whole
timeline:

- every visible object is a small glowing sphere in its displayed colour and
  size, moving along its track and growing in at its birth (like the
  viewer's birth fade);
- each visible dust volume becomes three stacks of see-through image slices,
  integrated from the voxels with the renderer's own opacity model, and fades
  in as the timeline reaches today (dust maps are present-day);
- a dark base plate carries distance rings around the Sun, the direction of
  the Galactic centre and a time readout that follows the animation.

It holds briefly at the start and for a few seconds at the present day.
**View this moment in AR** (figures with a timeline) opens a still of the
time you have scrubbed to instead, labelled with that time. The model is
centred on the figure, with most present-day objects inside a radius of
0.4 m and Galactic north up; objects that roam farther earlier in time are
hidden until they come into view. Pinch in AR to scale it.

For a published page, ship a prebuilt time-lapse next to it so Quick Look
opens a plain file instead of one generated in the browser:

.. code-block:: python

   figure = OvizFigure(bundle=bundle, ar_model="figure.usdz")
   figure.write_html("figure.html")
   # Save figure.usdz next to it: Share → "Save AR time-lapse (USDZ)",
   # or window.Oviz.viewer.arModel() in the browser.

Serve ``.usdz`` files as ``model/vnd.usdz+zip`` (GitHub Pages does).
Check a model with ``usdchecker --arkit figure.usdz`` and preview it with
``usdrecord`` (both ship with macOS).


Browser API
-----------

``window.Oviz.viewer`` (or ``window.Oviz.get(rootId)``) exposes:

``setTime(t)``, ``play()``, ``pause()``, ``time``
    Timeline control.
``setViewMode("3d" | "sky")``
    Fly between Galactic 3D and Sky.
``cameraAnchor``, ``setCameraAnchor(anchor)``
    The camera anchor: ``"lsr"``, ``"sun"``, ``"free"``, or
    ``{trace, index}`` for an object (``{kind: ...}`` objects are accepted
    too). Setting one flies the orbit centre onto it.
``screenshot({scale})``
    A PNG ``Blob`` of the WebGL layer.
``arModel({moment})``, ``viewInAr({moment})``
    The AR time-lapse (or, with ``moment: true``, a still of the current
    time) as a USDZ (``{blob, summary}``), and the AR Quick Look hand-off.
``viewLink({views, present})``, ``presentationLink()``
    Promises of shareable links: the current view, with the saved views
    unless ``views: false``; ``presentationLink()`` (or ``present: true``)
    opens presenting them.
``getState()``, ``applyState(state, {instant, duration_ms, easing, keepCamera})``
    Capture the whole viewer as a State, and restore one exactly (animated
    unless ``instant``; returns a promise that resolves on arrival).
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
  grid. Event-density volumes instead keep their events (time, x, y, z) and
  the smoothing recipe; the viewer rebuilds the density at each time sample,
  bit-identical to the classic runtime, keeps the samples as half-float
  textures, and builds the rest while the page is idle.
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
- ``sky/``: Aladin integration, the survey catalogue, layer stack and
  background picker, the sky lens, SIMBAD look-ups, and member stars.
- ``app/``: the viewer, timeline, state model, States, export, the public API.
- ``ui/``: the shell, layers panel, dock, inspector, palette, story, filter,
  lasso, recorder, the widget host and its widgets (Birth tree, Relative SFH,
  Notes), actions and axes.

Per-object GPU state is a bitfield: 1 = dimmed, 2 = hidden (the filter),
4 = replaced by member stars in Sky view, 8 = dimmed and 16 = hidden by a
lasso, 32 = dimmed by the Birth tree. The shaders draw the lasso from a
separate per-object pair: the weight a fade starts from and the level in the
selection. That lets a State change crossfade between selections; bits 8 and
16 serve the CPU's readers.
