Overview
========

Oviz turns tables of Galactic objects into one interactive HTML file:

1. pandas and Astropy hold the tables and coordinates.
2. galpy integrates each object's orbit over a time grid.
3. :func:`oviz.viewer.compile_scene_spec` packs the frames, volumes and member
   stars into compact binary blobs inside the page.
4. A WebGL2 viewer replays them on the GPU, and Aladin Lite supplies the
   registered sky backgrounds.

Readers need only a browser. Authors can save a sequence of views and present
them from the same file.


Data model
----------

Use :class:`oviz.Trace` for objects whose orbits should be integrated. Each row
needs heliocentric Galactic Cartesian ``x``, ``y`` and ``z`` in pc (x towards
the Galactic centre), heliocentric ``U``, ``V`` and ``W`` in km/s, a unique
``name`` and ``age_myr``. Use :class:`oviz.Layer` for other spatial data; a
layer with only positions sets ``assume_stationary=True``.

:class:`oviz.TraceCollection` and :class:`oviz.LayerCollection` group them,
and :class:`oviz.Scene3D` (an alias of :class:`oviz.Animate3D`) builds the
frames. ``make_plot`` returns an :class:`oviz.OvizFigure`:

.. code-block:: python

   import numpy as np
   from oviz import Scene3D, Trace, TraceCollection

   trace = Trace(clusters, data_name="Young clusters", color="#ff5a5a")
   scene = Scene3D(TraceCollection([trace]), figure_theme="dark")
   figure = scene.make_plot(time=np.arange(0, -61, -1.0), galactic_mode=True, enable_sky_panel=True)
   figure.write_html("clusters.html")

The time grid is in Myr, must be evenly spaced and must include 0. Positions
are shown in a frame centred on the Local Standard of Rest, which follows a
circular orbit through the axisymmetric part of the potential.

Other ``make_plot`` inputs add to the same figure:

- ``volumes``: 3D dust, emission or density cubes, from a FITS file or an
  array with bounds in pc.
- ``cluster_members_file``: one row per member star (a cluster-name column and
  ``l``/``b`` or ``ra``/``dec``); in Sky view each cluster marker crossfades
  into its members.
- ``spiral_arm_models``: published spiral arms (below).
- ``annotations``: labels, curves, arrows, bubbles and shells drawn into the
  figure (:mod:`oviz.annotations`); readers can draw and edit them too (K).
- ``viewer_mode``, ``camera_anchor``, ``actions``: how the figure opens.

``viewer="classic"`` writes the previous Three.js runtime instead
(:class:`oviz.threejs_figure.ThreeJSFigure`), which still offers Slides and
Paper exports. Older figures convert to the current viewer with
:func:`oviz.upgrade_html`.


Potentials and spiral arms
--------------------------

Orbits use galpy's MWPotential2014 unless the scene is given another
``potential``. :func:`oviz.spiral_models.khalil2025_potential` adds the m = 2
and m = 3 spiral modes of Khalil et al. (2025) as galpy
``SpiralArmsPotential`` terms; in the Galactic plane they reproduce the
authors' released SPIBACK potential exactly. Their bar, their own axisymmetric
background and the modes' radial cutoffs are left out.

``make_plot(spiral_arm_models=...)`` draws published arms. Each model is one
line trace holding all of its arms, turning at its pattern speeds through the
timeline. The traces start hidden, and Sky view leaves them out.

.. code-block:: python

   from oviz.spiral_models import CASTRO_GINARD2021_ARMS, KHALIL2025_ARMS, khalil2025_potential

   scene = Scene3D(TraceCollection([trace]), figure_theme="dark", potential=khalil2025_potential())
   figure = scene.make_plot(
       time=time_myr,
       galactic_mode=True,
       spiral_arm_models=(KHALIL2025_ARMS, CASTRO_GINARD2021_ARMS),
   )

``KHALIL2025_ARMS`` are the troughs of the two modes (m = 2: Crux-Scutum, Local
and Outer arms at 13.1 km/s/kpc; m = 3: Carina-Sagittarius and Perseus at
16.4 km/s/kpc). ``CASTRO_GINARD2021_ARMS`` are the Perseus, Local, Sagittarius
and Scutum segments fitted to young open clusters and masers, each turning at
its own pattern speed (17.8, 33.8, 26.1 and 49.8 km/s/kpc).


Views and presentation
----------------------

A saved view (a *State*) captures the whole viewer: camera, time, mode, layer
and volume styles, Sky layers, filters, selections, widgets and display
settings. Views form a story that can be presented, exported as a
present-only copy, or shared as links. :doc:`viewer` describes the controls
and the browser API.


Output
------

The figure is a single file that can be attached, opened from disk, or put on
any static host such as GitHub Pages. The viewer runtime is inlined, so the
3D scene works offline; Aladin Lite and the HiPS sky surveys load from the
network. Large figures stay compact: the October 1 figure (about 1,700
clusters over 61 frames, three volumes and 590,000 member stars) is 26 MB.


Testing
-------

.. code-block:: bash

   pytest -q tests
   node --test tests/viewer_js/*.test.mjs

When the viewer changes, also open a regenerated figure in a browser.
