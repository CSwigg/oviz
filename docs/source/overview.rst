Overview
========

Oviz joins four parts of an astronomical visualization workflow:

1. pandas and Astropy prepare tables and coordinates.
2. galpy integrates phase-space measurements through a Galactic potential.
3. Three.js renders points, tracks, labels, images, and volumes in 3D.
4. Aladin Lite supplies registered HiPS backgrounds in Sky mode.

The output is a browser-ready HTML figure. The reader does not need Python, and
the author can save a sequence of complete viewer States to guide the reader
through the data.


Data model
----------

Use :class:`oviz.Trace` for an orbiting sample. Each row needs Galactic
Cartesian ``x``, ``y``, and ``z`` in parsecs; ``U``, ``V``, and ``W`` in km/s;
``name``; and ``age_myr``. Use :class:`oviz.Layer` for a more general spatial
layer. A layer with only XYZ coordinates must set ``assume_stationary=True``.

:class:`oviz.TraceCollection` and :class:`oviz.LayerCollection` combine these
objects. :class:`oviz.Animate3D` and its alias :class:`oviz.Scene3D` build the
time frames and return a :class:`oviz.threejs_figure.ThreeJSFigure`.


Potentials and spiral arms
--------------------------

Orbits use galpy's MWPotential2014 unless the scene is given another
``potential``. :func:`oviz.spiral_models.khalil2025_potential` adds the two
spiral modes of Khalil et al. (2025): MWPotential2014 plus their m = 2 and
m = 3 spirals as galpy ``SpiralArmsPotential`` terms. In the Galactic plane
these reproduce the authors' released SPIBACK spiral potential exactly. The
model leaves out their bar, their own axisymmetric background, and the modes'
radial cutoffs. The default LSR frame keeps its circular orbit in the
axisymmetric part of any potential.

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
and Scutum segments fitted to young open clusters and masers. Each segment
turns at its own pattern speed: 17.8, 33.8, 26.1 and 49.8 km/s/kpc.


Viewer model
------------

The maintained renderer is the standalone Three.js viewer. Its spatial modes
are 3D and Sky. Sky mode keeps Oviz-rendered points registered with Aladin Lite
background imagery. Optional cluster-member tables can replace bulk cluster
markers with member stars in Sky mode. Optional volume layers can show 3D dust,
emission, or other scalar fields.

The timeline interpolates between generated frames. Controls, trace styles,
selections, labels, panels, camera settings, and Sky layers can all be captured
in a State.


States and presentation
-----------------------

States are ordered, named snapshots of the complete viewer. The original scene
is the implicit State 0. A destination State controls its transition duration,
easing, and whether its saved camera is followed or the current camera is kept.

The States drawer supports capture, update, rename, duplicate, reorder, delete,
preview, and HTML export. Standard exports remain editable. Present-only exports
show simple navigation controls and move through the saved sequence.

The same system is available through the browser API described in
:doc:`browser_api`.


HTML export
-----------

Large scene specifications can be gzip-compressed inside the output HTML. The
figure is a single file that can be attached, copied to a static host, or placed
on GitHub Pages. Three.js, Aladin Lite, and remote HiPS surveys are fetched at
runtime, so an internet connection is required for those assets.


Testing
-------

Run the maintained tests with:

.. code-block:: bash

   pytest -q tests

When runtime code changes, regenerate and inspect the canonical HTML in a
browser in addition to running unit and regression tests.
