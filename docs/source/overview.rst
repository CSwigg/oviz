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
