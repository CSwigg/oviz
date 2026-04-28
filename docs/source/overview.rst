Overview
========

Architecture
------------

The maintained `oviz` workflow is:

1. build one or more :class:`oviz.Trace` objects from cluster/member tables
2. collect them in a :class:`oviz.TraceCollection`
3. pass that collection into :class:`oviz.Animate3D`
4. export either a Three.js viewer or a legacy Plotly animation

Three.js is now the main interactive engine. Plotly output is still supported,
but primarily for backward compatibility with existing notebooks and older figure
generation paths.


Key modules
-----------

- ``oviz.traces``: data containers, validation, and filtering
- ``oviz.orbit_maker``: orbit integration and reference-frame handling
- ``oviz.viz``: high-level visualization orchestration and scene generation
- ``oviz.threejs_*``: Three.js HTML/runtime export pipeline
- ``oviz.app``: optional Dash wrapper with sky-panel integration


Testing
-------

Run the repository test suite with:

.. code-block:: bash

   pytest -q tests

The current tests exercise the maintained trace filtering, Three.js export, and
sky-panel helper paths.
