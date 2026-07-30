Python API
==========

Most figures need only the objects exported from :mod:`oviz`.


Scene builders
--------------

.. autoclass:: oviz.Animate3D
   :members: make_plot

.. autoclass:: oviz.Scene3D
   :members:


Traces and layers
-----------------

.. autoclass:: oviz.Trace
   :members: integrate_orbits, limit_age, limit_by_name, copy

.. autoclass:: oviz.Layer
   :members: from_dataframe, integrate_orbits, supports_orbit_tracing

.. autoclass:: oviz.TraceCollection
   :members: add_cluster, get_cluster, get_all_clusters, integrate_all_orbits,
             set_all_cluster_sizes, limit_all_cluster_ages,
             limit_all_cluster_names

.. autoclass:: oviz.LayerCollection
   :members: add_layer, get_layer, get_all_layers, integrate_all_layers,
             set_all_layer_sizes


Viewer profiles
---------------

.. automodule:: oviz.threejs_profiles
   :members:


HTML figure
-----------

.. autoclass:: oviz.threejs_figure.ThreeJSFigure
   :members: to_dict, to_html, write_html, show


States, actions, and authored content
-------------------------------------

.. automodule:: oviz.threejs_states
   :members:

.. automodule:: oviz.threejs_actions
   :members: normalize_threejs_actions

.. automodule:: oviz.threejs_deck
   :members: normalize_deck_object, normalize_deck_block, normalize_deck_spec

.. autoclass:: oviz.paper.Paper
   :members:


Orbit and appearance helpers
----------------------------

.. automodule:: oviz.orbit_maker
   :members:

.. automodule:: oviz.point_sizes
   :members:

.. automodule:: oviz.spiral_arms
   :members:


Lower-level scene export
------------------------

These functions are useful when an application already has an Oviz scene
specification. Most users should prefer :meth:`oviz.Animate3D.make_plot`.

.. automodule:: oviz.threejs_scene
   :members: build_threejs_scene_spec

.. automodule:: oviz.threejs_embed
   :members: render_threejs_html, threejs_data_url, threejs_iframe_html
