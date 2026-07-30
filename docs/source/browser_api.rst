Browser API
===========

Every initialized viewer registers itself under ``window.Oviz``. Give
``window.Oviz.get`` the root element ID written into the HTML:

.. code-block:: javascript

   const root = document.querySelector('[id^="oviz-three-"]');
   const viewer = window.Oviz.get(root.id);
   const states = viewer.states;

Wait for the root's ``states-ready`` event before calling the API when a page
has just loaded.


State navigation
----------------

``list()``
   Return ordered State summaries. Indices are one-based.

``goTo(idOrIndex)``
   Enter a State by stable ID or one-based index. ``0`` and ``"original"``
   select the implicit original scene. Returns a Promise.

``next()``, ``previous()``, ``original()``
   Navigate the sequence. Navigation methods return Promises and resolve with a
   boundary result when no further State exists.


State authoring
---------------

Authoring methods require Edit mode. ``add`` and ``update`` are unavailable
during an active transition.

``capture()``
   Return a complete snapshot of the currently rendered viewer.

``add(options = {})`` and ``quickAdd(options = {})``
   Capture a State. Options may include ``name``, ``id``, ``transition``,
   ``camera_behavior``, and a supplied ``snapshot``.

``update(idOrIndex, options = {})``
   Replace a State snapshot and optionally update its name, transition, or
   camera behavior.

``rename(idOrIndex, name)``
   Rename a State.

``duplicate(idOrIndex, options = {})``
   Copy a State and insert the copy after it.

``move(idOrIndex, destinationIndex)``
   Move a State to a one-based destination index.

``remove(idOrIndex)``
   Delete a State.

``setCameraBehavior(idOrIndex, "follow" | "keep")``
   Choose whether entering this State applies its camera or preserves the live
   incoming camera.

``setMode("edit" | "present")``
   Change the States drawer mode.

``setPreserveCamera(enabled)``
   Set the project-wide camera-preservation preference.


Export
------

``exportHtml(options = {})``
   Export an editable State-enabled HTML. Pass ``download: false`` to receive
   the HTML string instead of starting a download. ``filename`` sets the
   suggested file name.

``exportPresentOnlyHtml(options = {})``
   Export a presentation that starts at the first State and exposes only the
   presentation navigation surface.


Events
------

The viewer root dispatches ``CustomEvent`` objects. Read payloads from
``event.detail``.

State lifecycle events include ``states-ready``, ``states-changed``,
``states-preload-complete``, ``transition-start``, ``transition-progress``,
``transition-complete``, and ``transition-error``.

Action lifecycle events include ``action-start``, ``action-step-start``,
``action-progress``, ``action-complete``, ``action-cancel``, and
``action-error``.


Parent-frame commands
---------------------

An embedding page can call any States method with ``postMessage``:

.. code-block:: javascript

   iframe.contentWindow.postMessage({
     source: "oviz-command",
     rootId: "oviz-three-...",
     requestId: "request-1",
     command: "goTo",
     args: [2]
   }, "*");

The viewer replies to the sender with ``source: "oviz-result"``, the same
``requestId``, an ``ok`` boolean, and either ``result`` or ``error``. In a
production host, validate message origins rather than accepting every origin.


Slide authoring
---------------

When a figure includes the optional slide editor, ``viewer.deck`` exposes slide
and object operations. The main methods are ``list``, ``current``, ``goTo``,
``add``, ``rename``, ``duplicate``, ``move``, ``remove``, ``undo``, ``redo``,
``present``, ``exit``, and ``exportHtml``. ``viewer.deck.objects`` manages text
and shapes; ``viewer.deck.guides`` reads and updates alignment guides.
