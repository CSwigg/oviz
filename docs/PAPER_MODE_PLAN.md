# Oviz Paper Mode — Design Plan

Status: draft for review (feature/paper-mode branch)

## 1. What we are building

A reader opens a single exported Oviz HTML file and sees the full 3D figure with a
scrollable panel on the right containing a complete academic paper — abstract,
sections, equations, figures, references. As they read, the scene behind the text
flies through authored **States**: when the text discusses the spiral arms, the
camera is looking at the spiral arms; when it discusses cluster kinematics at
−20 Myr, the timeline has scrubbed there; when it discusses the dust map, the
volume is framed and lit. The paper's plots appear inline in the text, and the
"main figure" of the paper is effectively live behind it.

One file, no server: the paper text, its figures, and the bindings all embed in
the exported HTML exactly like states and deck slides do today.

## 2. Why this fits the existing architecture

Almost every hard piece already exists:

| Need | Existing machinery |
| --- | --- |
| Capture/restore full viewer state (camera, time, legend, volume, sky/3D mode) | States system: `states.items[].snapshot` (`oviz/threejs_states.py`, `threejs_runtime_states.py`) |
| Smooth authored transitions between states | `ovizGoToState(...)` + transition runner (`duration_ms`, easing, tangent-space fov interpolation) |
| Serialized navigation that never double-fires during a transition | `ovizQueuePresentationStateNavigation` queue with generation-based cancellation (`threejs_runtime_states.py:2545`) |
| Content bound to a state | Deck already binds `slide.state_id -> state` (`oviz/threejs_deck.py`) |
| Interactive state authoring | In-viewer States drawer: fly camera, Save State, name it, export self-contained HTML |
| Embedding large binary content | gzip+base64 payload pipeline + `deduplicate_state_assets` content-addressed data-URL dedup (`threejs_states.py:98`) |
| Panel UI patterns, theming, Zen mode | Shell panels (legend, controls drawers), CSS variable theme, `data-zen` gating |

Paper mode is therefore mostly **one new spec section + one new runtime module +
one authoring API**, not a new engine.

## 3. UX design

### Layout
- Right-side panel, default ~42% of viewport width (author-configurable
  `panel.width_fraction`, clamped 0.25–0.6), full height, translucent dark
  panel consistent with existing `--oviz-panel-*` theming. Scene canvas stays
  fullscreen behind; the left ~58% remains fully interactive 3D.
- Topbar gets a `Paper ▸` toggle (same pattern as `Text ▸` / `Slides ▸`).
  Zen mode hides the panel. `data-paper-open` attribute on the root for CSS.
- Panel chrome: paper title + authors header, thin **progress rail** on the
  panel's left edge showing one tick per anchor (clickable), reading progress
  indicator, close button.
- Mobile (phase 2): full-screen sheet with a floating "show figure" toggle
  rather than a split view.

### Scroll → state sync
- Authors place **anchors** on blocks (a paragraph, heading, or figure). Each
  anchor binds to a state by *name* (resolved to state id at build time) plus an
  optional per-anchor transition override (duration/easing).
- A scrollspy watches which anchor most recently crossed the **reading line**
  (~35% from the panel top, IntersectionObserver + fallback scroll handler).
  Crossing an anchor queues a state jump. Scrolling backwards re-activates
  earlier anchors symmetrically — the sync is bidirectional and idempotent.
- **Last-wins coalescing**: fast scrolling past five anchors must not play five
  transitions. New runtime primitive `ovizQueuePresentationStateJump(idOrIndex)`
  reuses the existing navigation queue's await-active-transition logic but
  collapses the pending queue to the latest requested target (the existing
  next/previous queue runs every queued step; the paper needs "go directly to
  wherever the reader is now").
- Anchored paragraphs get a subtle margin marker (dot/►). The active anchor is
  highlighted; clicking any anchor marker applies its state immediately and
  scrolls it to the reading line.

### Reader control and trust
- The scene stays interactive while reading. If the reader grabs the camera
  (canvas pointerdown/wheel) while sync is on, sync **pauses** automatically —
  we never fight the reader for the camera. A small pill in the panel footer
  shows `Sync: on / paused — Resume`. Crossing the *next* anchor (or clicking
  Resume / any anchor marker) resumes sync. Policy is `sync.resume_policy:
  "next-anchor" | "manual"` (default `next-anchor`).
- Keyboard: PageUp/PageDown/space scroll the paper when the panel has focus.
  The existing global arrow-key state navigation stays available and simply
  scrolls the paper to the matching anchor when it lands on an anchored state
  (keeping text and scene consistent in both directions).
- Wheel/touch events inside the panel must never zoom the scene: the panel is a
  sibling of the canvas, so canvas handlers don't see them, but the root
  capture-phase invalidation listeners do — that is harmless (it only triggers
  renders). No special casing needed beyond `overscroll-behavior: contain`.

### Content presentation
- Blocks render sanitized author HTML: headings, paragraphs, emphasis, lists,
  links (open in new tab), block quotes, tables, footnotes.
- **Figures**: embedded as data-URL assets with captions ("Figure 3: ...").
  Click opens a lightbox. Optional `live: true` figure blocks render a framed
  placeholder that says "shown live in the scene" with a `Show` affordance that
  applies the bound state — for the case where the oviz view *is* the figure.
- **Math**: KaTeX auto-render loaded from CDN (same policy as three.js, Aladin,
  reveal.js — this export already requires network at runtime). Graceful
  degradation: if the CDN fails, raw `$...$` TeX remains readable. Delimiters
  `$...$` and `$$...$$`.
- **Citations** (phase 2): `[@key]` markers become popover references from an
  embedded bibliography list.

## 4. Data model

New top-level scene-spec section, normalized by a new `oviz/threejs_paper.py`
(mirroring `threejs_states.py` / `threejs_deck.py` conventions):

```jsonc
"paper": {
  "schema_version": 1,
  "available": true,          // panel + toggle exist in this export
  "enabled": true,            // panel open on load
  "title": "…",
  "authors": ["…"],
  "panel": { "width_fraction": 0.42 },
  "sync": {
    "mode": "scroll",         // "scroll" | "manual"
    "reading_line": 0.35,
    "resume_policy": "next-anchor"
  },
  "sections": [
    {
      "id": "sec-results",
      "level": 1,             // 1=section, 2=subsection…
      "title_html": "4. Results",
      "blocks": [
        { "id": "b1", "type": "html", "html": "<p>…</p>",
          "anchor": { "state": "Spiral arms edge-on",   // name, resolved to state_id at build
                      "state_id": "state-…",             // filled by the resolver
                      "transition": { "duration_ms": 1800, "easing": "easeInOutCubic" },
                      "label": "Fig. set: arm geometry" } },
        { "id": "b2", "type": "figure", "asset": "<sha256>",
          "caption_html": "Figure 3: …", "width_px": 1400, "live": false },
        { "id": "b3", "type": "html", "html": "<p>…</p>" }
      ]
    }
  ],
  "assets": { "<sha256>": "data:image/jpeg;base64,…" },  // deduplicate_state_assets reuse
  "bibliography": []           // phase 2
}
```

Anchors are a property of blocks (not standalone) so the reading line maps
cleanly to DOM elements. `state_id: null` anchors are allowed (marker with no
scene change — e.g., to force "return to overview").

## 5. Authoring pipeline

The intended workflow keeps state authoring **visual** and text authoring
**textual**:

1. Build the figure as usual; open it; fly the camera / set time / volume /
   sky view; save named states with the existing States drawer; use the
   in-viewer Save State export (or author states in Python, as
   `tests/main_figure_chronos_july4_reveal_authoring.py` demonstrates).
2. Write or convert the paper to Markdown (pandoc handles LaTeX → Markdown;
   phase 2 adds a direct ar5iv/LaTeX importer). Annotate bindings inline:

   ```markdown
   ## Results

   The clusters trace two kinematically distinct arms…
   <!-- oviz-state: "Spiral arms edge-on" duration=1800 -->

   ![Figure 3: Radial velocity vs. age.](figures/fig3.png)
   ```

3. Python API (new `oviz/paper.py`):

   ```python
   from oviz.paper import Paper

   paper = Paper.from_markdown("paper.md", figures_dir="figures/",
                               title="…", authors=[…],
                               max_figure_width_px=1600, jpeg_quality=0.85)
   scene["paper"] = paper.to_spec(states=scene["states"])   # resolves names → ids,
                                                            # errors on unknown states
   ```

   - Markdown conversion: small internal converter for the required subset
     (headings, paragraphs, emphasis, lists, links, images, code, blockquote,
     `$math$` passthrough), with `Paper.from_html(...)` as the escape hatch for
     full control. Avoids a hard new dependency; `markdown-it-py` can be an
     optional accelerator later.
   - Figures: lazy `PIL` import (same pattern as the volume PNG-atlas encoder)
     to downscale/re-encode to JPEG/PNG data URLs; dedup via
     `deduplicate_state_assets`.
   - Binding resolution by state **name** (exact, then case-insensitive), with
     a clear error listing available names on miss.

## 6. Runtime

New `oviz/threejs_runtime_paper.py` → `THREEJS_PAPER_RUNTIME_JS`, injected like
the other modules; one `<aside class="oviz-three-paper-panel">` container added
to the shell markup; `Paper ▸` topbar button (emitted only when
`paper.available`, mirroring the deck button pattern).

Responsibilities:
1. Render sections/blocks into the panel (author-trusted HTML; assets resolved
   from `paper.assets`).
2. KaTeX loader (CDN, deferred, optional) + auto-render over the panel only.
3. Scrollspy → `ovizQueuePresentationStateJump(stateId)` (new, in the states
   runtime, ~30 lines: last-wins queue on top of
   `ovizWaitForStatesControllerReady` + `ovizGoToState`).
4. Sync pause/resume: listen for canvas `pointerdown`/`wheel` while sync is on
   → paused; resume per policy. Root dataset flags for tests:
   `data-paper-open`, `data-paper-sync`, `data-paper-active-anchor`.
5. Progress rail + active-anchor highlight + click-to-jump.
6. Lightbox for figure blocks.

Non-goals for the runtime: no iframe, no reveal.js dependency (the deck keeps
reveal; the paper panel is plain scrolling DOM — simpler, accessible,
print-friendly).

### Verified integration constraints (from the code, exact anchors)

These are the sharp edges the implementation must handle explicitly:

- **State application entry**: `ovizGoToState(idOrIndex)`
  (`threejs_runtime_states.py:2464`) accepts a state id string, 1-based index,
  or `"original"`, and returns a promise. All runtime modules share one IIFE
  scope, so the paper runtime calls it directly. The serialized presentation
  queue (`ovizQueuePresentationStateNavigation`, `:2545`) awaits
  `ovizWaitForStatesControllerReady()` and the in-flight
  `ovizStateTransition.promise` — the new jump primitive copies that skeleton
  with last-wins coalescing.
- **Pointerdown cancels transitions**: `initializeOvizStates` installs a
  capturing root `pointerdown` listener that cancels any in-flight transition
  unless the event is inside `.oviz-states-shell`
  (`threejs_runtime_states.py:3097-3101`). The paper panel must be added to
  that exemption, or touching its scrollbar aborts scene transitions.
- **Keyboard falls through**: `onKeyDown` only ignores editable targets via
  `keyboardTargetIsEditable` (`threejs_runtime_viewer.py:2498-2510`). A
  scrollable div is not covered — ArrowUp/Down would drive WASD-family motion
  and ArrowLeft/Right time-stepping while the reader scrolls. Add the paper
  panel to that predicate (and `stopPropagation` on the panel like the states
  drawer at `threejs_runtime_states.py:2942` and deck editor at
  `threejs_runtime_deck.py:2282`).
- **Wheel is safe by construction**: zoom is bound to the canvas element only
  (`canvas.addEventListener("wheel", ...)`, OrbitControls on
  `renderer.domElement`); a sibling panel scrolls without zooming the scene.
- **Scene recentering**: the canvas is never narrowed; the correct mechanism
  for "panel covers the right 42%, center the subject in the visible left" is
  `camera.view_offset`, which is already captured per state, interpolated
  during transitions, and reapplied on resize
  (`applyActionCameraViewOffset`, `threejs_runtime_actions.py:361-382`).
  Paper mode applies a panel-width-derived offset while open
  (`panel.apply_view_offset`, default on) and removes it when closed.
- **Presentation-mode CSS hides all chrome** via an opt-out `display:none`
  list (`threejs_figure.py:6613-6631`); the paper panel needs an explicit
  exemption there only if it should coexist with presentation mode — v1
  treats paper mode and presentation mode as mutually exclusive drivers.
- **Naming collision**: a theme preset named `paper` already exists (light
  theme, `threejs_figure.py:133`). The spec section keeps the name `"paper"`
  (different namespace), but UI copy should avoid bare "Paper theme"
  ambiguity.
- **Wiring points** (five-touchpoint pattern, mirroring the deck): normalizer
  import + call in `ThreeJSFigure.__init__` (`threejs_figure.py:22550-22554`),
  runtime concatenated beside `THREEJS_DECK_RUNTIME_JS` (`:22598`) or via a
  new `__PAPER_RUNTIME_JS__` placeholder in `render_threejs_html`
  (`threejs_embed.py:183-232`), panel DOM in `threejs_shell.py`, init chained
  after `initializeOvizDeck()` (`threejs_figure.py:22462-22466`).
- **Round-trip**: unknown top-level spec keys survive the in-viewer
  Save-State/export path unchanged, so paper content persists through
  re-exports without extra work; mirroring into `ovizExportStatesHtml` is
  only needed if papers become editable in-viewer (phase 2).
- **Public API**: `window.Oviz.get(rootId).states` and the
  `oviz-command`/`oviz-result` postMessage channel already exist
  (`threejs_runtime_states.py:3034-3079`) — paper mode should register a
  small `viewer.paper = {open, close, goToAnchor, syncEnabled}` alongside,
  which also gives embedding pages control for free.
- **CDN caution**: the deck hard-fails presentation when the reveal CDN is
  unreachable (`threejs_runtime_deck.py:2222-2232`). KaTeX must do the
  opposite: fail soft, leave raw TeX readable, never block the panel.

## 7. Feature breakdown (implementation order)

Each item is small, testable, and lands green:

1. **Schema + normalizer** — `oviz/threejs_paper.py` (`normalize_paper_spec`),
   wired into `ThreeJSFigure.__init__` beside states/deck normalization.
   Unit tests mirroring `test_threejs_states`-style coverage.
2. **Shell + toggle** — panel container, CSS (theme vars, scrollbar, widths,
   `data-zen`/`data-paper-open` rules), topbar button gated on
   `paper.available`. String assertions in renderer tests.
3. **Panel renderer** — sections/blocks/figures/lightbox DOM builder; assets
   resolution; KaTeX loader. `node --check` + renderer string tests.
4. **State jump primitive** — `ovizQueuePresentationStateJump` in the states
   runtime + diagnostics dataset; reuse in paper runtime. Tests alongside the
   existing presentation-queue tests.
5. **Scrollspy + sync lifecycle** — reading line, bidirectional activation,
   pause-on-interaction, resume policy, progress rail.
6. **Authoring API** — `oviz/paper.py` (`Paper.from_markdown/from_html`,
   figure ingestion, name→id binding resolution, `to_spec`). Pure-Python tests.
7. **Demo artifact** — `tests/main_figure_paper_demo.py`: july25 scene + a
   3-section sample paper with 4–5 anchored states + 2 figures →
   `tests/main_figure_paper_demo.html`, plus its artifact test (size budget,
   node --check, schema assertions). This is the file we iterate on visually.
8. **Docs** — README section + this plan kept current.

Phase 2 backlog (explicitly out of scope for the first cut): mobile sheet
layout, deep links (`#sec-results` URL hash ↔ anchor + state), LaTeX/ar5iv
importer, bibliography popovers, partial-state anchors (e.g. camera-only,
keep the reader's legend toggles), per-anchor "hold" (disable sync between two
anchors), reading-position persistence, print stylesheet.

## 8. Risks / open questions

- **Panel width vs. scene composition**: authors may want states framed for
  the visible left region. Mitigation: apply `camera.view_offset` while the
  panel is open (the runtime already supports view offsets) so state cameras
  center on the *visible* area; make it a `panel.apply_view_offset` option
  (default on). Needs a little experimentation in the demo.
- **Source format**: plan assumes Markdown-first authoring (pandoc from
  LaTeX). If the primary source is arXiv HTML/LaTeX, the importer moves up
  the priority list. → confirm preferred input.
- **CDN math**: KaTeX from CDN matches existing network posture, but papers
  read offline would show raw TeX. Acceptable? (Pre-rendering at build time is
  possible later via `katex` npm, at the cost of a node build dependency.)
- **HTML size**: text is negligible; figures dominate. Downscaling defaults
  (1600 px, JPEG q0.85) keep a 10-figure paper ≈ 1–3 MB on top of the scene
  payload. Assets dedup already handles repeats.
- **Two content systems** (deck + paper): they stay independent; both bind to
  states. A single export should enable at most one of deck/paper by default
  to avoid competing state drivers (normalizer warning if both enabled).
