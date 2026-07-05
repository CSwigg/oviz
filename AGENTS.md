# Agent Notes

- The main figure runner belongs in `tests/main_figure.py`.
- The generated main figure artifact belongs in `tests/main_figure.html`.
- Keep `scripts/main_figure.py` as a compatibility wrapper that delegates to `tests/main_figure.py`.
- Generate the stored main figure with the compact payload enabled so `tests/main_figure.html` stays under 100 MB.
- When the user asks to upload an Oviz figure, copy the produced HTML into
  `/Users/swiggumc/Desktop/astro_research/cam_website/oviz_figures`, then in
  `/Users/swiggumc/Desktop/astro_research/cam_website` git add, commit, and
  push the copied figure so it is available from the GitHub website.
- Before uploading Oviz HTML, check the artifact size and prefer compact or
  mobile-safe exports for large scenes. Uploaded figures should keep the
  desktop layout on laptop/desktop browsers while auto-detecting iPhone/mobile
  browsers and switching to mobile controls.
