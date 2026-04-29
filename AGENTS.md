# Agent Notes

- The main figure runner belongs in `tests/main_figure.py`.
- The generated main figure artifact belongs in `tests/main_figure.html`.
- Keep `scripts/main_figure.py` as a compatibility wrapper that delegates to `tests/main_figure.py`.
- Generate the stored main figure with the compact payload enabled so `tests/main_figure.html` stays under 100 MB.
