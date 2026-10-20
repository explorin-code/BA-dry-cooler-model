# Dry-cooler model (bachelor's thesis, TUM)

Finned-tube dry cooler with optional adiabatic precooling: three solvers (LMTD, NTU, Cell) and a comparison of three
standard coolers (P1-P3) over a Munich climate year. Entry points: `main.py solver|compare`, `gui.py`. The thesis text is in a
separate LaTeX project; thesis symbols are used throughout (`docs/notation.md`).

- Start with `docs/handover.md`; the other documents are indexed there.
- All inputs are in `src/parameters.py`; the tube counts are fixed per cooler (`src/inputs/cooler_sets.py`).
- Tests: `venv/bin/python -m pytest`. A changed pinned number is a model change.
