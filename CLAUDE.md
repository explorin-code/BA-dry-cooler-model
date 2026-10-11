# CLAUDE.md

Project context and instructions for Claude Code sessions in this repo. The details live in `docs/` (index below); this
file is what a session must know before touching anything. **Keep it and `docs/` up to date when something changes.**

## What this is

A model of a finned-tube dry cooler (water in the tubes, air across plate fins, optional adiabatic precooling pad) for a
bachelor's thesis, with two workflows (`README.md` has the table): the **solver comparison** (one cooler, LMTD vs. NTU
element field / table vs. Cell, plus solver benchmarks) and the **cooler comparison** (the standard coolers P1-P3 over a
Munich climate year with one solver). Entry points: `main.py solver|compare`, `gui.py` (Streamlit, same two modes).

## Standing rules (from the user)

- Keep the code clean and simple, no "vibe-coded mess": shared helpers instead of copies, one pattern for the same job.
- Thesis notation everywhere the user sees it (LaTeX symbols in GUI labels and plots, symbols in result tables); naming
  scheme in `docs/notation.md`.
- Never edit the thesis LaTeX files: give blocks to copy and paste.
- The colours are the TUM Corporate Design (`src/plot_style.py`, `TUM_*`).
- Answer questions about results with numbers from the current code, not from memory; say what was NOT run.
- AI-written formulas were not checked against the sources: `docs/review-needed.md` lists what is verified and what is not.
  Do not present unreviewed numbers as confirmed.

## Commands

```
venv/bin/pip install -r requirements.txt
venv/bin/python main.py solver [--solvers lmtd ntu_table ...] [--outputs cooler_results caching ...] [--cooler P3]
venv/bin/python main.py compare geometry|fan_control|... [--stride 24] [--solver lmtd]     # python -m comparison does the same
venv/bin/streamlit run gui.py                      # restart it after editing anything under src/ or comparison/
venv/bin/python -m pytest                          # tests (~40 s), see below
venv/bin/python -m src.build_climate_data          # rebuilds climate_data/ from the raw DWD download (raw_climate_data/, optional)
```

## Layout

```
main.py, gui.py        entry points (gui.py is the Streamlit page; its parts are in ui/: parameters.py reader, labels, presets, widgets, runner, results)
src/                   the model + the solver comparison
  parameters.py          ALL inputs (the GUI builds its forms from it: section headers `# --- Name [solver|compare|both]: ... ---`)
  cooler_sets.py         the standard cells P1-P3 (with fixed reference tube counts N_t_ref)
  dry_cooler_physics.py  geometry; heat_transfer_core.py correlations; pressure_drop.py every pressure drop; economics.py powers + water
  fluid_properties.py    CoolProp lookups (cached); precooling.py the pad; weather.py the climate year (from climate_data/)
  solvers/               package: LMTD, NTU (field / table), Cell (2D / 3D), registry (SOLVER_LABELS, run_solver), see its __init__ docstring
  scenario_pipeline.py   solver comparison: dry + precooled scenario, then the selected outputs (outputs.py, catalogue.py)
  benchmark_solvers.py, benchmark_annual.py   the benchmarks; profiling.py their instrumentation
  figures.py, quiet.py, formatting.py         shared helpers
comparison/            the cooler comparison (read __init__.py, then analyses.py = the catalogue; config.py = the Settings object)
tests/                 smoke tests    ui/   the GUI's parts    docs/   the details    climate_data/   derived weather    thesis_figures/   schematics + out/ (gitignored)
```

## Patterns to reuse

- **Catalogue + runner** (`src/catalogue.py`): both workflows offer their outputs / analyses as `Item(key, title, compute, plot, log)` lists
  (`src/outputs.py`, `comparison/analyses.py`) run by `run_items`. To add one: append an `Item`; the selection lists
  `SOLVER_OUTPUTS` / `COMPARE_ANALYSES`, the CLI and the GUI follow.
- **Settings objects, not module state**: the comparison reads `parameters.py` once per run into `config.Settings` (`cfg`, first argument of
  everything); the solvers get `SolverSettings`, the run `RunModes`. Never read `parameters` deep inside the model.
- **Fixed geometry**: the tube counts are the cells' `N_t_ref` (34 / 32 / 30), independent of every operating input.
- **Figures**: `figures.save_figure(fig, folder, name, title)`; solver figures go to `thesis_figures/out/solver/`, comparison figures to `out/ch4/`.

## Tests

`venv/bin/python -m pytest` (about 40 s, all of it): `tests/test_smoke.py` runs every workflow and benchmark on a reduced problem (figure and CSV output redirected to a temp folder) and
pins what must not move silently: the tube counts 34 / 32 / 30 (equal to `comparison.sizing.design_tubes`), the independence of the geometry from
the operating inputs, the fan / water numbers of the P2 27 degC case (386.66 W dry, 734.58 W / 10.8 g/s precooled), the profiling counts.
Run it after every refactoring step; a changed pinned number is a model change and needs a decision.

## The docs

| file | content |
|---|---|
| `docs/workflows.md` | the two workflows, run modes, the solver outputs and benchmarks (caching, resolution, annual), shared helpers |
| `docs/cooler-comparison.md` | the comparison: weather, cells, sizing, fan control, cost, results so far |
| `docs/gui-and-inputs.md` | inputs in mm / %, the GUI (layout, flow selectors), benchmark instrumentation |
| `docs/notation.md` | naming scheme (thesis symbols), open notation decisions |
| `docs/decisions.md` | decision log (thesis chapter 2 review) and the 2026-10 restructuring changes |
| `docs/review-needed.md` | what was AI-written and is NOT independently verified, per module; the solver notes (NTU, Cell, LMTD, convergence) |
| `docs/todo.md` | open points and known placeholders (hardcoded 0.8 pad efficiency, fixed pad pressure drop, ...) |
