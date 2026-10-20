# Dry-cooler model

Model of a finned-tube dry cooler (coolant in tubes, air across plate fins) with optional adiabatic precooling, three
solvers of increasing detail (LMTD < NTU element field < Cell method) and a comparison of standard cooler designs over a
Munich climate year. Written for a bachelor's thesis (TUM, Chair of Energy Systems); start with `docs/handover.md`.

## Run

```
venv/bin/pip install -r requirements.txt
venv/bin/python main.py solver              # one cooler, the selected solvers compared (default)
venv/bin/python main.py compare geometry    # the standard coolers P1-P3 raced against each other
venv/bin/streamlit run gui.py               # the same two workflows with forms
```

## Two workflows

| | Solver comparison | Cooler comparison |
|---|---|---|
| Question | how do the solvers differ for one cooler in one ambient state? | which cooler design operates best over a year? |
| Cooler | one standard cell `GEOMETRY_SET` P1/P2/P3 (every value overridable) in a module `N_TUBES` x `N_ROWS` x `HEIGHT_MM` | all `COMPARE_SETS`, scaled to the same water flow (more tubes, one tube length) |
| Solver | the selected ones side by side (`SOLVERS`: LMTD, NTU field / table, Cell 2D / 3D) | one, `COMPARE_SOLVER` |
| Ambient | `T_AIR_IN`, `PHI_AIR_PERCENT` | the Munich climate year (`climate_data/`), representative states |
| Extras | ambient + precooled scenario, solver benchmarks (caching, resolution, annual, accuracy) | geometry, fan control, water and cost, sensitivity, rows, tubes, wall thickness |

Both read `src/parameters.py` (the GUI builds its forms from it; section headers are tagged `[solver]`, `[compare]`
or `[both]`) and share the cooler cells (`src/cooler_sets.py`) and the climate data (`src/weather.py`).

## Layout

```
main.py, gui.py       entry points (command line / Streamlit; the GUI's parts are in ui/)
src/                  the model: geometry, heat transfer, solvers (src/solvers/), pressure drop, precooling, economics, weather,
                      the solver-comparison pipeline and its outputs, the solver benchmarks, plot style and figure helpers
comparison/           the cooler comparison (start with __init__.py and analyses.py)
thesis_figures/       chapter-3 schematics; the rendered figures go to out/ (gitignored)
climate_data/         derived hourly weather, DWD station München-Stadt (03379), see its README
docs/                 handover, workflows, cooler comparison, model notes, decisions, notation, GUI and inputs
tests/                `venv/bin/python -m pytest` (about 100 s)
```
