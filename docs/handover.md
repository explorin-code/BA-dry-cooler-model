# Handover

A model of a finned-tube dry cooler (water in the tubes, air across plate fins, optional adiabatic precooling pad) for the bachelor's thesis *Modeling and Parameter Identification of a Dry Cooler for Integration
with Absorption Systems* (TUM, Chair of Energy Systems; submitted October 2026). The thesis text is in the separate LaTeX project; this repository holds the model, the figures and the studies behind its chapters 3 and 4.

## Start

```
python3 -m venv venv && venv/bin/pip install -r requirements.txt     # developed with Python 3.14
venv/bin/python -m pytest                                            # about 100 s (4 min the first time: the NTU tables are built)
venv/bin/python main.py solver                                       # one cooler, LMTD / NTU / Cell side by side (a few seconds)
venv/bin/python main.py compare pinch_fan                            # the three coolers over the climate year (needs the NTU tables, see below)
venv/bin/streamlit run gui.py                                        # both workflows with forms
```

Reading order: `README.md`, this file, `docs/workflows.md` (how to run and what comes out), `docs/model-notes.md` (sources, assumptions, solver design), `docs/cooler-comparison.md`,
`docs/decisions.md`; `docs/notation.md` is the naming scheme (thesis symbols), `docs/gui-and-inputs.md` the GUI and the units of the inputs. All inputs are in `src/parameters.py` (the GUI builds its forms from it).

## Status

- Solvers: LMTD, NTU (element field and precomputed table) and Cell (2D / 3D). The NTU element field reproduces the published tables of its source (`python -m src.validation_ntu`); the numbers of all solvers are pinned
  in `tests/golden_solvers.json`. Nothing is validated against measurements of a real cooler.
- Tests: `venv/bin/python -m pytest` runs every workflow and benchmark on a reduced problem (figures and CSVs go to a temporary folder) and pins what must not move silently: the tube counts 34 / 32 / 30, the
  independence of the geometry from the operating inputs, the fan and water numbers of the P2 27 °C case (386.66 W dry, 734.58 W and 10.8 g/s precooled), the solver numbers and the profiling counts.
  A changed pinned number is a model change and needs a decision. Run the tests after every refactoring step.
- Open points and limitations: `docs/model-notes.md`, section 3. Not built: two-parameter maps and a global sensitivity; a construction-cost model; a decade of climate data (needs the DWD historical archives of
  station 03379, see `climate_data/README.md`).

## Data and caches

| what | where | note |
|---|---|---|
| climate year | `climate_data/` (tracked) | DWD station 03379 München-Stadt, 2025-10-09 .. 2026-10-08; rebuilt from the raw download by `python -m src.build_climate_data` (raw folder `raw_climate_data/`, not tracked) |
| NTU tables | `cache/ntu_tables/` (ignored) | built on first use, once per row count: 110-230 s each at 60 elements; `NTU_TABLE_VERSION` forces a rebuild; the first `compare` run builds all that it needs |
| stored results | `output/comparison/*.csv`, `thesis_figures/out/ch4_results/solver/*.pkl` (ignored) | what `--plots-only` and `figures.load_data` redraw from; delete to force a recomputation |
| figures | `thesis_figures/out/` (ignored) | the working output of every run; the figures in the thesis project are the copies of the submission |
| defense slides | `defense/` (ignored) | the presentation files, outside version control |

## Figures: how each one is made

Commands (`MPLBACKEND=Agg` to run without windows). The solver outputs are `python main.py solver --solvers lmtd ntu_field ntu_table cell_2d cell_3d --outputs <output>`
(`--solvers` and `--outputs` default to the selection in `parameters.py`); the analyses are `python main.py compare <analysis> [--plots-only]`. Typical run times with all five solvers:
`cooler_results`, `profile_grids`, `convergence`, `accuracy` a few seconds each; `caching` about 40 s; `resolution` about 6.5 min; `annual` minutes with the default stride (the whole year for every solver: about
105 min, see `docs/workflows.md`); `geometry` seconds; `pinch_fan` about 2 min; `variation` about 8 min plus the first table builds. `--plots-only` redraws in seconds from the stored data.

| figure file | made by | thesis |
|---|---|---|
| `solver/cooler_results_dry`, `cooler_results_precooled` | `solver` output `cooler_results` | |
| `solver/convergence_dry`, `convergence_precooled` | output `convergence` | |
| `solver/profile_grid_absolute` | output `profile_grids` | |
| `solver/profile_grid_normalised` | output `profile_grids` | Results |
| `solver/profile_collapsed_normalised`, `profile_results_table` | output `profile_grids` | Appendix |
| `solver/accuracy_sweep` | output `accuracy` | Results |
| `solver/accuracy_sweep_pinch` | output `accuracy` | |
| `solver/caching_times`, `caching_split` | output `caching` | Results |
| `solver/caching_speedup`, `caching_presentation_0 ... _4`, `caching_presentation_with3d_0 ... _5`, `caching_counts.tex` | output `caching` | |
| `solver/resolution_performance` | output `resolution` | Results |
| `solver/annual_deviation`, `annual_times` | output `annual` | Results |
| `geometry/geometry_print` | `compare geometry` | Results |
| `geometry/module_geometries`, `bank_full`, `geometries` | `compare geometry` | |
| `cooler_performance/pinch_year`, `missed_hours`, `controlled_fan`, `fan_cop`, `annual_costs` | `compare pinch_fan` | Results |
| `cooler_performance/pinch_year_clouds`, `monthly_cost` | `compare pinch_fan` | Appendix |
| `cooler_performance/pinch_speed`, `water_vs_hours`, `annual_costs_appendix.tex` | `compare pinch_fan` | |
| `cooler_performance/variation_costs` | `compare variation` | Results |
| `cooler_performance/variation_totals` | `compare variation` | Appendix |
| `cooler_performance/ambient_cases`, `airflow_sweep` | `compare ambient_cases`, `airflow_sweep` | |
| `cooler_performance/control_year`, `fan_year`, `year_scatter`, `fan_curves`, `cost_year`, `annual_energy` | `compare fan_control` | |
| `cooler_performance/sensitivity_lines`, `sensitivity_tornado` | `compare sensitivity` | |
| `cooler_performance/rows_analysis`, `tubes_analysis`, `wall_analysis` | `compare rows`, `tubes`, `wall` | |
| `ch3_methodology/*` | `python -m thesis_figures.chapter3_schematics` | Chapter 3 (as `Plot_*.png`) |

Reproducing the thesis figures: the solver and cooler-comparison figures above regenerate the submitted ones (checked pixel by pixel against the thesis project, apart from `fan_cop`, which now draws the total and the
pump-only COP in addition to the fan-only curves of the thesis). The raw results of the long runs are stored (`annual_results.pkl`: the whole-year run, all five solvers; `resolution_results.pkl`), so their
figures need no new run. The caching benchmark times depend on the machine.

## Conventions

- Code: shared helpers instead of copies, one pattern for the same job. Both workflows offer their outputs as `Item` catalogues (`src/catalogue.py`; add an `Item` and the CLI and the GUI follow). Settings are
  objects (`comparison.config.Settings`, `SolverSettings`, `RunModes`), never module state read deep inside the model.
- Fixed geometry: the tube counts are the cells' `N_t_ref` (34 / 32 / 30), independent of every operating input.
- Notation: thesis symbols everywhere a user sees them (`docs/notation.md`). Colours: TUM Corporate Design (`src/plot_style.py`).
- Figures: print size 6.9 in wide at 300 dpi, no titles, styles as contexts (`docs/workflows.md`, Figure style).
- Logging: `logging.getLogger(__name__)`, no `print`.
