# Handover

A model of a finned-tube dry cooler (water in the tubes, air across plate fins, optional adiabatic precooling pad) for
the bachelor's thesis *Modeling and Parameter Identification of a Dry Cooler for Integration
with Absorption Systems* (TUM, Chair of Energy Systems; submitted October 2026). The thesis text is in the separate
LaTeX project; this repository holds the model, the figures and the studies behind its chapters 3 and 4.

## Start

```
python3 -m venv venv && venv/bin/pip install -r requirements.txt     # developed with Python 3.14
venv/bin/python -m pytest                                            # about 100 s (4 min the first time: the NTU tables are built)
venv/bin/python main.py solver                                       # one cooler, LMTD / NTU / Cell side by side (a few seconds)
venv/bin/python main.py compare pinch_fan                            # the three coolers over the climate year (needs the NTU tables, see below)
venv/bin/streamlit run gui.py                                        # both workflows with forms
```

Reading order: `README.md`, this file, `docs/workflows.md` (how to run and what comes out), `docs/model-notes.md`
(sources, assumptions, solver design), `docs/cooler-comparison.md`,
`docs/decisions.md`; `docs/notation.md` is the naming scheme (thesis symbols), `docs/gui-and-inputs.md` the GUI and the
units of the inputs. All inputs are in `src/parameters.py` (the GUI builds its forms from it).

## How a run works

`python main.py solver` (one cooler, the solvers side by side):

1. `main.py` reads the options and calls `src/scenario_pipeline.run_scenarios`.
2. The inputs come from `src/parameters.py`: `get_geometry()` builds the cooler (a preset cell P1-P3 from
   `src/cooler_sets.py` plus any
   override), `get_operating_conditions(geo)` the air and water state, `get_solver_settings()` the tolerances and
   resolutions.
3. `solvers.solve_scenario` runs every selected solver on the dry operating point. All three iterate the same loop:
   properties at the
   mean temperatures (`fluid_properties`) -> overall heat transfer coefficient k (`heat_transfer_core.calc_overall_k`)
   -> heat flow and outlet
   temperatures (LMTD formula, NTU element field or table, Cell grid) -> repeat until both outlets change by less than
   `CONVERGENCE_THRESHOLD`.
4. If the most detailed solver misses the target outlet temperature, `precooling.precooling_decider` switches the pad on
   and
   `precooling.calc_precooler` builds the precooled operating point, which is solved the same way.
5. Pressure drops (`pressure_drop`), fan and pump power and pad water (`economics`) are computed from the results.
6. The selected outputs (the catalogue in `src/outputs.py`) draw the figures (`plot_*.py`), saved through
   `figures.save_figure`.

`python main.py compare <analysis>` (the coolers over the climate year): `comparison/runner.run_analyses` runs the
chosen items of the
catalogue in `comparison/analyses.py`. Each one computes a table of results, one row per operating point
(`evaluate.run_point`: geometry and state ->
solver -> pressure drops -> powers -> one result row), stores it as CSV, and plots from the table.

## Use the model from Python

```python
from src import solvers
from src.dry_cooler_physics import get_geometry
from src.operating_conditions import get_operating_conditions
from src.pressure_drop import calc_delta_p_air
from src.economics import calc_fan_power
from src.solver_settings import get_solver_settings

geo = get_geometry()                                        # the cell GEOMETRY_SET of parameters.py (P2) and its overrides
ops = get_operating_conditions(geo, T_a_i=30.0, phi=0.4)    # any operating input can be overridden here
settings = get_solver_settings()
for key in ("lmtd", "ntu_table", "cell_2d"):
    r = solvers.run_solver(key, ops, geo, settings, profile=False)
    print(f"{key:10s} Q = {r.Q_dot / 1e3:6.2f} kW   T_c,o = {r.T_c_o:5.2f} °C   pinch = {r.T_c_o - ops.T_a_i:5.2f} K")
dp = calc_delta_p_air(geo, ops, r.T_a_o, through_pad=False)
print(dp, "Pa;", calc_fan_power(ops, r.T_a_o, dp), "W")
```

`python -c` from the repository root works; the first `ntu_table` call builds and stores its table (about 20 s at 20
elements). A cooler
that is not one of the presets is a `CellSet` (`src/cooler_sets.py`) or overrides in `parameters.py`.

## Glossary

| term | meaning |
|---|---|
| P1, P2, P3 | the three standard coolers from the literature (P1 dense: Wang & Chi; P2 medium and P3 open: Kayansayan) |
| elementary cell, module | the cell holds the coil's tube and fin dimensions, pitches and materials; a module is a cell with N_t tubes per row, N_r rows and a tube length |
| rows N_r | tube rows in the air direction; every tube is its own N_r-pass circuit, so the rows are also the coolant passes |
| target, pinch | the coolant outlet should reach 25 °C; pinch = T_c,o - the ambient air temperature entering the whole cooler |
| design point, P95 | the ambient state 20 °C / 60 %; the 95th-percentile warm hour of the climate year |
| dry, precooled, wet, hybrid | dry = air straight to the coil; precooled (wet, hybrid) = the adiabatic pad cools the air towards the wet-bulb temperature first |
| fixed fan, adaptive fan | face velocity fixed at 3 m/s, or per hour the smallest speed that still meets the target; `limited` = missed even at the top speed |
| solver, family, reference | LMTD, NTU and Cell are the families; `ntu_field` / `ntu_table` and `cell_2d` / `cell_3d` the variants; the most detailed selected solver is the reference for deviations |
| NTU field / table | the element field solved in every iteration, or a precomputed table P(NTU, R) interpolated |
| Cell 2D / 3D | cells with local properties; 3D resolves all tubes of a row with air exchange between neighbours, 2D a single tube per row (identical for uniform inlet air) |
| cache no / cold / warm | CoolProp lookups not cached / cache empty and filling / cache filled |
| catalogue, `Item` | a list of outputs or analyses, each a compute, a plot and a log step (`src/catalogue.py`) |

## For the follow-up project (cooling towers, once-through water cooling)

The next project models cooling towers and running-water cooling and compares their performance with the dry cooler. The
thermodynamics differ (evaporation into the air stream, a driving force set by the wet-bulb temperature, no finned-tube geometry), so
this repository is a starting point, not a base class. What carries over and what does not:

- **Carries over as it is**
  - `fluid_properties.py` (CoolProp water and humid air, with the lookup cache) and the moist-air functions in `precooling.py`
    (`calc_Y`, `calc_h`, the cooling limit = wet-bulb temperature, humidity ratio and enthalpy).
  - `weather.py` and `climate_data/` (the Munich year, the representative states, wet-bulb temperature).
  - `economics.py` (pump and fan power from a pressure drop, pad water) and the cost columns of `comparison/operation.add_costs`.
  - The comparison framework: `comparison/runner.py`, the catalogue in `analyses.py`, the `Settings` object, the year loop and the fan
    control in `operation.py`, and the print-size figures in `plots_print.py`; the figure conventions in `src/plot_style.py`.
- **Specific to the finned-tube dry cooler (to be replaced by the tower's counterparts)**: `dry_cooler_physics.py` (geometry),
  `cooler_sets.py` and `comparison/sizing.py`, `heat_transfer_core.py` (correlations, fin efficiency, k), the air-side pressure drop
  of `pressure_drop.py` (Wang), and `src/solvers/`. A tower needs its own transfer model (typically an enthalpy-based
  Merkel / effectiveness-NTU description of the fill) and its own pressure drop, and replaces the geometry; once-through cooling
  has no fan and needs a water-side heat exchanger and a source-water temperature series instead of the air climate.
- **The seam for a comparison**: `comparison/evaluate.run_point(cfg, geo, ops, method, mode)` turns one geometry and one ambient state into one
  result row. A tower model that returns a row with the same columns -- `T_c_o`, `Q_kW`, `pinch`, `P_fan`, `P_pump`, `P_total`, `m_dot_ev`
  (evaporation; add blowdown for a tower), `target_met`, plus its own geometry and mass columns for the cost model -- can be fed to the
  existing year loops, the fan control and the cost and figure code. The common basis of a comparison is the one used here: 1 kg/s of water
  from 37 °C to the 25 °C target over the Munich climate year; `pinch` is defined against the AMBIENT DRY-BULB temperature for every
  cooler type (a tower's approach to the wet-bulb temperature is a different number and should be reported separately).
- **Not yet general**: the comparison's cooler list (`config.sets`, `cooler_sets.py`) and `operation.control_point`'s face-velocity search
  are written for the finned-tube coolers; the cost model prices aluminium and copper write-off only.

## Status

- Solvers: LMTD, NTU (element field and precomputed table) and Cell (2D / 3D). The numbers of all solvers are pinned in
  `tests/golden_solvers.json`. Nothing is validated against measurements of a real cooler.
- Tests: `venv/bin/python -m pytest` runs every workflow and benchmark on a reduced problem (figures and CSVs go to a
  temporary folder) and pins what must not move silently: the tube counts 34 / 32 / 30, the
  independence of the geometry from the operating inputs, the fan and water numbers of the P2 27 °C case (386.66 W dry,
  734.58 W and 10.8 g/s precooled), the solver numbers and the profiling counts.
  A changed pinned number is a model change and needs a decision. Run the tests after every refactoring step.
- Open points and limitations: `docs/model-notes.md`, section 3. Not built: two-parameter maps and a global sensitivity;
  a construction-cost model; a decade of climate data (needs the DWD historical archives of
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

Commands (`MPLBACKEND=Agg` to run without windows). The solver outputs are `python main.py solver --solvers lmtd
ntu_field ntu_table cell_2d cell_3d --outputs <output>`
(`--solvers` and `--outputs` default to the selection in `parameters.py`); the analyses are `python main.py compare
<analysis> [--plots-only]`. Typical run times with all five solvers:
`cooler_results`, `profile_grids`, `convergence`, `accuracy` a few seconds each; `caching` about 40 s; `resolution`
about 6.5 min; `annual` minutes with the default stride (the whole year for every solver: about
105 min, see `docs/workflows.md`); `geometry` seconds; `pinch_fan` about 2 min; `variation` about 8 min plus the first
table builds. `--plots-only` redraws in seconds from the stored data.

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

Reproducing the thesis figures: the solver and cooler-comparison figures above regenerate the submitted ones (checked
pixel by pixel against the thesis project, apart from `fan_cop`, which now draws the total and the
pump-only COP in addition to the fan-only curves of the thesis). The raw results of the long runs are stored
(`annual_results.pkl`: the whole-year run, all five solvers; `resolution_results.pkl`), so their
figures need no new run. The caching benchmark times depend on the machine.

## Conventions

- Code: shared helpers instead of copies, one pattern for the same job. Both workflows offer their outputs as `Item`
  catalogues (`src/catalogue.py`; add an `Item` and the CLI and the GUI follow). Settings are
  objects (`comparison.config.Settings`, `SolverSettings`, `RunModes`), never module state read deep inside the model.
- Fixed geometry: the tube counts are the cells' `N_t_ref` (34 / 32 / 30), independent of every operating input.
- Notation: thesis symbols everywhere a user sees them (`docs/notation.md`). Colours: TUM Corporate Design
  (`src/plot_style.py`).
- Figures: print size 6.9 in wide at 300 dpi, no titles, styles as contexts (`docs/workflows.md`, Figure style).
- Logging: `logging.getLogger(__name__)`, no `print`.
