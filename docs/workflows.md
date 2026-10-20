# The two workflows, their outputs and benchmarks

Everything runs in one of two workflows: `python main.py solver|compare` (default `solver`) and the two modes of the GUI
(radio button
at the top of the sidebar; own parameter groups and run button each).

- **Solver comparison** (`src/`): ONE cooler in one ambient state, the SELECTED solvers (`SOLVERS`: `lmtd`, `ntu_field`,
  `ntu_table`, `cell_2d`,
  `cell_3d`; the most detailed one is the reference for deviations and for the precooling decision), an ambient and a
  precooled
  scenario, and the solver benchmarks. No prices. The cooler is a standard cell (`GEOMETRY_SET` P1 / P2 / P3,
  `src/cooler_sets.py`); every cell value
  can be overridden (`None` = preset); module = `N_TUBES` x `N_ROWS` x `HEIGHT_MM`. Defaults: 1 kg/s, 2000 mm tubes, 6
  rows,
  `N_TUBES = None` = the cell's fixed reference count `CellSet.N_t_ref` (P1 34 / P2 32 / P3 30, see
  cooler-comparison.md), ambient 20 °C / 60 %.
- **Cooler comparison** (`comparison/`): the standard coolers `COMPARE_SETS` scaled to the same water flow and compared
  over the climate
  year with ONE solver (`COMPARE_SOLVER`), see cooler-comparison.md.
- **Shared**: `src/cooler_sets.py` (P1-P3), `src/weather.py` (the Munich climate year), coolant temperatures and target,
  the pad path and prices.
  The section headers of `parameters.py` carry the mode, `# --- Name [solver|compare|both]: ... ---`; the GUI builds its
  forms from them.

## Shared helpers

- `src/catalogue.py`: `Item(key, title, compute, plot, log)` and `run_items`. Both workflows offer their outputs /
  analyses as `Item` lists
  (`src/outputs.py`, `comparison/analyses.py`). To add one, append an `Item`; the selection lists, the CLI and the GUI
  follow.
- `src/figures.py`: `save_figure(fig, folder, name, title)` writes the PNG to `thesis_figures/out/<folder>/`. The
  folders follow the thesis chapters:
  `ch3_methodology/` (schematics), `ch4_results/solver/`, `ch4_results/cooler_performance/`, `ch4_results/geometry/`
  (`figures.CH3 / SOLVER / COOLER / GEOMETRY`).
  `save_data` / `load_data` keep raw results next to the figures (`*.pkl`), so a long run is never lost and figures can
  be redrawn without it.
- `src/plot_style.py`: the TUM Corporate Design colours (`TUM_*`; LMTD green, NTU blue, NTU table light blue, Cell
  orange, Cell 2D a tint; the comparison uses
  P1 blue, P2 orange, P3 green; `#0065BD` is also the Streamlit `primaryColor`), `print_context()` for the print-size
  graphs and `large_format_context()` /
  `@in_large_format` for the overview figures. Styles are contexts, never seaborn's global theme (the GUI shows both
  workflows in one process).
- `src/quiet.py` (`silence_all` for benchmarks, `quiet_solvers` for the studies), `src/formatting.py` (`fmt_ms`,
  `fmt_duration`, `machine_info`),
  `src/profiling.py` (benchmark instrumentation, see gui-and-inputs.md), `weather.hourly_states` (the operating points
  of a year).

### Figure style

Print-size graphs are 6.9 in wide at 300 dpi (`PRINT_WIDTH`, `PRINT_DPI`, `PRINT_STYLE`: 6-8 pt text, `PRINT_LW` 1.1 pt
lines, `PRINT_DASHED` for second cases and reference
lines), have no figure titles (definitions go in a footer or inside the graph), put legends inside the graphs where
there is room and draw the grid only along the
value axis in bar charts. `print_context(talk=True)` (solver figures) and `print_context()` (cooler-comparison figures)
reproduce the two looks of the thesis.

## Run modes

Dependencies: `requirements.txt`. scipy is needed for the NTU table mode and the annual benchmark, streamlit and pandas
only for the GUI. The NTU table cache
(`cache/`) and the rendered figures (`thesis_figures/out/`) are gitignored.

Solver-comparison settings live in `parameters.py` and can be overridden per run: `python main.py solver --solvers ...
[--no-]insight --outputs ...`
(`--outputs` with no value = summary only). `SOLVER_OUTPUTS` is the default selection (GUI: one multiselect). Outputs
run in catalogue order and their
results are in `results['outputs'][key]`.
- `INSIGHT_MODE`: per-iteration progress and full result boxes in the terminal (logging at DEBUG). Without it only the
  per-scenario summary (Q, outlets, pinch,
  deviation from Cell, iterations, wall time). Project code logs through `logging.getLogger(__name__)`, no `print` (only
  a few standalone `__main__` helpers print).
- `SOLVERS`: which of the five solvers run. The scenario results show one NTU and one Cell (`solvers.scenario_keys`: one
  solver per family; the variants give identical
  results and differ only in speed, which the benchmarks compare).

## Solver outputs

- `cooler_results`: per scenario (dry / precooled) a large-format figure, one row per solver: input / economics header,
  result box (deviation from the reference: pinch in K, Q in %) and its profile panel (T, k, dQ/dl along the coolant
  path; dotted = row changes). Identical axes in both figures (`plot_profiles.common_ticks`). LMTD's profile is the
  log-mean reconstruction
- `profile_grids`: profiles at the run's water temperature (solid) and `PROFILE_T_COOLANT_HOT` (90 °C, dashed):
  `profile_grid_absolute`, `profile_grid_normalised` (theta = (T - T_a,i)/(T_c,i - T_a,i), k / mean k,
  (dQ/dl)/(Q/l_tot)), `profile_collapsed_normalised` (all solvers in one graph per quantity), `profile_results_table`
  (the result boxes of both temperatures)
- `convergence`: outlet temperatures and iteration errors per scenario
- `caching`: caching benchmark + `caching_times`, `caching_split`, `caching_speedup`, the presentation sequences
  `caching_presentation_0 ... _4` and (with Cell 3D) `caching_presentation_with3d_0 ... _5`, and `caching_counts.tex`
- `resolution`: resolution sweep + `resolution_performance`
- `annual`: the annual benchmark + `annual_times`, `annual_deviation`
- `accuracy`: deviation of the solvers from the reference over the water inlet temperature (`ACCURACY_T_COOLANT`) +
  `accuracy_sweep` (both graphs) and `accuracy_sweep_pinch` (the pinch graph alone)

**Caching benchmark** (`benchmark_solvers.py`, fresh silent runs on the ambient operating point; every selected solver
variant in the cache modes no / cold / warm).
`caching_times`: run time per solver and mode (log axis, the solver colour in three shades). `caching_split`: one graph
per solver, one stacked bar per mode (CoolProp
lookups / cached lookups / k evaluation / rest of the solver, and the table build for NTU table). `caching_speedup`:
cold and warm cache against no cache (NTU table's warm
factor excludes the build). The presentation frames show the warm-cache times on one linear axis with the bars appearing
one by one (the last two are almost invisible, their
written times show them). The NTU table build is counted in the no-cache and cold-cache runs, never in the warm one; it
is always built fresh over the benchmark axis
(NTU_a <= 6). About 40 s with all five solvers. Instrumentation: `src/profiling.py`.

**Resolution benchmark**: sweep over `BENCHMARK_RESOLUTIONS` (normal runs use `CELL_N_SEGMENTS` / `NTU_N_ELEMENTS` =
20); `resolution_performance` (two graphs, all axes
logarithmic): the pinch error |dT_pinch(N) - dT_pinch(N_max)| against N (NTU and Cell lie on one curve; Cell 2D = 3D;
the table reproduces the field) and the wall time
against N (Cell 3D / 2D, NTU field, LMTD as a line, NTU table prebuilt, the table BUILD as a dashed series; the default
N dotted). The error falls as N^-2 (a straight line on
log-log, R^2 0.9997; semi-log gives 0.87); Cell 3D time ~ N^0.9, NTU field ~ N^0.8, the table build ~ N^2. Tables are
built fresh for every N up to `BENCHMARK_TABLE_MAX_N` (100;
0.3 s at N = 5, 6 s at 20, 36 s at 50, 126 s at 100). The prebuilt-table series is drawn only where the table is valid
(R axis reaches 0.1 N; the design point has R about 1.4,
so N >= 15), where a table solve takes about 1 ms (cold property cache; 0.04 ms warm) independent of N. Median of
`BENCHMARK_SWEEP_REPEATS` (3; 15 for the table solves) per point
with `gc.collect()` before each run. About 6.5 min with all five solvers (the builds take about 5 min). Raw results:
`resolution_results.pkl`.

**Annual benchmark** (`benchmark_annual.py`): the selected solvers over the 8760 hourly operating points of the climate
year (`src/weather.py` or `ANNUAL_WEATHER_CSV`): LMTD,
NTU field and NTU table every hour, Cell (2D / 3D) every `ANNUAL_CELL_STRIDE`-th hour (4) with its year total
extrapolated (hatched bar, "≈"; a solver over
`ANNUAL_TIME_BUDGET_S` (120 s; 0 = no limit) is stopped and extrapolated the same way). Only the inlet air (T, phi, p)
changes; dry operation, no pad, no fan control. The property cache
is on and emptied per solver; the hours come in the same seeded random order. Per hour the time, Q and T_c,o are stored.
`annual_times`: cumulative wall time over the points
solved (the axes follow the non-Cell solvers) next to the year totals on a log axis (NTU table prebuilt and including
its build); `annual_deviation`: pinch deviation in K of LMTD and NTU
from the most detailed Cell run against the ambient temperature (only if a Cell variant ran). The raw results are saved
first (`annual_results.pkl`).
The whole year for every solver:
`MPLBACKEND=Agg venv/bin/python main.py solver --solvers lmtd ntu_field ntu_table cell_2d cell_3d --outputs annual
--annual-cell-stride 1 --annual-time-budget 0`
(Cell 2D about 6 min, Cell 3D about 100 min; keep the machine awake and otherwise idle). The run stored with the thesis
(all five solvers, every hour, 6290 s in total): LMTD 5.4 s,
NTU field 43.0 s, NTU table 1.6 s (7.5 s with the 5.8 s build), Cell 2D 5.3 min, Cell 3D 98.6 min; deviation from Cell
3D over all 8760 h: LMTD up to 1.78 % in Q, NTU table up to 1.44 %.

**Accuracy sweep** (`benchmark_accuracy.py`): the selected solvers (one per family) against the reference over
`ACCURACY_T_COOLANT` at the run's geometry and ambient point: Q deviation [%]
and pinch deviation [K], two graphs, the run's own water temperature dotted. Needs two solvers; about 5 s with Cell 3D.
At the P2 cooler and the ambient point, from 30 to 90 °C water,
LMTD deviates +0.43 ... +0.68 % in Q and -0.03 ... -0.34 K in the pinch, NTU +0.06 ... +0.26 %: the solvers do not drift
apart at high water temperatures.
The reference uses the same correlations and element relation as NTU, so these are differences between thermal
calculations, not a validation.

## Chapter-3 schematics

`thesis_figures/chapter3_schematics.py` (standalone, nothing in the model imports it): `python -m
thesis_figures.chapter3_schematics [--format png] [--out folder] [--rows N] [--segments N]`
draws the discretisation from the actual geometry and the solvers' own path (`cell_path_order`), coloured by a real
model run: the tube circuit (6 rows x 10 segments, only the four
temperature symbols as text), `marching_order` (one iteration of `sweep_two_way` on one circuit, topology only),
`precooler_paths` (dry against precooled), the 3D staggered Cell cells
(`CELL_VERSIONS`), `cells_exploded_3rows` and `cells_bend` / `cells_bend_real` (U-bends at the back). Qualitative on
purpose: no numbers, the colour bar only says cold / hot.
Fins are drawn as a few plates per cell, not at the real pitch. Output: `thesis_figures/out/ch3_methodology/`.

Pinch = T_c_o - T_amb, the ambient air entering the whole cooler (pad and coil), also when precooled (it can then be
negative); `pinch_coil` is the value against the air entering the coil.
One definition: `run_scenario.calc_pinch(T_c_o, T_amb)` (solver workflow) and `comparison.evaluate.run_point`
(comparison). Shown per solver with the LMTD / NTU deviations from Cell.
