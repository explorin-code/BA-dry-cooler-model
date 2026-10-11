# The two workflows, their outputs and benchmarks

## Two workflows (restructured 2026-10; renamed 'Solver comparison' / 'Cooler comparison')

Everything runs in one of two workflows -- `python main.py solver|compare` (default `solver`) and the two modes of the GUI
(radio at the top of the sidebar, own parameter groups and run button each):
- **Solver comparison** (`src/`): ONE cooler in an ambient state, the SELECTED solvers (`SOLVERS`: lmtd, ntu_field, ntu_table, cell_2d,
  cell_3d -- only these run; the most detailed one is the reference for deviations and the precooling decision) compared, ambient +
  precooled scenario, solver benchmarks (cache/resolution compare the algorithm families of the selection, annual runs the selected
  variants). No prices here (costs are not part of a solver comparison). The cooler is a standard cell (`GEOMETRY_SET` P1/P2/P3 in parameters.py, `src/cooler_sets.py`),
  every cell value can be overridden (`None` = preset); module = `N_TUBES` x `N_ROWS` x `HEIGHT_MM`. Defaults (2026-10) = the comparison's rule: 1 kg/s, 2000 mm tubes, 6 rows, `N_TUBES = None` -> the cell's fixed reference tube count `CellSet.N_t_ref` (P1 34 / P2 32 / P3 30, derived by `comparison.sizing.design_tubes`, see cooler-comparison.md). Ambient 20 degC / 60 % (the Munich year's mean RH at 20 +-1 degC is 58 %), same as the comparison design point.
- **Cooler comparison** (`comparison/`): the standard coolers `COMPARE_SETS` scaled to the same water flow and raced against
  each other over the climate year with ONE solver (`COMPARE_SOLVER`, any of the five), see cooler-comparison.md.
- **Shared**: `src/cooler_sets.py` (P1-P3), `src/weather.py` (the Munich climate year, used by the solver annual benchmark
  AND the comparison), coolant temperatures / target, pad path, prices. `parameters.py` section headers carry the mode:
  `# --- Name [solver|compare|both]: ... ---` (the GUI builds its forms from that).

## Shared helpers (2026-10)

`src/figures.py` (`save_figure(fig, folder, name, title)`: PNG in `thesis_figures/out/<folder>/` + window title; solver-comparison figures now go to
`out/solver/` named after their title, e.g. `cooler_results_dry.png`, `caching_performance_log.png`; the cooler comparison to `out/ch4/`; the old
`KEEP_OPEN` global is gone -- the comparison runner closes the figures of a finished analysis unless `keep_open`), `src/quiet.py`
(`silence_all` for benchmarks, `quiet_solvers` for the studies; was two private copies), `src/formatting.py` (`fmt_ms`, `fmt_duration`,
`machine_info`), `weather.hourly_states` (the operating points of a year, used by the annual benchmark and the comparison's year analyses).

## Run modes

Dependencies: `requirements.txt` (`venv/bin/pip install -r requirements.txt`);
scipy is needed for NTU table mode and the annual benchmark, streamlit and
pandas only for the GUI. Rendered figures (`thesis_figures/out/`) and the NTU
table cache (`cache/`) are gitignored.

Solver-comparison settings in `parameters.py`, each overridable per run (`python main.py --solvers ... --[no-]insight --outputs ...`). The
outputs are ONE catalogue (2026-10, `src/outputs.py` + `src/catalogue.py`, the same `Item`/`run_items` pattern as `comparison/analyses.py`):
`SOLVER_OUTPUTS = [...]` = any of `cooler_results`, `convergence`, `caching`, `resolution`, `annual` (replaces the five booleans PLOT_RESULTS,
PLOT_CONVERGENCE, BENCHMARK_MODE, RESOLUTION_MODE, ANNUAL_MODE; `--outputs` with no value = summary only; GUI: one multiselect, as `COMPARE_ANALYSES`
for the cooler comparison). An output = compute / plot / log of a run context; to add one append an `Item` to `OUTPUTS`. They run in catalogue order,
results in `results['outputs'][key]`.
- `INSIGHT_MODE`: per-iteration progress + full result boxes in the
  terminal (logging at DEBUG). Without it only the per-scenario summary
  (INFO: Q, outlets, pinch, deviation from Cell, iterations, wall time).
  Project code logs via `logging.getLogger(__name__)`, no `print`
  (except the standalone `validation_ntu.py`).
- Colours (2026-10): the solver graphs use the TUM Corporate Design (Handbuch 2016, p. 20-21), `plot_style.py` `TUM_*`: LMTD = accent green
  #A2AD00, NTU = accent blue #64A0C8 (NTU table = the second accent blue #98C6EA), Cell = accent orange #E37222 (Cell 2D = a tint), caching time
  parts in the TUM blues/grey/ivory, result boxes tints of the family colour. The handbook prints Pantone 300 as RGB 0/101/189 but Hex #3070b3
  (inconsistent) -- #0065BD is used (also as Streamlit `primaryColor`, `.streamlit/config.toml`). The cooler-comparison figures use P1 = blue, P2 = orange, P3 = green (and
  TUM greys/blue for fixed-vs-controlled and water), Arial is not set.
- Figure / GUI-tab names (solver comparison): 'Cooler results - Dry|Precooled', 'Convergence - Dry|Precooled', 'Caching performance (linear|log)', 'Resolution performance', 'Annual performance' -- performance figures do not say which scenario they ran on (always the dry ambient point).
- `cooler_results`: per scenario a cooler-results figure (identical T / k / dQ-dL axes in the dry and precooled figure, `plot_profiles.common_ticks`; inputs, economics,
  each solver's result box next to its spatial profile).
- `convergence`: per scenario outlet temperatures + iteration errors.
- `SOLVERS`: which of the five solvers run (lmtd, ntu_field, ntu_table, cell_2d, cell_3d). The SCENARIO results (profiles, tables,
  convergence) show one NTU and one Cell -- the variants give identical results and differ only in speed, which is the benchmarks'
  business (`solvers.scenario_keys`: one solver per family; labels LMTD / NTU / Cell). The benchmarks below differentiate the variants.
- `caching`: caching benchmark (`benchmark_solvers.py`, fresh silent runs on the ambient operating point) + solver-performance
  figure (`plot_solver.py`), drawn twice -- on linear and on logarithmic time axes (the solvers span ~4 orders of magnitude; the log page shows plain bars without the time-part split, which would mislead there; the old 'rest of solver' panel was removed): every selected
  solver variant in the cache modes no / cold / warm, and each solver's time split into CoolProp / cached lookups / k / rest. NTU (table) additionally shows 'table build' (fresh build, never from disk): counted in the
  no-cache and cold-cache runs, assumed done in the warm-cache run. ~25 s with all five (Cell 3D without cache ~10 s, the table build ~6 s).
- `resolution`: resolution sweep over `BENCHMARK_RESOLUTIONS` (normal runs use `CELL_N_SEGMENTS` / `NTU_N_ELEMENTS`) + its figure (wall time, time per iteration, iterations, T_c,o):
  NTU (field only -- the table gives identical results, so the two bottom panels are labelled 'NTU (field/table)' and the top panels
  'NTU (field)'; the table's cost is the caching benchmark's business; `BENCHMARK_TABLE_MAX_N` removed 2026-10), Cell 2D / 3D, LMTD as a flat line.
  The pinch panel runs from 0 K with a zoom inset. ~30 s with all five.
- Cell 2D vs. 3D: with uniform inlet air identical to machine precision (all tubes behave the same, so the neighbour averaging of air
  has no effect); 2D is ~17x faster with a warm cache, ~2.5x on a cold single run (same CoolProp calls). 3D only matters for
  non-uniform inlet air / per-tube conditions.
- `annual` (`benchmark_annual.py`): annual benchmark -- the selected solvers over the 8760 hourly operating points
  of the climate year: LMTD, NTU field and NTU table over every hour, Cell (2D / 3D) on every `ANNUAL_CELL_STRIDE`-th hour (4) as the
  accuracy / cost comparison, its year total extrapolated (hatched part of the bar, "≈"; a solver over `ANNUAL_TIME_BUDGET_S` (120 s)
  is stopped and extrapolated the same way). With all five selected the figure has six bars: NTU (table) appears twice, table prebuilt
  and incl. build (one solid bar with the whole time). Only inlet air T/phi/p change (coolant inlet and flows from parameters.py), dry
  operation only. Weather: the Munich climate year of `src/weather.py` (the synthetic year was removed 2026-10) or `ANNUAL_WEATHER_CSV`.
  Property cache on, emptied per solver; hours in the same seeded random order. Log: year totals, LMTD / NTU deviation from Cell on the
  common hours. Figure: cumulative time, log-scale totals, the year's T_a,i and Q. Minutes.
Thesis figures (Chapter 3): `thesis_figures/chapter3_schematics.py`, run with
`python -m thesis_figures.chapter3_schematics [--format png] [--rows N]
[--segments N]` -> writes every
accepted version into `thesis_figures/out/`: the tube circuit (serpentine via
`cell_path_order`, coloured by local mean coolant/air temperature; always
6 rows x 10 segments, `CIRCUIT_ROWS`/`CIRCUIT_SEGMENTS`, only the four
temperature symbols as text), `precooler_paths` (dry vs. precooled: bypass or
wetted pad into the cooler, each cooler coloured by its own run) and the 3D
staggered Cell cells (3 tubes of row r-1, 2 of row r; air of two lower cells
merges into each upper cell and splits again), versions in `CELL_VERSIONS`
(exploded, 1 and 2 segments along the tube; the row gap is a parameter),
`marching_order` (one iteration of `sweep_two_way` on one circuit, 2D:
pass 1 numbered in `cell_path_order`, pass 2 row by row; topology only),
`cells_exploded_3rows` (own 3-row run: air inlet, inner and air outlet row,
3/2/3 tubes, 2 segments), plus `cells_bend` / `cells_bend_real`:
2+2 cells at the tube ends with the U-bends at the back, rows pulled apart /
touching. Qualitative on purpose -- no
numbers, colours from a real Cell run, colour bar only cold/hot. Fins are
drawn as a few plates per cell, continuous across a row, not in the gap
between pulled-apart rows (not the real pitch); in the model they
only enter via A and the fin efficiency in k. Nothing in the model imports
this file (could later be rendered automatically per run).
Pinch point = T_c_o - T_a_i (cold end), shown per solver with
LMTD/NTU deviations from Cell.
