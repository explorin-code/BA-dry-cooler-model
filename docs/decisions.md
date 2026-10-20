# Decisions and changes

What was decided and why, and which defects were found and corrected on the way. Numbers are those at the time of the
change.

## Modelling decisions

- **Solid conductivities**: Cu 380, Al 160, steel 50 W/(m K), VDI-Wärmeatlas (2019) D6 Tab. 8 (design values of DIN EN
  ISO 10456), instead of the pure-metal 401 / 237 (k dropped about 5 %, Q about 2 %).
- **No fin collar by default** (`FIN_COLLAR_THICKNESS_MM` = 0, `Geometry.delta_c` = 0, d_c = d). `delta_c` only enters
  Wang's collar diameter d_c (air friction factor, sigma); heat-transfer areas, fin efficiency
  and Nu_a keep using d.
- **Water use** uses the dry-air mass flow m_dot_a/(1+Y) (thesis Eq. 2.28).
- **Cooling limit**: the enthalpy is per kg DRY air (`calc_h` uses CoolProp `Hda`, VDI's h_1+X); T_K at 20 °C / 30 % is
  10.84 °C, matching CoolProp's wet-bulb temperature (it was 10.89 °C with `Hha`; precooled
  Cell Q 17.51 -> 17.54 kW, water 3.692 -> 3.665 g/s). c_W comes from CoolProp (liquid water at T_K); the formula is VDI
  M8 §2.2 Gl. (15).
- **Design values**: 37 °C inlet, 25 °C target and 0.28 kg/s are the reference values of the task; the solver defaults
  and the comparison use 1 kg/s modules (`M_COOLANT`).
- **Cell**: air averaging between neighbouring tubes follows from the staggered geometry and is a modelling assumption,
  not a correlation. Each fluid's properties are evaluated at its own local mean temperature.
- `DEFAULT_PHI_AIR_FALLBACK` = 0.30 is an assumption (used only when no humidity input is given).
- **Fin efficiency**: the code (VDI M1 Eqs. 7, 8, 12, 14 with the plate-fin factor phi) is the reference for thesis Eq.
  2.57.
- **Pinch** is taken against the AMBIENT air entering the whole cooler, also when precooled (the pad belongs to the
  cooler); `pinch_coil` keeps the coil-only value. Affects the precooled scenario of the solver
  workflow (log, result box, GUI table) and the `pinch` column of every precooled comparison row. Not affected:
  everything computed dry and all deviations between solvers.
- **One target tolerance**: the target is met when T_c,o <= 25 °C + `config.TARGET_TOL` (0.01 K), used by the row flag,
  the controller's search and the hour counts.
- **Operating modes of the annual economics**: fixed dry, adaptive dry and adaptive wet. Costs are computed at plot time
  from the stored columns (`operation.add_costs`), so prices and the cost model can change
  without new solves.
- **NTU table**: R = C_a/C_c grows with the element count (`ntu_table_R(n_elements)`, R_max = 0.1 N_e = the validity of
  Eq. 1); studies use 60 elements (R <= 6). The NTU_a axis runs 0.05 ... 60 (64 points,
  `NTU_TABLE_VERSION` = 2) because the adaptive fan at 0.4 m/s and coolers with more rows or area reach NTU_a far above
  6, where the table fell back to a field solve (about 100x slower per call: 26 % of the steps
  at N_r = 10). The solver benchmarks keep the first range (0.05 ... 6, 48 points, `NTU_TABLE_NTU_BENCHMARK`) so the
  thesis' build times keep their meaning. Pinned in `tests/golden_solvers.json`
  (the wide grid changed the table solves by at most 1.5e-6 relative, same iteration counts).
- **Parameter variation**: cost basis = the wet adaptive-fan year at every 6th hour; variants that miss the target for
  more than 500 h are dotted; `rows`, `tubes`, `wall` and `sensitivity` stay as separate analyses.
  A design-space optimisation mode (Latin hypercube over the cell variables, rank regression, Pareto) was tried and
  removed: its figures were not informative.
- **Figure folders** follow the thesis chapters (`figures.CH3 / SOLVER / COOLER / GEOMETRY`); print-size graphs share
  `plot_style.PRINT_*`; styles are contexts, never a global seaborn theme.

## Defects found and corrected

- **Coolant friction** used the total tube length (all passes) per pass, counting it N_r times too often (10.7 instead
  of 2.5 kPa, pump 3.35 instead of 0.79 W); now L = one pass (Towler 19.8, thesis Eq. 2.13).
- **Cell**: both fluids' properties were evaluated at the average of coolant and air inlet temperature (coolant about
  4-5 K too cold, k too low); the per-cell effectiveness used `calc_P` (Holman's empirical
  crossflow fit), which gave P about 17 % too low at Cell's per-cell R of about 19 -- now the same Cabezas-Gómez element
  relation as NTU (Holman is no longer used anywhere); Cell's dQ/dl used the cell's outlet
  temperatures and added up to only 89 % of Q (reported Q was right) -- all three profiles now integrate to their
  solver's Q within 0.05 %; Cell's grids became INTEGER when a temperature was given as a whole number
  (`T_COOLANT_IN = 80`), truncating every temperature (NTU already cast to float); the mean air temperature for k is now
  the cell mean (Cell k +0.2 %, Q +0.08 %).
- **Precooled operating point** rebuilt its flow areas from the geometry of `parameters.py` instead of the run's `geo`
  (harmless only when they were the same); `calc_precooler` takes `geo` and
  `calc_precooled_ops(ops, geo)` is split out.
- **Relaxation**: Cell with omega 0.2 stopped before convergence (Q +0.17 %) and its iteration count depended
  erratically on the resolution; the outer loop had no iteration cap and looped forever on divergence;
  LMTD's omega 0.2 stopped too early at the 1e-3 K threshold; the annual benchmark found limit cycles at the 0.01 K
  cache rounding (4 of 2000 LMTD hours and 9 % of the NTU hours never converged) -- fixed by the
  oscillation damping (model-notes.md, section 2).
- **Pure water below 0 °C**: Cell failed at the coldest hour with fixed 2.5 m/s (NTU / LMTD also put T_c,o below 0 °C);
  the lookup now clamps pure water to >= 0.01 °C and flags `freeze`.
- **Target tolerance**: the row flag used 1e-6 K and the hour counts 0.01 K, so the flag read True in only 2000-6000 of
  8400 hours; one tolerance now gives met + limited = 8760.

## Structure

- `src/solvers/` is a package (`results`, `common`, `field`, `relaxation`, `lmtd`, `ntu_field`, `ntu_table`, `ntu`,
  `cell`, `scenario`); `ScenarioResult.results` is a dict of the selected solvers (`.reference` = the most
  detailed), selected by the keys `cell_2d` / `cell_3d`, `ntu_field` / `ntu_table`.
- `cooler_sets.py` holds the cells P1-P3; `geometry.get_geometry` builds the geometry from the preset
  `GEOMETRY_SET` plus overrides (`CELL_PARAMETERS`, `get_cell_values`).
- `weather.py` holds the climate year and design points for both workflows; `benchmark_annual.py` uses it with the
  hour's measured station pressure.
- Pressure drops live in `pressure_drop.py` (`calc_delta_p_air` = bundle + pad of one state), powers and water in
  `economics.py` (`calc_fan_power(ops, T_a_o, Delta_p)`); both workflows call the same functions.
- `run_modes.py` / `scenario_pipeline.py` carry the solver-comparison flags only; the cooler analyses live in
  `comparison/`.
- Benchmark instrumentation is `src/profiling.py` (timing hooks active only inside `profiling.record()`), not patching
  of the model functions.
