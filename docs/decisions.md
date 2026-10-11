# Decisions and changes

The decision log: what was decided when and why (thesis chapter 2 review, the 2026-10 restructuring). Newest conventions of use are in `../CLAUDE.md`.

## Decisions recorded 2026-10 (from the thesis chapter 2 review)

- Solid conductivities: Cu 380, Al 160, steel 50 W/(m·K), VDI-Wärmeatlas
  (2019) D6, Tab. 8 (design values from DIN EN ISO 10456) -- all three
  checked against the table 2026-10 (was uncited pure-metal 401/237; k
  dropped ~5 %, Q ~2 %).
- No fin collar by default (`FIN_COLLAR_THICKNESS_MM` = 0 -> `Geometry.delta_c` = 0,
  d_c = d). The input exists (parameters.py, GUI) but delta_c only enters
  Wang's collar diameter d_c (air friction factor, sigma); heat-transfer areas, fin
  efficiency and Nu_a keep using d.
- Water usage uses the dry-air mass flow, m_dot_a/(1+Y) (thesis Eq. 2.28).
- Cooling limit enthalpy basis (fixed 2026-10): `calc_h` now uses CoolProp
  'Hda' (per kg dry air, VDI's h_1+X) instead of 'Hha' (per kg humid air);
  T_K at 20 °C / 30 %: 10.89 -> 10.84 °C, matching CoolProp's thermodynamic
  wet-bulb temperature. Precooled Cell Q 17.51 -> 17.54 kW, water
  3.692 -> 3.665 g/s.
- Cooling limit: c_W from CoolProp (liquid water at T_K) instead of the
  uncited 4186; the formula's source is VDI M8 §2.2 Gl. (15) (thesis
  updated from N4 to M8).
- Design constraints from the supervisor ("Konrad" in parameters.py:
  37 °C inlet, 0.28 kg/s, 25 °C target) are the reference values. The
  thesis Sec. 1.2 still says 2 kg/s -- to be corrected there.
- Air averaging between neighbouring tubes in Cell is a modelling
  assumption from the staggered geometry, not a correlation.
- `DEFAULT_PHI_AIR_FALLBACK = 0.30` is an assumption.
- Thesis Eq. 2.57 (fin efficiency) is cited/written incorrectly in the
  thesis; the code (VDI M1 Eqs. 7, 8, 12, 14 with the plate-fin factor
  phi) is the reference.

## Changes in `src/` made for the comparison and the restructuring (2026-10)

- `solvers.py`: the NTU table's R = C_a/C_c axis now grows with the element count
  (`ntu_table_R(n_elements)`, R_max = 0.1 N_e = Eq. (1) validity; 20 elements: 0..2 as before). The
  1 kg/s modules run at R ~ 2-5 (the old table stopped at 2 -> silent field-solve fallback). Studies
  use 60 elements (R <= 6, build ~50 s per N_r, then cached; beyond R 6 the field solve). Table range
  check now reads the interpolator's knots. Flow rates need no extra tables, only R/NTU inside range.
- `precooling.py`: `calc_precooled_ops(ops, geo)` split out of `calc_precooler`; **bug fixed**: the
  precooled operating point used to rebuild its flow areas from parameters.py's geometry (it did not
  pass `geo`) -- only harmless when the run used the parameters.py geometry. `calc_precooler` takes `geo`.
- pressure drops live in `pressure_drop.py` (`calc_delta_p_air` = bundle + pad of one state, was `scenario_pipeline.solve_air_side`), powers and water in
  `economics.py` (`calc_fan_power(ops, T_a_o, Delta_p)`); both workflows call the two functions (2026-10, no formula changed).
- `fluid_properties.py`: pure-water clamp, see above.
- `cooler_sets.py` (new): P1-P3 cells; `dry_cooler_physics.get_geometry` builds the geometry from the preset
  `GEOMETRY_SET` + overrides (`CELL_PARAMETERS`, `get_cell_values`) -- the old hand-made 9 mm / 17 tube geometry is gone,
  so the solver-comparison baseline changed (P2 cell: ambient Cell Q 11.7 kW instead of 12.45 kW; `W_O` 2.0 -> 2.5 m/s).
- `weather.py` (new, moved from the comparison): the climate year + design points; `benchmark_annual.py` uses it
  (synthetic year removed) and the hour's measured station pressure; table axis adapted.
- `solvers.py` (since 2026-10 the package `src/solvers/`): `SOLVER_LABELS`/`SOLVER_FAMILY`/`DETAIL_ORDER`, `run_solver(key, ...)`; `ScenarioResult.results` is a dict of the
  selected solvers (`.reference` = the most detailed; the old `.lmtd/.ntu/.cell` attributes are gone). `CELL_2D` was replaced by
  selecting `cell_2d` / `cell_3d`. Convergence, profile figures, benchmarks and the pipeline follow the selection.
- `run_modes.py`/`scenario_pipeline.py`: solver-comparison flags only (the cooler analyses live in `comparison/`).
