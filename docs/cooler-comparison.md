# Cooler comparison (`comparison/`)

`python -m comparison [analysis ...|all] [--plots-only] [--stride N] [--solver KEY] [--sets P1 P3] [--sensitivity-set P2]`,
`python main.py compare ...` or the GUI's comparison mode. No names = the analyses switched on in parameters.py
(`COMPARE_ANALYSES` = the default selection). CSVs in `output/comparison/` (gitignored), figures in `thesis_figures/out/ch4/` (axis labels in thesis
symbols via `plots_common.LABEL`; figures stay open for plt.show()/the GUI). Structure (read `comparison/__init__.py`, then
`analyses.py` -- the catalogue: every analysis = compute / plot / log): `config.py` (2026-10: `settings()` reads parameters.py ONCE per run into a frozen `Settings`
object `cfg` that is passed as the first argument to everything that needs a setting -- no module-level live view, no values frozen at import; the solver
tuning of the comparison therefore follows the GUI like the solver mode does), `evaluate.py` (one operating point = one result dict), `sizing.py`, `operation.py` (states,
year, face-velocity sweep, fan control, summaries), `fan_curve.py`, `sensitivity.py`/`rows.py`/`tubes.py`/`wall.py`,
`plots_common|design|operation.py`, `plot_geometry.py` (module drawings), `runner.py`. Analyses: geometry, ambient_cases,
airflow_sweep, fan_control (incl. the fixed-speed year), sensitivity, rows, tubes, wall. The solver is chosen once for all of them
(`COMPARE_SOLVER`; 'ntu_table' is the default, analyses that change the rows/cell per point use the table-free `ntu_field`;
the Cell variants are ~100x slower -- use the stride with year analyses). The sensitivity cell is `COMPARE_SENSITIVITY_SET`
(GUI: selector at the top of the comparison mode, CLI `--sensitivity-set`).

- **Weather** (`src/weather.py`): DWD station 03379 Muenchen-Stadt, derived files in `climate_data/` (README there; made by `python -m src.build_climate_data` from the raw download, which now lives outside the project as `raw_climate_data/`, gitignored and deletable);
  last 365 days of the file = 2025-10-09 .. 2026-10-08 (the file has no full calendar year), UTC,
  51 interpolated hours, RH capped at 99 % (100 % breaks the property cache), measured station
  pressure (~955 hPa) used as `p_a`. Representative cases: coldest, median, P95, P99, hottest,
  wettest_warm (highest wet-bulb), `design` (`COMPARE_DESIGN_T_AIR`/`_PHI_PERCENT`, default 20 degC / 60 %).
- **Elementary cells** (`src/cooler_sets.py`, from the Perplexity audit "Project-Source Geometry Audit
  for Plain Plate-Fin Dry-Cooler Cases", NOT re-read in the PDFs): P1 Wang & Chi 2000 sample 2,
  P2 Kayansayan 1994 coil 5, P3 Kayansayan coil 3. Only the elementary cell (d, s_f, delta_f, s_t,
  s_l, fin type, materials) is from the sources. **Assumed by me**: tube wall (0.5 mm; P3 1.5 mm),
  hence d_i; collar thickness = fin thickness for all three (the collar is the fin sheet; the sources give no value);
  P1's listed diameter 7.3 mm is over the collar (Wang's D_c), so the model's tube d = 7.3 - 2 delta_c = 7.07 mm
  (`cooler_sets._tube_od`); P2/P3's listed 9.52 / 16.3 mm are the tube d (d_c = 9.92 / 16.7); N_t 58/43/26 -> 61/43/26;
  aluminium fins/copper tubes for P1; every tube = its own N_r-pass circuit (sources: 4 parallel
  multipass circuits); densities for the masses (Al 2700, Cu 8960, steel 7850 -- cite before use).
  P3 is outside Wang's range (d, s_t, s_l), P1 too at N_r > 6.
  A fourth "industrial" cell (d 25.4 mm, s_t 63.5 mm, 400 fins/m, Kroger/Towler) was tried and
  dropped (circular-fin tubes, a different problem; Wang range exceeded).
- **Scaling to 1 kg/s**: tube length one fixed value `COMPARE_TUBE_LENGTH_MM` (2000, a design choice -- 1.5 m left
  too many summer hours above the target; l changes the behaviour, it is not a linear scale
  factor), N_t fixed per cell (below), rows 6 and tube length 2 m fixed for all; the tube count is the cell's FIXED reference count
  `CellSet.N_t_ref` = P1 34 / P2 32 / P3 30 (2026-10), used by BOTH workflows and independent of every operating input (coolant temperature,
  flow, ambient, solver choice, ...): only an explicit number (`COMPARE_N_TUBES`, `N_TUBES`, GUI field) changes it, for that run. Derivation
  (`comparison/sizing.design_tubes`, `python -m comparison.sizing` prints it next to N_t_ref -- refresh N_t_ref when the cells or the solver
  change): the SMALLEST N_t that brings 1 kg/s of 37 degC water to 25 degC at 20 degC / 60 %, dry, 3 m/s, 6 rows, 2 m (solver ntu_field; P2 gets 31
  with LMTD). Wall: P1/P2 0.5 mm, P3 1 mm (Re_c ~3800 at P3, below the 4000 mark but not laminar; P3 with 0.5 / 1 mm walls at 2.5 m/s would need 66
  tubes at laminar Re_c ~1700). T_c,o is not monotone in N_t because of the transition band. The earlier Re = 5000 rule (61/43/26 tubes) is only the
  reference of the tubes/wall studies (`COMPARE_RE_TUBE`). Duty and T_c,o are RESULTS, not targets. Module = elementary rows scaled by N_t; 2 kg/s = 2 modules.
- **Operation**: only TWO fan modes (no boost mode, removed 2026-10): fixed speed `COMPARE_W_FR_DESIGN` = 3 m/s (also the speed the tube
  count is sized for and the TOP of the variable range; from the user's reading -- face velocities above ~3 m/s risk droplet
  carry-over/noise; cite; Towler's typical 2.5), or variable speed between `COMPARE_W_FR_MIN` 0.4 m/s and the design speed. `operation.control_point`: variable-speed fan, per hour the SMALLEST face
  velocity in [0.4, 3] m/s with T_c,o <= 25 degC (brentq on w_fr, ~10 solves), dry first, then
  precooled, else precooled at 3 m/s ("limited" = target missed). Precooling pad: eta_B 0.8
  placeholder, 50 Pa. Freezing risk flag `freeze` = T_c,o < `COMPARE_T_FREEZE_RISK` (2 degC, assumption).
- **Results so far** (run of 2026-10-10 with the fixed counts 34 / 32 / 30, 6 rows, 2 m, 3 m/s design speed; not yet reviewed):
  design point 20 degC / 60 % dry: all three give Q = 50-51 kW and T_c,o ~ 24.9 degC by construction (P1 891 W fan, 23 kg Al + 38 kg Cu; P2 622 W,
  46 + 49 kg; P3 861 W, 55 + 155 kg). Fixed 3 m/s over the year: target met in 8380 / 8395 / 8358 of 8760 h (P1 / P2 / P3), 0 h with freezing risk.
  Fan control: fan 2568 / 1639 / 2031 kWh/a, pump 412 / 97 / 10 kWh/a (P1's few tubes make the water fast), water 67 / 85 / 105 m3/a, annual cost
  incl. the 20-year material write-off 1044 / 714 / 909 EUR/a (P2 cheapest); 365-402 h per year stay above 25 degC even at 3 m/s with the pad.
  Sensitivity (P2 base): rows (+8.3 K at 1 row) and face velocity (+8 K at 0.5 m/s) matter most, then the fin pitch; the wall hardly (< 0.6 K).
  Tube count (`tubes.py`): Q rises with N_t down to Re ~3000, then falls (transition band) while fan power and mass keep growing.
  Wall (`wall.py`, 0.3-2 mm): a thicker wall gives a few % more Q (N_t fixed) at a multiple of the copper mass; N_t kept vs. Re kept both compared.
  P3's wall (1 mm) leaves Re_c ~3800 at its 30 tubes -- below 4000, not laminar; with 0.5 / 1 mm and 2.5 m/s it would need 66 tubes (laminar).
- **Fan curve** (`fan_curve.py`, part of the fan_control analysis): the solved control collapses onto monotone curves
  w_fr(T_a,i) and P_fan(T_a,i) per set, separately for dry and precooled operation (the pad adds 50 Pa and
  humidity; T_a,i = air entering the cooler). PCHIP fit through binned medians of control.csv; schedule per
  hour needs only T_amb and wet-bulb (T_a,i = T - eta_B (T - T_wb), cooling limit ~ thermodynamic wet-bulb).
  Annual fan energy from the curves within -3 ... +2 % of the solved control, mean |dw| 0.01 m/s (single
  hours at the dry/precooled switch off by up to 2 m/s). Valid for the study module and weather only.
- **Material write-off** (2026-10): copper `COPPER_EUR_PER_KG` 12.4 / aluminium `ALUMINIUM_EUR_PER_KG` 2.7 EUR/kg (raw-metal LME, Oct 2026,
  14 300 / 3 134 USD/t at 1.15 USD/EUR -- no fabrication premium; sources in the search of 2026-10-10: Cu 6.49 USD/lb on 2.10., Al cash on 6.10.),
  straight-line over `LIFETIME_YEARS` = 20 (no interest/scrap), masses from `evaluate.calc_masses` (straight tube + fin only): `operation.material_writeoff`,
  column `cost_material`, part of `cost_total` in `annual_energy` (the old bottom row of hours was removed: identical heat load design).
- **Cost and water** (parameters.py prices): electricity 0.30 EUR/kWh (user), water 1.84 EUR/m3 net (Munich drinking
  water from 2026-07-01 incl. Wassercent, SWM press release 3.6.2026; gross 1.97). Wastewater fee 0 by default
  (evaporated pad water normally not charged; Munich sewage ~2 EUR/m3, 2023 sample bill, outdated). Water is tracked
  per hour (`m_dot_ev`, VDI M8 humidity-ratio increase x dry-air flow) and shown in `fan_year`,
  `annual_energy` (water + cost panels) and `cost_year` (cumulative and monthly cost, electricity and water separate).
  A design-space "optimization" mode (Latin hypercube over the cell variables, rank regression, Pareto) was tried and
  removed 2026-10 because its figures were not informative; the one-at-a-time sensitivity, rows/tubes/wall studies stay.
- **Cell vs NTU at cold ambient**: Cell ends 2-3 K warmer than NTU/LMTD (P2 coldest hour 3.6 vs 1.0 degC)
  -- local viscosity/k of cold water (Cell) vs. mean-temperature properties; not investigated further.
- **Cell failure at the coldest hour (fixed 2.5 m/s, P1/P2)**: NTU/LMTD also put T_c,o below 0 degC
  there -- the water would freeze, and CoolProp has no liquid water below 0 degC. Fixed 2026-10:
  `fluid_properties.get_fluid_properties` clamps pure water to >= 0.01 degC (`WATER_T_MIN`), Cell
  then converges; results below 0 degC are flagged (`freeze`) and not physical.
- **Open**: two-parameter maps and a global sensitivity are planned, not built. Sensitivity (`sensitivity.py`) is one-at-a-time on the P2 base module; the
  `s_t` sweep is noisy (N_t integer rounding changes the face area).
