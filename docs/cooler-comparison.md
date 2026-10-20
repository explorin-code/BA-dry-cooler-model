# Cooler comparison (`comparison/`)

```
python -m comparison [analysis ...|all] [--plots-only] [--stride N] [--solver KEY] [--sets P1 P3] [--sensitivity-set P2]
python main.py compare ...          # the same
```
or the GUI's comparison mode. No names = the analyses switched on in `parameters.py` (`COMPARE_ANALYSES`). CSVs go to
`output/comparison/` (gitignored; `--plots-only`
redraws from them), figures to `thesis_figures/out/ch4_results/cooler_performance/` (the module drawings to
`.../geometry/`); axis labels use the thesis symbols
(`plots_common.LABEL`); figures stay open for `plt.show()` / the GUI.

## Structure

Read `comparison/__init__.py`, then `analyses.py` -- the catalogue; every analysis is compute / plot / log.
`config.py`: `settings()` reads `parameters.py` ONCE per run into a frozen `Settings` object `cfg` that is passed as the
first argument to everything that needs a setting (no module-level view, no value
frozen at import), so the solver tuning of the comparison follows the GUI like the solver mode does. `evaluate.py` (one
operating point = one result dict), `sizing.py`, `operation.py` (states, year,
face-velocity sweep, fan control, summaries, costs), `fan_curve.py`, `variation.py`, `sensitivity.py` / `rows.py` /
`tubes.py` / `wall.py`, `plots_print.py` (print-size figures), `plots_common.py` /
`plots_design.py` / `plots_operation.py` (large-format figures), `plot_geometry.py` (module drawings), `runner.py`.

The solver is chosen once for all analyses (`COMPARE_SOLVER`; `ntu_table` is the default, analyses that change rows or
cell per point use the table-free `ntu_field`; the Cell variants are about 100x
slower, use the stride with year analyses). The sensitivity cell is `COMPARE_SENSITIVITY_SET` (GUI: selector at the top
of the comparison mode, CLI `--sensitivity-set`).

| analysis | content | figures (in `ch4_results/...`) |
|---|---|---|
| `geometry` | the three modules | `geometry/geometry_print`, `module_geometries`, `bank_full`, `geometries` |
| `ambient_cases` | the representative ambient states | `ambient_cases` |
| `airflow_sweep` | pinch against face velocity | `airflow_sweep` |
| `fan_control` | the year with fixed and adaptive fan | `control_year`, `fan_year`, `year_scatter`, `fan_curves`, `cost_year`, `annual_energy` |
| `pinch_fan` | the year at fixed speed, the controlled year, annual economics (print size, about 2 min) | `pinch_year`, `pinch_year_clouds`, `missed_hours`, `pinch_speed`, `controlled_fan`, `annual_costs`, `monthly_cost`, `water_vs_hours`, `fan_cop`, `annual_costs_appendix.tex` |
| `variation` | one-at-a-time geometry variation, pinch and annual cost (print size, about 8 min) | `variation_costs`, `variation_totals` |
| `sensitivity` | one-at-a-time literature ranges on the P2 module | `sensitivity_lines`, `sensitivity_tornado` |
| `rows`, `tubes`, `wall` | rows, tubes per row, wall thickness | `rows_analysis`, `tubes_analysis`, `wall_analysis` |

## Inputs

- **Weather** (`src/inputs/weather.py`): DWD station 03379 München-Stadt, derived files in `climate_data/` (README there; made
  by `python -m src.inputs.build_climate_data` from the raw download, kept outside the repo in
  `raw_climate_data/`). The year is the last 365 days of the file, 2025-10-09 .. 2026-10-08 (the file has no full
  calendar year), UTC, 51 interpolated hours, RH capped at 99 % (100 % breaks the
  property cache), the measured station pressure (about 955 hPa) used as `p_a`. Representative cases: coldest, median,
  P95, P99, hottest, `wettest_warm` (highest wet-bulb) and `design`
  (`COMPARE_DESIGN_T_AIR` / `_PHI_PERCENT`, 20 °C / 60 %).
- **Elementary cells** (`src/inputs/cooler_sets.py`): P1 Wang & Chi (2000) sample 2, P2 Kayansayan (1994) coil 5, P3 Kayansayan
  coil 3. Only the elementary cell (d, s_f, delta_f, s_t, s_l, fin type,
  materials) is from the sources. Assumptions: the tube wall (0.5 mm; P3 1 mm) and with it d_i; the collar thickness
  equals the fin thickness (the collar is the fin sheet turned up around the tube);
  P1's listed 7.3 mm is over the collar (Wang's D_c), so the tube d = 7.3 - 2 delta_c = 7.07 mm
  (`cooler_sets._tube_od`), while P2 / P3's 9.52 / 16.3 mm are the tube d (d_c = 9.92 / 16.7); aluminium
  fins on copper tubes for P1; every tube is its own N_r-pass circuit. P3 is outside Wang's range (d, s_t, s_l), P1 too
  at N_r > 6.
- **Scaling to 1 kg/s**: tube length one fixed value `COMPARE_TUBE_LENGTH_MM` (2000; 1.5 m left too many summer hours
  above the target; l changes the behaviour, it is not a linear scale factor), 6 rows. The tube
  count is the cell's FIXED reference count `CellSet.N_t_ref` = P1 34 / P2 32 / P3 30, used by BOTH workflows and
  independent of every operating input (coolant temperature, flow, ambient, solver choice, ...); only an
  explicit number (`COMPARE_N_TUBES`, `N_TUBES`, the GUI field) changes it, for that run. Derivation
  (`comparison/sizing.design_tubes`; `python -m comparison.sizing` prints it next to `N_t_ref` -- refresh `N_t_ref`
  when the cells or the solver change): the SMALLEST N_t that brings 1 kg/s of 37 °C water to 25 °C at 20 °C / 60 %,
  dry, 3 m/s, 6 rows, 2 m (solver `ntu_field`; P2 gets 31 with LMTD). Walls: P1 / P2 0.5 mm,
  P3 1 mm (Re_c about 3800 at P3, below the 4000 mark but not laminar; with 0.5 / 1 mm walls at 2.5 m/s P3 would need 66
  tubes at laminar Re_c about 1700). T_c,o is not monotone in N_t because of the transition
  band. The Re = 5000 rule (61 / 43 / 26 tubes, `COMPARE_RE_TUBE`) is only the reference of the tubes and wall studies.
  Duty and T_c,o are RESULTS, not targets. Module = elementary rows scaled by N_t;
  2 kg/s = 2 modules.
- **Operation**: two fan modes. Fixed speed `COMPARE_W_FR_DESIGN` = 3 m/s (the speed the tube count is sized for and the
  top of the variable range; above about 3 m/s droplet carry-over, noise and fin damage
  become a risk, Towler's typical value is 2.5), or variable speed between `COMPARE_W_FR_MIN` 0.4 m/s and the design
  speed. `operation.control_point`: per hour the SMALLEST face velocity in [0.4, 3] m/s with
  T_c,o <= 25 °C (brentq on w_fr, about 10 solves), dry first, then precooled, else precooled at 3 m/s (`limited` =
  target missed). Precooling pad: eta_B 0.8, 50 Pa. Freezing risk flag `freeze` =
  T_c,o < `COMPARE_T_FREEZE_RISK` (2 °C, an assumption).
- **Definitions**: `pinch` = T_c,o - the AMBIENT air entering the whole cooler (pad and coil), also when precooled (it
  can be negative: the pad cools the water below the dry-bulb ambient); `pinch_coil` =
  T_c,o - the air entering the coil (equal to `pinch` when dry; the coil's own approach). The target is met when T_c,o
  <= 25 °C + `config.TARGET_TOL` (0.01 K): one definition for the row flag `target_met`,
  the controller's search and the hour counts (the controlled hours land within about 5 mK of 25 °C), so met + limited =
  8760 for every cooler.

## Print-size analyses

**`geometry`** (`plot_geometry.plot_print`): one row per cooler -- the whole tube bank to one scale (all tubes, width b
and depth), the marked 3 x 3-tube detail magnified with ONE zoom window for all three
(s_t below, s_l on the right), the fin section (one size and scale for all: d, wall, s_f, delta_f); the header of each
row gives the tube count and the metal masses (`evaluate.calc_masses`: straight tube and fin
only, the same as in the cost model).

**`pinch_fan`** (`plots_print.py`): the year at the fixed design speed in both modes (`operation.year_modes`), the
face-velocity sweep and the controlled year (`year_control`) give
`pinch_year` (pinch of the whole cooler and T_c,o against T_amb, median per 2 K bin), `pinch_year_clouds` (every hour of
the year as a point, dry and precooled, one column per cooler; the lines in these two
are medians per 2 K bin, those of `controlled_fan` are means), `missed_hours` (T_c,o against T_amb per cooler at the
fixed design speed in the adopted operation -- dry, precooled only where dry misses -- the missed
hours in soft red and, as an inset, the duration curve of the excess over the target; legend with the hours and their
share of the year), `pinch_speed` (pinch against w_fr at the median, the design state = P80 and
P95), `controlled_fan` (w_fr linear and P_fan logarithmic against the AMBIENT temperature, T_c,o held at the target; dry
= a line, precooled = every hour as a point plus an average line per cooler; the bins where
most hours miss the target are dotted). The control switches from dry (up to 3 m/s at about 20 °C) to precooled, where
the pad lets the fan drop to a lower speed again. Results (P1 / P2 / P3): the dry outlet reaches
25 °C at T_amb about 20.2 / 20.3 / 20.1 °C, precooled at about 28.8 / 28.9 / 28.5 °C; pinch at the warmest hour 0.1 K
dry, about -11 K precooled. At fixed 3 m/s with the pad the same 375 / 361 / 399 hours miss
the target as the controlled fan is limited at (identical hour sets); dry only it would be 1719 / 1702 / 1752 h (up to
11.8 K).

**Fan curve** (`fan_curve.py`): the solved control collapses onto monotone curves w_fr(T_a,i) and P_fan(T_a,i) per
cooler, separately for dry and precooled operation (PCHIP through binned medians; the schedule
per hour needs only T_amb and the wet-bulb, T_a,i = T - eta_B (T - T_wb)). The annual fan energy from the curves is
within -3 ... +2 % of the solved control; valid for the study module and weather only.

**Annual economics** (`pinch_fan`): three operating modes -- the dry cooler at the fixed 3 m/s, the dry cooler with the
adaptive fan (smallest sufficient speed, no pad: `year_control(allow_pad=False)`) and the wet cooler
(the adaptive fan with the pad when dry is not enough). `annual_costs` (grouped by mode, coolers side by side): energy
(fan, pump) and cost (electricity, water, material write-off), with the pad water [m3/a]
as a small graph; `monthly_cost` (wet adaptive fan); `water_vs_hours` (hours above the target against the water used);
`fan_cop` (heat removed per electrical watt against T_amb, per cooler and mode: the total
Q/(P_fan+P_pump) darkest, the fan alone Q/P_fan and the pump alone Q/P_pump lighter; at -9 °C the dry fan COP is
8000-9000 (P2, P3) / 1700 (P1), the total 7600 / 2900 / 660 for P3 / P2 / P1; the pump COP is a flat
line because the pump power is constant). Stacked components use colours no cooler uses (pump light grey, pad water
purple, material teal). Numbers (P1 / P2 / P3): fixed 3 m/s dry 7959 / 5606 / 7737 kWh fan,
2538 / 1747 / 2428 EUR/a; adaptive dry 2723 / 1592 / 2102 kWh, 967 / 543 / 737 EUR/a; wet adaptive 2569 / 1639 / 2032
kWh, 67 / 85 / 105 m3 water, 1044 / 714 / 910 EUR/a. Hours above the target: 1719 / 1702 / 1752
(both dry modes) against 375 / 361 / 399 (wet). The dry adaptive fan is the cheapest mode but misses the target in about
1700 h; the pad buys those hours for about 170-370 EUR/a. `operation.fixed_year` (the
fixed-speed rule "dry if it meets the target, else the pad") is tested against `year_fixed`.

**`variation`** (`variation.py`, `plots_print.variation_figures`): one-at-a-time geometry variations of ALL coolers --
N_r (1-10), N_t (5-100), l (0.25-5 m), s_l (10-40 mm), s_f (1-10 mm), delta_f (0.05-0.4 mm); the
ranges are wider than the literature ones on purpose (`variation.RANGES`). Per variant: the pinch at the design state
(20 °C / 60 %, dry, fixed 3 m/s) and the wet adaptive-fan year (`operation.control_point` every 6th
hour, `variation.STRIDE`; within 1 % of the full year for the standard coolers, 957 / 660 / 849 EUR/a at the standard
points) -> fan / pump energy, pad water, hours at the limit, metal masses (CSV
`variation_variation.csv`; the costs are computed at plot time, `variation.add_costs`). Two figures: `variation_costs`
(pinch of P1-P3 and the cost of the sensitivity cooler split into fan / pump / water / material /
total) and `variation_totals` (the totals of P1-P3). Rows: pinch above cost per parameter; the cooler's own value is the
large marker; a variant that misses the target for more than `plots_print.MISSED_H` = 500 h
is dotted (a cooler that is too small looks cheap, so cost comparisons there are not like for like). The coolers are
solved in parallel processes (one per cooler; the NTU tables of all row counts are built or loaded
first in the main process). The wiggle of the N_t curves is the tube Reynolds number crossing the transition region
(Nu_c falls to the laminar 5.0 between Re_c about 3800 and 2300 for P3): a property of the
correlation, not an error. Not a cost optimisation: the missed hours are not priced.

## Cost model

Prices (`parameters.py`): electricity `ELECTRICITY_EUR_PER_KWH` 0.2726 (SWM M-Strom, net, valid from 15 Sep 2026; the
125.80 EUR/a basic charge is not included), water 1.84 EUR/m3 net (SWM M-Wasser, from 1 Jul 2026,
including the "Wassercent"; gross 1.97), no wastewater fee by default (evaporated pad water is normally not charged),
copper 13.10465 and aluminium 2.72593 EUR/kg (Westmetall LME cash of 9 Oct 2026: Cu 14,689.00 /
Al 3,055.50 USD/t at 1.12090 USD/EUR; raw metal, no fabrication premium), `LIFETIME_YEARS` = 20. The material write-off
is straight-line over the lifetime without interest or scrap, with the masses of
`evaluate.calc_masses` (`operation.material_writeoff`, column `cost_material`). Water is tracked per hour (`m_dot_ev`,
VDI M8 humidity-ratio increase times the dry-air flow).

All cost columns are computed in ONE place, `operation.add_costs` (from `annual_summary`), at plot time: changing a
price or the cost model (construction cost per area, annuity with interest, a fixed part, the fan with
its inverter) needs only an edit there or in `parameters.py` and `python main.py compare pinch_fan --plots-only`
(seconds, no new solves). The LaTeX block for the thesis appendix (equations and a table of
the parameters with their citations) is `annual_costs_appendix.tex` next to the figures; update it when a price or the
cost model changes.

## Reference results

Run of 12 Oct 2026 with the fixed counts 34 / 32 / 30, 6 rows, 2 m, 3 m/s design speed. Design point 20 °C / 60 % dry:
all three give Q = 50-51 kW and T_c,o about 24.9 °C by construction (P1 891 W fan, 23 kg Al + 38 kg Cu;
P2 622 W, 46 + 49 kg; P3 861 W, 55 + 155 kg). Fixed 3 m/s over the year: target met in 8385 / 8399 / 8361 of 8760 h (P1
/ P2 / P3), 0 h with freezing risk. Fan control: fan 2568 / 1639 / 2031 kWh/a, pump
412 / 97 / 10 kWh/a (P1's few tubes make the water fast), water 67 / 85 / 105 m3/a, annual cost including the 20-year
material write-off 1044 / 714 / 909 EUR/a (P2 cheapest); 361-399 h per year stay above
25 °C even at 3 m/s with the pad. Sensitivity (P2 base): rows (+8.3 K at 1 row) and face velocity (+8 K at 0.5 m/s)
matter most, then the fin pitch; the wall hardly (< 0.6 K). Tube count (`tubes.py`): Q rises with
N_t down to Re about 3000, then falls (transition band) while fan power and mass keep growing. Wall (`wall.py`, 0.3-2
mm): a thicker wall gives a few % more Q (N_t fixed) at a multiple of the copper mass.
The `s_t` sweep of the sensitivity study is noisy (N_t integer rounding changes the face area).

Planned but not built: two-parameter maps and a global sensitivity.
