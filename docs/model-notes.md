# Model notes: sources, assumptions and solver design

What each module is based on, which assumptions it makes, and how the three solvers work. Parameter names are those of
`src/parameters.py`; the thesis chapters 2-3 give the equations.

## 1. Sources and assumptions by module

**`fluid_properties.py`** -- CoolProp (Bell et al. 2014). `HAPropsSI`'s `Vha` / `Cha` are per kg *humid* air, not per kg
dry air; the
enthalpy used for the cooling limit is `Hda` (per kg dry air, VDI's h_1+X). Lookups are cached with T rounded to 0.01 K
(`CACHE_T_DECIMALS`); pure water is clamped to >= 0.01 °C (`WATER_T_MIN`: CoolProp has no liquid water below the melting
point; results below 0 °C are flagged `freeze` and not physical). Coolants: `fluid_properties.COOLANTS` (water, ethylene
/
propylene glycol 20-50 % as CoolProp INCOMP mixtures, checked 5-60 °C); the precooler water is always water.

**`precooling.py`** -- Cooling limit T_K from VDI M8 §2.2 Gl. (15), p. 1778 (12th ed.), root-found with brentq (c_W from
CoolProp);
the temperature / humidity interpolation is Gl. (13). Gl. (15) assumes a Lewis factor of 1 (Gl. 14), inherited from the
book.
T_K at 20 °C / 30 % is 10.84 °C, equal to CoolProp's thermodynamic wet-bulb temperature. The degree of humidification
`calc_eta_B` is the constant 0.8 (ASHRAE 2020, Ch. 41). Water use = dry-air mass flow m_a/(1+Y) times the humidity-ratio
increase (thesis Eq. 2.28).

**`heat_transfer_core.py`** (VDI-Wärmeatlas, 12th ed.)
- Air side: M1 §4 Gl. (16) (staggered, n >= 4 rows, C = 0.38) and Gl. (18) (n = 1-3 rows, C = 0.33 / 0.36); fin
  efficiency M1 Eqs. 7, 8,
  12, 14 with the plate-fin factor phi; `Geometry.sigma` = A_e/A_fr from the worked example on p. 1689, adapted to
  continuous plate
  fins (only the transverse gap is checked; the diagonal gap of staggered banks is ignored). The correlation was fitted
  for
  circular-finned tubes; applying it to continuous plate fins is a modelling choice. Its range (10^3 < Re_d < 10^5, 5 <=
  A/A_p0 <= 30) is
  not enforced; P1-P3 have A/A_p0 = 17.2 / 15.7 / 11.9, but geometry variations can leave it.
- Coolant side: G1 §3.2.2 Gl. (25) (constant heat flux; components Gl. 17, 18, 24) and §4.1 Gl. (26) (psi from Gl. 27).
  Nu is blended
  linearly between the laminar value at Re 2300 and the turbulent one at 4000; with water Re is about 3400, inside the
  band.
- Solid conductivities (VDI D6 Tab. 8, design values of DIN EN ISO 10456): Cu 380, Al 160, steel 50 W/(m K). Densities:
  Al 2700,
  Cu 8960 (VDI), steel 7850 kg/m³.

**`dry_cooler_physics.py`** -- The fin count along one pass is `N_f = int(l / s_f)` (whole fins). At the default fin
pitch this is
about 617 fins/m, above the 300-450 fins/m Kröger (Vol. 2, Sec. 8.1) gives as typical. There is no fin collar by default
(`FIN_COLLAR_THICKNESS_MM` = 0): `delta_c` only enters Wang's collar diameter d_c (friction factor, sigma), not the
heat-transfer
areas, fin efficiency or Nu_a.

**`pressure_drop.py`**
- Coolant: friction over ONE pass length plus 1.5 (2 N_p - 1)/N_p velocity heads per pass for contraction, expansion and
  bends (Towler
  3rd ed., Sec. 19.8, Eq. 19.20); Darcy f = 8 j_f (VDI L1.2). The laminar / turbulent blend of `calc_f_coolant` (linear
  across
  Re 2320-3000) is a choice of this work (the book gives only Re_crit = 2320). The viscosity correction uses mu_w = mu
  (no
  wall-temperature model); coolant properties are taken at the inlet temperature.
- Air: Wang, Chi & Chang (2000, Part II) plain-fin friction factor with the Kays & London core relation (`calc_f_wang`,
  `calc_delta_p_platefin_air`, `calc_delta_p_bundle_air`). Evaluated at the Cell result's T_a_o (the same reference as
  the precooling
  decider), Re at the mean air temperature (G_e fixed by continuity), density terms with rho_a_i / rho_a_o.
  `check_wang_range` warns
  when a geometry leaves Wang's tested range: P3 (d, s_t and s_l) and any cooler with more than 6 rows; P1 and P2 are
  inside.
- Pad: a fixed 50 Pa (ASHRAE 2020 Ch. 41: 34.8-74.6 Pa for a 300 mm pad), always added when precooling, and in dry
  operation only if
  `PAD_IN_DRY_AIR_PATH` (False = separate bypass inlet).

**`economics.py`** -- Power = m_dot dp / (eta rho) with eta_pump = 0.9 (Towler Sec. 20.7) and eta_fan = 0.7 (Towler Sec.
19.16). `calc_fan_power`'s
`fan_position` is `'forced'` (inlet density, default) or `'induced'` (outlet density, about +3 %); nothing passes the
latter. The fan
delivers exactly the requested face velocity: there is no fan-curve matching.

## 2. The solvers (`src/solvers/`)

Package layout: `results`, `common`, `field`, `relaxation`, `lmtd`, `ntu_field`, `ntu_table`, `ntu`, `cell`, `scenario`;
`run_solver(key, ...)`
and `SOLVER_LABELS` / `SOLVER_FAMILY` / `DETAIL_ORDER` (the most detailed selected solver is the reference for
deviations and for the
precooling decision). The numbers are pinned by `tests/test_solvers_golden.py`.

**Shared outer loop** (`relaxation.relax_lmtd_ntu`, used by LMTD and NTU; Cell has its own field loop).
- Convergence: all solvers stop when both outlet temperatures change by less than `CONVERGENCE_THRESHOLD` (1e-3 K) per
  iteration; NTU's
  inner field uses a tenth of it. Errors against a fully converged run: LMTD Q -1.4e-3 %, NTU about 1e-4 %, Cell about
  3e-3 %. Keep the
  threshold >= 1e-4 K: the property cache rounds T to 0.01 K, so a state on a rounding boundary flips between two
  property sets (a
  period-2 limit cycle, about 1.5e-5 K for Cell, up to 1.5e-3 K for LMTD / NTU at cold, humid hours).
- Oscillation damping: omega is halved when an oscillation stalls (sign flip without the amplitude shrinking below
  `OSC_STALL_RATIO` = 0.7
  of two steps back, twice in a row; floor `OMEGA_MIN` = 0.01). Converging oscillations are left alone. A cache rounding
  of 0.001 K would
  also remove the limit cycles but costs Cell about +31 % in annual runs.
- Relaxation: LMTD `CENTRAL_OMEGA` = 0.3 (stable up to 0.4; 0.5 diverges in winter at -10 °C, >= 0.7 everywhere); NTU
  `NTU_OMEGA` = 1.0
  (monotone convergence); Cell `CELL_OMEGA` = 1.0 (9 iterations at every resolution; with omega < 1 the
  change-per-iteration criterion
  stops before convergence, with 0.2 Q was +0.17 % off). The start rule is max(omega, min(0.5, 10 omega)); the loop is
  capped
  (`OUTER_MAX_ITER` = 500, `CELL_MAX_ITER`) and raises instead of looping forever.

**LMTD** -- The whole exchanger as pure counterflow, deliberately without a correction factor F (the study compares
three methods of
increasing complexity). VDI C1 Eq. (33) gives F of about 0.987 for this 6-pass arrangement, so LMTD is about 0.5 % above
NTU. If F is ever
added, use it as Q = F k A dT_lm, not through Eqs. (25) / (26), which would turn LMTD into an effectiveness formulation.

**NTU** -- The element-wise method of Cabezas-Gómez, Navarro & Saiz-Jabardo (2006, J. Heat Transfer 129, 282-290) for
the counter-cross-flow
arrangement G^c_(n,1): `calc_element_effectiveness` (Eq. 1), `calc_element_outlets` (Eqs. 3-4), `solve_ntu_field`,
`ntu_step`.
- Check against the paper: during development the field solver was run against Table 1 Eqs. (9), (14)-(16) (1-4 passes,
  both C_min sides, max
  error about 1e-5) and 102 points of Table 5 (5-10 passes, max error about 8e-5 = table rounding). The script is no
  longer part of the repository
  (`git show 3b6e647:src/validation_ntu.py`).
- Assumptions that are not from the paper: air passes row to row at the same tube / axial position without neighbour
  averaging
  (reproduces the closed forms); all tube circuits are identical, so one circuit is solved and broadcast to N_t. The
  paper assumes constant U and cp: k and cp
  are held constant within a field solve and updated between field solves at the mean temperatures, like LMTD (no local
  k as in Cell).
  Eq. (1) needs C_c^e << C_h^e (a warning above 0.1); the result changes by < 0.01 % between 5 and 50 elements.
- Field traversal (`field.sweep_two_way`, shared with Cell): pass 1 follows the coolant along `cell_path_order` with the
  previous
  iteration's air, pass 2 goes row 0 to last with the live air. Both solvers use the same element relation and
  connections and differ
  only in where Gamma / B come from (NTU: one constant pair per field solve; Cell: per cell from local k / cp). The
  paper's coolant-direction-only march is kept as
  `sweep_coolant_direction` (inactive; swap the function name in `solve_ntu_field` / the Cell loop to compare): same
  converged field,
  but two-way needs about 2x fewer sweeps at 6 rows (3.4x at 24).
- Warm start: each field solve starts from the previous outer iteration's field (the first from the inlet temperatures,
  as in the paper).
- Per-element relaxation inside a sweep is wrong (each pass is an exact sequential solve, so omega != 1 compounds along
  the path; > 1 diverges
  and is caught by `_field_changes`). Block over-relaxation of the whole field gains little (about 1.1 at 6 rows,
  geometry dependent).
- Two modes (`ntu_field`, `ntu_table`): the field mode solves the element field in every outer iteration. The table mode
  looks up the
  characteristic P(NTU_a, R) (NTU_a = kA_tot/C_a, R = C_a/C_c), computed once per (N_r, N_e) with the same field in
  dimensionless form (C_a = 1,
  C_c = 1/R, T_c,i = 1, T_a,i = 0, so P = T_a,o; the (eps, NTU) relations do not depend on the inlet temperatures, paper
  p. 283) and
  interpolated bicubically (`RectBivariateSpline`). Axes: NTU_a `NTU_TABLE_NTU` = geomspace(0.05, 60, 64); R from 0.05
  to 0.1 N_e (Eq. 1 validity),
  20 points per unit R (40 points up to R = 2 at 20 elements). Only N_r and N_e matter, every other geometry change
  reuses the table.
  Tables are stored in `cache/ntu_tables/*.npz` (gitignored; `NTU_TABLE_VERSION` = 2, bump it when the element field
  changes).
  Outside the range the field is solved with a warning (about 100x slower per call). The solver BENCHMARKS build their
  tables over
  the first range only (`NTU_TABLE_NTU_BENCHMARK` = geomspace(0.05, 6, 48)), so the build times reported in the thesis
  keep their meaning
  (a build over the wide axis costs about 3x as much; lookups and results differ by <= 1.5e-6 relative; pinned by
  `test_benchmark_table_is_the_first_table`). Study tables use 60 elements (R <= 6): about 110 s (N_r = 6) to 230 s (N_r
  = 10) per row count, once.
  Interpolation error against the field: max 1.3e-5 in P on the benchmark axis, about 1e-7 up to NTU_a = 55 on the wide
  one.

**Cell** -- The grid of cells with local properties: each fluid is evaluated at its own local mean temperature; k is the
mean of the local
values the cell used; dQ/dl is the cell's energy balance. The per-cell effectiveness is the same Cabezas-Gómez element
relation as NTU
(Eqs. 1 / 4, P1 = 2B/(2+B)). The mean air temperature for k (`calc_overall_k`, `calc_diagnostics`) is the cell mean for
Cell and the cooler
mean (T_a_i + T_a_o)/2 for LMTD / NTU (thesis Eq. mean_temperatures). Air averaging between neighbouring tubes in the
staggered
geometry is a modelling assumption, not a correlation; for uniform inlet air Cell 2D and 3D agree to machine precision
(2D is about 17x
faster with a warm cache). At cold ambient Cell ends 2-3 K warmer than NTU / LMTD (P2 coldest hour 3.6 against 1.0 °C):
local viscosity
and k of the cold water against mean-temperature properties; not investigated further.

**Profiles (`analysis.py`)** -- Cell's and NTU's profiles are reductions of the grids their solvers computed. LMTD's is
a reconstruction: a
closed-form exponential T(x) from its converged scalars (`calc_lmtd_profile`); the cooler-results figure marks it as
such. LMTD's dQ/dl =
k local-dT area-per-length follows Baehr (2016), Sec. 3.1. All three profiles integrate to their solver's Q within 0.05
%
(`tests/test_solvers_golden.py`).

**Typical result** (default run, `python main.py solver --solvers lmtd ntu_field cell_3d`, P2 cell, 20 °C / 60 %): Q =
48.99 (LMTD) /
48.80 (NTU) / 48.75 kW (Cell), pinch 5.28 / 5.32 / 5.34 K; precooled 60.33 / 60.11 / 60.02 kW. The solvers share their
correlations
and element relation, so these are differences between the thermal calculations, not a validation against measurements.

## 3. Known limitations

- Nothing is compared with measurements of a real cooler. The air side could be cross-checked against Wang's j-factor
  correlation (Part II gives
  one next to the friction factor) -- an independent check of k that has not been done.
- Wang's pressure-drop correlation is used outside its tested range by P3 and by every cooler with more than 6 rows or
  s_l above 27.5 mm;
  `check_wang_range` warns. VDI M1's A/A_p0 range is not checked.
- Fixed placeholders: `eta_B` = 0.8, the pad pressure drop 50 Pa, `DEFAULT_PHI_AIR_FALLBACK` = 0.30 (only used if
  neither humidity input is given),
  `mu_w = mu`, the fan efficiency 0.7 / pump 0.9, no fan curve.
- Costs are electricity, pad water and a straight-line write-off of the raw metal (`LIFETIME_YEARS` = 20); there is no
  construction cost model.
- `SolverResult.local_diagnostics` (per-cell Re / Nu grids of Cell) is stored but not drawn.
