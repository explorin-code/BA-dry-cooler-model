# CLAUDE.md

Project context and instructions for Claude Code sessions in this repo.

## Needs independent review

The items below were implemented with heavy AI assistance across an
extended session (Aug 2026) and have NOT yet been independently verified
against source material or reference data. Treat their numeric output as
unconfirmed until checked off.

- **`fluid_properties.py`** — humid-air property functions (`_hapropsSI_state`
  etc.). CoolProp's `HAPropsSI` key conventions (`Vha`/`Cha` = per kg
  *humid* air, not per kg dry air) were never independently verified
  against a reference psychrometric chart/table.
- **`precooling.py`** — `calc_cooling_limit` and `calc_precooler`'s
  temperature/humidity interpolation formulas. **Resolved 2026-08**:
  checked directly against the physical VDI Wärmeatlas (12th ed.) —
  `calc_cooling_limit`'s residual is an exact match to Chapter M8, §2.2,
  Gl. (15), p.1778; the interpolation formulas match Gl. (13), same page.
  Still open: the underlying Lewis-factor=1 assumption (book's Gl. 14)
  wasn't separately re-derived, just inherited as the book's own stated
  prerequisite; and neither formula has been checked against an actual
  measured wet-bulb reference point (e.g. 20°C/30% RH ambient) — the
  citations confirm the *formula*, not that this specific implementation
  reproduces real hardware behavior.
- **`pressure_drop.py`**:
  - Coolant side: the laminar/turbulent blend in `calc_f_coolant`
    (linear interpolation across Re 2320-3000) is the author's own choice,
    not from VDI -- the book gives only the critical Re=2320 and a
    qualitative "may still be laminar up to ~8000" caveat.
  - Air side: **replaced 2026-09** -- the VDI L1.4 bundle correlation was
    swapped for the Wang, Chi & Chang (Part II) plain-fin friction factor
    + Kays & London core relation (`calc_f_wang`,
    `calc_delta_p_platefin_air`, orchestrated by `calc_delta_p_bundle_air`).
    Formulas were transcribed with AI help and not yet checked against the
    papers. `w_e`/`A_e`/sigma all come from `Geometry.sigma` = A_e/A_fr. The old VDI code was deleted (git history only).
  - Wang is used slightly outside its tested range at the current
    geometry: s_l = 30 mm (Wang's P_l, range 12.4-27.5 mm); `check_wang_range` warns.
    At 1.5 mm fin spacing (`FIN_PITCH_MM` = 1.62, in range): bundle ΔP
    ~104 Pa (+50 Pa pad when the air passes it, see `PAD_IN_DRY_AIR_PATH`);
    fan ~152 W dry with bypass inlet, ~227 W through the pad (forced draft).
  - Post-processing basis (2026-09): air ΔP/fan power use the **Cell**
    result's T_a_o (same reference as the precooling decider).
    Re_Dc/f at the mean air temperature (G_e fixed by continuity),
    Kays & London density terms with rho_a_i/rho_a_o. Effect vs. the old
    inlet-only evaluation: ~+2.4% bundle ΔP.
- **`dry_cooler_physics.py`** — fin count `Geometry.N_f` = int(l/s_f): the
  original hardcoded "9 fins/inch" value was replaced with the value
  derived from the fin pitch after finding it numerically inconsistent
  with spacing/thickness. Confirm which one actually matches the physical
  hardware -- Kröger (Vol. 2, Sec. 8.1) gives 300-450 fins/m as typical;
  the current 617 fins/m is above that.
  **2026-09**: the fin spacing (then input `FIN_SPACING`; since 2026-10 the
  input is `FIN_PITCH_MM` = spacing + thickness, in mm) changed from 0.00023 to 0.0015 m (the old
  value gave a 0.35 mm pitch, ~72 FPI, and ~1440 Pa air ΔP). The new
  value gives ~15.7 FPI, still not the original 9 FPI (~0.0027 m) --
  confirm against the hardware.
- **`solvers.py`**:
  - Cell relaxation (**resolved 2026-10**): with the old shared omega = 0.2
    the iteration count depended erratically on n_segments (86 / 70 / 8 / 8
    at 5 / 10 / 20 / 50) and Cell stopped *before* converging (Q 12.715
    instead of the tight-tolerance fixed point 12.693 kW, +0.17 %): the
    per-cell max-raw-change < 1e-3 K criterion doesn't bound the error
    accumulated along the coolant path when steps are damped. Cell now has
    its own `CELL_OMEGA = 1.0`: 9 iterations at every resolution 2-100,
    wall time linear in n_segments. The former two-stage relaxation was
    collapsed into one stage -- results bit-identical. Remaining caveat:
    any change-per-iteration criterion stops early at omega < 1.
  - **Convergence criterion (2026-10, common to all three)**: every solver
    stops when both outlet temperatures change by less than
    `CONVERGENCE_THRESHOLD` (1e-3 K) per iteration -- LMTD/NTU's outer
    loop as before, Cell switched from the per-cell max change to the
    outlets. Per-solver overrides only as hardcoded constants at the top of
    solvers.py (`LMTD_THRESHOLD`, `NTU_THRESHOLD`, `CELL_THRESHOLD`,
    `NTU_FIELD_THRESHOLD`; None = shared; NTU's inner field defaults to
    1/10 of NTU's threshold). `CELL_MIN_ITER` removed. Errors vs. fully
    converged at 1e-3: LMTD Q -1.4e-3 % (damped, error ~3x the step), NTU
    ~1e-4 %, Cell ~3e-3 %; NTU-vs-Cell gap 0.3729 % converged vs. ~0.37 %.
    Floor: the property cache rounds T to 0.01 K
    (`fluid_properties.CACHE_T_DECIMALS` = 2). A state on a rounding boundary
    flips between two property sets: a period-2 limit cycle, ~1.5e-5 K for
    Cell, but up to ~1.5e-3 K for LMTD/NTU at cold, humid hours -- found by
    the annual benchmark: at the 1e-3 K threshold 4 of 2000 hours (LMTD) and
    9 % (NTU at omega 1) never converged. Fix (2026-10): the shared outer loop
    halves omega on a *stalled* oscillation (sign flip without the amplitude
    shrinking below `OSC_STALL_RATIO` = 0.7 of two steps back, twice in a row;
    floor `OMEGA_MIN` = 0.01). Converging oscillations are left alone: normal
    operating points are bit-identical. Synthetic year: 0 failures; damping
    active at 256 h (LMTD, up to 20 it) / 665 h (NTU, up to 14 it); Q vs. a
    0.001 K cache: max 2.2e-2 % (LMTD), 8e-3 % (NTU) = the rounding itself.
    Alternative not taken: 0.001 K rounding (also 0 failures) costs ~x1.6-3.2
    CoolProp evaluations and Cell ~+31 % in annual runs (hit rate 99.3 ->
    97.7 %). Before this fix, the old warm-start rule only dropped omega once
    to the requested value -- a no-op for NTU at omega 1.
  - `solve_it_LMTD`/`solve_it_NTU`'s dynamic omega switch (`omega_warm` ->
    requested `omega`, triggered by step-size growth or a double sign
    flip) -- as of the Aug 2026 structural refactor this logic lives in the
    `_relax_lmtd_ntu` driver both solvers call, not duplicated inline in
    each. Tuned empirically against this project's own settings, not
    derived from a stability proof -- should converge to the same fixed
    point regardless of the omega schedule, but the trigger thresholds
    (grow-by-any-amount, exactly 2 sign flips) are heuristic.
- **`economics.py`** — pump/fan power via `ṁ·ΔP/(η·ρ)` with fixed
  efficiencies (`eta_pump` 0.9, Towler Sec. 20.7; `eta_fan` 0.7, Towler
  Sec. 19.16). `calc_fan_power`'s
  `fan_position` ('forced' = inlet density, default; 'induced' = outlet
  density, ~+3%) -- actual fan arrangement not yet decided; nothing
  passes it yet.
- **`solvers.py` -- NTU** (**replaced 2026-10**): the old P-series NTU
  (VDI C1 series combination of a Holman per-row P, via Brunner) was removed -- it's in
  git history (last in commit `9c66769`). NTU is now the element-wise
  method of Cabezas-Gomez, Navarro & Saiz-Jabardo (2007), J. Heat Transfer
  129, 282-289, for the counter-cross-flow arrangement G^c_(n,1)
  (`calc_element_effectiveness` Eq. 1, `calc_element_outlets` Eqs. 3-4,
  `solve_ntu_field`, `_ntu_step`).
  - Validated: `python -m src.validation_ntu` runs `solve_ntu_field`
    against Table 1 Eqs. (9), (14)-(16) (1-4 passes, both C_min sides,
    max error ~1e-5) and 102 Table 5 points (5-10 passes, max error
    ~8e-5, i.e. table rounding). The closed-form transcription itself is
    checked against Table 6.
  - Implementation assumptions (not from the paper): air is passed
    row-to-row at the same tube/axial position with no neighbour
    averaging (reproduces the closed forms); all tube circuits are
    identical, so one circuit is solved and broadcast to N_t.
  - **Circle back if results look off**: the paper assumes constant U and
    cp. k/cp are held constant within each field solve and updated
    between field solves at the mean temperatures by `_relax_lmtd_ntu`
    (same as LMTD) -- no locally varying k like Cell.
  - **Shared field traversal (2026-10)**: NTU and Cell iterate their field
    with the same `sweep_two_way` (solvers.py section 2): pass 1 coolant
    along `cell_path_order` with last iteration's air, pass 2 air row 0 ->
    last with live air. Both use the same element relation
    `calc_element_outlets` (Eqs. 3/4) and connections (`_coolant_inlet`,
    `_air_inlet`); they differ only in where Gamma/B come from (NTU one
    constant pair per field solve, Cell per cell from local k/cp) and NTU's
    single tube column (no neighbour averaging). The paper's own
    coolant-direction-only march (Fig. 2a, item 1.7) is kept as the
    INACTIVE `sweep_coolant_direction` -- swap the function name in
    `solve_ntu_field` / `_relax_cell_grid` to compare. Same converged
    field (one-way reproduces the old NTU bit for bit; NTU result change
    from switching ~5e-6 %, Cell identical to 1e-14 K); field sweeps
    one-way/two-way ~2x at 6 rows, ~3.4x at 24, Cell iterations 19/9 at
    6 rows, 117/36 at 24. At the same per-sweep threshold one-way leaves a
    larger remaining error (slower convergence). Cost per sweep: the
    generic sweep is ~2.8x the old inline NTU loop (two element
    evaluations per element in two-way, function-call overhead) -- 84 vs.
    235 us/sweep. With omega 1 and the common 1e-3 K threshold (field 1e-4)
    NTU is ~13 ms warm (start of day 36 ms); Cell 8 iterations, ~390 ms. NTU's energy-balance check (Q_c vs. Q_a,
    1e-6 relative) was removed 2026-10: it only tripped on slightly
    unconverged fields; NTU's Q is now the coolant side, like Cell's.
  - **Relaxation study (2026-10)**: (a) -- applied: NTU now has its own
    `NTU_OMEGA` = 1.0 (`settings.ntu_omega`; solve_scenario's omega is
    LMTD-only) and the start rule is max(omega, min(0.5, 10 omega)), so
    LMTD is unchanged. LMTD's `CENTRAL_OMEGA` 0.2 -> 0.3 (2026-10): scan
    over ambient / winter -10 °C / precooled / hot 33 °C -- stable up to
    0.4, 0.5 diverges in winter, >= 0.7 everywhere (the undamped LMTD step
    overshoots, kA/C ~ 1.6). 0.2 stopped too early at the 1e-3 K threshold
    (hot point Q -0.1 %); 0.3: 13 instead of 17 iterations, errors
    <= 5e-3 %. Outer loop now capped (`OUTER_MAX_ITER` = 500, raises; it had
    no cap and looped forever on divergence -- a strongly diverging run can
    also stop earlier in CoolProp, e.g. coolant below 0 °C); Cell likewise raises at
    `CELL_MAX_ITER` instead of silently returning an unconverged field). Before: NTU never ran at
    `CENTRAL_OMEGA` = 0.2: `_relax_lmtd_ntu` starts at min(0.5, 10 omega)
    = 0.5 and only drops on a growing step / double sign flip, which NTU's
    monotone convergence never triggers -> effectively 0.5, 14 outer
    iterations. Fixed omega 1.0: 5 outer iterations (4 at 24 rows), 64 ->
    23 ms, Q unchanged (3 mW). LMTD shares the loop and keeps its damping.
    Net vs. start of day (warm cache): NTU 36 -> 23 ms (two-way sweep alone
    would have been 66 ms), Cell unchanged (~440 ms). (b) omega != 1 per element inside
    a sweep is wrong: each pass is an exact sequential solve along its
    fluid, so per-element relaxation compounds along the path (> 1
    diverges -- now caught by `_field_changes`; < 1 stops early).
    (c) Block over-relaxation (relax each fluid's whole field after its
    pass): best ~1.1 at 6 rows (72 -> 45 sweeps), ~1.3 at 24 rows (258 ->
    128); geometry-dependent, small absolute gain.
  - **Warm start (2026-10)**: each field solve starts from the previous
    outer iteration's field (`solve_ntu_field(..., T_start=...)`, passed by
    `_ntu_step` via field_out); the first one starts from the inlet
    temperatures as in the paper (validation always does). Sweeps per field
    9, 11, 11, 11, 11 -> 9, 9, 6, 3, 1 (53 -> 28); NTU ~13 -> ~7 ms warm,
    ~8.5 ms without cache; Q change 4e-5 %. Bigger gain for deeper coolers
    (24 rows 174 -> 80 ms). Reading for the thesis: one continuous field
    iteration in which k/cp are updated every few sweeps; with constant
    properties it reduces to the paper's single field solve.
  - **Two NTU modes (2026-10)**, `NTU_MODE` / `solve_it_NTU(mode=...)`:
    'field' (default) solves the element field in every outer iteration;
    'table' looks up the characteristic P(NTU_a, R) (NTU_a = kA_tot/C_a,
    R = C_a/C_c), computed once per (N_r, N_e) with the same element field in
    dimensionless form (C_a = 1, C_c = 1/R, T_c,i = 1, T_a,i = 0 -> P =
    T_a,o; paper p. 283: (eps, NTU) relations are independent of the inlet
    temperatures; Tables 4/5 are the same kind of table). Grid 48 NTU_a
    (geomspace 0.05-6) x 40 R (0.05-2, R/N_e <= 0.1), bicubic
    RectBivariateSpline (scipy), stored as cache/ntu_tables/*.npz
    (gitignored), outside the range -> field solve with a warning. Only N_r
    and N_e matter, every other geometry change reuses the table.
    Measured: build 5.9 s (once), load 4 ms; P interpolation error max
    1.3e-5, mean 8e-7 (300 random points); Q vs. field mode -3e-5 %.
    Per operating point: field 8.6 / 8.1 / 6.8 ms (no / cold / warm cache),
    table 1.6 / 1.3 / 0.04 ms (profile=False), 4.5 / 3.9 / 2.7 ms with
    profile=True (one field solve at the end for the grids -- the default,
    so plots work; pass profile=False for multi-run use, grids are None).
    Tipping point with an empty disk cache ~900 operating points (5.9 s /
    ~6.8 ms), with a stored table from the first point. A year (8760 h):
    ~14 s instead of ~75 s (estimate from the per-point times).
  - Elements (`NTU_N_ELEMENTS`, default 20) are independent of Cell's
    `CELL_N_SEGMENTS`. Eq. (1) needs C_c^e << C_h^e (warns above 0.1);
    project result changes <0.01% between 5 and 50 elements.
  - NTU lands close to LMTD (6 passes is near pure counterflow); current
    ambient results: LMTD 12.55, NTU 12.50, Cell 12.45 kW.
- **`solvers.py` -- Cell vs. NTU gap** (2026-10): Cell evaluated both
  fluids' properties at the average of coolant and air inlet temperature
  (coolant ~4-5 K too cold, k too low). **Fixed**: each fluid now at its
  own local mean temperature (`_cell_local_states`), and Cell reports the
  mean of the local k it used (previously a diagnostic at overall means).
  **Per-cell effectiveness replaced (2026-10)**: Cell used `calc_P`
  (Holman's empirical crossflow fit, via Brunner), which at Cell's
  per-cell R ~ 19 gave P ~17% too low (0.0100 vs. an exact fine-grid
  0.0120). Removed; Cell now uses the same Cabezas-Gomez element relation
  as NTU (Eqs. 1/4, P1 = 2B/(2+B)), so both cite one source. Holman is no
  longer used anywhere. Current results (ambient, after all 2026-10 fixes
  incl. D6 conductivities, cooling-limit fix and Cell's local T_a_m):
  Cell Q = 12.45 kW (NTU 12.50 = +0.37 %, LMTD 12.55 = +0.8 %); precooled:
  Cell 17.56 kW (NTU 17.73 = +0.96 %, LMTD 17.79 = +1.35 %).
  **Mean air temperature for k (2026-10)**: `calc_overall_k`/
  `calc_diagnostics` take T_a_m directly (thesis Eq. mean_temperatures):
  LMTD/NTU the cooler mean (T_a_i + T_a_o)/2, Cell the cell mean -- the same
  temperature Cell already used for the air properties. Previously Cell
  averaged the cooler inlet with the cell inlet for w_e_m; switching to the
  cell mean gave Cell k +0.2 %, Q +0.08 % (precooled +0.10 %).
  Remaining Cell-vs-NTU differences are Cell's own modelling choices --
  local k and staggered neighbour averaging of air -- both small here.
- **`solvers.py` -- LMTD**: deliberately kept as the simplest method,
  pure counterflow (no correction factor F). The study compares three
  methods of increasing complexity (LMTD < NTU element field < Cell);
  keep LMTD simple. Known deviation: VDI C1 Eq. (33) gives F ~ 0.987 for
  this gegensinnig 6-pass arrangement; LMTD is ~0.5% above NTU at the
  current operating point (VDI: pure-counterflow equations only valid from
  ~20 passes at NTU > 5). If F is ever added, use it only as
  Q = F*k*A*dT_lm, not via Eqs. (25)/(26), which would turn LMTD into an
  effectiveness (NTU) formulation.
- **`analysis.py`** (spatial profiles for the LMTD/NTU/Cell comparison
  plot) — Cell's (`calc_cell_profile`) and NTU's (`calc_ntu_profile`)
  profiles are genuine reductions of the grids their solvers computed.
  LMTD's is NOT: `calc_lmtd_profile` reconstructs a closed-form
  exponential T(x) profile purely from LMTD's converged scalars --
  `solve_it_LMTD` never computes this internally. Treat the LMTD curve as
  an interpretation of what LMTD *implies*. The local `dQ/dL = k × local
  ΔT × (area per unit length)` used for LMTD and Cell follows Baehr,
  Thermodynamik (2016), Sec. 3.1 (dQ = k dA ΔT); NTU's dQ/dL is the
  element heat duty divided by element length.
- **`pressure_drop.py`** (continued) — `calc_delta_p_coolant`: **fixed
  2026-10** -- friction used the total tube length (all passes) per pass,
  i.e. counted N_r times too often (10.7 instead of 2.5 kPa, pump 3.35
  instead of 0.79 W); now L = one pass, as in Towler 19.8 / thesis
  Eq. 2.13. Bend/contraction/expansion allowance 1.5*(2N_p-1)/N_p velocity
  heads per pass (Towler 19.8). Viscosity correction still uses
  `mu_w = mu` (no wall-temperature model).
  `calc_delta_p_coolant` evaluates coolant properties at the inlet
  temperature only (the air side now uses in/out/mean states).

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

## Naming conventions (renamed 2026-10, aligned with the thesis notation)

- Temperatures: `T_<medium>_<state>` in °C (`T_c_i`, `T_a_o`, `T_K`,
  `T_a_pc`, grids `T_c_grid`/`T_a_grid`, `Delta_T_hot`/`Delta_T_cold`).
  theta/vartheta dropped entirely (2026-10), in code, figures and thesis.
  Kelvin only where absolute temperatures are needed, marked `_K` /
  `T_kelvin` (`T_m_K`, `T_i_K` in `calc_w_e_m`). The Cabezas-Gomez element
  functions use T_hot/T_cold (paper: T_h/T_c) so T_c always means coolant.
- Indices: c coolant, a air, i in, o out, m mean, pc precooled, K cooling
  limit, da dry air, ev evaporation, f fin, p pipe, fr front/face,
  e narrowest cross-section. Humidity ratio: `Y` (ops.Y, Y_i, Y_pc, Y_K).
- **Notation scheme (option A, 2026-10)**: one letter per kind of quantity,
  the index says which; lengths lower case (ISO 80000-3 / DIN 1304),
  P = power, p = pressure, s = pitch, delta = thickness, d = diameter.
- Geometry: `d` (tube outer diameter, no index -- o = outflow), `d_i`,
  `d_c` (collar diameter, d + 2 delta_c), `delta_f`/`delta_c` (fin / collar
  thickness), `s_t`/`s_l`/`s_f` (transverse / longitudinal tube pitch, fin
  pitch; Wang's P_t/P_l/F_p/D_c), clear fin spacing has no symbol (inline
  s_f - delta_f), `N_t`/`N_r` (tubes per row / rows = passes; N_r kept as
  "rows", not N_l), `N_f` (fins per pass), `l` (one pass), `l_tot` (all
  passes), `b` (width), `A_f`/`A_p`/`A_p0` (fin / exposed pipe / bare pipe
  area per tube pass; A_p = pi d (l - N_f delta_f)), `A_tot`, `A_fr`,
  `A_e`, `sigma` (= A_e/A_fr), `A_cs_c`, `lambda_f`/`lambda_p`, `delta_p`
  (tube wall thickness), `b_f`/`l_f` (fin rectangle, staggered, VDI M1
  Eq. 14 -- moved from the fin-efficiency function into Geometry 2026-10).
- Flows/powers: `w_c`, `w_fr`, `m_dot_c`, `m_dot_a`, `V_dot_c`,
  `V_dot_a`, `p_c`/`p_a` (pressures), `Q_dot`, `C_dot_c`/`C_dot_a`
  (capacity rates), `P_pump`/`P_fan` (powers), `m_dot_ev`, `Delta_p_*`.
- Heat-transfer chain (renamed 2026-10, side index a = air, c = coolant):
  `w_e_m` (narrowest-gap velocity at the mean air temperature; VDI w_eT),
  `Re_a`/`Pr_a`/`Nu_a` (VDI M1, row factor `C_r`), `alpha_a` (VDI alpha_R),
  fin efficiency `eta_f` via `phi_f_0`/`phi_f` (from Geometry.b_f/l_f) (VDI phi -- index
  f so it is not the relative humidity), `X`, `alpha_a_eff` (incl. fin
  efficiency; VDI "scheinbar" alpha_S); `Re_c`/`Pr_c`/`Nu_c` (VDI G1,
  `calc_Nu_c_lam`/`calc_Nu_c_turb`, blend factor gamma, friction psi),
  `alpha_c` (VDI alpha_i), coolant-side area `A_c` (was A_i); k with
  1/k = 1/alpha_a_eff + A/A_c (1/alpha_c + delta_p/lambda_p),
  delta_p = Geometry.delta_p. Wang's friction Re is `Re_d_c` (length d_c).
  Diagnostics dict keys follow (`Re_a`, `Nu_c`, `alpha_c`, ...).
- Deliberately NOT renamed: `mu`, and the Cabezas-Gomez paper notation
  inside the NTU element functions (T_hot/T_cold, i.e. coolant/air;
  B there is the paper's B, not a width).
- Fan and pump are spelled out to keep f = fin and p = pipe unambiguous:
  `eta_fan`/`P_fan`, `eta_pump`/`P_pump` (thesis Eqs. 2.6-2.8 write
  eta_f/W_f and eta_p/W_p -- to be aligned there).

## Open notation decisions (thesis "Notation" chapter, 2026-10)

Keep in mind -- not yet decided / deliberately left for now:
- ~~T vs. theta~~ -- resolved: T everywhere (°C unless absolute K needed).
- ~~Index i = inflow and inner~~ -- mostly resolved: alpha_i -> alpha_c,
  A_i -> A_c; only d_i (inner diameter) keeps i.
- ~~B = cooler width vs. Cabezas-Gomez B~~ -- resolved: the width is now `b`.

- **Friction factors (thesis, 2026-10)**: the thesis writes ALL friction
  factors as f -- coolant f_c (Darcy: 64/Re, Blasius, blend), air f_a
  (Wang, Fanning). Not psi (so the notation's psi entry is obsolete). The
  code's `psi` (Gnielinski/Konakov in `calc_Nu_c_turb`) and `f_c`/`f` in
  pressure_drop.py are internal names; keep thesis and code in mind when
  writing.
- **Temperatures**: thesis uses T everywhere (°C unless absolute K is
  needed); vartheta dropped.

## Run modes

Dependencies: `requirements.txt` (`venv/bin/pip install -r requirements.txt`);
scipy is needed for NTU table mode and the annual benchmark, streamlit and
pandas only for the GUI. Rendered figures (`thesis_figures/out/`) and the NTU
table cache (`cache/`) are gitignored.

Independent toggles in `parameters.py`, each overridable per run with a
flag (`python main.py --[no-]insight --[no-]plots --[no-]convergence
--[no-]benchmark --[no-]resolution --[no-]cell-2d`):
- `INSIGHT_MODE`: per-iteration progress + full result boxes in the
  terminal (logging at DEBUG). Without it only the per-scenario summary
  (INFO: Q, outlets, pinch, deviation from Cell, iterations, wall time).
  Project code logs via `logging.getLogger(__name__)`, no `print`
  (except the standalone `validation_ntu.py`).
- `PLOT_RESULTS`: per scenario a cooler-results figure (inputs, economics,
  each solver's result box next to its spatial profile).
- `PLOT_CONVERGENCE`: per scenario outlet temperatures + iteration errors.
- `BENCHMARK_MODE`: caching benchmark (`benchmark_solvers.py`, fresh silent
  runs on the ambient operating point) + solver-performance figure
  (`plot_solver.py`): solvers compared per cache mode no/cold/warm, and
  each solver's time split. ~20 s (Cell without cache alone ~12 s).
- `RESOLUTION_MODE`: resolution sweep over `BENCHMARK_RESOLUTIONS` (normal
  runs use `CELL_N_SEGMENTS` / `NTU_N_ELEMENTS`) + its figure. ~15 s.
- `CELL_2D` (`--[no-]cell-2d`): Cell solves one representative tube
  instead of all N_t (output broadcast to all tubes). With uniform inlet
  air 2D and 3D are identical to machine precision (all tubes behave the
  same, so the neighbour averaging of air has no effect); 2D is ~17x
  faster with a warm cache, ~2.5x on a cold single run (same CoolProp
  calls). 3D only matters for non-uniform inlet air / per-tube conditions.
  Passed to the solvers via `SolverSettings.cell_2d`.
- `ANNUAL_MODE` (`--[no-]annual`, `benchmark_annual.py`): annual benchmark --
  LMTD, NTU (field), NTU (table), Cell (3D) and Cell (2D) over 8760 hourly operating points
  (only inlet air T/phi change; coolant inlet and flows from parameters.py),
  dry operation only. Weather: `ANNUAL_WEATHER_CSV` (T_air_C, phi_percent)
  or a synthetic Munich-like year (`synthetic_year`, seeded: seasonal
  cosine 9 +- 9.5 °C, daily cycle, AR(1) weather swings; an assumption, not
  measured data). Property cache on, emptied per solver; NTU table rebuilt
  (build time shown separately). A solver over `ANNUAL_TIME_BUDGET_S`
  (120 s) is stopped and extrapolated linearly; hours are processed in the
  same seeded random order for all solvers so a stopped run covers the whole
  year. Figure: cumulative time with the table build and the tipping point,
  log-scale totals, the year's T_a,i and Q. Several minutes.
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
(exploded, 1 and 2 segments along the tube; the row gap is a parameter), plus `cells_bend` / `cells_bend_real`:
2+2 cells at the tube ends with the U-bends at the back, rows pulled apart /
touching. Qualitative on purpose -- no
numbers, colours from a real Cell run, colour bar only cold/hot. Fins are
drawn as a few plates per cell, continuous across a row, not in the gap
between pulled-apart rows (not the real pitch); in the model they
only enter via A and the fin efficiency in k. Nothing in the model imports
this file (could later be rendered automatically per run).
Pinch point = T_c_o - T_a_i (cold end), shown per solver with
LMTD/NTU deviations from Cell.

## Inputs in mm / % and the GUI

- Geometry lengths are entered in **mm** in `parameters.py` (`*_MM`, e.g.
  `D_TUBE_OUTER_MM`, `FIN_PITCH_MM`, `HEIGHT_MM`) and converted to SI
  metres on read: `from_parameters` accepts `(CONSTANT_NAME, MM)` and
  divides by 1000. The model itself only ever sees metres.
- Relative humidity likewise in **%** (`PHI_AIR_PERCENT`, `(name, PERCENT)`,
  divided by 100); the model works with the fraction `ops.phi`, output shows %.
- Coolant: `COOLANT_TYPE` is any CoolProp fluid name (default 'Water');
  the GUI offers `fluid_properties.COOLANTS` as a dropdown (water, ethylene /
  propylene glycol 20-50 % as CoolProp INCOMP mixtures, checked 5-60 °C).
  Note: with water the coolant side is at Re ~3400, inside VDI G1's
  heat-transfer transition band 2300-4000 (Nu blended linearly between the
  laminar value at 2300 and the turbulent one at 4000; the friction factor
  uses its own 2320-3000 blend); 30 % glycol drops it to ~1700 (laminar),
  alpha_c 1584 -> 383 W/m²K, Q -37 %. Precooler water stays
  'Water' (precooling.calc_c_W) regardless of the coolant.
- `gui.py` (Streamlit, `venv/bin/streamlit run gui.py`) is the interactive
  alternative to `main.py`. It builds its form from `parameters.py`
  (values = defaults, inline comments = help, `# --- Group: ... ---` =
  groups), sets edited values on the `src.parameters` module for the
  session only (never writes the file), runs `run_scenarios`, and shows
  tables, log output and all figures. One run at a time across browser
  tabs (shared module state). Restart the Streamlit server after editing
  anything under `src/` -- it only reloads `gui.py` itself.
- The benchmark's temporary patches (uncached lookups, counters) always
  restore the originals captured at import, so an interrupted or
  overlapping run can't leave the property cache broken.

## Optional TODOs

- **Validate the air side against Wang's heat-transfer correlation.**
  Wang, Chi & Chang (Part II, the same paper as the friction factor in
  `pressure_drop.py`) also gives a j-factor correlation for plain plate
  fins. Comparing its air-side alpha with VDI M1 Gl. (16)/(18) is the best
  available independent check of k without measured data -- especially
  since M1 is used outside its range here (A/A_p0 = 37.4 vs. 5-30, see
  below). Deferred for thesis scope.

- **Cell: 3D vs. 2D** -- done as the `CELL_2D` toggle (see Run modes).
- **Solver toggles**: run only a chosen subset of LMTD/NTU/Cell (e.g. only
  NTU for fast annual runs).
- **Optimization mode**: compare/optimize geometries (pitches, fin
  spacing, rows, ...) against Q, pressure drop/fan power, etc.

## Known placeholders / not yet implemented

These are explicitly marked in the code (`TEMP:`, `not yet implemented`,
fixed-value comments) rather than silently wrong -- listed here so they're
easy to find in one place too.

- `precooling.calc_eta_B()` — hardcoded `0.8`. Meant to be empirical and
  velocity-dependent; the plan (stated early in the design discussion) was
  to implement several swappable approximations, never done.
- `operating_conditions.DEFAULT_PHI_AIR_FALLBACK` — hardcoded `0.30`.
  Needs a real Munich-climate design value (annual average vs. a summer
  design-day value was left as an open question).
- `pressure_drop.calc_delta_p_pad()` — fixed `50` Pa (ASHRAE 2020
  Ch. 41 range for a 300 mm pad: 34.8-74.6 Pa), always added when
  precooling; in dry operation only if `PAD_IN_DRY_AIR_PATH = True`
  (False = separate bypass inlet, as some designs have). Same value wet
  or dry for now. Needs a velocity-dependent correlation
  (e.g. S. He et al., per the pad source material).
- No fan-curve / actual-operating-point matching. `calc_fan_power` assumes
  the fan delivers exactly the design `w_o` regardless of back-pressure;
  the earlier design discussion concluded this needs an outer loop
  matching the fan's curve against the system's ΔP-vs-flow curve, not yet
  built.
- No dollar-cost layer. `economics.py` currently reports physical
  quantities only (W, kg/s) -- nothing converts these into actual
  cost figures (electricity/water price).
- `SolverResult.local_diagnostics` (per-cell `Re_a`/`Nu_a`/
  `Re_c`/`Nu_c` grids from Cell's diagnostic pass) is captured
  but not yet visualized anywhere.
- `Geometry.sigma` (= A_e/A_fr, `dry_cooler_physics.py`) -- **source resolved
  2026-10**: VDI M1, p. 1689, worked example ("Verengter
  Stroemungsquerschnitt"), circular-fin formula adapted to continuous plate
  fins (D = t_q). The same example confirms that M1's Re_d uses the
  narrowest-gap velocity with the temperature correction w_eT (here w_e_m), as
  implemented. Still open: only the transverse gap is checked (diagonal
  gap of staggered banks ignored; fine at current geometry).
- `calc_Nu_a` (VDI M1 Gl. 16/18) is used outside its stated range:
  A/A_p0 = 37.4 vs. 5-30 (Re_d = 1812 is inside 1e3-1e5). The correlation
  was fitted for circular-finned tubes; applying it to continuous plate
  fins is a modelling choice. Wang's j-factor correlation (same paper as
  the friction factor) would be the independent check -- deferred.
  `calc_Nu_a`'s row-count prefactor (0.33/0.36/0.38) and
  `calc_Nu_c_turb`/`calc_Nu_c_lam` (`heat_transfer_core.py`) were
  **resolved 2026-08** — all three confirmed exact matches to VDI
  Wärmeatlas Chapter M1 §4 (Gl. 16/18) and Chapter G1 §3.2.2/§4.1
  (Gl. 25/26) respectively; see their docstrings for exact citations.
