# Needs independent review

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
> **2026-10: `solvers.py` is now the package `src/solvers/`** (`results`, `common`, `field`, `relaxation`, `lmtd`, `ntu_field`, `ntu_table`, `ntu`,
> `cell`, `scenario`; numbers bit-identical, pinned by `tests/test_solvers_golden.py`). Names below that changed: `_relax_lmtd_ntu` -> `relaxation.relax_lmtd_ntu`,
> `_ntu_step` -> `ntu_field.ntu_step`, `_ntu_table_step` -> `ntu_table.ntu_table_step`, `_coolant_inlet` / `_air_inlet` -> `field.coolant_inlet` / `field.air_inlet`;
> `sweep_two_way` and the inactive `sweep_coolant_direction` are in `field.py`; the four per-solver threshold overrides (`LMTD_THRESHOLD`, ...) were all `None`
> and are gone -- one shared `CONVERGENCE_THRESHOLD`, NTU's inner field at 1/10 of it.

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
  - **Two NTU modes (2026-10)**, `solve_it_NTU(mode=...)` (chosen by the solver key `ntu_field` / `ntu_table`; the `NTU_MODE` parameter was removed 2026-10 -- it was dead):
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
