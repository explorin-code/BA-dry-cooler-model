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
    papers. `w_e`/`A_c`/sigma reuse `Geometry.Ao_Ae_ratio` (taken as
    Afr/Ac). The old VDI code is still in the commented-out block at the
    bottom of the file (unmaintained, has known syntax bugs).
  - Wang is used slightly outside its tested range at the current
    geometry: P_l = 30 mm (range 12.4-27.5 mm); `check_wang_range` warns.
    At FIN_SPACING = 1.5 mm (fin pitch 1.62 mm, in range): bundle ΔP
    ~104 Pa, +50 Pa pad = ~154 Pa total, fan ~225 W (forced draft).
  - Post-processing basis (2026-09): air ΔP/fan power use the **Cell**
    result's T_air_out (same reference as the precooling decider).
    Re_Dc/f at the mean air temperature (G_c fixed by continuity),
    Kays & London density terms with rho_in/rho_out. Effect vs. the old
    inlet-only evaluation: ~+2.4% bundle ΔP.
- **`dry_cooler_physics.py`** — `Geometry.fin_density`: the original
  hardcoded "9 fins/inch" value was replaced with `1/(a+s)` (derived from
  fin spacing/thickness) after finding it numerically inconsistent with
  those two fields. Confirm which one (the original constant, or the
  derived formula) actually matches the physical hardware.
  **2026-09**: `FIN_SPACING` changed from 0.00023 to 0.0015 m (the old
  value gave a 0.35 mm pitch, ~72 FPI, and ~1440 Pa air ΔP). The new
  value gives ~15.7 FPI, still not the original 9 FPI (~0.0027 m) --
  confirm against the hardware.
- **`solvers.py`**:
  - The cell method's two-stage relaxation (`omega_warm` + raw-residual
    convergence check in `_relax_cell_grid`). Empirically fast and stable
    at the settings currently used throughout the app, but a documented
    sensitivity issue exists: the same physical problem converged in ~8
    iterations at one `n_segments`/`omega` combination and ran for ~1000
    at another. The benchmark (2026-10) shows it again: 86 / 70 / 8 / 8
    outer iterations at 5 / 10 / 20 / 50 segments, so Cell is *slower* at
    5-10 segments than at 20. **Cause found (2026-10)**: stage 1 always
    takes 5 iterations; the difference is all in stage 2 (omega = 0.2).
    Worse, at omega = 0.2 Cell stops *before* converging: default settings
    give Q = 12.715 kW, the true fixed point (tight tolerance, any omega)
    is 12.693 kW (+0.17 % error, ~0.02 K at the coolant outlet). The
    per-cell max raw change < 1e-3 K criterion doesn't bound the error
    accumulated along the coolant path when steps are damped. With
    omega = 1 Cell converges correctly in 9 iterations at every resolution
    (time then linear in n_segments). **Fixed**: Cell now uses its own
    `CELL_OMEGA = 1.0` (LMTD/NTU keep `CENTRAL_OMEGA`). Benchmark after
    the fix: 9 iterations at every resolution 2-100, wall time linear in
    n_segments (85 ms at 2, ~0.7 s at 20, 2.5 s at 100). With omega = 1 the
    two-stage relaxation is effectively one stage (both stages at 1.0). Not a correctness bug, but a robustness gap worth
    understanding before relying on it outside today's settings.
  - `solve_it_LMTD`/`solve_it_NTU`'s dynamic omega switch (`omega_warm` ->
    requested `omega`, triggered by step-size growth or a double sign
    flip) -- as of the Aug 2026 structural refactor this logic lives in the
    `_relax_lmtd_ntu` driver both solvers call, not duplicated inline in
    each. Tuned empirically against this project's own settings, not
    derived from a stability proof -- should converge to the same fixed
    point regardless of the omega schedule, but the trigger thresholds
    (grow-by-any-amount, exactly 2 sign flips) are heuristic.
- **`economics.py`** — pump/fan power via `ṁ·ΔP/(η·ρ)` with fixed
  efficiencies (pump 0.9, fan 0.7, Towler2022). `calc_fan_power`'s
  `fan_position` ('forced' = inlet density, default; 'induced' = outlet
  density, ~+3%) -- actual fan arrangement not yet decided; nothing
  passes it yet.
- **`solvers.py` -- NTU** (**replaced 2026-10**): the old P-series NTU
  (VDI C1 series combination of a Holman per-row P, via Brunner) was removed -- it's in
  git history (last in commit `9c66769`). NTU is now the element-wise
  method of Cabezas-Gomez, Navarro & Saiz-Jabardo (2007), J. Heat Transfer
  129, 282-289, for the counter-cross-flow arrangement G^c_(n,1)
  (`calc_element_effectiveness` Eq. 1, `calc_element_outlets` Eqs. 3-4,
  `solve_ntu_field`, `_ntu_step`). Eq. (4) was checked at 300 dpi against
  the paper: its second term is `2B/(2+B)`, no Gamma (an AI-generated
  spec had an extra Gamma, which breaks the energy balance).
  - Validated: `python -m src.validation_ntu` runs `solve_ntu_field`
    against Table 1 Eqs. (9), (14)-(16) (1-4 passes, both C_min sides,
    max error ~1e-5) and 102 Table 5 points (5-10 passes, max error
    ~8e-5, i.e. table rounding). The closed-form transcription itself is
    checked against Table 6.
  - Implementation assumptions (not from the paper): air is passed
    row-to-row at the same tube/axial position with no neighbour
    averaging (reproduces the closed forms); all tube circuits are
    identical, so one circuit is solved and broadcast to n_tubes.
  - **Circle back if results look off**: the paper assumes constant U and
    cp. k/cp are held constant within each field solve and updated
    between field solves at the mean temperatures by `_relax_lmtd_ntu`
    (same as LMTD) -- no locally varying k like Cell.
  - Elements (`NTU_N_ELEMENTS`, default 20) are independent of Cell's
    `CELL_N_SEGMENTS`. Eq. (1) needs C_c^e << C_h^e (warns above 0.1);
    project result changes <0.01% between 5 and 50 elements.
  - First result (ambient): NTU Q = 12.75 kW, LMTD 12.81 kW, Cell 11.32
    kW (old P-series NTU: 12.55 kW). NTU is now close to LMTD (6 passes
    is near pure counterflow).
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
  longer used anywhere. Result (ambient, after the CELL_OMEGA fix): Cell
  Q = 12.69 kW (NTU 12.75 = +0.5 %, LMTD 12.81 = +0.9 %); precooled: Cell
  17.85 kW (NTU 18.05 = +1.1 %, LMTD 18.12 = +1.6 %).
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
  an interpretation of what LMTD *implies*. The `dQ/dL = k × local ΔT ×
  (area per unit length)` formula used for LMTD and Cell hasn't been
  checked against an independent reference (NTU's dQ/dL is the element
  heat duty divided by element length).
- **`pressure_drop.py`** (continued) — `calc_delta_p_coolant`'s
  viscosity correction uses `mu_w = mu` (no correction) since no
  wall-temperature model exists yet, and uses a fixed 2.5 velocity-head
  allowance per pass for bends (Towler Eq. 19.20).
  `calc_delta_p_coolant` evaluates coolant properties at the inlet
  temperature only (the air side now uses in/out/mean states).

## Run modes

Independent toggles in `parameters.py`, each overridable per run with a
flag (`python main.py --[no-]insight --[no-]plots --[no-]convergence
--[no-]benchmark --[no-]resolution`):
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
Pinch point = T_coolant_out - T_air_in (cold end), shown per solver with
LMTD/NTU deviations from Cell.

## Optional TODOs

- **Validate the air side against Wang's heat-transfer correlation.**
  Wang, Chi & Chang (Part II, the same paper as the friction factor in
  `pressure_drop.py`) also gives a j-factor correlation for plain plate
  fins. Comparing its air-side alpha with VDI M1 Gl. (16)/(18) is the best
  available independent check of k without measured data -- especially
  since M1 is used outside its range here (A/A_Go = 37.4 vs. 5-30, see
  below). Deferred for thesis scope.

- **Cell: 3D vs. 2D.** Cell solves all `n_tubes` tubes separately, although
  with uniform inlet air nearly all behave identically (NTU solves one
  circuit). That makes Cell up to ~n_tubes x more expensive than needed --
  note this in any performance comparison. Collapsing Cell to 2D (one
  representative tube) would remove that cost but also the ability to
  model non-uniform inlet air or different conditions per tube (edge
  tubes, wall effects). Undecided whether that capability is needed.
- **Solver toggles**: run only a chosen subset of LMTD/NTU/Cell (e.g. only
  NTU for fast annual runs).
- **Optimization mode**: compare/optimize geometries (pitches, fin
  spacing, rows, ...) against Q, pressure drop/fan power, etc.

- **Index naming (for the planned variable/index streamlining):**
  `alpha_i` (coolant-side heat transfer coefficient, "inside" per VDI M1)
  should probably become `alpha_c` (c = coolant), because the subscript
  `i` is already used for "in"/inlet elsewhere. Left as `alpha_i` for now,
  since the indices are not yet consistent across the codebase anyway.

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
- `SolverResult.local_diagnostics` (per-cell `Re_air`/`Nu_air`/
  `Re_coolant`/`Nu_coolant` grids from Cell's diagnostic pass) is captured
  but not yet visualized anywhere.
- `Geometry.Ao_Ae_ratio` (`dry_cooler_physics.py`) -- **source resolved
  2026-10**: VDI M1, p. 1689, worked example ("Verengter
  Stroemungsquerschnitt"), circular-fin formula adapted to continuous plate
  fins (D = t_q). The same example confirms that M1's Re_d uses the
  narrowest-gap velocity with the temperature correction w_eT, as
  implemented. Still open: only the transverse gap is checked (diagonal
  gap of staggered banks ignored; fine at current geometry).
- `calc_Nu_air` (VDI M1 Gl. 16/18) is used outside its stated range:
  A/A_Go = 37.4 vs. 5-30 (Re_d = 1812 is inside 1e3-1e5). The correlation
  was fitted for circular-finned tubes; applying it to continuous plate
  fins is a modelling choice. Wang's j-factor correlation (same paper as
  the friction factor) would be the independent check -- deferred.
  `calc_Nu_air`'s row-count prefactor (0.33/0.36/0.38) and
  `calc_Nu_turbulent`/`calc_Nu_laminar` (`heat_transfer_core.py`) were
  **resolved 2026-08** — all three confirmed exact matches to VDI
  Wärmeatlas Chapter M1 §4 (Gl. 16/18) and Chapter G1 §3.2.2/§4.1
  (Gl. 25/26) respectively; see their docstrings for exact citations.
