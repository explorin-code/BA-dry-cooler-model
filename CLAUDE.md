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
    geometry: s_2 = 30 mm (range 12.4-27.5 mm); `check_wang_range` warns.
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
    at another. Not a correctness bug, but a robustness gap worth
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
- **`analysis.py`** (spatial profile reconstruction for the LMTD/NTU/Cell
  comparison plot) — Cell's profile (`calc_cell_profile`) is a genuine
  reduction of data the solver already computed; LMTD's and NTU's are NOT.
  `calc_lmtd_profile` reconstructs a closed-form exponential T(x) profile
  purely from LMTD's converged scalars -- `solve_it_LMTD` never computes
  this internally. `calc_ntu_profile` reconstructs a row-by-row profile via
  a linear-shooting counterflow march (two trial guesses + one solved
  march, exact for the affine per-row map) and self-checks against NTU's
  own closed-form outlet values, printing a warning if they disagree by
  >0.5°C -- an earlier parallel-flow version of this failed that check and
  was replaced; the current counterflow version passes it, but the
  underlying `dQ/dL = k × local ΔT × (area per unit length)` formula (used
  for all three solvers) still hasn't been checked against an independent
  reference. Treat every curve in the profile-comparison window except
  Cell's as an interpretation of what LMTD/NTU *imply*, not something
  either solver actually solved for.
- **`pressure_drop.py`** (continued) — `calc_delta_p_coolant`'s
  viscosity correction uses `mu_w = mu` (no correction) since no
  wall-temperature model exists yet, and uses a fixed 2.5 velocity-head
  allowance per pass for bends (Towler Eq. 19.20).
  `calc_delta_p_coolant` evaluates coolant properties at the inlet
  temperature only (the air side now uses in/out/mean states).

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
  Ch. 41 range for a 300 mm pad: 34.8-74.6 Pa), added in both scenarios
  since the pad is always installed (wet or dry). Same value wet or dry
  for now. Needs a velocity-dependent correlation
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
- `Geometry.Ao_Ae_ratio` (`dry_cooler_physics.py`) still has a "source not
  recorded" TODO — checked thoroughly through VDI Wärmeatlas Chapter M1
  (fin efficiency + row-bundle Nu) in the 2026-08 citation pass and not
  found there. Best remaining leads, not yet checked: L1.5 (baffled
  bundles) or N4 (Kühltürme — cooling towers, thematically the closest
  match to this whole project). New lead (2026-09): the formula is the
  inverse of Kays & London's sigma = Ac/Afr for a plate-fin core; now
  uses `d_c` instead of `d_a` (identical while `h_collar = 0`).
  `calc_Nu_air`'s row-count prefactor (0.33/0.36/0.38) and
  `calc_Nu_turbulent`/`calc_Nu_laminar` (`heat_transfer_core.py`) were
  **resolved 2026-08** — all three confirmed exact matches to VDI
  Wärmeatlas Chapter M1 §4 (Gl. 16/18) and Chapter G1 §3.2.2/§4.1
  (Gl. 25/26) respectively; see their docstrings for exact citations.
