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
- **`precooling.py`** — `calc_cooling_limit` (VDI M8 2.2 wet-bulb energy
  balance, root-found via `brentq`) and `calc_precooler`'s
  temperature/humidity interpolation formulas. Verify against a known
  wet-bulb reference point (e.g. 20°C/30% RH ambient).
- **`pressure_drop.py`**:
  - Coolant side: the laminar/turbulent transition blend for `calc_zeta_coolant`
    (linear interpolation across Re 2320-10^4) is a Claude-suggested
    approximation, not something the VDI source itself specifies.
  - Air side: the full staggered-bundle correlation (`calc_zeta_air_staggered`
    and friends) was transcribed from a secondary/paraphrased description of
    VDI L1.4, not the original Wärmeatlas text -- cross-check every equation
    against the actual source before trusting results.
- **`dry_cooler_physics.py`** — `Geometry.fin_density`: the original
  hardcoded "9 fins/inch" value was replaced with `1/(a+s)` (derived from
  fin spacing/thickness) after finding it numerically inconsistent with
  those two fields. Confirm which one (the original constant, or the
  derived formula) actually matches the physical hardware.
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
    flip). Tuned empirically against this project's own settings, not
    derived from a stability proof -- should converge to the same fixed
    point regardless of the omega schedule, but the trigger thresholds
    (grow-by-any-amount, exactly 2 sign flips) are heuristic.
- **`economics.py`** — pump/fan power via `ṁ·v·ΔP`. Reports hydraulic
  power only (no pump/fan efficiency factor); confirm this is the
  intended basis before using the numbers for sizing/cost estimates.
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
- **`pressure_drop.py`** (continued) — `calc_delta_p_bundle_air`'s
  viscosity correction defaults `eta_wall` to the bulk value (`fz_t =
  1.0`, no correction at all) since no wall-temperature model exists yet;
  `calc_delta_p_total`'s `zeta_u` bend coefficient defaults to `1.0`, the
  rough midpoint of VDI's stated 0.5-1.5 range, not a value specific to
  this geometry's actual bend radius.

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
- `pressure_drop.calc_delta_p_pad()` — returns `None`. Needs a
  manufacturer correlation for the adiabatic pad's own pressure drop
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
- Pre-existing, predating this session's changes: `Geometry.Ao_Ae_ratio`,
  `calc_Nu_air`'s row-count prefactor (0.33/0.36/0.38), and
  `calc_Nu_turbulent`/`calc_Nu_laminar` all have "source not recorded"
  TODOs in `heat_transfer_core.py`/`dry_cooler_physics.py`.
