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
- **`solvers.py`** — the cell method's two-stage relaxation (`omega_warm`
  + raw-residual convergence check in `_relax_cell_grid`). Empirically
  fast and stable at the settings currently used throughout the app, but
  a documented sensitivity issue exists: the same physical problem
  converged in ~8 iterations at one `n_segments`/`omega` combination and
  ran for ~1000 at another. Not a correctness bug, but a robustness gap
  worth understanding before relying on it outside today's settings.
- **`economics.py`** — pump/fan power via `ṁ·v·ΔP`. Reports hydraulic
  power only (no pump/fan efficiency factor); confirm this is the
  intended basis before using the numbers for sizing/cost estimates.
