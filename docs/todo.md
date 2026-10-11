# Open points

## Optional TODOs

- **Validate the air side against Wang's heat-transfer correlation.**
  Wang, Chi & Chang (Part II, the same paper as the friction factor in
  `pressure_drop.py`) also gives a j-factor correlation for plain plate
  fins. Comparing its air-side alpha with VDI M1 Gl. (16)/(18) is the best
  available independent check of k without measured data -- especially
  since M1 is used outside its range here (A/A_p0 = 37.4 vs. 5-30, see
  below). Deferred for thesis scope.

- Done: Cell 3D vs. 2D (solvers `cell_2d` / `cell_3d`); running only a chosen subset of the solvers (`SOLVERS`).
- Dropped: a design-space optimization mode (tried and removed 2026-10, see cooler-comparison.md).
- Planned, not built: two-parameter maps and a global sensitivity (cooler comparison).
- Possible: calendar years and a decade of climate data -- needs the DWD historical archives of station 03379 (climate_data/README.md).

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
