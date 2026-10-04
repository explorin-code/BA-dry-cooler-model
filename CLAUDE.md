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
    papers. `w_e`/`A_e`/sigma reuse `Geometry.Afr_Ae_ratio` (taken as
    Afr/Ac). The old VDI code was deleted (git history only).
  - Wang is used slightly outside its tested range at the current
    geometry: P_l = 30 mm (range 12.4-27.5 mm); `check_wang_range` warns.
    At 1.5 mm fin spacing (`FIN_PITCH_MM` = 1.62, in range): bundle ΔP
    ~104 Pa (+50 Pa pad when the air passes it, see `PAD_IN_DRY_AIR_PATH`);
    fan ~152 W dry with bypass inlet, ~227 W through the pad (forced draft).
  - Post-processing basis (2026-09): air ΔP/fan power use the **Cell**
    result's theta_a_o (same reference as the precooling decider).
    Re_Dc/f at the mean air temperature (G_e fixed by continuity),
    Kays & London density terms with rho_a_i/rho_a_o. Effect vs. the old
    inlet-only evaluation: ~+2.4% bundle ΔP.
- **`dry_cooler_physics.py`** — fin density `Geometry.n_R` = 1/F_p: the
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
    collapsed into one stage (`CELL_THRESHOLD`/`CELL_MAX_ITER`/
    `CELL_MIN_ITER`) -- results bit-identical. Remaining caveat: the
    stopping criterion is per-cell; at omega < 1 it would stop early again.
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
    identical, so one circuit is solved and broadcast to N_t.
  - **Circle back if results look off**: the paper assumes constant U and
    cp. k/cp are held constant within each field solve and updated
    between field solves at the mean temperatures by `_relax_lmtd_ntu`
    (same as LMTD) -- no locally varying k like Cell.
  - Elements (`NTU_N_ELEMENTS`, default 20) are independent of Cell's
    `CELL_N_SEGMENTS`. Eq. (1) needs C_c^e << C_h^e (warns above 0.1);
    project result changes <0.01% between 5 and 50 elements.
  - NTU lands close to LMTD (6 passes is near pure counterflow); current
    ambient results: LMTD 12.55, NTU 12.50, Cell 12.44 kW.
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
  incl. D6 conductivities): Cell Q = 12.44 kW (NTU 12.50 = +0.45 %, LMTD
  12.55 = +0.9 %); precooled: Cell 17.51 kW (NTU 17.70 = +1.05 %, LMTD
  17.76 = +1.4 %).
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
  instead of 0.79 W); now L_p = one pass, as in Towler 19.8 / thesis
  Eq. 2.13. Bend/contraction/expansion allowance 1.5*(2N_p-1)/N_p velocity
  heads per pass (Towler 19.8). Viscosity correction still uses
  `mu_w = mu` (no wall-temperature model).
  `calc_delta_p_coolant` evaluates coolant properties at the inlet
  temperature only (the air side now uses in/out/mean states).

## Decisions recorded 2026-10 (from the thesis chapter 2 review)

- Solid conductivities: Cu 380, Al 160 W/(m·K), VDI-Wärmeatlas (2019) D6
  (was uncited pure-metal 401/237; k dropped ~5 %, Q ~2 %). Carbon steel 50
  not yet checked against D6.
- No fin collar is modelled (`h_collar = 0`, D_c = d) -- deliberate.
- Water usage uses the dry-air mass flow, m_dot_a/(1+Y) (thesis Eq. 2.28).
- Cooling limit: c_W from CoolProp (liquid water at theta_K) instead of the
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

- Temperatures in °C: `theta_<medium>_<state>` (`theta_c_i`, `theta_a_o`,
  `theta_K`, `theta_a_pc`); Kelvin: `T_...` (only `T_m`, `T_i` in
  `calc_w_e_T`, `T_kelvin` for CoolProp). Displayed output (figures,
  terminal) still writes T; if a theta is ever shown, use `\vartheta`.
- Indices: c coolant, a air, i in, o out, m mean, pc precooled, K cooling
  limit, da dry air, ev evaporation, f fin, p pipe, fr front/face,
  e narrowest cross-section. Humidity ratio: `Y` (ops.Y, Y_i, Y_pc, Y_K).
- Geometry: `d` (tube outer diameter), `d_i`, `D_c` (Wang collar
  diameter, = d), `delta_f`, `F_p` (fin pitch, input), `t_s` (clear fin
  spacing, derived), `N_r`/`N_t` (rows = passes / tubes per row), `L_p`
  (one pass), `l` (all passes), `A_f`/`A_p`/`A_p0` (fin / exposed pipe /
  bare pipe area per tube and row), `A_tot`, `A_fr`, `A_e`,
  `Afr_Ae_ratio`, `A_cs_c`, `lambda_f`/`lambda_p`.
- Flows/powers: `w_c`, `w_fr`, `m_dot_c`, `m_dot_a`, `V_dot_c`,
  `V_dot_a`, `p_c`/`p_a` (pressures), `Q_dot`, `C_dot_c`/`C_dot_a`
  (capacity rates), `W_pump`/`W_fan`, `m_dot_ev`, `Delta_P_*`.
- Deliberately NOT renamed: `alpha_i`/`alpha_R`/`alpha_S`, all `Re_*`/
  `Nu_*`, `mu`, and the Cabezas-Gomez paper notation inside the NTU
  element functions (T_h/T_c = hot/cold, i.e. coolant/air).
- Fan and pump are spelled out to keep f = fin and p = pipe unambiguous:
  `eta_fan`/`W_fan`, `eta_pump`/`W_pump` (thesis Eqs. 2.6-2.8 write
  eta_f/W_f and eta_p/W_p -- to be aligned there).

## Run modes

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
Pinch point = theta_c_o - theta_a_i (cold end), shown per solver with
LMTD/NTU deviations from Cell.

## Inputs in mm and the GUI

- Geometry lengths are entered in **mm** in `parameters.py` (`*_MM`, e.g.
  `D_TUBE_OUTER_MM`, `FIN_PITCH_MM`, `HEIGHT_MM`) and converted to SI
  metres on read: `from_parameters` accepts `(CONSTANT_NAME, MM)` and
  divides by 1000. The model itself only ever sees metres.
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
- `Geometry.Afr_Ae_ratio` (`dry_cooler_physics.py`) -- **source resolved
  2026-10**: VDI M1, p. 1689, worked example ("Verengter
  Stroemungsquerschnitt"), circular-fin formula adapted to continuous plate
  fins (D = t_q). The same example confirms that M1's Re_d uses the
  narrowest-gap velocity with the temperature correction w_eT, as
  implemented. Still open: only the transverse gap is checked (diagonal
  gap of staggered banks ignored; fine at current geometry).
- `calc_Nu_air` (VDI M1 Gl. 16/18) is used outside its stated range:
  A/A_p0 = 37.4 vs. 5-30 (Re_d = 1812 is inside 1e3-1e5). The correlation
  was fitted for circular-finned tubes; applying it to continuous plate
  fins is a modelling choice. Wang's j-factor correlation (same paper as
  the friction factor) would be the independent check -- deferred.
  `calc_Nu_air`'s row-count prefactor (0.33/0.36/0.38) and
  `calc_Nu_turbulent`/`calc_Nu_laminar` (`heat_transfer_core.py`) were
  **resolved 2026-08** — all three confirmed exact matches to VDI
  Wärmeatlas Chapter M1 §4 (Gl. 16/18) and Chapter G1 §3.2.2/§4.1
  (Gl. 25/26) respectively; see their docstrings for exact citations.
