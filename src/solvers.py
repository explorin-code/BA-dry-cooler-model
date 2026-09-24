"""
solvers.py
==========
Three solvers (LMTD, NTU, Cell) for the same dry-cooler problem, converging
to the same (T_coolant_out, T_air_out, Q). LMTD and NTU share an outer
relaxation loop (_relax_lmtd_ntu) -- see its docstring for why that's a
convenience, not a design coupling between the two. solve_scenario() runs
all three against one OperatingConditions and bundles the results -- pure
computation, no printing/plotting (see run_scenario.py for that).
"""

from dataclasses import dataclass
from math import log, pi, exp
import numpy as np

from src.operating_conditions import get_operating_conditions
from src.dry_cooler_physics import get_geometry
from src.fluid_properties import get_fluid_properties, get_air_properties
from src.heat_transfer_core import calc_overall_k, calc_diagnostics
from src.solver_settings import get_solver_settings


# =============================================================================
# Result types
# =============================================================================

@dataclass
class SolverResult:
    k: float
    dQ: float
    T_coolant_out: float
    T_air_out: float
    history_hot: list
    history_cold: list
    history_T_coolant: list
    history_T_air: list
    diagnostics: dict
    # AI-REVIEW: populated only by Cell (see solve_it_cell) -- these are the
    # solver's own final converged grids, not reconstructions. LMTD/NTU leave
    # these None; their spatial profiles are reconstructed on demand by
    # analysis.py instead, since they never computed spatial data at all.
    # See CLAUDE.md.
    T_c_grid: object = None
    T_a_grid: object = None
    k_grid: object = None
    dQdL_grid: object = None
    local_diagnostics: dict = None   # Re_air/Nu_air/Re_coolant/Nu_coolant grids -- captured, not yet plotted


@dataclass
class ScenarioResult:
    lmtd: SolverResult
    ntu: SolverResult
    cell: SolverResult


# =============================================================================
# Small math helpers
# =============================================================================

def calc_LMTD(dT_hot: float, dT_cold: float) -> float:
    """Log-mean temperature difference, with guards against the two usual
    singularities (temperature crossover, dT_hot == dT_cold)."""
    if dT_hot <= 0 or dT_cold <= 0:
        return 1e-5                            # temperature crossover -- not physical, clamp
    if abs(dT_hot - dT_cold) < 1e-5:
        return dT_hot                          # avoid 0/0 as dT_hot -> dT_cold
    return (dT_hot - dT_cold) / log(dT_hot / dT_cold)


def calc_P(NTU1, R1) -> float:
    """P1 correlation for single-pass effectiveness.
    Source (secondary): Brunner, p. 14, citing Holman."""
    return 1 - exp(((NTU1**0.22) / R1) * (exp(-R1 * NTU1**0.78) - 1))

# Alternate (inactive) P1 formula, VDI Heat Atlas p. 50: return 1 - exp((exp(-R1 * NTU1) - 1) / R1)


# =============================================================================
# Solvers 1 & 2: LMTD-based and NTU/P-based
# =============================================================================
# LMTD and NTU are two DISTINCT algorithms -- LMTD marches the log-mean
# temperature difference, NTU combines a per-row P-effectiveness correlation
# across n_rows. They share an outer relaxation loop (_relax_lmtd_ntu) only
# because that loop turned out identical between them, not because they're
# secretly the same method. Each solver's own physics lives in its own
# top-level step function (_lmtd_step / _ntu_step) below, readable on its
# own. If one solver ever needs a different relaxation scheme, give it its
# own loop -- the driver has no solver-specific knowledge and doesn't need
# touching.
# =============================================================================

def _relax_lmtd_ntu(step_fn, tag, omega, ops, geo, settings):
    """Shared outer iteration for LMTD/NTU: initial guesses, under-relaxation
    with the dynamic omega switch, and history bookkeeping. step_fn (see
    _lmtd_step/_ntu_step) supplies the one thing that actually differs
    between the two solvers -- how (k, dT_hot, dT_cold) become the raw
    (T_coolant_out, T_air_out, dQ) for this iteration."""
    dm_coolant = ops.m_dot_1
    dm_air = ops.m_dot_2

    dT_hot_it = settings.dT_hot_it_init
    dT_cold_it = settings.dT_cold_it_init

    T_coolant_out = ops.T_coolant_in - dT_hot_it
    T_air_out = ops.T_air_in + dT_cold_it

    threshold = settings.convergence_threshold # convergence threshold [K]
    diff_hot = 1.0                             # seeded above threshold so the loop runs at least once
    diff_cold = 1.0

    # AI-REVIEW: dynamic omega switch (warm -> requested) triggered by
    # magnitude-growth or double sign-flip, tuned empirically against this
    # project's own settings/geometry, not derived from a stability proof.
    # Should converge to the same fixed point regardless of the omega
    # schedule, but the trigger thresholds themselves are heuristic. See
    # CLAUDE.md.
    # Drops to the requested omega on the first step-size growth or 2nd sign
    # flip (checked below); the sign-flip counters only compare real
    # iterations against each other, never against the seeded first step.
    # NOTE: Cell's own warm-start (solve_it_cell) uses a 1.0 ceiling here,
    # not 0.5 -- unreconciled, not something this refactor changed.
    omega_warm = min(0.5, 10 * omega)
    omega_active = omega_warm
    sign_flips_hot = 0
    sign_flips_cold = 0

    history_hot = []
    history_cold = []
    history_T_coolant = []                     # T_coolant_out at each iteration
    history_T_air = []                         # T_air_out at each iteration

    # Fluid states from the FINAL iteration, used afterwards for the
    # final-iteration diagnostics (Pr/Re/Nu).
    coolant_state = None
    air_state = None

    # while either outlet temperature is still moving significantly
    while (abs(diff_hot) > threshold) or (abs(diff_cold) > threshold):

        dT_hot = dT_hot_it
        dT_cold = dT_cold_it

        # Dynamic mean temperatures for this iteration
        T_air_mean = (T_air_out + ops.T_air_in) / 2.0
        T_coolant_mean = (ops.T_coolant_in + T_coolant_out) / 2.0

        # Fluid states at the mean temperature
        coolant_state = get_fluid_properties(ops, T_coolant_mean, ops.P_coolant)
        air_state = get_air_properties(ops, T_air_mean, ops.P_air)

        # Overall k for this iteration
        k = calc_overall_k(
            geo=geo,
            ops=ops,
            coolant_state=coolant_state,
            air_state=air_state,
            T_air_out=T_air_out,
        )

        raw_T_coolant_out, raw_T_air_out, dQ = step_fn(
            k, dT_hot, dT_cold, coolant_state, air_state, dm_coolant, dm_air, geo, ops
        )

        # Under-relaxation: slow down the update (dynamically chosen omega)
        T_coolant_out = (1 - omega_active) * T_coolant_out + omega_active * raw_T_coolant_out
        T_air_out = (1 - omega_active) * T_air_out + omega_active * raw_T_air_out

        # New diffs for the loop condition
        new_dT_hot = ops.T_coolant_in - T_air_out
        new_dT_cold = T_coolant_out - ops.T_air_in

        diff_hot = new_dT_hot - dT_hot_it
        diff_cold = new_dT_cold - dT_cold_it

        dT_hot_it = new_dT_hot
        dT_cold_it = new_dT_cold

        # Live progress printout -- useful for spotting slow/diverging runs
        print(f"{tag:<7}Hot Error: {diff_hot:.5f} | Cold Error: {diff_cold:.5f} | T_coolant_out: {T_coolant_out:.2f} | T_air_out: {T_air_out:.2f} | omega: {omega_active}")

        # Drop to the requested omega for good at the first sign of trouble:
        # the step size growing instead of shrinking, or the error flipping
        # sign twice (oscillating back and forth), whichever trips first.
        if history_hot:
            if diff_hot * history_hot[-1] < 0:
                sign_flips_hot += 1
            if diff_cold * history_cold[-1] < 0:
                sign_flips_cold += 1

            if (abs(diff_hot) >= abs(history_hot[-1]) or abs(diff_cold) >= abs(history_cold[-1])
                    or sign_flips_hot >= 2 or sign_flips_cold >= 2):
                omega_active = omega

        # error and temperature appended together -> same index = same iteration
        history_hot.append(diff_hot)
        history_cold.append(diff_cold)
        history_T_coolant.append(T_coolant_out)
        history_T_air.append(T_air_out)

    # --- Final-iteration diagnostics ---------------------------------------
    # Pr / Re / Nu on both sides, evaluated at the converged states/T_air_out
    # (mirrors exactly what the last calc_overall_k call inside the loop used).
    diagnostics = calc_diagnostics(
        geo=geo,
        ops=ops,
        coolant_state=coolant_state,
        air_state=air_state,
        T_air_out=T_air_out,
    )

    return (k, dQ, T_coolant_out, T_air_out,
            history_hot, history_cold, history_T_coolant, history_T_air,
            diagnostics)


def _lmtd_step(k, dT_hot, dT_cold, coolant_state, air_state, dm_coolant, dm_air, geo, ops):
    """LMTD solver's physics core: log-mean temperature difference drives
    dQ across the whole array, which implies new outlet temperatures."""
    LMTD = calc_LMTD(dT_hot, dT_cold)
    A_tot = geo.A * geo.n_tubes * geo.n_rows   # total outer area across the whole array
    dQ = A_tot * k * LMTD

    dT_air_rise = dQ / (dm_air * air_state.cp)
    dT_coolant_drop = dQ / (dm_coolant * coolant_state.cp)

    raw_T_coolant_out = ops.T_coolant_in - dT_coolant_drop
    raw_T_air_out = ops.T_air_in + dT_air_rise
    return raw_T_coolant_out, raw_T_air_out, dQ


def _ntu_step(k, dT_hot, dT_cold, coolant_state, air_state, dm_coolant, dm_air, geo, ops):
    """NTU/P-effectiveness solver's physics core: a single row's P1
    correlation, combined (P-series) across n_rows into an overall P1tot/
    P2tot that maps inlet temperatures directly to outlet temperatures.
    NOTE: the row-combination below is VDI C1 §5.1's formula for n
    DISCRETE heat exchangers coupled in series, applied here per-row as an
    approximation -- it is not the book's separate, more complex closed-form
    for a single crossflow bundle with n rows (Table 5, p.50, which uses a
    recursive term), which is not what's implemented."""
    n = geo.n_rows
    A_run = geo.A * geo.n_tubes                # outer area of a single row (one tube pass)

    # Dimensionless groups
    W1 = dm_coolant * coolant_state.cp
    W2 = dm_air * air_state.cp
    R1 = W1 / W2
    NTU1 = (k * A_run) / W1

    P1 = calc_P(NTU1, R1)

    # --- P-series: combine per-row P1 into overall P1tot across n rows ---
    if abs(R1 - 1.0) < 1e-6:
        # Source: VDI Wärmeatlas, Chapter C1, §5.1, Gl. (47), p.54.
        # Confirmed 2026-08 against the 12th ed. (previously cited as
        # p.53, which was off by one page).
        P1tot = (n * P1) / (1 + (n - 1) * P1)
    else:
        # Source (secondary): Brunner, p. 14, Tab. 2.2, Item 6 (rearranged
        # VDI Wärmeatlas C1 §5.1 Gl. (46), p.54). Confirmed 2026-08: this
        # rearrangement is algebraically exact -- Brunner's X is the
        # reciprocal of the book's own retention term.
        X = (1 - P1 * R1) / (1 - P1)                # retention term
        P1tot = (X**n - 1) / (X**n - R1)            # overall system effectiveness

    P2tot = P1tot * R1

    raw_T_coolant_out = ops.T_coolant_in - P1tot * (ops.T_coolant_in - ops.T_air_in)
    raw_T_air_out = ops.T_air_in + P2tot * (ops.T_coolant_in - ops.T_air_in)

    dQ = W1 * (ops.T_coolant_in - raw_T_coolant_out)
    return raw_T_coolant_out, raw_T_air_out, dQ


def solve_it_LMTD(omega: float = None, ops=None, geo=None, settings=None):
    """omega: under-relaxation factor; use the same value as solve_it_NTU()
    for comparable steps. Delegates its outer loop to _relax_lmtd_ntu();
    see _lmtd_step() for the LMTD-specific physics."""
    if geo is None:
        geo = get_geometry()
    if ops is None:
        ops = get_operating_conditions(geo=geo)
    if settings is None:
        settings = get_solver_settings()
    if omega is None:
        omega = settings.central_omega

    return _relax_lmtd_ntu(_lmtd_step, "[LMTD]", omega, ops, geo, settings)


def solve_it_NTU(omega: float = None, ops=None, geo=None, settings=None):
    """omega: under-relaxation factor; use the same value as solve_it_LMTD()
    for comparable steps. Delegates its outer loop to _relax_lmtd_ntu();
    see _ntu_step() for the NTU-specific physics."""
    if geo is None:
        geo = get_geometry()
    if ops is None:
        ops = get_operating_conditions(geo=geo)
    if settings is None:
        settings = get_solver_settings()
    if omega is None:
        omega = settings.central_omega

    return _relax_lmtd_ntu(_ntu_step, "[NTU]", omega, ops, geo, settings)


# =============================================================================
# Solver 3: Cell-Method
# =============================================================================

def _solve_cell_step(r, t, s, T_c, T_a_for_inlet, ops, geo, n_segments, dm_coolant, dm_air):
    """Shared per-cell computation for both passes of _relax_cell_grid: the
    local NTU/P1 relation for one (row, tube, segment) cell, given the
    coolant/air inlet temperatures for that cell.

    T_a_for_inlet is which air grid to read the neighbor's inlet from --
    the two passes disagree on this and it must stay a parameter, not a
    hardcoded grid: Pass 1 passes T_a_old (this iteration's air side hasn't
    moved yet), Pass 2 passes the live T_a (already updated by Pass 1 this
    same iteration). That ordering is the staggered-mixing scheme itself,
    not an implementation detail."""
    T_c_in = get_coolant_inlet(r, t, s, n_segments, T_c, ops, geo)
    T_a_in = get_staggered_air_inlet(r, t, s, T_a_for_inlet, ops, geo)

    T_mean = (T_c_in + T_a_in) / 2.0
    props_c = get_fluid_properties(ops, T_mean, ops.P_coolant)
    props_a = get_air_properties(ops, T_mean, ops.P_air)
    k_local = calc_overall_k(geo, ops, props_c, props_a, T_a_in)

    A_cell = geo.A / n_segments
    W1 = dm_coolant * props_c.cp
    W2 = dm_air * props_a.cp
    R_loc = W1 / W2
    NTU_loc = (k_local * A_cell) / W1
    P1_loc = calc_P(NTU_loc, R_loc)

    return T_c_in, T_a_in, P1_loc, R_loc


def cell_path_order(n_rows, n_segments):
    """(r, s) coordinates in coolant flow order: r_last -> 0, serpentine
    (s direction alternates per row). This is the traversal _relax_cell_grid's
    Pass 1 below actually walks (per tube) AND the path analysis.py's
    calc_cell_profile reduces the converged grid along for the spatial
    profile plot -- one definition instead of two hand-kept-in-sync copies."""
    path = []
    for r in range(n_rows - 1, -1, -1):
        s_range = range(0, n_segments) if r % 2 == 0 else range(n_segments - 1, -1, -1)
        for s in s_range:
            path.append((r, s))
    return path


def _relax_cell_grid(T_c, T_a, omega, dT_hot_it, dT_cold_it,
                      geo, ops, n_segments, dm_coolant, dm_air,
                      threshold, max_iter, min_iter = 1):
    history_hot = []
    history_cold = []
    history_T_coolant = []
    history_T_air = []

    r_last = geo.n_rows - 1
    r_exit = 0
    s_exit = n_segments - 1

    for iteration in range(max_iter):
        T_c_old = T_c.copy()
        T_a_old = T_a.copy()
        max_raw_cold = 0.0
        max_raw_hot = 0.0

        # Pass 1: coolant direction -- r_last -> 0 (see cell_path_order)
        for r, s in cell_path_order(geo.n_rows, n_segments):
            for t in range(geo.n_tubes):
                T_c_in, T_a_in, P1_loc, R_loc = _solve_cell_step(
                    r, t, s, T_c, T_a_old, ops, geo, n_segments, dm_coolant, dm_air
                )

                raw_T_c = T_c_in - P1_loc * (T_c_in - T_a_in)
                max_raw_cold = max(max_raw_cold, abs(raw_T_c - T_c_old[r, t, s]))
                T_c[r, t, s] = (1 - omega) * T_c_old[r, t, s] + omega * raw_T_c

        # Pass 2: air direction -- 0 -> r_last
        for r in range(geo.n_rows):
            for t in range(geo.n_tubes):

                s_range = range(0, n_segments) if r % 2 == 0 else range(s_exit, -1, -1)
                for s in s_range:
                    T_c_in, T_a_in, P1_loc, R_loc = _solve_cell_step(
                        r, t, s, T_c, T_a, ops, geo, n_segments, dm_coolant, dm_air
                    )
                    P2_loc = P1_loc * R_loc

                    raw_T_a = T_a_in + P2_loc * (T_c_in - T_a_in)
                    max_raw_hot = max(max_raw_hot, abs(raw_T_a - T_a_old[r, t, s]))
                    T_a[r, t, s] = (1 - omega) * T_a_old[r, t, s] + omega * raw_T_a

        # Aggregate scalar outlet temps for this iteration, for the plot
        T_coolant_out = float(np.mean(T_c[r_exit, :, s_exit]))
        T_air_out = float(np.mean(T_a[r_last, :, :]))

        new_dT_hot = ops.T_coolant_in - T_air_out
        new_dT_cold = T_coolant_out - ops.T_air_in
        diff_hot = new_dT_hot - dT_hot_it
        diff_cold = new_dT_cold - dT_cold_it
        dT_hot_it = new_dT_hot
        dT_cold_it = new_dT_cold

        history_hot.append(diff_hot)
        history_cold.append(diff_cold)
        history_T_coolant.append(T_coolant_out)
        history_T_air.append(T_air_out)

        max_err_hot = np.max(np.abs(T_a - T_a_old))    # air grid -- "hot" side error
        max_err_cold = np.max(np.abs(T_c - T_c_old))   # coolant grid -- "cold" side error

        print(f"[Cell] Hot Error: {max_raw_hot:.5f} | Cold Error: {max_raw_cold:.5f} | "
              f"T_coolant_out: {T_coolant_out:.2f} | T_air_out: {T_air_out:.2f} | omega: {omega}")

        if iteration + 1 >= min_iter and max_raw_cold < threshold and max_raw_hot < threshold:
            break

    return (T_c, T_a, dT_hot_it, dT_cold_it, T_coolant_out, T_air_out,
            history_hot, history_cold, history_T_coolant, history_T_air)


def solve_it_cell(n_segments: int = None, omega: float = None, ops=None, geo=None, settings=None):
    """omega: same under-relaxation factor as the other solvers, applied per-cell."""
    # --- 1. Load static configuration -----------------------------------
    if geo is None:
        geo = get_geometry()
    if ops is None:
        ops = get_operating_conditions(geo=geo)
    if settings is None:
        settings = get_solver_settings()
    if omega is None:
        omega = settings.central_omega
    if n_segments is None:
        n_segments = settings.cell_n_segments

    dm_coolant = ops.m_dot_1 / geo.n_tubes
    dm_air = ops.m_dot_2 / (geo.n_tubes * n_segments)

    # --- 2. Grid initialization (row, tube, segment) ----------------------
    T_c = np.full((geo.n_rows, geo.n_tubes, n_segments), ops.T_coolant_in)
    T_a = np.full((geo.n_rows, geo.n_tubes, n_segments), ops.T_air_in)

    dT_hot_it = settings.dT_hot_it_init
    dT_cold_it = settings.dT_cold_it_init

    # AI-REVIEW: two-stage relaxation (omega_warm + raw-residual convergence
    # check in _relax_cell_grid). Empirically fast/stable at the settings
    # used throughout this app, but the same physical problem has been
    # observed converging in ~8 iterations at one n_segments/omega
    # combination and ~1000 at another -- a documented robustness gap, not
    # a known correctness bug. See CLAUDE.md.
    # Stage 1: coarse, fast propagation at a large (capped) omega -- gets
    # the grid close without chasing tight precision at an undamped step.
    # NOTE: LMTD/NTU's warm-start (_relax_lmtd_ntu) uses a 0.5 ceiling here,
    # not 1.0 -- unreconciled, not something this refactor changed.
    omega_warm = min(1.0, 10 * omega)
    (T_c, T_a, dT_hot_it, dT_cold_it, T_coolant_out, T_air_out,
     hist_hot_1, hist_cold_1, hist_Tc_1, hist_Ta_1) = _relax_cell_grid(
        T_c, T_a, omega_warm, dT_hot_it, dT_cold_it,
        geo, ops, n_segments, dm_coolant, dm_air,
        threshold=settings.cell_stage1_threshold,
        max_iter=settings.cell_stage1_max_iter,
        min_iter=settings.cell_stage1_min_iter,
    )

    # Stage 2: fine polish at the requested omega, warm-started from stage 1
    ## in the initial setup the temps dont't change anymore within 0.01 K, but left in here for completeness, as it doesnt take long
    (T_c, T_a, dT_hot_it, dT_cold_it, T_coolant_out, T_air_out,
     hist_hot_2, hist_cold_2, hist_Tc_2, hist_Ta_2) = _relax_cell_grid(
        T_c, T_a, omega, dT_hot_it, dT_cold_it,
        geo, ops, n_segments, dm_coolant, dm_air,
        threshold=settings.cell_stage2_threshold,
        max_iter=settings.cell_stage2_max_iter,
        min_iter=settings.cell_stage2_min_iter,
    )

    history_hot = hist_hot_1 + hist_hot_2
    history_cold = hist_cold_1 + hist_cold_2
    history_T_coolant = hist_Tc_1 + hist_Tc_2
    history_T_air = hist_Ta_1 + hist_Ta_2

    # --- 3. Final-iteration diagnostics -----------------------------------
    T_coolant_mean = (ops.T_coolant_in + T_coolant_out) / 2.0
    T_air_mean = (ops.T_air_in + T_air_out) / 2.0
    coolant_state = get_fluid_properties(ops, T_coolant_mean, ops.P_coolant)
    air_state = get_air_properties(ops, T_air_mean, ops.P_air)

    k = calc_overall_k(
        geo=geo,
        ops=ops,
        coolant_state=coolant_state,
        air_state=air_state,
        T_air_out=T_air_out,
    )
    dQ = ops.m_dot_1 * coolant_state.cp * (ops.T_coolant_in - T_coolant_out)

    diagnostics = calc_diagnostics(
        geo=geo,
        ops=ops,
        coolant_state=coolant_state,
        air_state=air_state,
        T_air_out=T_air_out,
    )

    # --- 4. Per-cell local diagnostics (final converged grid only) --------
    # AI-REVIEW: reconstruction/interpretation layer, computed once after
    # convergence -- NOT tracked live during iteration. Reuses the exact
    # same calc_overall_k/calc_diagnostics calls (and local-T_a_in-as-
    # T_air_out convention) _relax_cell_grid already uses per cell, just
    # re-evaluated once on the final T_c/T_a rather than discarded each
    # sweep. See CLAUDE.md.
    A_cell = geo.A / n_segments
    dL_cell = geo.height / n_segments

    k_grid = np.zeros((geo.n_rows, geo.n_tubes, n_segments))
    dQdL_grid = np.zeros((geo.n_rows, geo.n_tubes, n_segments))
    Re_air_grid = np.zeros((geo.n_rows, geo.n_tubes, n_segments))
    Nu_air_grid = np.zeros((geo.n_rows, geo.n_tubes, n_segments))
    Re_coolant_grid = np.zeros((geo.n_rows, geo.n_tubes, n_segments))
    Nu_coolant_grid = np.zeros((geo.n_rows, geo.n_tubes, n_segments))

    for r in range(geo.n_rows):
        for t in range(geo.n_tubes):
            for s in range(n_segments):
                T_c_in = get_coolant_inlet(r, t, s, n_segments, T_c, ops, geo)
                T_a_in = get_staggered_air_inlet(r, t, s, T_a, ops, geo)
                T_mean = (T_c_in + T_a_in) / 2.0
                props_c = get_fluid_properties(ops, T_mean, ops.P_coolant)
                props_a = get_air_properties(ops, T_mean, ops.P_air)

                k_cell = calc_overall_k(geo, ops, props_c, props_a, T_a_in)
                diag_cell = calc_diagnostics(geo, ops, props_c, props_a, T_a_in)

                k_grid[r, t, s] = k_cell
                dQdL_grid[r, t, s] = k_cell * A_cell * (T_c[r, t, s] - T_a[r, t, s]) / dL_cell
                Re_air_grid[r, t, s] = diag_cell['Re_air']
                Nu_air_grid[r, t, s] = diag_cell['Nu_air']
                Re_coolant_grid[r, t, s] = diag_cell['Re_coolant']
                Nu_coolant_grid[r, t, s] = diag_cell['Nu_coolant']

    return (k, dQ, T_coolant_out, T_air_out,
            history_hot, history_cold, history_T_coolant, history_T_air,
            diagnostics,
            T_c, T_a, k_grid, dQdL_grid,
            {'Re_air': Re_air_grid, 'Nu_air': Nu_air_grid,
             'Re_coolant': Re_coolant_grid, 'Nu_coolant': Nu_coolant_grid})


# =============================================================================
# Scenario orchestrator -- runs all three, no printing/plotting
# =============================================================================

def solve_scenario(ops, geo=None, omega: float = None, n_segments: int = None, settings=None) -> ScenarioResult:
    """Runs LMTD/NTU/Cell against the same ops/geo/omega. Pure computation:
    no printing, no plotting -- see run_scenario.plot_scenario() for that."""
    if geo is None:
        geo = get_geometry()
    if settings is None:
        settings = get_solver_settings()
    if omega is None:
        omega = settings.central_omega
    if n_segments is None:
        n_segments = settings.cell_n_segments

    lmtd = SolverResult(*solve_it_LMTD(omega=omega, ops=ops, geo=geo, settings=settings))
    ntu = SolverResult(*solve_it_NTU(omega=omega, ops=ops, geo=geo, settings=settings))
    cell = SolverResult(*solve_it_cell(n_segments=n_segments, omega=omega, ops=ops, geo=geo, settings=settings))

    return ScenarioResult(lmtd=lmtd, ntu=ntu, cell=cell)


def get_coolant_inlet(r, t, s, n_segments, T_c, ops, geo):
    # Coolant: logic depends on circuit layout -- for now simple serpentines in z-direction
    r_last = geo.n_rows - 1

    if r % 2 == 0:
        s_prev = s - 1

        if s == 0:
            if r == r_last:
                return ops.T_coolant_in
            else:
                return T_c[r+1, t, s]
        else:
            return T_c[r, t, s_prev]
    else:
        s_prev = s + 1

        if s == n_segments - 1:
            if r == r_last:
                return ops.T_coolant_in
            else:
                return T_c[r+1, t, s]
        else:
            return T_c[r, t, s_prev]

def get_staggered_air_inlet(r, t, s, T_a, ops, geo):
    # Air: Staggered mixing logic (Source: VDI C1, 3.1)

    if r == 0:
        return ops.T_air_in
    else:
        if r % 2 == 0:
            if t == 0:
                return T_a[r-1, t, s]
            else:
                return (T_a[r-1, t, s] + T_a[r-1, t-1, s]) / 2.0
        else:
            if t == geo.n_tubes - 1:
                return T_a[r-1, t, s]
            else:
                return (T_a[r-1, t, s] + T_a[r-1, t+1, s]) / 2.0
