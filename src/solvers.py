"""
solvers.py
==========
Three solvers (LMTD, NTU, Cell) for the same dry-cooler problem, converging
to the same (T_coolant_out, T_air_out, Q). See solve_it_LMTD for the shared
iteration structure. solve_scenario() runs all three against one
OperatingConditions and bundles the results -- pure computation, no
printing/plotting (see run_scenario.py for that).
"""

from dataclasses import dataclass
from math import log, pi, exp
import numpy as np

from src.operating_conditions import get_operating_conditions
from src.dry_cooler_physics import get_geometry
from src.fluid_properties import get_fluid_properties, get_air_properties
from src.heat_transfer_core import calc_overall_k, calc_diagnostics

# Initial guesses for the iteration (both solvers start from the same point)
dT_hot_it_init = 30
dT_cold_it_init = 30


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
    """P1 correlation for single-pass effectiveness (Brunner, p. 14, citing Holman)."""
    return 1 - exp(((NTU1**0.22) / R1) * (exp(-R1 * NTU1**0.78) - 1))

# Alternate (inactive) P1 formula, VDI Heat Atlas p. 50: return 1 - exp((exp(-R1 * NTU1) - 1) / R1)


# =============================================================================
# Solver 1: LMTD-based
# =============================================================================

def solve_it_LMTD(omega: float = 0.1, ops=None, geo=None):
    """omega: under-relaxation factor; use the same value as solve_it_NTU() for comparable steps."""
    # --- 1. Load static configuration -----------------------------------
    if geo is None:
        geo = get_geometry()
    if ops is None:
        ops = get_operating_conditions()

    # Static mass flows from OperatingConditions (already resolved from
    # whichever of w/m/V was specified, evaluated at inlet density).
    dm_coolant = ops.m_coolant
    dm_air = ops.m_o

    # --- 2. Initial guesses for the iterative loop -----------------------
    dT_hot_it = dT_hot_it_init
    dT_cold_it = dT_cold_it_init
    dT_hot = 5.0                               # forced difference to kickstart the while loop
    dT_cold = 5.0

    T_coolant_out = ops.T_coolant_in - dT_hot_it
    T_air_out = ops.T_air_in + dT_cold_it

    threshold = 1e-3                           # convergence threshold [K]
    diff_hot = 1.0                             # seeded above threshold so the loop runs at least once
    diff_cold = 1.0

    # Start aggressive, drop to the requested omega for good the first time
    # the step size grows instead of shrinking, or the error flips sign
    # twice (see the check inside the loop). Sign-flip counters start at 0
    # and only ever compare against real history entries, never the seeded
    # dT_hot/dT_cold above -- so the trivial sign flip on the very first
    # real iteration (seed -> first computed value) never counts as one.
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

        LMTD = calc_LMTD(dT_hot, dT_cold)
        A_tot = geo.A * geo.n_tubes * geo.n_rows   # total outer area across the whole array
        dQ = A_tot * k * LMTD

        # New temperature drops implied by the heat transferred
        dT_air_rise = dQ / (dm_air * air_state.cp)
        dT_coolant_drop = dQ / (dm_coolant * coolant_state.cp)

        raw_T_coolant_out = ops.T_coolant_in - dT_coolant_drop
        raw_T_air_out = ops.T_air_in + dT_air_rise

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
        print(f"[LMTD] Hot Error: {diff_hot:.5f} | Cold Error: {diff_cold:.5f} | T_coolant_out: {T_coolant_out:.2f} | T_air_out: {T_air_out:.2f} | omega: {omega_active}")

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

    # --- 3. Final-iteration diagnostics -----------------------------------
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


# =============================================================================
# Solver 2: NTU / P-based
# =============================================================================

def solve_it_NTU(omega: float = 0.1, ops=None, geo=None):
    """omega: under-relaxation factor; use the same value as solve_it_LMTD() for comparable steps. (Originally omega=0.5 here.)"""
    # --- 1. Load static configuration -----------------------------------
    if geo is None:
        geo = get_geometry()
    if ops is None:
        ops = get_operating_conditions()

    dm_coolant = ops.m_coolant
    dm_air = ops.m_o

    # --- 2. Initial guesses for the iterative loop -----------------------
    dT_hot_it = dT_hot_it_init
    dT_cold_it = dT_cold_it_init
    dT_hot = 5.0                               # forced difference to kickstart the while loop
    dT_cold = 5.0

    T_coolant_out = ops.T_coolant_in - dT_hot_it
    T_air_out = ops.T_air_in + dT_cold_it

    threshold = 1e-3                           # convergence threshold [K]
    diff_hot = 1.0                             # seeded above threshold so the loop runs at least once
    diff_cold = 1.0

    # Start aggressive, drop to the requested omega for good the first time
    # the step size grows instead of shrinking, or the error flips sign
    # twice (see the check inside the loop). Sign-flip counters start at 0
    # and only ever compare against real history entries, never the seeded
    # dT_hot/dT_cold above -- so the trivial sign flip on the very first
    # real iteration (seed -> first computed value) never counts as one.
    omega_warm = min(0.5, 10 * omega)
    omega_active = omega_warm
    sign_flips_hot = 0
    sign_flips_cold = 0

    history_hot = []
    history_cold = []
    history_T_coolant = []                     # T_coolant_out at each iteration
    history_T_air = []                         # T_air_out at each iteration

    coolant_state = None
    air_state = None

    # while either outlet temperature is still moving significantly
    while (abs(diff_hot) > threshold) or (abs(diff_cold) > threshold):
        n = geo.n_rows

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

        A_run = geo.A * geo.n_tubes                # outer area of a single row (one tube pass)

        # Dimensionless groups
        W1 = dm_coolant * coolant_state.cp
        W2 = dm_air * air_state.cp
        R1 = W1 / W2
        NTU1 = (k * A_run) / W1

        P1 = calc_P(NTU1, R1)

        # --- P-series: combine per-row P1 into overall P1tot across n rows ---
        if abs(R1 - 1.0) < 1e-6:
            # Source: VDI Heat Atlas, Section C1, p. 53, Eq. (47)
            P1tot = (n * P1) / (1 + (n - 1) * P1)
        else:
            # Source: Brunner, p. 14, Tab. 2.2, Item 6 (rearranged VDI Eq. 46)
            X = (1 - P1 * R1) / (1 - P1)                # retention term
            P1tot = (X**n - 1) / (X**n - R1)            # overall system effectiveness

        P2tot = P1tot * R1

        raw_T_coolant_out = ops.T_coolant_in - P1tot * (ops.T_coolant_in - ops.T_air_in)
        raw_T_air_out = ops.T_air_in + P2tot * (ops.T_coolant_in - ops.T_air_in)

        dQ = W1 * (ops.T_coolant_in - raw_T_coolant_out)

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
        print(f"[NTU]  Hot Error: {diff_hot:.5f} | Cold Error: {diff_cold:.5f} | T_coolant_out: {T_coolant_out:.2f} | T_air_out: {T_air_out:.2f} | omega: {omega_active}")

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

    # --- 3. Final-iteration diagnostics -----------------------------------
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


# =============================================================================
# Solver 3: Cell-Method
# =============================================================================

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

        # Pass 1: coolant direction -- r_last -> 0
        for r in range(r_last, -1, -1):
            for t in range(geo.n_tubes):

                s_range = range(0, n_segments, 1) if r % 2 == 0 else range(s_exit, -1, -1)

                for s in s_range:
                    T_c_in = get_coolant_inlet(r, t, s, n_segments, T_c, ops, geo)
                    T_a_in = get_staggered_air_inlet(r, t, s, T_a_old, ops, geo)

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

                    raw_T_c = T_c_in - P1_loc * (T_c_in - T_a_in)
                    max_raw_cold = max(max_raw_cold, abs(raw_T_c - T_c_old[r, t, s]))
                    T_c[r, t, s] = (1 - omega) * T_c_old[r, t, s] + omega * raw_T_c

        # Pass 2: air direction -- 0 -> r_last
        for r in range(geo.n_rows):
            for t in range(geo.n_tubes):

                s_range = range(0, n_segments) if r % 2 == 0 else range(s_exit, -1, -1)
                for s in s_range:
                    T_c_in = get_coolant_inlet(r, t, s, n_segments, T_c, ops, geo)
                    T_a_in = get_staggered_air_inlet(r, t, s, T_a, ops, geo)

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


def solve_it_cell(n_segments: int = 10, omega: float = 0.2, ops=None, geo=None):
    """omega: same under-relaxation factor as the other solvers, applied per-cell."""
    # --- 1. Load static configuration -----------------------------------
    if geo is None:
        geo = get_geometry()
    if ops is None:
        ops = get_operating_conditions()

    dm_coolant = ops.m_coolant / geo.n_tubes
    dm_air = ops.m_o / (geo.n_tubes * n_segments)

    # --- 2. Grid initialization (row, tube, segment) ----------------------
    T_c = np.full((geo.n_rows, geo.n_tubes, n_segments), ops.T_coolant_in)
    T_a = np.full((geo.n_rows, geo.n_tubes, n_segments), ops.T_air_in)

    dT_hot_it = dT_hot_it_init
    dT_cold_it = dT_cold_it_init

    # AI-REVIEW: two-stage relaxation (omega_warm + raw-residual convergence
    # check in _relax_cell_grid). Empirically fast/stable at the settings
    # used throughout this app, but the same physical problem has been
    # observed converging in ~8 iterations at one n_segments/omega
    # combination and ~1000 at another -- a documented robustness gap, not
    # a known correctness bug. See CLAUDE.md.
    # Stage 1: coarse, fast propagation at a large (capped) omega -- gets
    # the grid close without chasing tight precision at an undamped step.
    omega_warm = min(1.0, 10 * omega)
    (T_c, T_a, dT_hot_it, dT_cold_it, T_coolant_out, T_air_out,
     hist_hot_1, hist_cold_1, hist_Tc_1, hist_Ta_1) = _relax_cell_grid(
        T_c, T_a, omega_warm, dT_hot_it, dT_cold_it,
        geo, ops, n_segments, dm_coolant, dm_air,
        threshold=1e-1, max_iter=50, min_iter=3
    )

    # Stage 2: fine polish at the requested omega, warm-started from stage 1
    ## in the initial setup the temps dont't change anymore within 0.01 K, but left in here for completeness, as it doesnt take long
    (T_c, T_a, dT_hot_it, dT_cold_it, T_coolant_out, T_air_out,
     hist_hot_2, hist_cold_2, hist_Tc_2, hist_Ta_2) = _relax_cell_grid(
        T_c, T_a, omega, dT_hot_it, dT_cold_it,
        geo, ops, n_segments, dm_coolant, dm_air,
        threshold=1e-3, max_iter=1000, min_iter=3
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
    dQ = ops.m_coolant * coolant_state.cp * (ops.T_coolant_in - T_coolant_out)

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


# =============================================================================
# Scenario orchestrator -- runs all three, no printing/plotting
# =============================================================================

def solve_scenario(ops, geo=None, omega: float = 0.2, n_segments: int = 20) -> ScenarioResult:
    """Runs LMTD/NTU/Cell against the same ops/geo/omega. Pure computation:
    no printing, no plotting -- see run_scenario.plot_scenario() for that."""
    if geo is None:
        geo = get_geometry()

    lmtd = SolverResult(*solve_it_LMTD(omega=omega, ops=ops, geo=geo))
    ntu = SolverResult(*solve_it_NTU(omega=omega, ops=ops, geo=geo))
    cell = SolverResult(*solve_it_cell(n_segments=n_segments, omega=omega, ops=ops, geo=geo))

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
