"""
solvers.py
==========
Two independent iterative solvers for the same dry-cooler heat-transfer
problem -- LMTD-based and NTU/P-based -- both driving toward the same
converged (T_coolant_out, T_air_out, Q). Kept side by side so they can
be cross-checked against each other (see main.py's convergence plot).

Layout: small math helpers first (calc_LMTD, calc_P), then the two
solvers. Both solvers share the same overall structure -- load static
config, set up the iteration, loop with under-relaxation until the
outlet temperatures stop moving, then run final diagnostics -- so if
you're comparing them, read solve_it_LMTD first.
"""

from math import log, pi, exp
import numpy as np

from src.parameters import OperatingConditions
from src.dry_cooler_physics import get_geometry
from src.fluid_properties import get_fluid_properties
from src.heat_transfer_core import calc_overall_k, calc_diagnostics

# Initial guesses for the iteration (both solvers start from the same point)
dT_hot_it_init = 30
dT_cold_it_init = 30


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


# -----------------------------------------------------------------------
# Source: Brunner, p. 14 (citing Holman, "Heat Transfer", 10th ed., 2010)
# -----------------------------------------------------------------------
def calc_P(NTU1, R1) -> float:
    """P1 correlation: single-pass effectiveness as a function of NTU1
    and the capacity-rate ratio R1."""
    return 1 - exp(((NTU1**0.22) / R1) * (exp(-R1 * NTU1**0.78) - 1))

# ALTERNATE (inactive) -- not currently used, kept for reference.
# Source: VDI Heat Atlas, Section [not recorded], p. 50 [eq. # not recorded]
#     return 1 - exp((exp(-R1 * NTU1) - 1) / R1)


# =============================================================================
# Solver 1: LMTD-based
# =============================================================================

def solve_it_LMTD(omega: float = 0.1):
    """
    omega: centralized under-relaxation factor. Pass the SAME value into
    solve_it_NTU() when comparing the two solvers so the step size is
    consistent between them.
    """
    # --- 1. Load static configuration -----------------------------------
    geo = get_geometry()
    ops = OperatingConditions()

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

    history_hot = []
    history_cold = []
    history_T_coolant = []                     # T_coolant_out at each iteration
    history_T_air = []                         # T_air_out at each iteration

    # Fluid states from the FINAL iteration, used afterwards for the
    # final-iteration diagnostics (Pr/Re/Nu).
    coolant_state = None
    air_state = None

    # while either outlet temperature is still moving significantly
    while (abs(diff_hot) > threshold) or (diff_cold > threshold):

        dT_hot = dT_hot_it
        dT_cold = dT_cold_it

        # Dynamic mean temperatures for this iteration
        T_air_mean = (T_air_out + ops.T_air_in) / 2.0
        T_coolant_mean = (ops.T_coolant_in + T_coolant_out) / 2.0

        # Fluid states at the mean temperature
        coolant_state = get_fluid_properties(ops.coolant_type, T_coolant_mean, ops.P_coolant)
        air_state = get_fluid_properties(ops.air_type, T_air_mean, ops.P_air)

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

        # Under-relaxation: slow down the update (centralized omega)
        T_coolant_out = (1 - omega) * T_coolant_out + omega * raw_T_coolant_out
        T_air_out = (1 - omega) * T_air_out + omega * raw_T_air_out

        # New diffs for the loop condition
        new_dT_hot = ops.T_coolant_in - T_air_out
        new_dT_cold = T_coolant_out - ops.T_air_in

        diff_hot = new_dT_hot - dT_hot_it
        diff_cold = new_dT_cold - dT_cold_it

        dT_hot_it = new_dT_hot
        dT_cold_it = new_dT_cold

        # Live progress printout -- useful for spotting slow/diverging runs
        print(f"Hot Error: {diff_hot:.5f} | Cold Error: {diff_cold:.5f} | T_coolant_out: {T_coolant_out:.2f} | T_air_out: {T_air_out:.2f}")

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

def solve_it_NTU(omega: float = 0.1):
    """
    omega: centralized under-relaxation factor. Pass the SAME value into
    solve_it_LMTD() when comparing the two solvers so the step size is
    consistent between them. (Originally this solver used omega=0.5 --
    override at the call site if you want to reproduce that behavior.)
    """
    # --- 1. Load static configuration -----------------------------------
    geo = get_geometry()
    ops = OperatingConditions()

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

    history_hot = []
    history_cold = []
    history_T_coolant = []                     # T_coolant_out at each iteration
    history_T_air = []                         # T_air_out at each iteration

    coolant_state = None
    air_state = None

    # while either outlet temperature is still moving significantly
    while (abs(diff_hot) > threshold) or (diff_cold > threshold):
        n = geo.n_rows

        dT_hot = dT_hot_it
        dT_cold = dT_cold_it

        # Dynamic mean temperatures for this iteration
        T_air_mean = (T_air_out + ops.T_air_in) / 2.0
        T_coolant_mean = (ops.T_coolant_in + T_coolant_out) / 2.0

        # Fluid states at the mean temperature
        coolant_state = get_fluid_properties(ops.coolant_type, T_coolant_mean, ops.P_coolant)
        air_state = get_fluid_properties(ops.air_type, T_air_mean, ops.P_air)

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

        # Under-relaxation: slow down the update (centralized omega)
        T_coolant_out = (1 - omega) * T_coolant_out + omega * raw_T_coolant_out
        T_air_out = (1 - omega) * T_air_out + omega * raw_T_air_out

        # New diffs for the loop condition
        new_dT_hot = ops.T_coolant_in - T_air_out
        new_dT_cold = T_coolant_out - ops.T_air_in

        diff_hot = new_dT_hot - dT_hot_it
        diff_cold = new_dT_cold - dT_cold_it

        dT_hot_it = new_dT_hot
        dT_cold_it = new_dT_cold

        # Live progress printout -- useful for spotting slow/diverging runs
        print(f"Hot Error: {diff_hot:.5f} | Cold Error: {diff_cold:.5f} | T_coolant_out: {T_coolant_out:.2f} | T_air_out: {T_air_out:.2f}")

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

def solve_it_cell(n_segments: int = 10, omega: float = 0.1):
    """
    omega: SAME centralized under-relaxation factor as the other two
    solvers. Applied per-cell (blended against that cell's previous-outer-
    iteration value) so the step size stays comparable across all three
    solvers when plotted together.
    """
    omega = 1 # TEMP BECAUSE SLOW

    # --- 1. Load static configuration -----------------------------------
    geo = get_geometry()
    ops = OperatingConditions()

    # Per-cell static mass flows (mirrors dm_coolant/dm_air in the other two
    # solvers, but divided down to a single cell's share). Coolant splits
    # only across the n_tubes parallel tubes (same flow serially through all
    # n_segments of a tube). Air splits across the whole frontal face of a
    # row, i.e. across BOTH n_tubes (transverse) AND n_segments (along the
    # tube length).
    dm_coolant = ops.m_coolant / geo.n_tubes
    dm_air = ops.m_o / (geo.n_tubes * n_segments)

    # --- 2. Grid initialization (row, tube, segment) ----------------------
    # T_coolant / T_air: storing the outlet of each segment
    T_c = np.full((geo.n_rows, geo.n_tubes, n_segments), ops.T_coolant_in)
    T_a = np.full((geo.n_rows, geo.n_tubes, n_segments), ops.T_air_in)

    threshold = 1e-3                           # grid convergence threshold [K]
    max_iter = 1000

    # Aggregate (scalar) outlet history -- tracked in parallel with the
    # grid so this solver's convergence can be plotted on the same axes
    # as LMTD/NTU (which only ever track scalar outlet temps).
    dT_hot_it = dT_hot_it_init
    dT_cold_it = dT_cold_it_init
    history_hot = []
    history_cold = []
    history_T_coolant = []
    history_T_air = []

    # Coolant exits the array at the last row, at whichever end of the
    # serpentine that row's traversal finishes on (see get_coolant_inlet).
    r_last = geo.n_rows - 1
    r_exit = 0
    s_exit = n_segments - 1

    # while the grid hasn't settled (checked at the end of each full sweep)
    for iteration in range(max_iter):
        T_c_old = T_c.copy()
        T_a_old = T_a.copy()

        # Pass 1: coolant direction -- r_last -> 0
        for r in range(r_last, -1, -1):
            for t in range(geo.n_tubes):

                s_range = range(0, n_segments, 1) if r % 2 == 0 else range(s_exit, -1, -1)
            
                for s in s_range:

                    # Local inlet temperatures
                    # Coolant: logic depends on circuit layout -- for now simple serpentines in z-direction
                    T_c_in = get_coolant_inlet(r, t, s, n_segments, T_c, ops, geo)

                    # Air: Staggered mixing logic (Source: VDI C1, 3.1)
                    T_a_in = get_staggered_air_inlet(r, t, s, T_a_old, ops, geo) # air not yet updated

                    # Local properties (k) at local T_mean
                    T_mean = (T_c_in + T_a_in) / 2.0  ## TODO: this seems very different, check later
                    props_c = get_fluid_properties(ops.coolant_type, T_mean, ops.P_coolant)
                    props_a = get_fluid_properties(ops.air_type, T_mean, ops.P_air)

                    k_local = calc_overall_k(geo, ops, props_c, props_a, T_a_in)

                    # Local dimensionless groups
                    A_cell = geo.A / n_segments  # area per cell
                    W1 = dm_coolant * props_c.cp
                    W2 = dm_air * props_a.cp
                    R_loc = W1 / W2
                    NTU_loc = (k_local * A_cell) / W1

                    P1_loc = calc_P(NTU_loc, R_loc)
                    # P2_loc = P1_loc * R_loc # P2 => T_a not required in coolant pass

                    # Update temperatures (under-relaxed against this cell's
                    # previous-outer-iteration value, same omega as LMTD/NTU)
                    raw_T_c = T_c_in - P1_loc * (T_c_in - T_a_in)
                    # raw_T_a = T_a_in + P2_loc * (T_c_in - T_a_in)
                    T_c[r, t, s] = (1 - omega) * T_c_old[r, t, s] + omega * raw_T_c
                    # T_a[r, t, s] = (1 - omega) * T_a_old[r, t, s] + omega * raw_T_a

        # Pass 2: air direction -- 0 -> r_last
        for r in range(geo.n_rows):
            for t in range(geo.n_tubes):

                s_range = range(0, n_segments) if r % 2 == 0 else range (s_exit, -1, -1)
                for s in s_range:
                    T_c_in = get_coolant_inlet(r, t, s, n_segments, T_c, ops, geo)  # already fresh from pass 1
                    T_a_in = get_staggered_air_inlet(r, t, s, T_a, ops, geo)        # fresh within this pass

                    T_mean = (T_c_in + T_a_in) / 2.0
                    props_c = get_fluid_properties(ops.coolant_type, T_mean, ops.P_coolant)
                    props_a = get_fluid_properties(ops.air_type, T_mean, ops.P_air)
                    k_local = calc_overall_k(geo, ops, props_c, props_a, T_a_in)

                    A_cell = geo.A / n_segments
                    W1 = dm_coolant * props_c.cp
                    W2 = dm_air * props_a.cp
                    R_loc = W1 / W2
                    NTU_loc = (k_local * A_cell) / W1
                    P1_loc = calc_P(NTU_loc, R_loc)
                    P2_loc = P1_loc * R_loc

                    raw_T_a = T_a_in + P2_loc * (T_c_in - T_a_in)
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

        # Grid-based errors: max absolute change over the WHOLE array this
        # iteration (computed once, reused for both the printout and the
        # convergence check below).
        max_err_hot = np.max(np.abs(T_a - T_a_old))    # air grid -- "hot" side error
        max_err_cold = np.max(np.abs(T_c - T_c_old))   # coolant grid -- "cold" side error

        # Live progress printout -- same format as LMTD/NTU, but the error
        # terms are the max change over the whole grid (not a single scalar
        # diff), and the outlet temps are grid averages (as they will be at
        # convergence).
        print(f"Hot Error: {max_err_hot:.5f} | Cold Error: {max_err_cold:.5f} | "
              f"T_coolant_out: {T_coolant_out:.2f} | T_air_out: {T_air_out:.2f}")

        # Convergence check (grid-based -- the authoritative check, since it
        # only passes once every cell, not just the aggregate outlets, has
        # settled)
        if max_err_cold < threshold and max_err_hot < threshold:
            break

    # --- 3. Final-iteration diagnostics -----------------------------------
    # The cell method has no single "final iteration state" the way LMTD/NTU
    # do (k/Pr/Re/Nu vary cell to cell) -- so, for a diagnostics summary
    # comparable to the other two solvers, evaluate one representative
    # state at the overall (array-average) mean bulk temperatures.
    T_coolant_mean = (ops.T_coolant_in + T_coolant_out) / 2.0
    T_air_mean = (ops.T_air_in + T_air_out) / 2.0
    coolant_state = get_fluid_properties(ops.coolant_type, T_coolant_mean, ops.P_coolant)
    air_state = get_fluid_properties(ops.air_type, T_air_mean, ops.P_air)

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
    # Air: Staggered mixing logic (Source: VDI C1, 3

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