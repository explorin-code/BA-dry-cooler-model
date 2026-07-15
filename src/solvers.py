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