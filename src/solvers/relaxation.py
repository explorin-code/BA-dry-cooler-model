"""
relaxation.py -- the outer loop shared by LMTD and NTU: mean-temperature states -> k -> one step -> under-relaxed
outlet temperatures, with the oscillation damping against the property cache's rounding steps.
"""

import logging

from src.heat_transfer_core import calc_diagnostics, calc_overall_k
from src.solvers.common import mean_states
from src.solvers.results import SolverResult

logger = logging.getLogger(__name__)

# Oscillation damping of the shared LMTD/NTU outer loop (relax_lmtd_ntu)
OSC_STALL_RATIO = 0.7      # an oscillation shrinking slower than this per two steps counts as stalled
OMEGA_MIN = 0.01           # lower bound when halving omega on stalled oscillations


def relax_lmtd_ntu(step_fn, tag, omega, threshold, ops, geo, settings, pass_Delta_T=False):
    """Outer loop shared by LMTD and NTU: mean-temperature states -> k ->
    step_fn -> under-relaxed outlet temperatures, until both outlets settle.
    step_fn(k, coolant_state, air_state, geo, ops[, Delta_T_hot,
    Delta_T_cold]) returns the raw (T_c_o, T_a_o, Q_dot); the
    temperature differences are only passed if pass_Delta_T (LMTD)."""
    Delta_T_hot = settings.Delta_T_hot_init
    Delta_T_cold = settings.Delta_T_cold_init
    T_c_o = ops.T_c_i - Delta_T_hot
    T_a_o = ops.T_a_i + Delta_T_cold

    # Start with a larger omega (never below the requested one), drop to the
    # requested one for good at the first growing step or second sign flip.
    # Heuristic, tuned empirically.
    omega_active = max(omega, min(0.5, 10 * omega))
    sign_flips_hot = sign_flips_cold = 0
    stalls = 0                                     # consecutive non-shrinking oscillation steps
    diff_hot = diff_cold = 1.0                     # > threshold: run at least once
    history_hot, history_cold, history_T_c_o, history_T_a_o = [], [], [], []

    while abs(diff_hot) > threshold or abs(diff_cold) > threshold:
        if len(history_hot) >= settings.outer_max_iter:
            raise RuntimeError(f"{tag.strip()} did not converge in {settings.outer_max_iter} iterations "
                               f"(last changes {diff_hot:.2e} / {diff_cold:.2e} K) -- omega {omega} too large?")
        coolant_state, air_state = mean_states(ops, T_c_o, T_a_o)
        k = calc_overall_k(geo=geo, ops=ops, coolant_state=coolant_state, air_state=air_state,
                           T_a_m=(ops.T_a_i + T_a_o) / 2.0)

        step_args = (k, coolant_state, air_state, geo, ops)
        if pass_Delta_T:
            step_args += (Delta_T_hot, Delta_T_cold)
        raw_T_c_o, raw_T_a_o, Q_dot = step_fn(*step_args)
        T_c_o = (1 - omega_active) * T_c_o + omega_active * raw_T_c_o
        T_a_o = (1 - omega_active) * T_a_o + omega_active * raw_T_a_o

        new_Delta_T_hot = ops.T_c_i - T_a_o
        new_Delta_T_cold = T_c_o - ops.T_a_i
        diff_hot = new_Delta_T_hot - Delta_T_hot
        diff_cold = new_Delta_T_cold - Delta_T_cold
        Delta_T_hot, Delta_T_cold = new_Delta_T_hot, new_Delta_T_cold

        logger.debug("%-7sHot Error: %.5f | Cold Error: %.5f | T_coolant_out: %.2f | T_air_out: %.2f | omega: %s",
                     tag, diff_hot, diff_cold, T_c_o, T_a_o, omega_active)

        if history_hot:
            sign_flips_hot += diff_hot * history_hot[-1] < 0
            sign_flips_cold += diff_cold * history_cold[-1] < 0
            if omega_active > omega and (abs(diff_hot) >= abs(history_hot[-1])
                                         or abs(diff_cold) >= abs(history_cold[-1])
                                         or sign_flips_hot >= 2 or sign_flips_cold >= 2):
                omega_active = omega
            # Stalled oscillation (e.g. the property cache's rounding steps,
            # 0.01 K, leave a period-2 limit cycle): a sign flip whose amplitude
            # did not shrink below OSC_STALL_RATIO of the one two steps back.
            # Twice in a row -> halve omega (down to OMEGA_MIN). A converging
            # oscillation shrinks and is left alone, so normal runs keep their omega.
            elif len(history_hot) >= 2 and any(
                    d * prev < 0 and abs(d) > OSC_STALL_RATIO * abs(prev2)
                    for d, prev, prev2 in ((diff_hot, history_hot[-1], history_hot[-2]),
                                           (diff_cold, history_cold[-1], history_cold[-2]))):
                stalls += 1
                if stalls >= 2:
                    omega_active = max(omega_active / 2, OMEGA_MIN)
                    stalls = 0
            else:
                stalls = 0

        history_hot.append(diff_hot)
        history_cold.append(diff_cold)
        history_T_c_o.append(T_c_o)
        history_T_a_o.append(T_a_o)

    # Diagnostics at the states the last k was evaluated with.
    diagnostics = calc_diagnostics(geo=geo, ops=ops, coolant_state=coolant_state, air_state=air_state,
                                   T_a_m=(ops.T_a_i + T_a_o) / 2.0)
    return SolverResult(k, Q_dot, T_c_o, T_a_o,
                        history_hot, history_cold, history_T_c_o, history_T_a_o, diagnostics)
