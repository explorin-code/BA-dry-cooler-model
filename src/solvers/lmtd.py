"""
lmtd.py -- LMTD: the whole exchanger as pure counterflow (no correction factor F, deliberately the simplest of the three
methods -- see docs/model-notes.md).
"""

from math import log

from src.solvers.common import fill_defaults
from src.solvers.relaxation import relax_lmtd_ntu


def calc_lmtd(Delta_T_hot: float, Delta_T_cold: float) -> float:
    """Log-mean temperature difference, guarded against crossover and
    Delta_T_hot == Delta_T_cold.
    Source: Holman, Heat Transfer, 10th ed. (2010), Sec. 10-5."""
    if Delta_T_hot <= 0 or Delta_T_cold <= 0:
        return 1e-5                                # crossover -- not physical, clamp
    if abs(Delta_T_hot - Delta_T_cold) < 1e-5:
        return Delta_T_hot                     # avoid 0/0
    return (Delta_T_hot - Delta_T_cold) / log(Delta_T_hot / Delta_T_cold)


def _lmtd_step(k, coolant_state, air_state, geo, ops, Delta_T_hot, Delta_T_cold):
    """Q = k A dT_lm, then both outlets from the energy balances.
    Sources: Holman (2010), Sec. 10-5; Kröger, Air-Cooled Heat Exchangers
    and Cooling Towers, Vol. 2, Ch. 8."""
    Q_dot = geo.A_tot * k * calc_lmtd(Delta_T_hot, Delta_T_cold)
    raw_T_c_o = ops.T_c_i - Q_dot / (ops.m_dot_c * coolant_state.cp)
    raw_T_a_o = ops.T_a_i + Q_dot / (ops.m_dot_a * air_state.cp)
    return raw_T_c_o, raw_T_a_o, Q_dot


def solve_lmtd(omega: float = None, ops=None, geo=None, settings=None):
    """omega: under-relaxation factor (default settings.central_omega)."""
    ops, geo, settings = fill_defaults(ops, geo, settings)
    if omega is None:
        omega = settings.central_omega
    return relax_lmtd_ntu(_lmtd_step, "[LMTD]", omega, settings.convergence_threshold, ops, geo, settings, pass_Delta_T=True)
