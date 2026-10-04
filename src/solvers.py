"""
solvers.py
==========
Three solvers for the same dry-cooler problem:

  1. Results         SolverResult, ScenarioResult
  2. Shared helpers  outer relaxation loop (LMTD + NTU), coolant path and
                     element effectiveness (NTU + Cell), defaults, states
  3. LMTD            log-mean temperature difference, whole exchanger
  4. NTU             element-wise e-NTU field (Cabezas-Gomez et al. 2007)
  5. Cell            VDI cell method with local properties and k
  6. Scenario        solve_scenario() runs all three

Pure computation: progress is logged at DEBUG level (insight mode),
plotting lives in plot_profiles.py / plot_solver.py.
"""

from dataclasses import dataclass
from functools import partial
import logging
from math import log, exp, isclose
import time
import warnings
import numpy as np

from src.operating_conditions import get_operating_conditions
from src.dry_cooler_physics import get_geometry
from src.fluid_properties import get_fluid_properties, get_air_properties
from src.heat_transfer_core import calc_overall_k, calc_diagnostics
from src.solver_settings import get_solver_settings

logger = logging.getLogger(__name__)


# =============================================================================
# 1. Results
# =============================================================================

@dataclass
class SolverResult:
    k: float
    Q_dot: float
    theta_c_o: float
    theta_a_o: float
    history_hot: list                # per-iteration change of Delta_theta_hot [K]
    history_cold: list               # per-iteration change of Delta_theta_cold [K]
    history_theta_c_o: list
    history_theta_a_o: list
    diagnostics: dict                # Pr/Re/Nu/alpha at the overall mean temperatures
    # (row, tube, segment) grids -- Cell and NTU only, LMTD leaves them None
    theta_c_grid: object = None
    theta_a_grid: object = None
    k_grid: object = None
    dQdL_grid: object = None
    local_diagnostics: dict = None   # Cell only: Re/Nu grids, not yet plotted
    solve_time: float = None         # wall time [s], set by solve_scenario


@dataclass
class ScenarioResult:
    lmtd: SolverResult
    ntu: SolverResult
    cell: SolverResult


# =============================================================================
# 2. Shared helpers
# =============================================================================

def _defaults(ops, geo, settings):
    """Fills in whatever the caller didn't pass."""
    if geo is None:
        geo = get_geometry()
    if ops is None:
        ops = get_operating_conditions(geo=geo)
    if settings is None:
        settings = get_solver_settings()
    return ops, geo, settings


def _mean_states(ops, theta_c_o, theta_a_o):
    """Coolant and air states at their mean temperatures (inlet + outlet)/2."""
    coolant_state = get_fluid_properties(ops, (ops.theta_c_i + theta_c_o) / 2.0, ops.p_c)
    air_state = get_air_properties(ops, (ops.theta_a_i + theta_a_o) / 2.0, ops.p_a)
    return coolant_state, air_state


def cell_path_order(N_r, n_segments):
    """(row, segment) in coolant flow order: enters the last row, serpentine
    towards row 0, axial direction alternating per row. Used by NTU, Cell
    and analysis.py, so all three walk the same path."""
    path = []
    for r in range(N_r - 1, -1, -1):
        s_range = range(0, n_segments) if r % 2 == 0 else range(n_segments - 1, -1, -1)
        for s in s_range:
            path.append((r, s))
    return path


def calc_element_effectiveness(UA_e: float, C_c_e: float) -> float:
    """Element effectiveness Gamma = 1 - exp(-UA_e / C_c_e) of a small
    mixed(tube)-unmixed(air) cross-flow element, valid for C_c_e << C_h_e.
    Source: Cabezas-Gomez, Navarro & Saiz-Jabardo (2007), J. Heat Transfer
    129, p. 283, Eq. (1)."""
    return 1.0 - exp(-UA_e / C_c_e)


def _relax_lmtd_ntu(step_fn, tag, omega, ops, geo, settings, pass_Delta_theta=False):
    """Outer loop shared by LMTD and NTU: mean-temperature states -> k ->
    step_fn -> under-relaxed outlet temperatures, until both outlets settle.
    step_fn(k, coolant_state, air_state, geo, ops[, Delta_theta_hot,
    Delta_theta_cold]) returns the raw (theta_c_o, theta_a_o, Q_dot); the
    temperature differences are only passed if pass_Delta_theta (LMTD)."""
    Delta_theta_hot = settings.Delta_theta_hot_init
    Delta_theta_cold = settings.Delta_theta_cold_init
    theta_c_o = ops.theta_c_i - Delta_theta_hot
    theta_a_o = ops.theta_a_i + Delta_theta_cold

    # Start with a larger omega, drop to the requested one for good at the
    # first growing step or second sign flip. Heuristic, tuned empirically.
    omega_active = min(0.5, 10 * omega)
    sign_flips_hot = sign_flips_cold = 0
    diff_hot = diff_cold = 1.0                     # > threshold: run at least once
    history_hot, history_cold, history_theta_c_o, history_theta_a_o = [], [], [], []

    while abs(diff_hot) > settings.convergence_threshold or abs(diff_cold) > settings.convergence_threshold:
        coolant_state, air_state = _mean_states(ops, theta_c_o, theta_a_o)
        k = calc_overall_k(geo=geo, ops=ops, coolant_state=coolant_state, air_state=air_state,
                           theta_a_o=theta_a_o)

        step_args = (k, coolant_state, air_state, geo, ops)
        if pass_Delta_theta:
            step_args += (Delta_theta_hot, Delta_theta_cold)
        raw_theta_c_o, raw_theta_a_o, Q_dot = step_fn(*step_args)
        theta_c_o = (1 - omega_active) * theta_c_o + omega_active * raw_theta_c_o
        theta_a_o = (1 - omega_active) * theta_a_o + omega_active * raw_theta_a_o

        new_Delta_theta_hot = ops.theta_c_i - theta_a_o
        new_Delta_theta_cold = theta_c_o - ops.theta_a_i
        diff_hot = new_Delta_theta_hot - Delta_theta_hot
        diff_cold = new_Delta_theta_cold - Delta_theta_cold
        Delta_theta_hot, Delta_theta_cold = new_Delta_theta_hot, new_Delta_theta_cold

        logger.debug("%-7sHot Error: %.5f | Cold Error: %.5f | T_coolant_out: %.2f | T_air_out: %.2f | omega: %s",
                     tag, diff_hot, diff_cold, theta_c_o, theta_a_o, omega_active)

        if history_hot:
            sign_flips_hot += diff_hot * history_hot[-1] < 0
            sign_flips_cold += diff_cold * history_cold[-1] < 0
            if (abs(diff_hot) >= abs(history_hot[-1]) or abs(diff_cold) >= abs(history_cold[-1])
                    or sign_flips_hot >= 2 or sign_flips_cold >= 2):
                omega_active = omega

        history_hot.append(diff_hot)
        history_cold.append(diff_cold)
        history_theta_c_o.append(theta_c_o)
        history_theta_a_o.append(theta_a_o)

    # Diagnostics at the states the last k was evaluated with.
    diagnostics = calc_diagnostics(geo=geo, ops=ops, coolant_state=coolant_state, air_state=air_state,
                                   theta_a_o=theta_a_o)
    return SolverResult(k, Q_dot, theta_c_o, theta_a_o,
                        history_hot, history_cold, history_theta_c_o, history_theta_a_o, diagnostics)


# =============================================================================
# 3. LMTD
# =============================================================================
# Whole exchanger as pure counterflow (no correction factor F, deliberately
# the simplest of the three methods -- see CLAUDE.md).

def calc_LMTD(Delta_theta_hot: float, Delta_theta_cold: float) -> float:
    """Log-mean temperature difference, guarded against crossover and
    Delta_theta_hot == Delta_theta_cold.
    Source: Holman, Heat Transfer, 10th ed. (2010), Sec. 10-5."""
    if Delta_theta_hot <= 0 or Delta_theta_cold <= 0:
        return 1e-5                                # crossover -- not physical, clamp
    if abs(Delta_theta_hot - Delta_theta_cold) < 1e-5:
        return Delta_theta_hot                     # avoid 0/0
    return (Delta_theta_hot - Delta_theta_cold) / log(Delta_theta_hot / Delta_theta_cold)


def _lmtd_step(k, coolant_state, air_state, geo, ops, Delta_theta_hot, Delta_theta_cold):
    """Q = k A dT_lm, then both outlets from the energy balances.
    Sources: Holman (2010), Sec. 10-5; Kröger, Air-Cooled Heat Exchangers
    and Cooling Towers, Vol. 2, Ch. 8."""
    Q_dot = geo.A_tot * k * calc_LMTD(Delta_theta_hot, Delta_theta_cold)
    raw_theta_c_o = ops.theta_c_i - Q_dot / (ops.m_dot_c * coolant_state.cp)
    raw_theta_a_o = ops.theta_a_i + Q_dot / (ops.m_dot_a * air_state.cp)
    return raw_theta_c_o, raw_theta_a_o, Q_dot


def solve_it_LMTD(omega: float = None, ops=None, geo=None, settings=None):
    """omega: under-relaxation factor (default settings.central_omega)."""
    ops, geo, settings = _defaults(ops, geo, settings)
    if omega is None:
        omega = settings.central_omega
    return _relax_lmtd_ntu(_lmtd_step, "[LMTD]", omega, ops, geo, settings, pass_Delta_theta=True)


# =============================================================================
# 4. NTU -- element-wise e-NTU method for multipass counter-cross-flow
# =============================================================================
# Source: Cabezas-Gomez, Navarro & Saiz-Jabardo (2007), "Thermal Performance
# of Multipass Parallel and Counter-Cross-Flow Heat Exchangers", J. Heat
# Transfer 129, 282-289, DOI 10.1115/1.2430719. Arrangement G^c_(N_r,1).
# Paper notation inside the element functions: h = coolant (tube fluid,
# mixed), c = air (external fluid, unmixed).
# Constant U and cp within one field solve (as in the paper); the outer loop
# updates them between field solves -- revisit if results look off.

def calc_element_outlets(T_h_in: float, T_c_in: float, Gamma: float, B: float):
    """Element outlet temperatures (T_c_out, T_h_out), B = C_c Gamma / C_h.
    Source: Cabezas-Gomez et al. (2007), p. 283, Eqs. (3)-(4)."""
    T_c_out = ((B + 2.0 * (1.0 - Gamma)) / (2.0 + B)) * T_c_in + (2.0 * Gamma / (2.0 + B)) * T_h_in
    T_h_out = ((2.0 - B) / (2.0 + B)) * T_h_in + (2.0 * B / (2.0 + B)) * T_c_in   # no Gamma here
    return T_c_out, T_h_out


def solve_ntu_field(N_r: int, n_elements: int, C_h_circuit: float, C_c_element: float,
                    UA_element: float, theta_c_i: float, theta_a_i: float,
                    tol: float, max_iter: int):
    """Converged element field of ONE coolant circuit (all circuits are
    identical). Plain numbers only, so validation_ntu.py runs this exact code.
    Returns (theta_c, theta_a, Q_element, theta_c_o, theta_a_o, residuals);
    grids are [row, element] at the physical axial position.
    Assumption (not from the paper): air leaving element (r-1, e) enters
    (r, e), no averaging -- reproduces the paper's closed forms."""
    Gamma = calc_element_effectiveness(UA_element, C_c_element)
    B = C_c_element * Gamma / C_h_circuit                      # p. 283, below Eq. (4)
    if C_c_element / C_h_circuit > 0.1:                        # Eq. (1) needs C_c << C_h
        warnings.warn(f"NTU element field: C_c^e/C_h^e = {C_c_element / C_h_circuit:.3f} is not "
                      f"<< 1 -- increase n_elements (currently {n_elements}).")

    path = cell_path_order(N_r, n_elements)
    theta_c = np.full((N_r, n_elements), float(theta_c_i))
    theta_a = np.full((N_r, n_elements), float(theta_a_i))     # assumed initial air field
    Q_element = np.zeros((N_r, n_elements))

    residuals = []
    for _ in range(max_iter):
        theta_c_old, theta_a_old = theta_c.copy(), theta_a.copy()

        # March along the circuit (Fig. 2b). Row r-1's air is from the
        # previous sweep -- the lag the iteration converges away.
        T_h = theta_c_i
        for r, e in path:
            T_c_in = theta_a_i if r == 0 else theta_a_old[r - 1, e]
            T_c_out, T_h_out = calc_element_outlets(T_h, T_c_in, Gamma, B)
            theta_a[r, e] = T_c_out
            theta_c[r, e] = T_h_out
            Q_element[r, e] = C_h_circuit * (T_h - T_h_out)
            T_h = T_h_out

        # Own criterion (max-norm); the paper compares mean air outlets (Fig. 2a).
        residual = max(np.max(np.abs(theta_c - theta_c_old)), np.max(np.abs(theta_a - theta_a_old)))
        residuals.append(residual)
        if residual < tol:
            break
    else:
        raise RuntimeError(f"NTU element field did not converge in {max_iter} sweeps "
                           f"(last residual {residuals[-1]:.2e} K, tolerance {tol:.0e} K).")

    theta_c_o = T_h                                            # after the last element
    theta_a_o = float(np.mean(theta_a[N_r - 1, :]))            # mean air outlet (Fig. 2a)
    return theta_c, theta_a, Q_element, theta_c_o, theta_a_o, residuals


def _ntu_step(k, coolant_state, air_state, geo, ops, n_elements, field_tol, field_max_iter, field_out=None):
    """One field solve at the current k and cp; returns outlets and Q.
    Element allocation per Fig. 2(a), items 1.5-1.6. field_out (dict)
    receives the field for the final grids."""
    C_h = ops.m_dot_c * coolant_state.cp                       # paper h = coolant
    C_c = ops.m_dot_a * air_state.cp                           # paper c = air
    A_element = geo.A / n_elements                             # geo.A = one tube, one row
    if not isclose(geo.N_t * geo.N_r * n_elements * A_element, geo.A_tot, rel_tol=1e-9):
        raise ValueError("NTU element areas do not add up to geo.A_tot.")

    theta_c, theta_a, Q_element, theta_c_o, theta_a_o, residuals = solve_ntu_field(
        geo.N_r, n_elements, C_h / geo.N_t, C_c / (geo.N_t * n_elements), k * A_element,
        ops.theta_c_i, ops.theta_a_i, field_tol, field_max_iter)

    # Both sides' Q are equal by construction of Eqs. (3)-(4) -- checked, not assumed.
    Q_h = C_h * (ops.theta_c_i - theta_c_o)
    Q_c = C_c * (theta_a_o - ops.theta_a_i)
    if abs(Q_h - Q_c) > 1e-6 * max(abs(Q_h), 1.0):
        raise RuntimeError(f"NTU energy balance mismatch: Q_h = {Q_h:.3f} W, Q_c = {Q_c:.3f} W.")

    if field_out is not None:
        C_min, C_max = min(C_h, C_c), max(C_h, C_c)
        field_out.update(theta_c=theta_c, theta_a=theta_a, Q_element=Q_element, residuals=residuals,
                         delta_Q=Q_h - Q_c,
                         # Paper's whole-exchanger groups (p. 286 Eqs. 7-8, p. 289)
                         C_star=C_min / C_max, NTU=k * geo.A_tot / C_min,
                         P=(theta_a_o - ops.theta_a_i) / (ops.theta_c_i - ops.theta_a_i), R=C_c / C_h)
    return theta_c_o, theta_a_o, 0.5 * (Q_h + Q_c)


def solve_it_NTU(omega: float = None, ops=None, geo=None, settings=None, n_elements: int = None):
    """omega: under-relaxation of the outer property/k loop (default
    settings.central_omega); n_elements per tube and row (default
    settings.ntu_n_elements)."""
    ops, geo, settings = _defaults(ops, geo, settings)
    if omega is None:
        omega = settings.central_omega
    if n_elements is None:
        n_elements = settings.ntu_n_elements

    last_field = {}
    step_fn = partial(_ntu_step, n_elements=n_elements, field_tol=settings.ntu_field_threshold,
                      field_max_iter=settings.ntu_field_max_iter, field_out=last_field)
    result = _relax_lmtd_ntu(step_fn, "[NTU]", omega, ops, geo, settings)

    # Grids (row, tube, element) like Cell's: the one solved circuit, broadcast to all tubes.
    shape = (geo.N_r, geo.N_t, n_elements)
    widen = lambda grid: np.broadcast_to(grid[:, None, :], shape).copy()
    result.theta_c_grid = widen(last_field['theta_c'])
    result.theta_a_grid = widen(last_field['theta_a'])
    result.k_grid = np.full(shape, result.k)
    result.dQdL_grid = widen(last_field['Q_element'] / (geo.L_p / n_elements))
    return result


# =============================================================================
# 5. Cell -- VDI cell method with local properties and k
# =============================================================================
# Method: VDI-Wärmeatlas (2019), C1 §3.1. Grid [row, tube, segment]; each
# cell is one small mixed(coolant)-unmixed(air) cross-flow element.
# 3D: all N_t tubes. 2D (settings.cell_2d): one representative tube -- the
# grid has a single tube column, so no neighbour averaging of air; the
# result is broadcast to N_t tubes for output. The tube count of a run is
# always taken from the grid shape, not from geo.

def get_coolant_inlet(r, t, s, n_segments, theta_c, ops, geo):
    """Coolant entering cell (r, t, s): serpentine via return bends, entering
    at the last row (hot air exhaust), direction alternating per row.
    Source: Kröger, Air-Cooled Heat Exchangers and Cooling Towers, Vol. 1, Sec. 3.5."""
    forward = r % 2 == 0                                       # even rows run s = 0 -> end
    s_first = 0 if forward else n_segments - 1
    if s == s_first:                                           # row entry: from the U-bend
        return ops.theta_c_i if r == geo.N_r - 1 else theta_c[r + 1, t, s]
    return theta_c[r, t, s - 1 if forward else s + 1]


def get_staggered_air_inlet(r, t, s, theta_a, ops, geo):
    """Air entering cell (r, t, s). Modelling assumption from the staggered
    bank (Kröger, Vol. 1, Sec. 5.1): a tube sits in the gap between two
    tubes of the row upstream, so its air is their mean (edge tube: one)."""
    if r == 0:
        return ops.theta_a_i
    t_neighbour = t - 1 if r % 2 == 0 else t + 1
    if 0 <= t_neighbour < theta_a.shape[1]:
        return (theta_a[r - 1, t, s] + theta_a[r - 1, t_neighbour, s]) / 2.0
    return theta_a[r - 1, t, s]


def _cell_local_states(r, t, s, theta_c, theta_a_for_inlet, ops, geo, n_segments):
    """Inlet temperatures and fluid states of one cell; each fluid at its own
    local mean temperature (inlet + current outlet)/2."""
    theta_c_in = get_coolant_inlet(r, t, s, n_segments, theta_c, ops, geo)
    theta_a_in = get_staggered_air_inlet(r, t, s, theta_a_for_inlet, ops, geo)
    props_c = get_fluid_properties(ops, (theta_c_in + theta_c[r, t, s]) / 2.0, ops.p_c)
    props_a = get_air_properties(ops, (theta_a_in + theta_a_for_inlet[r, t, s]) / 2.0, ops.p_a)
    return theta_c_in, theta_a_in, props_c, props_a


def _solve_cell_step(r, t, s, theta_c, theta_a_for_inlet, ops, geo, n_segments, m_dot_c_tube, m_dot_a_cell):
    """Coolant-side effectiveness P1 = 2B/(2+B) and R = C_c/C_a of one cell
    (Cabezas-Gomez et al. 2007, Eqs. 1, 4); air side P2 = P1 R is Eq. (3).
    theta_a_for_inlet: pass 1 reads last iteration's air, pass 2 the live one."""
    theta_c_in, theta_a_in, props_c, props_a = _cell_local_states(
        r, t, s, theta_c, theta_a_for_inlet, ops, geo, n_segments)
    k_local = calc_overall_k(geo, ops, props_c, props_a, theta_a_in)

    C_dot_c = m_dot_c_tube * props_c.cp
    C_dot_a = m_dot_a_cell * props_a.cp
    Gamma = calc_element_effectiveness(k_local * (geo.A / n_segments), C_dot_a)
    B = C_dot_a * Gamma / C_dot_c
    return theta_c_in, theta_a_in, 2.0 * B / (2.0 + B), C_dot_c / C_dot_a


def _relax_cell_grid(theta_c, theta_a, omega, Delta_theta_hot, Delta_theta_cold,
                     geo, ops, n_segments, m_dot_c_tube, m_dot_a_cell, threshold, max_iter, min_iter=1):
    """Two passes per iteration (coolant along its path, then air row by row)
    until no cell changes by more than threshold."""
    history_hot, history_cold, history_theta_c_o, history_theta_a_o = [], [], [], []
    n_tubes = theta_c.shape[1]                                     # N_t (3D) or 1 (2D)
    cell = lambda grid, r, t, s: _solve_cell_step(r, t, s, theta_c, grid, ops, geo, n_segments,
                                                  m_dot_c_tube, m_dot_a_cell)

    for iteration in range(max_iter):
        theta_c_old, theta_a_old = theta_c.copy(), theta_a.copy()
        max_raw_cold = max_raw_hot = 0.0

        # Pass 1: coolant, along its flow path
        for r, s in cell_path_order(geo.N_r, n_segments):
            for t in range(n_tubes):
                theta_c_in, theta_a_in, P1, R = cell(theta_a_old, r, t, s)
                raw = theta_c_in - P1 * (theta_c_in - theta_a_in)
                max_raw_cold = max(max_raw_cold, abs(raw - theta_c_old[r, t, s]))
                theta_c[r, t, s] = (1 - omega) * theta_c_old[r, t, s] + omega * raw

        # Pass 2: air, row 0 -> last
        for r in range(geo.N_r):
            for t in range(n_tubes):
                for s in (range(n_segments) if r % 2 == 0 else range(n_segments - 1, -1, -1)):
                    theta_c_in, theta_a_in, P1, R = cell(theta_a, r, t, s)
                    raw = theta_a_in + P1 * R * (theta_c_in - theta_a_in)
                    max_raw_hot = max(max_raw_hot, abs(raw - theta_a_old[r, t, s]))
                    theta_a[r, t, s] = (1 - omega) * theta_a_old[r, t, s] + omega * raw

        theta_c_o = float(np.mean(theta_c[0, :, n_segments - 1]))   # coolant leaves row 0
        theta_a_o = float(np.mean(theta_a[geo.N_r - 1, :, :]))       # air leaves the last row

        new_Delta_theta_hot = ops.theta_c_i - theta_a_o
        new_Delta_theta_cold = theta_c_o - ops.theta_a_i
        diff_hot = new_Delta_theta_hot - Delta_theta_hot
        diff_cold = new_Delta_theta_cold - Delta_theta_cold
        Delta_theta_hot, Delta_theta_cold = new_Delta_theta_hot, new_Delta_theta_cold
        history_hot.append(diff_hot)
        history_cold.append(diff_cold)
        history_theta_c_o.append(theta_c_o)
        history_theta_a_o.append(theta_a_o)

        logger.debug("[Cell] Hot Error: %.5f | Cold Error: %.5f | T_coolant_out: %.2f | T_air_out: %.2f | omega: %s",
                     max_raw_hot, max_raw_cold, theta_c_o, theta_a_o, omega)
        if iteration + 1 >= min_iter and max_raw_cold < threshold and max_raw_hot < threshold:
            break

    return (theta_c, theta_a, theta_c_o, theta_a_o,
            history_hot, history_cold, history_theta_c_o, history_theta_a_o)


def _cell_local_grids(theta_c, theta_a, ops, geo, n_segments):
    """k, dQ/dL and Re/Nu per cell on the converged grids (same local states
    as during the iteration). dQ = k dA dT: Baehr, Thermodynamik (2016), Sec. 3.1."""
    shape = theta_c.shape
    k_grid, dQdL_grid = np.zeros(shape), np.zeros(shape)
    local = {key: np.zeros(shape) for key in ('Re_air', 'Nu_air', 'Re_coolant', 'Nu_coolant')}
    A_cell = geo.A / n_segments
    dL_cell = geo.L_p / n_segments

    for r in range(geo.N_r):
        for t in range(shape[1]):
            for s in range(n_segments):
                _, theta_a_in, props_c, props_a = _cell_local_states(r, t, s, theta_c, theta_a, ops, geo, n_segments)
                k_grid[r, t, s] = calc_overall_k(geo, ops, props_c, props_a, theta_a_in)
                dQdL_grid[r, t, s] = k_grid[r, t, s] * A_cell * (theta_c[r, t, s] - theta_a[r, t, s]) / dL_cell
                diag = calc_diagnostics(geo, ops, props_c, props_a, theta_a_in)
                for key in local:
                    local[key][r, t, s] = diag[key]
    return k_grid, dQdL_grid, local


def solve_it_cell(n_segments: int = None, omega: float = None, ops=None, geo=None, settings=None):
    """omega: per-cell relaxation (default settings.cell_omega = 1.0; with
    omega < 1 the per-cell stopping criterion fires before convergence --
    see CLAUDE.md). n_segments per tube and row (default settings.cell_n_segments).
    Single relaxation stage (the former coarse/fine stages were identical at omega = 1)."""
    ops, geo, settings = _defaults(ops, geo, settings)
    if omega is None:
        omega = settings.cell_omega
    if n_segments is None:
        n_segments = settings.cell_n_segments

    m_dot_c_tube = ops.m_dot_c / geo.N_t                       # per tube, also in 2D
    m_dot_a_cell = ops.m_dot_a / (geo.N_t * n_segments)
    n_tubes = 1 if settings.cell_2d else geo.N_t
    theta_c = np.full((geo.N_r, n_tubes, n_segments), ops.theta_c_i)
    theta_a = np.full((geo.N_r, n_tubes, n_segments), ops.theta_a_i)

    (theta_c, theta_a, theta_c_o, theta_a_o,
     history_hot, history_cold, history_theta_c_o, history_theta_a_o) = _relax_cell_grid(
        theta_c, theta_a, omega, settings.Delta_theta_hot_init, settings.Delta_theta_cold_init,
        geo, ops, n_segments, m_dot_c_tube, m_dot_a_cell,
        threshold=settings.cell_threshold, max_iter=settings.cell_max_iter, min_iter=settings.cell_min_iter)

    # Whole-exchanger Q and diagnostics at the overall mean temperatures; the
    # reported k is the mean of the local k Cell actually used (equal cell areas).
    coolant_state, air_state = _mean_states(ops, theta_c_o, theta_a_o)
    Q_dot = ops.m_dot_c * coolant_state.cp * (ops.theta_c_i - theta_c_o)
    diagnostics = calc_diagnostics(geo=geo, ops=ops, coolant_state=coolant_state, air_state=air_state,
                                   theta_a_o=theta_a_o)
    k_grid, dQdL_grid, local_diagnostics = _cell_local_grids(theta_c, theta_a, ops, geo, n_segments)

    if n_tubes == 1:                                           # 2D: broadcast to all tubes for output
        widen = lambda grid: np.broadcast_to(grid, (geo.N_r, geo.N_t, n_segments)).copy()
        theta_c, theta_a, k_grid, dQdL_grid = map(widen, (theta_c, theta_a, k_grid, dQdL_grid))
        local_diagnostics = {key: widen(grid) for key, grid in local_diagnostics.items()}

    return SolverResult(float(np.mean(k_grid)), Q_dot, theta_c_o, theta_a_o,
                        history_hot, history_cold, history_theta_c_o, history_theta_a_o, diagnostics,
                        theta_c_grid=theta_c, theta_a_grid=theta_a, k_grid=k_grid, dQdL_grid=dQdL_grid,
                        local_diagnostics=local_diagnostics)


# =============================================================================
# 6. Scenario
# =============================================================================

def solve_scenario(ops, geo=None, omega: float = None, n_segments: int = None, settings=None,
                   n_elements: int = None) -> ScenarioResult:
    """Runs LMTD, NTU and Cell on one operating point and times each. omega
    applies to LMTD/NTU, Cell uses settings.cell_omega. The solvers share the
    property cache -- use benchmark_solvers.py for fair timings."""
    ops, geo, settings = _defaults(ops, geo, settings)
    omega = settings.central_omega if omega is None else omega
    n_segments = settings.cell_n_segments if n_segments is None else n_segments
    n_elements = settings.ntu_n_elements if n_elements is None else n_elements

    def timed(solve_fn):
        t_start = time.perf_counter()
        result = solve_fn()
        result.solve_time = time.perf_counter() - t_start
        return result

    return ScenarioResult(
        lmtd=timed(lambda: solve_it_LMTD(omega=omega, ops=ops, geo=geo, settings=settings)),
        ntu=timed(lambda: solve_it_NTU(omega=omega, ops=ops, geo=geo, settings=settings, n_elements=n_elements)),
        cell=timed(lambda: solve_it_cell(n_segments=n_segments, ops=ops, geo=geo, settings=settings)),
    )
