"""
solvers.py
==========
Three solvers for the same dry-cooler problem:

  1. Results         SolverResult, ScenarioResult
  2. Shared helpers  outer relaxation loop (LMTD + NTU); coolant path,
                     element relation (Cabezas-Gomez Eqs. 1, 3, 4), field
                     connections and traversal (NTU + Cell); defaults, states
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
from pathlib import Path
import time
import warnings
import numpy as np

from src.operating_conditions import get_operating_conditions
from src.dry_cooler_physics import get_geometry
from src.fluid_properties import get_fluid_properties, get_air_properties
from src.heat_transfer_core import calc_overall_k, calc_diagnostics
from src.solver_settings import get_solver_settings

logger = logging.getLogger(__name__)

# Convergence: every solver stops when both outlet temperatures change by
# less than settings.convergence_threshold (parameters.CONVERGENCE_THRESHOLD)
# per iteration. Per-solver overrides [K] -- None = the shared value; set a
# number here only to test edge cases.
LMTD_THRESHOLD = None
NTU_THRESHOLD = None               # NTU's outer (property) loop
CELL_THRESHOLD = None
NTU_FIELD_THRESHOLD = None         # NTU's inner element field (max change per sweep);
                                   # None = 1/10 of NTU's threshold, so the field is
                                   # always more accurate than the criterion it feeds


# Oscillation damping of the shared LMTD/NTU outer loop (_relax_lmtd_ntu)
OSC_STALL_RATIO = 0.7      # an oscillation shrinking slower than this per two steps counts as stalled
OMEGA_MIN = 0.01           # lower bound when halving omega on stalled oscillations


def _threshold(override, settings):
    """A solver's convergence threshold: its override, else the shared one."""
    return override if override is not None else settings.convergence_threshold


# =============================================================================
# 1. Results
# =============================================================================

@dataclass
class SolverResult:
    k: float
    Q_dot: float
    T_c_o: float
    T_a_o: float
    history_hot: list                # per-iteration change of Delta_T_hot [K]
    history_cold: list               # per-iteration change of Delta_T_cold [K]
    history_T_c_o: list
    history_T_a_o: list
    diagnostics: dict                # Pr/Re/Nu/alpha at the overall mean temperatures
    # (row, tube, segment) grids -- Cell and NTU only, LMTD leaves them None
    T_c_grid: object = None
    T_a_grid: object = None
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


def _mean_states(ops, T_c_o, T_a_o):
    """Coolant and air states at their mean temperatures (inlet + outlet)/2."""
    coolant_state = get_fluid_properties(ops, (ops.T_c_i + T_c_o) / 2.0, ops.p_c)
    air_state = get_air_properties(ops, (ops.T_a_i + T_a_o) / 2.0, ops.p_a)
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


def calc_element_outlets(T_hot_in: float, T_cold_in: float, Gamma: float, B: float):
    """Element outlet temperatures (T_cold_out, T_hot_out), B = C_c Gamma / C_h.
    Source: Cabezas-Gomez et al. (2007), p. 283, Eqs. (3)-(4)."""
    T_cold_out = ((B + 2.0 * (1.0 - Gamma)) / (2.0 + B)) * T_cold_in + (2.0 * Gamma / (2.0 + B)) * T_hot_in
    T_hot_out = ((2.0 - B) / (2.0 + B)) * T_hot_in + (2.0 * B / (2.0 + B)) * T_cold_in
    return T_cold_out, T_hot_out


# --- Field connections and traversal (NTU + Cell) ----------------------------
# Grids are [row, tube, segment]; NTU solves one circuit, i.e. a single tube
# column. Both solvers pass an element function
#     element(r, t, s, T_c_in, T_a_in, T_c_cur, T_a_cur) -> (T_c_out, T_a_out)
# where T_c_cur/T_a_cur are the element's current outlet values (Cell needs
# them for its local mean temperatures, NTU ignores them). Which traversal a
# solver uses is just the function name called in its iteration loop:
# sweep_two_way (active) or sweep_coolant_direction (inactive, for comparison).

def _coolant_inlet(r, t, s, T_c, T_c_i):
    """Coolant entering element (r, t, s): serpentine via return bends,
    entering at the last row, direction alternating per row (cell_path_order)."""
    N_r, _, n_segments = T_c.shape
    forward = r % 2 == 0                                       # even rows run s = 0 -> end
    s_first = 0 if forward else n_segments - 1
    if s == s_first:                                           # row entry: from the U-bend
        return T_c_i if r == N_r - 1 else T_c[r + 1, t, s]
    return T_c[r, t, s - 1 if forward else s + 1]


def _coolant_predecessors(path):
    """{(r, s): (r, s) of the element upstream on the coolant path, None at
    the circuit inlet} -- the same serpentine as _coolant_inlet, looked up
    once per sweep instead of per element."""
    return dict(zip(path, [None] + path[:-1]))


def _air_inlet(r, t, s, T_a, T_a_i):
    """Air entering element (r, t, s): mean of the two tubes of the row
    upstream it sits between (staggered bank); edge tube, a single tube
    column (NTU, Cell 2D) or row 0: straight through."""
    if r == 0:
        return T_a_i
    t_neighbour = t - 1 if r % 2 == 0 else t + 1
    if 0 <= t_neighbour < T_a.shape[1]:
        return (T_a[r - 1, t, s] + T_a[r - 1, t_neighbour, s]) / 2.0
    return T_a[r - 1, t, s]


def get_coolant_inlet(r, t, s, n_segments, T_c, ops, geo):
    """Coolant entering cell (r, t, s) -- see _coolant_inlet.
    Source: Kröger, Air-Cooled Heat Exchangers and Cooling Towers, Vol. 1, Sec. 3.5."""
    return _coolant_inlet(r, t, s, T_c, ops.T_c_i)


def get_staggered_air_inlet(r, t, s, T_a, ops, geo):
    """Air entering cell (r, t, s) -- see _air_inlet. Modelling assumption
    from the staggered bank (Kröger, Vol. 1, Sec. 5.1)."""
    return _air_inlet(r, t, s, T_a, ops.T_a_i)


def _field_changes(T_c, T_a, T_c_old, T_a_old, omega):
    """Largest raw change per fluid over one iteration [K]. The stored value
    is old + omega (raw - old), so the raw change is the stored one / omega --
    computed once per sweep instead of per element. A diverged field (NaN/inf)
    raises: it would otherwise pass every 'change < tol' test. Divergence
    happens e.g. with omega != 1 per element: each pass is an exact
    sequential solve along its fluid, so per-element relaxation compounds
    along the whole path."""
    change_c = float(np.max(np.abs(T_c - T_c_old))) / omega
    change_a = float(np.max(np.abs(T_a - T_a_old))) / omega
    if not (np.isfinite(change_c) and np.isfinite(change_a)):
        raise RuntimeError("Field traversal diverged (non-finite temperatures) -- check omega.")
    return change_c, change_a


def sweep_two_way(T_c, T_a, element, T_c_i, T_a_i, path, omega=1.0):
    """ACTIVE traversal of NTU and Cell: one iteration over the field, each
    fluid swept in its own flow direction --
      pass 1: coolant along its path, air inlets from the previous iteration;
      pass 2: air row 0 -> last row, with the live air of the row upstream
              and the coolant from pass 1.
    Converges to the same field as sweep_coolant_direction (same element
    relation and connections) in fewer iterations: ~2x at 6 rows, ~3.5x at
    24 (see CLAUDE.md), but evaluates every element twice per iteration.
    Updates T_c/T_a in place (relaxed with omega) and returns the largest
    raw change (coolant, air) [K]."""
    T_c_old, T_a_old = T_c.copy(), T_a.copy()
    N_r, n_tubes, n_segments = T_c.shape
    pred = _coolant_predecessors(path)
    relax = omega != 1.0

    for r, s in path:                                          # pass 1: coolant
        p = pred[r, s]
        for t in range(n_tubes):
            T_c_in = T_c_i if p is None else T_c[p[0], t, p[1]]
            T_a_in = _air_inlet(r, t, s, T_a_old, T_a_i)
            raw, _ = element(r, t, s, T_c_in, T_a_in, T_c[r, t, s], T_a_old[r, t, s])
            T_c[r, t, s] = T_c_old[r, t, s] + omega * (raw - T_c_old[r, t, s]) if relax else raw

    for r in range(N_r):                                       # pass 2: air
        for t in range(n_tubes):
            for s in range(n_segments):
                p = pred[r, s]
                T_c_in = T_c_i if p is None else T_c[p[0], t, p[1]]
                T_a_in = _air_inlet(r, t, s, T_a, T_a_i)
                _, raw = element(r, t, s, T_c_in, T_a_in, T_c[r, t, s], T_a[r, t, s])
                T_a[r, t, s] = T_a_old[r, t, s] + omega * (raw - T_a_old[r, t, s]) if relax else raw
    return _field_changes(T_c, T_a, T_c_old, T_a_old, omega)


def sweep_coolant_direction(T_c, T_a, element, T_c_i, T_a_i, path, omega=1.0):
    """INACTIVE -- not called anywhere; kept for comparison with
    sweep_two_way (swap the function name in solve_ntu_field and/or
    _relax_cell_grid). One pass along the coolant path only, as in
    Cabezas-Gomez et al. (2007), Fig. 2(a), item 1.7: each element updates
    both outlets, with air inlets from the previous iteration (the row
    upstream is only reached later on the coolant path). Same converged
    field, but more iterations, and at the same stopping threshold a larger
    remaining error (convergence is slower).
    Same in-place update and return value as sweep_two_way."""
    T_c_old, T_a_old = T_c.copy(), T_a.copy()
    n_tubes = T_c.shape[1]
    pred = _coolant_predecessors(path)
    relax = omega != 1.0

    for r, s in path:
        p = pred[r, s]
        for t in range(n_tubes):
            T_c_in = T_c_i if p is None else T_c[p[0], t, p[1]]
            T_a_in = _air_inlet(r, t, s, T_a_old, T_a_i)
            raw_c, raw_a = element(r, t, s, T_c_in, T_a_in, T_c[r, t, s], T_a_old[r, t, s])
            T_c[r, t, s] = T_c_old[r, t, s] + omega * (raw_c - T_c_old[r, t, s]) if relax else raw_c
            T_a[r, t, s] = T_a_old[r, t, s] + omega * (raw_a - T_a_old[r, t, s]) if relax else raw_a
    return _field_changes(T_c, T_a, T_c_old, T_a_old, omega)


def _relax_lmtd_ntu(step_fn, tag, omega, threshold, ops, geo, settings, pass_Delta_T=False):
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
        coolant_state, air_state = _mean_states(ops, T_c_o, T_a_o)
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


# =============================================================================
# 3. LMTD
# =============================================================================
# Whole exchanger as pure counterflow (no correction factor F, deliberately
# the simplest of the three methods -- see CLAUDE.md).

def calc_LMTD(Delta_T_hot: float, Delta_T_cold: float) -> float:
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
    Q_dot = geo.A_tot * k * calc_LMTD(Delta_T_hot, Delta_T_cold)
    raw_T_c_o = ops.T_c_i - Q_dot / (ops.m_dot_c * coolant_state.cp)
    raw_T_a_o = ops.T_a_i + Q_dot / (ops.m_dot_a * air_state.cp)
    return raw_T_c_o, raw_T_a_o, Q_dot


def solve_it_LMTD(omega: float = None, ops=None, geo=None, settings=None):
    """omega: under-relaxation factor (default settings.central_omega)."""
    ops, geo, settings = _defaults(ops, geo, settings)
    if omega is None:
        omega = settings.central_omega
    return _relax_lmtd_ntu(_lmtd_step, "[LMTD]", omega, _threshold(LMTD_THRESHOLD, settings), ops, geo, settings, pass_Delta_T=True)


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
# The field is traversed with the shared sweep_two_way (section 2), the same
# as Cell; the paper's own coolant-direction-only march is the inactive
# sweep_coolant_direction (same converged field, more sweeps).

def solve_ntu_field(N_r: int, n_elements: int, C_h_circuit: float, C_c_element: float,
                    UA_element: float, T_c_i: float, T_a_i: float,
                    tol: float, max_iter: int, T_start=None):
    """Converged element field of ONE coolant circuit (all circuits are
    identical). Plain numbers only, so validation_ntu.py runs this exact code.
    Returns (T_c, T_a, Q_element, T_c_o, T_a_o, residuals);
    grids are [row, element] at the physical axial position.
    Assumption (not from the paper): air leaving element (r-1, e) enters
    (r, e), no averaging -- reproduces the paper's closed forms.
    Traversal: sweep_two_way (shared with Cell) instead of the paper's
    coolant-direction-only march (Fig. 2a, item 1.7, = sweep_coolant_direction):
    same converged field, fewer sweeps.
    T_start: optional (T_c, T_a) [row, element] grids to start from (warm
    start from the previous outer iteration); None = the paper's start from
    the inlet temperatures."""
    Gamma = calc_element_effectiveness(UA_element, C_c_element)
    B = C_c_element * Gamma / C_h_circuit                      # p. 283, below Eq. (4)
    if C_c_element / C_h_circuit > 0.1:                        # Eq. (1) needs C_c << C_h
        warnings.warn(f"NTU element field: C_c^e/C_h^e = {C_c_element / C_h_circuit:.3f} is not "
                      f"<< 1 -- increase n_elements (currently {n_elements}).")

    def element(r, t, s, T_c_in, T_a_in, T_c_cur, T_a_cur):   # constant Gamma, B
        T_a_out, T_c_out = calc_element_outlets(T_c_in, T_a_in, Gamma, B)
        return T_c_out, T_a_out

    path = cell_path_order(N_r, n_elements)
    if T_start is None:                                       # one circuit = one tube column
        T_c = np.full((N_r, 1, n_elements), float(T_c_i))
        T_a = np.full((N_r, 1, n_elements), float(T_a_i))     # assumed initial air field
    else:
        T_c = np.array(T_start[0], dtype=float)[:, None, :]
        T_a = np.array(T_start[1], dtype=float)[:, None, :]

    residuals = []
    for _ in range(max_iter):
        change_c, change_a = sweep_two_way(T_c, T_a, element, T_c_i, T_a_i, path)
        # Own criterion (max-norm); the paper compares mean air outlets (Fig. 2a).
        residuals.append(max(change_c, change_a))
        if residuals[-1] < tol:
            break
    else:
        raise RuntimeError(f"NTU element field did not converge in {max_iter} sweeps "
                           f"(last residual {residuals[-1]:.2e} K, tolerance {tol:.0e} K).")

    Q_element = np.array([[C_h_circuit * (_coolant_inlet(r, 0, e, T_c, T_c_i) - T_c[r, 0, e])
                           for e in range(n_elements)] for r in range(N_r)])
    r_last, e_last = path[-1]
    T_c_o = float(T_c[r_last, 0, e_last])                     # after the last element
    T_a_o = float(np.mean(T_a[N_r - 1, 0, :]))                # mean air outlet (Fig. 2a)
    return T_c[:, 0, :], T_a[:, 0, :], Q_element, T_c_o, T_a_o, residuals


def _ntu_step(k, coolant_state, air_state, geo, ops, n_elements, field_tol, field_max_iter, field_out=None):
    """One field solve at the current k and cp; returns outlets and Q.
    Element allocation per Fig. 2(a), items 1.5-1.6. field_out (dict)
    receives the field for the final grids; a field already in it (from the
    previous outer iteration) is the warm start of this solve -- k/cp change
    little between outer iterations, so it is nearly converged already
    (sweeps per field 9, 11, 11, 11, 11 -> 9, 9, 6, 3, 1; NTU ~13 -> ~7 ms)."""
    C_h = ops.m_dot_c * coolant_state.cp                       # paper h = coolant
    C_c = ops.m_dot_a * air_state.cp                           # paper c = air
    A_element = geo.A / n_elements                             # geo.A = one tube, one row
    if not isclose(geo.N_t * geo.N_r * n_elements * A_element, geo.A_tot, rel_tol=1e-9):
        raise ValueError("NTU element areas do not add up to geo.A_tot.")

    T_c, T_a, Q_element, T_c_o, T_a_o, residuals = solve_ntu_field(
        geo.N_r, n_elements, C_h / geo.N_t, C_c / (geo.N_t * n_elements), k * A_element,
        ops.T_c_i, ops.T_a_i, field_tol, field_max_iter,
        T_start=(field_out['T_c'], field_out['T_a']) if field_out and 'T_c' in field_out else None)

    Q_dot = C_h * (ops.T_c_i - T_c_o)                          # coolant side, as Cell

    if field_out is not None:
        C_min, C_max = min(C_h, C_c), max(C_h, C_c)
        field_out.update(T_c=T_c, T_a=T_a, Q_element=Q_element, residuals=residuals,
                         # Paper's whole-exchanger groups (p. 286 Eqs. 7-8, p. 289)
                         C_star=C_min / C_max, NTU=k * geo.A_tot / C_min,
                         P=(T_a_o - ops.T_a_i) / (ops.T_c_i - ops.T_a_i), R=C_c / C_h)
    return T_c_o, T_a_o, Q_dot


# --- NTU operating modes -------------------------------------------------------
# 'field' (default): the element field is solved in every outer iteration --
#     the internal temperature field (NTU profile) comes for free.
# 'table' (multi-run, constant geometry): the arrangement's characteristic
#     P(NTU_a, R) is independent of the inlet temperatures (Cabezas-Gomez et
#     al. 2007, p. 283) and of everything in the geometry except N_r and the
#     element count. It is computed once with the same element field (as the
#     paper's Tables 4/5), stored on disk and interpolated; the outer loop then
#     only looks it up. No internal field unless profile=True (one field solve
#     at the end). Outside the table range the field is solved instead.
NTU_TABLE_VERSION = 1                         # bump when the element field changes -> tables rebuilt
NTU_TABLE_NTU = np.geomspace(0.05, 6.0, 48)   # NTU_a = kA/C_a grid, denser at small NTU (curvature)
NTU_TABLE_R = np.linspace(0.05, 2.0, 40)      # R = C_a/C_c grid (R/N_e <= 0.1 at N_e = 20, Eq. 1)
NTU_TABLE_TOL = 1e-8                          # field tolerance per table point (dimensionless, dT_max = 1)
NTU_TABLE_DIR = Path(__file__).resolve().parent.parent / 'cache' / 'ntu_tables'
_ntu_tables = {}                              # in-memory: (N_r, n_elements) -> interpolator


def build_ntu_table(N_r: int, n_elements: int) -> np.ndarray:
    """P(NTU_a, R) on the NTU_TABLE_NTU x NTU_TABLE_R grid: one dimensionless
    field solve per point (C_a = 1, C_c = 1/R, kA = NTU_a, T_c,i = 1,
    T_a,i = 0, so P = T_a,o) -- the paper's procedure, Fig. 2(a). Each point
    warm-starts from its neighbour in R."""
    P = np.empty((len(NTU_TABLE_NTU), len(NTU_TABLE_R)))
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')                        # element-size warning at large R
        for i, NTU_a in enumerate(NTU_TABLE_NTU):
            start = None
            for j, R in enumerate(NTU_TABLE_R):
                T_c, T_a, _, _, P[i, j], _ = solve_ntu_field(
                    N_r, n_elements, 1.0 / R, 1.0 / n_elements, NTU_a / (N_r * n_elements),
                    1.0, 0.0, NTU_TABLE_TOL, 100000, T_start=start)
                start = (T_c, T_a)
    return P


def get_ntu_table(N_r: int, n_elements: int):
    """Interpolator P(NTU_a, R) for this row and element count: from memory,
    else from disk (NTU_TABLE_DIR), else built and stored."""
    key = (N_r, n_elements)
    if key not in _ntu_tables:
        from scipy.interpolate import RectBivariateSpline
        path = NTU_TABLE_DIR / f"P_Nr{N_r}_Ne{n_elements}_v{NTU_TABLE_VERSION}.npz"
        P = None
        if path.exists():
            data = np.load(path)
            if np.array_equal(data['NTU'], NTU_TABLE_NTU) and np.array_equal(data['R'], NTU_TABLE_R):
                P = data['P']
        if P is None:
            logger.info("[NTU] building P(NTU, R) table for N_r = %d, N_e = %d (once, stored in %s)",
                        N_r, n_elements, path.parent)
            P = build_ntu_table(N_r, n_elements)
            path.parent.mkdir(parents=True, exist_ok=True)
            np.savez(path, P=P, NTU=NTU_TABLE_NTU, R=NTU_TABLE_R)
        _ntu_tables[key] = RectBivariateSpline(NTU_TABLE_NTU, NTU_TABLE_R, P, kx=3, ky=3)
    return _ntu_tables[key]


def _ntu_table_step(k, coolant_state, air_state, geo, ops, table, field_step):
    """One outer step with the table: NTU_a, R -> P -> Q -> both outlets.
    Outside the table range: field_step (an ordinary field solve)."""
    C_c = ops.m_dot_c * coolant_state.cp
    C_a = ops.m_dot_a * air_state.cp
    NTU_a, R = k * geo.A_tot / C_a, C_a / C_c
    if not (NTU_TABLE_NTU[0] <= NTU_a <= NTU_TABLE_NTU[-1] and NTU_TABLE_R[0] <= R <= NTU_TABLE_R[-1]):
        logger.warning("[NTU] NTU_a = %.3f, R = %.3f outside the table -- solving the field", NTU_a, R)
        return field_step(k, coolant_state, air_state, geo, ops)
    Q_dot = float(table.ev(NTU_a, R)) * C_a * (ops.T_c_i - ops.T_a_i)
    return ops.T_c_i - Q_dot / C_c, ops.T_a_i + Q_dot / C_a, Q_dot


def solve_it_NTU(omega: float = None, ops=None, geo=None, settings=None, n_elements: int = None,
                 mode: str = None, profile: bool = True):
    """omega: relaxation of the outer property/k loop (default
    settings.ntu_omega = 1.0, undamped); n_elements per tube and row
    (default settings.ntu_n_elements); mode 'field' or 'table' (default
    settings.ntu_mode, see above). profile: in 'table' mode, solve the field
    once at the converged state for the grids (plots) -- pass False for
    multi-run use, the grids are then None."""
    ops, geo, settings = _defaults(ops, geo, settings)
    if omega is None:
        omega = settings.ntu_omega
    if n_elements is None:
        n_elements = settings.ntu_n_elements
    if mode is None:
        mode = settings.ntu_mode

    last_field = {}
    threshold = _threshold(NTU_THRESHOLD, settings)
    field_tol = NTU_FIELD_THRESHOLD if NTU_FIELD_THRESHOLD is not None else threshold / 10
    field_step = partial(_ntu_step, n_elements=n_elements, field_tol=field_tol,
                         field_max_iter=settings.ntu_field_max_iter, field_out=last_field)
    if mode == 'field':
        step_fn = field_step
    elif mode == 'table':
        step_fn = partial(_ntu_table_step, table=get_ntu_table(geo.N_r, n_elements), field_step=field_step)
    else:
        raise ValueError(f"Unknown NTU mode {mode!r} -- 'field' or 'table'.")
    result = _relax_lmtd_ntu(step_fn, "[NTU]", omega, threshold, ops, geo, settings)

    if mode == 'table' and profile:                            # one field solve for the grids
        field_step(result.k, *_mean_states(ops, result.T_c_o, result.T_a_o), geo, ops)
    if 'T_c' in last_field:
        # Grids (row, tube, element) like Cell's: the one solved circuit, broadcast to all tubes.
        shape = (geo.N_r, geo.N_t, n_elements)
        widen = lambda grid: np.broadcast_to(grid[:, None, :], shape).copy()
        result.T_c_grid = widen(last_field['T_c'])
        result.T_a_grid = widen(last_field['T_a'])
        result.k_grid = np.full(shape, result.k)
        result.dQdL_grid = widen(last_field['Q_element'] / (geo.l / n_elements))
    return result


# =============================================================================
# 5. Cell -- VDI cell method with local properties and k
# =============================================================================
# Method: VDI-Wärmeatlas (2019), C1 §3.1. Grid [row, tube, segment]; each
# cell is one small mixed(coolant)-unmixed(air) cross-flow element with the
# Cabezas-Gomez element relation (as NTU) but local k and cp, traversed with
# the shared sweep_two_way (section 2). Coolant/air connections: _coolant_inlet,
# _air_inlet (staggered neighbour averaging).
# 3D: all N_t tubes. 2D (settings.cell_2d): one representative tube -- the
# grid has a single tube column, so no neighbour averaging of air; the
# result is broadcast to N_t tubes for output. The tube count of a run is
# always taken from the grid shape, not from geo.

def _cell_local_states(T_c_in, T_a_in, T_c_cur, T_a_cur, ops):
    """Local mean air temperature and fluid states of one cell; each fluid at
    its own local mean temperature (inlet + current outlet)/2 (thesis Eq.
    mean_temperatures, applied per cell)."""
    T_a_m = (T_a_in + T_a_cur) / 2.0
    props_c = get_fluid_properties(ops, (T_c_in + T_c_cur) / 2.0, ops.p_c)
    props_a = get_air_properties(ops, T_a_m, ops.p_a)
    return T_a_m, props_c, props_a


def _cell_element(r, t, s, T_c_in, T_a_in, T_c_cur, T_a_cur,
                  ops, geo, n_segments, m_dot_c_tube, m_dot_a_cell):
    """Outlets (T_c_out, T_a_out) of one cell: the same Cabezas-Gomez element
    relation as NTU (Eqs. 1, 3, 4), with Gamma and B from the cell's local k
    and cp instead of one value for the whole field."""
    T_a_m, props_c, props_a = _cell_local_states(T_c_in, T_a_in, T_c_cur, T_a_cur, ops)
    k_local = calc_overall_k(geo, ops, props_c, props_a, T_a_m)

    C_dot_c = m_dot_c_tube * props_c.cp
    C_dot_a = m_dot_a_cell * props_a.cp
    Gamma = calc_element_effectiveness(k_local * (geo.A / n_segments), C_dot_a)
    B = C_dot_a * Gamma / C_dot_c
    T_a_out, T_c_out = calc_element_outlets(T_c_in, T_a_in, Gamma, B)
    return T_c_out, T_a_out


def _relax_cell_grid(T_c, T_a, omega, Delta_T_hot, Delta_T_cold,
                     geo, ops, n_segments, m_dot_c_tube, m_dot_a_cell, threshold, max_iter):
    """Iterates the cell field with sweep_two_way (shared with NTU) until both
    outlet temperatures change by less than threshold per iteration -- the
    same criterion as LMTD/NTU's outer loop."""
    history_hot, history_cold, history_T_c_o, history_T_a_o = [], [], [], []
    path = cell_path_order(geo.N_r, n_segments)
    element = partial(_cell_element, ops=ops, geo=geo, n_segments=n_segments,
                      m_dot_c_tube=m_dot_c_tube, m_dot_a_cell=m_dot_a_cell)

    for iteration in range(max_iter):
        change_c, change_a = sweep_two_way(T_c, T_a, element, ops.T_c_i, ops.T_a_i, path, omega)

        T_c_o = float(np.mean(T_c[0, :, n_segments - 1]))   # coolant leaves row 0
        T_a_o = float(np.mean(T_a[geo.N_r - 1, :, :]))       # air leaves the last row

        new_Delta_T_hot = ops.T_c_i - T_a_o
        new_Delta_T_cold = T_c_o - ops.T_a_i
        diff_hot = new_Delta_T_hot - Delta_T_hot
        diff_cold = new_Delta_T_cold - Delta_T_cold
        Delta_T_hot, Delta_T_cold = new_Delta_T_hot, new_Delta_T_cold
        history_hot.append(diff_hot)
        history_cold.append(diff_cold)
        history_T_c_o.append(T_c_o)
        history_T_a_o.append(T_a_o)

        logger.debug("[Cell] Air change: %.5f | Coolant change: %.5f | T_coolant_out: %.2f | T_air_out: %.2f | omega: %s",
                     change_a, change_c, T_c_o, T_a_o, omega)
        if abs(diff_hot) < threshold and abs(diff_cold) < threshold:
            break
    else:
        raise RuntimeError(f"[Cell] did not converge in {max_iter} iterations "
                           f"(last changes {diff_hot:.2e} / {diff_cold:.2e} K, threshold {threshold:.0e} K) -- "
                           f"cap too low, or threshold below the property-cache floor (~1.5e-5 K)?")

    return (T_c, T_a, T_c_o, T_a_o,
            history_hot, history_cold, history_T_c_o, history_T_a_o)


def _cell_local_grids(T_c, T_a, ops, geo, n_segments):
    """k, dQ/dL and Re/Nu per cell on the converged grids (same local states
    as during the iteration). dQ = k dA dT: Baehr, Thermodynamik (2016), Sec. 3.1."""
    shape = T_c.shape
    k_grid, dQdL_grid = np.zeros(shape), np.zeros(shape)
    local = {key: np.zeros(shape) for key in ('Re_a', 'Nu_a', 'Re_c', 'Nu_c')}
    A_cell = geo.A / n_segments
    dL_cell = geo.l / n_segments

    for r in range(geo.N_r):
        for t in range(shape[1]):
            for s in range(n_segments):
                T_a_m, props_c, props_a = _cell_local_states(
                    _coolant_inlet(r, t, s, T_c, ops.T_c_i), _air_inlet(r, t, s, T_a, ops.T_a_i),
                    T_c[r, t, s], T_a[r, t, s], ops)
                k_grid[r, t, s] = calc_overall_k(geo, ops, props_c, props_a, T_a_m)
                dQdL_grid[r, t, s] = k_grid[r, t, s] * A_cell * (T_c[r, t, s] - T_a[r, t, s]) / dL_cell
                diag = calc_diagnostics(geo, ops, props_c, props_a, T_a_m)
                for key in local:
                    local[key][r, t, s] = diag[key]
    return k_grid, dQdL_grid, local


def solve_it_cell(n_segments: int = None, omega: float = None, ops=None, geo=None, settings=None):
    """omega: per-cell relaxation (default settings.cell_omega = 1.0; with
    omega < 1 the change-per-iteration criterion fires before convergence --
    see CLAUDE.md). n_segments per tube and row (default settings.cell_n_segments)."""
    ops, geo, settings = _defaults(ops, geo, settings)
    if omega is None:
        omega = settings.cell_omega
    if n_segments is None:
        n_segments = settings.cell_n_segments

    m_dot_c_tube = ops.m_dot_c / geo.N_t                       # per tube, also in 2D
    m_dot_a_cell = ops.m_dot_a / (geo.N_t * n_segments)
    n_tubes = 1 if settings.cell_2d else geo.N_t
    T_c = np.full((geo.N_r, n_tubes, n_segments), ops.T_c_i)
    T_a = np.full((geo.N_r, n_tubes, n_segments), ops.T_a_i)

    (T_c, T_a, T_c_o, T_a_o,
     history_hot, history_cold, history_T_c_o, history_T_a_o) = _relax_cell_grid(
        T_c, T_a, omega, settings.Delta_T_hot_init, settings.Delta_T_cold_init,
        geo, ops, n_segments, m_dot_c_tube, m_dot_a_cell,
        threshold=_threshold(CELL_THRESHOLD, settings), max_iter=settings.cell_max_iter)

    # Whole-exchanger Q and diagnostics at the overall mean temperatures; the
    # reported k is the mean of the local k Cell actually used (equal cell areas).
    coolant_state, air_state = _mean_states(ops, T_c_o, T_a_o)
    Q_dot = ops.m_dot_c * coolant_state.cp * (ops.T_c_i - T_c_o)
    diagnostics = calc_diagnostics(geo=geo, ops=ops, coolant_state=coolant_state, air_state=air_state,
                                   T_a_m=(ops.T_a_i + T_a_o) / 2.0)
    k_grid, dQdL_grid, local_diagnostics = _cell_local_grids(T_c, T_a, ops, geo, n_segments)

    if n_tubes == 1:                                           # 2D: broadcast to all tubes for output
        widen = lambda grid: np.broadcast_to(grid, (geo.N_r, geo.N_t, n_segments)).copy()
        T_c, T_a, k_grid, dQdL_grid = map(widen, (T_c, T_a, k_grid, dQdL_grid))
        local_diagnostics = {key: widen(grid) for key, grid in local_diagnostics.items()}

    return SolverResult(float(np.mean(k_grid)), Q_dot, T_c_o, T_a_o,
                        history_hot, history_cold, history_T_c_o, history_T_a_o, diagnostics,
                        T_c_grid=T_c, T_a_grid=T_a, k_grid=k_grid, dQdL_grid=dQdL_grid,
                        local_diagnostics=local_diagnostics)


# =============================================================================
# 6. Scenario
# =============================================================================

def solve_scenario(ops, geo=None, omega: float = None, n_segments: int = None, settings=None,
                   n_elements: int = None) -> ScenarioResult:
    """Runs LMTD, NTU and Cell on one operating point and times each. omega
    applies to LMTD only; NTU uses settings.ntu_omega, Cell settings.cell_omega. The solvers share the
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
        ntu=timed(lambda: solve_it_NTU(ops=ops, geo=geo, settings=settings, n_elements=n_elements)),
        cell=timed(lambda: solve_it_cell(n_segments=n_segments, ops=ops, geo=geo, settings=settings)),
    )
