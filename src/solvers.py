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
from functools import partial
import logging
from math import log, pi, exp, isclose
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
    # AI-REVIEW: populated by Cell and NTU (see solve_it_cell/solve_it_NTU) --
    # these are the solver's own final converged grids, not reconstructions.
    # LMTD leaves these None; its spatial profile is reconstructed on demand
    # by analysis.py instead, since it never computes spatial data at all.
    # See CLAUDE.md.
    T_c_grid: object = None
    T_a_grid: object = None
    k_grid: object = None
    dQdL_grid: object = None
    local_diagnostics: dict = None   # Re_air/Nu_air/Re_coolant/Nu_coolant grids -- captured, not yet plotted
    solve_time: float = None         # wall time of the solve [s], set by solve_scenario


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


# =============================================================================
# Solvers 1 & 2: LMTD and NTU
# =============================================================================
# LMTD and NTU are two DISTINCT algorithms -- LMTD marches the log-mean
# temperature difference over the whole array, NTU solves an element-wise
# temperature field (Cabezas-Gomez et al. 2007, see _ntu_step). They share
# an outer relaxation loop (_relax_lmtd_ntu) only because "update fluid
# properties -> k -> new outlet temperatures" is identical between them,
# not because they're secretly the same method. Each solver's own physics
# lives in its own step function (_lmtd_step / _ntu_step) below. If one
# solver ever needs a different relaxation scheme, give it its own loop --
# the driver has no solver-specific knowledge and doesn't need touching.
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

        # Live progress (insight mode) -- useful for spotting slow/diverging runs.
        # %-style args, so nothing is formatted when insight mode is off.
        logger.debug("%-7sHot Error: %.5f | Cold Error: %.5f | T_coolant_out: %.2f | T_air_out: %.2f | omega: %s",
                     tag, diff_hot, diff_cold, T_coolant_out, T_air_out, omega_active)

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


# -----------------------------------------------------------------------------
# NTU: element-wise e-NTU method for multipass counter-cross-flow
# -----------------------------------------------------------------------------
# Source for this whole block: Cabezas-Gomez, Navarro & Saiz-Jabardo (2007),
# "Thermal Performance of Multipass Parallel and Counter-Cross-Flow Heat
# Exchangers", J. Heat Transfer 129, pp. 282-289, DOI 10.1115/1.2430719.
#
# This dry cooler is the paper's G^c_(n,1) arrangement (p. 284, Fig. 3 and
# Section 3: superscript c = counter-cross-flow, first subscript = number of
# tube-fluid passes, second = tube rows per pass), with n = geo.n_rows.
#
# Paper notation -> project:
#   h  = hot tube-side fluid ("Fluid A", mixed)  = coolant
#   c  = cold external fluid ("Fluid B", unmixed) = air
#   N_r = geo.n_rows, N_t = N_c = geo.n_tubes (every tube is one circuit),
#   N_e = n_elements (elements per tube per row)
#   U = k (calc_overall_k, outer-area based), A = geo.A_total
# Source assumptions (p. 283, Section 2): steady state, no heat loss to the
# surroundings, no heat sources, tube fluid mixed, external fluid unmixed,
# constant properties and heat-transfer coefficients, no phase change.
#
# NOTE (circle back if results look off): the source assumes constant U and
# cp. Here they are held constant within one field solve, and updated
# between field solves from the mean temperatures by the shared outer loop
# (_relax_lmtd_ntu) -- the same k-update LMTD uses. A locally varying k
# (as Cell does) would be a departure from the source and is not done.
# -----------------------------------------------------------------------------

def calc_element_effectiveness(UA_e: float, C_c_e: float) -> float:
    """Local element thermal effectiveness Gamma^e [-] of one small
    mixed(tube)-unmixed(external) cross-flow element. Valid for
    C_c^e << C_h^e, i.e. sufficiently small elements.
    Source: Cabezas-Gomez, Navarro & Saiz-Jabardo (2007),
    J. Heat Transfer 129, p. 283, Eq. (1)."""
    return 1.0 - exp(-UA_e / C_c_e)


def calc_element_outlets(T_h_in: float, T_c_in: float, Gamma: float, B: float):
    """Element outlet temperatures (T_c_out, T_h_out) -- paper notation,
    h = coolant, c = air. B = C_c^e * Gamma^e / C_h^e.
    Derived in the source by combining Eq. (1) with the element-mean tube
    temperature T_h^e = 0.5*(T_h,i^e + T_h,o^e) (Eq. (2)) and both fluids'
    energy balances; the two coefficients of each equation sum to 1 and
    C_c*(T_c_out - T_c_in) == C_h*(T_h_in - T_h_out) holds exactly.
    Source: Cabezas-Gomez, Navarro & Saiz-Jabardo (2007),
    J. Heat Transfer 129, p. 283, Eqs. (3)-(4)."""
    # Source: Cabezas-Gomez, Navarro & Saiz-Jabardo (2007),
    # J. Heat Transfer 129, p. 283, Eq. (3).
    T_c_out = ((B + 2.0 * (1.0 - Gamma)) / (2.0 + B)) * T_c_in + (2.0 * Gamma / (2.0 + B)) * T_h_in
    # Source: Cabezas-Gomez, Navarro & Saiz-Jabardo (2007),
    # J. Heat Transfer 129, p. 283, Eq. (4). (No Gamma in the second term.)
    T_h_out = ((2.0 - B) / (2.0 + B)) * T_h_in + (2.0 * B / (2.0 + B)) * T_c_in
    return T_c_out, T_h_out


def solve_ntu_field(n_rows: int, n_elements: int, C_h_circuit: float, C_c_element: float,
                    UA_element: float, T_coolant_in: float, T_air_in: float,
                    tol: float, max_iter: int):
    """Element temperature field of ONE coolant circuit through the G^c_(n,1)
    bundle, constant U and cp. Takes plain numbers only (no geo/ops), so the
    validation (validation_ntu.py) runs exactly this code.

    Returns (T_coolant, T_air, Q_element, T_coolant_out, T_air_out, residuals):
      T_coolant[r, e] -- coolant temperature leaving element (r, e) [degC]
      T_air[r, e]     -- air temperature leaving element (r, e) [degC]
      Q_element[r, e] -- heat transferred in element (r, e) [W]
    Index e is the PHYSICAL axial position along the tube (same in every
    row), not the flow-order index.

    Implementation assumption (not directly specified by source):
    every tube circuit sees the same uniform inlet air and, with the air
    unmixed position by position, behaves identically -- so one circuit is
    solved and represents all n_tubes. The air leaving element (r-1, e)
    enters element (r, e) at the same tube/axial position, without any
    averaging with neighbouring positions (this mapping of the source's
    "external fluid unmixed" to the row/element grid is the project's, not
    a statement of the paper; it reproduces the paper's Table 1 closed forms,
    see validation_ntu.py)."""
    Gamma = calc_element_effectiveness(UA_element, C_c_element)
    # Source: Cabezas-Gomez, Navarro & Saiz-Jabardo (2007),
    # J. Heat Transfer 129, p. 283, text below Eq. (4): B = C_c^e Gamma^e / C_h^e.
    B = C_c_element * Gamma / C_h_circuit

    # Source: p. 283, paragraph before Eq. (1): Eq. (1) requires C_c^e << C_h^e.
    if C_c_element / C_h_circuit > 0.1:
        warnings.warn(f"NTU element field: C_c^e/C_h^e = {C_c_element / C_h_circuit:.3f} is not "
                      f"<< 1 -- increase n_elements (currently {n_elements}).")

    # Coolant flow order: enters the last air-side row, serpentine towards
    # row 0, axial direction alternating per row (U-bend at the same end) --
    # the same routing the Cell solver uses (cell_path_order /
    # get_coolant_inlet), reused here unchanged.
    path = cell_path_order(n_rows, n_elements)

    # Initial field: the source's counter-cross-flow procedure ASSUMES the
    # external-fluid inlet temperatures and iterates (p. 284, Section 2,
    # second-to-last bullet; Fig. 2, p. 283). Assumed here: air unheated
    # everywhere, coolant at its inlet temperature everywhere.
    T_coolant = np.full((n_rows, n_elements), float(T_coolant_in))
    T_air = np.full((n_rows, n_elements), float(T_air_in))
    Q_element = np.zeros((n_rows, n_elements))

    residuals = []
    for _ in range(max_iter):
        T_coolant_old = T_coolant.copy()
        T_air_old = T_air.copy()

        # Source: Fig. 2(b), p. 283 -- march element by element along the
        # tube-fluid circuit, computing both outlet temperatures per element.
        T_h = T_coolant_in
        for r, e in path:
            # Air inlet: row 0 sees the ambient inlet; row r sees the air that
            # left row r-1 at the same position. Rows are visited r_last -> 0,
            # so row r-1 is always taken from the previous sweep -- this lag
            # is what the outer iteration converges away.
            T_c_in = T_air_in if r == 0 else T_air_old[r - 1, e]
            T_c_out, T_h_out = calc_element_outlets(T_h, T_c_in, Gamma, B)

            T_air[r, e] = T_c_out
            T_coolant[r, e] = T_h_out
            Q_element[r, e] = C_h_circuit * (T_h - T_h_out)
            T_h = T_h_out

        # Numerical implementation choice (not a source correlation):
        # max-norm change of both fields. The source's own criterion compares
        # successive MEAN external-fluid outlet temperatures,
        # |T_new - T| / T < tolerance (Fig. 2(a), p. 283).
        residual = max(np.max(np.abs(T_coolant - T_coolant_old)),
                       np.max(np.abs(T_air - T_air_old)))
        residuals.append(residual)
        if residual < tol:
            break
    else:
        raise RuntimeError(f"NTU element field did not converge in {max_iter} sweeps "
                           f"(last residual {residuals[-1]:.2e} K, tolerance {tol:.0e} K).")

    T_coolant_out = T_h                                  # leaves the last element of the circuit
    # Source: Fig. 2(a), item 1.8, p. 283 -- mean external-fluid outlet temperature.
    T_air_out = float(np.mean(T_air[n_rows - 1, :]))

    return T_coolant, T_air, Q_element, T_coolant_out, T_air_out, residuals


def _ntu_step(k, dT_hot, dT_cold, coolant_state, air_state, dm_coolant, dm_air, geo, ops,
              n_elements, field_tol, field_max_iter, field_out=None):
    """NTU solver's physics core: element quantities from the current k and
    cp, one converged element field (solve_ntu_field), outlet temperatures
    and Q. dT_hot/dT_cold are unused (shared step_fn signature).
    field_out (dict, optional) receives the field for the final grids."""
    # Whole-exchanger capacity rates -- paper h = coolant, c = air.
    C_h = dm_coolant * coolant_state.cp
    C_c = dm_air * air_state.cp

    # Element allocation -- Source: Cabezas-Gomez, Navarro & Saiz-Jabardo
    # (2007), J. Heat Transfer 129, p. 283, Fig. 2(a), items 1.5-1.6:
    # (UA)^e = UA/(N_e N_t N_r), C_h^e = C_h/N_c, C_c^e = C_c/(N_e N_t).
    # geo.A is the outer area of one tube in one row, so A_e = geo.A / N_e.
    A_element = geo.A / n_elements
    if not isclose(geo.n_tubes * geo.n_rows * n_elements * A_element, geo.A_total, rel_tol=1e-9):
        raise ValueError("NTU element areas do not add up to geo.A_total.")
    UA_element = k * A_element
    C_h_circuit = C_h / geo.n_tubes
    C_c_element = C_c / (geo.n_tubes * n_elements)

    T_coolant, T_air, Q_element, T_coolant_out, T_air_out, residuals = solve_ntu_field(
        geo.n_rows, n_elements, C_h_circuit, C_c_element, UA_element,
        ops.T_coolant_in, ops.T_air_in, field_tol, field_max_iter,
    )

    # Heat duty from both sides -- equal by construction of Eqs. (3)-(4);
    # a mismatch would mean a bug, so it is checked, not hidden.
    Q_h = C_h * (ops.T_coolant_in - T_coolant_out)
    Q_c = C_c * (T_air_out - ops.T_air_in)
    delta_Q = Q_h - Q_c
    if abs(delta_Q) > 1e-6 * max(abs(Q_h), 1.0):
        raise RuntimeError(f"NTU energy balance mismatch: Q_h = {Q_h:.3f} W, Q_c = {Q_c:.3f} W.")
    dQ = 0.5 * (Q_h + Q_c)

    if field_out is not None:
        C_min, C_max = min(C_h, C_c), max(C_h, C_c)
        field_out.update(
            T_coolant=T_coolant, T_air=T_air, Q_element=Q_element, residuals=residuals,
            delta_Q=delta_Q,
            # Whole-exchanger groups for reporting/validation -- Source: p. 289,
            # Nomenclature (C* = C_min/C_max, NTU = UA/C_min); p. 286, Eqs. (7)-(8)
            # (P = air-side temperature effectiveness, R = C_c/C_h).
            C_star=C_min / C_max,
            NTU=k * geo.A_total / C_min,
            P=(T_air_out - ops.T_air_in) / (ops.T_coolant_in - ops.T_air_in),
            R=C_c / C_h,
        )
    return T_coolant_out, T_air_out, dQ


def solve_it_NTU(omega: float = None, ops=None, geo=None, settings=None, n_elements: int = None):
    """omega: under-relaxation factor of the outer property/k loop; use the
    same value as solve_it_LMTD() for comparable steps. n_elements: elements
    per tube per row (default: settings.ntu_n_elements). Delegates its outer
    loop to _relax_lmtd_ntu(); see _ntu_step()/solve_ntu_field() for the
    NTU-specific physics."""
    if geo is None:
        geo = get_geometry()
    if ops is None:
        ops = get_operating_conditions(geo=geo)
    if settings is None:
        settings = get_solver_settings()
    if omega is None:
        omega = settings.central_omega
    if n_elements is None:
        n_elements = settings.ntu_n_elements

    last_field = {}
    step_fn = partial(_ntu_step, n_elements=n_elements, field_tol=settings.ntu_field_threshold,
                      field_max_iter=settings.ntu_field_max_iter, field_out=last_field)

    (k, dQ, T_coolant_out, T_air_out,
     history_hot, history_cold, history_T_coolant, history_T_air,
     diagnostics) = _relax_lmtd_ntu(step_fn, "[NTU]", omega, ops, geo, settings)

    # --- Final grids, (row, tube, element) like Cell's ---------------------
    # From the last field solve, i.e. the same k returned above. All circuits
    # are identical (see solve_ntu_field), so the one solved circuit is
    # broadcast across n_tubes. The reported outlet temperatures are the
    # outer loop's under-relaxed values and agree with this field to within
    # the outer convergence threshold.
    shape = (geo.n_rows, geo.n_tubes, n_elements)
    T_c_grid = np.broadcast_to(last_field['T_coolant'][:, None, :], shape).copy()   # coolant
    T_a_grid = np.broadcast_to(last_field['T_air'][:, None, :], shape).copy()       # air
    k_grid = np.full(shape, k)
    dL_element = geo.height / n_elements
    dQdL_grid = np.broadcast_to((last_field['Q_element'] / dL_element)[:, None, :], shape).copy()

    return (k, dQ, T_coolant_out, T_air_out,
            history_hot, history_cold, history_T_coolant, history_T_air,
            diagnostics,
            T_c_grid, T_a_grid, k_grid, dQdL_grid,
            None)


# =============================================================================
# Solver 3: Cell-Method
# =============================================================================

def _cell_local_states(r, t, s, T_c, T_a_for_inlet, ops, geo, n_segments):
    """Inlet temperatures and fluid states of one (row, tube, segment) cell.
    Each fluid's properties are evaluated at its OWN local mean temperature,
    (inlet + outlet)/2 across the cell -- coolant at the coolant's, air at
    the air's -- the bulk-temperature convention the tube-side (VDI G1) and
    air-side (VDI M1) correlations are written for, and the same convention
    LMTD/NTU use for the whole exchanger. The outlets are the current
    iterate's values in T_c / T_a_for_inlet, so they settle together with
    the field. (Previously both fluids were evaluated at the average of
    coolant and air inlet temperatures, i.e. the coolant ~4-5 K too cold.)"""
    T_c_in = get_coolant_inlet(r, t, s, n_segments, T_c, ops, geo)
    T_a_in = get_staggered_air_inlet(r, t, s, T_a_for_inlet, ops, geo)

    T_c_mean = (T_c_in + T_c[r, t, s]) / 2.0
    T_a_mean = (T_a_in + T_a_for_inlet[r, t, s]) / 2.0
    props_c = get_fluid_properties(ops, T_c_mean, ops.P_coolant)
    props_a = get_air_properties(ops, T_a_mean, ops.P_air)
    return T_c_in, T_a_in, props_c, props_a


def _solve_cell_step(r, t, s, T_c, T_a_for_inlet, ops, geo, n_segments, dm_coolant, dm_air):
    """Shared per-cell computation for both passes of _relax_cell_grid: the
    local element effectiveness P1 for one (row, tube, segment) cell, given the
    coolant/air inlet temperatures for that cell.

    T_a_for_inlet is which air grid to read the neighbor's inlet from --
    the two passes disagree on this and it must stay a parameter, not a
    hardcoded grid: Pass 1 passes T_a_old (this iteration's air side hasn't
    moved yet), Pass 2 passes the live T_a (already updated by Pass 1 this
    same iteration). That ordering is the staggered-mixing scheme itself,
    not an implementation detail."""
    T_c_in, T_a_in, props_c, props_a = _cell_local_states(
        r, t, s, T_c, T_a_for_inlet, ops, geo, n_segments
    )
    k_local = calc_overall_k(geo, ops, props_c, props_a, T_a_in)

    A_cell = geo.A / n_segments
    W1 = dm_coolant * props_c.cp
    W2 = dm_air * props_a.cp
    R_loc = W1 / W2

    # Per-cell effectiveness: the same element relation the NTU solver uses
    # (one cell = one small mixed(coolant)-unmixed(air) cross-flow element).
    # Source: Cabezas-Gomez, Navarro & Saiz-Jabardo (2007),
    # J. Heat Transfer 129, p. 283, Eqs. (1) and (4):
    #   Gamma = 1 - exp(-UA_cell/W2),  B = W2*Gamma/W1,
    #   coolant side P1 = (T_c,in - T_c,out)/(T_c,in - T_a,in) = 2B/(2+B).
    # The air side then follows as P2 = P1*R_loc = 2*Gamma/(2+B), which is
    # exactly Eq. (3) -- so Pass 2's P1*R_loc stays energy-consistent.
    Gamma = calc_element_effectiveness(k_local * A_cell, W2)
    B = W2 * Gamma / W1
    P1_loc = 2.0 * B / (2.0 + B)

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

        logger.debug("[Cell] Hot Error: %.5f | Cold Error: %.5f | T_coolant_out: %.2f | T_air_out: %.2f | omega: %s",
                     max_raw_hot, max_raw_cold, T_coolant_out, T_air_out, omega)

        if iteration + 1 >= min_iter and max_raw_cold < threshold and max_raw_hot < threshold:
            break

    return (T_c, T_a, dT_hot_it, dT_cold_it, T_coolant_out, T_air_out,
            history_hot, history_cold, history_T_coolant, history_T_air)


def solve_it_cell(n_segments: int = None, omega: float = None, ops=None, geo=None, settings=None):
    """omega: per-cell relaxation factor, default settings.cell_omega (1.0).
    Deliberately not the LMTD/NTU central_omega: at 0.2 the per-cell
    max-raw-change criterion fired before convergence (errors accumulate
    along the coolant path while each damped step is tiny) -- Q came out
    0.17 % high, and the iteration count depended erratically on
    n_segments. At 1.0 Cell converges to the same fixed point a tight
    tolerance gives, in ~9 iterations at every resolution tested (2-100)."""
    # --- 1. Load static configuration -----------------------------------
    if geo is None:
        geo = get_geometry()
    if ops is None:
        ops = get_operating_conditions(geo=geo)
    if settings is None:
        settings = get_solver_settings()
    if omega is None:
        omega = settings.cell_omega
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
    # Whole-exchanger Pr/Re/Nu at the overall mean temperatures, for the
    # output box (comparable with LMTD/NTU). The reported k is NOT taken from
    # here -- it's the mean of the local k Cell actually used, see section 4.
    T_coolant_mean = (ops.T_coolant_in + T_coolant_out) / 2.0
    T_air_mean = (ops.T_air_in + T_air_out) / 2.0
    coolant_state = get_fluid_properties(ops, T_coolant_mean, ops.P_coolant)
    air_state = get_air_properties(ops, T_air_mean, ops.P_air)

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
                T_c_in, T_a_in, props_c, props_a = _cell_local_states(
                    r, t, s, T_c, T_a, ops, geo, n_segments
                )

                k_cell = calc_overall_k(geo, ops, props_c, props_a, T_a_in)
                diag_cell = calc_diagnostics(geo, ops, props_c, props_a, T_a_in)

                k_grid[r, t, s] = k_cell
                dQdL_grid[r, t, s] = k_cell * A_cell * (T_c[r, t, s] - T_a[r, t, s]) / dL_cell
                Re_air_grid[r, t, s] = diag_cell['Re_air']
                Nu_air_grid[r, t, s] = diag_cell['Nu_air']
                Re_coolant_grid[r, t, s] = diag_cell['Re_coolant']
                Nu_coolant_grid[r, t, s] = diag_cell['Nu_coolant']

    # Reported k: area-weighted mean of the local k (all cells have equal
    # area A_cell, so a plain mean), i.e. the k that UA_total / A_total Cell
    # actually used.
    k = float(np.mean(k_grid))

    return (k, dQ, T_coolant_out, T_air_out,
            history_hot, history_cold, history_T_coolant, history_T_air,
            diagnostics,
            T_c, T_a, k_grid, dQdL_grid,
            {'Re_air': Re_air_grid, 'Nu_air': Nu_air_grid,
             'Re_coolant': Re_coolant_grid, 'Nu_coolant': Nu_coolant_grid})


# =============================================================================
# Scenario orchestrator -- runs all three, no printing/plotting
# =============================================================================

def solve_scenario(ops, geo=None, omega: float = None, n_segments: int = None, settings=None,
                   n_elements: int = None) -> ScenarioResult:
    """Runs LMTD/NTU/Cell against the same ops/geo (omega applies to LMTD/NTU;
    Cell uses settings.cell_omega) and records each
    solver's wall time (SolverResult.solve_time). No plotting; progress is
    only logged at DEBUG level (insight mode). Note: the solvers share the
    property cache, so later solvers profit from earlier ones' lookups --
    use benchmark_solvers.py for fair timings."""
    if geo is None:
        geo = get_geometry()
    if settings is None:
        settings = get_solver_settings()
    if omega is None:
        omega = settings.central_omega
    if n_segments is None:
        n_segments = settings.cell_n_segments
    if n_elements is None:
        n_elements = settings.ntu_n_elements

    def timed(solve_fn):
        t_start = time.perf_counter()
        result = SolverResult(*solve_fn())
        result.solve_time = time.perf_counter() - t_start
        return result

    lmtd = timed(lambda: solve_it_LMTD(omega=omega, ops=ops, geo=geo, settings=settings))
    ntu = timed(lambda: solve_it_NTU(omega=omega, ops=ops, geo=geo, settings=settings, n_elements=n_elements))
    cell = timed(lambda: solve_it_cell(n_segments=n_segments, ops=ops, geo=geo, settings=settings))   # own omega

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
