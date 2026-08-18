"""
analysis.py
============
Reconstructs a comparable spatial profile (temperature, k, dQ/dL along the
coolant's flow path) for each solver from its already-converged output.
Cell's profile is a genuine reduction of data it already computed. LMTD's
and NTU's are post-hoc reconstructions from their own converged scalars --
NOT something either solver actually tracks internally. See CLAUDE.md.
"""

from dataclasses import dataclass
import numpy as np

from src.fluid_properties import get_fluid_properties, get_air_properties
from src.solvers import calc_P


@dataclass
class ProfileData:
    x_frac: object      # 0 (coolant inlet) to 1 (coolant outlet)
    T_coolant: object
    T_air: object
    k: object
    dQdL: object


# AI-REVIEW: closed-form reconstruction of LMTD's implied exponential T(x)
# profile, derived from the log-mean assumption itself -- not something
# solve_it_LMTD computes. See CLAUDE.md.
def calc_lmtd_profile(lmtd_result, ops, geo, n_points: int = 200) -> ProfileData:
    dT_hot = ops.T_coolant_in - lmtd_result.T_air_out
    dT_cold = lmtd_result.T_coolant_out - ops.T_air_in

    T_coolant_mean = (ops.T_coolant_in + lmtd_result.T_coolant_out) / 2.0
    T_air_mean = (ops.T_air_in + lmtd_result.T_air_out) / 2.0
    coolant_state = get_fluid_properties(ops, T_coolant_mean, ops.P_coolant)
    air_state = get_air_properties(ops, T_air_mean, ops.P_air)
    W_coolant = ops.m_coolant * coolant_state.cp
    W_air = ops.m_o * air_state.cp

    x_frac = np.linspace(0.0, 1.0, n_points)
    ratio = dT_cold / dT_hot
    delta_T = dT_hot * ratio**x_frac
    Q = dT_hot * (1 - ratio**x_frac) / (1 / W_coolant - 1 / W_air)

    T_coolant = ops.T_coolant_in - Q / W_coolant
    T_air = lmtd_result.T_air_out - Q / W_air
    k_profile = np.full(n_points, lmtd_result.k)

    A_total = geo.A * geo.n_tubes * geo.n_rows
    dQdL = lmtd_result.k * delta_T * (A_total / geo.l)

    return ProfileData(x_frac=x_frac, T_coolant=T_coolant, T_air=T_air, k=k_profile, dQdL=dQdL)


# AI-REVIEW: row-by-row reconstruction by cascading NTU's own (constant)
# per-row P1 forward from the known inlet conditions -- not something
# solve_it_NTU computes; it only ever evaluates the closed-form combination.
# Self-checked against that closed form below. See CLAUDE.md.
def calc_ntu_profile(ntu_result, ops, geo) -> ProfileData:
    n = geo.n_rows
    T_coolant_mean = (ops.T_coolant_in + ntu_result.T_coolant_out) / 2.0
    T_air_mean = (ops.T_air_in + ntu_result.T_air_out) / 2.0
    coolant_state = get_fluid_properties(ops, T_coolant_mean, ops.P_coolant)
    air_state = get_air_properties(ops, T_air_mean, ops.P_air)

    W1 = ops.m_coolant * coolant_state.cp
    W2 = ops.m_o * air_state.cp
    R1 = W1 / W2
    A_run = geo.A * geo.n_tubes
    NTU1 = (ntu_result.k * A_run) / W1
    P1 = calc_P(NTU1, R1)
    P2 = P1 * R1

    def march(T_a0_guess):
        """Marches coolant-direction (j=0 coolant inlet -> j=n coolant outlet),
        given a guessed air temperature at the coolant-inlet boundary (T_a[0] --
        physically air's OUTLET there, since this is counterflow). T_a[j+1] is
        solved from row j's own air-side relation rather than assumed known,
        which is what makes this a genuine counterflow march instead of the
        parallel-flow one from before."""
        T_c = np.zeros(n + 1)
        T_a = np.zeros(n + 1)
        T_c[0] = ops.T_coolant_in
        T_a[0] = T_a0_guess
        for j in range(n):
            T_a[j + 1] = (T_a[j] - P2 * T_c[j]) / (1 - P2)
            T_c[j + 1] = T_c[j] - P1 * (T_c[j] - T_a[j + 1])
        return T_c, T_a

    # AI-REVIEW: linear shooting -- exact for this affine per-row map, not an
    # approximation. Two trial guesses for the unknown air-outlet boundary
    # pin down the guess -> T_a[n] relationship exactly, then a third march
    # with the solved guess gives the real profile. See CLAUDE.md.
    _, T_a_trial0 = march(0.0)
    _, T_a_trial1 = march(1.0)
    b0, b1 = T_a_trial0[-1], T_a_trial1[-1]
    T_a0_solved = (ops.T_air_in - b0) / (b1 - b0)
    T_c, T_a = march(T_a0_solved)

    if abs(T_c[-1] - ntu_result.T_coolant_out) > 0.5:
        print(f"[analysis] WARNING: NTU row-recursion reconstruction "
              f"(T_coolant_out={T_c[-1]:.2f}) doesn't match the closed-form "
              f"result ({ntu_result.T_coolant_out:.2f}). See CLAUDE.md.")

    x_frac = np.linspace(0.0, 1.0, n + 1)
    k_profile = np.full(n + 1, ntu_result.k)
    dQdL = ntu_result.k * (T_c - T_a) * (A_run / geo.height)   # both arrays now share the same boundary indexing

    return ProfileData(x_frac=x_frac, T_coolant=T_c, T_air=T_a, k=k_profile, dQdL=dQdL)


def _coolant_path_order(n_rows, n_segments):
    """(r, s) coordinates in coolant flow order -- mirrors _relax_cell_grid's
    own Pass 1 traversal exactly, so this is provably the same path the
    solver itself walks, not a separate guess at it."""
    path = []
    for r in range(n_rows - 1, -1, -1):
        s_range = range(0, n_segments) if r % 2 == 0 else range(n_segments - 1, -1, -1)
        for s in s_range:
            path.append((r, s))
    return path


# Cell's profile is a genuine reduction, not a reconstruction -- see
# solve_it_cell's own diagnostic pass for where this data actually comes from.
def calc_cell_profile(cell_result, geo, n_segments: int) -> ProfileData:
    path = _coolant_path_order(geo.n_rows, n_segments)
    n_points = len(path)

    T_coolant = np.array([cell_result.T_c_grid[r, :, s].mean() for r, s in path])
    T_air = np.array([cell_result.T_a_grid[r, :, s].mean() for r, s in path])
    k = np.array([cell_result.k_grid[r, :, s].mean() for r, s in path])
    # dQdL is extensive (per-tube rates from parallel tubes add, not average) --
    # sum across tubes here, unlike the mean() used for the intensive quantities
    # above, so this matches LMTD/NTU's dQdL (both already scaled by n_tubes).
    dQdL = np.array([cell_result.dQdL_grid[r, :, s].sum() for r, s in path])

    x_frac = np.linspace(0.0, 1.0, n_points)

    return ProfileData(x_frac=x_frac, T_coolant=T_coolant, T_air=T_air, k=k, dQdL=dQdL)
