"""
profiles.py
============
Reconstructs a comparable spatial profile (temperature, k, dQ/dL along the
coolant's flow path) for each solver from its already-converged output.
Cell's and NTU's profiles are genuine reductions of the element/cell grids
they already computed. LMTD's is a post-hoc reconstruction from its own
converged scalars -- NOT something it tracks internally. See docs/model-notes.md.
"""

from dataclasses import dataclass
import numpy as np

from src.thermo.fluid_properties import get_fluid_properties, get_air_properties
from src.solvers import cell_path_order


@dataclass
class ProfileData:
    x_frac: object      # 0 (coolant inlet) to 1 (coolant outlet)
    T_c: object
    T_a: object
    k: object
    dQdL: object


def calc_lmtd_profile(lmtd_result, ops, geo, n_points: int = 200) -> ProfileData:
    """The exponential T(x) profile implied by LMTD's log-mean assumption, reconstructed in closed form from its converged
    scalars -- not something solve_lmtd computes (docs/model-notes.md)."""
    Delta_T_hot = ops.T_c_i - lmtd_result.T_a_o
    Delta_T_cold = lmtd_result.T_c_o - ops.T_a_i

    T_c_m = (ops.T_c_i + lmtd_result.T_c_o) / 2.0
    T_a_m = (ops.T_a_i + lmtd_result.T_a_o) / 2.0
    coolant_state = get_fluid_properties(ops, T_c_m, ops.p_c)
    air_state = get_air_properties(ops, T_a_m, ops.p_a)
    C_dot_c = ops.m_dot_c * coolant_state.cp
    C_dot_a = ops.m_dot_a * air_state.cp

    x_frac = np.linspace(0.0, 1.0, n_points)
    ratio = Delta_T_cold / Delta_T_hot
    delta_T = Delta_T_hot * ratio**x_frac
    Q = Delta_T_hot * (1 - ratio**x_frac) / (1 / C_dot_c - 1 / C_dot_a)

    T_c = ops.T_c_i - Q / C_dot_c
    T_a = lmtd_result.T_a_o - Q / C_dot_a
    k_profile = np.full(n_points, lmtd_result.k)

    # Local dQ = k dA (T_B - T_A): Baehr, Thermodynamik, 16th ed. (2016), Sec. 3.1.
    dQdL = lmtd_result.k * delta_T * (geo.A_tot / geo.l_tot)

    return ProfileData(x_frac=x_frac, T_c=T_c, T_a=T_a, k=k_profile, dQdL=dQdL)


def calc_ntu_profile(ntu_result, geo, n_elements: int) -> ProfileData:
    """NTU's profile: a reduction of its own converged element field along the coolant path the solver marches
    (cell_path_order). Same reduction as calc_cell_profile; all tubes are identical for NTU, so the mean / sum over tubes is exact."""
    path = cell_path_order(geo.N_r, n_elements)
    n_points = len(path)

    T_c = np.array([ntu_result.T_c_grid[r, :, e].mean() for r, e in path])
    T_a = np.array([ntu_result.T_a_grid[r, :, e].mean() for r, e in path])
    k = np.array([ntu_result.k_grid[r, :, e].mean() for r, e in path])
    dQdL = np.array([ntu_result.dQdL_grid[r, :, e].sum() for r, e in path])   # extensive: sum over tubes

    x_frac = np.linspace(0.0, 1.0, n_points)

    return ProfileData(x_frac=x_frac, T_c=T_c, T_a=T_a, k=k, dQdL=dQdL)


def calc_cell_profile(cell_result, geo, n_segments: int) -> ProfileData:
    """Cell's profile: a reduction of its converged grids (not a reconstruction), along solvers.cell_path_order -- the same
    path the shared field traversal (sweep_two_way, pass 1) walks."""
    path = cell_path_order(geo.N_r, n_segments)
    n_points = len(path)

    T_c = np.array([cell_result.T_c_grid[r, :, s].mean() for r, s in path])
    T_a = np.array([cell_result.T_a_grid[r, :, s].mean() for r, s in path])
    k = np.array([cell_result.k_grid[r, :, s].mean() for r, s in path])
    # dQdL is extensive (per-tube rates from parallel tubes add, not average) --
    # sum across tubes here, unlike the mean() used for the intensive quantities
    # above, so this matches LMTD/NTU's dQdL (both already scaled by N_t).
    dQdL = np.array([cell_result.dQdL_grid[r, :, s].sum() for r, s in path])

    x_frac = np.linspace(0.0, 1.0, n_points)

    return ProfileData(x_frac=x_frac, T_c=T_c, T_a=T_a, k=k, dQdL=dQdL)
