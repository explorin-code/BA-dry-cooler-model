"""
ntu_table.py -- NTU 'table' mode: the arrangement's characteristic P(NTU_a, R) computed once with the element field
(like the paper's Tables 4 / 5), stored on disk (cache/ntu_tables/, gitignored), interpolated in the outer loop.

P is independent of the inlet temperatures (Cabezas-Gomez et al. 2006, p. 283) and of everything in the geometry
except N_r and the element count. Outside the table range the field is solved instead (ntu.py).
"""

import logging
from pathlib import Path
import warnings

import numpy as np

from src.solvers.ntu_field import solve_ntu_field

logger = logging.getLogger(__name__)

NTU_TABLE_VERSION = 2                         # bump when the element field changes -> tables rebuilt
NTU_TABLE_NTU = np.geomspace(0.05, 60.0, 64)  # NTU_a = kA/C_a grid, denser at small NTU (curvature); the top covers the slow fan of the adaptive control
NTU_TABLE_NTU_BENCHMARK = np.geomspace(0.05, 6.0, 48)   # the range of the first table (NTU_a <= 6): the SOLVER BENCHMARKS build their tables over this one, so their
                                              # build times stay those reported in the thesis (a build over NTU_TABLE_NTU costs ~3x as much; the lookups do not differ)
NTU_TABLE_R_PER_ELEMENT = 0.1                 # Eq. (1) validity: R / N_e <= 0.1 -> R grid reaches 0.1 N_e
NTU_TABLE_R_DENSITY = 20                      # R grid points per unit R (40 points up to R = 2 at N_e = 20)
NTU_TABLE_TOL = 1e-8                          # field tolerance per table point (dimensionless, dT_max = 1)
NTU_TABLE_DIR = Path(__file__).resolve().parent.parent.parent / 'cache' / 'ntu_tables'
tables = {}                                   # in memory: (N_r, n_elements) -> interpolator


def ntu_table_R(n_elements: int) -> np.ndarray:
    """R = C_a/C_c axis of the table: 0.05 up to the largest R that Eq. (1)
    allows for this element count (R/N_e <= 0.1) -- 0.05..2 at N_e = 20.
    More elements = a larger flow-rate range (air capacity vs. coolant)."""
    R_max = NTU_TABLE_R_PER_ELEMENT * n_elements
    return np.linspace(0.05, R_max, round(NTU_TABLE_R_DENSITY * R_max))


def build_ntu_table(N_r: int, n_elements: int, axis: np.ndarray = NTU_TABLE_NTU) -> np.ndarray:
    """P(NTU_a, R) on the `axis` (default NTU_TABLE_NTU) x ntu_table_R(n_elements) grid: one dimensionless
    field solve per point (C_a = 1, C_c = 1/R, kA = NTU_a, T_c,i = 1,
    T_a,i = 0, so P = T_a,o) -- the paper's procedure, Fig. 2(a). Each point
    warm-starts from its neighbour in R."""
    R_axis = ntu_table_R(n_elements)
    P = np.empty((len(axis), len(R_axis)))
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')                        # element-size warning at large R
        for i, NTU_a in enumerate(axis):
            start = None
            for j, R in enumerate(R_axis):
                T_c, T_a, _, _, P[i, j], _ = solve_ntu_field(
                    N_r, n_elements, 1.0 / R, 1.0 / n_elements, NTU_a / (N_r * n_elements),
                    1.0, 0.0, NTU_TABLE_TOL, 100000, T_start=start)
                start = (T_c, T_a)
    return P


def get_ntu_table(N_r: int, n_elements: int):
    """Interpolator P(NTU_a, R) for this row and element count: from memory,
    else from disk (NTU_TABLE_DIR), else built and stored."""
    key = (N_r, n_elements)
    if key not in tables:
        path = NTU_TABLE_DIR / f"P_Nr{N_r}_Ne{n_elements}_v{NTU_TABLE_VERSION}.npz"
        P = None
        if path.exists():
            data = np.load(path)
            if np.array_equal(data['NTU'], NTU_TABLE_NTU) and np.array_equal(data['R'], ntu_table_R(n_elements)):
                P = data['P']
        if P is None:
            logger.info("[NTU] building P(NTU, R) table for N_r = %d, N_e = %d (once, stored in %s)",
                        N_r, n_elements, path.parent)
            P = build_ntu_table(N_r, n_elements)
            path.parent.mkdir(parents=True, exist_ok=True)
            np.savez(path, P=P, NTU=NTU_TABLE_NTU, R=ntu_table_R(n_elements))
        install_table(N_r, n_elements, P)
    return tables[key]


def install_table(N_r: int, n_elements: int, P: np.ndarray, axis: np.ndarray = NTU_TABLE_NTU):
    """Puts a table (as built by build_ntu_table over the same `axis`) into the in-memory tables of this run (benchmarks, which build their own)."""
    from scipy.interpolate import RectBivariateSpline
    tables[(N_r, n_elements)] = RectBivariateSpline(axis, ntu_table_R(n_elements), P, kx=3, ky=3)


def ntu_table_step(k, coolant_state, air_state, geo, ops, table, field_step):
    """One outer step with the table: NTU_a, R -> P -> Q -> both outlets.
    Outside the table range: field_step (an ordinary field solve)."""
    C_c = ops.m_dot_c * coolant_state.cp
    C_a = ops.m_dot_a * air_state.cp
    NTU_a, R = k * geo.A_tot / C_a, C_a / C_c
    (NTU_lo, NTU_hi), (R_lo, R_hi) = (axis[[0, -1]] for axis in table.get_knots())
    if not (NTU_lo <= NTU_a <= NTU_hi and R_lo <= R <= R_hi):
        logger.warning("[NTU] NTU_a = %.3f, R = %.3f outside the table -- solving the field", NTU_a, R)
        return field_step(k, coolant_state, air_state, geo, ops)
    Q_dot = float(table.ev(NTU_a, R)) * C_a * (ops.T_c_i - ops.T_a_i)
    return ops.T_c_i - Q_dot / C_c, ops.T_a_i + Q_dot / C_a, Q_dot
