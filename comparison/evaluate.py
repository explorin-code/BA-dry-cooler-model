"""
evaluate.py
============
One operating point -> one flat result row (dict), for every analysis of the cooler comparison.
Builds Geometry/OperatingConditions directly (no parameters.py module state)
and reuses the project's solvers, pressure drops and economics.

Every function takes the run's settings `cfg` (config.Settings) first.

Methods: the five solvers of src/solvers/ -- 'lmtd', 'ntu_field', 'ntu_table' (fast, any flow rate inside the
table range -- otherwise it falls back to the field solve), 'cell_2d', 'cell_3d' (reference). `method=None` means
the comparison's solver (cfg.solver()).
Operating modes: 'dry' (ambient air) or 'precooled' (after the adiabatic
pad, forced). 'auto' = dry unless the chosen method misses the target, then
precooled.
"""

import logging
import warnings

from src import solvers
from src.inputs.cooler_sets import CellSet
from src.inputs.geometry import Geometry
from src.inputs.operating_conditions import OperatingConditions
from src.auxiliaries.pressure_drop import calc_delta_p_air, calc_delta_p_coolant, check_wang_range
from src.thermo.fluid_properties import get_fluid_properties
from src.auxiliaries.economics import calc_fan_power, calc_pump_power, calc_water_usage
from src.thermo.precooling import calc_precooled_ops

from comparison import config

logger = logging.getLogger(__name__)

def make_geometry(gs: CellSet, N_t: int, l: float) -> Geometry:
    return Geometry(d=gs.d_mm / 1e3, delta_f=gs.delta_f_mm / 1e3, s_f=gs.s_f_mm / 1e3,
                    d_i=gs.d_i_mm / 1e3, N_t=int(N_t), N_r=gs.N_r, s_t=gs.s_t_mm / 1e3,
                    s_l=gs.s_l_mm / 1e3, l=l, pipe_material=gs.pipe_material,
                    fin_material=gs.fin_material, delta_c=gs.delta_c_mm / 1e3)


def make_ops(cfg, geo, T_a_i, phi, p_a, w_fr=None, m_dot_c=None, T_c_i=None) -> OperatingConditions:
    """Operating point of a module: ambient air (T_a_i, phi, p_a) at face velocity w_fr (default: the design
    velocity), water m_dot_c / T_c_i (default: the comparison's settings)."""
    w_fr = cfg.w_fr_design if w_fr is None else w_fr
    m_dot_c = cfg.m_dot_c if m_dot_c is None else m_dot_c
    T_c_i = cfg.T_c_i if T_c_i is None else T_c_i
    return OperatingConditions(T_c_i=T_c_i, T_a_i=T_a_i, w_c=0.0, m_dot_c=m_dot_c, V_dot_c=0.0,
                               w_fr=w_fr, m_dot_a=0.0, V_dot_a=0.0, p_a=p_a, phi=phi, geo=geo)


def solve(cfg, method, ops, geo):
    """Runs one solver of src/solvers/ (SOLVER_LABELS): 'lmtd', 'ntu_field', 'ntu_table', 'cell_2d', 'cell_3d'."""
    return solvers.run_solver(method, ops, geo, cfg.solver_settings, profile=False)


def _range_flags(geo):
    """Wang friction-factor range (pressure_drop.WANG_RANGE) and VDI M1
    A/A_p0 range 5-30 (heat transfer)."""
    with warnings.catch_warnings(record=True) as caught:
        warnings.filterwarnings('always', message='Wang correlation')
        check_wang_range(geo)
    ratio = geo.A / geo.A_p0
    return len(caught) == 0, 5.0 <= ratio <= 30.0, ratio


def run_point(cfg, geo, ops, method=None, mode='dry') -> dict:
    """Solves one operating point in the given mode ('dry' or 'precooled')
    and returns the result row. For 'precooled' `ops` is the AMBIENT point;
    the pad is applied here."""
    ops_amb = ops
    if mode == 'precooled':
        ops = calc_precooled_ops(ops_amb, geo)
    elif mode != 'dry':
        raise ValueError(f"mode must be 'dry' or 'precooled', got {mode!r}")

    method = cfg.solver() if method is None else method
    res = solve(cfg, method, ops, geo)
    coolant = get_fluid_properties(ops, ops.T_c_i, ops.p_c)
    Delta_p_c = calc_delta_p_coolant(geo, ops, coolant)
    P_pump = calc_pump_power(ops, coolant, Delta_p_c)
    through_pad = (mode == 'precooled') or cfg.pad_in_dry_path
    Delta_p_a = calc_delta_p_air(geo, ops, res.T_a_o, through_pad)
    P_fan = calc_fan_power(ops, res.T_a_o, Delta_p_a)
    wang_ok, vdi_ok, area_ratio = _range_flags(geo)

    # pinch = of the WHOLE cooler: T_c,o - ambient air (the pad is part of it; negative values are possible when precooled);
    # pinch_coil = of the coil alone: T_c,o - the air entering the coil (= pinch when dry)
    return dict(
        method=method, mode=mode,
        T_a_amb=ops_amb.T_a_i, phi_amb=ops_amb.phi, T_a_i=ops.T_a_i,
        Q_kW=res.Q_dot / 1e3, T_c_o=res.T_c_o, T_a_o=res.T_a_o,
        pinch=res.T_c_o - ops_amb.T_a_i, pinch_coil=res.T_c_o - ops.T_a_i, k=res.k, Re_c=res.diagnostics['Re_c'],
        w_fr=ops.w_fr, m_dot_a=ops.m_dot_a, V_dot_a=ops.V_dot_a, w_c=ops.w_c,
        Delta_p_a=Delta_p_a, Delta_p_c=Delta_p_c,
        P_fan=P_fan, P_pump=P_pump, P_total=P_fan + P_pump,
        m_dot_ev=calc_water_usage(ops_amb, ops) if mode == 'precooled' else 0.0,
        target_met=res.T_c_o <= cfg.T_c_target + config.TARGET_TOL,
        wang_ok=wang_ok, vdi_ok=vdi_ok, A_ratio=area_ratio,
        freeze=res.T_c_o < cfg.T_freeze_risk,   # water at risk of freezing; below 0 °C the numbers are not physical
    )


def run_auto(cfg, geo, ops, method=None) -> dict:
    """Dry if it meets the target, otherwise the precooled result."""
    dry = run_point(cfg, geo, ops, method, 'dry')
    if dry['target_met']:
        return dry
    return run_point(cfg, geo, ops, method, 'precooled')


def calc_masses(geo) -> dict:
    """Material masses [kg] of the fins and the straight tube passes
    (no return bends, headers, collars, frame): fin plate = (s_t s_l - pi d^2/4)
    thick delta_f per fin and tube pass; tube wall = pi/4 (d^2 - d_i^2) l per pass.
    Densities: config.DENSITY (assumed)."""
    passes = geo.N_t * geo.N_r
    m_fin = (geo.A_f / 2) * geo.delta_f * passes * config.DENSITY[geo.fin_material]
    m_tube = (3.141592653589793 / 4) * (geo.d**2 - geo.d_i**2) * geo.l * passes * config.DENSITY[geo.pipe_material]
    by_material = {}
    for material, m in ((geo.fin_material, m_fin), (geo.pipe_material, m_tube)):
        by_material[material] = by_material.get(material, 0.0) + m
    return dict(m_fin=m_fin, m_tube=m_tube, m_total=m_fin + m_tube,
                m_Al=by_material.get('Aluminum', 0.0), m_Cu=by_material.get('Copper', 0.0))


def geometry_summary(geo) -> dict:
    """Footprint and material numbers of one module."""
    depth = geo.N_r * geo.s_l
    return dict(N_t=geo.N_t, N_r=geo.N_r, l=geo.l, width=geo.b, height_l=geo.l, depth=depth,
                A_fr=geo.A_fr, volume=geo.A_fr * depth, A_tot=geo.A_tot,
                N_f=geo.N_f, tubes_total=geo.N_t * geo.N_r, **calc_masses(geo))
