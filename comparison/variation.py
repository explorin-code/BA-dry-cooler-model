"""
variation.py -- one-at-a-time geometry variations with the annual economics of the adaptive fan.

Each parameter of RANGES is varied alone (the cooler's own value is always one of the points); for every variant and cooler:
  * the pinch (T_c,o - ambient) at the design state, dry, fixed face velocity (what the correlations are good for),
  * the climate year with the adaptive fan and the pad (operation.control_point, every STRIDE-th hour): fan and pump energy, pad water,
    hours at the limit (target missed), and the metal masses for the write-off (prices at plot time: add_costs),
  * whether the variant is inside the validity ranges of Wang (pressure drop) and VDI M1 (heat transfer), evaluate._range_flags.
The ranges are wider than the literature ranges of config.RANGES on purpose (the variants outside are marked in the figure).

Needs a solver with a table or a fast field (the comparison's 'ntu_table': the tables of every row count in RANGES['N_r'] are built / loaded once up front,
then the coolers are solved in parallel processes). STRIDE = 6 stays within 1 % of the full year for the standard coolers.
"""

import dataclasses
from concurrent.futures import ProcessPoolExecutor
from itertools import repeat

import numpy as np
import pandas as pd

from comparison import evaluate as ev, operation, sizing
from src.inputs import weather
from src.quiet import quiet_solvers

STRIDE = 6                                    # every n-th hour of the year (a coarser cfg.year_stride wins)
RANGES = {                                    # parameter: (low, high, points); N_r: every integer from low to high
    'N_r': (1, 10, 10), 'N_t': (5, 100, 12), 'l': (0.25, 5.0, 10),
    's_l_mm': (10.0, 40.0, 10), 's_f_mm': (1.0, 10.0, 10), 'delta_f_mm': (0.05, 0.4, 10)}


def _values(param, base, ranges):
    low, high, n = ranges[param]
    if param == 'N_r':
        values = range(int(low), int(high) + 1)
    elif param == 'N_t':
        values = np.unique(np.round(np.linspace(low, high, n)).astype(int)).tolist()
    else:
        values = np.round(np.linspace(low, high, n), 3).tolist()
    return sorted(set(values) | {base[param]})


def _is_geometry(gs):
    """False when the tubes would touch, the fins overlap or no bore is left."""
    return not (gs.d_mm > min(gs.s_t_mm, float(np.hypot(gs.s_t_mm / 2, gs.s_l_mm))) - 0.5 or gs.s_f_mm <= gs.delta_f_mm or gs.d_i_mm <= 0.5)


def _variant(base_set, base_geo, param, value):
    """(GeometrySet, N_t, l) with `param` set to `value`."""
    gs, N_t, l = base_set, base_geo.N_t, base_geo.l
    if param == 'N_r':
        gs = dataclasses.replace(base_set, N_r=int(value))
    elif param == 'N_t':
        N_t = int(value)
    elif param == 'l':
        l = float(value)
    else:
        gs = dataclasses.replace(base_set, **{param: float(value)})
    return gs, N_t, l


def _run_set(cfg, key, stride, ranges) -> pd.DataFrame:
    """All variants of cooler `key` (runs in a worker process)."""
    method, pt = cfg.solver(), cfg.design_point()
    states = list(weather.hourly_states(weather.load_year(), stride))
    base_set, base_geo = cfg.cell(key), sizing.get_sized(cfg, key)
    base = dict(N_r=base_set.N_r, N_t=base_geo.N_t, l=base_geo.l, s_l_mm=base_set.s_l_mm, s_f_mm=base_set.s_f_mm, delta_f_mm=base_set.delta_f_mm)
    rows = []
    with quiet_solvers():
        for param in ranges:
            for value in _values(param, base, ranges):
                gs, N_t, l = _variant(base_set, base_geo, param, value)
                if not _is_geometry(gs):
                    continue
                geo = ev.make_geometry(gs, N_t, l)
                try:
                    design = ev.run_point(cfg, geo, ev.make_ops(cfg, geo, pt['T_a_i'], pt['phi'], pt['p_a']), method, 'dry')
                    year = pd.DataFrame([dict(limited=(c := operation.control_point(cfg, geo, p, method))['limited'], **c['row']) for _, p in states])
                except RuntimeError:                                  # a solver that does not converge on an extreme variant
                    continue
                masses = ev.calc_masses(geo)
                rows.append(dict(set=key, param=param, value=float(value), standard=bool(np.isclose(value, base[param])),
                                 T_c_o=design['T_c_o'], pinch=design['pinch'], wang_ok=design['wang_ok'], vdi_ok=design['vdi_ok'],
                                 E_fan_kWh=year.P_fan.sum() * stride / 1e3, E_pump_kWh=year.P_pump.sum() * stride / 1e3,
                                 water_m3=year.m_dot_ev.sum() * 3600 * stride / 1e3, h_limited=year.limited.sum() * stride,
                                 w_fr_mean=year.w_fr.mean(), m_Cu=masses['m_Cu'], m_Al=masses['m_Al']))
    return pd.DataFrame(rows)


def _prepare_tables(cfg, ranges):
    """Builds / loads the NTU table of every row count once, in this process (the workers then only read the files)."""
    if cfg.solver() == 'ntu_table':
        from src.solvers.ntu_table import get_ntu_table
        rows = range(int(ranges['N_r'][0]), int(ranges['N_r'][1]) + 1) if 'N_r' in ranges else ()
        for n_r in set(rows) | {cfg.cell(k).N_r for k in cfg.sets}:
            get_ntu_table(n_r, cfg.solver_settings.ntu_n_elements)


def run(cfg, stride=None, ranges=None) -> pd.DataFrame:
    """The variation table: one row per cooler, parameter and value (see the module docstring). ranges: like RANGES (a subset for tests)."""
    stride = max(STRIDE, cfg.year_stride) if stride is None else stride
    ranges = RANGES if ranges is None else ranges
    keys = list(cfg.sets)
    with quiet_solvers():
        _prepare_tables(cfg, ranges)
    if len(keys) == 1:
        return _run_set(cfg, keys[0], stride, ranges)
    with ProcessPoolExecutor(len(keys)) as pool:                              # one process per cooler (the module is re-imported there: ranges go as an argument)
        return pd.concat(pool.map(_run_set, repeat(cfg), keys, repeat(stride), repeat(ranges)), ignore_index=True)


def add_costs(cfg, df) -> pd.DataFrame:
    """Annual cost columns (prices of parameters.py): cost_fan, cost_pump, cost_water, cost_material (write-off), cost_total."""
    df = df.copy()
    df['cost_fan'] = df.E_fan_kWh * cfg.electricity_eur_per_kwh
    df['cost_pump'] = df.E_pump_kWh * cfg.electricity_eur_per_kwh
    df['cost_water'] = df.water_m3 * cfg.water_eur_per_m3
    df['cost_material'] = [operation.writeoff(cfg, r) for _, r in df[['m_Cu', 'm_Al']].iterrows()]
    df['cost_total'] = df.cost_fan + df.cost_pump + df.cost_water + df.cost_material
    return df
