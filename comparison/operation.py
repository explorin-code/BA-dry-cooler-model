"""
operation.py -- how the compared modules operate: in representative ambient states, with a fixed fan over the
climate year, with the face velocity swept, and with the variable-speed fan controlled hour by hour.

Every function takes the run's settings `cfg` (config.Settings) first and solves with the comparison's solver
(cfg.solver()) unless `method` is given.

  ambient_cases     the representative states (coldest hour, median, P95, ..., design), dry and precooled
  year_fixed        the climate year at ONE fixed face velocity: dry if the target is met, else precooled
  airflow_sweep     face velocity min..max in every representative state, dry and precooled
  year_modes        the climate year at ONE fixed face velocity, EVERY hour in both modes (dry and precooled)
  year_control      the climate year with the variable-speed fan: per hour the SMALLEST face velocity that
                    reaches the target (dry first, then the pad, else the maximum = target missed)
  annual_summary    energy, water, cost and hour counts per module of a (controlled) year
"""

import logging
import time

import numpy as np
import pandas as pd
from scipy.optimize import brentq

from comparison import evaluate as ev, sizing
from src.inputs import weather

logger = logging.getLogger(__name__)

N_SWEEP = 25


def _method(cfg, method):
    return cfg.solver() if method is None else method


def _warm_up(cfg, method, geo):
    """First call builds the NTU table (cached on disk) -- keep it out of the timing of the year."""
    ev.solve(cfg, method, ev.make_ops(cfg, geo, 20.0, 0.5, 95500.0), geo)


# --- representative states ------------------------------------------------------------------------------------

def ambient_cases(cfg, method=None) -> pd.DataFrame:
    """Every module in every representative ambient state, dry and precooled."""
    method = _method(cfg, method)
    rows = []
    for key in cfg.sets:
        geo = sizing.get_sized(cfg, key)
        for case, pt in cfg.design_points().items():
            ops = ev.make_ops(cfg, geo, pt['T_a_i'], pt['phi'], pt['p_a'])
            for mode in ('dry', 'precooled'):
                rows.append(dict(set=key, case=case, **ev.run_point(cfg, geo, ops, method, mode)))
    return pd.DataFrame(rows)


def airflow_sweep(cfg, method=None, n=N_SWEEP) -> pd.DataFrame:
    """Face velocity w_fr_min..w_fr_design in every representative state, dry and precooled."""
    method = _method(cfg, method)
    w_grid = np.linspace(cfg.w_fr_min, cfg.w_fr_design, n)
    rows = []
    for key in cfg.sets:
        geo = sizing.get_sized(cfg, key)
        for case, pt in cfg.design_points().items():
            for w in w_grid:
                ops = ev.make_ops(cfg, geo, pt['T_a_i'], pt['phi'], pt['p_a'], w_fr=float(w))
                for mode in ('dry', 'precooled'):
                    rows.append(dict(set=key, case=case, **ev.run_point(cfg, geo, ops, method, mode)))
    return pd.DataFrame(rows)


# --- the climate year -------------------------------------------------------------------------------------------

def year_fixed(cfg, w_fr=None, stride=1, method=None) -> pd.DataFrame:
    """Hourly operation at the FIXED face velocity w_fr (default: the design velocity); every hour dry if that
    reaches the target, else precooled. One row per module and hour (stride > 1: every n-th hour)."""
    method = _method(cfg, method)
    rows = []
    for key in cfg.sets:
        geo = sizing.get_sized(cfg, key)
        _warm_up(cfg, method, geo)
        t0 = time.perf_counter()
        for ts, pt in weather.hourly_states(weather.load_year(), stride):
            ops = ev.make_ops(cfg, geo, pt['T_a_i'], pt['phi'], pt['p_a'], w_fr=w_fr)
            rows.append(dict(set=key, hour=ts, **ev.run_auto(cfg, geo, ops, method)))
        logger.info("[Year, fixed] %s: %.1f s", key, time.perf_counter() - t0)
    return pd.DataFrame(rows)


def year_modes(cfg, w_fr=None, stride=1, method=None) -> pd.DataFrame:
    """Hourly operation at the FIXED face velocity w_fr (default: the design velocity) in BOTH modes, dry and precooled, for every hour
    (year_fixed picks one of the two): one row per module, hour and mode -- the pinch / outlet temperature against the ambient temperature."""
    method = _method(cfg, method)
    rows = []
    for key in cfg.sets:
        geo = sizing.get_sized(cfg, key)
        _warm_up(cfg, method, geo)
        t0 = time.perf_counter()
        for ts, pt in weather.hourly_states(weather.load_year(), stride):
            ops = ev.make_ops(cfg, geo, pt['T_a_i'], pt['phi'], pt['p_a'], w_fr=w_fr)
            for mode in ('dry', 'precooled'):
                rows.append(dict(set=key, hour=ts, **ev.run_point(cfg, geo, ops, method, mode)))
        logger.info("[Year, both modes] %s: %.1f s", key, time.perf_counter() - t0)
    return pd.DataFrame(rows)


def fixed_year(df) -> pd.DataFrame:
    """The fixed-speed year of the operator's rule (what year_fixed solves): of the year_modes rows `df`, dry in the hours in which dry meets the
    target, otherwise precooled. One row per module and hour."""
    dry, wet = df[df['mode'] == 'dry'], df[df['mode'] == 'precooled']
    missed = dry[~dry.target_met][['set', 'hour']]
    return pd.concat([dry[dry.target_met], wet.merge(missed, on=['set', 'hour'])]).sort_values(['set', 'hour']).reset_index(drop=True)


def target_crossing(df, mode, target, bin_width=2.0):
    """Ambient temperature [degC] at which the median water outlet of `mode` in df (year_modes rows of ONE module) reaches `target`
    (interpolated over bins of bin_width; NaN if it does not within the year)."""
    d = df[df['mode'] == mode]
    g = d.groupby((d.T_a_amb // bin_width).astype(int).rename('bin'))[['T_a_amb', 'T_c_o']].median().sort_values('T_a_amb')
    below, above = g[g.T_c_o <= target], g[g.T_c_o > target]
    if below.empty or above.empty:
        return float('nan')
    lo, hi = below.iloc[-1], above.iloc[0]
    return float(lo.T_a_amb + (target - lo.T_c_o) / (hi.T_c_o - lo.T_c_o) * (hi.T_a_amb - lo.T_a_amb))


def _min_flow(cfg, geo, pt, mode, method):
    """(w_fr, row) with the smallest face velocity in [w_fr_min, w_fr_design] that reaches the target in `mode`;
    None if even w_fr_design misses."""
    def at(w):
        ops = ev.make_ops(cfg, geo, pt['T_a_i'], pt['phi'], pt['p_a'], w_fr=w)
        return ev.run_point(cfg, geo, ops, method, mode)

    low = at(cfg.w_fr_min)
    if low['target_met']:
        return cfg.w_fr_min, low                       # the slowest fan setting already overcools
    high = at(cfg.w_fr_design)
    if not high['target_met']:
        return None
    if high['T_c_o'] >= cfg.T_c_target:                 # met only within the tolerance: full speed is the smallest sufficient setting
        return cfg.w_fr_design, high
    w = brentq(lambda w: at(w)['T_c_o'] - cfg.T_c_target, cfg.w_fr_min, cfg.w_fr_design, xtol=1e-3)
    return w, at(w)


def control_point(cfg, geo, pt, method=None, allow_pad=True) -> dict:
    """The variable-speed fan in one ambient state pt: dry at the smallest sufficient face velocity, else the
    pad likewise, else the pad at full speed (limited = target missed). allow_pad=False: no pad at all -- dry at the smallest
    sufficient velocity, else dry at full speed (limited)."""
    method = _method(cfg, method)
    modes = ('dry', 'precooled') if allow_pad else ('dry',)
    for mode in modes:
        found = _min_flow(cfg, geo, pt, mode, method)
        if found:
            return dict(row=found[1], limited=False)
    ops = ev.make_ops(cfg, geo, pt['T_a_i'], pt['phi'], pt['p_a'], w_fr=cfg.w_fr_design)
    return dict(row=ev.run_point(cfg, geo, ops, method, 'precooled' if allow_pad else 'dry'), limited=True)


def year_control(cfg, stride=1, method=None, allow_pad=True) -> pd.DataFrame:
    """The climate year with the variable-speed fan, one row per module and hour (see control_point; allow_pad=False: the dry cooler without a pad)."""
    method = _method(cfg, method)
    rows = []
    for key in cfg.sets:
        geo = sizing.get_sized(cfg, key)
        _warm_up(cfg, method, geo)
        t0 = time.perf_counter()
        for ts, pt in weather.hourly_states(weather.load_year(), stride):
            c = control_point(cfg, geo, pt, method, allow_pad)
            rows.append(dict(set=key, hour=ts, limited=c['limited'], **c['row']))
        logger.info("[Year, control%s] %s: %.1f s", "" if allow_pad else ", dry", key, time.perf_counter() - t0)
    return pd.DataFrame(rows)


# --- summaries --------------------------------------------------------------------------------------------------

def writeoff(cfg, masses) -> float:
    """Yearly write-off of the copper and aluminium of a module [EUR/a]: raw-metal masses (ev.calc_masses) x prices,
    straight-line over cfg.lifetime_years."""
    return (masses['m_Cu'] * cfg.copper_eur_per_kg + masses['m_Al'] * cfg.aluminium_eur_per_kg) / cfg.lifetime_years


def material_writeoff(cfg, key) -> float:
    """writeoff of the sized module `key`."""
    return writeoff(cfg, ev.calc_masses(sizing.get_sized(cfg, key)))


def add_costs(cfg, df):
    """Annual cost columns of a table indexed by set: electricity (fan + pump), water, material write-off
    (prices: parameters.py); cost_total is their sum."""
    df = df.copy()
    df['cost_el'] = (df.E_fan_kWh + df.E_pump_kWh) * cfg.electricity_eur_per_kwh
    df['cost_water'] = df.water_m3 * cfg.water_eur_per_m3
    df['cost_material'] = [material_writeoff(cfg, k) for k in df.index]
    df['cost_total'] = df.cost_el + df.cost_water + df.cost_material
    return df


def annual_summary(cfg, df, hours_per_row=1.0) -> pd.DataFrame:
    """Per module: fan/pump energy [kWh], water [m3], cost [EUR], hours precooled / at the limit / meeting
    the target / at freezing risk. hours_per_row = stride of the year the rows come from."""
    g = df.groupby('set')
    return pd.DataFrame(dict(
        E_fan_kWh=g['P_fan'].sum() * hours_per_row / 1e3,
        E_pump_kWh=g['P_pump'].sum() * hours_per_row / 1e3,
        water_m3=g['m_dot_ev'].sum() * 3600 * hours_per_row / 1e3,
        h_precooled=g['mode'].apply(lambda m: (m == 'precooled').sum() * hours_per_row),
        h_limited=g['limited'].sum() * hours_per_row if 'limited' in df else np.nan,
        h_freeze=g['freeze'].sum() * hours_per_row,
        h_target_met=g['target_met'].sum() * hours_per_row,
        w_fr_mean=g['w_fr'].mean())).pipe(lambda d: add_costs(cfg, d))
