"""
sensitivity.py -- one-at-a-time sensitivities around a base module
(a sized geometry set) over the literature ranges of config.RANGES.

Constraints while one parameter moves (stated in the thesis):
  * s_t: the module width stays (N_t = round(b / s_t)), so the face area and
    the air mass flow stay; the water velocity follows N_t.
  * d, s_f, delta_f, s_l, N_r: N_t and l stay.
  * w_fr: the air flow (operating parameter).
Two bases of comparison: 'w_fr' = face velocity held at the design value,
'P_fan' = face velocity adjusted so the fan power equals the base module's
(fair for geometries that change the air-side pressure drop).
"""

import dataclasses
import logging

import numpy as np
import pandas as pd
from scipy.optimize import brentq

from comparison import config, evaluate as ev, sizing

logger = logging.getLogger(__name__)

N_POINTS = 9
PARAMS = ['wall_mm', 'd_mm', 's_f_mm', 'delta_f_mm', 's_t_mm', 's_l_mm', 'N_r', 'w_fr']


def _variant(base_set, base_geo, param, value):
    """(GeometrySet, Geometry) with `param` set to `value`."""
    if param == 'w_fr':
        return base_set, base_geo
    kw = {param: int(round(value)) if param == 'N_r' else float(value)}
    gs = dataclasses.replace(base_set, **kw)
    N_t = base_geo.N_t
    if param == 's_t_mm':
        N_t = max(1, round(base_geo.b / (gs.s_t_mm / 1e3)))
    return gs, ev.make_geometry(gs, N_t, base_geo.l)


def _row(cfg, geo, pt, w_fr, basis, param, value):
    ops = ev.make_ops(cfg, geo, pt['T_a_i'], pt['phi'], pt['p_a'], w_fr=w_fr)
    return dict(param=param, value=value, basis=basis, **ev.geometry_summary(geo),
                **ev.run_point(cfg, geo, ops, cfg.solver(table_ok=False), 'dry'))


def run(cfg, set_key=None, n=None):
    """One-at-a-time sensitivities of the base module (`set_key`, default the configured sensitivity cooler) over config.RANGES with `n` values per parameter. Every variant is solved at the design face velocity and at the face velocity that gives the base fan power: one row per variant and basis."""
    n = N_POINTS if n is None else n
    set_key = cfg.sensitivity_set if set_key is None else set_key
    pt = cfg.design_point()
    base_set = cfg.cell(set_key)
    base_geo = sizing.get_sized(cfg, set_key)
    P_fan_base = _row(cfg, base_geo, pt, cfg.w_fr_design, 'w_fr', 'base', np.nan)['P_fan']
    rows = [dict(_row(cfg, base_geo, pt, cfg.w_fr_design, 'w_fr', 'base', np.nan), set=set_key)]
    base_values = dict(wall_mm=base_set.wall_mm, d_mm=base_set.d_mm, s_f_mm=base_set.s_f_mm, delta_f_mm=base_set.delta_f_mm,
                       s_t_mm=base_set.s_t_mm, s_l_mm=base_set.s_l_mm, N_r=base_set.N_r, w_fr=cfg.w_fr_design)

    for param in PARAMS:
        lo, hi, unit, _ = config.RANGES[param]
        values = np.unique(np.round(np.linspace(lo, hi, n), 3)) if param != 'N_r' else np.arange(lo, hi + 1)
        values = np.unique(np.append(values, base_values[param]))
        for v in values:
            gs, geo = _variant(base_set, base_geo, param, v)
            w_fr_fixed = float(v) if param == 'w_fr' else cfg.w_fr_design
            try:
                rows.append(dict(_row(cfg, geo, pt, w_fr_fixed, 'w_fr', param, v), set=set_key))
            except RuntimeError as err:               # a solver that does not converge on an extreme variant
                logger.warning("[Sens] %s = %s skipped: %s", param, v, err)
                continue
            if param != 'w_fr':           # same fan power as the base module
                try:
                    f = lambda w: _row(cfg, geo, pt, w, 'P_fan', param, v)['P_fan'] - P_fan_base
                    w = brentq(f, 0.05, 12.0, xtol=1e-3)
                    rows.append(dict(_row(cfg, geo, pt, w, 'P_fan', param, v), set=set_key))
                except (ValueError, RuntimeError):
                    logger.warning("[Sens] %s = %s: no face velocity gives the base fan power", param, v)
        logger.info("[Sens] %s done", param)
    return pd.DataFrame(rows)


if __name__ == '__main__':
    import time
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    t = time.time()
    df = run(config.settings())
    print(len(df), "rows in", round(time.time() - t), "s")
    pd.set_option('display.width', 250)
    print(df[df.basis == 'w_fr'][['param', 'value', 'Q_kW', 'T_c_o', 'P_fan', 'Delta_p_a', 'w_c', 'wang_ok', 'vdi_ok']].round(3).to_string())
