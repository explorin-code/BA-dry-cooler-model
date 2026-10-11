"""
tubes.py -- what happens if the number of tubes per row N_t is not the one that gives
the design tube Reynolds number? N_t is varied at fixed water flow (1 kg/s), tube length
and rows, so the water velocity / Reynolds number changes with it AND the face area
(hence the air mass flow at fixed face velocity) and the material mass grow with N_t.
Design point (dry) and P95 (precooled) per elementary cell.
"""

import logging

import numpy as np
import pandas as pd

from comparison import config, evaluate as ev, sizing

logger = logging.getLogger(__name__)

FACTORS = np.array([0.3, 0.4, 0.5, 0.65, 0.8, 1.0, 1.25, 1.5, 2.0, 2.5])      # N_t / N_t(Re_design)


def run(cfg, method=None) -> pd.DataFrame:
    method = cfg.solver(table_ok=False) if method is None else method
    points = cfg.design_points()
    cases = [('design', 'dry'), ('P95', 'precooled')]
    rows = []
    for key, gs in cfg.sets.items():
        base = sizing.get_sized(cfg, key)
        for N_t in sorted({max(1, round(base.N_t * f)) for f in FACTORS}):
            geo = ev.make_geometry(gs, N_t, base.l)
            for case, mode in cases:
                pt = points[case]
                ops = ev.make_ops(cfg, geo, pt['T_a_i'], pt['phi'], pt['p_a'])
                rows.append(dict(set=key, case=case, N_t_base=base.N_t, **ev.geometry_summary(geo),
                                 **ev.run_point(cfg, geo, ops, method, mode)))
        logger.info("[Tubes] %s done", key)
    df = pd.DataFrame(rows)
    df['Q_per_P_fan'] = df.Q_kW * 1e3 / df.P_fan
    df['kg_per_kW'] = df.m_total / df.Q_kW
    return df


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    pd.set_option('display.width', 250)
    d = run(config.settings())
    print(d[d.case == 'design'][['set', 'N_t', 'Re_c', 'w_c', 'A_fr', 'Q_kW', 'T_c_o', 'P_fan', 'P_pump', 'm_total', 'kg_per_kW']].round(2).to_string())
