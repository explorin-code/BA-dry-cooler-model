"""
wall.py -- tube wall thickness delta_p (outer diameter fixed, d_i = d - 2 delta_p).
Two ways to compare:
  'N_t fixed' : the base tube count -> a thicker wall makes the bore smaller, the
                water faster (higher Re, alpha_c), pressure drop up, conduction
                resistance up, copper mass up.
  'Re fixed'  : N_t is rescaled to keep the tube Reynolds number -> a thicker wall
                needs FEWER tubes (smaller bore at the same Re carries less water each,
                so MORE tubes -- see output), face area and air flow change with it.
The wall was an assumption for P1, P3 (0.5 mm) and P4 (2 mm).
"""

import dataclasses
import logging

import numpy as np
import pandas as pd

from comparison import config, evaluate as ev, sizing

logger = logging.getLogger(__name__)

WALLS = np.array([0.3, 0.5, 0.75, 1.0, 1.5, 2.0])


def run(cfg, method=None) -> pd.DataFrame:
    method = cfg.solver(table_ok=False) if method is None else method
    pt = cfg.design_point()
    rows = []
    for key, gs in cfg.sets.items():
        base_N_t = sizing.get_sized(cfg, key).N_t
        for wall in sorted(set(WALLS) | {gs.wall_mm}):
            gv = dataclasses.replace(gs, wall_mm=float(wall))
            if 2 * wall >= gv.d_mm - 1.0:
                continue
            for basis in ('N_t fixed', 'Re fixed'):
                N_t = base_N_t if basis == 'N_t fixed' else sizing.tubes_for_flow(cfg, gv, pt)
                geo = ev.make_geometry(gv, N_t, cfg.l_tube)
                ops = ev.make_ops(cfg, geo, pt['T_a_i'], pt['phi'], pt['p_a'])
                rows.append(dict(set=key, case='design', basis=basis, wall_mm=float(wall), wall_base=gs.wall_mm,
                                 **ev.geometry_summary(geo), **ev.run_point(cfg, geo, ops, method, 'dry')))
        logger.info("[Wall] %s done", key)
    df = pd.DataFrame(rows)
    df['kg_per_kW'] = df.m_total / df.Q_kW
    return df


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    pd.set_option('display.width', 250)
    d = run(config.settings())
    print(d[['set', 'basis', 'wall_mm', 'N_t', 'Re_c', 'Q_kW', 'T_c_o', 'P_pump', 'm_Cu', 'kg_per_kW']].round(2).to_string())
