"""
rows.py -- how many tube rows N_r? Sweep N_r = 1..N_R_MAX for each elementary cell
(tube length and water flow fixed, N_t from the tube Reynolds number) and
report, at the design point (dry, 20 °C / 60 %) and at the P95 summer point
(precooled), duty, outlet temperature, powers, depth and material masses.
Criteria for the choice (thesis text): reaches the 25 °C target at the design
point, fan power and material grow ~linearly with N_r while the duty gains
flatten, Wang's range (N_r <= 6).
"""

import logging

import pandas as pd

from comparison import config, evaluate as ev, sizing

logger = logging.getLogger(__name__)

N_R_MAX = 6                   # upper limit chosen for the thesis (Wang's tested range; Kroger 3-6 typical)


def run(cfg, method=None) -> pd.DataFrame:
    method = cfg.solver(table_ok=False) if method is None else method
    points = cfg.design_points()
    cases = [('design', 'dry'), ('P95', 'precooled')]
    rows = []
    for key in cfg.sets:
        for n in range(1, N_R_MAX + 1):
            geo = sizing.get_sized(cfg, key, N_r=n)
            for case, mode in cases:
                pt = points[case]
                ops = ev.make_ops(cfg, geo, pt['T_a_i'], pt['phi'], pt['p_a'])
                rows.append(dict(set=key, case=case, **ev.geometry_summary(geo),
                                 **ev.run_point(cfg, geo, ops, method, mode)))
        logger.info("[Rows] %s done", key)
    df = pd.DataFrame(rows)
    df['Q_per_P_fan'] = df.Q_kW * 1e3 / df.P_fan
    df['kg_per_kW'] = df.m_total / df.Q_kW
    return df


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    pd.set_option('display.width', 250)
    df = run(config.settings())
    d = df[df.case == 'design']
    print(d[['set', 'N_r', 'N_t', 'Q_kW', 'T_c_o', 'pinch', 'P_fan', 'P_pump', 'Delta_p_a', 'depth', 'm_Al', 'm_Cu', 'kg_per_kW', 'wang_ok']].round(2).to_string())
