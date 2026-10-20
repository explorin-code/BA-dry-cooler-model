"""
analyses.py -- the catalogue of cooler-comparison analyses.

Every analysis is one `Item` (src/catalogue.py): `compute(cfg)` -> {name: DataFrame}, `plot(data, cfg)` -> figures,
`log(data, cfg)` -> a few summary lines for the terminal / GUI log (cfg = config.Settings of the run). runner.py runs any selection of them the same way
for the command line and the GUI; the default selection is parameters.COMPARE_ANALYSES.

To add an analysis: write its compute/plot/log and append an Item to ANALYSES.
"""

import pandas as pd

from comparison import evaluate as ev, operation, plots_design, plots_operation, plots_print, rows, sensitivity, sizing, tubes, variation, wall
from src.catalogue import Item
from src.plots.plot_style import in_large_format, large_format_context


# --- geometry ---------------------------------------------------------------------------------------------------

def _compute_geometry(cfg):
    pt = cfg.design_point()
    rows = []
    for key, cell in cfg.sets.items():
        geo = sizing.get_sized(cfg, key)
        ops = ev.make_ops(cfg, geo, pt['T_a_i'], pt['phi'], pt['p_a'])
        rows.append(dict(set=key, label=cell.label, **ev.geometry_summary(geo), **ev.run_point(cfg, geo, ops, cfg.solver(), 'dry')))
    return {'modules': pd.DataFrame(rows)}


def _plot_geometry(data, cfg):
    """The print-size drawing (plot_geometry.plot_print) and the large-format overview figures."""
    from comparison import plot_geometry
    plot_geometry.plot_print(cfg)
    with large_format_context():
        plot_geometry.plot(cfg)
        plot_geometry.plot_bank(cfg)
        plots_design.geometries(data['modules'], cfg)


def _log_geometry(data, cfg):
    return [f"{r['set']}: N_t = {r.N_t}, N_r = {r.N_r}, l = {r.l:.2f} m, {r.width:.2f} m wide, {r.depth * 1e3:.0f} mm deep | "
            f"Q = {r.Q_kW:.1f} kW, T_c,o = {r.T_c_o:.1f} °C | P_fan = {r.P_fan:.0f} W | m_Al = {r.m_Al:.0f} kg, m_Cu = {r.m_Cu:.0f} kg"
            for _, r in data['modules'].iterrows()]


# --- operation --------------------------------------------------------------------------------------------------

def _log_year(df, w_fr):
    return [f"{key} at {w_fr:.1f} m/s: target met in {int(d.target_met.sum())} of {len(d)} h "
            f"({int((d['mode'] == 'precooled').sum())} h with the pad), fan {d.P_fan.mean():.0f} W mean, "
            f"{int(d.freeze.sum())} h with freezing risk" for key, d in df.groupby('set')]


def _compute_fan_control(cfg):
    return {'control': operation.year_control(cfg, cfg.year_stride),
            'fixed': operation.year_fixed(cfg, w_fr=cfg.w_fr_design, stride=cfg.year_stride)}


def _log_control(data, cfg):
    s = operation.annual_summary(cfg, data['control'], hours_per_row=float(cfg.year_stride))
    return _log_year(data['fixed'], cfg.w_fr_design) + [f"{key}: fan {r.E_fan_kWh:.0f} kWh/a, pump {r.E_pump_kWh:.0f} kWh/a, water {r.water_m3:.0f} m3/a "
            f"| cost {r.cost_total:.0f} EUR/a (electricity {r.cost_el:.0f}, water {r.cost_water:.0f}, material {r.cost_material:.0f}) "
            f"| {r.h_precooled:.0f} h precooled, {r.h_limited:.0f} h at the top fan speed without reaching the target"
            for key, r in s.iterrows()]


def _plot_fan_control(data, cfg):
    plots_operation.control_year(data['control'], cfg)
    plots_operation.fan_year(data['control'], cfg)
    plots_operation.annual_energy(data['fixed'], data['control'], cfg)
    plots_operation.cost_year(data['control'], cfg)
    plots_operation.year_scatter(data['fixed'], cfg)
    if len(data['control']) >= 8000 * data['control']['set'].nunique():      # the fan curves need the full year
        from src.inputs import weather
        plots_operation.fan_curves(data['control'], weather.load_year(), cfg)


def _compute_pinch_fan(cfg):
    return {'year': operation.year_modes(cfg, cfg.w_fr_design, cfg.year_stride),
            'sweep': operation.airflow_sweep(cfg),
            'control': operation.year_control(cfg, cfg.year_stride),
            'control_dry': operation.year_control(cfg, cfg.year_stride, allow_pad=False)}


def _plot_pinch_fan(data, cfg):
    plots_print.pinch_year(data['year'], cfg)
    plots_print.missed_hours(data['year'], cfg)
    plots_print.pinch_year_clouds(data['year'], cfg)
    plots_print.annual_costs(data['year'], data['control'], data['control_dry'], cfg)
    plots_print.monthly_cost(data['control'], cfg)
    plots_print.water_vs_hours(data['year'], data['control'], data['control_dry'], cfg)
    plots_print.pinch_speed(data['sweep'], cfg)
    plots_print.controlled_fan(data['control'], cfg)
    plots_print.fan_cop(data['control'], cfg)


def _log_pinch_fan(data, cfg):
    lines = []
    for key in cfg.sets:
        d = data['year'][data['year'].set == key]
        limited = int(data['control'][data['control'].set == key].limited.sum() * cfg.year_stride)
        lines.append(f"{key} at {cfg.w_fr_design:g} m/s: the median outlet reaches {cfg.T_c_target:g} °C at T_amb = "
                     f"{operation.target_crossing(d, 'dry', cfg.T_c_target):.1f} °C dry, {operation.target_crossing(d, 'precooled', cfg.T_c_target):.1f} °C precooled; "
                     f"pinch at the warmest hour {d[d['mode'] == 'dry'].pinch.min():.1f} K dry / {d[d['mode'] == 'precooled'].pinch.min():.1f} K precooled; "
                     f"controlled fan misses the target in {limited} h")
    return lines


# --- design analyses ---------------------------------------------------------------------------------------------

def _log_variation(data, cfg):
    """Per parameter: the cheapest variant of the sensitivity cooler that misses the target for at most plots_print.MISSED_H hours, against its standard value."""
    df = variation.add_costs(cfg, data['variation'])
    df = df[df.set == cfg.sensitivity_set]
    lines = [f"{cfg.sensitivity_set}: standard cost {df[df.standard].cost_total.mean():.0f} EUR/a; cheapest variant per parameter with at most "
             f"{plots_print.MISSED_H} h above {cfg.T_c_target:g} °C:"]
    for param, d in df.groupby('param', sort=False):
        best = d[d.h_limited <= plots_print.MISSED_H].nsmallest(1, 'cost_total').iloc[0]
        lines.append(f"  {param:<10} = {best.value:g}: {best.cost_total:.0f} EUR/a, pinch {best.pinch:.1f} K, {best.h_limited:.0f} h missed")
    return lines


def _log_sensitivity(data, cfg):
    df = data['sensitivity']
    base = df[df.param == 'base'].iloc[0]
    lines = [f"base module {cfg.sensitivity_set}: Q = {base.Q_kW:.1f} kW, T_c,o = {base.T_c_o:.1f} °C; "
             f"change of T_c,o over the literature range:"]
    for p, d in df[(df.param != 'base') & (df.basis == 'w_fr')].groupby('param'):
        d = d.sort_values('value')
        lines.append(f"  {p:<10} {d.value.iloc[0]:5.1f} .. {d.value.iloc[-1]:5.1f} -> T_c,o "
                     f"{d.T_c_o.iloc[0] - base.T_c_o:+.1f} K .. {d.T_c_o.iloc[-1] - base.T_c_o:+.1f} K")
    return lines


def _log_reach(df, name, cfg):
    """First row/tube/wall value at which a set reaches the target at the design point."""
    lines = []
    for key, d in df[df.case == 'design'].groupby('set'):
        ok = d[d.T_c_o <= cfg.T_c_target]
        lines.append(f"{key}: " + (f"target reached from {name} = {ok[name].min():g}" if len(ok) else f"target not reached up to {name} = {d[name].max():g}"))
    return lines


ANALYSES = {a.key: a for a in (
    Item('geometry', "Geometry overview", _compute_geometry,
         _plot_geometry, _log_geometry),
    Item('ambient_cases', "Ambient states",
         lambda cfg: {'cases': operation.ambient_cases(cfg)},
         in_large_format(lambda d, cfg: plots_operation.ambient_cases(d['cases'], cfg))),
    Item('airflow_sweep', "Face-velocity sweep",
         lambda cfg: {'sweep': operation.airflow_sweep(cfg)},
         in_large_format(lambda d, cfg: plots_operation.airflow_sweep(d['sweep'], cfg))),
    Item('fan_control', "Year with fan control",
         _compute_fan_control, in_large_format(_plot_fan_control), _log_control),
    Item('pinch_fan', "Pinch, fan and annual economics (print size)", _compute_pinch_fan, _plot_pinch_fan, _log_pinch_fan),
    Item('variation', "Parameter variation: pinch and annual cost (print size)", lambda cfg: {'variation': variation.run(cfg)},
         lambda d, cfg: plots_print.variation_figures(d['variation'], cfg), _log_variation),
    Item('sensitivity', "Sensitivity", lambda cfg: {'sensitivity': sensitivity.run(cfg)},
         in_large_format(lambda d, cfg: plots_design.sensitivity(d['sensitivity'], cfg)), _log_sensitivity),
    Item('rows', "Rows", lambda cfg: {'rows': rows.run(cfg)}, in_large_format(lambda d, cfg: plots_design.rows_analysis(d['rows'], cfg)),
         lambda d, cfg: _log_reach(d['rows'], 'N_r', cfg)),
    Item('tubes', "Tubes per row", lambda cfg: {'tubes': tubes.run(cfg)},
         in_large_format(lambda d, cfg: plots_design.tubes_analysis(d['tubes'], cfg))),
    Item('wall', "Wall thickness", lambda cfg: {'wall': wall.run(cfg)}, in_large_format(lambda d, cfg: plots_design.wall_analysis(d['wall'], cfg))),
)}
