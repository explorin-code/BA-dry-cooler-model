"""
plots_operation.py -- figures of the operation analyses: ambient states, face velocity, the climate year with a fixed
and with the controlled fan, fan curves, energy, water and cost.
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from comparison.plots_common import LABEL, SET_COLOR, mode_color, save_figure, target_line, year_figure
from src.plot_style import TUM_ACCENT_BLUE, TUM_GREY_50, TUM_IVORY, lighten_color


def ambient_cases(df, cfg):
    cell = df
    order = ['coldest', 'median', 'design', 'wettest_warm', 'P95', 'P99', 'hottest']
    fig, axes = plt.subplots(1, 3, figsize=(24, 7))
    for ax, col in zip(axes, ('T_c_o', 'pinch', 'Q_kW')):
        for j, (mode, hatch) in enumerate((('dry', ''), ('precooled', '//'))):
            d = cell[cell['mode'] == mode]
            for i, key in enumerate(cfg.sets):
                v = [d[(d.case == c) & (d.set == key)][col].iloc[0] for c in order]
                n = len(cfg.sets)
                x = np.arange(len(order)) + (i - (n - 1) / 2) * 0.1 + (j - 0.5) * 0.45
                ax.bar(x, v, width=0.1, color=mode_color(key, mode), edgecolor='white', label=f'{key} {mode}')
        ax.set_xticks(range(len(order))); ax.set_xticklabels([o.replace('_', ' ') for o in order], rotation=30)
        ax.set_ylabel(LABEL[col])
        if col == 'T_c_o':
            target_line(ax, cfg)
    axes[0].legend(fontsize=9, ncol=2)
    fig.tight_layout()
    save_figure(fig, 'ambient_cases')


def year_scatter(df, cfg, w_fr=None):
    """T_c,o against ambient temperature over the year at fixed face velocity (dry dark, precooled bright)."""
    fig, axes = plt.subplots(1, len(cfg.sets), figsize=(8 * len(cfg.sets), 7), sharey='row', squeeze=False)
    axes = axes[0]
    for ax, key in zip(axes, cfg.sets):
        d = df[df.set == key]
        for mode in ('dry', 'precooled'):
            m = d[d['mode'] == mode]
            ax.scatter(m.T_a_amb, m.T_c_o, s=5, color=mode_color(key, mode), alpha=0.6, label=mode)
        target_line(ax, cfg)
        ax.set_title(key); ax.set_xlabel(LABEL['T_a_amb'])
    axes[0].set_ylabel(LABEL['T_c_o']); axes[0].legend(markerscale=3)
    v = cfg.w_fr_design if w_fr is None else w_fr
    fig.suptitle(rf'Year at fixed $w_{{fr}}$ = {v:g} m/s')
    fig.tight_layout()
    save_figure(fig, 'year_scatter')


def airflow_sweep(df, cfg):
    cases = ['design', 'P95', 'hottest', 'wettest_warm']
    fig, axes = plt.subplots(2, len(cases), figsize=(26, 11), sharex='row')
    for j, case in enumerate(cases):
        for key in cfg.sets:
            for mode, ls in (('dry', '-'), ('precooled', '--')):
                d = df[(df.set == key) & (df.case == case) & (df['mode'] == mode)].sort_values('w_fr')
                axes[0, j].plot(d.w_fr, d.T_c_o, ls, color=mode_color(key, mode), label=f'{key} {mode}')
                axes[1, j].plot(d.P_fan, d.Q_kW, ls, color=mode_color(key, mode))
        target_line(axes[0, j], cfg); axes[0, j].set_title(case.replace('_', ' '))
        axes[0, j].set_xlabel(LABEL['w_fr']); axes[1, j].set_xlabel(LABEL['P_fan'])
        axes[1, j].set_xscale('log')
    axes[0, 0].set_ylabel(LABEL['T_c_o']); axes[1, 0].set_ylabel(LABEL['Q_kW'])
    axes[0, 0].legend(fontsize=9)
    fig.tight_layout()
    save_figure(fig, 'airflow_sweep')


def control_year(df, cfg):
    """Required face velocity / fan power of the controlled fan vs ambient (dry dark, precooled bright)."""
    fig, axes = plt.subplots(1, 2, figsize=(20, 7))
    for key in cfg.sets:
        for mode in ('dry', 'precooled'):
            d = df[(df.set == key) & (df['mode'] == mode)]
            axes[0].scatter(d.T_a_amb, d.w_fr, s=4, color=mode_color(key, mode), alpha=0.6, label=f'{key} {mode}')
            axes[1].scatter(d.T_a_amb, d.P_fan, s=4, color=mode_color(key, mode), alpha=0.6)
    axes[0].set_ylabel(r'required $w_{fr}$ [m/s]'); axes[1].set_ylabel(LABEL['P_fan']); axes[1].set_yscale('log')
    for ax in axes:
        ax.set_xlabel(LABEL['T_a_amb'])
    axes[0].legend(markerscale=4, ncol=2, fontsize=11)
    fig.suptitle(rf'Fan control (bright: with pad), target $T_{{c,o}}$ = {cfg.T_c_target:.0f} °C')
    fig.tight_layout()
    save_figure(fig, 'control_year')


def annual_energy(year, control, cfg):
    """Annual operation with a fixed and with a controlled fan: fan energy, water, and the yearly cost
    (electricity + water + write-off of the copper and aluminium over the service life; prices from config)."""
    from comparison import operation
    keys = list(cfg.sets)
    ctrl = operation.annual_summary(cfg, control)
    g = year.groupby('set')
    fixed = pd.DataFrame(dict(E_fan_kWh=g.P_fan.sum() / 1e3, E_pump_kWh=g.P_pump.sum() / 1e3,
                              water_m3=g.m_dot_ev.sum() * 3.6)).pipe(lambda d: operation.add_costs(cfg, d))
    x = np.arange(len(keys))
    cols = [SET_COLOR[k] for k in keys]
    fig, axes = plt.subplots(1, 3, figsize=(22, 7))

    def pair(ax, f, c, ylabel):
        ax.bar(x - 0.2, f[keys], 0.4, color=TUM_GREY_50, label=rf'fixed $w_{{fr}}$ = {cfg.w_fr_design:g} m/s')
        ax.bar(x + 0.2, c[keys], 0.4, color=cols, label='variable speed')
        ax.set_ylabel(ylabel); ax.set_xticks(x); ax.set_xticklabels(keys)

    pair(axes[0], fixed.E_fan_kWh, ctrl.E_fan_kWh, r'fan energy $\sum P_{fan}\,\Delta t$ [kWh/year]')
    axes[0].legend(loc='upper center', bbox_to_anchor=(0.5, -0.08), ncol=2)
    pair(axes[1], fixed.water_m3, ctrl.water_m3, r'water use for the pad [m³/year]')
    # stacked cost: electricity, water, material write-off
    ax = axes[2]
    for off, df, alpha in ((-0.2, fixed, 0.55), (0.2, ctrl, 1.0)):
        el, wa, ma = (df.loc[keys, c] for c in ('cost_el', 'cost_water', 'cost_material'))
        ax.bar(x + off, el, 0.4, color=[TUM_GREY_50 if off < 0 else c for c in cols], alpha=alpha)
        ax.bar(x + off, wa, 0.4, bottom=el, color=TUM_ACCENT_BLUE, alpha=alpha, label='water' if off > 0 else None)
        ax.bar(x + off, ma, 0.4, bottom=el + wa, color=TUM_IVORY, edgecolor=TUM_GREY_50, alpha=alpha,
               label=f'copper + aluminium, {cfg.lifetime_years:g} years' if off > 0 else None)
        for xi, k in zip(x, keys):
            ax.text(xi + off, df.loc[k, 'cost_total'] + 10, f"{df.loc[k, 'cost_total']:.0f}", ha='center', fontsize=11)
    ax.set_ylabel('annual cost [€/year]'); ax.set_xticks(x); ax.set_xticklabels(keys)
    ax.set_title('electricity + water + material')
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.08), ncol=2)
    fig.suptitle('Annual operation: fixed air flow (grey) vs. fan control')
    fig.tight_layout()
    save_figure(fig, 'annual_energy')


def fan_curves(control, year, cfg):
    """Fan curves w_fr(T_a,i), P_fan(T_a,i) per set (dry dark, precooled bright), duration curve of the
    fan power, and the check of the curve-based schedule against the solved control."""
    from comparison import fan_curve
    from src import weather
    year = weather.add_wet_bulb(year)
    curves = {k: fan_curve.FanCurve(cfg, control, k) for k in cfg.sets}
    sched = {k: fan_curve.schedule(curves[k], year) for k in cfg.sets}

    fig, axes = plt.subplots(2, 2, figsize=(18, 13))
    T = np.linspace(-11, 25, 400)
    for col, (quantity, ylabel) in enumerate((('w_fr', LABEL['w_fr']), ('P_fan', LABEL['P_fan']))):
        ax = axes[0, col]
        for k, c in curves.items():
            for mode, br in (('dry', c.dry), ('precooled', c.pre)):
                d = control[(control.set == k) & (control['mode'] == mode)]
                ax.scatter(d.T_a_i, d[quantity], s=2, color=mode_color(k, mode), alpha=0.25)
                Tb = T[(T >= br.T_min) & (T <= br.T_max)]
                label = f'{k} {mode}' + (f' (to {c.T_dry_max:.1f} °C)' if mode == 'dry' else '')
                ax.plot(Tb, getattr(br, quantity)(Tb), color=mode_color(k, mode), lw=2.4, label=label)
        ax.set_xlabel(r'cooler inlet air $T_{a,i}$ [°C]'); ax.set_ylabel(ylabel)
        if quantity == 'w_fr':
            ax.axhline(cfg.w_fr_design, color='k', ls=':', lw=1); ax.axhline(cfg.w_fr_min, color='k', ls=':', lw=1)
            ax.legend(fontsize=9, ncol=2)
            ax.set_title('required face velocity')
        else:
            ax.set_yscale('log'); ax.set_title('fan power')
    ax = axes[1, 0]
    for k, s_ in sched.items():
        ax.plot(np.arange(len(s_)) / len(s_) * 100, np.sort(s_.P_fan.values)[::-1], color=SET_COLOR[k], lw=2, label=k)
    ax.set_yscale('log'); ax.set_xlabel('share of the year [%]'); ax.set_ylabel(LABEL['P_fan'])
    ax.set_title('fan power, duration curve'); ax.legend()
    ax = axes[1, 1]
    rows = []
    for k in cfg.sets:
        v = fan_curve.validate(curves[k], control, k, year)
        rows.append((k, v['E_curve_kWh'], v['E_control_kWh'], v['w_err_mean']))
    x = np.arange(len(rows))
    ax.bar(x - 0.2, [r[1] for r in rows], 0.4, color=[SET_COLOR[r[0]] for r in rows], label='from the fan curve')
    ax.bar(x + 0.2, [r[2] for r in rows], 0.4, color='0.6', label='solved hour by hour')
    ax.set_xticks(x); ax.set_xticklabels([r[0] for r in rows]); ax.set_ylabel('fan energy [kWh/year]'); ax.legend()
    ax.set_title('annual fan energy: curve vs. solved')
    fig.suptitle('Fan curves')
    fig.tight_layout()
    save_figure(fig, 'fan_curves')


def fan_year(control, cfg):
    """The year of the controlled fan: ambient temperature, face velocity and fan power per set,
    daily mean (left) and daily maximum (right) in separate columns."""
    c = control.copy()
    c['hour'] = pd.to_datetime(c['hour'])
    first = c[c.set == list(cfg.sets)[0]].set_index('hour')
    per_set = {k: c[c.set == k].set_index('hour') for k in cfg.sets}
    rows = [
        (LABEL['T_a_amb'], 'ambient temperature', [(first.T_a_amb, '0.25', 'ambient')], {}),
        (LABEL['w_fr'], 'fan speed (face velocity)', [(d.w_fr, SET_COLOR[k], k) for k, d in per_set.items()],
         {'hlines': [cfg.w_fr_min, cfg.w_fr_design], 'legend': True}),
        (LABEL['P_fan'], 'fan power', [(d.P_fan, SET_COLOR[k], k) for k, d in per_set.items()], {'logy': True}),
        (LABEL['m_dot_ev'], 'water use of the pad', [(d.m_dot_ev * 1000, SET_COLOR[k], k) for k, d in per_set.items()], {}),
    ]
    year_figure((('mean', 'daily mean'), ('max', 'daily maximum')), rows,
                 'Fan control over the year 2025-10-09 … 2026-10-08', 'fan_year', cfg, figsize=(24, 17))


def cost_year(control, cfg):
    """Operating cost over the year with the fan control, electricity (fan + pump) and water separate:
    cumulative curves (left) and monthly cost (right). Prices from config."""
    c = control.copy()
    c['hour'] = pd.to_datetime(c['hour'])
    price_w = cfg.water_eur_per_m3
    fig, axes = plt.subplots(1, 2, figsize=(22, 8))
    months = None
    width = 0.8 / len(cfg.sets)
    for i, key in enumerate(cfg.sets):
        d = c[c.set == key].set_index('hour').sort_index()
        el = (d.P_fan + d.P_pump) / 1e3 * cfg.electricity_eur_per_kwh                  # EUR per hourly row
        wa = d.m_dot_ev * 3.6 * price_w                                                     # kg/s -> m3/h -> EUR
        axes[0].plot(d.index, el.cumsum(), color=SET_COLOR[key], lw=2.2, label=f'{key} electricity (fan + pump)')
        axes[0].plot(d.index, wa.cumsum(), color=SET_COLOR[key], lw=2.2, ls='--', label=f'{key} water')
        m_el, m_wa = el.resample('MS').sum(), wa.resample('MS').sum()
        months = m_el.index
        x = np.arange(len(months)) + (i - (len(cfg.sets) - 1) / 2) * width
        axes[1].bar(x, m_el.values, width, color=SET_COLOR[key], label=f'{key} electricity' if i == 0 or True else None)
        axes[1].bar(x, m_wa.values, width, bottom=m_el.values, color=lighten_color(SET_COLOR[key], 0.55),
                    label=f'{key} water')
    axes[0].set_ylabel('cumulative cost [€]'); axes[0].legend(fontsize=10, ncol=2)
    axes[0].set_title('cumulative cost (water dashed)')
    axes[1].set_xticks(np.arange(len(months))); axes[1].set_xticklabels([m.strftime('%b %y') for m in months], rotation=45)
    axes[1].set_ylabel('cost per month [€]'); axes[1].legend(fontsize=9, ncol=3)
    axes[1].set_title('monthly cost (water bright)')
    fig.suptitle('Operating cost with fan control')
    fig.tight_layout()
    save_figure(fig, 'cost_year')


