"""
plots_common.py -- what all comparison figures share: colours, thesis-symbol labels, saving.

Figures use the thesis symbols (mathtext); the CSV columns of the analyses keep the code names and LABEL maps
one to the other. `save_figure` writes the PNG (src/figures.py); the figure stays open, so main.py (plt.show) and the
GUI (collects open figures) can display it.
"""

import matplotlib.pyplot as plt
import seaborn as sns

from src import figures
from src.plot_style import TUM_BLUE, TUM_GREEN, TUM_ORANGE, lighten_color

SET_COLOR = {'P1': TUM_BLUE, 'P2': TUM_ORANGE, 'P3': TUM_GREEN}      # TUM Corporate Design
sns.set_theme(style='whitegrid', context='talk')

# CSV column / parameter name -> thesis symbol with unit
LABEL = {
    'Q_kW': r'$\dot{Q}$ [kW]', 'T_c_o': r'$T_{c,o}$ [°C]', 'T_a_o': r'$T_{a,o}$ [°C]',
    'pinch': r'$\Delta T_{pinch}$ [K]', 'P_fan': r'$P_{fan}$ [W]', 'P_pump': r'$P_{pump}$ [W]',
    'A_fr': r'$A_{fr}$ [m²]', 'A_tot': r'$A_{tot}$ [m²]', 'depth': r'depth $N_r\,s_l$ [m]',
    'm_total': r'$m_{Al+Cu}$ [kg]', 'm_Al': r'$m_{Al}$ (fins) [kg]', 'm_Cu': r'$m_{Cu}$ (tubes) [kg]',
    'kg_per_kW': r'$m_{Al+Cu}/\dot{Q}$ [kg/kW]', 'N_t': r'$N_t$ [–]', 'N_r': r'$N_r$ [–]',
    'w_c': r'$w_c$ [m/s]', 'w_fr': r'$w_{fr}$ [m/s]', 'Re_c': r'$Re_c$ [–]', 'T_a_amb': r'ambient $T_{a,i}$ [°C]',
    'm_dot_ev': r'$\dot{m}_{water}$ [g/s]',
    'wall_mm': r'$\delta_p$ [mm]', 'd_mm': r'$d$ [mm]', 's_f_mm': r'$s_f$ [mm]', 'delta_f_mm': r'$\delta_f$ [mm]',
    's_t_mm': r'$s_t$ [mm]', 's_l_mm': r'$s_l$ [mm]',
}
SYMBOL = {k: v.split(' [')[0] for k, v in LABEL.items()}      # symbol without unit
DELTA_UNIT = {'Q_kW': 'kW', 'T_c_o': 'K'}

def save_figure(fig, name):
    """PNG in thesis_figures/out/ch4/, window title = the name (figures stay open; the runner closes them)."""
    figures.save_figure(fig, 'ch4', name, title=name.replace('_', ' ').capitalize())


def mode_color(key, mode):
    """Set colour; precooled operation in a brighter shade of it."""
    return SET_COLOR[key] if mode == 'dry' else lighten_color(SET_COLOR[key], 0.5)


def daily_line(ax, series, color, label=None, stat='mean', lw=1.6):
    """Daily mean or daily maximum of an hourly series (DatetimeIndex) as a line."""
    ax.plot(series.resample('D').agg(stat).index, series.resample('D').agg(stat).values, color=color, lw=lw, label=label)


def year_figure(columns, rows, suptitle, name, cfg, figsize=(24, 13)):
    """rows: [(ylabel, title, [(series, colour, label)], {'logy': bool, 'target': bool, 'hlines': [..]})];
    one column per statistic in `columns` (('mean', 'daily mean'), ('max', 'daily maximum'))."""
    fig, axes = plt.subplots(len(rows), len(columns), figsize=figsize, sharex=True, sharey='row')
    for i, (ylabel, title, lines, opt) in enumerate(rows):
        for j, (stat, stat_name) in enumerate(columns):
            ax = axes[i, j]
            for series, color, label in lines:
                daily_line(ax, series, color, label if i == len(rows) - 1 or opt.get('legend') else None, stat)
            if opt.get('target'):
                target_line(ax, cfg)
            for h in opt.get('hlines', []):
                ax.axhline(h, color='k', ls=':', lw=1)
            if opt.get('logy'):
                ax.set_yscale('log')
            if j == 0:
                ax.set_ylabel(ylabel)
            ax.set_title(f'{title}: {stat_name}')
    for ax in axes.flat:
        if ax.get_legend_handles_labels()[0]:
            ax.legend(fontsize=11, ncol=3)
    fig.suptitle(suptitle)
    fig.tight_layout()
    save_figure(fig, name)


def target_line(ax, cfg, horizontal=True):
    label = rf'target $T_{{c,o}}$ = {cfg.T_c_target:.0f} °C'
    (ax.axhline if horizontal else ax.axvline)(cfg.T_c_target, color='k', ls='--', lw=1, label=label)


