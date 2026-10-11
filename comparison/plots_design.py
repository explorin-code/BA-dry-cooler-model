"""
plots_design.py -- figures of the design analyses: module data, sensitivity, rows, tubes, wall thickness.
"""

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import seaborn as sns

from src.plot_style import TUM_BLUE, TUM_ORANGE
from comparison.plots_common import DELTA_UNIT, LABEL, SET_COLOR, SYMBOL, save_figure, target_line


def geometries(df, cfg):
    cell = df
    panels = ['A_fr', 'depth', 'A_tot', 'P_fan', 'm_Al', 'm_Cu']
    fig, axes = plt.subplots(2, 3, figsize=(16, 9))
    for ax, col in zip(axes.flat, panels):
        sns.barplot(data=cell, x='set', y=col, hue='set', palette=SET_COLOR, ax=ax, legend=False)
        ax.set_xlabel(''); ax.set_ylabel(LABEL[col])
    fig.suptitle(rf"Modules for $\dot{{m}}_c$ = {cfg.m_dot_c:g} kg/s")
    fig.tight_layout()
    save_figure(fig, 'geometries')


def sensitivity(df, cfg):
    base = df[df.param == 'base'].iloc[0]
    params = [p for p in df.param.unique() if p != 'base']
    fig, axes = plt.subplots(2, 4, figsize=(22, 10))
    for ax, p in zip(axes.flat, params):
        for basis, style in (('w_fr', '-'), ('P_fan', '--')):
            d = df[(df.param == p) & (df.basis == basis)].sort_values('value')
            if d.empty:
                continue
            ax.plot(d.value, d.Q_kW, style, marker='o', ms=4, color=TUM_BLUE,
                    label=r'$w_{fr}$ held' if basis == 'w_fr' else r'$P_{fan}$ held')
            bad = d[~(d.wang_ok & d.vdi_ok)]
            ax.plot(bad.value, bad.Q_kW, 'x', color='r', ms=8)
        ax.set_xlabel(LABEL[p]); ax.set_ylabel(LABEL['Q_kW'])
    axes.flat[0].legend(fontsize=11)
    for ax in axes.flat[len(params):]:
        ax.axis('off')
    fig.suptitle(r"Sensitivity of $\dot{Q}$ (red ×: outside validity range)")
    fig.tight_layout()
    save_figure(fig, 'sensitivity_lines')

    # Tornado: each parameter alone is moved from the base value to the lowest / highest
    # literature value (all others at base) -> change of the result.
    fig, axes = plt.subplots(1, 2, figsize=(17, 6.5))
    for ax, col in zip(axes, ('Q_kW', 'T_c_o')):
        spans = []
        for p in params:
            d = df[(df.param == p) & (df.basis == 'w_fr')].sort_values('value')
            spans.append((p, d[col].iloc[0] - base[col], d[col].iloc[-1] - base[col]))
        spans.sort(key=lambda s: -max(abs(s[1]), abs(s[2])))
        for i, (p, lo, hi) in enumerate(spans):
            ax.barh(i, lo, color=TUM_BLUE); ax.barh(i, hi, color=TUM_ORANGE)
        ax.set_yticks(range(len(spans))); ax.set_yticklabels([SYMBOL[s[0]] for s in spans]); ax.invert_yaxis()
        ax.axvline(0, color='k', lw=1)
        ax.set_xlabel(rf"change of {SYMBOL[col]} vs. base [{DELTA_UNIT[col]}]")
    axes[0].legend(handles=[mpatches.Patch(color=TUM_BLUE, label='parameter at lowest literature value'),
                            mpatches.Patch(color=TUM_ORANGE, label='parameter at highest literature value')],
                   fontsize=11, loc='lower left')
    fig.suptitle(rf"Base module {cfg.sensitivity_set}: $\dot{{Q}}$ = {base.Q_kW:.1f} kW, $T_{{c,o}}$ = {base.T_c_o:.1f} °C")
    fig.tight_layout()
    save_figure(fig, 'sensitivity_tornado')


def rows_analysis(df, cfg):
    """T_c,o, Q, powers, masses vs. number of rows for each elementary cell."""
    d = df[df.case == 'design']
    fig, axes = plt.subplots(2, 3, figsize=(21, 11))
    panels = [('T_c_o', None), ('Q_kW', None), ('P_fan', None), ('m_Al', 'm_Cu'), ('kg_per_kW', None), ('P_pump', None)]
    for ax, (col, col2) in zip(axes.flat, panels):
        for key in cfg.sets:
            x = d[d.set == key]
            ax.plot(x.N_r, x[col], '-o', ms=3, color=SET_COLOR[key], label=key)
            if col2:
                ax.plot(x.N_r, x[col2], '--s', ms=3, color=SET_COLOR[key])
            if col == 'T_c_o':
                cross = x[x.T_c_o <= cfg.T_c_target].head(1)
                if len(cross):
                    ax.plot(cross.N_r, cross[col], '*', ms=16, color=SET_COLOR[key], mec='k')
        ax.set_xlabel(LABEL['N_r']); ax.set_ylabel(LABEL[col] if not col2 else r'$m$ [kg]: solid Al fins, dashed Cu tubes')
        if col == 'T_c_o':
            target_line(ax, cfg)
    axes[0, 0].legend()
    fig.suptitle(r"Row count (star: first $N_r$ reaching the target)")
    fig.tight_layout()
    save_figure(fig, 'rows_analysis')


def tubes_analysis(df, cfg):
    """Effect of the tubes per row N_t (water Reynolds number and face area change with it)."""
    d = df[df.case == 'design']
    fig, axes = plt.subplots(2, 3, figsize=(21, 11))
    panels = ['T_c_o', 'Q_kW', 'Re_c', 'P_fan', 'P_pump', 'kg_per_kW']
    for ax, col in zip(axes.flat, panels):
        for key in cfg.sets:
            x = d[d.set == key].sort_values('N_t')
            ax.plot(x.N_t, x[col], '-o', ms=4, color=SET_COLOR[key], label=key)
            base = x[x.N_t == x.N_t_base]
            ax.plot(base.N_t, base[col], '*', ms=16, color=SET_COLOR[key], mec='k')
        ax.set_xlabel(LABEL['N_t']); ax.set_ylabel(LABEL.get(col, r'$Re_c$ [–]' if col == 'Re_c' else col))
        if col == 'T_c_o':
            target_line(ax, cfg)
        if col == 'Re_c':
            ax.axhspan(0, 4000, color='0.88', zorder=0)
            ax.axhline(2300, color='k', lw=0.8, ls=':')
    axes[0, 0].legend()
    fig.suptitle(r"Tubes per row (star: the sized module, grey: $Re_c$ < 4000)")
    fig.tight_layout()
    save_figure(fig, 'tubes_analysis')


def wall_analysis(df, cfg):
    """Tube wall thickness: N_t fixed (water speeds up as d_i shrinks) vs. Re fixed (N_t rescaled)."""
    d = df[df.case == 'design']
    fig, axes = plt.subplots(2, 4, figsize=(26, 11))
    panels = ['T_c_o', 'Q_kW', 'Re_c', 'P_pump', 'N_t', 'm_Cu', 'kg_per_kW', 'P_fan']
    for ax, col in zip(axes.flat, panels):
        for key in cfg.sets:
            for basis, ls in (('N_t fixed', '-'), ('Re fixed', '--')):
                x = d[(d.set == key) & (d.basis == basis)].sort_values('wall_mm')
                ax.plot(x.wall_mm, x[col], ls, marker='o', ms=3, color=SET_COLOR[key],
                        label=f'{key}, {basis}' if col == 'T_c_o' else None)
                base = x[x.wall_mm == x.wall_base]
                ax.plot(base.wall_mm, base[col], '*', ms=14, color=SET_COLOR[key], mec='k')
        ax.set_xlabel(LABEL['wall_mm']); ax.set_ylabel(LABEL.get(col, r'$Re_c$ [–]' if col == 'Re_c' else col))
        if col == 'T_c_o':
            target_line(ax, cfg)
        if col == 'Re_c':
            ax.axhspan(0, 4000, color='0.88', zorder=0)
    axes[0, 0].legend(fontsize=8, ncol=2)
    fig.suptitle(r"Tube wall thickness $\delta_p$ (solid: $N_t$ fixed, dashed: $Re_c$ kept)")
    fig.tight_layout()
    save_figure(fig, 'wall_analysis')


