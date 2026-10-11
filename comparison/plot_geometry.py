"""
plot_geometry.py -- scaled drawings of the three sized modules (P1-P3).

Row 1: cross-section normal to the tubes (air flows left -> right), same scale
       for all sets, shows d, s_t, s_l and the staggering.
Row 2: section along a tube, same scale for all sets, shows fin pitch s_f,
       fin thickness delta_f, tube d / d_i (fin thickness is drawn to scale).
Row 3: whole module, front view (width b x tube length l) and side view
       (depth N_r s_l x l), common scale in metres.
"""

import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Circle, Rectangle

from comparison import evaluate as ev, sizing
from comparison.plots_common import SET_COLOR, save_figure

FIN, TUBE, WATER = '#9aa5b1', '#b87333', '#a8d0f0'


def _dim(ax, p0, p1, text, off=(0, 0), **kw):
    ax.annotate('', xy=p0, xytext=p1, arrowprops=dict(arrowstyle='<->', lw=1.2, color='k'))
    ax.text((p0[0] + p1[0]) / 2 + off[0], (p0[1] + p1[1]) / 2 + off[1], text, fontsize=11, **kw)


def _pitch(ax, a, b, at, start, arrow, horizontal=True):
    """Pitch between positions a and b, measured along one axis: dotted extension lines from `start` to the
    dimension line at `at`, arrows pointing at the two lines from OUTSIDE (readable for tiny pitches)."""
    kw = dict(arrowstyle='-|>', lw=1.0, color='k', mutation_scale=12)
    for v, sign in ((a, -1), (b, 1)):
        if horizontal:
            ax.plot([v, v], [start(v), at], color='k', lw=0.6, ls=':')
            ax.annotate('', xy=(v, at), xytext=(v + sign * arrow, at), arrowprops=kw)
        else:
            ax.plot([start(v), at], [v, v], color='k', lw=0.6, ls=':')
            ax.annotate('', xy=(at, v), xytext=(at, v + sign * arrow), arrowprops=kw)


def _cross_section(ax, gs, n_across=4):
    """Staggered bank, N_r rows; rows alternate by s_t/2."""
    for r in range(gs.N_r):
        for t in range(n_across):
            x = r * gs.s_l_mm
            y = t * gs.s_t_mm + (gs.s_t_mm / 2 if r % 2 else 0)
            ax.add_patch(Circle((x, y), gs.d_mm / 2, fc=TUBE, ec='k', lw=0.6))
            ax.add_patch(Circle((x, y), gs.d_i_mm / 2, fc=WATER, ec='k', lw=0.4))
    ax.set_xlim(-60, 215); ax.set_ylim(-55, 170); ax.set_aspect('equal')        # one scale, P3 (the largest) fills it
    y_l, x_t = -28.0, -28.0
    _pitch(ax, 0, gs.s_l_mm, y_l, lambda v: (gs.s_t_mm / 2 if v else 0) - gs.d_mm / 2, 8)
    ax.text(gs.s_l_mm / 2, y_l - 8, rf'$s_l$ = {gs.s_l_mm:.4g} mm', fontsize=11, ha='center', va='top')
    _pitch(ax, 0, gs.s_t_mm, x_t, lambda v: -gs.d_mm / 2, 8, horizontal=False)
    ax.text(x_t - 6, gs.s_t_mm / 2, rf'$s_t$ = {gs.s_t_mm:.4g} mm', fontsize=11, ha='right', va='center', rotation=90)
    ax.annotate('air', xy=(210, 160), xytext=(180, 160), arrowprops=dict(arrowstyle='->', lw=1.5, color='k'), va='center', fontsize=11)
    ax.set_title(f"{gs.key}: $d$ = {gs.d_mm:.3g} mm, $d_i$ = {gs.d_i_mm:.3g} mm", fontsize=13)
    ax.set_xlabel('[mm]')


def _fin_section(ax, gs, length_mm=18.0, half=12.0):
    """Cut along the tube axis through the tube centre: tube wall, water, fins."""
    d, d_i, p, t = gs.d_mm, gs.d_i_mm, gs.s_f_mm, gs.delta_f_mm
    ax.add_patch(Rectangle((0, -d / 2), length_mm, d, fc=TUBE, ec='none'))
    ax.add_patch(Rectangle((0, -d_i / 2), length_mm, d_i, fc=WATER, ec='none'))
    x = 0.0
    while x + t <= length_mm:
        ax.add_patch(Rectangle((x, -half), t, 2 * half, fc=FIN, ec='k', lw=0.3))
        x += p
    ax.set_xlim(-11, length_mm + 14); ax.set_ylim(-half - 3.5, half + 5); ax.set_aspect('equal')
    _pitch(ax, 0, p, half + 3, lambda v: half, 1.6)
    ax.text(p + 2.2, half + 3, rf'$s_f$ = {p:.3g} mm', fontsize=11, ha='left', va='center')
    x_f = x - p                                                           # left edge of the last fin
    _pitch(ax, x_f, x_f + t, -half - 1.8, lambda v: -half, 1.2)
    ax.text(x_f + t + 2.2, -half - 1.8, rf'$\delta_f$ = {t:.3g} mm', fontsize=11, ha='left', va='center')
    for y in (-d / 2, d / 2):
        ax.plot([-3.2, 0], [y, y], color='k', lw=0.6, ls=':')
    ax.annotate('', xy=(-2.6, -d / 2), xytext=(-2.6, d / 2), arrowprops=dict(arrowstyle='<->', lw=1.2, color='k'))
    ax.text(-3.6, 0, rf'$d$ = {d:.3g} mm', fontsize=11, ha='right', va='center', rotation=90)
    x_w = length_mm + 1.8                                                  # wall thickness at the right end of the tube
    _pitch(ax, d_i / 2, d / 2, x_w, lambda v: length_mm, 1.0, horizontal=False)
    ax.text(x_w + 1.6, d / 2 + 1.0, rf'$\delta_p$ = {(d - d_i) / 2:.3g} mm', fontsize=11, ha='left', va='bottom')
    ax.axis('off')



def _module(ax, geo, key, w_axis, l_axis, w_fr):
    """Front view (b x l) and side view (depth x l) to one scale."""
    ax.add_patch(Rectangle((0, 0), geo.b, geo.l, fc=SET_COLOR[key], alpha=0.4, ec='k'))
    x_side = geo.b + 0.35
    depth = geo.N_r * geo.s_l
    ax.add_patch(Rectangle((x_side, 0), depth, geo.l, fc=SET_COLOR[key], alpha=0.7, ec='k'))
    m = ev.calc_masses(geo)
    ax.text(geo.b / 2, geo.l / 2, rf'$N_t$ = {geo.N_t}' '\n' rf'$A_{{fr}}$ = {geo.A_fr:.2f} m²' '\n\n' rf'$\sigma = A_e/A_{{fr}}$ = {geo.sigma:.2f}' '\n'
            rf'$w_e = w_{{fr}}/\sigma$ = {w_fr / geo.sigma:.1f} m/s' '\n\n'
            rf'$m_{{Al}}$ = {m["m_Al"]:.0f} kg' '\n' rf'$m_{{Cu}}$ = {m["m_Cu"]:.0f} kg', ha='center', va='center', fontsize=12)
    _dim(ax, (-0.15, 0), (-0.15, geo.l), rf'$l$ = {geo.l:.2f} m', off=(-0.05, 0), ha='right', va='center', rotation=90)
    y_dim = geo.l + 0.1
    _dim(ax, (0, y_dim), (geo.b, y_dim), '')
    ax.text(geo.b / 2, y_dim + 0.06, rf'$b$ = {geo.b:.2f} m', ha='center', fontsize=11)
    _pitch(ax, x_side, x_side + depth, y_dim, lambda v: geo.l, 0.1)
    ax.text(x_side + depth / 2, y_dim + 0.06, rf'{depth * 1e3:.0f} mm', ha='center', fontsize=11)
    ax.annotate('air', xy=(x_side + depth + 0.2, geo.l / 2), xytext=(x_side - 0.2, geo.l / 2),
                arrowprops=dict(arrowstyle='->', lw=1.5, color='k'), va='center', fontsize=11)
    ax.set_xlim(-0.6, w_axis); ax.set_ylim(-0.3, l_axis); ax.set_aspect('equal')
    ax.set_title(key, fontsize=13)
    ax.axis('off')


def plot(cfg):
    sns.set_theme(style='white', context='talk')
    n = len(cfg.sets)
    fig, axes = plt.subplots(3, n, figsize=(5.6 * n, 15), gridspec_kw=dict(height_ratios=[1.3, 0.7, 1.5]), squeeze=False)
    geos = {k: sizing.get_sized(cfg, k) for k in cfg.sets}
    w_axis = max(g.b for g in geos.values()) + 0.35 + 0.7
    l_axis = max(g.l for g in geos.values()) + 0.7
    for j, (key, gs) in enumerate(cfg.sets.items()):
        _cross_section(axes[0, j], gs)
        _fin_section(axes[1, j], gs)
        _module(axes[2, j], geos[key], key, w_axis, l_axis, cfg.w_fr_design)
        axes[1, j].set_title(rf"fin section, $N_f$ = {geos[key].N_f}", fontsize=13)
    fig.suptitle(f"Modules for {cfg.m_dot_c:g} kg/s: tube bank, fin section, module size")
    fig.tight_layout()
    save_figure(fig, 'module_geometries')


def plot_bank(cfg):
    """Complete tube bank of each module (all N_t tubes per row, all N_r rows), same scale.
    Air flows left -> right; fins (not drawn) span the whole bank."""
    sns.set_theme(style='white', context='talk')
    geos = {k: sizing.get_sized(cfg, k) for k in cfg.sets}
    fig, axes = plt.subplots(1, len(cfg.sets), figsize=(5.6 * len(cfg.sets), 15), sharey=True, squeeze=False)
    for ax, (key, geo) in zip(axes[0], geos.items()):
        gs = cfg.sets[key]
        for r in range(geo.N_r):
            for t in range(geo.N_t):
                x = (r + 0.5) * gs.s_l_mm
                y = (t + 0.5) * gs.s_t_mm + (gs.s_t_mm / 2 if r % 2 else 0)
                ax.add_patch(Circle((x, y), gs.d_mm / 2, fc=TUBE, ec='k', lw=0.2))
                ax.add_patch(Circle((x, y), gs.d_i_mm / 2, fc=WATER, ec='none'))
        ax.add_patch(Rectangle((0, 0), geo.N_r * gs.s_l_mm, geo.b * 1e3, fill=False, ec='k', lw=1))
        ax.set_xlim(-20, 380); ax.set_ylim(-20, 1400); ax.set_aspect('equal')
        ax.set_title(f"{key}: $N_t$ = {geo.N_t}, $N_r$ = {geo.N_r}",
                     fontsize=13)
        ax.set_xlabel('depth, air →  [mm]')
    axes[0, 0].set_ylabel('width [mm]')
    fig.suptitle('Tube banks (same scale)')
    fig.tight_layout()
    save_figure(fig, 'bank_full')


if __name__ == '__main__':
    from comparison.config import settings
    plot(settings())
    plot_bank(settings())
