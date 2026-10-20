"""
plot_geometry.py -- scaled drawings of the three sized modules (P1-P3).

Row 1: cross-section normal to the tubes (air flows left -> right), same scale
       for all sets, shows d, s_t, s_l and the staggering.
Row 2: section along a tube, same scale for all sets, shows fin pitch s_f,
       fin thickness delta_f, tube d / d_i (fin thickness is drawn to scale).
Row 3: whole module, front view (width b x tube length l) and side view
       (depth N_r s_l x l), common scale in metres.

plot_print: the print-size version (one figure, 6.9 in wide): per cooler a row with the WHOLE tube bank (all tubes, to one scale, with its
       width and depth), the marked detail of 3 x 3 tubes magnified (ONE zoom for all coolers: s_t, s_l) and the fin section (one size and one scale
       for all: d, wall, fin pitch and thickness).
"""

import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Circle, Rectangle

from comparison import evaluate as ev, sizing
from comparison.plots_common import SET_COLOR, save_figure
from src.plots.figures import GEOMETRY
from src.plots.plot_style import PRINT_DPI, PRINT_WIDTH, TUM_GREY_50, TUM_GREY_80, print_context

FIN, TUBE, WATER = '#9aa5b1', '#b87333', '#a8d0f0'


def _dim(ax, p0, p1, text, off=(0, 0), **kw):
    ax.annotate('', xy=p0, xytext=p1, arrowprops=dict(arrowstyle='<->', lw=1.2, color='k'))
    ax.text((p0[0] + p1[0]) / 2 + off[0], (p0[1] + p1[1]) / 2 + off[1], text, fontsize=11, **kw)


def _pitch(ax, a, b, at, start, arrow, horizontal=True, scale=12, lw=1.0):
    """Pitch between positions a and b, measured along one axis: dotted extension lines from `start` to the
    dimension line at `at`, arrows pointing at the two lines from OUTSIDE (readable for tiny pitches: no crossed arrowheads).
    scale / lw: arrowhead size and line width (smaller for the print-size figure)."""
    kw = dict(arrowstyle='-|>', lw=lw, color='k', mutation_scale=scale)
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


def _fin_section(ax, gs, length_mm=18.0, half=12.0, fs=11, scale=12):
    """Cut along the tube axis through the tube centre: tube wall, water, fins. Same axis limits (mm) for every set: the sets are drawn at
    ONE scale when the axes have one size. fs: font size, scale: arrowhead size (smaller for the print-size figure)."""
    d, d_i, p, t = gs.d_mm, gs.d_i_mm, gs.s_f_mm, gs.delta_f_mm
    ax.add_patch(Rectangle((0, -d / 2), length_mm, d, fc=TUBE, ec='none'))
    ax.add_patch(Rectangle((0, -d_i / 2), length_mm, d_i, fc=WATER, ec='none'))
    x = 0.0
    while x + t <= length_mm:
        ax.add_patch(Rectangle((x, -half), t, 2 * half, fc=FIN, ec='k', lw=0.3))
        x += p
    ax.set_xlim(-11, length_mm + 14); ax.set_ylim(-half - 3.5, half + 5); ax.set_aspect('equal')
    _pitch(ax, 0, p, half + 3, lambda v: half, 1.6, scale=scale, lw=fs / 11)
    ax.text(p + 2.2, half + 3, rf'$s_f$ = {p:.3g} mm', fontsize=fs, ha='left', va='center')
    x_f = x - p                                                           # left edge of the last fin
    _pitch(ax, x_f, x_f + t, -half - 1.8, lambda v: -half, 1.2, scale=scale, lw=fs / 11)
    ax.text(x_f + t + 2.2, -half - 1.8, rf'$\delta_f$ = {t:.3g} mm', fontsize=fs, ha='left', va='center')
    for y in (-d / 2, d / 2):
        ax.plot([-3.2, 0], [y, y], color='k', lw=0.6, ls=':')
    ax.annotate('', xy=(-2.6, -d / 2), xytext=(-2.6, d / 2), arrowprops=dict(arrowstyle='<->', lw=1.2 * fs / 11, color='k'))
    ax.text(-3.6, 0, rf'$d$ = {d:.3g} mm', fontsize=fs, ha='right', va='center', rotation=90)
    x_w = length_mm + 1.8                                                  # wall thickness at the right end of the tube
    _pitch(ax, d_i / 2, d / 2, x_w, lambda v: length_mm, 1.0, horizontal=False, scale=scale, lw=fs / 11)
    ax.text(x_w + 1.6, d / 2 + 1.0, rf'$\delta_p$ = {(d - d_i) / 2:.3g} mm', fontsize=fs, ha='left', va='bottom')
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
    with sns.axes_style('white'), sns.plotting_context('talk'):                   # not the global theme
        return _plot(cfg)


def _plot(cfg):
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
    save_figure(fig, 'module_geometries', folder=GEOMETRY)


def plot_bank(cfg):
    """Complete tube bank of each module (all N_t tubes per row, all N_r rows), same scale.
    Air flows left -> right; fins (not drawn) span the whole bank."""
    with sns.axes_style('white'), sns.plotting_context('talk'):
        return _plot_bank(cfg)


def _plot_bank(cfg):
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
    save_figure(fig, 'bank_full', folder=GEOMETRY)


# --- print size ---------------------------------------------------------------------------------------------------------

PRINT_FS = 6.5                              # text size [pt]
BANK_XLIM = (-70.0, 1215.0)                 # mm: one horizontal scale for the three banks (P3, the widest, fills it; room on the left for the depth)
ZOOM_XLIM, ZOOM_YLIM = (-23.0, 146.0), (-52.0, 92.0)       # mm: ONE zoom window for all coolers (P3's 3 x 3 tubes, the region frame at -s_t/2 = -20 mm and the dimensions fit)
ZOOM_ROWS, ZOOM_ACROSS = 3, 3               # the detail: tubes across x rows


def _zoom_x(gs, r, t):
    """x position [mm] of tube t of row r in the detail (rows alternate by s_t / 2)."""
    return t * gs.s_t_mm + (gs.s_t_mm / 2 if r % 2 else 0.0)


def _bank_strip(ax, key, geo, gs, i0):
    """The whole tube bank of one cooler in mm (tubes along x = the width b, rows along y = the depth), with the box that marks the detail
    (ZOOM_ACROSS tubes from tube i0, ZOOM_ROWS rows); b and the depth are dimensioned. Returns the box (x0, x1, y0, y1)."""
    width, depth = geo.b * 1e3, geo.N_r * gs.s_l_mm
    ax.add_patch(Rectangle((0, 0), width, depth, fc=SET_COLOR[key], alpha=0.12, ec=SET_COLOR[key], lw=0.6))
    for r in range(geo.N_r):
        for t in range(geo.N_t):
            ax.add_patch(Circle(((t + 0.5) * gs.s_t_mm + (gs.s_t_mm / 2 if r % 2 else 0.0), (r + 0.5) * gs.s_l_mm), gs.d_mm / 2, fc=TUBE, ec='none'))
    box = (i0 * gs.s_t_mm, (i0 + ZOOM_ACROSS + 0.5) * gs.s_t_mm, 0.0, ZOOM_ROWS * gs.s_l_mm)
    ax.add_patch(Rectangle((box[0], box[2]), box[1] - box[0], box[3] - box[2], fc='none', ec='k', lw=0.7))
    ax.set_xlim(*BANK_XLIM); ax.set_ylim(-50, 250); ax.set_aspect('equal'); ax.axis('off')
    kw = dict(arrowstyle='<->', lw=0.7, color=TUM_GREY_80, shrinkA=0, shrinkB=0)
    ax.annotate('', xy=(0, -22), xytext=(width, -22), arrowprops=kw)
    ax.text(width / 2, -30, rf"$b$ = {width:.0f} mm", fontsize=PRINT_FS, ha='center', va='top', color=TUM_GREY_80)
    ax.annotate('', xy=(-16, 0), xytext=(-16, depth), arrowprops=kw)                       # the depth (thickness) of the bank, on the left:
    ax.text(-26, depth / 2, f"{depth:.0f} mm", fontsize=PRINT_FS, ha='right', va='center', color=TUM_GREY_80, rotation=90)   # the leader lines run to the right
    ax.text(0, 250, key, fontsize=8, fontweight='bold', color=SET_COLOR[key], ha='left', va='top')
    m = ev.calc_masses(geo)
    ax.text(95, 250, rf"{geo.N_t} tubes $\times$ {geo.N_r} rows;  $m_{{Al}}$ = {m['m_Al']:.0f} kg,  $m_{{Cu}}$ = {m['m_Cu']:.0f} kg", fontsize=PRINT_FS, ha='left', va='top',
            color=TUM_GREY_80)                                                       # left: the leader lines run to the right; masses as in the cost model (fins + straight tubes)
    return box


def _zoom(ax, gs):
    """The magnified detail: ZOOM_ACROSS x ZOOM_ROWS tubes (tube and water), s_t below, s_l on the RIGHT; the same mm window for every cooler.
    Returns the region (x0, x1, y0, y1) in the window that corresponds to the box on the bank."""
    d, d_i, s_t, s_l = gs.d_mm, gs.d_i_mm, gs.s_t_mm, gs.s_l_mm
    for r in range(ZOOM_ROWS):
        for t in range(ZOOM_ACROSS):
            ax.add_patch(Circle((_zoom_x(gs, r, t), r * s_l), d / 2, fc=TUBE, ec='k', lw=0.3))
            ax.add_patch(Circle((_zoom_x(gs, r, t), r * s_l), d_i / 2, fc=WATER, ec='k', lw=0.25))
    region = (-0.5 * s_t, (ZOOM_ACROSS + 0.0) * s_t, -0.5 * s_l, (ZOOM_ROWS - 0.5) * s_l)
    ax.add_patch(Rectangle((region[0], region[2]), region[1] - region[0], region[3] - region[2], fc='none', ec='k', lw=0.5))
    ax.set_xlim(*ZOOM_XLIM); ax.set_ylim(*ZOOM_YLIM); ax.set_aspect('equal'); ax.axis('off')
    y_dim = region[2] - 9.0
    _pitch(ax, 0, s_t, y_dim, lambda v: -d / 2, 4.5, scale=7, lw=0.6)
    ax.text(s_t / 2, y_dim - 7.5, rf"$s_t$ = {s_t:.3g} mm", fontsize=PRINT_FS, ha='center', va='top')
    x_dim = region[1] + 8.0
    _pitch(ax, 0, s_l, x_dim, lambda v: _zoom_x(gs, round(v / s_l), ZOOM_ACROSS - 1) + d / 2, 4.5, horizontal=False, scale=7, lw=0.6)
    ax.text(x_dim + 6.0, s_l / 2, rf"$s_l$ = {s_l:.3g} mm", fontsize=PRINT_FS, ha='left', va='center', rotation=90)
    return region


def plot_print(cfg):
    """One row per cooler: the whole tube bank (all to one scale) | the magnified detail (one zoom) | the fin section (one size and scale)."""
    keys = list(cfg.sets)
    with print_context():
        fig = plt.figure(figsize=(PRINT_WIDTH, 3.6))
        grid = fig.add_gridspec(len(keys), 3, width_ratios=[3.0, 1.2, 1.0], hspace=0.12, wspace=0.07, left=0.005, right=0.995, top=0.995, bottom=0.01)
        for row, key in enumerate(keys):
            geo, gs = sizing.get_sized(cfg, key), cfg.sets[key]
            ax_bank, ax_zoom, ax_fin = (fig.add_subplot(grid[row, c]) for c in range(3))
            box = _bank_strip(ax_bank, key, geo, gs, i0=geo.N_t // 2)
            region = _zoom(ax_zoom, gs)
            _fin_section(ax_fin, gs, fs=PRINT_FS, scale=7)
            for (bx, by), (zx, zy) in (((box[1], box[3]), (region[0], region[3])), ((box[1], box[2]), (region[0], region[2]))):
                ax_zoom.annotate('', xy=(bx, by), xycoords=ax_bank.transData, xytext=(zx, zy), textcoords='data', annotation_clip=False,
                                 arrowprops=dict(arrowstyle='-', lw=0.5, color=TUM_GREY_50, shrinkA=0, shrinkB=0), zorder=5)
        save_figure(fig, 'geometry_print', dpi=PRINT_DPI, folder=GEOMETRY)
    return fig


if __name__ == '__main__':
    from comparison.config import settings
    plot(settings())
    plot_bank(settings())
