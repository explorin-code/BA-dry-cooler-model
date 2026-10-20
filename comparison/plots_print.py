"""
plots_print.py -- the print-size figures of the cooler comparison (6.9 in wide, 300 dpi, plot_style.print_context; no figure titles,
one small corner label for the constant that matters, every graph with its own legend). The large-format figures are in
plots_operation.py / plots_design.py.

  pinch_year       pinch of the whole cooler and water outlet temperature against the ambient temperature, fixed face velocity
  pinch_speed      pinch against the face velocity at three ambient states (median, the design state, P95), dry and precooled
  pinch_year_clouds  the appendix version of pinch_year: every hour of the year as a point (dry and precooled), one column per cooler, pinch above T_c,o
  annual_costs     energy (fan, pump) and cost (electricity, water, material) per year for three operating modes (the groups): dry at the fixed speed / dry
                   with the adaptive fan / wet (the adaptive fan with the pad); the coolers side by side in each group; the pad water as a small graph
  monthly_cost     the cost per month of the adaptive fan (electricity, water)
  water_vs_hours   what the pad water buys: hours above the target against the water used, per operating mode
  missed_hours     the year at the fixed face velocity: T_c,o against the ambient temperature per cooler, the hours that miss the target highlighted,
                   with the duration curve of the missed part (dry / precooled) as an inset and the hours counted
  fan_cop          the heat removed per fan watt (Q / P_fan) of the adaptive fan against the ambient temperature
  variation_costs / variation_totals   pinch (design state) and annual cost of the adaptive fan against six geometry parameters (N_r, N_t, l / s_l, s_f, delta_f),
                   P1-P3; the variants that miss the target for more than MISSED_H hours dotted, those outside the Wang / VDI ranges with open markers; the first
                   with the cost of ONE cooler split into fan / pump / water / material (main text), the second with the totals of all three (appendix)
  controlled_fan   face velocity and fan power of the variable-speed fan against the ambient temperature (outlet held at the target; dry = line, precooled = every hour as a point + an average line)

Pinch = T_c,o - the ambient air that enters the WHOLE cooler (evaluate.run_point); dry curves solid, precooled dashed, P1-P3 in their colours.
"""

import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.legend_handler import HandlerTuple
from matplotlib.patches import Patch
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

from comparison import operation, variation
from comparison.config import TARGET_TOL
from comparison.plots_common import SET_COLOR, save_figure
from src import weather
from src.plot_style import PRINT_DASHED, PRINT_DPI, PRINT_LW, PRINT_WIDTH, TUM_GREY_20, TUM_GREY_50, TUM_GREY_80, lighten_color, print_context

HEIGHT = 2.5                                  # one row of graphs [in]
BIN = 2.0                                     # temperature bins of the year curves [K]
MODE_LS = {'dry': '-', 'precooled': PRINT_DASHED}


def binned(x, y, lo=5, hi=95, min_n=5, center=np.median):
    """(bin centres, centre (median, or e.g. np.mean), low percentile, high percentile) of y in BIN-wide bins of x (bins with fewer than min_n points are dropped)."""
    edges = np.arange(np.floor(np.min(x) / BIN) * BIN, np.max(x) + BIN, BIN)
    index = np.digitize(x, edges) - 1
    out = [(edges[i] + BIN / 2, center(y[index == i]), np.percentile(y[index == i], lo), np.percentile(y[index == i], hi))
           for i in range(len(edges) - 1) if (index == i).sum() >= min_n]
    return tuple(np.array(column) for column in zip(*out)) if out else (np.array([]),) * 4


def legend_handles(sets):
    return ([Line2D([], [], color=SET_COLOR[k], lw=PRINT_LW, label=k) for k in sets]
            + [Line2D([], [], color='k', lw=PRINT_LW, label="dry"), Line2D([], [], color='k', lw=PRINT_LW, ls=PRINT_DASHED, label="precooled")])


def _tight_layout(fig, **kwargs):
    """fig.tight_layout for figures with an inset graph: tight_layout leaves the inset out (expected), so its warning is not shown."""
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', UserWarning)
        fig.tight_layout(**kwargs)


def corner(ax, text, loc='lower left'):
    """The one constant that matters, small in a corner of the graph."""
    x, ha = (0.03, 'left') if 'left' in loc else (0.97, 'right')
    y, va = (0.04, 'bottom') if 'lower' in loc else (0.96, 'top')
    ax.text(x, y, text, transform=ax.transAxes, ha=ha, va=va, fontsize=7, color=TUM_GREY_80,
            bbox=dict(facecolor='white', edgecolor='none', alpha=0.85, pad=1.5))


def _finish(fig, name):
    fig.set_size_inches(PRINT_WIDTH, HEIGHT)
    fig.tight_layout(w_pad=1.0)
    save_figure(fig, name, dpi=PRINT_DPI)


def pinch_year(df, cfg):
    """df: operation.year_modes. Left: pinch against T_amb, right: T_c,o against T_amb (target dotted); median per bin."""
    with print_context():
        fig, axes = plt.subplots(1, 2)
        for key in cfg.sets:
            for mode, ls in MODE_LS.items():
                d = df[(df.set == key) & (df['mode'] == mode)]
                for ax, col in zip(axes, ('pinch', 'T_c_o')):
                    x, m, *_ = binned(d.T_a_amb.values, d[col].values)
                    ax.plot(x, m, color=SET_COLOR[key], lw=PRINT_LW, ls=ls)
        axes[1].axhline(cfg.T_c_target, color=TUM_GREY_50, lw=0.9, ls=':', zorder=0)
        axes[1].annotate(f"target {cfg.T_c_target:g} °C", (0.02, cfg.T_c_target), xycoords=('axes fraction', 'data'),
                         xytext=(0, 3), textcoords='offset points', fontsize=6, color=TUM_GREY_50)
        for ax in axes:
            ax.set_xlabel(r"$T_{amb}$ [°C]")
        axes[0].set_ylabel(r"$\Delta T_{pinch}$ [K]")
        axes[1].set_ylabel(r"$T_{c,o}$ [°C]")
        axes[0].legend(handles=legend_handles(cfg.sets), loc='upper right', ncol=2)
        corner(axes[0], rf"$w_{{fr}}$ = {cfg.w_fr_design:g} m/s")
        _finish(fig, 'pinch_year')
    return fig


def _state_titles(cfg, sweep):
    """Panel titles of the three ambient states: 'median (11 °C)', 'P80 (20 °C)' (the design state: its percentile in the year), 'P95 (28 °C)'."""
    T = weather.load_year()['T_air_C'].values
    names = {'median': "median", 'design': f"P{100 * (T <= cfg.design_T_air).mean():.0f}", 'P95': "P95"}
    return {case: f"{name} ({sweep[sweep.case == case].T_a_amb.iloc[0]:.0f} °C)" for case, name in names.items()}


def pinch_speed(sweep, cfg):
    """sweep: operation.airflow_sweep. One panel per ambient state: pinch against the face velocity, dry and precooled."""
    titles = _state_titles(cfg, sweep)
    with print_context():
        fig, axes = plt.subplots(1, len(titles), sharey=True)
        for ax, (case, title) in zip(axes, titles.items()):
            for key in cfg.sets:
                for mode, ls in MODE_LS.items():
                    d = sweep[(sweep.case == case) & (sweep.set == key) & (sweep['mode'] == mode)].sort_values('w_fr')
                    ax.plot(d.w_fr, d.pinch, color=SET_COLOR[key], lw=PRINT_LW, ls=ls)
            ax.set_title(title)
            ax.set_xlabel(r"$w_{fr}$ [m/s]")
        axes[0].set_ylabel(r"$\Delta T_{pinch}$ [K]")
        axes[0].legend(handles=legend_handles(cfg.sets), loc='upper right', ncol=2)
        _finish(fig, 'pinch_speed')
    return fig


POINTS = dict(s=1.6, alpha=0.30, linewidths=0, rasterized=True)       # the hours of a scatter: small, transparent, light on the page
HALO_LW = PRINT_LW + 0.8                                               # a white solid line under a dashed one keeps it visible over the points


def _average_line(ax, data, y, color, mode):
    """The average of `y` (a column of `data`: control rows of ONE cooler and mode) per ambient-temperature bin as a line; precooled lines get a white
    outline, the bins in which most hours miss the target (`limited`) are dotted."""
    x, m, *_ = binned(data.T_a_amb.values, y.values, center=np.mean)
    missed = binned(data.T_a_amb.values, data.limited.values.astype(float), center=np.mean)[1] > 0.5
    first = int(np.argmax(missed)) if missed.any() else len(x)
    parts = [(slice(0, first), MODE_LS[mode])] + ([(slice(max(first - 1, 0), None), ':')] if first < len(x) else [])
    for part, ls in parts:
        if mode == 'precooled':
            ax.plot(x[part], m[part], color='white', lw=HALO_LW, zorder=3)
        ax.plot(x[part], m[part], color=color, lw=PRINT_LW, ls=ls, zorder=4)


COP_SHADES = (('pump', 0.7, TUM_GREY_20, r"pump: $\dot{Q}\,/\,P_{pump}$"),                   # (column of P, lightening of the cooler colour, legend grey, legend label)
              ('fan', 0.45, TUM_GREY_50, r"fan: $\dot{Q}\,/\,P_{fan}$"),
              ('total', 0.0, TUM_GREY_80, r"total: $\dot{Q}\,/\,(P_{fan}+P_{pump})$"))


def fan_cop(control, cfg):
    """The heat removed per electrical watt against the ambient temperature, one line per cooler and mode (average per bin; dotted where most hours
    miss the target) in three shades: Q / (P_fan + P_pump) (the total COP, dark), Q / P_fan (the fan 'COP') and Q / P_pump (the pump alone, lightest).
    Huge in the cold, where the fan hardly runs -- there the pump dominates the total --, and a few tens near the limit."""
    power = dict(fan=lambda d: d.P_fan, total=lambda d: d.P_fan + d.P_pump, pump=lambda d: d.P_pump)
    with print_context():
        fig, ax = plt.subplots(figsize=(PRINT_WIDTH * 0.62, 2.7))
        for key in cfg.sets:
            for mode in MODE_LS:
                d = control[(control.set == key) & (control['mode'] == mode)]
                if not d.empty:
                    for name, light, _, _ in COP_SHADES:                      # the lightest first: the dark lines lie on top
                        _average_line(ax, d, d.Q_kW * 1e3 / power[name](d), lighten_color(SET_COLOR[key], light), mode)
        ax.set_yscale('log')
        ax.set_ylim(top=1.2e6)                                                    # room for the legend above the pump lines
        ax.set_xlabel(r"$T_{amb}$ [°C]")
        ax.set_ylabel(r"$\dot{Q}\,/\,P$ [–]")
        sets, (dry, precooled) = legend_handles(cfg.sets)[:-2], legend_handles(cfg.sets)[-2:]
        shades = {name: Line2D([], [], color=grey, lw=PRINT_LW, label=label) for name, _, grey, label in COP_SHADES}
        missed = Line2D([], [], color='k', lw=PRINT_LW, ls=':', label="target missed")
        ax.legend(handles=sets + [dry, precooled, missed, shades['total'], shades['fan'], shades['pump']], loc='upper right', ncol=3)
        corner(ax, rf"adaptive fan, $T_{{c,o}}$ = {cfg.T_c_target:g} °C", 'lower left')
        fig.tight_layout()
        save_figure(fig, 'fan_cop', dpi=PRINT_DPI)
    return fig


def controlled_fan(control, cfg):
    """control: operation.year_control. Face velocity (linear) and fan power (logarithmic) against the AMBIENT temperature. Dry hours: a line (the dry
    control is a function of T_amb). Precooled hours: every hour as a point -- their spread is real (the pad's effect depends on the humidity) --
    plus one average line per cooler. The hours at the top speed that still miss the target (`limited`) are included (the fan runs at full
    speed); the averages of the bins in which most hours miss the target are dotted."""
    with print_context():
        fig, axes = plt.subplots(1, 2)
        for ax, col, label in ((axes[0], 'w_fr', r"$w_{fr}$ [m/s]"), (axes[1], 'P_fan', r"$P_{fan}$ [W]")):
            for key in cfg.sets:
                for mode in ('precooled', 'dry'):
                    d = control[(control.set == key) & (control['mode'] == mode)]
                    if d.empty:
                        continue
                    if mode == 'precooled':
                        ax.scatter(d.T_a_amb, d[col], color=SET_COLOR[key], **POINTS)
                    _average_line(ax, d, d[col], SET_COLOR[key], mode)
            ax.set_xlabel(r"$T_{amb}$ [°C]")
            ax.set_ylabel(label)
        axes[1].set_yscale('log')
        axes[0].legend(handles=[Line2D([], [], color=SET_COLOR[k], lw=PRINT_LW, label=k) for k in cfg.sets]
                       + [Line2D([], [], color='k', lw=PRINT_LW, label="dry"),
                          Line2D([], [], color='k', lw=PRINT_LW, ls=PRINT_DASHED, label="precooled: average"),
                          Line2D([], [], color='0.5', marker='o', ls='none', ms=2.5, alpha=0.6, label="precooled: hours"),
                          Line2D([], [], color='k', lw=PRINT_LW, ls=':', label="target missed")],
                       loc='upper left', ncol=2, columnspacing=0.8)
        corner(axes[0], rf"$T_{{c,o}}$ = {cfg.T_c_target:g} °C", 'lower right')
        _finish(fig, 'controlled_fan')
    return fig


MISSED = '#e88e8e'                            # the hours that miss the target: a soft red


def missed_hours(df, cfg):
    """df: operation.year_modes (every hour in both modes at the fixed design speed). One panel per cooler: T_c,o against T_amb in the adopted operation
    (dry; precooled only for the hours in which dry misses the target), the missed hours in red; in the lower right an inset with the duration curve of the
    excess T_c,o - target for dry and for precooled operation (only the part above zero), its legend giving the hours and their share of the year."""
    tol = TARGET_TOL
    with print_context():
        fig, axes = plt.subplots(1, len(cfg.sets), sharey=True)
        for ax, key in zip(axes, cfg.sets):
            d = df[df.set == key]
            dry, wet = d[d['mode'] == 'dry'], d[d['mode'] == 'precooled']
            dry_missed = dry[dry.T_c_o > cfg.T_c_target + tol]
            adopted = pd.concat([dry[dry.T_c_o <= cfg.T_c_target + tol], wet[wet.hour.isin(dry_missed.hour)]])
            met, miss = adopted[adopted.T_c_o <= cfg.T_c_target + tol], adopted[adopted.T_c_o > cfg.T_c_target + tol]
            ax.scatter(met.T_a_amb, met.T_c_o, s=1.6, color=lighten_color(SET_COLOR[key], 0.45), linewidths=0, rasterized=True)
            ax.scatter(miss.T_a_amb, miss.T_c_o, s=2.2, color=MISSED, linewidths=0, rasterized=True)
            ax.axhline(cfg.T_c_target, color=TUM_GREY_80, lw=0.9, ls='--', zorder=0)
            ax.set_title(key, color=SET_COLOR[key], fontweight='bold')
            ax.set_xlabel(r"$T_{amb}$ [°C]")

            ins = inset_axes(ax, width="100%", height="100%", bbox_to_anchor=(0.62, 0.075, 0.355, 0.34), bbox_transform=ax.transAxes, loc='lower left', borderpad=0)
            n = len(d) // 2                                                      # hours of the year
            for data, ls, name in ((dry, '-', "dry"), (wet, PRINT_DASHED, "precooled")):
                excess = np.sort((data.T_c_o - cfg.T_c_target).values)[::-1]
                hours = int((excess > tol).sum())
                ins.plot(np.arange(1, hours + 1), excess[:hours], color=SET_COLOR[key], lw=PRINT_LW, ls=ls,
                         label=f"{name} {hours} h ({100 * hours / n:.0f} %)")
            ins.axhline(0, color=TUM_GREY_80, lw=0.7, ls='--', zorder=0)
            ins.set_xlim(left=0); ins.set_ylim(bottom=0)
            ins.set_xticks([0, 1000], ["0 h", "1000 h"])               # the units in the tick labels, no axis titles
            ins.set_yticks([0, 5, 10], ["0 K", "+5 K", "+10 K"])
            ins.tick_params(labelsize=5, pad=1, length=1.5)
            ins.tick_params(axis='x', labelrotation=0)
            ins.grid(True, lw=0.3)
            ins.legend(loc='lower right', bbox_to_anchor=(1.0, 1.03), fontsize=5.5, borderpad=0.25, labelspacing=0.15, handlelength=1.2, handletextpad=0.3)
        axes[0].set_ylabel(r"$T_{c,o}$ [°C]")
        axes[0].legend(handles=[Line2D([], [], color='0.65', marker='o', ls='none', ms=2.5, label="target met"),
                                Line2D([], [], color=MISSED, marker='o', ls='none', ms=2.5, label="target missed")], loc='upper left')
        corner(axes[-1], rf"$w_{{fr}}$ = {cfg.w_fr_design:g} m/s", 'upper left')
        fig.set_size_inches(PRINT_WIDTH, 2.7)
        _tight_layout(fig, w_pad=0.8)
        save_figure(fig, 'missed_hours', dpi=PRINT_DPI)
    return fig


def pinch_year_clouds(df, cfg):
    """df: operation.year_modes. The appendix version of pinch_year: one column per cooler, the pinch (top) and the water outlet temperature (bottom)
    against T_amb with EVERY hour of the year as a point, dry and precooled (the clouds that pinch_year reduces to median lines), plus the median line per
    2 K bin of each mode. The dry hours lie on a thin band, the precooled ones spread with the humidity."""
    with print_context():
        fig, axes = plt.subplots(2, len(cfg.sets), sharex=True, sharey='row')
        for col, key in enumerate(cfg.sets):
            for row, (name, label) in enumerate((('pinch', r"$\Delta T_{pinch}$ [K]"), ('T_c_o', r"$T_{c,o}$ [°C]"))):
                ax = axes[row, col]
                for mode, shade in (('precooled', 0.3), ('dry', 0.55)):
                    d = df[(df.set == key) & (df['mode'] == mode)]
                    ax.scatter(d.T_a_amb, d[name], s=1.3, color=lighten_color(SET_COLOR[key], shade), linewidths=0, rasterized=True)
                for mode, ls in MODE_LS.items():
                    d = df[(df.set == key) & (df['mode'] == mode)]
                    x, m, *_ = binned(d.T_a_amb.values, d[name].values)
                    ax.plot(x, m, color='white', lw=HALO_LW, zorder=3)
                    ax.plot(x, m, color=SET_COLOR[key], lw=PRINT_LW, ls=ls, zorder=4)
                if name == 'T_c_o':
                    ax.axhline(cfg.T_c_target, color=TUM_GREY_50, lw=0.9, ls=':', zorder=0)
                if row == 0:
                    ax.set_title(key, color=SET_COLOR[key], fontweight='bold')
                if row == 1:
                    ax.set_xlabel(r"$T_{amb}$ [°C]")
                if col == 0:
                    ax.set_ylabel(label)
        axes[0, 0].legend(handles=[Line2D([], [], color='k', lw=PRINT_LW, label="dry: median"),
                                   Line2D([], [], color='k', lw=PRINT_LW, ls=PRINT_DASHED, label="precooled: median"),
                                   Line2D([], [], color='0.6', marker='o', ls='none', ms=2.5, label="every hour")], loc='upper right')
        corner(axes[1, 0], rf"$w_{{fr}}$ = {cfg.w_fr_design:g} m/s", 'lower right')
        fig.set_size_inches(PRINT_WIDTH, 4.6)
        fig.tight_layout(w_pad=0.8, h_pad=0.8)
        save_figure(fig, 'pinch_year_clouds', dpi=PRINT_DPI)
    return fig


# --- the annual economics -------------------------------------------------------------------------------------------------------------------

SCENARIOS = (("dry, fixed 3 m/s", 0.6), ("dry, adaptive fan", 0.3), ("wet, adaptive fan", 0.0))      # (name, lightening of the cooler colour)
STACK = dict(pump=TUM_GREY_20, water='#7B4B94', material='#2A9D8F')        # purple / teal: nothing stacked may look like a cooler's colour (blue, orange, green)


def annual_summaries(year, control, control_dry, cfg) -> dict:
    """{scenario name: operation.annual_summary} of the three operating modes of one year: the dry cooler at the fixed design speed, the dry cooler with the
    adaptive fan (smallest sufficient speed, no pad), and the adaptive fan with the pad when dry is not enough ('wet').
    year = operation.year_modes rows, control / control_dry = operation.year_control rows with / without the pad."""
    frames = {SCENARIOS[0][0]: year[year['mode'] == 'dry'], SCENARIOS[1][0]: control_dry, SCENARIOS[2][0]: control}
    return {name: operation.annual_summary(cfg, frames[name], hours_per_row=float(cfg.year_stride)) for name, _ in SCENARIOS}


def annual_costs(year, control, control_dry, cfg):
    """Two graphs grouped by operating mode (dry fixed / dry adaptive / wet adaptive), the coolers side by side in their colours inside every group: the
    energy (fan, pump stacked) and the cost (electricity, pad water, material write-off stacked). The free upper right of the energy graph holds a small graph
    of the pad water of the wet mode (its bars in the coolers' colours, no labels of their own). The legends are inside the graphs, the cooler legend at the
    upper right of the energy graph. The costs come from operation.annual_summary / add_costs (prices: parameters.py),
    so a change of the prices or of the cost model needs only a redraw (`--plots-only`)."""
    summaries = annual_summaries(year, control, control_dry, cfg)
    names = [name for name, _ in SCENARIOS]
    w = 0.26
    patch = lambda color, label=None: Patch(facecolor=color, edgecolor='black', linewidth=0.3, label=label)
    cooler_swatch = tuple(patch(SET_COLOR[k]) for k in cfg.sets)                     # 'the coolers' colours' as one legend symbol
    with print_context():
        fig, axes = plt.subplots(1, 2)
        for i, name in enumerate(names):
            for j, key in enumerate(cfg.sets):
                x, r, color = i + (j - (len(cfg.sets) - 1) / 2) * w, summaries[name].loc[key], SET_COLOR[key]
                kw = dict(width=w * 0.94, edgecolor='black', linewidth=0.3, zorder=3)
                axes[0].bar(x, r.E_fan_kWh, color=color, **kw)
                axes[0].bar(x, r.E_pump_kWh, bottom=r.E_fan_kWh, color=STACK['pump'], **kw)
                axes[1].bar(x, r.cost_el, color=color, **kw)
                axes[1].bar(x, r.cost_water, bottom=r.cost_el, color=STACK['water'], **kw)
                axes[1].bar(x, r.cost_material, bottom=r.cost_el + r.cost_water, color=STACK['material'], **kw)
        for ax, label in zip(axes, (r"energy [kWh/a]", r"cost [€/a]")):
            ax.set_xticks(range(len(names)), [n.replace(", ", ",\n") for n in names])
            ax.set_ylabel(label)
            ax.grid(False)
            ax.grid(True, axis='y')
        axes[1].set_ylim(top=axes[1].get_ylim()[1] * 1.12)

        # the pad water of the wet mode: a small graph in the free upper right of the energy graph
        inset = inset_axes(axes[0], width="100%", height="100%", bbox_to_anchor=(0.40, 0.50, 0.36, 0.46), bbox_transform=axes[0].transAxes,
                           loc='lower left', borderpad=0)
        water = [summaries[names[-1]].loc[k].water_m3 for k in cfg.sets]
        inset.bar(range(len(water)), water, 0.7, color=[SET_COLOR[k] for k in cfg.sets], edgecolor='black', linewidth=0.3, zorder=3)
        for i, v in enumerate(water):
            inset.annotate(f"{v:.0f}", (i, v), xytext=(0, 1.5), textcoords='offset points', ha='center', fontsize=5.5, color=TUM_GREY_80)
        inset.set_xticks([])                                                      # the bars carry the coolers' colours (legend of the graph)
        inset.set_title("wet mode: pad water [m³/a]", fontsize=5.8, pad=2)
        inset.set_ylim(0, max(water) * 1.25)
        inset.tick_params(labelsize=5, pad=1, length=1.5)
        inset.grid(False)
        inset.grid(True, axis='y', lw=0.3)

        # the legends, inside the graphs (in the free space above the adaptive modes)
        axes[0].legend(handles=[patch(SET_COLOR[k], k) for k in cfg.sets] + [patch(STACK['pump'], "pump")], loc='upper right', ncol=1, fontsize=6, borderpad=0.4)
        axes[1].legend(handles=[cooler_swatch, patch(STACK['water']), patch(STACK['material'])],
                       labels=["electricity", "pad water", f"material ({cfg.lifetime_years:g} a)"], handler_map={tuple: HandlerTuple(ndivide=None, pad=0)},
                       loc='upper right', ncol=1, fontsize=6, borderpad=0.4, handlelength=2.4)
        fig.set_size_inches(PRINT_WIDTH, 2.7)
        _tight_layout(fig, w_pad=0.9)
        save_figure(fig, 'annual_costs', dpi=PRINT_DPI)
    return fig


def monthly_cost(control, cfg):
    """The cost per month of the adaptive fan, one panel per cooler: electricity (fan + pump) and pad water stacked."""
    with print_context():
        fig, axes = plt.subplots(1, len(cfg.sets), sharey=True)
        for ax, key in zip(axes, cfg.sets):
            d = control[control.set == key].copy()
            d['month'] = d.hour.dt.to_period('M')
            g = d.groupby('month').agg(fan=('P_fan', 'sum'), pump=('P_pump', 'sum'), water=('m_dot_ev', 'sum'))
            electricity = (g.fan + g.pump) / 1e3 * cfg.year_stride * cfg.electricity_eur_per_kwh
            water = g.water * 3600 / 1e3 * cfg.year_stride * cfg.water_eur_per_m3
            x = np.arange(len(g))
            ax.bar(x, electricity, 0.8, color=SET_COLOR[key], edgecolor='black', linewidth=0.3, zorder=3)
            ax.bar(x, water, 0.8, bottom=electricity, color=STACK['water'], edgecolor='black', linewidth=0.3, zorder=3)
            ax.set_xticks(x, [m.strftime('%b')[0] for m in g.index])
            ax.set_title(key, color=SET_COLOR[key], fontweight='bold')
            ax.grid(False)
            ax.grid(True, axis='y')
        axes[0].set_ylabel("cost [€/month]")
        swatch = lambda color: Patch(facecolor=color, edgecolor='black', linewidth=0.3)
        axes[0].legend(handles=[tuple(swatch(SET_COLOR[k]) for k in cfg.sets), swatch(STACK['water'])], labels=["electricity", "pad water"],
                       handler_map={tuple: HandlerTuple(ndivide=None, pad=0)}, loc='upper left', handlelength=2.4)       # the same legend as annual_costs
        corner(axes[-1], "adaptive fan", 'upper left')
        fig.set_size_inches(PRINT_WIDTH, 2.4)
        fig.tight_layout(w_pad=0.9)
        save_figure(fig, 'monthly_cost', dpi=PRINT_DPI)
    return fig


def water_vs_hours(year, control, control_dry, cfg):
    """What the pad water buys: the hours above the target against the water used per year, one marker per operating mode, the three modes of one
    cooler joined. Both dry modes miss the target for ~1700 h without water; the pad cuts that to ~370 h for 67-105 m3 of water."""
    summaries = annual_summaries(year, control, control_dry, cfg)
    marks = {SCENARIOS[0][0]: 's', SCENARIOS[1][0]: '^', SCENARIOS[2][0]: 'o'}
    n_hours = len(year[(year.set == cfg.set_keys[0]) & (year['mode'] == 'dry')]) * cfg.year_stride
    with print_context():
        fig, ax = plt.subplots(figsize=(PRINT_WIDTH * 0.62, 2.7))
        for key in cfg.sets:
            points = [(summaries[name].loc[key].water_m3, n_hours - summaries[name].loc[key].h_target_met) for name, _ in SCENARIOS]
            ax.plot(*zip(*points), color=SET_COLOR[key], lw=0.6, alpha=0.6, zorder=2)
            for (water, missed), (name, _) in zip(points, SCENARIOS):
                ax.scatter(water, missed, marker=marks[name], s=22, color=SET_COLOR[key], edgecolor='black', linewidth=0.3, zorder=3)
        ax.set_xlabel(r"pad water [m³/a]")
        ax.set_ylabel(f"hours above {cfg.T_c_target:g} °C [h/a]")
        ax.legend(handles=[Line2D([], [], marker=m, ls='none', color='0.5', markeredgecolor='black', markeredgewidth=0.3, ms=4, label=name) for name, m in marks.items()]
                  + [Line2D([], [], color=SET_COLOR[k], lw=PRINT_LW, label=k) for k in cfg.sets], loc='upper right')
        fig.tight_layout()
        save_figure(fig, 'water_vs_hours', dpi=PRINT_DPI)
    return fig


# --- parameter variations ---------------------------------------------------------------------------------------------

MISSED_H = 500                                # a variant that misses the target for more hours than this is dotted
VARIATION_BLOCKS = ((('N_r', r"$N_r$ [–]"), ('N_t', r"$N_t$ [–]"), ('l', r"$l$ [m]")),
                    (('s_l_mm', r"$s_l$ [mm]"), ('s_f_mm', r"$s_f$ [mm]"), ('delta_f_mm', r"$\delta_f$ [mm]")))
COST_LINES = (('cost_fan', "fan", '#C2185B'), ('cost_pump', "pump", TUM_GREY_50), ('cost_water', "water", STACK['water']),
              ('cost_material', "material", STACK['material']), ('cost_total', "total", 'black'))
MISSED_LS = (0, (1, 1.6))


def _variation_curve(ax, d, y, color, lw=PRINT_LW):
    """The variants of one cooler as a line with markers: solid where the target is met in at most MISSED_H hours, dotted (faint markers) beyond; open markers
    outside the correlation ranges (d.outside, only where the cooler's own value is inside them)."""
    x, y, missed, outside = d.value.to_numpy(), np.asarray(y), d.missed.to_numpy(), d.outside.to_numpy()
    for i in range(len(x) - 1):
        ax.plot(x[i:i + 2], y[i:i + 2], color=color, lw=lw, ls=MISSED_LS if missed[i] or missed[i + 1] else '-', zorder=2)
    ax.plot(x[~outside & ~missed], y[~outside & ~missed], ls='none', color=color, marker='o', ms=1.8)
    ax.plot(x[outside & ~missed], y[outside & ~missed], ls='none', marker='o', ms=2.2, mfc='white', mec=color, mew=0.5)
    ax.plot(x[missed], y[missed], ls='none', marker='o', ms=1.4, color=color, alpha=0.45)
    std = d.standard.to_numpy()
    ax.scatter(x[std], y[std], s=26, color=color, edgecolor='black', linewidth=0.7, zorder=5)


def _variation_figure(df, cfg, split, name):
    """Two blocks of 3 columns (parameters), each a pinch row above a cost row. split: the cost row shows the cost of the sensitivity cooler by item
    (fan, pump, water, material, total), else the totals of every cooler."""
    pt = cfg.design_point()
    with print_context():
        fig, grid = plt.subplots(5, 3, figsize=(PRINT_WIDTH, 3.95), gridspec_kw=dict(
            height_ratios=[1, 1, 0.05, 1, 1], hspace=0.34, wspace=0.05, left=0.085, right=0.995, top=0.99, bottom=0.105))
        for ax in grid[2]:
            ax.remove()                                                       # the spacer row between the two blocks
        axes = grid[[0, 1, 3, 4]]
        for b, block in enumerate(VARIATION_BLOCKS):
            for j, (param, xlabel) in enumerate(block):
                ax_pinch, ax_cost = axes[2 * b, j], axes[2 * b + 1, j]
                for key in cfg.sets:
                    d = df[(df.set == key) & (df.param == param)].sort_values('value')
                    _variation_curve(ax_pinch, d, d.pinch, SET_COLOR[key])
                    if not split:
                        _variation_curve(ax_cost, d, d.cost_total, SET_COLOR[key])
                if split:
                    d = df[(df.set == cfg.sensitivity_set) & (df.param == param)].sort_values('value')
                    for column, _, color in COST_LINES:
                        _variation_curve(ax_cost, d, d[column], color, lw=PRINT_LW * (1.4 if column == 'cost_total' else 1))
                ax_pinch.set_ylim(0, 10)
                ax_pinch.tick_params(labelbottom=False)                      # the cost graph below has the same x axis
                ax_cost.set_ylim(0, 1500)
                ax_cost.set_yticks([0, 500, 1000, 1500])
                ax_cost.set_xlabel(xlabel)
            axes[2 * b, 0].set_ylabel(r"$\Delta T_{pinch}$ [K]")
            axes[2 * b + 1, 0].set_ylabel(f"{cfg.sensitivity_set} cost [€/a]" if split else "cost [€/a]")
            for row in (2 * b, 2 * b + 1):
                for ax in axes[row, 1:]:
                    ax.sharey(axes[row, 0])
                    ax.tick_params(labelleft=False)
        kw = dict(fontsize=6.2, borderpad=0.3, labelspacing=0.18, handlelength=1.5, borderaxespad=0.3)
        axes[0, 0].legend(handles=[Line2D([], [], color=SET_COLOR[k], lw=PRINT_LW, label=k) for k in cfg.sets]
                          + [Line2D([], [], color='0.6', marker='o', ls='none', ms=2.6, mfc='white', mec='0.4', mew=0.5, label="outside range"),
                             Line2D([], [], color='0.45', lw=PRINT_LW, ls=MISSED_LS, label=f"> {MISSED_H} h above {cfg.T_c_target:g} °C")],
                          loc='lower left', **kw)
        if split:
            axes[1, 0].legend(handles=[Line2D([], [], color=c, lw=PRINT_LW, label=label) for _, label, c in COST_LINES], loc='upper left', **kw)
        corner(axes[0, 2], rf"{pt['T_a_i']:g} °C, {pt['phi'] * 100:g} %, {cfg.w_fr_design:g} m/s, dry", 'upper right')
        save_figure(fig, name, dpi=PRINT_DPI)
    return fig


def variation_figures(df, cfg):
    """Both variation figures from the variation table (variation.run): the cost split of one cooler and the totals of all."""
    df = variation.add_costs(cfg, df)
    inside = {key: bool((g.wang_ok & g.vdi_ok).all()) for key, g in df[df.standard].groupby('set')}
    df['outside'] = [inside[r.set] and not (r.wang_ok and r.vdi_ok) for r in df.itertuples()]      # leaving the ranges the cooler itself is inside
    df['missed'] = df.h_limited > MISSED_H
    return _variation_figure(df, cfg, True, 'variation_costs'), _variation_figure(df, cfg, False, 'variation_totals')
