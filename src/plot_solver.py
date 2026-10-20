"""
plot_solver.py
===============
Solver figures -- how the solvers behave, as opposed to what the cooler
does (plot_profiles.plot_cooler_results):
  plot_convergence        -- outlet temperatures and iteration errors per
                             iteration (output 'convergence'), per scenario
  plot_caching            -- caching benchmark (output 'caching'): run time per
                             solver and cache mode, where the time goes, the
                             speed-up of the cache, and a presentation sequence
  plot_resolution         -- resolution benchmark (output 'resolution'): pinch error
                             and wall time against the number of elements / segments
  plot_annual             -- annual benchmark (output 'annual'): cumulative time and
                             year totals, pinch deviation from Cell over the year
"""

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch

from src.run_scenario import draw_convergence
from src.solvers import SOLVER_FAMILY
from src import solvers
from src.benchmark_solvers import CACHE_MODES, PART_COLORS, TIME_PARTS, resolution_labels
from src.figures import save_figure, slug, SOLVER
from src.formatting import fmt_compact, fmt_duration
from matplotlib.ticker import NullFormatter

from src.plot_style import PRINT_DASHED, PRINT_DPI, PRINT_LW, PRINT_WIDTH, TUM_GREY_50, in_large_format, lighten_color, print_context, solver_color


@in_large_format
def plot_convergence(result, ops, label: str, settings):
    """Returns the figure (does NOT call plt.show())."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(17, 7.5))
    draw_convergence(ax1, ax2, result, ops, settings.central_omega, settings.cell_omega)
    resolutions = resolution_labels(settings)
    used = [f"{SOLVER_FAMILY[key]}: {resolutions[key]}" for key in result.results if resolutions[key]]
    fig.suptitle(f"Convergence - {label}\n" + (f"({', '.join(used)} per tube and row)" if used else ""),
                 fontsize=15, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    save_figure(fig, SOLVER, slug(f"Convergence - {label}"), f"Convergence - {label}")
    return fig


# --- caching benchmark -------------------------------------------------------------------------------------

SHADE = {'no cache': 0.65, 'cold cache': 0.35, 'warm cache': 0.0}      # lighten factor of the solver colour per cache mode
SLIDE_ORDER = ('cell_2d', 'ntu_field', 'lmtd', 'ntu_table')            # the presentation graph: bars appear in this order
SLIDE_ORDER_WITH_3D = ('cell_3d',) + SLIDE_ORDER                       # the same with Cell (3D) in front (when it was benchmarked)


def _name(key):
    """'NTU (field)' -> two lines for an axis tick."""
    return solvers.SOLVER_LABELS[key].replace(" (", "\n(")


def _grid_along_values(ax):
    """Grid lines only perpendicular to the value (y) axis: the kind of solver on x needs none."""
    ax.grid(False)
    ax.grid(True, axis='y')


def _mode_patches(modes):
    return [Patch(facecolor=lighten_color('#555555', SHADE[m]), edgecolor='black', linewidth=0.4, label=m) for m in modes]


def plot_caching(cache_results):
    """The caching benchmark's figures (output 'caching'), each a graph for itself at print size:
      caching_times    run time of one solve per solver and cache mode (log axis; the solvers span 5 orders of magnitude)
      caching_split    where the time goes: one graph per solver, one bar per cache mode, split into the time parts
      caching_speedup  speed-up of the cold and the warm cache against no cache
    plus the presentation sequences (plot_caching_presentation: caching_presentation_0 ... _4 and, if Cell (3D) was benchmarked,
    caching_presentation_with3d_0 ... _5). Returns the figures."""
    keys = list(cache_results)
    total = lambda key, mode: cache_results[key][mode]['total'] * 1e3
    figures = []

    with print_context(talk=True):
        fig, ax = plt.subplots(figsize=(PRINT_WIDTH, 2.6))
        width = 0.29
        for i, key in enumerate(keys):
            for j, mode in enumerate(CACHE_MODES):
                value, x = total(key, mode), i + (j - 1) * width
                ax.bar(x, value, width * 0.94, color=lighten_color(solver_color(key), SHADE[mode]), edgecolor='black',
                       linewidth=0.4, zorder=3)
                ax.annotate(fmt_compact(value), (x, value), xytext=(0, 2), textcoords='offset points', ha='center', fontsize=5.5)
        ax.set_yscale('log')
        ax.set_ylim(min(total(k, m) for k in keys for m in CACHE_MODES) / 8, max(total(k, m) for k in keys for m in CACHE_MODES) * 60)
        ax.set_xticks(range(len(keys)), [_name(k) for k in keys])
        ax.set_ylabel("wall time per run [ms]")
        _grid_along_values(ax)
        ax.legend(handles=_mode_patches(CACHE_MODES), loc='upper left', ncol=3)
        fig.tight_layout()
        save_figure(fig, SOLVER, 'caching_times', "Caching times", dpi=PRINT_DPI)
        figures.append(fig)

        fig, axes = plt.subplots(1, len(keys), figsize=(PRINT_WIDTH, 2.5), sharey=True, squeeze=False)
        axes = axes[0]
        used = [(part, color) for part, color in zip(TIME_PARTS, PART_COLORS)
                if any(cache_results[k][m]['split'].get(part, 0.0) > 0 for k in keys for m in CACHE_MODES)]
        for ax, key in zip(axes, keys):
            for j, mode in enumerate(CACHE_MODES):
                split = cache_results[key][mode]['split']
                whole = sum(split.values())
                bottom = 0.0
                for part, color in used:
                    share = split.get(part, 0.0) / whole * 100
                    ax.bar(j, share, 0.8, bottom=bottom, color=color, edgecolor='black', linewidth=0.3, zorder=3)
                    bottom += share
            ax.set_xticks(range(len(CACHE_MODES)), [m.replace(" cache", "") for m in CACHE_MODES])
            ax.set_title(solvers.SOLVER_LABELS[key], color=solver_color(key), fontweight='bold')
            ax.set_ylim(0, 100)
            _grid_along_values(ax)
        axes[0].set_ylabel("share of run time [%]")
        fig.legend(handles=[Patch(facecolor=c, edgecolor='black', linewidth=0.3, label=p) for p, c in used],
                   loc='lower center', ncol=len(used), frameon=False, bbox_to_anchor=(0.5, 0.0))
        fig.text(0.5, 0.1, "property-cache mode", ha='center', fontsize=7)
        fig.tight_layout(rect=[0, 0.14, 1, 1], w_pad=0.5)
        save_figure(fig, SOLVER, 'caching_split', "Caching time split", dpi=PRINT_DPI)
        figures.append(fig)

        fig, ax = plt.subplots(figsize=(PRINT_WIDTH, 2.3))
        width = 0.36
        modes = ('cold cache', 'warm cache')
        for i, key in enumerate(keys):
            for j, mode in enumerate(modes):
                factor, x = total(key, 'no cache') / total(key, mode), i + (j - 0.5) * width
                ax.bar(x, factor, width * 0.95, color=lighten_color(solver_color(key), SHADE[mode]), edgecolor='black',
                       linewidth=0.4, zorder=3)
                ax.annotate(f"{fmt_compact(factor)}×", (x, factor), xytext=(0, 2), textcoords='offset points', ha='center', fontsize=6.5)
        ax.axhline(1, color=TUM_GREY_50, lw=0.8, zorder=4)
        ax.set_yscale('log')
        ax.set_ylim(0.5, max(total(k, 'no cache') / total(k, 'warm cache') for k in keys) * 20)
        ax.set_xticks(range(len(keys)), [_name(k) for k in keys])
        ax.set_ylabel("speed-up vs. no cache [–]")
        _grid_along_values(ax)
        ax.legend(handles=_mode_patches(modes), loc='upper left', ncol=2)
        fig.tight_layout()
        save_figure(fig, SOLVER, 'caching_speedup', "Caching speed-up", dpi=PRINT_DPI)
        figures.append(fig)

    figures.append(plot_caching_presentation(cache_results))
    if 'cell_3d' in cache_results:
        figures.append(plot_caching_presentation(cache_results, order=SLIDE_ORDER_WITH_3D, name='caching_presentation_with3d'))
    return figures


def plot_caching_presentation(cache_results, order=SLIDE_ORDER, y_max: float = None, name: str = 'caching_presentation'):
    """A sequence of slides frames on one absolute (linear) scale: the warm-cache run times of Cell (2D), NTU (field),
    LMTD and NTU (table), the bars appearing one by one in this order (caching_presentation_0 ... _n: identical axes, to
    be switched on the slide). The last bars are invisible on this scale -- only their written times show that they
    have appeared. order: the solver keys left to right (those not benchmarked are skipped); y_max: top of the axis
    [ms], default 10 % above the tallest bar; name: file name prefix. Returns the last frame (the others are saved and closed)."""
    order = [key for key in order if key in cache_results]
    if y_max is None:
        y_max = 1.1 * max(cache_results[key]['warm cache']['total'] * 1e3 for key in order)
    slide = {'font.size': 11, 'axes.labelsize': 11, 'xtick.labelsize': 11, 'ytick.labelsize': 10}
    fig = None
    with print_context(slide, talk=True):
        for n in range(len(order) + 1):
            if fig is not None:
                plt.close(fig)
            fig, ax = plt.subplots(figsize=(6.0, 3.6))
            for i, key in enumerate(order[:n]):
                value = cache_results[key]['warm cache']['total'] * 1e3
                ax.bar(i, value, 0.6, color=solver_color(key), edgecolor='black', linewidth=0.4 if value > 0.003 * y_max else 0.0, zorder=3)      # a bar below 0.3 % of the axis: no outline (it would show as a hairline)
                ax.annotate(f"{fmt_compact(value)} ms", (i, value), xytext=(0, 4), textcoords='offset points', ha='center',
                            fontsize=12, fontweight='bold')
            ax.set_xlim(-0.6, len(order) - 0.4)
            ax.set_ylim(0, y_max)
            ax.set_xticks(range(len(order)), [solvers.SOLVER_LABELS[k] for k in order])
            ax.set_ylabel("wall time per run, warm cache [ms]")
            _grid_along_values(ax)
            fig.tight_layout()
            save_figure(fig, SOLVER, f'{name}_{n}', "Caching presentation", dpi=PRINT_DPI)
    return fig


# --- resolution benchmark ----------------------------------------------------------------------------------

def plot_resolution(results, settings):
    """The resolution benchmark (output 'resolution'), two graphs at print size, all axes logarithmic: the pinch error against the
    finest swept resolution (NTU and Cell) and the wall time of every swept solver -- LMTD (no resolution) as a line, the NTU table
    build as a dashed series, the NTU table solve with the table prebuilt only where the table covers the operating point (below,
    the solve falls back to the field solve); the default resolution(s) of the settings dotted in both. Returns the figure."""
    n_ref = max(n for n, *_ in next(rows for key, rows in results.items() if key in ('ntu_field', 'cell_2d', 'cell_3d') and rows))
    default_n = sorted({settings.ntu_n_elements, settings.cell_n_segments})

    def error(key):
        rows = results[key]
        ref = next(p for n, _, _, p in rows if n == n_ref)
        return [n for n, *_ in rows if n != n_ref], [abs(p - ref) for n, _, _, p in rows if n != n_ref]

    def style(ax):
        ticks = [2, 5, 10, 20, 50, 100]
        ax.set_xscale('log')                                      # first: it resets the tick locator and formatter
        ax.set_yscale('log')
        ax.set_xticks(ticks, [str(t) for t in ticks])
        ax.xaxis.set_minor_formatter(NullFormatter())
        ax.set_xlabel("$N$ (elements / segments)")
        for n in default_n:
            ax.axvline(n, color=TUM_GREY_50, lw=0.9, ls=':', zorder=0, label=f"default $N$ = {n}")

    ntu_label = f"NTU ({results['ntu_variants'].replace('/', ' / ')})"
    cell_key = 'cell_2d' if results.get('cell_2d') else 'cell_3d'
    cell_label = "Cell (2D / 3D)" if results.get('cell_2d') and results.get('cell_3d') else solvers.SOLVER_LABELS[cell_key]
    with print_context(talk=True):
        fig, (ax_e, ax_t) = plt.subplots(1, 2, figsize=(PRINT_WIDTH, 2.5))
        for key, label in (('ntu_field', ntu_label), (cell_key, cell_label)):
            if results.get(key):
                n, e = error(key)
                ax_e.plot(n, e, marker='o', ms=2.5, lw=PRINT_LW, color=solver_color(key), label=label)
        style(ax_e)
        ax_e.set_ylabel(r"$\Delta T_{pinch}$ error [K]")
        ax_e.legend(loc='upper right')

        series = []                                               # (n, time [ms]) of everything drawn, for the axis limits
        for key in ('cell_3d', 'cell_2d', 'ntu_field'):
            if results.get(key):
                n, t = zip(*[(n, t * 1e3) for n, t, *_ in results[key]])
                ax_t.plot(n, t, marker='o', ms=2.5, lw=PRINT_LW, color=solver_color(key), label=solvers.SOLVER_LABELS[key])
                series += t
        if 'lmtd' in results:
            ax_t.axhline(results['lmtd'][0] * 1e3, color=solver_color('lmtd'), lw=PRINT_LW, ls=PRINT_DASHED, label="LMTD")
            series.append(results['lmtd'][0] * 1e3)
        covered = [(n, t * 1e3) for n, t, *_ in results.get('ntu_table', []) if n in results['table_covers']]
        if covered:
            ax_t.plot(*zip(*covered), marker='o', ms=2.5, lw=PRINT_LW, color=solver_color('ntu_table'), label="NTU (table, prebuilt)")
            series += [t for _, t in covered]
        if results.get('table_build'):
            n, t = zip(*[(n, t * 1e3) for n, t in results['table_build']])
            ax_t.plot(n, t, marker='s', ms=2.5, lw=PRINT_LW, ls=PRINT_DASHED, color=solver_color('ntu_table'), label="NTU table build")
            series += t
        style(ax_t)
        ax_t.set_ylim(min(series) / 3, max(series) * 60)          # room above the curves for the legend
        ax_t.set_ylabel("wall time [ms]")
        ax_t.legend(loc='upper left', ncol=2)
        fig.tight_layout(w_pad=1.2)
        save_figure(fig, SOLVER, 'resolution_performance', "Resolution performance", dpi=PRINT_DPI)
    return fig


# --- annual benchmark --------------------------------------------------------------------------------------

def plot_annual(results):
    """The annual benchmark's figures (output 'annual'), print size, each a graph pair / graph for itself:
      annual_times      cumulative wall time over the operating points solved (the axes follow the non-Cell solvers: Cell runs
                        out of the image, which shows how steep it is) | the year total per solver (log axis; the NTU table with
                        and without its build; what was extrapolated hatched)
      annual_deviation  pinch deviation of the other solvers from the most detailed Cell run against the ambient temperature
                        (only if a Cell variant ran)
    Returns the figures."""
    sol = results['solvers']
    keys = {label: key for key, label in solvers.SOLVER_LABELS.items()}
    color = lambda name: solver_color(keys[name])
    n_hours = len(results['T_air'])

    def cumulative(r):
        times = r['times'][r['hours']][:r['n_done']]
        return np.arange(1, r['n_done'] + 1), np.cumsum(times)

    figures_ = []
    with print_context(talk=True):
        fig, (ax_cum, ax_tot) = plt.subplots(1, 2, figsize=(PRINT_WIDTH, 2.7), gridspec_kw={'width_ratios': [1, 1.15]})
        steady = [n for n in sol if not n.startswith('Cell')] or list(sol)
        y_top = 0.0
        for name in steady + [n for n in sol if n.startswith('Cell')]:
            x, y = cumulative(sol[name])
            offset = results['t_build'] if name == 'NTU (table)' else 0.0
            if name in steady:
                y_top = max(y_top, offset + y[-1])
            ax_cum.plot(np.concatenate([[0], x]), offset + np.concatenate([[0.0], y]), lw=PRINT_LW, color=color(name),
                        ls=PRINT_DASHED if offset else '-', label=name + (" incl. build" if offset else ""))
        ax_cum.set_xlim(0, n_hours)
        ax_cum.set_ylim(0, 1.05 * y_top)
        ax_cum.set_xlabel("hours solved")
        ax_cum.set_ylabel("cumulative time [s]")
        ax_cum.legend(loc='upper left')

        rows = []
        for name, r in sol.items():
            if name == 'NTU (table)':
                rows += [("NTU (table)\nprebuilt", name, r, 0.0), ("NTU (table)\nincl. build", name, r, results['t_build'])]
            else:
                rows.append((name.replace(" (", "\n("), name, r, 0.0))
        for i, (label, name, r, offset) in enumerate(rows):
            total = r['total'] + offset
            ax_tot.barh(i, total, 0.66, color=lighten_color(color(name), 0.45 if offset else 0.0), edgecolor='black',
                        linewidth=0.4, hatch=None if r['complete'] else '////', zorder=3)
            ax_tot.annotate((" " if r['complete'] else " ≈") + fmt_duration(total), (total, i), va='center', fontsize=6.5)
        ax_tot.set_xscale('log')
        ax_tot.set_xlim(min(r['total'] + o for _, _, r, o in rows) / 5, max(r['total'] + o for _, _, r, o in rows) * 14)
        ax_tot.set_yticks(range(len(rows)), [row[0] for row in rows])
        ax_tot.invert_yaxis()
        ax_tot.set_xlabel(f"time for the year ({n_hours} hours) [s]")
        ax_tot.grid(False)
        ax_tot.grid(True, axis='x')
        if any(not r['complete'] for _, _, r, _ in rows):
            ax_tot.legend(handles=[Patch(facecolor='white', edgecolor='black', hatch='////', linewidth=0.4, label="extrapolated")],
                          loc='upper right')
        fig.tight_layout(w_pad=1.0)
        save_figure(fig, SOLVER, 'annual_times', "Annual times", dpi=PRINT_DPI)
        figures_.append(fig)

        reference = next((n for n in ('Cell (3D)', 'Cell (2D)') if n in sol), None)
        ntu = 'NTU (table)' if 'NTU (table)' in sol else 'NTU (field)'                 # the two give the same results: one is drawn
        others = [n for n in ('LMTD', ntu) if n in sol]
        if reference and others:
            fig, ax = plt.subplots(figsize=(PRINT_WIDTH * 0.62, 2.6))
            ref = sol[reference]['T_c_o']
            for name in others:
                T_c_o = sol[name]['T_c_o']
                ok = ~np.isnan(T_c_o) & ~np.isnan(ref)
                ax.plot(results['T_air'][ok], T_c_o[ok] - ref[ok], '.', ms=1.6, color=color(name), alpha=0.6, label=name, rasterized=True)
            ax.axhline(0.0, color=color(reference), lw=PRINT_LW, label=f"{reference} (reference)")
            ax.set_xlabel(r"$T_{a,i}$ [°C]")
            ax.set_ylabel(r"$\Delta T_{pinch}$ deviation [K]")
            ax.legend(loc='best', markerscale=5)
            fig.tight_layout()
            save_figure(fig, SOLVER, 'annual_deviation', "Annual deviation", dpi=PRINT_DPI)
            figures_.append(fig)
    return figures_
