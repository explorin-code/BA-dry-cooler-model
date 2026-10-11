"""
benchmark_solvers.py
=====================
Solver benchmarks (outputs 'caching' / 'resolution'; main.py --outputs caching resolution): the computational cost of
the selected solvers (SOLVERS: LMTD, NTU field / table, Cell 2D / 3D) on one operating point. Runs its own fresh,
silent solves -- never the normal run's, whose property cache is already warm and whose terminal output itself
costs time. Two parts:

1. Caching (run_cache_benchmark): each solver at its configured resolution in three property-cache modes -- no
   cache (every lookup calls CoolProp), cold cache (cache on, starting empty: a normal single run) and warm
   cache (pure algorithm cost) -- with the time split into lookup (CoolProp), lookup (cached), k evaluation and the
   rest of the solver, plus work counts. NTU (table) additionally pays the table build in the no-cache and
   cold-cache runs (a first run needs the table); in the warm-cache runs the table is assumed to be built already.
   Drawn into the solver-performance figure (plot_solver.py).
2. Resolution (run_resolution_benchmark): NTU (n_elements) and Cell (n_segments) over BENCHMARK_RESOLUTIONS, cache on
   and emptied before every run (= a normal run): wall time, outer iterations, time per iteration and the pinch point. NTU runs
   as the field solver only (the table gives the same results; its build is part of the caching benchmark). LMTD has no resolution.

The measurement uses the hooks of src/profiling.py (`@timed` on the property lookups, CoolProp calls and the k
evaluation, active only inside `profiling.record()`) and fluid_properties.caching() / clear_cache() for the cache modes -- no
function is patched. The hooks' own overhead (~2-4 % of Cell's time, also outside benchmarks) is included in the totals.
"""


from contextlib import contextmanager
import logging
import statistics
import time

import dataclasses

import matplotlib.pyplot as plt

import src.fluid_properties as fluid_properties
import src.solvers as solvers
from src.solvers import ntu_table
from src import profiling
from src.figures import save_figure
from src.formatting import fmt_ms, machine_info
from src.quiet import silence_all
from src.plot_style import (COLOR_LMTD, TUM_ACCENT_BLUE, TUM_BLUE, TUM_BLUE_540, TUM_GREY_20, TUM_IVORY,
                            family_color, solver_color)

logger = logging.getLogger(__name__)

CACHE_MODES = ('no cache', 'cold cache', 'warm cache')
TIME_PARTS = ('lookup (CoolProp)', 'lookup (cached)', 'k evaluation', 'rest of solver', 'table build')
PART_COLORS = (TUM_BLUE_540, TUM_BLUE, TUM_ACCENT_BLUE, TUM_GREY_20, TUM_IVORY)      # TUM blues, grey, ivory


# =============================================================================
# Measurement helpers
# =============================================================================

@contextmanager
def _cache_mode(mode: str, fn):
    """Puts the property cache into `mode` before a measured run of fn(): 'no cache' switches the cache off,
    'cold cache' starts it empty, 'warm cache' fills it with one unmeasured run of fn()."""
    if mode == 'no cache':
        with fluid_properties.caching(False):
            yield
        return
    fluid_properties.clear_cache()
    if mode == 'warm cache':
        fn()
    yield


def _measured_run(fn):
    """One counted run of fn(): total wall time, time split and counts (src/profiling.py hooks)."""
    with profiling.record() as rec:
        t_start = time.perf_counter()
        output = fn()
        total = time.perf_counter() - t_start
    t_lookups = rec.time['lookup']
    return {
        'total': total,
        'split': {
            'lookup (CoolProp)': rec.time['coolprop'],
            'lookup (cached)': max(t_lookups - rec.time['coolprop'], 0.0),
            'k evaluation': rec.time['k'],
            'rest of solver': max(total - t_lookups - rec.time['k'], 0.0),
        },
        'iterations': len(output.history_hot),
        'k evaluations': rec.calls['k'],
        'property lookups': rec.calls['lookup'],
        'CoolProp calls': rec.calls['coolprop'],
        'Q [kW]': output.Q_dot / 1000,
    }


@contextmanager
def preserved_ntu_tables():
    """The in-memory NTU tables as they were, restored afterwards (a benchmark installs freshly built ones)."""
    saved = dict(ntu_table.tables)
    try:
        yield
    finally:
        ntu_table.tables.clear()
        ntu_table.tables.update(saved)


def build_fresh_table(N_r: int, n_elements: int) -> float:
    """Builds the NTU table P(NTU, R) from scratch (never from disk), installs it for this run and returns the
    build time [s]."""
    t_start = time.perf_counter()
    ntu_table.install_table(N_r, n_elements, ntu_table.build_ntu_table(N_r, n_elements))
    return time.perf_counter() - t_start


def _solver_runs(ops, geo, settings, resolution=None):
    """One zero-argument callable per solver key; NTU elements / Cell segments at `resolution` (None -> the
    settings' own)."""
    if resolution is not None:
        settings = dataclasses.replace(settings, ntu_n_elements=resolution, cell_n_segments=resolution)
    return {key: (lambda key=key: solvers.run_solver(key, ops, geo, settings, profile=False))
            for key in solvers.SOLVER_LABELS}


# =============================================================================
# 1. Caching
# =============================================================================

def run_cache_benchmark(ops, geo, settings, modes):
    """results[solver key][cache mode] = _measured_run dict (median of benchmark_repeats runs by total time; 'no
    cache' runs once -- it is slow and its spread is small relative to its size). NTU (table): the build time is
    added to 'no cache' and 'cold cache' (split part 'table build'), not to 'warm cache'."""
    t_benchmark = time.perf_counter()
    results = {}
    with silence_all(), preserved_ntu_tables():
        for key, fn in _solver_runs(ops, geo, settings).items():
            if key not in modes.solvers:
                continue
            if key == 'ntu_table':
                t_build = build_fresh_table(geo.N_r, settings.ntu_n_elements)       # installed: runs below find it
            results[key] = {}
            for mode in CACHE_MODES:
                runs = []
                for _ in range(1 if mode == 'no cache' else modes.benchmark_repeats):
                    with _cache_mode(mode, fn):
                        runs.append(_measured_run(fn))
                runs.sort(key=lambda r: r['total'])
                run = runs[len(runs) // 2]
                if key == 'ntu_table':
                    build = 0.0 if mode == 'warm cache' else t_build
                    run['total'] += build
                    run['split']['table build'] = build
                results[key][mode] = run
    fluid_properties.clear_cache()

    logger.info("[Benchmark: caching] %s (took %.1f s)", machine_info(), time.perf_counter() - t_benchmark)
    logger.info("  %-12s %-11s %10s %6s %9s %10s %9s %8s",
                "", "", "time", "iter", "k evals", "lookups", "CoolProp", "Q [kW]")
    for key, per_mode in results.items():
        for mode, r in per_mode.items():
            logger.info("  %-12s %-11s %7.1f ms %6d %9d %10d %9d %8.2f", solvers.SOLVER_LABELS[key], mode,
                        r['total'] * 1e3, r['iterations'], r['k evaluations'], r['property lookups'],
                        r['CoolProp calls'], r['Q [kW]'])
    return results


def resolution_labels(settings) -> dict:
    """Resolution each solver runs at (LMTD has none), by solver key."""
    elements, segments = f"{settings.ntu_n_elements} elements", f"{settings.cell_n_segments} segments"
    return {'lmtd': "", 'ntu_field': elements, 'ntu_table': elements, 'cell_2d': segments, 'cell_3d': segments}


def draw_cache_panel(ax, key: str, cache_results, resolution: str = "", show_modes: bool = True, log: bool = False):
    """One solver: a horizontal bar per cache mode, length = wall time on this solver's own ms axis (the solvers
    differ by ~4 orders of magnitude, so they can't share one), split into the time parts. log: logarithmic axis,
    one plain bar in the solver's colour per cache mode -- the stacked parts would be misleading there (a part's
    share of the length is not its share of the time)."""
    per_mode, name = cache_results[key], solvers.SOLVER_LABELS[key]
    totals = [r['total'] * 1e3 for r in per_mode.values()]
    x_min = min(totals) / 30                                              # log axis: where the bars start
    for row, mode in enumerate(CACHE_MODES):
        left = 0.0
        if log:
            left = per_mode[mode]['total'] * 1e3
            ax.barh(row, left - x_min, left=x_min, color=solver_color(key), edgecolor='black', linewidth=0.5)
        else:
            for part, color in zip(TIME_PARTS, PART_COLORS):
                width = per_mode[mode]['split'].get(part, 0.0) * 1e3
                ax.barh(row, width, left=left, color=color, edgecolor='black', linewidth=0.5)
                left += width
        ax.annotate(f" {fmt_ms(per_mode[mode]['total'] * 1e3)}", (left, row), va='center', fontsize=10)
    ax.set_yticks(range(len(CACHE_MODES)), CACHE_MODES if show_modes else [])
    ax.invert_yaxis()
    if log:
        ax.set_xscale('log')
        ax.set_xlim(x_min, max(totals) * 6)                               # room for the labels
    else:
        ax.set_xlim(0, max(totals) * 1.35)
    ax.set_xlabel("Wall time [ms]" + (" (log)" if log else ""))
    counts = per_mode['cold cache']
    ax.set_title(f"{name}{f' ({resolution})' if resolution else ''} — caching\n"
                 f"{counts['iterations']} it, {counts['k evaluations']} k evals, "
                 f"{counts['property lookups']} lookups", fontsize=12, color=family_color(key),
                 fontweight="bold")


MODE_TITLES = {
    'no cache': "No cache: every lookup calls CoolProp",
    'cold cache': "Cold cache: single normal run",
    'warm cache': "Warm cache: repeated runs (default case)",
}


def _set_time_axis(ax, values, log: bool):
    """The wall-time axis of a bar panel: linear from 0, or logarithmic around the values."""
    if log:
        ax.set_yscale('log')
        ax.set_ylim(min(values) / 3, max(values) * 4)
    else:
        ax.set_ylim(0, max(values) * 1.15)
    ax.set_ylabel("Wall time per run [ms]" + (" (log)" if log else ""))


def draw_mode_comparison(ax, cache_results, mode: str, resolutions=None, log: bool = False):
    """The solvers in one cache mode, on ONE linear axis. Values are
    labelled, since LMTD's bar is often too small to see."""
    keys = list(cache_results)
    names = [solvers.SOLVER_LABELS[k] for k in keys]
    values = [cache_results[k][mode]['total'] * 1e3 for k in keys]
    bars = ax.bar(names, values, color=[solver_color(k) for k in keys], edgecolor='black', linewidth=0.6)
    for bar, value in zip(bars, values):
        ax.annotate(fmt_ms(value), (bar.get_x() + bar.get_width() / 2, value),
                    textcoords="offset points", xytext=(0, 4), ha='center', fontsize=11)
    if resolutions:
        ax.set_xticks(range(len(names)), [f"{n}\n({resolutions[k]})" if resolutions[k] else n for n, k in zip(names, keys)])
        ax.tick_params(axis='x', labelsize=9)
    _set_time_axis(ax, values, log)
    ax.set_title(MODE_TITLES[mode])


# =============================================================================
# 2. Resolution
# =============================================================================

# series of the resolution benchmark: key -> (label, colour, line style)
RESOLUTION_SERIES = {
    'cell_3d': ("Cell (3D)", solver_color('cell_3d'), '-'),
    'cell_2d': ("Cell (2D)", solver_color('cell_2d'), '-'),
    'ntu_field': ("NTU (field)", solver_color('ntu_field'), '-'),      # also stands for NTU (table): same results
}


def run_resolution_benchmark(ops, geo, settings, modes):
    """results[series key] = list of (n, time [s], iterations, pinch [K]) over modes.benchmark_resolutions for the
    selected NTU / Cell variants (series: RESOLUTION_SERIES), results['lmtd'] = (time [s], iterations, pinch [K]).
    NTU runs once, as the field solver: the table gives the same results (its P(NTU, R) lookup is the field solved
    once and interpolated), and its build cost is the caching benchmark's business; results['ntu_variants'] names which
    of the two were selected. Cache on and emptied before every run (= a normal run)."""
    t_benchmark = time.perf_counter()
    swept = [k for k in ('ntu_field', 'cell_2d', 'cell_3d')
             if k in modes.solvers or (k == 'ntu_field' and 'ntu_table' in modes.solvers)]
    results = {key: [] for key in swept}
    results['ntu_variants'] = '/'.join(v for v in ('field', 'table') if f'ntu_{v}' in modes.solvers)
    with silence_all():
        for n in modes.benchmark_resolutions:
            runs = _solver_runs(ops, geo, settings, n)
            for key in swept:
                times = []
                for _ in range(modes.benchmark_sweep_repeats):
                    fluid_properties.clear_cache()
                    t_start = time.perf_counter()
                    output = runs[key]()
                    times.append(time.perf_counter() - t_start)
                t = statistics.median(times)
                results[key].append((n, t, len(output.history_hot), output.T_c_o - ops.T_a_i))
        if 'lmtd' in modes.solvers:
            fluid_properties.clear_cache()
            t_start = time.perf_counter()
            output = _solver_runs(ops, geo, settings)['lmtd']()
            results['lmtd'] = (time.perf_counter() - t_start, len(output.history_hot), output.T_c_o - ops.T_a_i)
    fluid_properties.clear_cache()

    logger.info("[Benchmark: resolution] %s (took %.1f s)", machine_info(), time.perf_counter() - t_benchmark)
    for key, (name, _, _) in RESOLUTION_SERIES.items():
        if results.get(key):
            logger.info("  %-28s %s", name, "   ".join(f"n={n}: {t * 1e3:.0f} ms/{it} it" for n, t, it, _ in results[key]))
    return results


def _pinch_inset(ax, results):
    """The pinch panel runs from 0 K (the true scale); the solvers differ by ~1 % of it, so an inset shows the details."""
    values = [p for key in RESOLUTION_SERIES for *_, p in results.get(key, [])]
    values += [results['lmtd'][2]] if 'lmtd' in results else []
    ax.set_ylim(0, max(values) * 1.15)
    inset = ax.inset_axes([0.5, 0.2, 0.47, 0.42])
    for key, (name, color, line) in RESOLUTION_SERIES.items():
        if results.get(key):
            n, _, _, pinch = zip(*results[key])
            inset.plot(n, pinch, marker='o', ms=3, color=color, linestyle=line)
    if 'lmtd' in results:
        inset.axhline(results['lmtd'][2], color=COLOR_LMTD, linestyle='--')
    inset.set_title("zoom", fontsize=9)
    inset.tick_params(labelsize=8)
    inset.set_xlim(left=0)


def plot_resolution(results):
    """2x2, all linear, every swept solver in every panel: wall time, time per outer iteration, iterations to
    convergence, pinch point (full scale from 0 K, with an inset zoom -- the solvers differ by ~1 % of it)."""
    fig, ((ax_t, ax_ti), (ax_it, ax_P)) = plt.subplots(2, 2, figsize=(17, 11))

    iterations = []
    for key, (name, color, line) in RESOLUTION_SERIES.items():
        if not results.get(key):
            continue
        n, t, it, pinch = zip(*results[key])
        style = dict(marker='o', color=color, linestyle=line)
        ax_t.plot(n, [x * 1e3 for x in t], label=name, **style)
        ax_ti.plot(n, [x * 1e3 / i for x, i in zip(t, it)], label=name, **style)
        # the result panels name both NTU variants (the table gives the same values; only the field is drawn)
        name_result = f"NTU ({results['ntu_variants']})" if key == 'ntu_field' else name
        ax_it.plot(n, it, label=name_result, **style)
        ax_P.plot(n, pinch, label=name_result, **style)
        iterations += list(it)

    # LMTD has no resolution: a dashed reference line in every panel.
    if 'lmtd' in results:
        lmtd_t, lmtd_it, lmtd_P = results['lmtd']
        iterations.append(lmtd_it)
        for ax, value in ((ax_t, lmtd_t * 1e3), (ax_ti, lmtd_t * 1e3 / lmtd_it), (ax_it, lmtd_it), (ax_P, lmtd_P)):
            ax.axhline(value, color=COLOR_LMTD, linestyle='--', label="LMTD (no resolution)")

    for ax, title, ylabel in ((ax_t, "Wall time", "Wall time [ms]"),
                              (ax_ti, "Time per outer iteration", "Time per iteration [ms]"),
                              (ax_it, "Outer iterations to convergence", "Iterations"),
                              (ax_P, r"Result: pinch point $\Delta T_{pinch} = T_{c,o} - T_{a,i}$", r"$\Delta T_{pinch}$ [K]")):
        ax.set_title(title, fontsize=13)
        ax.set_ylabel(ylabel)
        ax.set_xlabel("Cell segments / NTU elements per tube and row")
        ax.set_xlim(left=0)
        ax.legend(fontsize=10)
    for ax in (ax_t, ax_ti):
        ax.set_ylim(bottom=0)
    _pinch_inset(ax_P, results)
    ax_P.legend(fontsize=10, loc='lower left')
    # Headroom above the (flat) top curve, so it doesn't sit on the frame.
    ax_it.set_ylim(0, max(iterations) * 1.3)

    fig.suptitle(f"Resolution performance\n"
                 f"(cache on, emptied before every run = normal single run; {machine_info()})",
                 fontsize=14, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    save_figure(fig, 'solver', 'resolution_performance', "Resolution performance")
    return fig
