"""
benchmark_solvers.py
=====================
Benchmark mode (BENCHMARK_MODE / main.py --benchmark): measures the
computational cost of LMTD, NTU and Cell on one operating point. Runs its
own fresh, silent solves -- never the normal run's, whose property cache is
already warm and whose terminal output itself costs time. Two parts:

1. Caching (run_cache_benchmark): each solver at its configured resolution
   in three property-cache modes -- no cache (every lookup calls CoolProp),
   cold cache (cache on, starting empty: a normal single run) and warm
   cache (pure algorithm cost) -- with the time split into lookup
   (CoolProp), lookup (cached), k evaluation and the rest of the solver,
   plus work counts. Drawn by draw_mode_comparison()/draw_cache_panel() into the
   solver-performance figure (plot_solver.py). Toggle: BENCHMARK_MODE.
2. Resolution (run_resolution_benchmark), toggle RESOLUTION_MODE: NTU (n_elements) and Cell
   (n_segments) over BENCHMARK_RESOLUTIONS, cache on and emptied before
   every run (= a normal run): wall time, outer iterations, time per
   iteration and Q. Own figure (plot_resolution). LMTD has no resolution.

No solver code is changed: counting/timing wrappers are patched into the
solvers module and CoolProp only for the duration of a measurement. Their
own overhead (~2-3 % of Cell's time) is included in the totals.

Note: Cell solves all n_tubes tubes separately although with uniform inlet
air they behave nearly identically (NTU solves one circuit) -- up to
~n_tubes x more work than necessary. Kept as is; see CLAUDE.md.
"""

from contextlib import contextmanager, redirect_stdout
import io
import logging
import platform
import statistics
import time
import warnings

import CoolProp
import CoolProp.CoolProp as CP
import matplotlib.pyplot as plt

import src.fluid_properties as fluid_properties
import src.solvers as solvers
from src.plot_style import COLOR_LMTD, COLOR_NTU, COLOR_CELL

logger = logging.getLogger(__name__)

SOLVER_COLORS = {'LMTD': COLOR_LMTD, 'NTU': COLOR_NTU, 'Cell': COLOR_CELL}
CACHE_MODES = ('no cache', 'cold cache', 'warm cache')
TIME_PARTS = ('lookup (CoolProp)', 'lookup (cached)', 'k evaluation', 'rest of solver')
PART_COLORS = ('#b2182b', '#ef8a62', '#67a9cf', '#d9d9d9')


def machine_info() -> str:
    return (f"{platform.processor() or platform.machine()}, Python {platform.python_version()}, "
            f"CoolProp {CoolProp.__version__}")


# =============================================================================
# Measurement helpers
# =============================================================================

@contextmanager
def _quiet():
    """Silences logging, warnings and stdout -- output would distort timings."""
    logging.disable(logging.CRITICAL)
    try:
        with warnings.catch_warnings(), redirect_stdout(io.StringIO()):
            warnings.simplefilter("ignore")
            yield
    finally:
        logging.disable(logging.NOTSET)


def _clear_property_cache():
    fluid_properties._get_fluid_properties_cached.cache_clear()
    fluid_properties._get_air_properties_cached.cache_clear()


@contextmanager
def _cache_mode(mode: str, fn):
    """Puts the property cache into `mode` before a measured run of fn():
    'no cache' bypasses the lru_cache, 'cold cache' empties it, 'warm cache'
    fills it with one unmeasured run of fn()."""
    if mode == 'no cache':
        cached_c = fluid_properties._get_fluid_properties_cached
        cached_a = fluid_properties._get_air_properties_cached
        fluid_properties._get_fluid_properties_cached = cached_c.__wrapped__
        fluid_properties._get_air_properties_cached = cached_a.__wrapped__
        try:
            yield
        finally:
            fluid_properties._get_fluid_properties_cached = cached_c
            fluid_properties._get_air_properties_cached = cached_a
    else:
        _clear_property_cache()
        if mode == 'warm cache':
            fn()
        yield


class _Counter:
    """Wraps a function, counting calls and accumulating their wall time."""
    def __init__(self, fn):
        self.fn = fn
        self.calls = 0
        self.time = 0.0

    def __call__(self, *args, **kwargs):
        t_start = time.perf_counter()
        try:
            return self.fn(*args, **kwargs)
        finally:
            self.time += time.perf_counter() - t_start
            self.calls += 1


@contextmanager
def _counting():
    """Temporarily wraps what the solvers call by name (k, property lookups)
    and CoolProp's own entry points, so time and calls can be attributed."""
    targets = [(solvers, 'calc_overall_k'), (solvers, 'get_fluid_properties'),
               (solvers, 'get_air_properties'), (CP, 'PropsSI'), (CP, 'HAPropsSI')]
    originals = [(module, name, getattr(module, name)) for module, name in targets]
    counters = {}
    for module, name, fn in originals:
        counters[name] = _Counter(fn)
        setattr(module, name, counters[name])
    try:
        yield counters
    finally:
        for module, name, fn in originals:
            setattr(module, name, fn)


def _measured_run(fn):
    """One counted run of fn(): total wall time, time split and counts."""
    with _counting() as c:
        t_start = time.perf_counter()
        output = fn()
        total = time.perf_counter() - t_start
    t_coolprop = c['PropsSI'].time + c['HAPropsSI'].time
    t_lookups = c['get_fluid_properties'].time + c['get_air_properties'].time
    t_k = c['calc_overall_k'].time
    return {
        'total': total,
        'split': {
            'lookup (CoolProp)': t_coolprop,
            'lookup (cached)': max(t_lookups - t_coolprop, 0.0),
            'k evaluation': t_k,
            'rest of solver': max(total - t_lookups - t_k, 0.0),
        },
        'iterations': len(output[4]),
        'k evaluations': c['calc_overall_k'].calls,
        'property lookups': c['get_fluid_properties'].calls + c['get_air_properties'].calls,
        'CoolProp calls': c['PropsSI'].calls + c['HAPropsSI'].calls,
        'Q [kW]': output[1] / 1000,
    }


def _solver_runs(ops, geo, settings, resolution=None):
    """One zero-argument callable per solver; NTU/Cell at `resolution`
    (None -> settings.ntu_n_elements / settings.cell_n_segments)."""
    return {
        'LMTD': lambda: solvers.solve_it_LMTD(ops=ops, geo=geo, settings=settings),
        'NTU': lambda: solvers.solve_it_NTU(ops=ops, geo=geo, settings=settings, n_elements=resolution),
        'Cell': lambda: solvers.solve_it_cell(n_segments=resolution, ops=ops, geo=geo, settings=settings),
    }


# =============================================================================
# 1. Caching
# =============================================================================

def run_cache_benchmark(ops, geo, settings, modes):
    """results[solver][cache mode] = _measured_run dict (median of
    benchmark_repeats runs by total time; 'no cache' runs once -- it is slow
    and its spread is small relative to its size)."""
    t_benchmark = time.perf_counter()
    results = {}
    with _quiet():
        for name, fn in _solver_runs(ops, geo, settings).items():
            results[name] = {}
            for mode in CACHE_MODES:
                runs = []
                for _ in range(1 if mode == 'no cache' else modes.benchmark_repeats):
                    with _cache_mode(mode, fn):
                        runs.append(_measured_run(fn))
                runs.sort(key=lambda r: r['total'])
                results[name][mode] = runs[len(runs) // 2]
    _clear_property_cache()

    logger.info("[Benchmark: caching] %s (took %.1f s)", machine_info(), time.perf_counter() - t_benchmark)
    logger.info("  %-5s %-11s %10s %6s %9s %10s %9s %8s",
                "", "", "time", "iter", "k evals", "lookups", "CoolProp", "Q [kW]")
    for name, per_mode in results.items():
        for mode, r in per_mode.items():
            logger.info("  %-5s %-11s %7.1f ms %6d %9d %10d %9d %8.2f", name, mode, r['total'] * 1e3,
                        r['iterations'], r['k evaluations'], r['property lookups'], r['CoolProp calls'],
                        r['Q [kW]'])
    return results


def resolution_labels(settings) -> dict:
    """Resolution each solver ran at in the caching benchmark (LMTD has none)."""
    return {'LMTD': "", 'NTU': f"{settings.ntu_n_elements} elements",
            'Cell': f"{settings.cell_n_segments} segments"}


def draw_cache_panel(ax, name: str, cache_results, show_legend: bool = False, resolution: str = ""):
    """One solver: a horizontal bar per cache mode, length = wall time on
    this solver's own LINEAR ms axis (the solvers differ by ~4 orders of
    magnitude, so they can't share one), split into the time parts."""
    per_mode = cache_results[name]
    for row, mode in enumerate(CACHE_MODES):
        left = 0.0
        for part, color in zip(TIME_PARTS, PART_COLORS):
            width = per_mode[mode]['split'][part] * 1e3
            ax.barh(row, width, left=left, color=color, edgecolor='black', linewidth=0.5,
                    label=part if row == 0 else None)
            left += width
        ax.annotate(f" {per_mode[mode]['total'] * 1e3:.1f} ms", (left, row), va='center', fontsize=10)
    ax.set_yticks(range(len(CACHE_MODES)), CACHE_MODES)
    ax.invert_yaxis()
    ax.set_xlim(0, max(r['total'] for r in per_mode.values()) * 1e3 * 1.35)   # room for the labels
    ax.set_xlabel("Wall time [ms]")
    counts = per_mode['cold cache']
    ax.set_title(f"{name}{f' ({resolution})' if resolution else ''} — caching\n"
                 f"{counts['iterations']} it, {counts['k evaluations']} k evals, "
                 f"{counts['property lookups']} lookups", fontsize=12, color=SOLVER_COLORS[name],
                 fontweight="bold")
    if show_legend:
        ax.legend(fontsize=9, loc='lower right', framealpha=1.0)


MODE_TITLES = {
    'no cache': "No cache: every lookup calls CoolProp",
    'cold cache': "Cold cache: single normal run",
    'warm cache': "Warm cache: repeated runs (default case)",
}


def draw_mode_comparison(ax, cache_results, mode: str, resolutions=None):
    """All three solvers in one cache mode, on ONE linear axis. Values are
    labelled, since LMTD's bar is often too small to see."""
    names = list(cache_results)
    values = [cache_results[n][mode]['total'] * 1e3 for n in names]
    bars = ax.bar(names, values, color=[SOLVER_COLORS[n] for n in names], edgecolor='black', linewidth=0.6)
    for bar, value in zip(bars, values):
        ax.annotate(f"{value:.1f} ms", (bar.get_x() + bar.get_width() / 2, value),
                    textcoords="offset points", xytext=(0, 4), ha='center', fontsize=11)
    if resolutions:
        ax.set_xticks(range(len(names)), [f"{n}\n({resolutions[n]})" if resolutions[n] else n for n in names])
        ax.tick_params(axis='x', labelsize=10)
    ax.set_ylim(0, max(values) * 1.15)
    ax.set_ylabel("Wall time per run [ms]")
    ax.set_title(MODE_TITLES[mode])


def draw_algorithm_comparison(ax, cache_results, resolutions=None):
    """The solvers' own work only (warm cache run), one linear axis: "rest
    of solver" = everything except property lookups AND k evaluation --
    k is left out too because, like the lookups, it scales with how often
    a solver evaluates local states (Cell: once per cell)."""
    names = list(cache_results)
    rest = [cache_results[n]['warm cache']['split']['rest of solver'] * 1e3 for n in names]
    bars = ax.bar(names, rest, color=[SOLVER_COLORS[n] for n in names], edgecolor='black', linewidth=0.6)
    for bar, value in zip(bars, rest):
        ax.annotate(f"{value:.2f} ms" if value < 1 else f"{value:.1f} ms", (bar.get_x() + bar.get_width() / 2, value),
                    textcoords="offset points", xytext=(0, 4), ha='center', fontsize=11)
    if resolutions:
        ax.set_xticks(range(len(names)), [f"{n}\n({resolutions[n]})" if resolutions[n] else n for n in names])
        ax.tick_params(axis='x', labelsize=10)
    ax.set_ylim(0, max(rest) * 1.15)
    ax.set_ylabel("Wall time per run [ms]")
    ax.set_title("Solver work only: rest of solver\n(without lookups and k evaluation)")


# =============================================================================
# 2. Resolution
# =============================================================================

def run_resolution_benchmark(ops, geo, settings, modes):
    """results['NTU'/'Cell'] = list of (n, time [s], iterations, Q [kW]) over
    modes.benchmark_resolutions; results['LMTD'] = (time [s], iterations, Q [kW]).
    Cache on and emptied before every run (= a normal single run)."""
    t_benchmark = time.perf_counter()
    results = {'NTU': [], 'Cell': []}
    with _quiet():
        for n in modes.benchmark_resolutions:
            runs = _solver_runs(ops, geo, settings, n)
            for name in ('NTU', 'Cell'):
                times = []
                for _ in range(modes.benchmark_sweep_repeats):
                    _clear_property_cache()
                    t_start = time.perf_counter()
                    output = runs[name]()
                    times.append(time.perf_counter() - t_start)
                results[name].append((n, statistics.median(times), len(output[4]), output[1] / 1000))
        _clear_property_cache()
        t_start = time.perf_counter()
        output = _solver_runs(ops, geo, settings)['LMTD']()
        results['LMTD'] = (time.perf_counter() - t_start, len(output[4]), output[1] / 1000)
    _clear_property_cache()

    logger.info("[Benchmark: resolution] %s (took %.1f s)", machine_info(), time.perf_counter() - t_benchmark)
    for name in ('NTU', 'Cell'):
        logger.info("  %-4s %s", name,
                    "   ".join(f"n={n}: {t * 1e3:.0f} ms/{it} it" for n, t, it, _ in results[name]))
    return results


def plot_resolution(results, label: str):
    """2x2, all linear, NTU and Cell together in every panel: wall time,
    time per outer iteration, iterations to convergence, Q."""
    fig, ((ax_t, ax_ti), (ax_it, ax_q)) = plt.subplots(2, 2, figsize=(17, 11))
    fig.canvas.manager.set_window_title(f"{label} -- Solver performance vs. resolution")

    for name in ('Cell', 'NTU'):
        n, t, it, q = zip(*results[name])
        style = dict(marker='o', color=SOLVER_COLORS[name], label=name)
        ax_t.plot(n, [x * 1e3 for x in t], **style)
        ax_ti.plot(n, [x * 1e3 / i for x, i in zip(t, it)], **style)
        ax_it.plot(n, it, **style)
        ax_q.plot(n, q, **style)

    # LMTD has no resolution: a dashed reference line in every panel.
    lmtd_t, lmtd_it, lmtd_q = results['LMTD']
    for ax, value in ((ax_t, lmtd_t * 1e3), (ax_ti, lmtd_t * 1e3 / lmtd_it), (ax_it, lmtd_it), (ax_q, lmtd_q)):
        ax.axhline(value, color=COLOR_LMTD, linestyle='--', label="LMTD (no resolution)")

    for ax, title, ylabel in ((ax_t, "Wall time", "Wall time [ms]"),
                              (ax_ti, "Time per outer iteration", "Time per iteration [ms]"),
                              (ax_it, "Outer iterations to convergence", "Iterations"),
                              (ax_q, "Result Q", "Q [kW]")):
        ax.set_title(title, fontsize=13)
        ax.set_ylabel(ylabel)
        ax.set_xlabel("Cell segments / NTU elements per tube and row")
        ax.set_xlim(left=0)
        ax.legend(fontsize=10)
    for ax in (ax_t, ax_ti):
        ax.set_ylim(bottom=0)
    # Headroom above the (flat) top curve, so it doesn't sit on the frame.
    ax_it.set_ylim(0, max([lmtd_it] + [it for name in ('Cell', 'NTU') for _, _, it, _ in results[name]]) * 1.3)

    fig.suptitle(f"Solver Performance vs. Resolution — {label}\n"
                 f"(cache on, emptied before every run = normal single run; {machine_info()})",
                 fontsize=14, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    return fig
