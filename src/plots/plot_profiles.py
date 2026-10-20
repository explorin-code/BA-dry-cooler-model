"""
plot_profiles.py
=================
The profile figures. Every panel stands on its own (title, axis labels, legend), so single panels can be cut out of a page.

  plot_cooler_results      one scenario: inputs and economics on top, then one row per solver: its result box and its
                           profile (temperatures, k and dQ/dL on three y-axes along the coolant's flow path)
  plot_profile_grid        two water temperatures (the run's own solid, the second dashed): one row per solver, columns
                           temperatures | k | dQ/dL, absolute or normalised (theta, k/k_mean, (dQ/dL)/(Q/L))
  plot_profile_collapsed   the normalised grid with all solvers in one graph per quantity (appendix: the lines coincide)
  plot_result_table        the solvers' result boxes of both water temperatures

Axes: k spans at least +-5 % around its mean and dQ/dL starts at zero, so the figures show how LITTLE k varies and how
large dQ/dL is. Cell's and NTU's profiles are genuinely computed; LMTD's is a reconstruction from its scalars (the
log-mean profile, profiles.py, docs/model-notes.md). The path axis has no labels: only the row changes are marked.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D

from src.solvers import SOLVER_FAMILY, ScenarioResult
from src.solvers.profiles import calc_lmtd_profile, calc_ntu_profile, calc_cell_profile
from src.plots.plot_style import (COLOR_CELL, COLOR_LMTD, COLOR_NTU, TUM_BLUE, TUM_BLUE_540, TUM_GREY_50, TUM_GREY_80,
                            family_color, lighten_color, nice_tick_sets, draw_text_grid, in_large_format, in_solver_print_style, print_context, PRINT_WIDTH, PRINT_DPI, PRINT_LW, PRINT_DASHED)
from src.plots.figures import save_figure, slug, SOLVER
from src.workflow.run_scenario import input_condition_cells, economics_cells, output_condition_cells

COLOR_K = TUM_GREY_80      # dark gray -- neutral, not tied to hot/cold identity
COLOR_DQDL = TUM_BLUE_540  # dark blue -- distinct from the LMTD/NTU/Cell colours
K_MIN_WINDOW = 0.05        # the k axis spans at least +-5 % around the mean, so a 3 % drift does not look like a collapse
LW = PRINT_LW
CASE_LS = ('-', PRINT_DASHED)   # the run's own water temperature solid, the second one dashed

BOX_FACECOLORS = {family: lighten_color(color, 0.85) for family, color in
                  (("LMTD", COLOR_LMTD), ("NTU", COLOR_NTU), ("Cell", COLOR_CELL))}      # by algorithm family


def _profile(key, r, ops, geo, n_segments, n_elements):
    """The spatial profile of one solver: LMTD's is reconstructed from its scalars, NTU's and Cell's reduce the
    grids they computed."""
    family = SOLVER_FAMILY[key]
    if family == 'LMTD':
        return calc_lmtd_profile(r, ops, geo)
    if family == 'NTU':
        return calc_ntu_profile(r, geo, n_elements)
    return calc_cell_profile(r, geo, n_segments)


def _path_axis(ax, n_rows):
    """x axis along the coolant path, without labels: faint lines where the coolant changes rows."""
    ax.set_xlim(0, 1)
    ax.set_xticks([])
    for i in range(1, n_rows):
        ax.axvline(i / n_rows, color=TUM_GREY_50, lw=0.8, ls=':', alpha=0.8, zorder=0)


def _range(arrays, include=()):
    values = np.concatenate([np.asarray(a).ravel() for a in arrays] + [np.asarray(include, dtype=float)])
    return float(values.min()), float(values.max())


def common_ticks(scenarios, geo, n_segments: int, n_elements: int):
    """Axis ticks (T, k, dQ/dL) that span the profiles of ALL given scenarios [(result, ops), ...], so figures (dry and
    precooled, or the two water temperatures) have identical coordinate systems."""
    profiles = [_profile(k, result.results[k], ops, geo, n_segments, n_elements)
                for result, ops in scenarios for k in result.results]
    k_all = np.concatenate([p.k for p in profiles])
    k_mean = float(k_all.mean())
    return tuple(nice_tick_sets([
        _range([a for p in profiles for a in (p.T_c, p.T_a)]),
        _range([k_all], include=[(1 - K_MIN_WINDOW) * k_mean, (1 + K_MIN_WINDOW) * k_mean]),
        _range([p.dQdL for p in profiles], include=[0.0])]))


def _plot_solver_panel(ax_T, profile, key, ticks, n_rows, title):
    """One solver's panel: temperature on the primary axis, k and dQ/dL on two offset twin axes. Each axis's ylim spans
    exactly its own nice ticks, so the gridlines of all three coincide. Self-contained: title, x label, one legend."""
    T_ticks, k_ticks, dQdL_ticks = ticks
    color = family_color(key)
    color_air = lighten_color(color, 0.55)

    ax_T.plot(profile.x_frac, profile.T_c, color=color, linewidth=2.2, label=r"$T_c$ (coolant)")
    ax_T.plot(profile.x_frac, profile.T_a, color=color_air, linewidth=2.2, linestyle='--', label=r"$T_a$ (air)")
    ax_T.set_ylabel("Temperature [°C]", color=color)
    ax_T.set_ylim(T_ticks[0], T_ticks[-1])
    ax_T.set_yticks(T_ticks)
    ax_T.tick_params(axis='y', labelcolor=color)
    _path_axis(ax_T, n_rows)

    ax_k = ax_T.twinx()
    ax_k.grid(False)   # only ax_T should draw gridlines -- twin axes' own grid was rendering on top of the legend
    ax_k.plot(profile.x_frac, profile.k, color=COLOR_K, linewidth=1.6, linestyle=':', label=r"$k$")
    ax_k.set_ylabel(r"$k$ [W/m²K]", color=COLOR_K, labelpad=4)
    ax_k.set_ylim(k_ticks[0], k_ticks[-1])
    ax_k.set_yticks(k_ticks)
    ax_k.tick_params(axis='y', labelcolor=COLOR_K)

    ax_dQdL = ax_T.twinx()
    ax_dQdL.grid(False)
    ax_dQdL.spines['right'].set_position(('outward', 78))
    ax_dQdL.plot(profile.x_frac, profile.dQdL, color=COLOR_DQDL, linewidth=1.6, linestyle='-.', label=r"$\mathrm{d}\dot{Q}/\mathrm{d}l$")
    ax_dQdL.set_ylabel(r"$\mathrm{d}\dot{Q}/\mathrm{d}l$ [W/m]", color=COLOR_DQDL)
    ax_dQdL.set_ylim(dQdL_ticks[0], dQdL_ticks[-1])
    ax_dQdL.set_yticks(dQdL_ticks)
    ax_dQdL.tick_params(axis='y', labelcolor=COLOR_DQDL)

    ax_T.set_title(title, color=color, fontweight="bold")

    # Legend goes on ax_dQdL (created last, so its render pass draws on top of both other axes' lines).
    lines = ax_T.get_lines()[:2] + ax_k.get_lines() + ax_dQdL.get_lines()
    handles = lines + [Line2D([], [], color=TUM_GREY_50, lw=0.8, ls=':')]
    labels = [l.get_label() for l in lines] + ["row change"]
    ax_dQdL.legend(handles, labels, fontsize=8, loc='upper right', framealpha=1.0, facecolor='white', edgecolor='gray')


def _solver_title(key):
    return SOLVER_FAMILY[key]


@in_large_format
def plot_cooler_results(result: ScenarioResult, ops, geo, label: str, n_segments: int, n_elements: int,
                        P_pump=None, P_fan=None, m_dot_ev=None, ticks=None, T_amb=None):
    """Returns the figure (does NOT call plt.show()): one row per solver that ran (result box | profile). The
    boxes of the other solvers show their deviation from the reference (most detailed) solver.
    P_pump/P_fan/m_dot_ev: economics, None for whichever aren't computed. T_amb: the ambient air temperature the pinch of the
    whole cooler is taken against (default: the air entering the coil, i.e. the dry case). ticks: from common_ticks, for identical
    axes in several figures; default = the ranges of this figure."""
    keys = list(result.results)
    profiles = {k: _profile(k, result.results[k], ops, geo, n_segments, n_elements) for k in keys}
    if ticks is None:
        ticks = common_ticks([(result, ops)], geo, n_segments, n_elements)

    fig = plt.figure(figsize=(21, 3.2 + 3.6 * len(keys)))
    grid = GridSpec(1 + len(keys), 2, figure=fig, height_ratios=[0.6] + [1] * len(keys), width_ratios=[1.2, 1.45])

    ax_header = fig.add_subplot(grid[0, :])
    ax_header.axis("off")
    draw_text_grid(ax_header, 0.5, 0.66, input_condition_cells(ops, geo), "whitesmoke", "gray")
    draw_text_grid(ax_header, 0.5, 0.08, economics_cells(P_pump, P_fan, m_dot_ev),
                   lighten_color(TUM_BLUE, 0.9), TUM_BLUE)

    for row, key in enumerate(keys, start=1):
        family = SOLVER_FAMILY[key]
        ax_text = fig.add_subplot(grid[row, 0])
        ax_text.axis("off")
        is_reference = key == result.reference_key
        draw_text_grid(ax_text, 0.5, 0.5,
                       output_condition_cells(family, result.results[key], ops.T_a_i if T_amb is None else T_amb,
                                              reference=None if is_reference else result.reference,
                                              reference_label=SOLVER_FAMILY[result.reference_key]),
                       BOX_FACECOLORS[family], family_color(key))
        _plot_solver_panel(fig.add_subplot(grid[row, 1]), profiles[key], key, ticks, geo.N_r, _solver_title(key))

    fig.suptitle(f"Cooler results - {label}", fontsize=15, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 0.95, 0.95])
    save_figure(fig, SOLVER, slug(f"Cooler results - {label}"), f"Cooler results - {label}")
    return fig


# --- two water temperatures: grids, collapsed grid, result table --------------------------------------------

ROW_HEIGHT = 1.2                       # one graph row [in]: with the labels about 2:1 (width : height) per graph
FOOTNOTE = (r"$\theta = (T - T_{a,i})/(T_{c,i} - T_{a,i})$;   $\bar{k}$: mean of $k$ along the coolant path;   "
            r"$\dot{Q}/l_\mathrm{tot}$: mean heat flow per path length")

GRID_YLABELS = {False: (r"$T$ [°C]", r"$k$ [W/m²K]", r"$\mathrm{d}\dot{Q}/\mathrm{d}l$ [kW/m]"),
                True: (r"$\theta$ [–]", r"$k\,/\,\bar{k}$ [–]", r"$(\mathrm{d}\dot{Q}/\mathrm{d}l)\,/\,(\dot{Q}/l_\mathrm{tot})$ [–]")}
GRID_TITLES = {False: ("temperatures", "coefficient k", r"heat flow $\mathrm{d}\dot{Q}/\mathrm{d}l$"),
               True: ("normalised T", "normalised k", r"normalised $\mathrm{d}\dot{Q}/\mathrm{d}l$")}


def _case_curves(cases, geo, n_segments, n_elements, normalised):
    """{(case index, solver key): (x, T_c, T_a, k, dQ/dL)}, absolute [°C, °C, W/m2K, kW/m] or normalised
    (theta = (T - T_a,i)/(T_c,i - T_a,i); k / its mean along the path; (dQ/dL) / (Q/L))."""
    curves = {}
    for n, (ops, result) in enumerate(cases):
        for key, r in result.results.items():
            p = _profile(key, r, ops, geo, n_segments, n_elements)
            if normalised:
                theta = lambda T: (T - ops.T_a_i) / (ops.T_c_i - ops.T_a_i)
                curves[n, key] = (p.x_frac, theta(p.T_c), theta(p.T_a), p.k / p.k.mean(), p.dQdL / (r.Q_dot / geo.l_tot))
            else:
                curves[n, key] = (p.x_frac, p.T_c, p.T_a, p.k, p.dQdL / 1e3)
    return curves


def _centred_ticks(arrays):
    """Five ticks symmetric around 1 (1 in the middle of the axis) that span the data: the half range is the smallest
    round value that covers the largest deviation from 1 plus a margin."""
    deviation = 1.05 * max(float(np.abs(np.asarray(a) - 1.0).max()) for a in arrays)
    half = next((h for h in (0.01, 0.02, 0.025, 0.05, 0.1, 0.2, 0.25, 0.5, 1.0) if h >= deviation), 1.0)
    return 1.0 + half * np.array([-1.0, -0.5, 0.0, 0.5, 1.0])


def _column_ticks(curves, normalised):
    """Ticks of the three columns (temperature, k, dQ/dL), shared by all graphs of a column: from the data of both cases.
    Absolute: k at least +-5 % around its mean, dQ/dL from zero. Normalised: theta from 0 to 1, k/k_mean and
    (dQ/dL)/(Q/L) symmetric around 1."""
    values = list(curves.values())
    if normalised:
        return [nice_tick_sets([(0.0, 1.0)])[0], _centred_ticks([v[3] for v in values]), _centred_ticks([v[4] for v in values])]
    k_all = np.concatenate([v[3] for v in values])
    k_mean = float(k_all.mean())
    ranges = [_range([a for v in values for a in (v[1], v[2])]),
              _range([k_all], include=[(1 - K_MIN_WINDOW) * k_mean, (1 + K_MIN_WINDOW) * k_mean]),
              _range([v[4] for v in values], include=[0.0])]
    return [nice_tick_sets([r])[0] for r in ranges]


def _style_graph(ax, ticks, ylabel, title, n_rows, color=None):
    ax.set_ylim(ticks[0], ticks[-1])
    ax.set_yticks(ticks)
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontweight='bold' if color else 'normal', color=color or 'black')
    _path_axis(ax, n_rows)


def _case_handles(cases):
    return [Line2D([], [], color='k', lw=LW, ls=CASE_LS[n], label=f"{ops.T_c_i:.0f} °C water")
            for n, (ops, _) in enumerate(cases)]


def _check_cases(cases):
    if len(cases) != 2:
        raise ValueError("The profile grids compare exactly two water temperatures.")


def _footer(fig, height, margin_in):
    """The definition of the normalised quantities as one line at the bottom of the figure (not in the axes);
    returns the bottom margin to leave free for it (and what else sits there), `margin_in` inches, as figure fraction."""
    fig.text(0.5, 0.05 / height, FOOTNOTE, ha='center', va='bottom', fontsize=7, color=TUM_GREY_80)
    return margin_in / height


def _legend(ax, handles):
    ax.legend(handles=handles, loc='best', ncol=2 if len(handles) > 3 else 1)


def plot_profile_grid(cases, geo, n_segments: int, n_elements: int, normalised: bool):
    """cases = [(ops, ScenarioResult), ...] (two water temperatures). One row per solver, columns temperatures | k | dQ/dL;
    both temperatures in every graph (the second dashed), every graph with its own legend and title; the graphs of a
    column share their axes. Print size (print_context)."""
    _check_cases(cases)
    curves = _case_curves(cases, geo, n_segments, n_elements, normalised)
    ticks = _column_ticks(curves, normalised)
    keys = list(cases[0][1].results)
    with print_context(talk=True):
        height = ROW_HEIGHT * len(keys) + (0.3 if normalised else 0.0)
        fig, axes = plt.subplots(len(keys), 3, figsize=(PRINT_WIDTH, height), squeeze=False)
        for row, key in enumerate(keys):
            color = family_color(key)
            for n in range(len(cases)):
                x, T_c, T_a, k, dQdL = curves[n, key]
                axes[row, 0].plot(x, T_c, color=color, lw=LW, ls=CASE_LS[n])
                axes[row, 0].plot(x, T_a, color=lighten_color(color, 0.55), lw=LW, ls=CASE_LS[n])
                axes[row, 1].plot(x, k, color=color, lw=LW, ls=CASE_LS[n])
                axes[row, 2].plot(x, dQdL, color=color, lw=LW, ls=CASE_LS[n])
            medium = [Line2D([], [], color=color, lw=LW, label=r"$T_c$"),
                      Line2D([], [], color=lighten_color(color, 0.55), lw=LW, label=r"$T_a$")]
            for col in range(3):
                _style_graph(axes[row, col], ticks[col], GRID_YLABELS[normalised][col],
                             f"{_solver_title(key)}: {GRID_TITLES[normalised][col]}", geo.N_r, color)
                _legend(axes[row, col], (medium if col == 0 else []) + _case_handles(cases))
        bottom = _footer(fig, height, 0.3) if normalised else 0.0
        fig.tight_layout(rect=[0, bottom, 1, 1], h_pad=0.8, w_pad=0.8)
        kind = "normalised" if normalised else "absolute"
        save_figure(fig, SOLVER, f'profile_grid_{kind}', f"Profile grid ({kind})", dpi=PRINT_DPI)
    return fig


def plot_profile_collapsed(cases, geo, n_segments: int, n_elements: int):
    """The normalised quantities with all solvers in one graph per quantity. Slightly transparent lines: where solvers
    coincide the colours blend (that is the point -- the normalised profiles hardly depend on the water temperature).
    The legends and the definitions of the normalised quantities sit inside the graphs (in their empty corners). Print size (print_context)."""
    _check_cases(cases)
    curves = _case_curves(cases, geo, n_segments, n_elements, True)
    ticks = _column_ticks(curves, True)
    keys = list(cases[0][1].results)
    with print_context(talk=True):
        fig, axes = plt.subplots(1, 3, figsize=(PRINT_WIDTH, 1.9))
        for (n, key), (x, T_c, T_a, k, dQdL) in curves.items():
            color = family_color(key)
            kw = dict(lw=LW, ls=CASE_LS[n], alpha=0.8)
            axes[0].plot(x, T_c, color=color, **kw)
            axes[0].plot(x, T_a, color=lighten_color(color, 0.5), **kw)
            axes[1].plot(x, k, color=color, **kw)
            axes[2].plot(x, dQdL, color=color, **kw)
        for col, ax in enumerate(axes):
            title = GRID_TITLES[True][col]
            _style_graph(ax, ticks[col], GRID_YLABELS[True][col], title[0].upper() + title[1:], geo.N_r)
        # legends and definitions inside the graphs (the empty corners): the two temperatures | the solvers | the two water temperatures
        legend_kw = dict(loc='upper right', fontsize=6.5, borderpad=0.3, labelspacing=0.2, handlelength=1.5, borderaxespad=0.3)
        axes[0].legend(handles=[Line2D([], [], color='k', lw=LW, label=r"$T_c$"),
                                Line2D([], [], color=lighten_color('k', 0.5), lw=LW, label=r"$T_a$")], **legend_kw)
        axes[1].legend(handles=[Line2D([], [], color=family_color(key), lw=LW, label=_solver_title(key)) for key in keys], **legend_kw)
        axes[2].legend(handles=_case_handles(cases), **legend_kw)
        for ax, text in zip(axes, (r"$\theta = \dfrac{T - T_{a,i}}{T_{c,i} - T_{a,i}}$", r"$\bar{k}$: mean of $k$" + "\nalong the coolant path",
                                   r"$\dot{Q}/l_\mathrm{tot}$: mean heat flow" + "\nper path length")):
            ax.text(0.04, 0.05, text, transform=ax.transAxes, ha='left', va='bottom', fontsize=6.5, color=TUM_GREY_80,
                    bbox=dict(facecolor='white', edgecolor='none', alpha=0.85, pad=1.2))
        fig.tight_layout(w_pad=0.8)
        save_figure(fig, SOLVER, 'profile_collapsed_normalised', "Profile graphs (normalised, solvers combined)", dpi=PRINT_DPI)
    return fig


TABLE_FONT, TABLE_LINE_SPACING, TABLE_PAD_PT = 7.0, 1.45, 10.0     # result-table text [pt, line height / font size, box padding]


@in_solver_print_style
def plot_result_table(cases, geo):
    """The result boxes of every solver at both water temperatures for the appendix: page width (PRINT_WIDTH), the
    boxes stacked without gaps -- per water temperature the inputs, then per solver outlets, pinch and Q (deviation
    from the reference solver in K / %) and the dimensionless numbers. Each box is as high as its text."""
    _check_cases(cases)
    keys = list(cases[0][1].results)
    boxes = []                               # (blocks, face colour, edge colour, extra gap above [in])
    for n, (ops, result) in enumerate(cases):
        boxes.append((input_condition_cells(ops, geo), "whitesmoke", "gray", 0.0 if n == 0 else 0.25))
        for key in keys:
            family = SOLVER_FAMILY[key]
            reference = None if key == result.reference_key else result.reference
            blocks = output_condition_cells(family, result.results[key], ops.T_a_i, reference, SOLVER_FAMILY[result.reference_key])
            boxes.append((blocks, BOX_FACECOLORS[family], family_color(key), 0.0))
    heights = [(sum(len(block) for block in blocks) * TABLE_FONT * TABLE_LINE_SPACING + TABLE_PAD_PT * 1.6) / 72
               for blocks, *_ in boxes]
    gaps = [gap + 0.04 for *_, gap in boxes]
    total = sum(heights) + sum(gaps) + 0.05
    fig = plt.figure(figsize=(PRINT_WIDTH, total))
    top = total - 0.05
    for (blocks, face, edge, _), h, gap in zip(boxes, heights, gaps):
        top -= gap
        ax = fig.add_axes([0.0, (top - h) / total, 1.0, h / total])
        ax.axis("off")
        draw_text_grid(ax, 0.5, 0.5, blocks, face, edge, fontsize=TABLE_FONT, line_spacing=TABLE_LINE_SPACING, pad_pt=TABLE_PAD_PT)
        top -= h
    save_figure(fig, SOLVER, 'profile_results_table', "Profile results table", dpi=PRINT_DPI)
    return fig
