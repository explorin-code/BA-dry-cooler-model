"""
chapter3_schematics.py
=======================
Thesis Chapter 3 figures: schematics of the discretisation, drawn from the
actual geometry and solvers.py's own path / neighbour logic, coloured by a
real model run (ambient operating point from parameters.py).

  circuit             side view of one tube circuit: rows x segments, the
                      coolant's serpentine path (solvers.cell_path_order)
                      through U-bends, each cell coloured by its local mean
                      coolant (tube) and air (background) temperature; its own
                      run at CIRCUIT_ROWS x CIRCUIT_SEGMENTS (6 x 10), text
                      only the four inlet/outlet temperature symbols.
  precooler_paths     dry operation vs. adiabatic precooling: ambient air
                      enters through the bypass or through the wetted pad
                      (damper picks one, PAD_IN_DRY_AIR_PATH = False); the
                      cooler is each scenario's own circuit, one colour scale.
  cells_bend(_real)   3D view of the tube ends at the U-bends: 2 cells of
                      row r above, 2 of row r-1 below (segment 0), bends at
                      the back joining the same tube index (get_coolant_inlet);
                      rows pulled apart / touching (real geometry, real d).
  cells_*             3D view of staggered Cell cells (3 tubes of row r-1,
                      2 of row r between them; versions in CELL_VERSIONS):
                      coolant directions alternating per row, air from two
                      lower cells merging into each upper cell and splitting
                      again towards the next row (get_staggered_air_inlet).

Qualitative on purpose: no numbers on the figures; the colours come from the
run, the colour bar only says cold / hot.

Standalone -- nothing in the model imports this file:
    python -m thesis_figures.chapter3_schematics [--format png] [--out folder]
                                                 [--rows N] [--segments N]
Writes every accepted version into thesis_figures/out/ (default).
"""

import argparse
from types import SimpleNamespace
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from matplotlib.transforms import Affine2D
from matplotlib.patches import FancyArrowPatch, Rectangle, Wedge
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

import src.parameters as parameters
from src.dry_cooler_physics import get_geometry
from src.operating_conditions import get_operating_conditions
from src.precooling import calc_precooler
from src.solver_settings import get_solver_settings
from src.solvers import get_coolant_inlet, get_staggered_air_inlet, solve_it_cell, solve_it_NTU

COLOR_OUTLINE = "#333333"
CMAP = "RdYlBu_r"                # one temperature scale for coolant and air
TUBE_WIDTH = 0.38                # tube band height as a fraction of the row


def _arrow(ax, start, end, color, lw=1.4, style="-|>", mutation_scale=10, **kw):
    ax.add_patch(FancyArrowPatch(start, end, arrowstyle=style, color=color, lw=lw,
                                 mutation_scale=mutation_scale, shrinkA=0, shrinkB=0, **kw))


# =============================================================================
# Side view: one tube circuit, rows x segments, coloured by a real run
# =============================================================================

def cell_inlet_temperatures(result, ops, geo, solver: str, tube: int):
    """Per-cell inlet temperatures of one tube circuit, [row, segment], for
    coolant and air (outlets are the result grids). Inlets come from the solvers' own
    neighbour rules: coolant along the serpentine (get_coolant_inlet); air
    straight from the row upstream (NTU) or via the staggered-bank average
    (Cell, get_staggered_air_inlet -- identical tubes in 2D)."""
    T_c, T_a = result.T_c_grid, result.T_a_grid
    N_r, _, n_segments = T_c.shape
    in_c = np.empty((N_r, n_segments))
    in_a = np.empty((N_r, n_segments))
    for r in range(N_r):
        for s in range(n_segments):
            c_in = get_coolant_inlet(r, tube, s, n_segments, T_c, ops, geo)
            if solver == "cell":
                a_in = get_staggered_air_inlet(r, tube, s, T_a, ops, geo)
            else:
                a_in = ops.T_a_i if r == 0 else T_a[r - 1, tube, s]
            in_c[r, s], in_a[r, s] = c_in, a_in
    return in_c, in_a


def cell_mean_temperatures(result, ops, geo, solver: str, tube: int):
    """Per-cell mean temperatures (inlet + outlet)/2, [row, segment]."""
    in_c, in_a = cell_inlet_temperatures(result, ops, geo, solver, tube)
    return ((in_c + result.T_c_grid[:, tube, :]) / 2.0,
            (in_a + result.T_a_grid[:, tube, :]) / 2.0)


def _temperature_colors(ops, t_range=None):
    """Shared colour scale from air inlet (cold) to coolant inlet (hot), or
    t_range = (cold, hot) when several scenarios share one scale."""
    vmin, vmax = t_range if t_range else (ops.T_a_i, ops.T_c_i)
    norm = Normalize(vmin=vmin, vmax=vmax)
    cmap = plt.get_cmap(CMAP)
    return norm, cmap, (lambda value: cmap(norm(value)))


def add_qualitative_colorbar(fig, mappable, ax, **kw):
    """Colour bar without numbers: the colours come from a real run, the
    figure stays qualitative."""
    cbar = fig.colorbar(mappable, ax=ax, **kw)
    cbar.set_ticks([mappable.norm.vmin, mappable.norm.vmax], labels=["cold", "hot"])
    cbar.ax.tick_params(length=0, labelsize=9)
    cbar.set_label("temperature", fontsize=9)
    return cbar


def draw_circuit(ax, result, ops, geo, solver: str = "cell", annotate: bool = True, t_range=None,
                 origin=None):
    """Cells at (x = segment, y = row). Cell background = local mean air
    temperature, tube = local mean coolant temperature, one shared colour
    scale. Air enters below row 0; coolant enters row N_r-1. Returns the
    ScalarMappable for a colourbar. annotate = False drops the temperature
    symbols; origin = (x, y) draws it shifted into a larger figure (no axis
    setup) -- used by the precooler figure."""
    n_patches, n_texts = len(ax.patches), len(ax.texts)
    N_r, N_t, n_segments = result.T_c_grid.shape
    tube = N_t // 2                                       # a middle tube (all identical with uniform inlet air)
    mean_c, mean_a = cell_mean_temperatures(result, ops, geo, solver, tube)
    T_c, T_a = result.T_c_grid[:, tube, :], result.T_a_grid[:, tube, :]

    norm, cmap, color = _temperature_colors(ops, t_range)
    half = TUBE_WIDTH / 2

    # Air: one rectangle per cell
    for r in range(N_r):
        for s in range(n_segments):
            ax.add_patch(Rectangle((s, r), 1, 1, facecolor=color(mean_a[r, s]), edgecolor="white", lw=0.5))

    # Tube: one band piece per cell, U-bends outside the grid
    for r in range(N_r):
        y = r + 0.5
        forward = r % 2 == 0
        for s in range(n_segments):
            ax.add_patch(Rectangle((s, y - half), 1, TUBE_WIDTH, facecolor=color(mean_c[r, s]),
                                   edgecolor="none", zorder=2))
        ax.add_patch(Rectangle((0, y - half), n_segments, TUBE_WIDTH, fill=False,
                               edgecolor=COLOR_OUTLINE, lw=0.8, zorder=3))
        for s in range(1, n_segments, max(1, n_segments // 4)):   # flow direction chevrons
            x = s + 0.5
            dx = 0.22 if forward else -0.22
            _arrow(ax, (x - dx, y), (x + dx, y), "white", lw=1.2, mutation_scale=9, zorder=4)
        if r > 0:                                         # bend from row r to row r-1 at this row's end
            x_end = n_segments if forward else 0
            s_end = n_segments - 1 if forward else 0
            theta1, theta2 = (-90, 90) if forward else (90, 270)
            ax.add_patch(Wedge((x_end, y - 0.5), 0.5 + half, theta1, theta2, width=TUBE_WIDTH,
                               facecolor=color(T_c[r, s_end]), edgecolor=COLOR_OUTLINE, lw=0.8, zorder=3))

    # Coolant in / out
    x_in = n_segments if (N_r - 1) % 2 == 1 else 0
    side = 1 if x_in else -1
    y_in = N_r - 0.5
    ax.add_patch(Rectangle((x_in, y_in - half), side * 1.0, TUBE_WIDTH, facecolor=color(ops.T_c_i),
                           edgecolor=COLOR_OUTLINE, lw=0.8, zorder=3))
    _arrow(ax, (x_in + side * 1.9, y_in), (x_in + side * 1.05, y_in), COLOR_OUTLINE, lw=1.6)
    ax.text(x_in + side * 2.0, y_in, "$T_{c,i}$",
            ha="left" if side > 0 else "right", va="center", fontsize=12)
    x_out = n_segments                                    # row 0 runs forward
    ax.add_patch(Rectangle((x_out, 0.5 - half), 1.0, TUBE_WIDTH, facecolor=color(result.T_c_o),
                           edgecolor=COLOR_OUTLINE, lw=0.8, zorder=3))
    _arrow(ax, (x_out + 1.05, 0.5), (x_out + 1.9, 0.5), COLOR_OUTLINE, lw=1.6)
    ax.text(x_out + 2.0, 0.5, "$T_{c,o}$",
            ha="left", va="center", fontsize=12)

    # Air arrows: inlet uniform, outlet coloured by each column's outlet temperature
    for s in range(n_segments):
        _arrow(ax, (s + 0.5, -0.85), (s + 0.5, -0.1), color(ops.T_a_i), lw=1.6)
        _arrow(ax, (s + 0.5, N_r + 0.1), (s + 0.5, N_r + 0.75), color(T_a[N_r - 1, s]), lw=1.6)
    ax.text(n_segments / 2, -0.95, "$T_{a,i}$",
            ha="center", va="top", fontsize=12)
    ax.text(n_segments / 2, N_r + 0.85, "$T_{a,o}$", ha="center", va="bottom", fontsize=12)

    if not annotate:
        for text in ax.texts[n_texts:]:
            text.remove()
    if origin is not None:                                # shift everything just drawn
        shift = Affine2D().translate(*origin) + ax.transData
        for artist in ax.patches[n_patches:] + ax.texts[n_texts:]:
            artist.set_transform(shift)
        return ScalarMappable(norm=norm, cmap=cmap)

    ax.set_xlim(-3.2, n_segments + 4.8)
    ax.set_ylim(-1.5, N_r + 1.4)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    return ScalarMappable(norm=norm, cmap=cmap)


# =============================================================================
# 3D view: staggered cells -- how a cell's inlet air is built
# =============================================================================
# Cell 3D's staggered-bank rule (get_staggered_air_inlet): three tubes of
# row r-1 below, two of row r between them above; each upper cell's inlet
# air is the mean of the two lower outlets, and its outlet feeds two tubes of
# the next row. Axes: x across the bank (s_t), y along the tubes (coolant),
# z = air flow. Schematic proportions.

def _box_edges(ax, x0, x1, y0, y1, z0, z1, **kw):
    """The 12 edges of an axis-aligned box."""
    corners = lambda x, y, z: (x, y, z)
    for a, b in [((x0, y0, z0), (x1, y0, z0)), ((x0, y1, z0), (x1, y1, z0)),
                 ((x0, y0, z1), (x1, y0, z1)), ((x0, y1, z1), (x1, y1, z1)),
                 ((x0, y0, z0), (x0, y1, z0)), ((x1, y0, z0), (x1, y1, z0)),
                 ((x0, y0, z1), (x0, y1, z1)), ((x1, y0, z1), (x1, y1, z1)),
                 ((x0, y0, z0), (x0, y0, z1)), ((x1, y0, z0), (x1, y0, z1)),
                 ((x0, y1, z0), (x0, y1, z1)), ((x1, y1, z0), (x1, y1, z1))]:
        ax.plot(*zip(corners(*a), corners(*b)), **kw)


def _tube(ax, x_c, z_c, y0, y1, radius, color, n=40):
    """Cylinder along y."""
    theta = np.linspace(0, 2 * np.pi, n)
    T, Y = np.meshgrid(theta, np.array([y0, y1]))
    ax.plot_surface(x_c + radius * np.cos(T), Y, z_c + radius * np.sin(T), color=color,
                    shade=True, linewidth=0, antialiased=True)


def draw_staggered_cells_3d(ax, result, ops, geo, solver: str = "cell", gap: float = 0.9,
                            depth: int = 1, fins_per_cell: int = None):
    """Three tubes of row r-1 (t..t+2) below, two of row r (t, t+1) between
    them above, r odd so that (r, t) sits between t and t+1 upstream (as in
    get_staggered_air_inlet); depth = cells per tube along the tube
    (consecutive segments). Each upper cell's inlet air merges the outlets
    of the two cells below it; each upper outlet splits towards the two tubes
    of the next row. gap = distance between the rows / s_l (0 = real
    geometry; > 0 pulls the rows apart to show the mixing). Fins: plates
    continuous across each row's tubes, inside the cells only, a few per
    cell (not the real pitch).
    Colours from the run: tubes by the cells' mean coolant temperature, air
    arrows by the air leaving / entering each cell."""
    N_r, N_t, n_segments = result.T_c_grid.shape
    r = N_r // 2 if (N_r // 2) % 2 == 1 else N_r // 2 + 1          # odd row with a row below it
    t = N_t // 2 - 1
    s0 = n_segments // 2 - depth // 2
    # y index j (0 = front): the same physical segment in both rows, so air
    # passes straight up. Upper row (odd) flows to the back with s falling,
    # lower row (even) to the front with s rising -> front cell = highest s.
    seg = [s0 + depth - 1 - j for j in range(depth)]
    norm, cmap, color = _temperature_colors(ops)
    if fins_per_cell is None:
        fins_per_cell = 4 if depth == 1 else 2

    def temps(row, tube, s):
        in_c, in_a = cell_inlet_temperatures(result, ops, geo, solver, tube)
        return ((in_c[row, s] + result.T_c_grid[row, tube, s]) / 2,   # mean coolant
                in_a[row, s], result.T_a_grid[row, tube, s])         # air in, air out

    s_t, s_l = geo.s_t * 1000, geo.s_l * 1000
    L, R = 1.2 * s_t, geo.d / 2 * 1000 * 1.4                      # schematic cell length / tube radius
    G = gap * s_l
    z_low, z_top = 0.0, s_l + G                                    # bottom of each row's cells
    y_len = depth * L
    lw_air = 3.0

    # Fins: plates normal to the tubes, continuous across the tubes of a row,
    # contained in the cells (not drawn in the gap between pulled-apart rows)
    for x0, x1, z0 in ((-0.5 * s_t, 2.5 * s_t, z_low), (0.0, 2.0 * s_t, z_top)):
        for f in range(depth * fins_per_cell):
            yf = (f + 0.5) * L / fins_per_cell
            ax.add_collection3d(Poly3DCollection([[(x0, yf, z0), (x1, yf, z0), (x1, yf, z0 + s_l), (x0, yf, z0 + s_l)]],
                                                 facecolor="#b0b0b0", edgecolor="#8a8a8a", lw=0.5, alpha=0.16))

    def tube_row(row, x_centres, z0, into_page, labels):
        zc = z0 + s_l / 2
        for k, x_c in enumerate(x_centres):
            means = []
            for j in range(depth):
                mean_c, _, _ = temps(row, t + k, seg[j])
                means.append(mean_c)
                _tube(ax, x_c, zc, j * L, (j + 1) * L, R, color(mean_c))
                _box_edges(ax, x_c - s_t / 2, x_c + s_t / 2, j * L, (j + 1) * L, z0, z0 + s_l,
                           color="#555555", lw=0.9, ls="--")
            # Coolant arrows in front of and behind the tube, in its flow direction
            c_first, c_last = (means[0], means[-1]) if into_page else (means[-1], means[0])
            if into_page:
                ax.quiver(x_c, -0.6 * L, zc, 0, 0.55 * L, 0, color=color(c_first), lw=3.5, arrow_length_ratio=0.35)
                ax.quiver(x_c, y_len + 0.05 * L, zc, 0, 0.55 * L, 0, color=color(c_last), lw=3.5,
                          arrow_length_ratio=0.35)
            else:
                ax.quiver(x_c, y_len + 0.6 * L, zc, 0, -0.55 * L, 0, color=color(c_first), lw=3.5,
                          arrow_length_ratio=0.35)
                ax.quiver(x_c, -0.05 * L, zc, 0, -0.55 * L, 0, color=color(c_last), lw=3.5, arrow_length_ratio=0.35)
            ax.text(x_c - s_t / 2 + 0.06 * s_t, 0, z0 + 0.06 * s_l, labels[k], fontsize=10, color="#333333",
                    zorder=20)

    tube_row(r - 1, [0.0, s_t, 2 * s_t], z_low, False, [f"$(r-1,\\ t{'+' + str(k) if k else ''})$" for k in range(3)])
    tube_row(r, [0.5 * s_t, 1.5 * s_t], z_top, True, ["$(r,\\ t)$", "$(r,\\ t+1)$"])

    # Air, per segment: straight in below, merge into each upper cell, split above
    z_from = z_low + s_l + 0.05 * max(G, 0.3 * s_l)
    z_mix = z_top - 0.05 * max(G, 0.3 * s_l)
    reach = max(G, 0.55 * s_l)                                     # length of the split arrows above
    for j, s in enumerate(seg):
        y_air = (j + 0.5) * L
        alpha = 1.0 if j == 0 else 0.3                              # back segments faded: front one tells the story
        outs = []
        for k in range(3):
            _, a_in, a_out = temps(r - 1, t + k, s)
            ax.quiver(k * s_t, y_air, z_low - 0.65 * s_l, 0, 0, 0.6 * s_l, color=color(a_in), lw=lw_air,
                      arrow_length_ratio=0.3, alpha=alpha)
            outs.append(a_out)
        if G > 0:                                                  # merging arrows need the gap
            for k, (x_from, a_out_low) in enumerate(zip((0.0, s_t, 2 * s_t), outs)):
                for x_to in ((k - 0.5) * s_t, (k + 0.5) * s_t):
                    ax.quiver(x_from, y_air, z_from, (x_to - x_from) * 0.92, 0, z_mix - z_from,
                              color=color(a_out_low), lw=lw_air, arrow_length_ratio=0.18, alpha=alpha)
        for k in range(2):
            _, a_in, a_out = temps(r, t + k, s)
            x_c = (k + 0.5) * s_t
            ax.scatter([x_c], [y_air], [z_top], s=40, color=color(a_in), edgecolor="#333333", zorder=10,
                       alpha=alpha)
            for dx in (-0.5 * s_t, 0.5 * s_t):
                ax.quiver(x_c, y_air, z_top + s_l + 0.05 * reach, dx * 0.92, 0, 0.88 * reach, color=color(a_out),
                          lw=lw_air, arrow_length_ratio=0.18, alpha=alpha)

    if depth > 1:                                                  # segment indices along the tube
        for j, s in enumerate(seg):
            ax.text(2.5 * s_t + 0.1 * s_t, (j + 0.5) * L, z_low, f"$s{'+' + str(s - s0) if s > s0 else ''}$",
                    fontsize=10, color="#333333")

    z_max = z_top + s_l + 1.1 * reach
    ax.set_xlim(-1.0 * s_t, 3.0 * s_t)
    ax.set_ylim(-0.7 * L, y_len + 0.7 * L)
    ax.set_zlim(z_low - 1.0 * s_l, z_max)
    ax.set_box_aspect((4.0 * s_t, y_len + 1.4 * L, z_max - z_low + 1.0 * s_l), zoom=1.1)
    ax.view_init(elev=14 if depth == 1 else 18, azim=-62)        # a bit steeper with depth: segments separate
    ax.set_axis_off()
    return ScalarMappable(norm=norm, cmap=cmap)


def plot_staggered_cells_3d(result, ops, geo, solver: str = "cell", gap: float = 0.9, depth: int = 1):
    fig = plt.figure(figsize=(10 + 1.5 * (depth - 1), 8 if gap > 0.5 else 7))
    if fig.canvas.manager:
        fig.canvas.manager.set_window_title("Staggered cells (3D)")
    ax = fig.add_axes((-0.12, -0.08, 1.05, 1.12), projection="3d")   # 3D axes leave wide margins
    mappable = draw_staggered_cells_3d(ax, result, ops, geo, solver, gap=gap, depth=depth)
    add_qualitative_colorbar(fig, mappable, None, cax=fig.add_axes((0.9, 0.25, 0.02, 0.5)))
    return fig


def _bend(ax, a, b, y0, radius, color, n_arc=40, n_circ=30):
    """U-bend tube from point a = (x, z) to b = (x, z), both at y = y0,
    bulging towards +y: a half circle in the plane of a-b and the y axis."""
    a = np.array([a[0], y0, a[1]])
    b = np.array([b[0], y0, b[1]])
    centre, rho = (a + b) / 2, np.linalg.norm(a - b) / 2
    u = (a - centre) / rho                                         # in-plane, towards a
    yhat = np.array([0.0, 1.0, 0.0])
    binormal = np.cross(u, yhat)                                   # normal of the bend plane
    phi = np.linspace(0, np.pi, n_arc)
    theta = np.linspace(0, 2 * np.pi, n_circ)
    P, T = np.meshgrid(phi, theta, indexing="ij")
    radial = np.cos(P)[..., None] * u + np.sin(P)[..., None] * yhat   # unit vector centre -> arc point
    pts = (centre + rho * radial + radius * (np.cos(T)[..., None] * radial
                                             + np.sin(T)[..., None] * binormal))
    for i in range(0, n_arc - 1, 2):                               # short pieces: mplot3d sorts per surface
        piece = pts[i:i + 3]
        ax.plot_surface(piece[..., 0], piece[..., 1], piece[..., 2], color=color, shade=True, linewidth=0,
                        antialiased=True)


def draw_bend_cells_3d(ax, result, ops, geo, solver: str = "cell", gap: float = 0.9, fins_per_cell: int = 4):
    """Tube ends at the U-bends: two cells of row r (tubes t, t+1) above, two
    of row r-1 below, all at the bend end of the pass (segment 0). The
    coolant runs to the back in row r, through the U-bend into the same
    tube index of row r-1 (get_coolant_inlet) -- diagonal in the staggered
    bank -- and back to the front. Air as in draw_staggered_cells_3d: two
    lower cells merge into each upper cell (the right upper cell's second
    source, tube t+2, lies outside the picture), outlets split again.
    gap = 0: real geometry, rows touching -- the merge is drawn across the
    interface inside the cells.
    Colours from the run."""
    N_r, N_t, n_segments = result.T_c_grid.shape
    r = N_r // 2 if (N_r // 2) % 2 == 1 else N_r // 2 + 1          # odd row with a row below it
    t, s = N_t // 2 - 1, 0                                         # segment 0 = bend end of both rows
    norm, cmap, color = _temperature_colors(ops)

    def temps(row, tube):
        in_c, in_a = cell_inlet_temperatures(result, ops, geo, solver, tube)
        return ((in_c[row, s] + result.T_c_grid[row, tube, s]) / 2,   # mean coolant
                in_a[row, s], result.T_a_grid[row, tube, s])         # air in, air out

    s_t, s_l = geo.s_t * 1000, geo.s_l * 1000
    L = 1.2 * s_t                                                  # schematic cell length
    R = geo.d / 2 * 1000 * (1.4 if gap > 0 else 1.0)               # real tube radius when the rows touch
    G = gap * s_l
    z_low, z_top = 0.0, s_l + G
    lw_air = 3.0
    lower_x, upper_x = [0.0, s_t], [0.5 * s_t, 1.5 * s_t]
    stub = 0.25 * L                                                # straight piece between fins and bend

    # Fins: continuous across each row's tubes, inside the cells only
    for x0, x1, z0 in ((-0.5 * s_t, 1.5 * s_t, z_low), (0.0, 2.0 * s_t, z_top)):
        for f in range(fins_per_cell):
            yf = (f + 0.5) * L / fins_per_cell
            ax.add_collection3d(Poly3DCollection([[(x0, yf, z0), (x1, yf, z0), (x1, yf, z0 + s_l), (x0, yf, z0 + s_l)]],
                                                 facecolor="#b0b0b0", edgecolor="#8a8a8a", lw=0.5, alpha=0.16))

    for k in range(2):
        mean_up, _, _ = temps(r, t + k)
        mean_low, _, _ = temps(r - 1, t + k)
        x_up, x_low = upper_x[k], lower_x[k]
        z_up_c, z_low_c = z_top + s_l / 2, z_low + s_l / 2
        # Upper cell: coolant from the front into the back
        _tube(ax, x_up, z_up_c, 0, L, R, color(mean_up))
        _box_edges(ax, x_up - s_t / 2, x_up + s_t / 2, 0, L, z_top, z_top + s_l, color="#555555", lw=0.9, ls="--")
        ax.quiver(x_up, -0.6 * L, z_up_c, 0, 0.55 * L, 0, color=color(mean_up), lw=3.5, arrow_length_ratio=0.35)
        # Lower cell: coolant from the bend back to the front
        _tube(ax, x_low, z_low_c, 0, L, R, color(mean_low))
        _box_edges(ax, x_low - s_t / 2, x_low + s_t / 2, 0, L, z_low, z_low + s_l, color="#555555", lw=0.9, ls="--")
        ax.quiver(x_low, -0.05 * L, z_low_c, 0, -0.55 * L, 0, color=color(mean_low), lw=3.5, arrow_length_ratio=0.35)
        # Straight stubs out of the fin block and the U-bend, at the coolant
        # temperature leaving the upper row (= entering the lower row)
        T_bend = result.T_c_grid[r, t + k, s]
        _tube(ax, x_up, z_up_c, L, L + stub, R, color(T_bend))
        _tube(ax, x_low, z_low_c, L, L + stub, R, color(T_bend))
        _bend(ax, (x_up, z_up_c), (x_low, z_low_c), L + stub, R, color(T_bend))
        ax.text(x_up - s_t / 2 + 0.06 * s_t, 0, z_top + 0.06 * s_l, f"$(r,\\ t{'+1' if k else ''})$",
                fontsize=10, color="#333333", zorder=20)
        ax.text(x_low - s_t / 2 + 0.06 * s_t, 0, z_low + 0.06 * s_l, f"$(r-1,\\ t{'+1' if k else ''})$",
                fontsize=10, color="#333333", zorder=20)

    # Air in the middle of the cell length
    y_air = 0.5 * L
    reach = max(G, 0.55 * s_l)                                     # length of the split arrows above
    if G > 0:                                                      # merge in the gap between the rows
        z_from, z_mix = z_low + s_l + 0.05 * G, z_top - 0.05 * G
    else:                                                          # rows touching: merge across the interface
        z_from, z_mix = z_low + 0.8 * s_l, z_top + 0.2 * s_l
    lower_out = {}
    for k, x_c in enumerate(lower_x):
        _, a_in, a_out = temps(r - 1, t + k)
        ax.quiver(x_c, y_air, z_low - 0.65 * s_l, 0, 0, 0.6 * s_l, color=color(a_in), lw=lw_air, arrow_length_ratio=0.3)
        lower_out[k] = a_out
        for x_to in (x_c - 0.5 * s_t, x_c + 0.5 * s_t):
            ax.quiver(x_c, y_air, z_from, (x_to - x_c) * 0.92, 0, z_mix - z_from, color=color(a_out),
                      lw=lw_air, arrow_length_ratio=0.18)
    # Second source of the right upper cell: tube t+2 of row r-1, outside the picture
    _, _, a_out_t2 = temps(r - 1, t + 2)
    ax.quiver(2 * s_t, y_air, z_from, -0.5 * s_t * 0.92, 0, z_mix - z_from, color=color(a_out_t2),
              lw=lw_air, arrow_length_ratio=0.18, alpha=0.45)
    for k, x_c in enumerate(upper_x):
        _, a_in, a_out = temps(r, t + k)
        ax.scatter([x_c], [y_air], [z_mix if G == 0 else z_top], s=40, color=color(a_in), edgecolor="#333333", zorder=10)
        for dx in (-0.5 * s_t, 0.5 * s_t):
            ax.quiver(x_c, y_air, z_top + s_l + 0.05 * reach, dx * 0.92, 0, 0.88 * reach, color=color(a_out),
                      lw=lw_air, arrow_length_ratio=0.18)

    rho = np.hypot(0.5 * s_t, s_l + G) / 2
    y_max = L + stub + rho + R
    z_max = z_top + s_l + 1.1 * reach
    ax.set_xlim(-1.0 * s_t, 2.5 * s_t)
    ax.set_ylim(-0.7 * L, y_max)
    ax.set_zlim(z_low - 1.0 * s_l, z_max)
    ax.set_box_aspect((3.5 * s_t, y_max + 0.7 * L, z_max - z_low + 1.0 * s_l), zoom=1.1)
    ax.view_init(elev=16, azim=-35)                               # from the side: bends visibly at the back
    ax.set_axis_off()
    return ScalarMappable(norm=norm, cmap=cmap)


def plot_bend_cells_3d(result, ops, geo, solver: str = "cell", gap: float = 0.9):
    fig = plt.figure(figsize=(10, 8))
    if fig.canvas.manager:
        fig.canvas.manager.set_window_title("Cells at the U-bends (3D)")
    ax = fig.add_axes((-0.02, -0.08, 1.05, 1.12), projection="3d")   # 3D axes leave wide margins
    mappable = draw_bend_cells_3d(ax, result, ops, geo, solver, gap=gap)
    add_qualitative_colorbar(fig, mappable, None, cax=fig.add_axes((0.9, 0.25, 0.02, 0.5)))
    return fig


# =============================================================================
# Precooler: two air paths into the dry cooler
# =============================================================================
# Two panels, dry operation and adiabatic precooling, on one colour scale:
# ambient air enters either directly (bypass) or through the wetted pad
# (cooled towards the cooling limit) -- a damper picks one path, as in the
# model (PAD_IN_DRY_AIR_PATH = False: separate bypass inlet). The cooler on
# top is the circuit figure of each scenario's own Cell run.

COLOR_WATER = "#2b6cb0"
COLOR_WALL = "#444444"


def _channel(ax, x0, x1, y0, y1, active):
    """Duct walls; inactive ducts are greyed out."""
    if not active:
        ax.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0, facecolor="#ececec", edgecolor="none", zorder=0))
    for x in (x0, x1):
        ax.plot([x, x], [y0, y1], color=COLOR_WALL, lw=1.6, zorder=3)


def _damper(ax, x0, x1, y, is_open):
    """Damper at the duct entry: blades parallel to the flow when open,
    across the duct when closed."""
    xm, half = (x0 + x1) / 2, (x1 - x0) / 2
    if is_open:
        for x in (xm - half / 2, xm + half / 2):
            ax.plot([x, x], [y - 0.35, y + 0.35], color=COLOR_WALL, lw=2.2, zorder=4)
            ax.plot([x], [y], marker="o", ms=3.5, color=COLOR_WALL, zorder=5)
    else:
        ax.plot([x0, x1], [y, y], color=COLOR_WALL, lw=3.0, zorder=4)
        for x in (xm - half / 2, xm + half / 2):
            ax.plot([x], [y], marker="o", ms=3.5, color=COLOR_WALL, zorder=5)


def _pad(ax, x0, x1, y0, y1, wet):
    """Evaporative pad with a water distributor on top (droplets if wet)."""
    ax.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0, facecolor="#cfe3ef" if wet else "#e4e4e4",
                           edgecolor=COLOR_WALL, lw=1.0, hatch="////", zorder=2))
    y_pipe = y1 + 0.35
    ax.plot([x0 - 1.2, x1 - 0.2], [y_pipe, y_pipe], color=COLOR_WATER if wet else "#aaaaaa", lw=3.0,
            solid_capstyle="butt", zorder=3)
    if wet:
        for x in np.linspace(x0 + 0.5, x1 - 0.5, 4):
            ax.plot([x], [y_pipe - 0.22], marker="v", ms=5, color=COLOR_WATER, zorder=4)


def draw_precooler_panel(ax, ops_ambient, result, ops, geo, precooled: bool, t_range):
    """Ducts, damper, pad and the cooler (circuit of this scenario's run)."""
    norm, cmap, color = _temperature_colors(ops, t_range)
    pad_x, bypass_x = (0.5, 4.5), (5.5, 9.5)
    y_entry, y_top = 2.6, 9.2
    pad_y = (4.6, 6.8)

    # Ambient air below both ducts (faded in front of the closed damper)
    for (x0, x1), is_open in ((pad_x, precooled), (bypass_x, not precooled)):
        xm = (x0 + x1) / 2
        _arrow(ax, (xm, 0.4), (xm, y_entry - 0.55), color(ops_ambient.T_a_i), lw=3, mutation_scale=16,
               alpha=1.0 if is_open else 0.3)
    ax.text(5.0, 0.0, "$T_{a,i}$", ha="center", va="top", fontsize=12)

    # Pad duct and bypass duct; the damper opens one of them
    _channel(ax, *pad_x, y_entry, y_top, active=precooled)
    _channel(ax, *bypass_x, y_entry, y_top, active=not precooled)
    _damper(ax, *pad_x, y_entry, is_open=precooled)
    _damper(ax, *bypass_x, y_entry, is_open=not precooled)
    _pad(ax, pad_x[0] + 0.25, pad_x[1] - 0.25, *pad_y, wet=precooled)
    if precooled:
        xm = sum(pad_x) / 2
        _arrow(ax, (xm, y_entry + 0.5), (xm, pad_y[0] - 0.15), color(ops_ambient.T_a_i), lw=3, mutation_scale=16)
        _arrow(ax, (xm, pad_y[1] + 0.75), (xm, y_top - 0.1), color(ops.T_a_i), lw=3, mutation_scale=16)
        ax.text(pad_x[0] - 0.25, (pad_y[1] + y_top) / 2 + 0.4, "$T_{a,pc}$", ha="right", va="center", fontsize=12)
    else:
        xm = sum(bypass_x) / 2
        _arrow(ax, (xm, y_entry + 0.5), (xm, y_top - 0.1), color(ops.T_a_i), lw=3, mutation_scale=16)

    # Plenum below the cooler, filled with the air the cooler receives
    ax.add_patch(Rectangle((0.5, y_top), 9.0, 1.0, facecolor=color(ops.T_a_i), edgecolor="none", zorder=1))
    for x in (0.5, 9.5):
        ax.plot([x, x], [y_top, y_top + 1.0], color=COLOR_WALL, lw=1.6, zorder=3)

    # The dry cooler: circuit of this scenario, unlabelled, on the shared scale
    draw_circuit(ax, result, ops, geo, "cell", annotate=False, t_range=t_range, origin=(0.0, y_top + 2.0))
    N_r = result.T_c_grid.shape[0]

    ax.set_xlim(-2.0, 12.5)
    ax.set_ylim(-0.9, y_top + 2.0 + N_r + 1.2)
    ax.set_aspect("equal")
    ax.axis("off")
    return ScalarMappable(norm=norm, cmap=cmap)


def run_precooler_scenarios(n_segments: int = None, n_rows: int = None):
    """Ambient and precooled Cell runs at the circuit size. The precooled
    state comes from precooling.calc_precooler (VDI M8); if the decider
    would not engage at this operating point it is forced, for the figure."""
    result, ops, geo = run_model("cell", n_segments or CIRCUIT_SEGMENTS, n_rows or CIRCUIT_ROWS)
    ops_pc, engaged = calc_precooler(ops, SimpleNamespace(cell=result))
    if not engaged:
        ops_pc, _ = calc_precooler(ops, SimpleNamespace(cell=SimpleNamespace(T_c_o=float("inf"))))
    result_pc = solve_it_cell(n_segments=result.T_c_grid.shape[2], ops=ops_pc, geo=geo,
                              settings=get_solver_settings())
    return (result, ops), (result_pc, ops_pc), geo


def plot_precooler_paths(ambient, precooled, geo):
    (result, ops), (result_pc, ops_pc) = ambient, precooled
    t_range = (ops_pc.T_a_i, ops.T_c_i)                   # coldest air .. hottest coolant
    fig, axes = plt.subplots(1, 2, figsize=(12, 8))
    if fig.canvas.manager:
        fig.canvas.manager.set_window_title("Precooler air paths")
    draw_precooler_panel(axes[0], ops, result, ops, geo, precooled=False, t_range=t_range)
    mappable = draw_precooler_panel(axes[1], ops, result_pc, ops_pc, geo, precooled=True, t_range=t_range)
    axes[0].set_title("dry operation", fontsize=12)
    axes[1].set_title("adiabatic precooling", fontsize=12)
    fig.subplots_adjust(left=0.02, right=0.88, bottom=0.01, top=0.95, wspace=0.05)
    add_qualitative_colorbar(fig, mappable, None, cax=fig.add_axes((0.9, 0.3, 0.018, 0.4)))
    return fig


# =============================================================================
# Figures
# =============================================================================

# The circuit rendering: smallest size that still shows the counterflow
# effect (coolant leaving close to the air inlet temperature). Its own run,
# independent of N_ROWS / CELL_N_SEGMENTS; --rows / --segments override it.
CIRCUIT_ROWS = 6
CIRCUIT_SEGMENTS = 10

# Accepted versions of the staggered-cell figure: name -> (gap / s_l, cells along the tube)
CELL_VERSIONS = {
    "cells_exploded": (0.9, 1),        # rows pulled apart, merge/split arrows clearly visible
    "cells_exploded_2seg": (0.9, 2),   # exploded + a second segment along the tube (2x2 over 3x2)
}

def run_model(solver: str = "cell", n_segments: int = None, n_rows: int = None):
    """Ambient operating point from parameters.py, solved by Cell or NTU at
    n_segments (default: CELL_N_SEGMENTS / NTU_N_ELEMENTS); n_rows overrides
    N_ROWS for this run only. Returns (result, ops, geo)."""
    original_rows = parameters.N_ROWS
    if n_rows is not None:
        parameters.N_ROWS = n_rows
    try:
        geo = get_geometry()
        ops = get_operating_conditions(geo)
    finally:
        parameters.N_ROWS = original_rows
    settings = get_solver_settings()
    if solver == "cell":
        result = solve_it_cell(n_segments=n_segments, ops=ops, geo=geo, settings=settings)
    else:
        result = solve_it_NTU(n_elements=n_segments, ops=ops, geo=geo, settings=settings)
    return result, ops, geo


def plot_flow_schematic(result, ops, geo, solver: str = "cell"):
    """Temperature-coloured tube circuit of a solved result."""
    N_r, _, n_segments = result.T_c_grid.shape
    unit = 0.5                                                     # inches per cell, both directions
    fig, ax = plt.subplots(figsize=(unit * (n_segments + 8.0) + 1.6, unit * (N_r + 2.9) + 1.3))
    if fig.canvas.manager:
        fig.canvas.manager.set_window_title("Flow path schematic")
    mappable = draw_circuit(ax, result, ops, geo, solver)
    add_qualitative_colorbar(fig, mappable, ax, fraction=0.025, pad=0.01, shrink=0.75)
    fig.tight_layout()
    return fig


def render_all(out_dir: Path, fmt: str = "pdf", solver: str = "cell", n_segments: int = None,
               n_rows: int = None):
    """Renders every accepted figure into out_dir; returns the written paths."""
    out_dir.mkdir(parents=True, exist_ok=True)
    run = run_model(solver, n_segments, n_rows)
    cell_run = run if solver == "cell" else run_model("cell", n_segments, n_rows)   # neighbour mixing: Cell only
    circuit_run = run_model(solver, n_segments or CIRCUIT_SEGMENTS, n_rows or CIRCUIT_ROWS)
    figures = {"circuit": plot_flow_schematic(*circuit_run, solver=solver)}
    figures["precooler_paths"] = plot_precooler_paths(*run_precooler_scenarios(n_segments, n_rows))
    for name, (gap, depth) in CELL_VERSIONS.items():
        figures[name] = plot_staggered_cells_3d(*cell_run, solver="cell", gap=gap, depth=depth)
    figures["cells_bend"] = plot_bend_cells_3d(*cell_run, solver="cell", gap=0.9)
    figures["cells_bend_real"] = plot_bend_cells_3d(*cell_run, solver="cell", gap=0.0)
    paths = []
    for name, fig in figures.items():
        path = out_dir / f"{name}.{fmt}"
        fig.savefig(path, dpi=200)
        plt.close(fig)
        paths.append(path)
    return paths


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Chapter 3 schematics (qualitative), coloured by a model run.")
    parser.add_argument("--solver", choices=("cell", "ntu"), default="cell",
                        help="run for the circuit figure (the staggered-cell figures always use Cell 3D)")
    parser.add_argument("--segments", type=int, default=None,
                        help="segments / elements per pass (default: CELL_N_SEGMENTS / NTU_N_ELEMENTS)")
    parser.add_argument("--rows", type=int, default=None, help="tube rows for this rendering (default: N_ROWS)")
    parser.add_argument("--out", default=str(Path(__file__).parent / "out"), help="output folder")
    parser.add_argument("--format", default="pdf", help="file format, e.g. pdf, png, svg")
    args = parser.parse_args()
    for written in render_all(Path(args.out), args.format, args.solver, args.segments, args.rows):
        print(written)
