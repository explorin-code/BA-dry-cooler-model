"""
plot_style.py
==============
Shared plotting style/helpers for the convergence figure (run_scenario.py)
and the spatial-profile figure (plot_profiles.py): one color per solver,
a color-lightening helper, and generic (domain-independent) round-number
axis-tick math.
"""

import numpy as np
import matplotlib.colors as mc

COLOR_LMTD = "#1b7837"  # green
COLOR_NTU = "#762a83"   # purple
COLOR_CELL = "#e08214"  # orange


def lighten_color(color, factor=0.5):
    """Blend `color` toward white. factor=0 -> unchanged, factor=1 -> white."""
    r, g, b = mc.to_rgb(color)
    return (r + (1 - r) * factor, g + (1 - g) * factor, b + (1 - b) * factor)


_NICE_FRACTIONS = [1, 2, 2.5, 5, 10]


def _nice_step(raw_step):
    """Rounds raw_step up to the nearest 'nice' number (1/2/2.5/5/10 x
    10^n) -- so tick labels are whole/round numbers, not arbitrary fractions."""
    if raw_step <= 0:
        return 1.0
    exponent = np.floor(np.log10(raw_step))
    fraction = raw_step / 10**exponent
    for f in _NICE_FRACTIONS:
        if fraction <= f:
            return f * 10**exponent
    return 10 * 10**exponent


def _next_nice_step(step):
    """The next larger step in the same 1/2/2.5/5/10 x 10^n sequence."""
    exponent = np.floor(np.log10(step))
    fraction = round(step / 10**exponent, 6)
    idx = _NICE_FRACTIONS.index(fraction) if fraction in _NICE_FRACTIONS else 0
    if idx + 1 < len(_NICE_FRACTIONS):
        return _NICE_FRACTIONS[idx + 1] * 10**exponent
    return _NICE_FRACTIONS[0] * 10**(exponent + 1)


def _nice_ticks(data_lo, data_hi, n_intervals: int = 4):
    """n_intervals+1 evenly spaced round-number ticks covering [data_lo,
    data_hi]. Using the same n_intervals for every axis in a panel makes
    unrelated quantities' (e.g. temperature, k, dQ/dL) gridlines align,
    even though their actual values differ."""
    span = data_hi - data_lo
    if span <= 0:
        span = abs(data_hi) if data_hi != 0 else 1.0
    step = _nice_step(span / n_intervals)
    for _ in range(4):   # a couple of nice-step bumps is always enough in practice
        nice_lo = np.floor(data_lo / step) * step
        ticks = nice_lo + step * np.arange(n_intervals + 1)
        if ticks[-1] >= data_hi - 1e-9:
            return ticks
        step = _next_nice_step(step)
    return ticks


# =============================================================================
# Column-aligned text boxes (cooler-results figure)
# =============================================================================
# Mathtext subscripts ($T_{c,o}$) are drawn in a different font than the
# surrounding monospace text, so aligning columns by padding strings with
# spaces breaks. draw_text_grid instead places every cell as its own text at
# measured column positions. Positions are in points relative to an anchor
# in axes coordinates, so the layout survives figure resizing.

def _text_width_pt(ax, s: str, fontsize: float) -> float:
    """Rendered width of s in points (exact, incl. spaces and mathtext),
    measured with a temporary text on the figure's renderer."""
    if not s:
        return 0.0
    fig = ax.figure
    probe = ax.text(0, 0, s, fontsize=fontsize, family='monospace')
    width_px = probe.get_window_extent(renderer=fig.canvas.get_renderer()).width
    probe.remove()
    return width_px * 72 / fig.dpi


def draw_text_grid(ax, x_anchor, y_anchor, blocks, facecolor, edgecolor, fontsize=11,
                   label_gap_pt=4.0, pair_gap_pt=16.0, line_spacing=1.6, pad_pt=10.0):
    """Draws a rounded box of column-aligned text centered at (x_anchor,
    y_anchor) in axes coordinates. blocks: list of blocks stacked
    vertically; each block is a list of rows; a row is either a list of
    cells or a plain string spanning the whole box (centered). Cell rows
    are [row label, name, value, name, value, ...]: a small gap follows
    each name ("T ="), a large one each value and the row label.
    fontsize is the maximum: the box shrinks its font to fit the axes'
    width, and refits whenever the window is resized."""
    from matplotlib.patches import FancyBboxPatch
    from matplotlib.transforms import Affine2D, ScaledTranslation

    fig = ax.figure
    transform = (Affine2D().scale(1 / 72) + fig.dpi_scale_trans
                 + ScaledTranslation(x_anchor, y_anchor, ax.transAxes))
    artists = []

    def layout(fs):
        """Column positions per block and the total box width, at font size fs."""
        scale = fs / fontsize
        layouts, total_width = [], 0.0
        for block in blocks:
            n_cols = max((len(row) for row in block if not isinstance(row, str)), default=0)
            widths = [0.0] * n_cols
            for row in block:
                if not isinstance(row, str):
                    for j, cell in enumerate(row):
                        widths[j] = max(widths[j], _text_width_pt(ax, cell, fs))
            gaps = [(label_gap_pt if j % 2 == 1 else pair_gap_pt) * scale for j in range(n_cols)]
            x_cols = [sum(widths[:j]) + sum(gaps[:j]) for j in range(n_cols)]
            block_width = (x_cols[-1] + widths[-1]) if n_cols else 0.0
            span_width = max((_text_width_pt(ax, row, fs) for row in block if isinstance(row, str)),
                             default=0.0)
            total_width = max(total_width, block_width, span_width)
            layouts.append(x_cols)
        return layouts, total_width

    def draw(*_):
        for artist in artists:
            artist.remove()
        artists.clear()

        # Fit the font to the axes' current width (never larger than fontsize).
        _, natural_width = layout(fontsize)
        available = ax.get_window_extent().width * 72 / fig.dpi - 2 * pad_pt
        fs = fontsize * min(1.0, available / natural_width) if natural_width > 0 else fontsize
        layouts, total_width = layout(fs)
        pad = pad_pt * fs / fontsize

        line_height = fs * line_spacing
        height = sum(len(block) for block in blocks) * line_height
        x_left, y_top = -total_width / 2, height / 2
        row_index = 0
        for block, x_cols in zip(blocks, layouts):
            for row in block:
                y = y_top - (row_index + 0.72) * line_height          # text baseline
                if isinstance(row, str):
                    artists.append(ax.text(0, y, row, transform=transform, ha='center', va='baseline',
                                           fontsize=fs, family='monospace', clip_on=False))
                else:
                    for x, cell in zip(x_cols, row):
                        artists.append(ax.text(x_left + x, y, cell, transform=transform, ha='left',
                                               va='baseline', fontsize=fs, family='monospace', clip_on=False))
                row_index += 1
        box = FancyBboxPatch((x_left - pad, -height / 2 - pad / 2), total_width + 2 * pad, height + pad,
                             boxstyle="round,pad=0,rounding_size=6", transform=transform,
                             facecolor=facecolor, edgecolor=edgecolor, alpha=0.95, clip_on=False, zorder=0)
        ax.add_patch(box)
        artists.append(box)

    draw()
    fig.canvas.mpl_connect('resize_event', draw)
