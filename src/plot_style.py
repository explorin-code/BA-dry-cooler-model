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
