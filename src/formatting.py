"""
formatting.py
==============
Text shared by the benchmark logs and figures: time labels and the machine description.
"""

import platform

import CoolProp


def fmt_ms(value_ms: float) -> str:
    """A time given in ms: two decimals below 10 ms (LMTD, NTU table, ...), one above."""
    return f"{value_ms:.2f} ms" if value_ms < 10 else f"{value_ms:.1f} ms"


def fmt_duration(seconds: float) -> str:
    """A duration in the unit that fits it: ms, s, min or h."""
    if seconds < 1:
        return f"{seconds * 1e3:.0f} ms"
    if seconds < 120:
        return f"{seconds:.1f} s"
    if seconds < 7200:
        return f"{seconds / 60:.1f} min"
    return f"{seconds / 3600:.1f} h"


def machine_info() -> str:
    return (f"{platform.processor() or platform.machine()}, Python {platform.python_version()}, "
            f"CoolProp {CoolProp.__version__}")
