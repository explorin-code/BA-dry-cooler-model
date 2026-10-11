"""
profiling.py
=============
Optional instrumentation of the model's hot spots, for the caching benchmark: functions marked `@timed('name')` count
their calls and add up their wall time into the active `Recorder`; without an active recorder they just call through
(one `is None` test). Nothing is patched at run time.

  'lookup'   fluid_properties.get_fluid_properties / get_air_properties (incl. the CoolProp calls of a cache miss)
  'coolprop' every CoolProp call of the model (PropsSI / HAPropsSI)
  'k'        heat_transfer_core.calc_overall_k

    with profiling.record() as rec:
        run()
    rec.calls['k'], rec.time['coolprop'], ...
"""

import functools
import time
from collections import defaultdict
from contextlib import contextmanager

_active = None


class Recorder:
    def __init__(self):
        self.calls = defaultdict(int)
        self.time = defaultdict(float)


def timed(name):
    """Decorator: count the calls and the wall time of the function under `name` while a recorder is active."""
    def decorator(fn):
        @functools.wraps(fn)
        def wrapper(*args, **kwargs):
            recorder = _active
            if recorder is None:
                return fn(*args, **kwargs)
            t_start = time.perf_counter()
            try:
                return fn(*args, **kwargs)
            finally:
                recorder.time[name] += time.perf_counter() - t_start
                recorder.calls[name] += 1
        return wrapper
    return decorator


@contextmanager
def record():
    """Activates a fresh Recorder for the block and yields it."""
    global _active
    previous, _active = _active, Recorder()
    try:
        yield _active
    finally:
        _active = previous
