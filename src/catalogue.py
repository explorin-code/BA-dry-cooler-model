"""
catalogue.py
=============
The one pattern both workflows use to offer "things to run": a catalogue of `Item`s, each with a compute, a plot and a
log step, run in catalogue order by `run_items`.

  solver comparison   src/outputs.py        cooler results, convergence, caching / resolution / annual benchmarks
  cooler comparison   comparison/analyses.py   geometry, fan control, sensitivity, ...

The steps take the same extra arguments `*args` (the solver comparison passes its run context, the cooler
comparison nothing). Which items run is a list of keys (parameters.SOLVER_OUTPUTS / COMPARE_ANALYSES, the command
line, the GUI's multiselects).
"""

import logging
import time
from dataclasses import dataclass
from typing import Callable


@dataclass(frozen=True)
class Item:
    key: str
    title: str
    compute: Callable                      # (*args) -> data (any object; the cooler comparison: {name: DataFrame})
    plot: Callable                         # (data, *args) -> None, builds figures
    log: Callable = lambda data, *args: []     # (data, *args) -> a few summary lines


def run_items(catalogue: dict, keys, args=(), logger=None, get_data=None, after=None) -> dict:
    """Runs the items `keys` of `catalogue` in catalogue order: compute -> plot -> summary lines, with the progress
    logged. get_data(item) replaces `item.compute(*args)` (the cooler comparison reads / stores CSVs there); after(item)
    is called when an item is done (the cooler comparison closes its figures there). Returns {key: data}."""
    logger = logger or logging.getLogger(__name__)
    unknown = set(keys) - set(catalogue)
    if unknown:
        raise ValueError(f"Unknown item(s) {sorted(unknown)} -- one of {list(catalogue)}")
    results = {}
    for key in (k for k in catalogue if k in keys):
        item = catalogue[key]
        t_start = time.perf_counter()
        logger.info("[%s]", item.title)
        data = get_data(item) if get_data else item.compute(*args)
        item.plot(data, *args)
        for line in item.log(data, *args):
            logger.info("  %s", line)
        logger.info("  (done in %.1f s)", time.perf_counter() - t_start)
        results[key] = data
        if after:
            after(item)
    return results
