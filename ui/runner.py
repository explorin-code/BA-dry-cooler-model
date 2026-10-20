"""
runner.py -- running a mode with the form's values: the values go onto the parameters module (everything else back to the
file's defaults), the model's log is captured, only one run at a time per server.
"""

import logging
import threading

import matplotlib.pyplot as plt
import streamlit as st

import src.parameters as parameters
from src.inputs.run_modes import get_run_modes
from src.workflow.scenario_pipeline import run_scenarios


class _Capture(logging.Handler):
    """Collects the model's log output (summary / insight lines, warnings)."""
    def __init__(self):
        super().__init__()
        self.lines = []

    def emit(self, record):
        self.lines.append(self.format(record))


@st.cache_resource
def _model_lock():
    """One lock for the whole server: all browser tabs share the model's modules (parameter overrides, the
    property-cache switch and the instrumentation recorder), so only one run may execute at a time."""
    return threading.Lock()


def run_model(mode, values, defaults):
    """Runs the model under the global lock (the model state is shared), waiting if another browser tab is still running."""
    lock = _model_lock()
    if not lock.acquire(blocking=False):
        st.warning("Another run (other browser tab) is still in progress -- waiting for it ...")
        lock.acquire()
    try:
        return _run(mode, values, defaults)
    finally:
        lock.release()


def _run(mode, values, defaults):
    """Sets the GUI values on the parameters module (everything else back to the file's defaults), runs the
    mode, returns (results, log text, [(title, figure)])."""
    for name, value in defaults.items():
        setattr(parameters, name, value)
    for name, value in values.items():
        setattr(parameters, name, value)

    capture = _Capture()
    capture.setFormatter(logging.Formatter("%(message)s"))
    src_logger, warn_logger = logging.getLogger("src"), logging.getLogger("py.warnings")
    insight = mode == "solver" and get_run_modes().insight_mode
    src_logger.setLevel(logging.DEBUG if insight else logging.INFO)
    logging.captureWarnings(True)
    for logger in (src_logger, warn_logger):
        logger.addHandler(capture)

    plt.close("all")
    try:
        if mode == "solver":
            results = run_scenarios(get_run_modes())
        else:
            from comparison.runner import run_analyses
            results = run_analyses(keep_open=True)
    finally:
        for logger in (src_logger, warn_logger):
            logger.removeHandler(capture)

    figures = []
    for number in plt.get_fignums():
        fig = plt.figure(number)
        figures.append((fig.canvas.manager.get_window_title() if fig.canvas.manager else f"Figure {number}", fig))
    return results, "\n".join(capture.lines), figures
