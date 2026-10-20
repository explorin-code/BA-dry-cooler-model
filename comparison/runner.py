"""
runner.py -- runs a selection of analyses (analyses.py): compute -> CSV -> figures -> summary lines.

One code path for `python -m comparison`, `python main.py compare` and the GUI. Results are written to
output/comparison/<analysis>_<name>.csv so figures can be restyled without recomputing (plots_only).
Figures are saved to thesis_figures/out/ch4_results/cooler_performance/ (the module drawings to .../geometry/); with keep_open they also stay open for plt.show() / the GUI.
Progress (src/catalogue.run_items) and the summary lines go to the logger "src.comparison" (the GUI shows the "src" loggers).
"""

import logging
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from comparison import config
import src.parameters as parameters
from comparison.analyses import ANALYSES
from src.catalogue import run_items
from src.quiet import quiet_solvers

DATA_DIR = Path(__file__).resolve().parent.parent / 'output' / 'comparison'
logger = logging.getLogger('src.comparison')
logger.setLevel(logging.INFO)


def _write(key, data):
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    for name, df in data.items():
        df.to_csv(DATA_DIR / f'{key}_{name}.csv', index=False)


def _read(key):
    files = sorted(DATA_DIR.glob(f'{key}_*.csv'))
    if not files:
        raise FileNotFoundError(f"No stored results for '{key}' in {DATA_DIR} -- run it once without plots_only.")
    data = {}
    for f in files:
        df = pd.read_csv(f)
        if 'hour' in df:
            df['hour'] = pd.to_datetime(df['hour'])
        data[f.stem[len(key) + 1:]] = df
    return data


def run_analyses(keys=None, plots_only=False, keep_open=False) -> dict:
    """Runs the analyses `keys` (default: parameters.COMPARE_ANALYSES) -- results from the stored CSVs with plots_only.
    Returns {key: {name: DataFrame}}."""
    keys = list(parameters.COMPARE_ANALYSES if keys is None else keys)
    cfg = config.settings()                              # parameters.py as it is now, for this whole run
    if not keys:
        logger.info("No analysis selected (COMPARE_ANALYSES in parameters.py).")

    def get_data(analysis):
        if plots_only:
            return _read(analysis.key)
        data = analysis.compute(cfg)
        _write(analysis.key, data)
        return data

    def after(_):
        if not keep_open:                                # figures are saved by the plots; free the memory
            plt.close('all')

    with quiet_solvers():
        return run_items(ANALYSES, keys, args=(cfg,), logger=logger, get_data=get_data, after=after)
