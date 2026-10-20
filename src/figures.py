"""
figures.py
===========
Saving figures, shared by both workflows: a PNG in thesis_figures/out/<folder>/ and the window title (the GUI shows
every open figure under its title; main.py's plt.show() names the windows with it).

The folders under thesis_figures/out/ follow the thesis chapters (names defined here, never spelled out elsewhere):

  CH3       ch3_methodology/                   Chapter 3, methodology and models: the schematics (thesis_figures/chapter3_schematics.py)
  SOLVER    ch4_results/solver/                Chapter 4, solver comparison: profiles, caching, resolution, accuracy, annual (+ raw results, tables)
  COOLER    ch4_results/cooler_performance/    Chapter 4, cooler comparison: pinch, fan, sensitivities, cost, ambient states
  GEOMETRY  ch4_results/geometry/              Chapter 4, the three modules: the drawings and the size / mass comparison
save_data / load_data keep the raw results of the long benchmarks (annual, resolution) next to their figures.
"""

import pickle
import re
from pathlib import Path

OUTPUT_DIR = Path(__file__).resolve().parent.parent / 'thesis_figures' / 'out'
CH3 = 'ch3_methodology'
SOLVER = 'ch4_results/solver'
COOLER = 'ch4_results/cooler_performance'
GEOMETRY = 'ch4_results/geometry'


def slug(title: str) -> str:
    """'Caching performance (log)' -> 'caching_performance_log'."""
    return re.sub(r'[^a-z0-9]+', '_', title.lower()).strip('_')


def save_text(folder: str, name: str, text: str):
    """Writes OUTPUT_DIR/folder/name (a text file, e.g. a LaTeX table to paste into the thesis) and returns its path."""
    out = OUTPUT_DIR / folder
    out.mkdir(parents=True, exist_ok=True)
    (out / name).write_text(text)
    return out / name


def save_data(folder: str, name: str, data) -> str:
    """Pickles the raw results of a long benchmark next to its figures (OUTPUT_DIR/folder/name), so the figures can be redrawn
    without the run (see load_data); returns the path."""
    out = OUTPUT_DIR / folder
    out.mkdir(parents=True, exist_ok=True)
    with open(out / name, 'wb') as f:
        pickle.dump(data, f)
    return str(out / name)


def load_data(folder: str, name: str):
    """The results saved by save_data."""
    with open(OUTPUT_DIR / folder / name, 'rb') as f:
        return pickle.load(f)


def save_figure(fig, folder: str, name: str, title: str = None, dpi: int = 200):
    """Writes OUTPUT_DIR/folder/name.png and titles the figure's window (default: the name)."""
    out = OUTPUT_DIR / folder
    out.mkdir(parents=True, exist_ok=True)
    fig.savefig(out / f'{name}.png', dpi=dpi, bbox_inches='tight')
    if fig.canvas.manager is not None:
        fig.canvas.manager.set_window_title(title or name)
