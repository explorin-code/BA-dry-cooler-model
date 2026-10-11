"""
figures.py
===========
Saving figures, shared by both workflows: a PNG in thesis_figures/out/<folder>/ and the window title (the GUI shows
every open figure under its title; main.py's plt.show() names the windows with it).

  solver comparison   folder 'solver'  (named after the figure title, e.g. cooler_results_dry.png)
  cooler comparison   folder 'ch4'
"""

import re
from pathlib import Path

OUTPUT_DIR = Path(__file__).resolve().parent.parent / 'thesis_figures' / 'out'


def slug(title: str) -> str:
    """'Caching performance (log)' -> 'caching_performance_log'."""
    return re.sub(r'[^a-z0-9]+', '_', title.lower()).strip('_')


def save_figure(fig, folder: str, name: str, title: str = None, dpi: int = 200):
    """Writes OUTPUT_DIR/folder/name.png and titles the figure's window (default: the name)."""
    out = OUTPUT_DIR / folder
    out.mkdir(parents=True, exist_ok=True)
    fig.savefig(out / f'{name}.png', dpi=dpi, bbox_inches='tight')
    if fig.canvas.manager is not None:
        fig.canvas.manager.set_window_title(title or name)
