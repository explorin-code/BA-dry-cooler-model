"""
comparison -- the cooler comparison: the standard coolers (src/cooler_sets.py) raced against each other.

Every cooler cell is scaled to the same water flow (sizing.py: a fixed tube count per cell, N_t_ref in cooler_sets.py,
with common rows and one common tube length) and then analysed in representative ambient states and over the climate year (src/weather.py).

  config.py        `settings()` = parameters.py read once per run into a Settings object (`cfg`, the first argument of
                   every function that needs it), and the literature ranges
  evaluate.py      one operating point -> one result row (solver, pressure drops, fan/pump power, water)
  sizing.py        the sized module of a cell (fixed tube count; design_tubes re-derives it)
  operation.py     ambient states, fixed-speed year, face-velocity sweep, variable-speed fan control, summaries
  fan_curve.py     the controlled fan as a function of the cooler inlet air temperature
  sensitivity.py, rows.py, tubes.py, wall.py   design analyses around the sized modules
  analyses.py      the catalogue: compute / plot / log of every analysis      <- start reading here
  runner.py        runs a selection of them (CSV cache, figures, summary lines)
  plots_*.py, plot_geometry.py   the figures (thesis symbols)

Entry points: `python -m comparison`, `python main.py compare`, the GUI's comparison mode.
"""
