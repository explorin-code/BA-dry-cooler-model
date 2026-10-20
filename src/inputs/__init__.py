"""
inputs -- what a run is given.

  param_loader.py         builds the settings dataclasses from the constants of parameters.py (units: mm and % in, SI inside)
  solver_settings.py      tolerances, relaxation factors, resolutions of the solvers
  run_modes.py            what the solver comparison runs and produces (solvers, outputs, benchmark options)
  operating_conditions.py the air and water state of an operating point
  geometry.py             the cooler's geometry (Geometry) and its materials; get_geometry() builds it from a cell plus overrides
  cooler_sets.py          the standard elementary cells P1-P3 and their sources
  weather.py              the Munich climate year, wet-bulb temperature and the representative states
  build_climate_data.py   rebuilds climate_data/ from the raw DWD download
"""
