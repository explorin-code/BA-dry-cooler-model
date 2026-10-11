"""
ui -- the parts of the Streamlit GUI (gui.py is the page itself):

  parameter_file.py  parameters.py -> the list of parameters with group, mode tag, default and help text
  labels.py          display labels (thesis symbols), modes, group order, the "exactly one nonzero" flow triples
  presets.py         what a parameter means "by default" on the chosen standard cooler
  widgets.py         one input widget per parameter, the sidebar form
  runner.py          running a mode with the form's values (log capture, one run at a time)
  results.py         the result tables and the figure tabs
"""
