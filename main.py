"""
main.py
=======
Runs the dry-cooler scenario twice: once on ambient conditions, once with
precooling.calc_precooler() applied. Both windows are shown together so
precooling's effect (or lack of one, while the decider is still a stub)
is visible side by side.
"""

import matplotlib.pyplot as plt

from src.operating_conditions import get_operating_conditions
from src.precooling import calc_precooler
from src.run_scenario import run_scenario

ops_ambient = get_operating_conditions()
run_scenario(ops_ambient, label="Ambient (no precooling)")

ops_final, was_precooled = calc_precooler(ops_ambient)
if was_precooled:
    run_scenario(ops_final, label="Precooled")

plt.show()