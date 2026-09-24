"""
main.py
=======
Entry point. Runs the dry-cooler scenario(s) via scenario_pipeline.run_scenarios()
and shows the resulting figures together. See scenario_pipeline.py for what
actually gets solved/plotted.
"""

import matplotlib.pyplot as plt

from src.scenario_pipeline import run_scenarios


def run():
    run_scenarios()
    plt.show()


if __name__ == "__main__":
    run()
