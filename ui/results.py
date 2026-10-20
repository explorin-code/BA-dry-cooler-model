"""
results.py -- what the GUI shows after a run: the result tables (solver mode), the summary log (comparison mode), and one tab per figure.
"""

import pandas as pd
import streamlit as st

from src.solvers import SOLVER_FAMILY
from ui.labels import label


def markdown_table(df, decimals=2):
    """A markdown table (the column names may contain LaTeX -- st.dataframe cannot render it)."""
    cell = lambda v: "–" if v is None or v != v else (f"{v:.{decimals}f}" if isinstance(v, float) else str(v))
    lines = ["| " + " | ".join(df.columns) + " |", "|" + "---|" * len(df.columns)]
    lines += ["| " + " | ".join(cell(v) for v in row) + " |" for row in df.itertuples(index=False)]
    return "\n".join(lines)


def summary_table(results):
    """One row per scenario and solver: Q, outlet temperatures, pinch, k, iterations and time."""
    rows = []
    for key in ("ambient", "precooled"):
        sc = results.get(key)
        if not sc:
            continue
        for key, r in sc["result"].results.items():
            rows.append({"Scenario": sc["label"], "Solver": SOLVER_FAMILY[key],
                         r"$\dot{Q}$ [kW]": r.Q_dot / 1000, r"$T_{c,o}$ [°C]": r.T_c_o, r"$T_{a,o}$ [°C]": r.T_a_o,
                         r"$\Delta T_{pinch}$ [K]": r.T_c_o - results["ambient"]["ops"].T_a_i, r"$k$ [W/m²K]": r.k,
                         "Iterations": len(r.history_hot), "Time [ms]": round((r.solve_time or 0) * 1e3, 1)})
    return pd.DataFrame(rows)


def economics_table(results):
    rows = [{"Scenario": sc["label"], r"$P_{pump}$ [W]": sc["P_pump"], r"$P_{fan}$ [W]": sc["P_fan"],
             r"$\dot{m}_{ev}$ [g/s]": None if sc["m_dot_ev"] is None else sc["m_dot_ev"] * 1000}
            for sc in (results.get(k) for k in ("ambient", "precooled")) if sc]
    return pd.DataFrame(rows)


def show_output(mode, output):
    """Draws a finished run on the page: the changed parameters, the tables, the log and the figures of the selected mode."""
    results, log_text, figures, changed = output
    if changed:
        st.caption("Changed vs. parameters.py: " + "; ".join(f"{label(n)} = {v!r}" for n, v in changed.items()))
    if mode == "solver":
        st.subheader("Results")
        st.markdown(markdown_table(summary_table(results)))
        st.subheader("Economics")
        st.markdown(markdown_table(economics_table(results)))
    else:
        st.subheader("Summary")
        st.code(log_text or "(no output)", language=None)
    if mode == "solver":
        with st.expander("Terminal output", expanded=False):
            st.code(log_text or "(no output)", language=None)
    if figures:
        for tab, (title, fig) in zip(st.tabs([title for title, _ in figures]), figures):
            with tab:
                st.pyplot(fig, width="stretch")
    else:
        st.caption("No figures -- switch on figures / benchmarks / analyses in the sidebar.")
