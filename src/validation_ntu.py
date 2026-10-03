"""
validation_ntu.py
==================
Validates the NTU solver's element-field core (solvers.solve_ntu_field)
against the published results of its source, using exactly the same
function the solver runs -- in dimensionless form (T_h,in = 1, T_c,in = 0,
C_h = 1), constant U and cp, one tube per row (N_t = 1).

Source: Cabezas-Gomez, Navarro & Saiz-Jabardo (2007), "Thermal Performance
of Multipass Parallel and Counter-Cross-Flow Heat Exchangers", J. Heat
Transfer 129, pp. 282-289, DOI 10.1115/1.2430719.

Checks:
  1. Closed-form transcription: Table 1 Eq. (16) against the paper's own
     Table 6 (four-pass column) -- catches transcription errors before
     the closed forms are used as references.
  2. Closed forms: G_(1,1) and G^c_(2,1)..G^c_(4,1), Table 1 Eqs. (9),
     (14)-(16), p. 285, over a grid of C*, NTU and both C_min sides.
  3. Regression: Table 5, pp. 287-288 (5-10 passes), used only at its
     listed (R, NTU) points, never interpolated.

Run:  python -m src.validation_ntu
"""

from math import exp
import sys

from src.solvers import solve_ntu_field


# =============================================================================
# Closed-form references -- Table 1, p. 285
# =============================================================================
# Table 1 notation: Fluid A mixed (= tube fluid = coolant), Fluid B unmixed
# (= external fluid = air); C*_A = C_A/C_B (may exceed 1), NTU_A = UA/C_A,
# eps_A = tube-side effectiveness (T_h,i - T_h,o)/(T_h,i - T_c,i).

def eps_A_G11(NTU_A: float, C_A: float) -> float:
    """Pure cross-flow, one row. Source: Cabezas-Gomez, Navarro &
    Saiz-Jabardo (2007), J. Heat Transfer 129, p. 285, Table 1, Eq. (9)."""
    return 1 - exp(-(1 - exp(-NTU_A * C_A)) / C_A)


def eps_A_Gc21(NTU_A: float, C_A: float) -> float:
    """Two-pass counter-cross-flow. Source: Cabezas-Gomez, Navarro &
    Saiz-Jabardo (2007), J. Heat Transfer 129, p. 285, Table 1, Eq. (14)."""
    K = 1 - exp(-NTU_A * (C_A / 2))
    return 1 - 1 / (K / 2 + (1 - K / 2) * exp(2 * K / C_A))


def eps_A_Gc31(NTU_A: float, C_A: float) -> float:
    """Three-pass counter-cross-flow. Source: Cabezas-Gomez, Navarro &
    Saiz-Jabardo (2007), J. Heat Transfer 129, p. 285, Table 1, Eq. (15)."""
    K = 1 - exp(-NTU_A * (C_A / 3))
    return 1 - 1 / ((1 - K / 2)**2 * exp(3 * K / C_A)
                    + (K * (1 - K / 4) - (1 - K / 2) * K**2 / C_A) * exp(K / C_A))


def eps_A_Gc41(NTU_A: float, C_A: float) -> float:
    """Four-pass counter-cross-flow. Source: Cabezas-Gomez, Navarro &
    Saiz-Jabardo (2007), J. Heat Transfer 129, p. 285, Table 1, Eq. (16)."""
    K = 1 - exp(-NTU_A * C_A / 4)
    return 1 - 1 / ((K / 2) * (1 - K / 2 + K**2 / 4)
                    + K * (1 - K / 2) * (1 - 2 * (K / C_A) * (1 - K / 2)) * exp(2 * K / C_A)
                    + (1 - K / 2)**3 * exp(4 * K / C_A))


CLOSED_FORMS = {1: eps_A_G11, 2: eps_A_Gc21, 3: eps_A_Gc31, 4: eps_A_Gc41}


def closed_form_eps(n_rows: int, C_star: float, NTU: float, C_min_side: str) -> float:
    """Exchanger effectiveness eps = q/q_max (based on C_min) from the
    Table 1 closed forms. C_star = C_min/C_max, NTU = UA/C_min (p. 289,
    Nomenclature); C_min_side = 'air' or 'tube'. Uses the Table 1 footnote
    conversions eps_B = eps_A C*_A, NTU_B = NTU_A C*_A, C*_A = 1/C*_B."""
    if C_min_side == 'tube':
        C_A, NTU_A = C_star, NTU                     # C_A = C_min
        return CLOSED_FORMS[n_rows](NTU_A, C_A)
    C_A, NTU_A = 1 / C_star, NTU * C_star            # C_A = C_max
    return CLOSED_FORMS[n_rows](NTU_A, C_A) * C_A    # eps_B, B = air = C_min


# =============================================================================
# Numerical model in the same dimensionless setting
# =============================================================================

def numerical_eps_and_P(n_rows: int, C_h: float, C_c: float, UA: float, n_elements: int):
    """Runs solvers.solve_ntu_field for one circuit (N_t = 1) and returns
    (eps, P): eps = q/(C_min (T_h,i - T_c,i)), P = air-side temperature
    effectiveness (p. 286, Eq. (7))."""
    _, _, _, T_h_out, T_c_out, _ = solve_ntu_field(
        n_rows, n_elements,
        C_h_circuit=C_h,
        C_c_element=C_c / n_elements,                # C_c^e = C_c/(N_e N_t), N_t = 1
        UA_element=UA / (n_rows * n_elements),       # (UA)^e = UA/(N_e N_t N_r)
        T_coolant_in=1.0, T_air_in=0.0,
        tol=1e-11, max_iter=100000,
    )
    q = C_h * (1.0 - T_h_out)
    return q / min(C_h, C_c), T_c_out


# =============================================================================
# Reference data
# =============================================================================

# Table 6, p. 288, "Four pass counter cross-flow" columns: (C*, NTU) ->
# (eps for C_min = C_air, eps for C_min = C_t). Transcription check only.
TABLE6_FOUR_PASS = {
    (0.2, 2): (0.8293, 0.8294),
    (0.2, 12): (0.9993, 0.9997),
    (1.0, 2): (0.6578, 0.6578),
    (1.0, 12): (0.8537, 0.8537),
}

# Table 5, pp. 287-288: P for counter-cross-flow, (R, NTU) -> P for
# [5, 6, 7, 8, 9, 10] passes. R = C_c/C_h (Eq. (8)), P = air-side
# temperature effectiveness (Eq. (7)), NTU = UA/C_min (Nomenclature; the
# R = 3 rows only match with C_min, not C_air).
TABLE5_PASSES = [5, 6, 7, 8, 9, 10]
TABLE5 = {
    (0.2, 2.0):  [0.8302, 0.8306, 0.8310, 0.8311, 0.8312, 0.8313],
    (0.2, 5.0):  [0.9832, 0.9839, 0.9842, 0.9845, 0.9847, 0.9848],
    (0.2, 10.0): [0.9992, 0.9994, 0.9995, 0.9996, 0.9996, 0.9996],
    (0.5, 2.0):  [0.7710, 0.7720, 0.7727, 0.7732, 0.7735, 0.7737],
    (0.5, 5.0):  [0.9477, 0.9506, 0.9524, 0.9535, 0.9543, 0.9549],
    (0.5, 10.0): [0.9882, 0.9913, 0.9930, 0.9940, 0.9946, 0.9950],
    (1.0, 1.0):  [0.4984, 0.4988, 0.4992, 0.4993, 0.4995, 0.4996],
    (1.0, 2.0):  [0.6609, 0.6626, 0.6637, 0.6644, 0.6648, 0.6652],
    (1.0, 5.0):  [0.8123, 0.8183, 0.8221, 0.8246, 0.8264, 0.8277],
    (1.0, 10.0): [0.8661, 0.8772, 0.8847, 0.8898, 0.8935, 0.8963],
    (1.5, 2.0):  [0.4902, 0.4911, 0.4916, 0.4920, 0.4922, 0.4924],
    (1.5, 5.0):  [0.6092, 0.6120, 0.6137, 0.6148, 0.6156, 0.6162],
    (1.5, 10.0): [0.6465, 0.6501, 0.6523, 0.6537, 0.6548, 0.6555],
    (3.0, 2.0):  [0.2683, 0.2685, 0.2687, 0.2688, 0.2689, 0.2689],
    (3.0, 5.0):  [0.3239, 0.3243, 0.3246, 0.3247, 0.3249, 0.3249],
    (3.0, 10.0): [0.3326, 0.3327, 0.3328, 0.3329, 0.3329, 0.3329],
    (7.0, 1.0):  [0.0875, 0.0875, 0.0875, 0.0875, 0.0875, 0.0875],
}


# =============================================================================
# Runner
# =============================================================================

def run_validation(n_elements: int = 100, tol_closed: float = 1e-4, tol_table: float = 1.5e-4) -> bool:
    """Returns True if every check passes. tol_table allows for the
    tables' 4-decimal rounding; n_elements = 100 keeps C_c^e/C_h^e small
    for every listed R (largest R = 7)."""
    ok = True

    print("1) Closed-form transcription vs. Table 6 (four passes)")
    for (C_star, NTU), refs in TABLE6_FOUR_PASS.items():
        for side, ref in zip(('air', 'tube'), refs):
            eps = closed_form_eps(4, C_star, NTU, side)
            passed = abs(eps - ref) <= 5e-5 + 1e-9
            ok &= passed
            print(f"   C*={C_star:<4} NTU={NTU:<3} C_min={side:<4}  Eq.(16)={eps:.5f}  Table 6={ref:.4f}  "
                  f"{'ok' if passed else 'FAIL'}")

    print(f"\n2) Numerical model vs. Table 1 closed forms (N_e = {n_elements})")
    worst = 0.0
    for n_rows in (1, 2, 3, 4):
        for C_star in (0.2, 0.5, 1.0):
            for NTU in (0.5, 2.0, 5.0):
                for side in ('air', 'tube'):
                    C_h, C_c = (1.0, 1.0 / C_star) if side == 'tube' else (1.0 / C_star, 1.0)
                    eps_num, _ = numerical_eps_and_P(n_rows, C_h, C_c, NTU * min(C_h, C_c), n_elements)
                    err = abs(eps_num - closed_form_eps(n_rows, C_star, NTU, side))
                    worst = max(worst, err)
                    if err > tol_closed:
                        ok = False
                        print(f"   FAIL n={n_rows} C*={C_star} NTU={NTU} C_min={side}: |d eps| = {err:.2e}")
    print(f"   {4 * 3 * 3 * 2} cases, max |eps_num - eps_closed| = {worst:.2e}  "
          f"({'ok' if worst <= tol_closed else 'FAIL'}, tol {tol_closed:.0e})")

    print(f"\n3) Numerical model vs. Table 5 regression points (N_e = {n_elements})")
    worst = 0.0
    for (R, NTU), refs in TABLE5.items():
        C_h, C_c = 1.0, R                            # R = C_c/C_h
        for n_rows, ref in zip(TABLE5_PASSES, refs):
            _, P = numerical_eps_and_P(n_rows, C_h, C_c, NTU * min(C_h, C_c), n_elements)
            err = abs(P - ref)
            worst = max(worst, err)
            if err > tol_table:
                ok = False
                print(f"   FAIL n={n_rows} R={R} NTU={NTU}: P = {P:.5f}, Table 5 = {ref:.4f}")
    print(f"   {len(TABLE5) * len(TABLE5_PASSES)} points, max |P_num - P_table| = {worst:.2e}  "
          f"({'ok' if worst <= tol_table else 'FAIL'}, tol {tol_table:.1e})")

    print(f"\nNTU validation {'PASSED' if ok else 'FAILED'}")
    return ok


if __name__ == "__main__":
    sys.exit(0 if run_validation() else 1)
