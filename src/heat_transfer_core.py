"""
heat_transfer_core.py
======================
All heat-transfer correlations (Nu/Re/alpha on both sides, fin
efficiency) plus the two orchestrator functions that combine them into
an overall k-value: calc_overall_k() (used inside the iteration loop)
and calc_diagnostics() (used once after convergence to report the
final Re/Nu/Pr/alpha numbers).

Layout: air-side helpers (including the inactive-but-ready in-line fin
variant), then coolant-side helpers, then the two orchestrators.
"""

from math import log10, log, tanh


# =============================================================================
# 1. AIR-SIDE HELPERS
# =============================================================================

def calc_w_e(w_fr: float, Afr_Ae_ratio: float) -> float:
    """Effective (minimum free cross-section) air velocity [m/s],
    w_e = w_0 * A_0/A_e. Source: VDI-Wärmeatlas (2019), Chapter M1, worked example p. 1689."""
    return w_fr * Afr_Ae_ratio


def calc_w_e_T(w_fr: float, Afr_Ae_ratio: float, theta_m: float, theta_i: float) -> float:
    """Effective air velocity, corrected for thermal expansion between
    inlet and mean bulk temperature [m/s]: w_eT = w_e * T_m / T_i (Kelvin).
    Source: VDI-Wärmeatlas (2019), Chapter M1, worked example p. 1689."""
    w_e = calc_w_e(w_fr, Afr_Ae_ratio)
    T_m = theta_m + 273.15                 # mean air temperature [K]
    T_i = theta_i + 273.15                 # inlet air temperature [K]
    return w_e * (T_m / T_i)


def calc_Re_air(d: float, w_e_T: float, rho_air: float, mu_air: float) -> float:
    """Air-side Reynolds number, based on diameter d (tube outer diameter
    d for heat transfer, fin-collar diameter D_c for Wang's friction factor).
    Source: VDI-Wärmeatlas (2019), Chapter M1, worked example p. 1689 (Re_d = d w_eT rho / eta); Wang, Chi & Chang (2000),
    Part II, Nomenclature (Re_Dc)."""
    return (d * w_e_T * rho_air) / mu_air


# -----------------------------------------------------------------------
# Source: VDI Wärmeatlas, Chapter M1, §4, Gl. (16) [staggered, n>=4 rows,
# C=0.38] and Gl. (18) [staggered, n=1-3 rows: C=0.33 for n=2, C=0.36 for
# n=3]. Confirmed 2026-08 against the 12th ed. -- the row-count prefactor
# is the book's own C, not an invented correction.
# Book's stated validity range: 10^3 < Re_d < 10^5, 5 <= A/A_p0 <= 30
# (not currently enforced/warned on here).
# -----------------------------------------------------------------------
def calc_Nu_air(A_ratio: float, Pr_air: float, d: float, n: int, w_e_T: float, rho_air: float, mu_air: float) -> float:
    """Air-side Nusselt number for the fin-tube bundle."""
    Re_air = calc_Re_air(d, w_e_T, rho_air, mu_air)
    C: float

    # staggered
    if n <= 2: C = 0.33
    elif n == 3: C = 0.36
    else: C = 0.38

    # inline:
    # if n <= 1: C = 0.20
    # else: C = 0.22

    return C * (Re_air**0.6) * (A_ratio**(-0.15)) * (Pr_air**(1/3))


def calc_alpha_R(Nu_air: float, lambda_air: float, d: float) -> float:
    """Air-side heat transfer coefficient, referred to the outer (finned)
    surface [W/m²K]: alpha_R = Nu_d lambda / d. Source: VDI-Wärmeatlas (2019), Chapter M1, worked example p. 1689."""
    return (Nu_air * lambda_air) / d


# -----------------------------------------------------------------------
# STATUS: INACTIVE -- not currently wired into calc_overall_k(), which
# uses calc_fin_efficiency_staggered() below instead. Kept here, ready
# to swap in (same signature) if the tube layout ever changes from
# staggered to in-line rows.
#
# Source: VDI Heat Atlas, Section M1, p. 1687, Eq. (13)   [in-line rows]
# -----------------------------------------------------------------------
def calc_fin_efficiency_inline(P_t: float, P_l: float, d: float, alpha_R: float, lambda_f: float, delta_f: float) -> float:
    # Determine bR and lR such that lR >= bR
    bR = min(P_t, P_l)
    lR = max(P_t, P_l)

    # Source: VDI M1, p. 1687, Eq. (13)
    phi_0 = 1.28 * (bR / d) * ((lR / bR) - 0.2)**0.5

    # Source: VDI M1, p. 1687, Eq. (12)
    phi = (phi_0 - 1) * (1 + 0.35 * log(phi_0))

    # Source: VDI M1, p. 1687, Eq. (7) [eta_f = tanh(X)/X] and Eq. (8) [X]
    X = phi * d / 2 * ((2 * alpha_R) / (lambda_f * delta_f))**0.5
    return tanh(X) / X


# -----------------------------------------------------------------------
# Source: VDI Heat Atlas, Section M1, p. 1687, Eq. (14)   [staggered rows]
# (book notation s_1/s_2 == P_t/P_l here)
# (shares the phi_0 -> phi -> X -> tanh(X)/X chain with Eq. (12), see
# calc_fin_efficiency_inline() above for the Eq. (13) in-line variant)
# -----------------------------------------------------------------------
def calc_fin_efficiency_staggered(P_t: float, P_l: float, d: float, alpha_R: float, lambda_f: float, delta_f: float) -> float:
    if P_l >= P_t / 2:
        bR = P_t
    else:
        bR = 2 * P_l

    lR = (P_l**2 + (P_t / 2)**2)**0.5

    # Source: VDI M1, p. 1687, Eq. (14)
    phi_0 = 1.27 * (bR / d) * ((lR / bR) - 0.3)**0.5

    # Source: VDI M1, p. 1687, Eq. (12)
    phi = (phi_0 - 1) * (1 + 0.35 * log(phi_0))

    # Source: VDI M1, p. 1687, Eq. (7) [eta_f = tanh(X)/X] and Eq. (8) [X]
    X = phi * d / 2 * ((2 * alpha_R) / (lambda_f * delta_f))**0.5
    return tanh(X) / X


def calc_alpha_S(alpha_R: float, eta_f: float, A: float, A_f: float) -> float:
    """Air-side heat transfer coefficient corrected for fin efficiency,
    referred to the total outer surface [W/m²K]:
    alpha_S = alpha_R * [1 - (1 - eta_f) * A_f/A].
    Source: VDI-Wärmeatlas (2019), Chapter M1 (and worked example p. 1689)."""
    return alpha_R * (1 - (1 - eta_f) * (A_f / A))


# =============================================================================
# 2. COOLANT-SIDE HELPERS
# =============================================================================

# -----------------------------------------------------------------------
# Source: VDI Wärmeatlas, Chapter G1, §4.1, Gl. (26); friction factor psi
# from Gl. (27) ("Hanakov" relation, per the book's own attribution).
# Confirmed 2026-08 against the 12th ed. -- exact match, entrance-length
# term included.
# -----------------------------------------------------------------------
def calc_Nu_turbulent(Re: float, Pr: float, d_i: float, l: float) -> float:
    psi = (1.8 * log10(Re) - 1.5)**(-2)
    return (((psi / 8) * (Re - 1000) * Pr) / (1 + 12.7 * (psi / 8)**0.5 * (Pr**(2/3) - 1))) * (1 + (d_i / l)**(2/3))


# -----------------------------------------------------------------------
# Source: VDI Wärmeatlas, Chapter G1, §3.2.2, Gl. (25); components
# Gl. (17) [Nu_mq1=4.364], Gl. (18) [Nu_mq2], Gl. (24) [Nu_mq3]. Confirmed
# 2026-08 against the 12th ed. -- this is the CONSTANT-HEAT-FLUX boundary
# condition ("konstante Wärmestromdichte"), not the constant-wall-temperature
# variant (which uses 0.7/Nu1=3.66 instead of this formula's 0.6/4.364) --
# the two look superficially similar but are genuinely different formulas.
# -----------------------------------------------------------------------
def calc_Nu_laminar(Re: float, Pr: float, d_i: float, l: float) -> float:
    Nu_mq1 = 4.364
    Nu_mq2 = 1.953 * (Re * Pr * (d_i / l))**(1/3)
    Nu_mq3 = 0.924 * Pr**(1/3) * (Re * (d_i / l))**0.5
    return (Nu_mq1**3 + 0.6**3 + (Nu_mq2 - 0.6)**3 + Nu_mq3**3)**(1/3)


def calc_Re_coolant(w: float, d_i: float, rho_coolant: float, cool_mu: float) -> float:
    """Coolant-side (tube) Reynolds number Re = w d_i rho / eta.
    Source: VDI-Wärmeatlas (2019), Chapter G1."""
    return (d_i * w * rho_coolant) / cool_mu


def calc_Nu_coolant(w: float, d_i: float, l: float, rho_coolant: float, cool_mu: float, cool_Pr: float):
    """Returns (Nu, Re) for the coolant (tube) side, handling the
    laminar / transitional / turbulent blend.
    Source: VDI Wärmeatlas, Chapter G1, §4.2, Gl. (29)/(30) -- the 2300-4000
    band and linear blend are the book's own transition treatment (its
    section heading literally reads "... bei 2300 < Re < 4*10^3"), not an
    invented approximation. The book's own endpoint formulas (Gl. 34/37)
    are algebraically identical to evaluating calc_Nu_laminar(2300, ...)/
    calc_Nu_turbulent(4000, ...) directly, exactly as done below.
    NOTE: this transition band (2300-4000) is not the same one
    pressure_drop.calc_f_coolant uses for its own laminar/turbulent
    blend (2320-3000) -- confirmed genuinely different sources (this one is
    from the book; that one is the author's choice, see its own docstring),
    not an inconsistency to reconcile."""
    Re = calc_Re_coolant(w, d_i, rho_coolant, cool_mu)

    if Re < 2300:
        Nu = calc_Nu_laminar(Re, cool_Pr, d_i, l)
    elif 2300 <= Re < 4000:
        gamma = (Re - 2300) / (4000 - 2300)
        Nu = (1 - gamma) * calc_Nu_laminar(2300, cool_Pr, d_i, l) + gamma * calc_Nu_turbulent(4000, cool_Pr, d_i, l)
    else:
        Nu = calc_Nu_turbulent(Re, cool_Pr, d_i, l)

    return Nu, Re


def calc_alpha_i(w: float, d_i: float, l: float, rho_coolant: float, cool_mu: float, cool_Pr: float, cool_lambda: float) -> float:
    """Coolant-side heat transfer coefficient alpha_i = Nu lambda / d_i
    [W/m²K]. Source: VDI-Wärmeatlas (2019), Chapter G1."""
    Nu, Re = calc_Nu_coolant(w, d_i, l, rho_coolant, cool_mu, cool_Pr)
    return (Nu * cool_lambda) / d_i


# =============================================================================
# 3. ORCHESTRATORS -- the main API surface used by the solvers
# =============================================================================

def _compute_air_side(geo, ops, air_state, theta_a_o: float):
    """Air-side intermediates shared by calc_overall_k and calc_diagnostics:
    thermally-corrected effective velocity, Nusselt number, and the
    resulting heat transfer coefficient referred to the bare tube surface."""
    theta_a_m = (ops.theta_a_i + theta_a_o) / 2.0

    w_e_T = calc_w_e_T(
        w_fr=ops.w_fr,
        Afr_Ae_ratio=geo.Afr_Ae_ratio,
        theta_m=theta_a_m,
        theta_i=ops.theta_a_i,
    )

    Nu_air = calc_Nu_air(
        A_ratio=(geo.A / geo.A_p0),
        Pr_air=air_state.Pr,
        d=geo.d,
        n=geo.N_r,
        w_e_T=w_e_T,
        rho_air=air_state.rho,
        mu_air=air_state.mu,
    )

    alpha_R = calc_alpha_R(Nu_air, air_state.lambda_, geo.d)
    return w_e_T, Nu_air, alpha_R


def calc_overall_k(geo, ops, coolant_state, air_state, theta_a_o: float) -> float:
    """Overall heat transfer coefficient k [W/m²K], referred to the outer
    (air-side) surface: 1/k = 1/alpha_S + A/A_i * (1/alpha_i + (d - d_i)/(2 lambda_G)).
    Source: VDI-Wärmeatlas (2019), Chapter M1; VDI-Wärmeatlas (2019), Chapter M1, worked example p. 1689.
    Called once per solver iteration with the current guess for theta_a_o."""

    # --- Air side (outer) ---------------------------------------------
    w_e_T, Nu_air, alpha_R = _compute_air_side(geo, ops, air_state, theta_a_o)
    eta_f = calc_fin_efficiency_staggered(geo.P_t, geo.P_l, geo.d, alpha_R, geo.lambda_f, geo.delta_f)
    alpha_S = calc_alpha_S(alpha_R, eta_f, geo.A, geo.A_f)

    # --- Coolant side (inner) ------------------------------------------
    alpha_i = calc_alpha_i(
        w=ops.w_c,
        d_i=geo.d_i,
        l=geo.l,
        rho_coolant=coolant_state.rho,
        cool_mu=coolant_state.mu,
        cool_Pr=coolant_state.Pr,
        cool_lambda=coolant_state.lambda_,
    )

    # --- Combine into overall k -----------------------------------------
    k_inv = (1 / alpha_S) + (geo.A / geo.A_i) * ((1 / alpha_i) + (geo.d - geo.d_i) / (2 * geo.lambda_p))

    return k_inv ** (-1)


def calc_diagnostics(geo, ops, coolant_state, air_state, theta_a_o: float) -> dict:
    """
    Recomputes the dimensionless groups (Pr, Re, Nu) for both sides using
    the SAME states/theta_a_o that were fed into calc_overall_k for a given
    iteration. Intended to be called once more after convergence, using the
    final converged states, to report final Re/Nu/Pr alongside the outlet
    temperatures and heat transfer rate.
    """
    # --- Air side --------------------------------------------------------
    w_e_T, Nu_air, alpha_R = _compute_air_side(geo, ops, air_state, theta_a_o)
    Re_air = calc_Re_air(geo.d, w_e_T, air_state.rho, air_state.mu)

    # --- Coolant side ------------------------------------------------------
    Nu_coolant, Re_coolant = calc_Nu_coolant(
        w=ops.w_c,
        d_i=geo.d_i,
        l=geo.l,
        rho_coolant=coolant_state.rho,
        cool_mu=coolant_state.mu,
        cool_Pr=coolant_state.Pr,
    )
    alpha_i = (Nu_coolant * coolant_state.lambda_) / geo.d_i

    # NOTE on naming: all three follow VDI Chapter M1's own notation (as in
    # its worked example and k equation): alpha_i = coolant side, inside
    # the tube; alpha_R / alpha_S = air side, bare vs. fin-efficiency-
    # corrected outer coefficient.
    return {
        'Pr_air': air_state.Pr,
        'Re_air': Re_air,
        'Nu_air': Nu_air,
        'alpha_R': alpha_R,
        'Pr_coolant': coolant_state.Pr,
        'Re_coolant': Re_coolant,
        'Nu_coolant': Nu_coolant,
        'alpha_i': alpha_i,
    }



