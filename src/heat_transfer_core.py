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

def calc_w_e(w_fr: float, sigma: float) -> float:
    """Effective (minimum free cross-section) air velocity [m/s],
    w_e = w_0 * A_0/A_e. Source: VDI-Wärmeatlas (2019), Chapter M1, worked example p. 1689."""
    return w_fr / sigma


def calc_w_e_m(w_fr: float, sigma: float, T_m: float, T_i: float) -> float:
    """Effective air velocity, corrected for thermal expansion between
    inlet and mean bulk temperature [m/s] (VDI M1: w_eT): w_e_m = w_e * T_m / T_i (Kelvin).
    Source: VDI-Wärmeatlas (2019), Chapter M1, worked example p. 1689."""
    w_e = calc_w_e(w_fr, sigma)
    T_m_K = T_m + 273.15               # mean air temperature [K]
    T_i_K = T_i + 273.15               # inlet air temperature [K]
    return w_e * (T_m_K / T_i_K)


def calc_Re_a(d: float, w_e_m: float, rho_a: float, mu_a: float) -> float:
    """Air-side Reynolds number, based on diameter d (tube outer diameter
    d for heat transfer, fin-collar diameter d_c for Wang's friction factor).
    Source: VDI-Wärmeatlas (2019), Chapter M1, worked example p. 1689 (Re_d = d w_eT rho / eta; here Re_a with w_e_m); Wang, Chi & Chang (2000),
    Part II, Nomenclature (Re_d_c)."""
    return (d * w_e_m * rho_a) / mu_a


# -----------------------------------------------------------------------
# Source: VDI Wärmeatlas, Chapter M1, §4, Gl. (16) [staggered, n>=4 rows,
# C=0.38] and Gl. (18) [staggered, n=1-3 rows: C=0.33 for n=2, C=0.36 for
# n=3]. Confirmed 2026-08 against the 12th ed. -- the row-count prefactor
# is the book's own C (here C_r), not an invented correction.
# Book's stated validity range: 10^3 < Re_d < 10^5, 5 <= A/A_p0 <= 30
# (not currently enforced/warned on here).
# -----------------------------------------------------------------------
def calc_Nu_a(A_ratio: float, Pr_a: float, d: float, n: int, w_e_m: float, rho_a: float, mu_a: float) -> float:
    """Air-side Nusselt number for the fin-tube bundle."""
    Re_a = calc_Re_a(d, w_e_m, rho_a, mu_a)
    C_r: float

    # staggered
    if n <= 2: C_r = 0.33
    elif n == 3: C_r = 0.36
    else: C_r = 0.38

    # inline:
    # if n <= 1: C_r = 0.20
    # else: C_r = 0.22

    return C_r * (Re_a**0.6) * (A_ratio**(-0.15)) * (Pr_a**(1/3))


def calc_alpha_a(Nu_a: float, lambda_a: float, d: float) -> float:
    """Air-side heat transfer coefficient, referred to the outer (finned)
    surface [W/m²K] (VDI M1: alpha_R): alpha_a = Nu_a lambda_a / d. Source: VDI-Wärmeatlas (2019), Chapter M1, worked example p. 1689."""
    return (Nu_a * lambda_a) / d


# -----------------------------------------------------------------------
# STATUS: INACTIVE -- not currently wired into calc_overall_k(), which
# uses calc_fin_efficiency_staggered() below instead. Kept here, ready
# to swap in (same signature) if the tube layout ever changes from
# staggered to in-line rows; the fin rectangle then comes from the
# in-line rule b_f = min(s_t, s_l), l_f = max(s_t, s_l) (Geometry would
# need a layout switch for that).
#
# Source: VDI Heat Atlas, Section M1, p. 1687, Eq. (13)   [in-line rows]
# -----------------------------------------------------------------------
def calc_fin_efficiency_inline(b_f: float, l_f: float, d: float, alpha_a: float, lambda_f: float, delta_f: float) -> float:
    # Source: VDI M1, p. 1687, Eq. (13)
    phi_f_0 = 1.28 * (b_f / d) * ((l_f / b_f) - 0.2)**0.5

    # Source: VDI M1, p. 1687, Eq. (12)
    phi_f = (phi_f_0 - 1) * (1 + 0.35 * log(phi_f_0))

    # Source: VDI M1, p. 1687, Eq. (7) [eta_f = tanh(X)/X] and Eq. (8) [X]
    X = phi_f * d / 2 * ((2 * alpha_a) / (lambda_f * delta_f))**0.5
    return tanh(X) / X


# -----------------------------------------------------------------------
# Source: VDI Heat Atlas, Section M1, p. 1687, Eq. (14)   [staggered rows]
# (fin rectangle b_f, l_f from Geometry.b_f / Geometry.l_f)
# (shares the phi_f_0 -> phi_f -> X -> tanh(X)/X chain with Eq. (12), see
# calc_fin_efficiency_inline() above for the Eq. (13) in-line variant)
# -----------------------------------------------------------------------
def calc_fin_efficiency_staggered(b_f: float, l_f: float, d: float, alpha_a: float, lambda_f: float, delta_f: float) -> float:
    # Source: VDI M1, p. 1687, Eq. (14)
    phi_f_0 = 1.27 * (b_f / d) * ((l_f / b_f) - 0.3)**0.5

    # Source: VDI M1, p. 1687, Eq. (12)
    phi_f = (phi_f_0 - 1) * (1 + 0.35 * log(phi_f_0))

    # Source: VDI M1, p. 1687, Eq. (7) [eta_f = tanh(X)/X] and Eq. (8) [X]
    X = phi_f * d / 2 * ((2 * alpha_a) / (lambda_f * delta_f))**0.5
    return tanh(X) / X


def calc_alpha_a_eff(alpha_a: float, eta_f: float, A: float, A_f: float) -> float:
    """Effective air-side heat transfer coefficient, corrected for fin
    efficiency and referred to the total outer surface [W/m²K] (VDI M1:
    "scheinbarer" alpha_S):
    alpha_a_eff = alpha_a * [1 - (1 - eta_f) * A_f/A].
    Source: VDI-Wärmeatlas (2019), Chapter M1 (and worked example p. 1689)."""
    return alpha_a * (1 - (1 - eta_f) * (A_f / A))


# =============================================================================
# 2. COOLANT-SIDE HELPERS
# =============================================================================

# -----------------------------------------------------------------------
# Source: VDI Wärmeatlas, Chapter G1, §4.1, Gl. (26); friction factor psi
# from Gl. (27) ("Hanakov" relation, per the book's own attribution).
# Confirmed 2026-08 against the 12th ed. -- exact match, entrance-length
# term included.
# -----------------------------------------------------------------------
def calc_Nu_c_turb(Re: float, Pr: float, d_i: float, l: float) -> float:
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
def calc_Nu_c_lam(Re: float, Pr: float, d_i: float, l: float) -> float:
    Nu_mq1 = 4.364
    Nu_mq2 = 1.953 * (Re * Pr * (d_i / l))**(1/3)
    Nu_mq3 = 0.924 * Pr**(1/3) * (Re * (d_i / l))**0.5
    return (Nu_mq1**3 + 0.6**3 + (Nu_mq2 - 0.6)**3 + Nu_mq3**3)**(1/3)


def calc_Re_c(w: float, d_i: float, rho_c: float, mu_c: float) -> float:
    """Coolant-side (tube) Reynolds number Re = w d_i rho / eta.
    Source: VDI-Wärmeatlas (2019), Chapter G1."""
    return (d_i * w * rho_c) / mu_c


def calc_Nu_c(w: float, d_i: float, l: float, rho_c: float, mu_c: float, Pr_c: float):
    """Returns (Nu, Re) for the coolant (tube) side, handling the
    laminar / transitional / turbulent blend.
    Source: VDI Wärmeatlas, Chapter G1, §4.2, Gl. (29)/(30) -- the 2300-4000
    band and linear blend are the book's own transition treatment (its
    section heading literally reads "... bei 2300 < Re < 4*10^3"), not an
    invented approximation. The book's own endpoint formulas (Gl. 34/37)
    are algebraically identical to evaluating calc_Nu_c_lam(2300, ...)/
    calc_Nu_c_turb(4000, ...) directly, exactly as done below.
    NOTE: this transition band (2300-4000) is not the same one
    pressure_drop.calc_f_coolant uses for its own laminar/turbulent
    blend (2320-3000) -- confirmed genuinely different sources (this one is
    from the book; that one is the author's choice, see its own docstring),
    not an inconsistency to reconcile."""
    Re = calc_Re_c(w, d_i, rho_c, mu_c)

    if Re < 2300:
        Nu = calc_Nu_c_lam(Re, Pr_c, d_i, l)
    elif 2300 <= Re < 4000:
        gamma = (Re - 2300) / (4000 - 2300)
        Nu = (1 - gamma) * calc_Nu_c_lam(2300, Pr_c, d_i, l) + gamma * calc_Nu_c_turb(4000, Pr_c, d_i, l)
    else:
        Nu = calc_Nu_c_turb(Re, Pr_c, d_i, l)

    return Nu, Re


def calc_alpha_c(w: float, d_i: float, l: float, rho_c: float, mu_c: float, Pr_c: float, lambda_c: float) -> float:
    """Coolant-side heat transfer coefficient alpha_c = Nu lambda / d_i
    [W/m²K]. Source: VDI-Wärmeatlas (2019), Chapter G1."""
    Nu, Re = calc_Nu_c(w, d_i, l, rho_c, mu_c, Pr_c)
    return (Nu * lambda_c) / d_i


# =============================================================================
# 3. ORCHESTRATORS -- the main API surface used by the solvers
# =============================================================================

def _compute_air_side(geo, ops, air_state, T_a_m: float):
    """Air-side intermediates shared by calc_overall_k and calc_diagnostics:
    thermally-corrected effective velocity, Nusselt number, and the
    resulting heat transfer coefficient referred to the bare tube surface."""

    w_e_m = calc_w_e_m(
        w_fr=ops.w_fr,
        sigma=geo.sigma,
        T_m=T_a_m,
        T_i=ops.T_a_i,
    )

    Nu_a = calc_Nu_a(
        A_ratio=(geo.A / geo.A_p0),
        Pr_a=air_state.Pr,
        d=geo.d,
        n=geo.N_r,
        w_e_m=w_e_m,
        rho_a=air_state.rho,
        mu_a=air_state.mu,
    )

    alpha_a = calc_alpha_a(Nu_a, air_state.lambda_, geo.d)
    return w_e_m, Nu_a, alpha_a


def calc_overall_k(geo, ops, coolant_state, air_state, T_a_m: float) -> float:
    """Overall heat transfer coefficient k [W/m²K], referred to the outer
    (air-side) surface: 1/k = 1/alpha_a_eff + A/A_c * (1/alpha_c + delta_p/lambda_p), with the
    tube wall thickness delta_p = Geometry.delta_p.
    Source: VDI-Wärmeatlas (2019), Chapter M1; VDI-Wärmeatlas (2019), Chapter M1, worked example p. 1689.
    T_a_m: mean air temperature (thesis Eq. mean_temperatures) -- of the
    whole cooler for LMTD/NTU, of the cell for Cell; sets the velocity
    correction w_e_m (the air properties come in via air_state)."""

    # --- Air side (outer) ---------------------------------------------
    w_e_m, Nu_a, alpha_a = _compute_air_side(geo, ops, air_state, T_a_m)
    eta_f = calc_fin_efficiency_staggered(geo.b_f, geo.l_f, geo.d, alpha_a, geo.lambda_f, geo.delta_f)
    alpha_a_eff = calc_alpha_a_eff(alpha_a, eta_f, geo.A, geo.A_f)

    # --- Coolant side (inner) ------------------------------------------
    alpha_c = calc_alpha_c(
        w=ops.w_c,
        d_i=geo.d_i,
        l=geo.l_tot,
        rho_c=coolant_state.rho,
        mu_c=coolant_state.mu,
        Pr_c=coolant_state.Pr,
        lambda_c=coolant_state.lambda_,
    )

    # --- Combine into overall k -----------------------------------------
    k_inv = (1 / alpha_a_eff) + (geo.A / geo.A_c) * ((1 / alpha_c) + geo.delta_p / geo.lambda_p)

    return k_inv ** (-1)


def calc_diagnostics(geo, ops, coolant_state, air_state, T_a_m: float) -> dict:
    """
    Recomputes the dimensionless groups (Pr, Re, Nu) for both sides using
    the SAME states/T_a_m that were fed into calc_overall_k for a given
    iteration. Intended to be called once more after convergence, using the
    final converged states, to report final Re/Nu/Pr alongside the outlet
    temperatures and heat transfer rate.
    """
    # --- Air side --------------------------------------------------------
    w_e_m, Nu_a, alpha_a = _compute_air_side(geo, ops, air_state, T_a_m)
    Re_a = calc_Re_a(geo.d, w_e_m, air_state.rho, air_state.mu)

    # --- Coolant side ------------------------------------------------------
    Nu_c, Re_c = calc_Nu_c(
        w=ops.w_c,
        d_i=geo.d_i,
        l=geo.l_tot,
        rho_c=coolant_state.rho,
        mu_c=coolant_state.mu,
        Pr_c=coolant_state.Pr,
    )
    alpha_c = (Nu_c * coolant_state.lambda_) / geo.d_i

    # Naming (thesis notation): index a = air, c = coolant. VDI M1 writes
    # alpha_i (coolant), alpha_R (air, from Nu) and alpha_S (air, incl. fin
    # efficiency) -- here alpha_c, alpha_a, alpha_a_eff.
    return {
        'Pr_a': air_state.Pr,
        'Re_a': Re_a,
        'Nu_a': Nu_a,
        'alpha_a': alpha_a,
        'Pr_c': coolant_state.Pr,
        'Re_c': Re_c,
        'Nu_c': Nu_c,
        'alpha_c': alpha_c,
    }



