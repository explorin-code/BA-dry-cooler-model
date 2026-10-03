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

def calc_w_e(w_f: float, Ao_Ae_ratio: float) -> float:
    """Effective (minimum free cross-section) air velocity [m/s]."""
    return w_f * Ao_Ae_ratio


def calc_w_e_T(w_f: float, Ao_Ae_ratio: float, T_mean: float, T_in: float) -> float:
    """Effective air velocity, corrected for thermal expansion between
    inlet and mean bulk temperature [m/s]."""
    w_e = calc_w_e(w_f, Ao_Ae_ratio)
    return w_e * ((273.15 + T_mean) / (273.15 + T_in))


def calc_Re_air(d: float, w_e_T: float, rho_air: float, mu_air: float) -> float:
    """Air-side Reynolds number, based on diameter d (tube outer diameter
    d_a for heat transfer, fin-collar diameter d_c for Wang's friction factor)."""
    return (d * w_e_T * rho_air) / mu_air


# -----------------------------------------------------------------------
# Source: VDI Wärmeatlas, Chapter M1, §4, Gl. (16) [staggered, n>=4 rows,
# C=0.38] and Gl. (18) [staggered, n=1-3 rows: C=0.33 for n=2, C=0.36 for
# n=3]. Confirmed 2026-08 against the 12th ed. -- the row-count prefactor
# is the book's own C, not an invented correction.
# Book's stated validity range: 10^3 < Re_d < 10^5, 5 <= A/A_Go <= 30
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
    surface [W/m²K]."""
    return (Nu_air * lambda_air) / d


# -----------------------------------------------------------------------
# STATUS: INACTIVE -- not currently wired into calc_overall_k(), which
# uses calc_fin_efficiency_staggered() below instead. Kept here, ready
# to swap in (same signature) if the tube layout ever changes from
# staggered to in-line rows.
#
# Source: VDI Heat Atlas, Section M1, p. 1687, Eq. (13)   [in-line rows]
# -----------------------------------------------------------------------
def calc_fin_efficiency_inline(P_t: float, P_l: float, d_a: float, alpha_R: float, lambda_f: float, delta_R: float) -> float:
    # Determine bR and lR such that lR >= bR
    bR = min(P_t, P_l)
    lR = max(P_t, P_l)

    # Source: VDI M1, p. 1687, Eq. (13)
    phi_0 = 1.28 * (bR / d_a) * ((lR / bR) - 0.2)**0.5

    # Source: VDI M1, p. 1687, Eq. (12)
    phi = (phi_0 - 1) * (1 + 0.35 * log(phi_0))

    # Source: VDI M1, p. 1687, Eq. (7) [eta_R = tanh(X)/X] and Eq. (8) [X]
    X = phi * d_a / 2 * ((2 * alpha_R) / (lambda_f * delta_R))**0.5
    return tanh(X) / X


# -----------------------------------------------------------------------
# Source: VDI Heat Atlas, Section M1, p. 1687, Eq. (14)   [staggered rows]
# (book notation s_1/s_2 == P_t/P_l here)
# (shares the phi_0 -> phi -> X -> tanh(X)/X chain with Eq. (12), see
# calc_fin_efficiency_inline() above for the Eq. (13) in-line variant)
# -----------------------------------------------------------------------
def calc_fin_efficiency_staggered(P_t: float, P_l: float, d_a: float, alpha_R: float, lambda_f: float, delta_R: float) -> float:
    if P_l >= P_t / 2:
        bR = P_t
    else:
        bR = 2 * P_l

    lR = (P_l**2 + (P_t / 2)**2)**0.5

    # Source: VDI M1, p. 1687, Eq. (14)
    phi_0 = 1.27 * (bR / d_a) * ((lR / bR) - 0.3)**0.5

    # Source: VDI M1, p. 1687, Eq. (12)
    phi = (phi_0 - 1) * (1 + 0.35 * log(phi_0))

    # Source: VDI M1, p. 1687, Eq. (7) [eta_R = tanh(X)/X] and Eq. (8) [X]
    X = phi * d_a / 2 * ((2 * alpha_R) / (lambda_f * delta_R))**0.5
    return tanh(X) / X


def calc_alpha_S(alpha_R: float, eta_R: float, A: float, A_R: float) -> float:
    """Air-side heat transfer coefficient corrected for fin efficiency,
    referred to the total outer surface [W/m²K]."""
    return alpha_R * (1 - (1 - eta_R) * (A_R / A))


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
    """Coolant-side (tube) Reynolds number."""
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
    """Coolant-side heat transfer coefficient [W/m²K]."""
    Nu, Re = calc_Nu_coolant(w, d_i, l, rho_coolant, cool_mu, cool_Pr)
    return (Nu * cool_lambda) / d_i


# =============================================================================
# 3. ORCHESTRATORS -- the main API surface used by the solvers
# =============================================================================

def _compute_air_side(geo, ops, air_state, T_air_out: float):
    """Air-side intermediates shared by calc_overall_k and calc_diagnostics:
    thermally-corrected effective velocity, Nusselt number, and the
    resulting heat transfer coefficient referred to the bare tube surface."""
    T_air_mean = (ops.T_air_in + T_air_out) / 2.0

    w_e_T = calc_w_e_T(
        w_f=ops.w_f,
        Ao_Ae_ratio=geo.Ao_Ae_ratio,
        T_mean=T_air_mean,
        T_in=ops.T_air_in,
    )

    Nu_air = calc_Nu_air(
        A_ratio=(geo.A / geo.A_Go),
        Pr_air=air_state.Pr,
        d=geo.d_a,
        n=geo.n_rows,
        w_e_T=w_e_T,
        rho_air=air_state.rho,
        mu_air=air_state.mu,
    )

    alpha_R = calc_alpha_R(Nu_air, air_state.lambda_, geo.d_a)
    return w_e_T, Nu_air, alpha_R


def calc_overall_k(geo, ops, coolant_state, air_state, T_air_out: float) -> float:
    """Overall heat transfer coefficient k [W/m²K], referred to the outer
    (air-side) surface. Called once per solver iteration with the
    current guess for T_air_out."""

    # --- Air side (outer) ---------------------------------------------
    w_e_T, Nu_air, alpha_R = _compute_air_side(geo, ops, air_state, T_air_out)
    eta_R = calc_fin_efficiency_staggered(geo.P_t, geo.P_l, geo.d_a, alpha_R, geo.lambda_f, geo.delta_R)
    alpha_S = calc_alpha_S(alpha_R, eta_R, geo.A, geo.A_R)

    # --- Coolant side (inner) ------------------------------------------
    alpha_i = calc_alpha_i(
        w=ops.u_i,
        d_i=geo.d_i,
        l=geo.l,
        rho_coolant=coolant_state.rho,
        cool_mu=coolant_state.mu,
        cool_Pr=coolant_state.Pr,
        cool_lambda=coolant_state.lambda_,
    )

    # --- Combine into overall k -----------------------------------------
    k_inv = (1 / alpha_S) + (geo.A / geo.A_i) * ((1 / alpha_i) + (geo.d_a - geo.d_i) / (2 * geo.lambda_p))

    return k_inv ** (-1)


def calc_diagnostics(geo, ops, coolant_state, air_state, T_air_out: float) -> dict:
    """
    Recomputes the dimensionless groups (Pr, Re, Nu) for both sides using
    the SAME states/T_air_out that were fed into calc_overall_k for a given
    iteration. Intended to be called once more after convergence, using the
    final converged states, to report final Re/Nu/Pr alongside the outlet
    temperatures and heat transfer rate.
    """
    # --- Air side --------------------------------------------------------
    w_e_T, Nu_air, alpha_R = _compute_air_side(geo, ops, air_state, T_air_out)
    Re_air = calc_Re_air(geo.d_a, w_e_T, air_state.rho, air_state.mu)

    # --- Coolant side ------------------------------------------------------
    Nu_coolant, Re_coolant = calc_Nu_coolant(
        w=ops.u_i,
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

# NOTE: an earlier capacity-flow-ratio helper chain (calc_heat_cap_flow_*,
# calc_cap_flow_ratio) was removed here as dead code (unreferenced; NTU's
# solver computes R1=W1/W2 inline instead). It used geo.d_i -- a coolant-side
# dimension -- for both sides, which would have been a bug if ever revived.

