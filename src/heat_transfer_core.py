"""
heat_transfer_core.py
======================
All heat-transfer correlations (Nu/Re/alpha on both sides, fin
efficiency) plus the two orchestrator functions that combine them into
an overall k-value: calc_overall_k() (used inside the iteration loop)
and calc_diagnostics() (used once after convergence to report the
final Re/Nu/Pr/alpha numbers).

Layout: air-side helpers (including the inactive-but-ready in-line fin
variant), then coolant-side helpers, then the two orchestrators, then a
trash pile of currently-unused functions at the bottom.
"""

from math import log10, log, tanh


# =============================================================================
# 1. AIR-SIDE HELPERS
# =============================================================================

def calc_w_e(w_o: float, Ao_Ae_ratio: float) -> float:
    """Effective (minimum free cross-section) air velocity [m/s]."""
    return w_o * Ao_Ae_ratio


def calc_w_eT(w_o: float, Ao_Ae_ratio: float, T_mean: float, T_in: float) -> float:
    """Effective air velocity, corrected for thermal expansion between
    inlet and mean bulk temperature [m/s]."""
    w_e = calc_w_e(w_o, Ao_Ae_ratio)
    return w_e * ((273.15 + T_mean) / (273.15 + T_in))


def calc_Re_air(d: float, w_eT: float, rho_air: float, eta_air: float) -> float:
    """Air-side Reynolds number, based on tube outer diameter."""
    return (d * w_eT * rho_air) / eta_air


# -----------------------------------------------------------------------
# Source: [not given in original notes -- TODO: find & fill in]
# (row-count-dependent prefactor: 0.33 for n<=2 rows, 0.36 for n==3, 0.38
# for n>=4 rows -- reads like a VDI-style row-correction but unconfirmed)
# -----------------------------------------------------------------------
def calc_Nu_air(A_ratio: float, Pr_air: float, d: float, n: int, w_eT: float, rho_air: float, eta_air: float) -> float:
    """Air-side Nusselt number for the fin-tube bundle."""
    Re_air = calc_Re_air(d, w_eT, rho_air, eta_air)
    if n <= 2:
        return 0.33 * (Re_air**0.6) * (A_ratio**(-0.15)) * (Pr_air**(1/3))
    elif n == 3:
        return 0.36 * (Re_air**0.6) * (A_ratio**(-0.15)) * (Pr_air**(1/3))
    else:  # n >= 4
        return 0.38 * (Re_air**0.6) * (A_ratio**(-0.15)) * (Pr_air**(1/3))


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
def calc_fin_efficiency_inline(t_q: float, t_l: float, d: float, alpha_R: float, lambda_R: float, s: float) -> float:
    # Determine bR and lR such that lR >= bR
    bR = min(t_q, t_l)
    lR = max(t_q, t_l)

    # Source: VDI M1, p. 1687, Eq. (13)
    phi_0 = 1.28 * (bR / d) * ((lR / bR) - 0.2)**0.5

    # Source: VDI M1, p. 1687, Eq. (12)
    phi = (phi_0 - 1) * (1 + 0.35 * log(phi_0))

    X = phi * d / 2 * ((2 * alpha_R) / (lambda_R * s))**0.5
    return tanh(X) / X


# -----------------------------------------------------------------------
# Source: VDI Heat Atlas, Section M1, p. 1687, Eq. (14)   [staggered rows]
# (shares the phi_0 -> phi -> X -> tanh(X)/X chain with Eq. (12), see
# calc_fin_efficiency_inline() above for the Eq. (13) in-line variant)
# -----------------------------------------------------------------------
def calc_fin_efficiency_staggered(t_q: float, t_l: float, d: float, alpha_R: float, lambda_R: float, s: float) -> float:
    if t_l >= t_q / 2:
        bR = t_q
    else:
        bR = 2 * t_l

    lR = (t_l**2 + (t_q / 2)**2)**0.5

    # Source: VDI M1, p. 1687, Eq. (14)
    phi_0 = 1.27 * (bR / d) * ((lR / bR) - 0.3)**0.5

    # Source: VDI M1, p. 1687, Eq. (12)
    phi = (phi_0 - 1) * (1 + 0.35 * log(phi_0))

    X = phi * d / 2 * ((2 * alpha_R) / (lambda_R * s))**0.5
    return tanh(X) / X


def calc_alpha_S(alpha_R: float, eta_R: float, A: float, A_R: float) -> float:
    """Air-side heat transfer coefficient corrected for fin efficiency,
    referred to the total outer surface [W/m²K]."""
    return alpha_R * (1 - (1 - eta_R) * (A_R / A))


# =============================================================================
# 2. COOLANT-SIDE HELPERS
# =============================================================================

# -----------------------------------------------------------------------
# Source: [not given in original notes -- TODO: find & fill in]
# (Gnielinski-style turbulent correlation, incl. entrance-length term)
# -----------------------------------------------------------------------
def calc_Nu_turbulent(Re: float, Pr: float, d_i: float, l: float) -> float:
    psi = (1.8 * log10(Re) - 1.5)**(-2)
    return (((psi / 8) * (Re - 1000) * Pr) / (1 + 12.7 * (psi / 8)**0.5 * (Pr**(2/3) - 1))) * (1 + (d_i / l)**(2/3))


# -----------------------------------------------------------------------
# Source: [not given in original notes -- TODO: find & fill in]
# (combined laminar correlation -- constant term + two entrance-effect
# terms, cubic-mean blend; classic Gnielinski/Hausen-style form)
# -----------------------------------------------------------------------
def calc_Nu_laminar(Re: float, Pr: float, d_i: float, l: float) -> float:
    Nu_mq1 = 4.364
    Nu_mq2 = 1.953 * (Re * Pr * (d_i / l))**(1/3)
    Nu_mq3 = 0.924 * Pr**(1/3) * (Re * (d_i / l))**0.5
    return (Nu_mq1**3 + 0.6**3 + (Nu_mq2 - 0.6)**3 + Nu_mq3**3)**(1/3)


def calc_Re_coolant(w: float, d_i: float, rho_coolant: float, cool_eta: float) -> float:
    """Coolant-side (tube) Reynolds number."""
    return (d_i * w * rho_coolant) / cool_eta


def calc_Nu_coolant(w: float, d_i: float, l: float, rho_coolant: float, cool_eta: float, cool_Pr: float):
    """Returns (Nu, Re) for the coolant (tube) side, handling the
    laminar / transitional / turbulent blend."""
    Re = calc_Re_coolant(w, d_i, rho_coolant, cool_eta)

    if Re < 2300:
        Nu = calc_Nu_laminar(Re, cool_Pr, d_i, l)
    elif 2300 <= Re < 4000:
        gamma = (Re - 2300) / (4000 - 2300)
        Nu = (1 - gamma) * calc_Nu_laminar(2300, cool_Pr, d_i, l) + gamma * calc_Nu_turbulent(4000, cool_Pr, d_i, l)
    else:
        Nu = calc_Nu_turbulent(Re, cool_Pr, d_i, l)

    return Nu, Re


def calc_alpha_i(w: float, d_i: float, l: float, rho_coolant: float, cool_eta: float, cool_Pr: float, cool_lambda: float) -> float:
    """Coolant-side heat transfer coefficient [W/m²K]."""
    Nu, Re = calc_Nu_coolant(w, d_i, l, rho_coolant, cool_eta, cool_Pr)
    return (Nu * cool_lambda) / d_i


# =============================================================================
# 3. ORCHESTRATORS -- the main API surface used by the solvers
# =============================================================================

def calc_overall_k(geo, ops, coolant_state, air_state, T_air_out: float) -> float:
    """Overall heat transfer coefficient k [W/m²K], referred to the outer
    (air-side) surface. Called once per solver iteration with the
    current guess for T_air_out."""

    # --- Air side (outer) ---------------------------------------------
    T_air_mean = (ops.T_air_in + T_air_out) / 2.0

    w_eT = calc_w_eT(
        w_o=ops.w_o,
        Ao_Ae_ratio=geo.Ao_Ae_ratio,
        T_mean=T_air_mean,
        T_in=ops.T_air_in,
    )

    Nu_air = calc_Nu_air(
        A_ratio=(geo.A / geo.A_Go),
        Pr_air=air_state.Pr,
        d=geo.d,
        n=geo.n_rows,
        w_eT=w_eT,
        rho_air=air_state.rho,
        eta_air=air_state.eta,
    )

    alpha_R = calc_alpha_R(Nu_air, air_state.lambda_, geo.d)
    eta_R = calc_fin_efficiency_staggered(geo.t_q, geo.t_l, geo.d, alpha_R, geo.lambda_R, geo.s)
    alpha_S = calc_alpha_S(alpha_R, eta_R, geo.A, geo.A_R)

    # --- Coolant side (inner) ------------------------------------------
    alpha_i = calc_alpha_i(
        w=ops.w_coolant,
        d_i=geo.d_i,
        l=geo.l,
        rho_coolant=coolant_state.rho,
        cool_eta=coolant_state.eta,
        cool_Pr=coolant_state.Pr,
        cool_lambda=coolant_state.lambda_,
    )

    # --- Combine into overall k -----------------------------------------
    k_inv = (1 / alpha_S) + (geo.A / geo.A_i) * ((1 / alpha_i) + (geo.d - geo.d_i) / (2 * geo.lambda_R))

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
    T_air_mean = (ops.T_air_in + T_air_out) / 2.0
    w_eT = calc_w_eT(
        w_o=ops.w_o,
        Ao_Ae_ratio=geo.Ao_Ae_ratio,
        T_mean=T_air_mean,
        T_in=ops.T_air_in,
    )
    Re_air = calc_Re_air(geo.d, w_eT, air_state.rho, air_state.eta)
    Nu_air = calc_Nu_air(
        A_ratio=(geo.A / geo.A_Go),
        Pr_air=air_state.Pr,
        d=geo.d,
        n=geo.n_rows,
        w_eT=w_eT,
        rho_air=air_state.rho,
        eta_air=air_state.eta,
    )
    alpha_R = calc_alpha_R(Nu_air, air_state.lambda_, geo.d)

    # --- Coolant side ------------------------------------------------------
    Nu_coolant, Re_coolant = calc_Nu_coolant(
        w=ops.w_coolant,
        d_i=geo.d_i,
        l=geo.l,
        rho_coolant=coolant_state.rho,
        cool_eta=coolant_state.eta,
        cool_Pr=coolant_state.Pr,
    )
    alpha_i = (Nu_coolant * coolant_state.lambda_) / geo.d_i

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


# =============================================================================
# TRASH PILE -- unused code, kept for reference / possible reuse
# =============================================================================
# Everything below is not called anywhere in this codebase (checked July
# 2026). Each block says why it might still be worth keeping.
# =============================================================================

# -----------------------------------------------------------------------
# calc_heat_cap_flow_coolant / calc_heat_cap_flow_air / calc_cap_flow_ratio:
# a dead chain -- calc_cap_flow_ratio() is the only caller of the other
# two, and nothing calls calc_cap_flow_ratio() itself. Looks like an
# early attempt at a capacity-rate ratio (R = W_min/W_max) before the
# NTU solver ended up computing R1 = W1/W2 inline instead (see
# solvers.py, solve_it_NTU). Note these use `geo.d_i` for BOTH sides,
# which is a coolant-side dimension -- almost certainly a bug if this
# were ever revived for the air side.
#
# Probably safe to delete, but flagging rather than removing outright
# since "capacity flow ratio" as a named, reusable helper could still be
# handy if the NTU solver's inline R1/R2 computation ever gets factored
# out.
# -----------------------------------------------------------------------
# def calc_heat_cap_flow_coolant(geo, ops, coolant_state):
#     return coolant_state.rho * geo.d_i * ops.w_coolant
#
# def calc_heat_cap_flow_air(geo, ops, air_state):
#     return air_state.rho * geo.d_i * ops.w_o
#
# def calc_cap_flow_ratio(geo, ops, air_state, coolant_state):
#     return calc_heat_cap_flow_air(geo, ops, air_state) / calc_heat_cap_flow_coolant(geo, ops, coolant_state)