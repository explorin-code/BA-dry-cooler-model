"""
pressure_drop.py
=================
Coolant-side (tube) and air-side (finned bundle) pressure drop. Decoupled
from the thermal solve -- only used for fan/pump sizing (see economics.py).
Coolant source: VDI Heat Atlas, Section L1.2, p. 1355-1356 (straight pipe);
Section L1.3, p. 1369 (bends). Air source: Wang, Chi & Chang, Part II
(plate-fin friction factor) + Kays & London core relation.
"""

import warnings
from math import log

from src.heat_transfer_core import calc_Re_coolant, calc_Re_air

# =============================================================================
# Coolant side
# =============================================================================
# Friction coefficients f_c: VDI L1.2, p. 1356, Gl. (4)/(5).
# Pressure-drop equation, bend allowance, viscosity correction: all Towler2022,
# Ch. 19, p. 851, Eq. (19.18)-(19.20).
# Laminar/turbulent blend across Re 2320-3000: author's choice.

def calc_f_c_lam(Re: float) -> float:
    """Laminar f_c (Darcy lambda). Source: VDI L1.2, p. 1356, Gl. (4)."""
    return 64 / Re


def calc_f_c_turb(Re: float) -> float:
    """Turbulent f_c (Darcy lambda), Blasius, valid 3e3<=Re<=1e5.
    Source: VDI L1.2, p. 1356, Gl. (5)."""
    return 0.3164 / (Re**0.25)


def calc_f_coolant(Re: float) -> float:
    """f_c for straight tube flow. Laminar/turbulent branches: VDI L1.2,
    p. 1356, Gl. (4)/(5). Blend across 2320-3000: author's choice."""
    if Re < 2320:
        return calc_f_c_lam(Re)
    elif Re < 3000:
        gamma = (Re - 2320) / (3000 - 2320)
        return (1 - gamma) * calc_f_c_lam(Re) + gamma * calc_f_c_turb(Re)
    else:
        return calc_f_c_turb(Re)


def calc_visc_corr_c(mu_bulk: float, mu_wall: float, Re: float) -> float:
    """Viscosity correction (mu/mu_w)^-m for non-isothermal flow.
    Source: Towler2022, Ch. 19, p. 851, Eq. (19.19).
    NOTE: this m-threshold (2100) is not the same as calc_f_coolant's own
    laminar/turbulent blend band (2320-3000) -- different sources, not
    reconciled. Currently inactive (mu_wall = mu_bulk, see caller)."""
    m = 0.25 if Re < 2100 else 0.14
    return (mu_bulk / mu_wall)**(-m)


def calc_delta_p_coolant(geo, ops, coolant_state) -> float:
    """Total coolant-side pressure drop [Pa]: N_p = N_r passes, each with
    straight-pipe friction over ONE pass length L_p (= geo.L_p) plus
    1.5*(2*N_p - 1)/N_p velocity heads per pass for contraction (0.5) and
    expansion (1.0) at every pass and the N_p - 1 180-degree bends (1.5 each).
    Source: Towler, Chemical Engineering Design, 3rd ed., Sec. 19.8
    (Eq. 19.20); friction factor as Darcy f = 8 j_f (VDI L1.2).
    Fixed 2026-10: previously used geo.l (all passes) as the per-pass
    length, counting friction N_r times too often (10.7 vs. 2.5 kPa)."""
    w_c = ops.w_c
    rho = coolant_state.rho
    mu = coolant_state.mu
    mu_w = mu  # TEMP: no wall-temperature model yet -- visc_corr_c = 1.0 until then

    Re = calc_Re_coolant(w_c, geo.d_i, rho, mu)
    f_c = calc_f_coolant(Re)
    visc_corr_c = calc_visc_corr_c(mu, mu_w, Re)

    N_p = geo.N_r
    L_p = geo.L_p                                   # length of ONE tube pass
    K_pass = 1.5 * (2 * N_p - 1) / N_p                 # = 0.5 + 1.0 + 1.5*(N_p - 1)/N_p
    return N_p * (f_c * (L_p / geo.d_i) * visc_corr_c + K_pass) * (rho * w_c**2 / 2)
    

# =============================================================================
# Air side
# =============================================================================
# Plate-fin friction factor: Wang, Chi & Chang, Part II (2000), Eq. (12)-(15).
# Pressure-drop relation: Kays & London core equation as used by Wang & Chi,
# Part I, Eq. (17). Wang's f was reduced from data with this same relation
# and no separate entrance/exit coefficients (Kc/Ke) -- do not add them here.
# Minimum-gap velocity/area reuse Geometry.Afr_Ae_ratio (= Afr/Ac = 1/sigma).
# =============================================================================

# Wang Part II applicability (p. 2700, Conclusions), SI units.
WANG_RANGE = {
    'N_r': (1, 6),
    'd':    (6.35e-3, 12.7e-3),
    'F_p':    (1.19e-3, 8.7e-3),
    'P_t':    (17.7e-3, 31.75e-3),
    'P_l':    (12.4e-3, 27.5e-3),
}


def check_wang_range(geo) -> None:
    """Warns for every geometry value outside Wang's tested range."""
    for name, (lo, hi) in WANG_RANGE.items():
        value = getattr(geo, name)
        if not lo <= value <= hi:
            warnings.warn(f"Wang correlation: {name} = {value:g} outside tested range [{lo:g}, {hi:g}]")


def calc_wang_F1(P_t: float, P_l: float, F_p: float, D_c: float, N_r: int) -> float:
    """Source: Wang et al., Part II, p. 2699, Eq. (13):
    F1 = -0.764 + 0.739*(Pt/Pl) + 0.177*(Fp/Dc) - 0.00758/N."""
    return -0.764 + 0.739 * (P_t / P_l) + 0.177 * (F_p / D_c) - 0.00758 / N_r


def calc_wang_F2(Re_Dc: float) -> float:
    """Source: Wang et al., Part II, p. 2699, Eq. (14): F2 = -15.689 + 64.021/ln(ReDc)."""
    return -15.689 + 64.021 / log(Re_Dc)


def calc_wang_F3(Re_Dc: float) -> float:
    """Source: Wang et al., Part II, p. 2699, Eq. (15): F3 = 1.696 - 15.695/ln(ReDc)."""
    return 1.696 - 15.695 / log(Re_Dc)


def calc_f_wang(Re_Dc: float, P_t: float, P_l: float, F_p: float, D_c: float, N_r: int) -> float:
    """Wang plain-fin friction factor f [-]. Tube pitches P_t/P_l, fin
    pitch F_p and collar diameter D_c follow Wang's own notation. Re_Dc is based on D_c and the minimum-gap
    velocity w_e (heat_transfer_core.calc_Re_air with d = D_c).
    Source: Wang et al., Part II, p. 2699, Eq. (12):
    f = 0.0267 * ReDc^F1 * (Pt/Pl)^F2 * (Fp/Dc)^F3.
    Applicability: see WANG_RANGE."""
    F1 = calc_wang_F1(P_t, P_l, F_p, D_c, N_r)
    F2 = calc_wang_F2(Re_Dc)
    F3 = calc_wang_F3(Re_Dc)
    return 0.0267 * Re_Dc**F1 * (P_t / P_l)**F2 * (F_p / D_c)**F3


def calc_G_e(m_dot_air: float, A_e: float) -> float:
    """Mass flux at the minimum flow area [kg/m^2s].
    Source: Wang & Chi, Part I, Nomenclature, p. 2682: Gc = rho*Vmax."""
    return m_dot_air / A_e


def calc_delta_p_platefin_air(f: float, A_tot: float, A_e: float, G_e: float, sigma: float,
                              rho_a_i: float, rho_a_o: float = None) -> float:
    """Air-side pressure drop across a continuous plate-fin bundle [Pa].
    rho_a_o defaults to rho_a_i (no acceleration term) if not given.
    Source: Wang & Chi, Part I, p. 2687, Eq. (17), citing Kays and
    London, Compact Heat Exchangers, 3rd ed., 1984."""
    if rho_a_o is None:
        rho_a_o = rho_a_i
    rho_a_m = 2.0 * rho_a_i * rho_a_o / (rho_a_i + rho_a_o)  # mean of specific volumes
    return (G_e**2 / (2.0 * rho_a_i)) * (f * (A_tot / A_e) * (rho_a_i / rho_a_m)
                                        + (1.0 + sigma**2) * (rho_a_i / rho_a_o - 1.0))


def calc_delta_p_bundle_air(geo, ops, air_in, air_out, air_mean) -> float:
    """Air-side pressure drop across the finned-tube bundle [Pa], Wang
    plain-fin correlation. air_in/air_out: states at the bundle inlet/
    outlet (Kays & London density terms). air_mean: state at the mean air
    temperature, used for Re_Dc (and thus f)."""
    check_wang_range(geo)

    A_e = geo.A_e  # TODO check other ratio from WangII source
    G_e = calc_G_e(ops.m_dot_a, A_e)
    sigma = 1.0 / geo.Afr_Ae_ratio

    # G_e is fixed by continuity, so the min-gap velocity at the mean
    # state is G_e / rho_mean.
    w_e_m = G_e / air_mean.rho
    Re_Dc = calc_Re_air(geo.D_c, w_e_m, air_mean.rho, air_mean.mu)
    f = calc_f_wang(Re_Dc, geo.P_t, geo.P_l, geo.F_p, geo.D_c, geo.N_r)

    return calc_delta_p_platefin_air(f, geo.A_tot, A_e, G_e, sigma, air_in.rho, air_out.rho)


def calc_delta_p_pad() -> float:
    """Pressure drop across the adiabatic pre-cooler pad [Pa]."""
    return 50  # 300 mm rigid-media pad: 34.8-74.6 Pa (ASHRAE Handbook 2020, Ch. 41); dry operation: see PAD_IN_DRY_AIR_PATH


def calc_delta_p_air_total(Delta_P_bundle: float, Delta_P_pad: float = None) -> float:
    """Total air-side pressure drop [Pa]: bundle plus pad, if the air passes
    the pad (always when precooling; in dry operation see
    parameters.PAD_IN_DRY_AIR_PATH). Delta_P_pad None = no pad in the path."""
    if Delta_P_pad is None:
        return Delta_P_bundle
    return Delta_P_bundle + Delta_P_pad