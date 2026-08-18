"""
pressure_drop.py
=================
Coolant-side (tube) and air-side (finned bundle) pressure drop. Decoupled
from the thermal solve -- only used for fan/pump sizing (see economics.py).
Coolant source: VDI Heat Atlas, Section L1.2, p. 1355-1356 (straight pipe);
Section L1.3, p. 1369 (bends). Air source: VDI Heat Atlas, Section L1.4
(single bundle), Section L1.5 (array-specific).
"""

from math import log10, pi, exp, sqrt

from src.heat_transfer_core import calc_Re_coolant

ZETA_U_DEFAULT = 1.0  # bend resistance coefficient -- VDI L1.3 gives 0.5-1.5 for 180deg bends


# =============================================================================
# Coolant side
# =============================================================================

# AI-REVIEW: the 2320-10^4 linear transition blend below is a
# Claude-suggested approximation, not something the VDI source itself
# specifies. See CLAUDE.md.
def calc_zeta_coolant(Re: float) -> float:
    """Friction factor (Widerstandsbeiwert) for straight tube flow.
    Source: VDI Section L1.2, p. 1356, Gl. (4) [laminar] & Gl. (6)
    [turbulent, Hanakov/Konakov]. Linearly blended across the gap between
    the two correlations' stated validity ranges (2320-10^4), same
    approach as calc_Nu_coolant's laminar/turbulent blend."""
    if Re < 2320:
        return 64 / Re
    elif Re < 1e4:
        gamma = (Re - 2320) / (1e4 - 2320)
        zeta_lam = 64 / 2320
        zeta_turb = (1.8 * log10(1e4) - 1.5)**(-2)
        return (1 - gamma) * zeta_lam + gamma * zeta_turb
    else:
        return (1.8 * log10(Re) - 1.5)**(-2)


def calc_delta_p_straight(zeta: float, l: float, d_i: float, rho: float, u_i: float) -> float:
    """Straight-pipe pressure drop [Pa]. Source: VDI L1.2, p. 1355, Gl. (1)."""
    return zeta * (l / d_i) * (rho * u_i**2 / 2)


def calc_delta_p_bend(zeta_u: float, rho: float, u_i: float) -> float:
    """Single 180° bend pressure drop [Pa]. Source: VDI L1.3, p. 1369, Abb. 16."""
    return zeta_u * (rho * u_i**2 / 2)


def calc_delta_p_total(geo, ops, coolant_state, zeta_u: float = ZETA_U_DEFAULT) -> float:
    """Total coolant-side pressure drop [Pa]: straight pipe plus one 180°
    bend per row-to-row transition in the serpentine."""
    u_i = ops.w_coolant
    rho = coolant_state.rho
    eta = coolant_state.eta

    Re = calc_Re_coolant(u_i, geo.d_i, rho, eta)
    zeta = calc_zeta_coolant(Re)

    delta_p_straight = calc_delta_p_straight(zeta, geo.l, geo.d_i, rho, u_i)

    n_bends = geo.n_rows - 1
    delta_p_bend = calc_delta_p_bend(zeta_u, rho, u_i)

    return delta_p_straight + n_bends * delta_p_bend


# =============================================================================
# Air side
# =============================================================================
# AI-REVIEW: this whole section was transcribed from a secondary/paraphrased
# description of VDI L1.4, not the original Waermeatlas text -- cross-check
# every equation below against the actual source before trusting results.
# See CLAUDE.md.
# =============================================================================

def calc_air_pitch_ratios(geo):
    """Transverse/longitudinal/diagonal pitch ratios [-].
    Source: VDI L1.4, p. 1376, Gl. (6), (7), (8)."""
    a = geo.t_q / geo.d
    b = geo.t_l / geo.d
    c = sqrt((a / 2)**2 + b**2)
    return a, b, c


def calc_narrow_gap_velocity_staggered(w_o: float, a: float, b: float, c: float, n_rows: int):
    """Narrowest-gap velocity we [m/s], row-resistance count nW [-], and
    whether the diagonal (vs. transverse) gap is narrowest, for a
    staggered bundle. Source: VDI L1.4, p. 1376, Gl. (2)-(5)."""
    diagonal_narrowest = b < 0.5 * sqrt(2 * a + 1)
    if diagonal_narrowest:
        we = w_o * (a / (2 * (c - 1)))
        nW = n_rows - 1
    else:
        we = w_o * (a / (a - 1))
        nW = n_rows
    return we, nW, diagonal_narrowest


def calc_narrow_gap_velocity_inline(w_o: float, a: float, n_rows: int):
    """STATUS: INACTIVE -- not currently wired into calc_delta_p_bundle_air(),
    which uses the staggered functions instead. Kept here, ready to swap in
    (same shape) if the tube layout ever changes from staggered to in-line.
    Source: VDI L1.4, p. 1376, Gl. (2)-(3)."""
    we = w_o * (a / (a - 1))
    nW = n_rows
    return we, nW


def calc_Re_air_bundle(we: float, d: float, nu_air: float) -> float:
    """Reynolds number for the bundle pressure-drop correlation, based on
    the narrowest-gap velocity we (NOT the same Re as the heat-transfer
    correlations, which use the thermally-corrected w_eT -- see
    heat_transfer_core.calc_Re_air for that one).
    Source: VDI L1.4, p. 1380, Gl. (20)."""
    return (we * d) / nu_air


def calc_zeta_air_staggered(Re: float, a: float, b: float, c: float, diagonal_narrowest: bool) -> float:
    """Ideal friction coefficient for a staggered bundle.
    Source: VDI L1.4, p. 1378-1379, Gl. (16)-(22)."""
    Fv = 1 - exp(-(Re + 200) / 1000)

    if diagonal_narrowest:
        fa_l = (280 * pi * (b**0.5 - 0.6)**2 + 0.75) / ((4 * a * b - pi) * c**1.6)
    else:
        fa_l = (280 * pi * (b**0.5 - 0.6)**2 + 0.75) / ((4 * a * b - pi) * a**1.6)
    xi_l = fa_l / Re

    fa_t = 2.5 + 1.2 / (a - 0.85)**1.08 + 0.4 * (b/a - 1)**3 - 0.01 * (a/b - 1)**3
    xi_t = fa_t / (Re**0.25)

    return xi_l + xi_t * Fv


def calc_zeta_air_inline(Re: float, a: float, b: float) -> float:
    """STATUS: INACTIVE -- see calc_narrow_gap_velocity_inline() above.
    Ideal friction coefficient for an in-line bundle.
    Source: VDI L1.4, p. 1377-1378, Gl. (10)-(15)."""
    Ff = 1 - exp(-(Re + 1000) / 2000)

    fa_l = (280 * pi * (b**0.5 - 0.6)**2 + 0.75) / ((4 * a * b - pi) * a**1.6)
    xi_l = fa_l / Re

    fa_t = 0.22 + 1.2 * (1 - 0.94/b)**0.6 / (a - 0.85)**1.3 + 10**(0.47 * (b/a) - 1.5) + 0.03 * (a - 1) * (b - 1)
    xi_t = fa_t / (Re**(0.1 * b / a))

    return xi_l + xi_t * Ff


def calc_viscosity_correction(eta_wall: float, eta_bulk: float) -> float:
    """Near-wall viscosity correction for non-isothermal flow.
    Source: VDI L1.4, p. 1378, Gl. (24)."""
    return (eta_wall / eta_bulk)**0.14


def calc_delta_p_bundle_air(geo, ops, air_state, eta_wall: float = None) -> float:
    """Air-side pressure drop across the finned-tube bundle [Pa], staggered
    arrangement. eta_wall defaults to the bulk viscosity (fz_t = 1.0, no
    correction) until wall temperature is modeled elsewhere.
    Source: VDI L1.4/L1.5 (see individual functions for exact citations)."""
    if eta_wall is None:
        eta_wall = air_state.eta  # TEMP: no wall-temperature model yet -- fz_t = 1.0 until then

    a, b, c = calc_air_pitch_ratios(geo)
    we, nW, diagonal_narrowest = calc_narrow_gap_velocity_staggered(ops.w_o, a, b, c, geo.n_rows)

    nu_air = air_state.eta / air_state.rho
    Re = calc_Re_air_bundle(we, geo.d, nu_air)

    xi_ideal = calc_zeta_air_staggered(Re, a, b, c, diagonal_narrowest)
    fz_t = calc_viscosity_correction(eta_wall, air_state.eta)

    return xi_ideal * fz_t * nW * (air_state.rho * we**2 / 2)


def calc_delta_p_pad() -> float:
    """Pressure drop across the adiabatic pre-cooler pad [Pa]."""
    return None  ## not yet implemented -- needs a manufacturer correlation (e.g. S. He et al.)


def calc_delta_p_air_total(delta_p_bundle: float, delta_p_pad: float = None) -> float:
    """Total air-side pressure drop [Pa]: bundle plus pad, if the pad's
    contribution is known -- otherwise just the bundle (not None), so
    fan power can still be reported with the caveat that it's a lower
    bound while the pad is active but unmodeled.
    Source: DuEPublico thesis ref, p. 14."""
    if delta_p_pad is None:
        return delta_p_bundle
    return delta_p_bundle + delta_p_pad