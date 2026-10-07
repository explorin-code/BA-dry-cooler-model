import dataclasses
import logging

import CoolProp.CoolProp as CP
from scipy.optimize import brentq

import src.parameters as parameters

logger = logging.getLogger(__name__)


def calc_eta_B ():
    """Return degree of humidification depending on geometry and operational conditions.
    Placeholder: VDI Wärmeatlas M8 §2.2 confirms there is no universal
    formula for eta_B -- it's inherently apparatus-specific empirical data
    (the book discusses several competing empirical forms, e.g. its
    Gl. (16): eta_B = 1 - exp(-K_Y*M_dot_W/M_dot_L), and states none
    generalizes across hardware). If real pad test data becomes available,
    Gl. (16)'s form is the simplest one to fit against it."""
    # Assumption: 0.8, within ASHRAE Handbook 2020, Ch. 41: rigid-media pads
    # reach 75-95 %, random-media pads about 80 %.
    return 0.8

def precooling_decider(cell_result):
    """Precool if the already-solved ambient scenario's Cell T_c_o
    is above the target. Takes the Cell SolverResult directly -- no
    re-solving, so this can't drift from what was actually reported."""
    target = parameters.T_COOLANT_TARGET_OUT
    decision = cell_result.T_c_o > target
    comparison = ">" if decision else "<="
    logger.info("[Decider] Cell T_coolant_out = %.2f °C %s target = %.2f °C -> %s",
                cell_result.T_c_o, comparison, target,
                'precooling engaged' if decision else 'no precooling needed')
    return decision

def calc_c_W(T, P):
    """Specific heat capacity of LIQUID water [J/kg-K] at T [°C] -- the
    c_W of the cooling-limit equation (c_W * T_K = enthalpy of the
    evaporating water). Clamped to >= 0.01 °C: CoolProp's water has no ice.
    Source: CoolProp (Bell et al. 2014), IAPWS-95 water."""
    return CP.PropsSI('C', 'T', max(T, 0.01) + 273.15, 'P', P, 'Water')


def calc_Y(T_wb_guess, p_a, phi):
    """Humidity ratio [kg water/kg dry air] at the given temperature,
    pressure and relative humidity. Source: CoolProp HAPropsSI
    (Bell et al. 2014)."""
    return CP.HAPropsSI('W', 'T', T_wb_guess + 273.15, 'P', p_a, 'R', phi)

def calc_h(T, P, phi):
    """Moist-air enthalpy per kg DRY air [J/kg dry air] (VDI's h_1+X, the
    basis of M8 Gl. (15), whose denominator Y_i - Y_K is also per kg dry
    air). Source: CoolProp HAPropsSI (Bell et al. 2014), key 'Hda'.
    Fixed 2026-10: was 'Hha' (per kg humid air), which put T_K ~0.05 K too
    high (20 °C / 30 %: 10.89 instead of 10.84 °C = CoolProp's own
    thermodynamic wet-bulb temperature)."""
    return CP.HAPropsSI('Hda', 'T', T + 273.15, 'P', P, 'R', phi)

# Source: VDI Wärmeatlas, Chapter M8, §2.2, Gl. (15), p.1778 -- confirmed
# 2026-08 against the 12th ed., exact match (residual below is literally
# (h_LE-h_LK)/(Y_i-Y_K) - c_W*T_K = 0, i.e. Gl.(15) root-found via
# brentq). The formula assumes a Lewis factor of 1 (book's Gl. 14,
# alpha/(beta_Y*c_p(1+Y)) = 1) -- not separately checked here, inherited
# from the book's own stated prerequisite for Gl.(15)/(13) to hold.
def calc_cooling_limit(ops):
    """Return the cooling limit of the dry cooler depending on geometry and operational conditions."""
    if ops.phi >= 0.999:
        # already (essentially) saturated -- no evaporative cooling potential,
        # and the wet-bulb energy balance is singular right at 100% RH
        return ops.T_a_i

    Y_i = ops.Y
    h_LE = calc_h(ops.T_a_i, ops.p_a, ops.phi)

    def residual(T_K_guess):
        h_LK = calc_h(T_K_guess, ops.p_a, 1.0)
        Y_K = calc_Y(T_K_guess, ops.p_a, 1.0)

        return (h_LE - h_LK) / (Y_i - Y_K) - calc_c_W(T_K_guess, ops.p_a) * T_K_guess

    T_dp = CP.HAPropsSI('D', 'T', ops.T_a_i + 273.15, 'P', ops.p_a, 'W', Y_i) - 273.15

    try:
        return brentq(residual, T_dp + 0.1, ops.T_a_i)
    except ValueError:
        # bracket collapsed (near-saturated air) or failed to bracket a root
        # -- fall back to "no cooling potential" rather than crashing
        return ops.T_a_i

def calc_precooler(ops, ambient_result):
    """returns (ops, was_precooled). ops unchanged and was_precooled=False
    if the decider says no precooling is needed. ambient_result is the
    already-solved ScenarioResult for `ops` -- the decider reads its Cell
    result rather than re-solving."""
    if precooling_decider(ambient_result.cell):
        eta_B = calc_eta_B()

        T_K = calc_cooling_limit(ops)
        T_a_i = ops.T_a_i

        Y_i = ops.Y
        Y_K = calc_Y(T_K, ops.p_a, 1.0)

        # Source: VDI Wärmeatlas, Chapter M8, §2.2, Gl. (13), p.1778 --
        # confirmed 2026-08 against the 12th ed. Gl.(13) states
        # eta_B = (Y_pc-Y_i)/(Y_K-Y_i) ~= (T_a_i-T_a_pc)/(T_a_i-T_K);
        # both lines below are that same equation solved for Y_pc and
        # T_a_pc respectively.
        T_a_pc = T_a_i - eta_B * (T_a_i - T_K)
        Y_pc = Y_i + eta_B * (Y_K - Y_i)

        ops = dataclasses.replace(
            ops,
            T_a_i=T_a_pc,
            Y=Y_pc,
            phi=None,
            w_fr=ops.w_fr,
            m_dot_a=0.0,
            V_dot_a=0.0,
            w_c=0.0,
            V_dot_c=0.0,
        )
        return ops, True

    return ops, False
