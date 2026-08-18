import dataclasses

import src.param_loader as param_loader
import src.solvers as solvers
import src.fluid_properties as fluid_props
import CoolProp.CoolProp as CP
from scipy.optimize import brentq

TARGET_T_COOLANT_OUT = 25.0  # °C -- Konrad's target output temp

def calc_eta_B ():
    """Return degree of humidification depending on geometry and operational conditions."""
    return 0.8 ## fixed value for now, usually empirical and velocity dependent

def precooling_decider(ops):
    """Run the cell method once against ambient ops; precool if the
    resulting T_coolant_out is above the target."""
    _, _, T_coolant_out, *_ = solvers.solve_it_cell(ops=ops)
    return T_coolant_out > TARGET_T_COOLANT_OUT

def calc_X(T_wb_guess, P_air, phi_air):
    """return the humidity ratio of the air at the wet bulb temperature and given pressure"""
    return CP.HAPropsSI('W', 'T', T_wb_guess + 273.15, 'P', P_air, 'R', phi_air)

def calc_h(T, P, phi_air):
    """return the air enthalpy. given conditions"""
    return CP.HAPropsSI('Hha', 'T', T + 273.15, 'P', P, 'R', phi_air)

# AI-REVIEW: VDI M8 2.2 wet-bulb energy balance, root-found via brentq --
# not yet checked against a known wet-bulb reference point (e.g. 20degC/30%
# RH ambient). See CLAUDE.md.
def calc_cooling_limit(ops):
    """Return the cooling limit of the dry cooler depending on geometry and operational conditions."""
    if ops.phi_air >= 0.999:
        # already (essentially) saturated -- no evaporative cooling potential,
        # and the wet-bulb energy balance is singular right at 100% RH
        return ops.T_air_in

    X_E = ops.X_air
    c_water = 4186.0 # specific heat of water [J/kg-K] from VDI
    h_LE = calc_h(ops.T_air_in, ops.P_air, ops.phi_air)

    def residual(theta_K_guess):
        h_LK = calc_h(theta_K_guess, ops.P_air, 1.0)
        X_K = calc_X(theta_K_guess, ops.P_air, 1.0)

        return (h_LE - h_LK) / (X_E - X_K) - c_water * theta_K_guess

    T_dewpoint = CP.HAPropsSI('D', 'T', ops.T_air_in + 273.15, 'P', ops.P_air, 'W', X_E) - 273.15

    try:
        return brentq(residual, T_dewpoint + 0.1, ops.T_air_in)
    except ValueError:
        # bracket collapsed (near-saturated air) or failed to bracket a root
        # -- fall back to "no cooling potential" rather than crashing
        return ops.T_air_in

def calc_precooler(ops):
    """returns (ops, was_precooled). ops unchanged and was_precooled=False
    if the decider says no precooling is needed."""
    if precooling_decider(ops):
        eta_B = calc_eta_B()

        theta_K = calc_cooling_limit(ops)
        theta_LE = ops.T_air_in

        X_E = ops.X_air
        X_K = calc_X(theta_K, ops.P_air, 1.0)

        # AI-REVIEW: temperature/humidity interpolation formulas below not
        # yet independently verified. See CLAUDE.md.
        theta_LA = theta_LE - eta_B * (theta_LE - theta_K)
        X_A = X_E + eta_B * (X_K - X_E)

        ops = dataclasses.replace(
            ops,
            T_air_in=theta_LA,
            X_air=X_A,
            phi_air=None,
            w_o=ops.w_o,
            m_o=0.0,
            V_o=0.0,
            w_coolant=0.0,
            V_coolant=0.0,
        )
        return ops, True

    return ops, False