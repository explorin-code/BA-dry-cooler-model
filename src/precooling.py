import dataclasses

import CoolProp.CoolProp as CP
from scipy.optimize import brentq

import src.parameters as parameters


def calc_eta_B ():
    """Return degree of humidification depending on geometry and operational conditions.
    Placeholder: VDI Wärmeatlas M8 §2.2 confirms there is no universal
    formula for eta_B -- it's inherently apparatus-specific empirical data
    (the book discusses several competing empirical forms, e.g. its
    Gl. (16): eta_B = 1 - exp(-K_Y*M_dot_W/M_dot_L), and states none
    generalizes across hardware). If real pad test data becomes available,
    Gl. (16)'s form is the simplest one to fit against it."""
    return 0.8 ## fixed value for now, usually empirical and velocity dependent

def precooling_decider(cell_result):
    """Precool if the already-solved ambient scenario's Cell T_coolant_out
    is above the target. Takes the Cell SolverResult directly -- no
    re-solving, so this can't drift from what was actually reported."""
    target = parameters.T_COOLANT_TARGET_OUT
    decision = cell_result.T_coolant_out > target
    comparison = ">" if decision else "<="
    print(f"[Decider] Cell T_coolant_out = {cell_result.T_coolant_out:.2f} °C {comparison} "
          f"target = {target:.2f} °C -> "
          f"{'precooling engaged' if decision else 'no precooling needed'}")
    return decision

def calc_X(T_wb_guess, P_air, phi):
    """return the humidity ratio of the air at the wet bulb temperature and given pressure"""
    return CP.HAPropsSI('W', 'T', T_wb_guess + 273.15, 'P', P_air, 'R', phi)

def calc_h(T, P, phi):
    """return the air enthalpy. given conditions"""
    return CP.HAPropsSI('Hha', 'T', T + 273.15, 'P', P, 'R', phi)

# Source: VDI Wärmeatlas, Chapter M8, §2.2, Gl. (15), p.1778 -- confirmed
# 2026-08 against the 12th ed., exact match (residual below is literally
# (h_LE-h_LK)/(X_E-X_K) - c_water*theta_K = 0, i.e. Gl.(15) root-found via
# brentq). The formula assumes a Lewis factor of 1 (book's Gl. 14,
# alpha/(beta_Y*c_p(1+Y)) = 1) -- not separately checked here, inherited
# from the book's own stated prerequisite for Gl.(15)/(13) to hold.
def calc_cooling_limit(ops):
    """Return the cooling limit of the dry cooler depending on geometry and operational conditions."""
    if ops.phi >= 0.999:
        # already (essentially) saturated -- no evaporative cooling potential,
        # and the wet-bulb energy balance is singular right at 100% RH
        return ops.T_air_in

    X_E = ops.X
    c_water = 4186.0 # specific heat of water [J/kg-K] from VDI
    h_LE = calc_h(ops.T_air_in, ops.P_air, ops.phi)

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

def calc_precooler(ops, ambient_result):
    """returns (ops, was_precooled). ops unchanged and was_precooled=False
    if the decider says no precooling is needed. ambient_result is the
    already-solved ScenarioResult for `ops` -- the decider reads its Cell
    result rather than re-solving."""
    if precooling_decider(ambient_result.cell):
        eta_B = calc_eta_B()

        theta_K = calc_cooling_limit(ops)
        theta_LE = ops.T_air_in

        X_E = ops.X
        X_K = calc_X(theta_K, ops.P_air, 1.0)

        # Source: VDI Wärmeatlas, Chapter M8, §2.2, Gl. (13), p.1778 --
        # confirmed 2026-08 against the 12th ed. Gl.(13) states
        # eta_B = (X_A-X_E)/(X_K-X_E) ~= (theta_LE-theta_LA)/(theta_LE-theta_K);
        # both lines below are that same equation solved for X_A and
        # theta_LA respectively.
        theta_LA = theta_LE - eta_B * (theta_LE - theta_K)
        X_A = X_E + eta_B * (X_K - X_E)

        ops = dataclasses.replace(
            ops,
            T_air_in=theta_LA,
            X=X_A,
            phi=None,
            w_f=ops.w_f,
            m_dot_2=0.0,
            V_o=0.0,
            u_i=0.0,
            V_coolant=0.0,
        )
        return ops, True

    return ops, False
