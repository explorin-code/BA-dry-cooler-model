from parameters import OperatingConditions
from dry_cooler_physics import get_geometry
from fluid_properties import get_fluid_properties
from heat_transfer_core import calc_overall_k
from math import log, pi

def solve_it_LMTD():
    # 1. Load static configurations
    geo = get_geometry()
    ops = OperatingConditions() 

    # 2. Calculate static mass flows (evaluated at the inlet!)
    rho_coolant_in = get_fluid_properties(ops.coolant_type, ops.T_coolant_in, ops.P_coolant).rho
    rho_air_in = get_fluid_properties(ops.air_type, ops.T_air_in, ops.P_air).rho
    
    # Assuming a 1-pass tube arrangement for the inner flow area
    A_flow_coolant = geo.n_tubes * (pi / 4) * (geo.d_i ** 2)

    dm_coolant = rho_coolant_in * ops.w_coolant * A_flow_coolant
    dm_air = rho_air_in * ops.w_o * geo.inflow_cross_section

    # 3. Initial guesses for the iterative loop
    dT_hot_it = 15.0
    dT_cold_it = 15.0
    dT_hot = 5.0    # Forced difference to kickstart the while loop
    dT_cold = 5.0
    
    T_coolant_out = ops.T_coolant_in - dT_hot_it    
    T_air_out = ops.T_air_in + dT_cold_it

    # Convergence threshold (e.g., 0.001 K)
    threshold = 1e-3 

    # Initialize these with values larger than the threshold
    diff_hot = 1.0 
    diff_cold = 1.0

    # Store history for return
    history_hot = []
    history_cold = []

    # while either is still moving significantly
    while (abs(diff_hot) > threshold) or ((diff_cold) > threshold):
        
        dT_hot = dT_hot_it
        dT_cold = dT_cold_it
        
        # Calculate dynamic mean temperatures for this loop
        T_air_mean = (T_air_out + ops.T_air_in) / 2.0
        T_coolant_mean = (ops.T_coolant_in + T_coolant_out) / 2.0

        # Get dynamic fluid states at the mean temperature
        coolant_state = get_fluid_properties(ops.coolant_type, T_coolant_mean, ops.P_coolant)
        air_state = get_fluid_properties(ops.air_type, T_air_mean, ops.P_air)

        # Pass everything to the orchestrator
        k = calc_overall_k(
            geo=geo, 
            ops=ops, 
            coolant_state=coolant_state, 
            air_state=air_state, 
            T_air_out=T_air_out
        )

        # Log Mean Temperature Difference
        LMTD = calc_LMTD(dT_hot, dT_cold)
        
        # Total Area (Make sure n_rows is defined in your Geometry class!)
        A_tot = geo.A * geo.n_tubes * geo.n_rows

        # Total Heat Transfer
        dQ = A_tot * k * LMTD

        # Calculate new temperature drops based on the heat transferred
        dT_air_rise = dQ / (dm_air * air_state.cp)               
        dT_coolant_drop = dQ / (dm_coolant * coolant_state.cp)   

        # calculate the raw exiting temperatures for the next iteration
        raw_T_coolant_out = ops.T_coolant_in - dT_coolant_drop
        raw_T_air_out = ops.T_air_in + dT_air_rise

        # UNDER-RELAXATION: Slow down the update!
        # omega = 0.1 means "take 10% of the new guess, keep 90% of the old guess"
        omega = 0.1 
        T_coolant_out = (1 - omega) * T_coolant_out + omega * raw_T_coolant_out
        T_air_out = (1 - omega) * T_air_out + omega * raw_T_air_out

        # Calculate new diffs for the loop condition
        new_dT_hot = ops.T_coolant_in - T_air_out
        new_dT_cold = T_coolant_out - ops.T_air_in
        
        diff_hot = new_dT_hot - dT_hot_it
        diff_cold = new_dT_cold - dT_cold_it
        
        dT_hot_it = new_dT_hot
        dT_cold_it = new_dT_cold

        # Add this to see what the solver is actually doing!
        print(f"Hot Error: {diff_hot:.5f} | Cold Error: {diff_cold:.5f} | T_coolant_out: {T_coolant_out:.2f} | T_air_out: {T_air_out:.2f}")
        
        history_hot.append(diff_hot)
        history_cold.append(diff_cold)

    return k, dQ, T_coolant_out, T_air_out, history_hot, history_cold

def calc_LMTD(dT_hot: float, dT_cold: float) -> float:
    """Calculates LMTD and safely catches mathematical singularities."""
    # Prevent math domain error from temperature crossovers
    if dT_hot <= 0 or dT_cold <= 0:
        return 1e-5 
        
    # Prevent division by zero if differences are identical
    if abs(dT_hot - dT_cold) < 1e-5:
        return dT_hot
        
    return (dT_hot - dT_cold) / log(dT_hot / dT_cold)