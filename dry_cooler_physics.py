from math import log10, pi, ln, tanh
from CoolProp.CoolProp import PropsSI

import fluid_properties

# Cooler topologies, flow arrangements, nusselt correlations, fin and tube geometries

def calc_Nu_turbulent(Re: float, Pr: float, d_i: float, L: float):

    psi = (1.8 * log10(Re) - 1.5)**(-2)                                                             # Correction factor for turbulent flow (G1 27)
    Nu_mT = (((psi/8)*(Re - 1000)*Pr)/(1 + 12.7*(psi/8)**(1/2)*(Pr**(2/3) - 1)))*(1+(d_i/L)**(2/3)) # Gnielinski equation for turbulent flow with correction factor (G1 26)

    return Nu_mT

def calc_Nu_laminar(Re: float, Pr: float, d_i: float, L: float):

    Nu_mq1 = 4.364                                                      # Nusselt number for fully developed laminar flow with constant heat flux (G1 17)
    Nu_mq2 = 1.953 * (Re * Pr * (d_i/L))**(1/3)                         # Nusselt number for developing laminar flow with constant heat flux (G1 18)
    Nu_mq3 = 0.924 * Pr**(1/3) * (Re * (d_i/L))**(1/2)                  # Nusselt number for developing laminar flow with constant heat flux and Prandtl number correction (G1 24)

    Nu_mL = (Nu_mq1**3 + 0.6**3 + (Nu_mq2-0.6)**3 + Nu_mq3**3)**(1/3)   # Combined Nusselt number for laminar flow with constant heat flux (G1 25)

    return Nu_mL

def calc_alpha_i(w: float, d_i: float, L: float, coolant: str, T_coolant: float, P_coolant: float):
   
    rho, c_p, lambda_coolant, eta, Pr = get_coolant_properties(P_coolant, T_coolant, coolant)

    Re = (w * d_i) / eta                                              # Reynolds number
    Nu                                                                # Nusselt number

    if Re < 2300:                                                                                           # Laminar flow (G1, 3.1)
        Nu = calc_Nu_laminar(Re, Pr, d_i, L)
       
    elif 2300 <= Re < 4000:                                                                                 # Transitional flow (G1, 4.2)
        gamma = (Re - 2300) / (4000 - 2300)                                                                 # Transition parameter (G1 30)
        Nu = (1-gamma) * calc_Nu_laminar(2300, Pr, d_i, L) + gamma * calc_Nu_turbulent(4000, Pr, d_i, L)    # Linear interpolation between laminar and turbulent Nusselt numbers (G1 29)
        
    else:                                                                                                   # turbulent flow (G1 4.1)
        Nu = calc_Nu_turbulent(Re, Pr, d_i, L)

    alpha_i = (Nu * lambda_coolant) / d_i                                                                     # Convective heat transfer coefficient [W/m²-K]

    return alpha_i


def get_geometry():                     # exemplary values for geometry from M1 4
    D = 0.056       # Outer diameter of the fin [m]
    d = 0.0254      # Inner diameter of the tube [m]
    s = 0.0004      # Fin thickness [m]
    a = 0.00242     # Fin spacing [m]
    d_i = 0.0254    # Inner diameter of the tube [m]
    
    material = 'Aluminum'  # Fin material
    lambda_R = PropsSI('L', 'T', 293.15, 'P', 101325, material) # Thermal conductivity of the fin material at room temperature [W/m-K]

    fin_density = 9*2.54*100   # Number of fins per meter
    t_R = 0.00282   # fin spacing [m] (from fin left edge to left edge of next fin)
    t_q = 0.06      # tube spacing [m] (from tube center to tube center)
    inflow_cross_section = 1 # Cross-sectional area for the inflow [m²]

    w_o = 2.0  # Inflow velocity [m/s]
    
    coolant = 'Water'  # Coolant fluid
    v_coolant = 0.5  # Cooling velocity [m/s]

    n_tubes = 17

    width = n_tubes * t_q                  # Width of the cooler [m]
    height = inflow_cross_section / width  # Height of the cooler [m]

    fins_per_pipe = 348

    A_R = 2 * (pi/4) * (D**2 - d**2) * fins_per_pipe  # Outer surface area of a segment of the tube between two rips [m²]
    A_G = (fins_per_pipe + 1) * pi * d * a
    
    A = A_R * A_G
    A_i = height * d_i * pi                 # Inner surface area of one tube [m²]

    Ao_Ae_ratio = (t_q*(a+s))/((t_q-d)*a + (t_q-D)*s)  # Ratio of the total cross-sectional area to the minimum cross-sectional area for the airflow between the fins and tubes [-]

    return D, d, s, a, d_i, material, lambda_R, fin_density, t_R, t_q, inflow_cross_section, v_coolant, w_o, A_i, Ao_Ae_ratio



def calc_fin_efficiency(D: str, d: str, alpha_R: float, lambda_R: float, s: float):         # fluchtende Anordnung
    phi = (D/d - 1) * (1 + 0.35 * ln(D/d))
    X = phi * d/2 * ((2 * alpha_R) / (lambda_R * s))**(1/2)
    eta_R = tanh(X)/X

    return eta_R
