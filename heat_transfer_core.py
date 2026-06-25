from math import log10, log, tanh

# ==========================================
# 1. PURE MATH WORKERS
# ==========================================

def calc_w_e(w_o: float, Ao_Ae_ratio: float) -> float:
    return w_o * Ao_Ae_ratio

def calc_w_eT(w_o: float, Ao_Ae_ratio: float, T_mean: float, T_in: float) -> float:
    w_e = calc_w_e(w_o, Ao_Ae_ratio)
    return w_e * ((273.15 + T_mean) / (273.15 + T_in))

def calc_fin_efficiency(D: float, d: float, alpha_R: float, lambda_R: float, s: float) -> float:
    phi = (D/d - 1) * (1 + 0.35 * log(D/d))
    X = phi * d/2 * ((2 * alpha_R) / (lambda_R * s))**(0.5)
    return tanh(X)/X

def calc_Nu_air(A_ratio: float, Pr_air: float, d: float, w_eT: float, rho_air: float, eta_air: float) -> float:
    Re_air = (d * w_eT * rho_air) / eta_air
    return 0.22 * (Re_air**0.6) * (A_ratio**(-0.15)) * (Pr_air**(1/3))

def calc_alpha_R(Nu_air: float, lambda_air: float, d: float) -> float: 
    return (Nu_air * lambda_air) / d

def calc_alpha_S(alpha_R: float, eta_R: float, A: float, A_R: float) -> float:
    return alpha_R * (1 - (1 - eta_R) * (A_R / A))

def calc_Nu_turbulent(Re: float, Pr: float, d_i: float, l: float) -> float:
    psi = (1.8 * log10(Re) - 1.5)**(-2)                                                             
    return (((psi / 8) * (Re - 1000) * Pr) / (1 + 12.7 * (psi / 8)**(0.5) * (Pr**(2/3) - 1))) * (1 + (d_i / l)**(2/3)) 

def calc_Nu_laminar(Re: float, Pr: float, d_i: float, l: float) -> float:
    Nu_mq1 = 4.364                                                      
    Nu_mq2 = 1.953 * (Re * Pr * (d_i / l))**(1/3)                         
    Nu_mq3 = 0.924 * Pr**(1/3) * (Re * (d_i / l))**(0.5)                  
    return (Nu_mq1**3 + 0.6**3 + (Nu_mq2 - 0.6)**3 + Nu_mq3**3)**(1/3)   

def calc_alpha_i(w: float, d_i: float, l: float, cool_eta: float, cool_Pr: float, cool_lambda: float) -> float:
    Re = (w * d_i) / cool_eta                                              

    if Re < 2300:                                                                                           
        Nu = calc_Nu_laminar(Re, cool_Pr, d_i, l)
    elif 2300 <= Re < 4000:                                                                                 
        gamma = (Re - 2300) / (4000 - 2300)                                                                 
        Nu = (1 - gamma) * calc_Nu_laminar(2300, cool_Pr, d_i, l) + gamma * calc_Nu_turbulent(4000, cool_Pr, d_i, l)    
    else:                                                                                                   
        Nu = calc_Nu_turbulent(Re, cool_Pr, d_i, l)

    return (Nu * cool_lambda) / d_i                                                                     

# ==========================================
# 2. THE MAIN API ENDPOINT
# ==========================================

def calc_overall_k(geo, ops, coolant_state, air_state, T_air_out: float) -> float:

    # ---------------------------------------------------------
    # 1. AIR SIDE (Outer)
    # ---------------------------------------------------------
    T_air_mean = (ops.T_air_in + T_air_out) / 2.0
    
    # Corrected velocity using the ops data natively
    w_eT = calc_w_eT(
        w_o=ops.w_o, 
        Ao_Ae_ratio=geo.Ao_Ae_ratio, 
        T_mean=T_air_mean, 
        T_in=ops.T_air_in
    )
    
    Nu_air = calc_Nu_air(
        A_ratio=(geo.A / geo.A_Go), 
        Pr_air=air_state.Pr, 
        d=geo.d, 
        w_eT=w_eT, 
        rho_air=air_state.rho, 
        eta_air=air_state.eta
    )
    
    alpha_R = calc_alpha_R(Nu_air, air_state.lambda_, geo.d)
    eta_R = calc_fin_efficiency(geo.D, geo.d, alpha_R, geo.lambda_R, geo.s)
    alpha_S = calc_alpha_S(alpha_R, eta_R, geo.A, geo.A_R)

    # ---------------------------------------------------------
    # 2. COOLANT SIDE (Inner)
    # ---------------------------------------------------------
    alpha_i = calc_alpha_i(
        w=ops.w_coolant, 
        d_i=geo.d_i, 
        l=geo.l, 
        cool_eta=coolant_state.eta, 
        cool_Pr=coolant_state.Pr, 
        cool_lambda=coolant_state.lambda_
    )
    
    # ---------------------------------------------------------
    # 3. OVERALL K-VALUE
    # ---------------------------------------------------------
    k_inv = (1 / alpha_S) + (geo.A / geo.A_i) * ((1 / alpha_i) + (geo.d - geo.d_i) / (2 * geo.lambda_R))
    
    return k_inv ** (-1)