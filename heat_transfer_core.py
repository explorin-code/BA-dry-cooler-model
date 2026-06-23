# basic heat balances

def calc_alpha_S(alpha_R: float, eta_R: float, A: float, A_R: float):
    alpha_S = alpha_R * (1 - (1 - eta_R) * A_R/A)
    return alpha_S

def calc_w_e(w_o, Ao_Ae_ratio):
    w_e = w_o * Ao_Ae_ratio
    return w_e

def calc_w_eT(w_o, Ao_Ae_ratio, T_upper, T_lower):
    w_eT = calc_w_e(w_o, Ao_Ae_ratio) * ((273.15 + 0.5 * (T_upper - T_lower)) / (273.15 + T_lower))
    return w_eT

def calc_alpha_R(Nu: float, lambda_air: float, d: float): # mean alpha of fin and pipes
    alpha_R = (Nu * lambda_air) / d
    return alpha_R

def calc_heat_transfer_coeff(alpha_S: float, A: float, A_i: float, alpha_i: float, d: float, d_i: float, lambda_R: float):
    k = (1/alpha_S + (A/A_i)*(1/alpha_i + (d-d_i)/(2*lambda_R)))**(-1)
    return k
