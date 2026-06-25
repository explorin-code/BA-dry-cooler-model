# implementation of different solvers (cell method, correction factors, etc.)
from math import ln

def solve_it_LMTD(T_coolant_in: float, T_air_in: float):

    dT_in = 10
    dT_out = 10

    T_coolant_out = T_air_in + 10          # first assumption about exiting coolant temp
    T_air_out = T_coolant_in + 10           # first assumption about exiting air temp

    LMTD: float
    LMTD_it: float




    return

def calc_LMTD(dT_in: float, dT_out):

    LMTD = (dT_in-dT_out)/ln(dT_in/dT_out)

    return LMTD





