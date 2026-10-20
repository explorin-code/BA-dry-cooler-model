"""
fan_curve.py -- the variable-speed fan as a FUNCTION of the cooler inlet air temperature.

The hour-by-hour control (operation.control_point) needs a root-finding per hour; its result
collapses onto one monotone curve w_fr(T_a,i) per module, T_a,i = air entering the cooler
(ambient when dry, after the pad when precooled; humidity and station pressure only add a small
scatter). This module fits that curve (monotone PCHIP through binned medians of control.csv) and
turns it into a per-hour fan schedule without solving anything:

    T_amb, T_wb -> dry if T_amb <= T_dry_max (the curve reaches W_FR_DESIGN there), else the pad:
    T_a,i = T_amb - eta_B (T_amb - T_wb)  ->  w_fr = curve(T_a,i)  ->  P_fan = power_curve(T_a,i)

Approximations: the cooling limit is taken as the thermodynamic wet-bulb temperature (VDI M8
Gl. 15 differs by < 0.05 K, see docs/model-notes.md), eta_B from precooling.calc_eta_B, the curve is
valid for the weather/pressure range of the study year. `validate` compares the schedule with the
full control run.
"""

import numpy as np
import pandas as pd
from scipy.interpolate import PchipInterpolator

from src.inputs import weather
from src.thermo.precooling import calc_eta_B

BIN_K = 0.5


class _Branch:
    """w_fr(T_a,i) and P_fan(T_a,i) of one operating mode (dry or precooled): PCHIP through
    binned medians of the solved hours, monotone, valid between T_min and T_max."""
    def __init__(self, d: pd.DataFrame, w_min: float, w_max: float):
        self.w_min, self.w_max = w_min, w_max
        d = d.copy()
        d['bin'] = (d.T_a_i / BIN_K).round() * BIN_K
        g = d.groupby('bin').agg(w=('w_fr', 'median'), P=('P_fan', 'median'), n=('w_fr', 'size'))
        g = g[g.n >= 3]
        self.T = g.index.values
        self.w = np.maximum.accumulate(g.w.values)
        self.P = np.maximum.accumulate(g.P.values)
        self._w = PchipInterpolator(self.T, self.w, extrapolate=False)
        self._P = PchipInterpolator(self.T, self.P, extrapolate=False)
        self.T_min, self.T_max = self.T[0], self.T[-1]

    def w_fr(self, T):
        T = np.asarray(T, float)
        w = np.where(T < self.T_min, self.w_min, self._w(np.clip(T, self.T_min, self.T_max)))
        return np.clip(w, self.w_min, self.w_max)

    def P_fan(self, T):
        T = np.asarray(T, float)
        return np.where(T < self.T_min, self.P[0], self._P(np.clip(T, self.T_min, self.T_max)))


class FanCurve:
    """Fan curves of one module: separate branches for dry and precooled operation (the pad adds its
    pressure drop and the pad air is more humid, so P_fan is not one function of T_a,i)."""
    def __init__(self, cfg, control: pd.DataFrame, set_key: str):
        self.w_max = cfg.w_fr_design
        d = control[(control.set == set_key) & (~control.limited) & (control.w_fr > cfg.w_fr_min + 1e-3)]
        self.dry = _Branch(d[d['mode'] == 'dry'], cfg.w_fr_min, cfg.w_fr_design)
        self.pre = _Branch(d[d['mode'] == 'precooled'], cfg.w_fr_min, cfg.w_fr_design)
        # end of the dry range: T_a,i at which the dry curve reaches the top fan speed
        self.T_dry_max = float(np.interp(self.w_max, self.dry.w, self.dry.T))


def schedule(curve: FanCurve, year: pd.DataFrame) -> pd.DataFrame:
    """Fan schedule over a weather year: per hour mode, T_a,i, w_fr, P_fan, limited."""
    y = year if 'T_wb_C' in year else weather.add_wet_bulb(year)
    eta_B = calc_eta_B()
    T_amb, T_wb = y['T_air_C'].values, y['T_wb_C'].values
    dry = T_amb <= curve.T_dry_max
    T_pc = T_amb - eta_B * (T_amb - T_wb)
    T_a_i = np.where(dry, T_amb, T_pc)
    w = np.where(dry, curve.dry.w_fr(T_amb), curve.pre.w_fr(T_pc))
    P = np.where(dry, curve.dry.P_fan(T_amb), curve.pre.P_fan(T_pc))
    limited = ~dry & (T_pc > curve.pre.T_max)             # the pad is not enough even at the top fan speed
    w = np.where(limited, curve.w_max, w)
    P = np.where(limited, curve.pre.P[-1], P)
    return pd.DataFrame(dict(T_amb=T_amb, T_a_i=T_a_i, precooled=~dry, w_fr=w, P_fan=P, limited=limited), index=y.index)


def validate(curve: FanCurve, control: pd.DataFrame, set_key: str, year: pd.DataFrame) -> dict:
    """Schedule vs. the solved control run: annual energy and per-hour velocity error."""
    s = schedule(curve, year)
    c = control[control.set == set_key].copy()
    c['hour'] = pd.to_datetime(c['hour'])
    c = c.set_index('hour').reindex(s.index)
    err = (s.w_fr - c.w_fr).abs()
    return dict(E_curve_kWh=s.P_fan.sum() / 1e3, E_control_kWh=c.P_fan.sum() / 1e3,
                w_err_mean=err.mean(), w_err_p95=err.quantile(0.95), w_err_max=err.max())
