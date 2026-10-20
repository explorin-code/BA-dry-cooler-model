"""
weather.py
===========
The climate data of the project: the hourly weather of the DWD station 03379 München-Stadt, as the derived files in
climate_data/ (src/inputs/build_climate_data.py makes them from the raw DWD download) -> the past year as hourly (T, phi, p),
plus representative design points (coldest / median / hottest hour ...) picked from it. Used by BOTH workflows: the
annual solver benchmark (benchmark_annual.py) and the cooler comparison (comparison/).

Timestamps are UTC (DWD convention) and are never converted: only the values are used, not clock times.
Gaps in the hourly series were filled by linear interpolation when the files were built. Relative humidity is capped at
RH_MAX when a year is loaded.
"""

from functools import lru_cache
from pathlib import Path
import logging

import pandas as pd
from CoolProp.HumidAirProp import HAPropsSI

logger = logging.getLogger(__name__)

CLIMATE_DIR = Path(__file__).resolve().parent.parent.parent / 'climate_data'
N_HOURS = 8760
RH_MAX = 99.0     # [%] cap: at 100 % the property cache's rounding of Y (6 digits) gives phi > 1 in CoolProp


@lru_cache(maxsize=None)
def load_year() -> pd.DataFrame:
    """The past year (climate_data/munich_past_year_*.csv, 8760 hours) as a gap-free hourly DataFrame, index = UTC
    timestamp, columns T_air_C [°C], phi_percent [%] (capped at RH_MAX), p_Pa [Pa] (station pressure). The result is
    cached -- do not modify it. Attribute `n_filled` = number of interpolated hours."""
    path = next(CLIMATE_DIR.glob('munich_past_year_[0-9]*.csv'))
    df = pd.read_csv(path, index_col='utc', parse_dates=True)
    df.attrs['n_filled'] = int(df.pop('interpolated').sum())
    df['phi_percent'] = df['phi_percent'].clip(upper=RH_MAX)
    logger.info("[Weather] %s: %s .. %s, %d hours (%d interpolated)", path.name, df.index[0], df.index[-1], len(df), df.attrs['n_filled'])
    return df


def climate_year(csv_path=None) -> pd.DataFrame:
    """The year used by the annual runs: the DWD Munich year, or -- if csv_path is given -- a CSV with the
    columns T_air_C, phi_percent (and optionally p_Pa, default station mean), 8760 hourly rows."""
    if csv_path is None:
        return load_year()
    df = pd.read_csv(csv_path)
    if len(df) != N_HOURS:
        raise ValueError(f"{csv_path}: expected {N_HOURS} hourly rows, found {len(df)}.")
    if 'p_Pa' not in df:
        df['p_Pa'] = load_year()['p_Pa'].mean()
    return df[['T_air_C', 'phi_percent', 'p_Pa']]


def hourly_states(df: pd.DataFrame, stride: int = 1):
    """The operating points of a weather year: [(timestamp, {'T_a_i': [°C], 'phi': [0-1], 'p_a': [Pa]}), ...], every
    stride-th hour -- what the annual benchmark and the cooler comparison's year analyses feed to the solvers."""
    return [(ts, dict(T_a_i=float(row['T_air_C']), phi=float(row['phi_percent']) / 100.0, p_a=float(row['p_Pa'])))
            for ts, row in df.iloc[::stride].iterrows()]


def add_wet_bulb(df: pd.DataFrame) -> pd.DataFrame:
    """Adds T_wb_C (thermodynamic wet-bulb, CoolProp) to the frame."""
    T_wb = [HAPropsSI('Twb', 'T', T + 273.15, 'R', phi / 100.0, 'P', p) - 273.15
            for T, phi, p in zip(df['T_air_C'], df['phi_percent'], df['p_Pa'])]
    out = df.copy()
    out['T_wb_C'] = T_wb
    return out


def design_points(df: pd.DataFrame, design_T=20.0, design_phi=0.50) -> dict:
    """name -> dict(T_a_i, phi, p_a, note). Real hours (coldest, hottest,
    most humid-warm = highest wet-bulb) plus percentile/median points whose
    humidity and pressure are the medians of the hours within +-0.5 K of
    that temperature. 'design' is the design ambient (default 20 °C / 50 %, mean station pressure)."""
    if 'T_wb_C' not in df:
        df = add_wet_bulb(df)
    p_mean = float(df['p_Pa'].mean())

    def real_hour(label, idx, note):
        r = df.loc[idx]
        return dict(T_a_i=float(r['T_air_C']), phi=float(r['phi_percent']) / 100, p_a=float(r['p_Pa']),
                    note=f"{note} ({idx:%Y-%m-%d %H}h UTC)")

    def at_percentile(q, note):
        T = float(df['T_air_C'].quantile(q))
        near = df[(df['T_air_C'] - T).abs() <= 0.5]
        return dict(T_a_i=T, phi=float(near['phi_percent'].median()) / 100, p_a=float(near['p_Pa'].median()),
                    note=f"{note} (T at the {q * 100:g} % quantile, median phi of nearby hours)")

    return {
        'coldest':   real_hour('coldest', df['T_air_C'].idxmin(), "coldest hour"),
        'median':    at_percentile(0.50, "median"),
        'P95':       at_percentile(0.95, "95 % quantile"),
        'P99':       at_percentile(0.99, "99 % quantile"),
        'hottest':   real_hour('hottest', df['T_air_C'].idxmax(), "hottest hour"),
        'wettest_warm': real_hour('wettest', df['T_wb_C'].idxmax(), "highest wet-bulb hour"),
        'design':    dict(T_a_i=design_T, phi=design_phi, p_a=p_mean,
                          note=f"design ambient: {design_T:g} °C / {design_phi * 100:g} % RH"),
    }


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    year = add_wet_bulb(load_year())
    print(year.describe().round(1))
    for name, pt in design_points(year).items():
        print(f"{name:13s} T={pt['T_a_i']:5.1f} °C  phi={pt['phi'] * 100:4.0f} %  p={pt['p_a'] / 100:6.1f} hPa  {pt['note']}")
