"""
build_climate_data.py
======================
Derives the project's climate files (climate_data/) from the raw DWD download (raw_climate_data/, not part of the
project -- it can be deleted once the files exist; to rebuild, put the DWD "stundenwerte_TF_<station>_akt" folder
back and run `python -m src.inputs.build_climate_data`).

  munich_past_year_<from>_to_<to>.csv     the last 8760 hours (365 days) of the data -- what every annual run uses
  munich_all_available_<from>_to_<to>.csv every hour of the download (gaps filled, flagged)
  munich_monthly_statistics.csv           per calendar month of all available hours: temperature, humidity, pressure
  munich_past_year_design_points.csv      the representative states of the past year (coldest, median, P95, P99, hottest, wettest warm hour)

Columns of the hourly files: utc, T_air_C [°C], phi_percent [%] (as measured, not capped), p_Pa [Pa] (station
pressure), interpolated (1 = the hour is missing in the download and was filled by linear interpolation).
"""

import logging
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent.parent
RAW_DIR = ROOT / 'raw_climate_data'
CLIMATE_DIR = ROOT / 'climate_data'
N_HOURS = 8760

logger = logging.getLogger(__name__)


def read_raw(raw_dir=RAW_DIR) -> pd.DataFrame:
    """The DWD file as a gap-free hourly frame (T_air_C, phi_percent, p_Pa, interpolated), index = UTC timestamps."""
    path = next(Path(raw_dir).glob('stundenwerte_TF_*/produkt_tf_stunde_*.txt'))
    raw = pd.read_csv(path, sep=';', skipinitialspace=True, na_values=['-999'])
    raw.index = pd.to_datetime(raw['MESS_DATUM'].astype(str), format='%Y%m%d%H')
    df = pd.DataFrame({'T_air_C': raw['TT_STD'], 'phi_percent': raw['RF_STD'], 'p_Pa': raw['P_STD'] * 100.0})
    df = df.reindex(pd.date_range(df.index[0], df.index[-1], freq='h'))
    df['interpolated'] = df['T_air_C'].isna().astype(int)
    df[['T_air_C', 'phi_percent', 'p_Pa']] = df[['T_air_C', 'phi_percent', 'p_Pa']].interpolate(limit_direction='both')
    df.index.name = 'utc'
    return df


def monthly_statistics(df: pd.DataFrame) -> pd.DataFrame:
    """Per calendar month: temperature (mean / min / max), hours below 0 °C and above 25 / 30 °C, mean humidity and pressure, interpolated hours."""
    g = df.groupby(df.index.to_period('M'))
    stats = pd.DataFrame({
        'hours': g.size(),
        'T_mean_C': g['T_air_C'].mean(), 'T_min_C': g['T_air_C'].min(), 'T_max_C': g['T_air_C'].max(),
        'hours_below_0C': g['T_air_C'].apply(lambda T: int((T < 0).sum())),
        'hours_above_25C': g['T_air_C'].apply(lambda T: int((T > 25).sum())),
        'hours_above_30C': g['T_air_C'].apply(lambda T: int((T > 30).sum())),
        'phi_mean_percent': g['phi_percent'].mean(), 'p_mean_hPa': g['p_Pa'].mean() / 100.0,
        'interpolated_hours': g['interpolated'].sum()})
    stats.index = stats.index.astype(str)
    stats.index.name = 'month'
    return stats.round(2)


def build(raw_dir=RAW_DIR, out_dir=CLIMATE_DIR):
    """Rebuilds the derived files of climate_data/ from the raw DWD download: the past year and all available hours, the monthly statistics and the design points. Returns the hourly table."""
    from src.inputs import weather                                        # design points use the project's definitions
    df = read_raw(raw_dir)
    out_dir.mkdir(exist_ok=True)
    for old in out_dir.glob('munich_*.csv'):
        old.unlink()
    span = lambda d: f"{d.index[0]:%Y-%m-%d}_to_{d.index[-1]:%Y-%m-%d}"
    past_year = df.iloc[-N_HOURS:]
    files = {f'munich_past_year_{span(past_year)}.csv': past_year,
             f'munich_all_available_{span(df)}.csv': df}
    for name, frame in files.items():
        frame.round({'T_air_C': 3, 'phi_percent': 3, 'p_Pa': 1}).to_csv(out_dir / name)
    monthly_statistics(df).to_csv(out_dir / 'munich_monthly_statistics.csv')
    points = weather.design_points(weather.add_wet_bulb(past_year.assign(phi_percent=past_year.phi_percent.clip(upper=weather.RH_MAX))))
    points.pop('design')                                           # a setting (COMPARE_DESIGN_*), not a property of the climate
    pd.DataFrame.from_dict(points, orient='index').rename_axis('case').to_csv(out_dir / 'munich_past_year_design_points.csv')
    for name in [*files, 'munich_monthly_statistics.csv', 'munich_past_year_design_points.csv']:
        logger.info("written: %s", out_dir / name)
    return df


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    build()
