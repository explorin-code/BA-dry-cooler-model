# climate_data -- the project's weather

Hourly weather of the DWD station **03379 München-Stadt** (48.163 °N, 11.543 °E, 515 m; times in **UTC**, never
converted). Source: Deutscher Wetterdienst, Climate Data Center, hourly station observations (`stundenwerte_TF_*`, the
"recent" package). The files are derived from the raw download by `python -m src.inputs.build_climate_data` (the raw folder
`raw_climate_data/` is not part of the project and can be deleted).

| file | content |
|---|---|
| `munich_past_year_2025-10-09_to_2026-10-08.csv` | the last 8760 hours (365 days) of the download -- the year every annual run uses (`src/inputs/weather.py`) |
| `munich_all_available_2025-04-07_to_2026-10-08.csv` | every hour of the download (about 18 months) |
| `munich_monthly_statistics.csv` | per calendar month of all available hours: temperature (mean / min / max), hours below 0 / above 25 / above 30 °C, mean humidity and pressure |
| `munich_past_year_design_points.csv` | the representative states of the past year: coldest, median, P95, P99, hottest, wettest warm hour |

Hourly columns: `utc`, `T_air_C` [°C], `phi_percent` [%] (as measured; the model caps it at 99 %), `p_Pa` [Pa] (station
pressure, ~955 hPa), `interpolated` (1 = the hour is missing in the DWD download and was filled by linear interpolation:
134 hours in all, 51 in the past year).
