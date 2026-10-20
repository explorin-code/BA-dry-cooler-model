# Inputs, GUI and instrumentation

## Inputs in mm / % and the GUI

- Geometry lengths are entered in **mm** in `parameters.py` (`*_MM`, e.g.
  `D_TUBE_OUTER_MM`, `FIN_PITCH_MM`, `HEIGHT_MM`) and converted to SI
  metres on read: `from_parameters` accepts `(CONSTANT_NAME, MM)` and
  divides by 1000. The model itself only ever sees metres.
- Relative humidity likewise in **%** (`PHI_AIR_PERCENT`, `(name, PERCENT)`,
  divided by 100); the model works with the fraction `ops.phi`, output shows %.
- Coolant: `COOLANT_TYPE` is any CoolProp fluid name (default 'Water');
  the GUI offers `fluid_properties.COOLANTS` as a dropdown (water, ethylene /
  propylene glycol 20-50 % as CoolProp INCOMP mixtures, checked 5-60 °C).
  Note: with water the coolant side is at Re ~3400, inside VDI G1's
  heat-transfer transition band 2300-4000 (Nu blended linearly between the
  laminar value at 2300 and the turbulent one at 4000; the friction factor
  uses its own 2320-3000 blend); 30 % glycol drops it to ~1700 (laminar),
  alpha_c 1584 -> 383 W/m²K, Q -37 %. Precooler water stays
  'Water' (precooling.calc_c_W) regardless of the coolant.
- GUI layout: groups in `gui.GROUP_ORDER` per mode (solver: Run, Operating conditions, Coolant, Precooling, Cooler, Output settings, Solver
  tuning; compare: Analyses, Comparison, Coolant, Precooling, Prices; the first ones open, the rest collapsed), booleans as toggles, the outputs /
  analyses as ONE multiselect each (from the catalogues), the 'exactly one nonzero' flow triples (M/W/V_COOLANT, W_O/M_O/V_O) as a 'given as'
  selectbox ABOVE the form (so the number field follows at once; its start value = the equivalent of parameters.py's default on the preset cooler,
  `flow_equivalents`) + one number, `N_TUBES` pre-filled with the preset's reference count, `COMPARE_N_TUBES` 0 = each cooler's own, the annual weather
  a selectbox over the full-year files of `climate_data/`.
- `gui.py` (Streamlit, `venv/bin/streamlit run gui.py`; the page is ~80 lines, its parts are the package `ui/`: `parameter_file` (reads parameters.py), `labels`, `presets`,
  `widgets`, `runner`, `results`; tests: `tests/test_gui_parts.py` + the AppTest runs in `test_smoke.py`) is the interactive
  alternative to `main.py`. It builds its form from `parameters.py`
  (values = defaults, inline comments = help, `# --- Group: ... ---` =
  groups), sets edited values on the `src.parameters` module for the
  session only (never writes the file), runs `run_scenarios`, and shows
  tables, log output and all figures. One run at a time across browser
  tabs (shared module state). Restart the Streamlit server after editing
  anything under `src/` -- it only reloads `gui.py` itself.
- Benchmark instrumentation: `src/profiling.py` -- `@timed('lookup' | 'coolprop' | 'k')`
  hooks on `get_fluid_properties` / `get_air_properties`, the CoolProp calls of `fluid_properties` and `calc_overall_k`, active only inside
  `profiling.record()`; the cache modes use `fluid_properties.caching(False)` / `clear_cache()` (a context manager that restores in `finally`).
  The hooks cost ~2-4 % of Cell's warm time in every run (measured: LMTD 0.11 -> 0.12 ms, NTU 4.6 -> 4.7, Cell 2D 20.2 -> 21.2, Cell 3D 648 -> 675 ms);
  call counts are identical to the patched version (LMTD 19 / 38, Cell 2D 2040 / 4082, Cell 3D 65280 / 130562 k / lookups).

## Notes behind the parameter comments

The comments in `parameters.py` are the GUI's help texts, so they are kept to one line; the longer notes:
- `W_O` 2.5 m/s is a normal operating speed (Towler's typical value); the cooler's design / top fan speed is 3 m/s (`COMPARE_W_FR_DESIGN`). The lower speed also
  triggers precooling at 27 degC, which is handy for testing.
- `N_TUBES = None`: the cell's `N_t_ref` -- 1 kg/s of 37 degC water to 25 degC at 20 degC / 60 %, 3 m/s -- fixed, not following the operating inputs.
- `CENTRAL_OMEGA` 0.3 (LMTD): stable up to 0.4, diverges from 0.5 (winter -10 degC) / 0.7 (everywhere); the loop is shared with NTU and Cell, which use
  their own `NTU_OMEGA` / `CELL_OMEGA` (model-notes.md, section 2). `CELL_OMEGA` 0.2 stopped ~0.02 K before convergence (Q +0.17 %).
- `CONVERGENCE_THRESHOLD`: both outlet temperatures change by less than it per iteration (NTU's inner field: 1/10 of it); keep >= ~1e-4 -- the property cache
  rounds T to 0.01 K, which leaves Cell a ~1.5e-5 K limit cycle.
- `PAD_IN_DRY_AIR_PATH`: precooled operation always passes the pad; in dry operation True counts its pressure drop, False = separate bypass inlet.
- `COPPER_EUR_PER_KG` / `ALUMINIUM_EUR_PER_KG`: LME cash prices of 9 Oct 2026 converted to EUR (cooler-comparison.md, cost model); no premium for tubes / fin stock.
- `COMPARE_W_FR_DESIGN`: above ~3 m/s droplet carry-over, noise and fin damage; Towler's typical face velocity is 2.5 m/s.
- `COMPARE_RE_TUBE`: > 4000 is out of the transition band; only the tubes / wall studies use it (the tube counts are the cells' `N_t_ref`).
