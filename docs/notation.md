# Notation

## Naming conventions (aligned with the thesis notation)

- Temperatures: `T_<medium>_<state>` in °C (`T_c_i`, `T_a_o`, `T_K`,
  `T_a_pc`, grids `T_c_grid`/`T_a_grid`, `Delta_T_hot`/`Delta_T_cold`).
  theta / vartheta are not used in the code or the thesis equations; EXCEPTION: the normalised profile figures
  (`profile_grid_normalised`, `profile_collapsed_normalised`) label their axis theta = (T - T_a,i)/(T_c,i - T_a,i); it
  is listed in the
  thesis Notation chapter as 'normalised temperature (figures)'.
  Kelvin only where absolute temperatures are needed, marked `_K` /
  `T_kelvin` (`T_m_K`, `T_i_K` in `calc_w_e_m`). The Cabezas-Gomez element
  functions use T_hot/T_cold (paper: T_h/T_c) so T_c always means coolant.
- Indices: c coolant, a air, i in, o out, m mean, pc precooled, K cooling
  limit, da dry air, ev evaporation, f fin, p pipe, fr front/face,
  e narrowest cross-section. Humidity ratio: `Y` (ops.Y, Y_i, Y_pc, Y_K).
- **Notation scheme**: one letter per kind of quantity,
  the index says which; lengths lower case (ISO 80000-3 / DIN 1304),
  P = power, p = pressure, s = pitch, delta = thickness, d = diameter.
- Geometry: `d` (tube outer diameter, no index -- o = outflow), `d_i`,
  `d_c` (collar diameter, d + 2 delta_c), `delta_f`/`delta_c` (fin / collar
  thickness), `s_t`/`s_l`/`s_f` (transverse / longitudinal tube pitch, fin
  pitch; Wang's P_t/P_l/F_p/D_c), clear fin spacing has no symbol (inline
  s_f - delta_f), `N_t`/`N_r` (tubes per row / rows = passes; N_r kept as
  "rows", not N_l), `N_f` (fins per pass), `l` (one pass), `l_tot` (all
  passes), `b` (width), `A_f`/`A_p`/`A_p0` (fin / exposed pipe / bare pipe
  area per tube pass; A_p = pi d (l - N_f delta_f)), `A_tot`, `A_fr`,
  `A_e`, `sigma` (= A_e/A_fr), `A_cs_c`, `lambda_f`/`lambda_p`, `delta_p`
  (tube wall thickness), `b_f`/`l_f` (fin rectangle, staggered, VDI M1
  Eq. 14, part of Geometry).
- Flows/powers: `w_c`, `w_fr`, `m_dot_c`, `m_dot_a`, `V_dot_c`,
  `V_dot_a`, `p_c`/`p_a` (pressures), `Q_dot`, `C_dot_c`/`C_dot_a`
  (capacity rates), `P_pump`/`P_fan` (powers), `m_dot_ev`, `Delta_p_*`.
- Heat-transfer chain (side index a = air, c = coolant):
  `w_e_m` (narrowest-gap velocity at the mean air temperature; VDI w_eT),
  `Re_a`/`Pr_a`/`Nu_a` (VDI M1, row factor `C_r`), `alpha_a` (VDI alpha_R),
  fin efficiency `eta_f` via `phi_f_0`/`phi_f` (from Geometry.b_f/l_f) (VDI phi -- index
  f so it is not the relative humidity), `X`, `alpha_a_eff` (incl. fin
  efficiency; VDI "scheinbar" alpha_S); `Re_c`/`Pr_c`/`Nu_c` (VDI G1,
  `calc_Nu_c_lam`/`calc_Nu_c_turb`, blend factor gamma, friction psi),
  `alpha_c` (VDI alpha_i), coolant-side area `A_c` (was A_i); k with
  1/k = 1/alpha_a_eff + A/A_c (1/alpha_c + delta_p/lambda_p),
  delta_p = Geometry.delta_p. Wang's friction Re is `Re_d_c` (length d_c).
  Diagnostics dict keys follow (`Re_a`, `Nu_c`, `alpha_c`, ...).
- Deliberately NOT renamed: `mu`, and the Cabezas-Gomez paper notation
  inside the NTU element functions (T_hot/T_cold, i.e. coolant/air;
  B there is the paper's B, not a width).
- Fan and pump are spelled out to keep f = fin and p = pipe unambiguous:
  `eta_fan`/`P_fan`, `eta_pump`/`P_pump` (thesis Eqs. 2.6-2.8 write
  eta_f/W_f and eta_p/W_p).

## Decisions

- Temperatures are T everywhere (°C unless absolute K is needed), no vartheta; the one exception is theta in the
  normalised profile figures.
- Index i is inflow, except in `d_i` (inner diameter); the coolant-side quantities carry c (`alpha_c`, `A_c`). The
  cooler width is `b` (not B, which is Cabezas-Gómez's B in the NTU element functions).
- Friction factors are all written f in the thesis (coolant f_c: Darcy 64/Re, Blasius, blend; air f_a: Wang, Fanning).
  The code's `psi` (Gnielinski / Konakov in `calc_Nu_c_turb`) and `f_c` / `f` in
  `pressure_drop.py` are internal names.

## Figure labels aligned with the thesis Notation
- Result boxes (`src/workflow/run_scenario.py`): `N_r`, `N_t`, `l x b` (was n_rows, n_tubes, H x W), `A_fr` (was A_front),
  `P_tot` (was P_total),
  `m_dot_ev` (was m_dot_water; also `ui/results.py`, `comparison/plots_common.LABEL`).
- Heat-flow profiles: `d Q_dot / d l` and `Q_dot / l_tot` (was dQ/dL, Q/L; `L` is the service life in the cost
  equations).
- Cost model symbols of the appendix (thesis `Appendix.tex`): C (annual cost), c (price), E (energy), m (mass), V
  (volume), L (service life),
  t (time), summation index j (not h = enthalpy); indices amb, pad, pinch, Al, Cu.
- Thesis symbols that exist only there (not in the code): beta_Y (mass-transfer coefficient), gamma (laminar/turbulent
  blend), F_1..F_3 (Wang
  exponents), C* (heat capacity ratio), Me (Merkel number), z (tower height coordinate), primes ' / '' (thermodynamic /
  saturated state).
