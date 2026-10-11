# Notation

## Naming conventions (renamed 2026-10, aligned with the thesis notation)

- Temperatures: `T_<medium>_<state>` in °C (`T_c_i`, `T_a_o`, `T_K`,
  `T_a_pc`, grids `T_c_grid`/`T_a_grid`, `Delta_T_hot`/`Delta_T_cold`).
  theta/vartheta dropped entirely (2026-10), in code, figures and thesis.
  Kelvin only where absolute temperatures are needed, marked `_K` /
  `T_kelvin` (`T_m_K`, `T_i_K` in `calc_w_e_m`). The Cabezas-Gomez element
  functions use T_hot/T_cold (paper: T_h/T_c) so T_c always means coolant.
- Indices: c coolant, a air, i in, o out, m mean, pc precooled, K cooling
  limit, da dry air, ev evaporation, f fin, p pipe, fr front/face,
  e narrowest cross-section. Humidity ratio: `Y` (ops.Y, Y_i, Y_pc, Y_K).
- **Notation scheme (option A, 2026-10)**: one letter per kind of quantity,
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
  Eq. 14 -- moved from the fin-efficiency function into Geometry 2026-10).
- Flows/powers: `w_c`, `w_fr`, `m_dot_c`, `m_dot_a`, `V_dot_c`,
  `V_dot_a`, `p_c`/`p_a` (pressures), `Q_dot`, `C_dot_c`/`C_dot_a`
  (capacity rates), `P_pump`/`P_fan` (powers), `m_dot_ev`, `Delta_p_*`.
- Heat-transfer chain (renamed 2026-10, side index a = air, c = coolant):
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
  eta_f/W_f and eta_p/W_p -- to be aligned there).

## Open notation decisions (thesis "Notation" chapter, 2026-10)

Keep in mind -- not yet decided / deliberately left for now:
- ~~T vs. theta~~ -- resolved: T everywhere (°C unless absolute K needed).
- ~~Index i = inflow and inner~~ -- mostly resolved: alpha_i -> alpha_c,
  A_i -> A_c; only d_i (inner diameter) keeps i.
- ~~B = cooler width vs. Cabezas-Gomez B~~ -- resolved: the width is now `b`.

- **Friction factors (thesis, 2026-10)**: the thesis writes ALL friction
  factors as f -- coolant f_c (Darcy: 64/Re, Blasius, blend), air f_a
  (Wang, Fanning). Not psi (so the notation's psi entry is obsolete). The
  code's `psi` (Gnielinski/Konakov in `calc_Nu_c_turb`) and `f_c`/`f` in
  pressure_drop.py are internal names; keep thesis and code in mind when
  writing.
- **Temperatures**: thesis uses T everywhere (°C unless absolute K is
  needed); vartheta dropped.
