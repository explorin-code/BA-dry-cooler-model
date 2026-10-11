"""
cooler_sets.py
===============
The standard "elementary cells" P1-P3 of the dry cooler: tube and fin dimensions, pitches and materials
of three plain plate-fin-and-tube coils taken from the literature. Only the elementary cell is from the
sources; the module around it (tubes per row, tube length, rows) is chosen by the run:

  * single-cooler runs take a cell as the default geometry (parameters.GEOMETRY_SET) and may override
    any value (dry_cooler_physics.get_geometry);
  * cooler comparisons size every cell to the same water flow and race them against each other
    (comparison/).

Source of the values: "Project-Source Geometry Audit for Plain Plate-Fin Dry-Cooler Cases" (Perplexity
report, 2026-10-09) -- NOT re-read in the original PDFs. What the sources do not give is an ASSUMPTION
and marked as such in `notes`: tube wall (0.5 mm; P3 1 mm), collar thickness (= fin thickness: the collar is
the fin sheet turned up around the tube), materials of P1.

P1: the diameter listed by Wang & Chi is the OUTER DIAMETER OVER THE COLLAR (Wang's D_c), so the tube's own outer diameter
is d = listed - 2 delta_c (`_tube_od`); P2/P3 (Kayansayan): the listed value is taken as the tube outer diameter d, the
collar adds on top (d_c = d + 2 delta_c).
"""

from dataclasses import dataclass


@dataclass(frozen=True)
class CellSet:
    key: str
    label: str
    source: str
    d_mm: float                       # tube outer diameter
    delta_f_mm: float                 # fin thickness
    s_f_mm: float                     # fin pitch (spacing + thickness)
    s_t_mm: float                     # transverse tube pitch
    s_l_mm: float                     # longitudinal tube pitch
    wall_mm: float = 0.5              # tube wall thickness -> d_i = d - 2 wall (not reported: assumption)
    delta_c_mm: float = 0.0           # fin-collar thickness (only enters Wang's d_c)
    pipe_material: str = 'Copper'
    fin_material: str = 'Aluminum'
    N_r_source: int = 0               # rows of the source coil (information only)
    N_r: int = 6                      # rows of the module built from the cell (comparison runs set this)
    N_t_ref: int = 0                  # tubes per row of the reference module (single-cooler runs, see below)
    notes: str = ""

    @property
    def d_i_mm(self) -> float:
        return self.d_mm - 2 * self.wall_mm


def _tube_od(listed_mm: float, delta_c_mm: float) -> float:
    """Tube outer diameter from the listed diameter over the collar (d_c = d + 2 delta_c)."""
    return round(listed_mm - 2 * delta_c_mm, 4)


# N_t_ref: the tube count a single-cooler (solver comparison) run uses by default, FIXED so the geometry does not change with the
# operating inputs (coolant temperature, flow, ...): the smallest count that brings 1 kg/s water of 37 degC to 25 degC at 20 degC / 60 %,
# dry, 3 m/s, 6 rows, 2 m tubes (comparison.sizing.design_tubes, 2026-10; the cooler comparison itself re-sizes live).
CELL_SETS = {
    'P1': CellSet(
        'P1', "P1 dense (Wang sample 2)", "Wang & Chi (2000), Part I, Table 2, sample 2",
        d_mm=_tube_od(7.3, 0.115), delta_f_mm=0.115, s_f_mm=1.22, s_t_mm=21.0, s_l_mm=12.4, delta_c_mm=0.115, N_r_source=2, N_t_ref=34,
        notes="lab air-side coil; listed d 7.3 mm is over the collar; tube ID, fin material and circuit not reported (aluminium fins, 0.5 mm wall assumed)"),
    'P2': CellSet(
        'P2', "P2 medium (Kayansayan type 5)", "Kayansayan (1994), Table 1, coil 5",
        d_mm=9.52, delta_f_mm=0.2, s_f_mm=1000 / 312, s_t_mm=30.0, s_l_mm=26.0, delta_c_mm=0.2, N_r_source=4, N_t_ref=32,
        notes="16 tubes/row, 4 parallel circuits in the source; collar (= fin thickness) assumed"),
    'P3': CellSet(
        'P3', "P3 open (Kayansayan type 3)", "Kayansayan (1994), Table 1, coil 3",
        d_mm=16.3, delta_f_mm=0.2, s_f_mm=1000 / 238, s_t_mm=40.0, s_l_mm=34.67, wall_mm=1.0, delta_c_mm=0.2, N_r_source=4, N_t_ref=30,
        notes="12 tubes/row; d, s_t, s_l above Wang's validity range (friction factor extrapolated); wall 1 mm and collar (= fin thickness) assumed"),
}

# An industrial reference cell (d 25.4 mm, s_t 63.5 mm, 400 fins/m; Kroger / Towler) was tried and dropped
# (individually finned circular-fin tubes, outside Wang's range -- a different problem).
