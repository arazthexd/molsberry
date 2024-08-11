CUBY4_MAIN_TEMPLATE = """
geometry: {geometry}
cuby_threads: {cuby_threads}"""

CUBY4_INTERACTION_ADDON = """
job: interaction
molecule_a:
  selection: {selection_a}
  {optional_a}
molecule_b:
  selection: {selection_b}
  {optional_b}"""

CUBY4_OPTIMIZE_ADDON = """
job: optimize
maxcycles: {max_cycles}
restart_file: {restart_file}"""

CUBY4_MOPAC_MAIN_ADDON = """
interface: mopac
method: {method}
mopac_exe: {mopac_exe}
mopac_mozyme: {mopac_mozyme}
mopac_corrections: {mopac_corrections}
solvent_model: {solvent_model}"""

CUBY4_MOPAC_MOLSPEC_ADDON = """
mopac_setcharge:
  {mopac_setcharge}
mopac_setpi: [{mopac_setpi}]
charge: {charge}"""

CUBY4_QMMM_MAIN_ADDON = """
interface: qmmm
qmmm_core: {qmmm_core}
qmmm_embedding: {qmmm_embedding}
gradient_on_point_charges: {grad_on_points}
calculation_qm:
  gradient_on_point_charges: {grad_on_points}
  {qm_config}
calculation_mm:
  {mm_config}
calculation_qmregion_mm:
  {qmregion_mm_config}"""

CUBY4_AMBER_MAIN_ADDON = """
interface: amber
amber_amberhome: {amber_home}
{topology_config}"""