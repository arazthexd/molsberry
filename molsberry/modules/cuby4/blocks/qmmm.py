from copy import deepcopy, copy
from typing import Dict, List, Any, Type, Tuple

from rdkit import Chem
from openmm import unit
import parmed

from ....core import *
from ...mopac import MOPACInputMolRep
from ...openmm import OpenMMInputMolRep
from ...rdkit import RDKitMolRep
from ...parmed import ParmedMolRep

from ..configs import (
    Cuby4MOPACInterfaceConfig as C4MIC, 
    Cuby4EnergyJobConfig, 
    Cuby4MergedConfig,
    MOPACMozymeConfig,
    Cuby4OptimizeJobConfig,
    Cuby4Config,
    Cuby4InterfaceConfig,
    Cuby4JobConfig,
    Cuby4AMBERInterfaceConfig
)
from ..interface import Cuby4Interface

from .general import (
    Cuby4GeneralBlock,
    Cuby4GeneralEnergyCalculator,
    Cuby4GeneralEnergyOptimizer
)

#####################################################################
##                           QMMM Blocks                           ##
#####################################################################

class Cuby4QMMMEnergyCalculator(Cuby4GeneralEnergyCalculator):
    name = "c4qmmmspc"
    display_name = "Cuby4-QMMM Single Point Calculator"
    inputs = [
        ("qm_region", MoleculeData, ParmedMolRep, False),
        ("nonqm_region", MoleculeData, ParmedMolRep, False)
    ]
    outputs = [
        ("energy", NumericData, FloatRep, False),
    ]
    batch_groups = [("qm_region", "nonqm_region")]

    def __init__(self, 
                 interface_config,
                 debug: bool = False, 
                 save_output: bool = False,
                 work_dir: str = ".",
                 cuby4_exe: str = "auto") -> None:
        
        interface_config = deepcopy(interface_config)

        super().__init__(
            interface_config=interface_config,
            debug=debug,
            save_output=save_output,
            work_dir=work_dir,
            cuby4_exe=cuby4_exe
        )

    def generate_job_config(self, input_dict: Dict[str, ParmedMolRep]) \
        -> Cuby4JobConfig:

        struct_qm = input_dict["qm_region"].content
        struct_nonqm = input_dict["nonqm_region"].content
        struct_whole: parmed.Structure = struct_qm + struct_nonqm

        # geometry, core
        geometry = generate_path_in_dir(6, self.base_dir, ".pdb")
        struct_whole.save(geometry)

        core = f"1-{len(struct_qm.atoms)}"
        charge_whole = round(sum([atom.charge for atom in struct_whole.atoms]))

        config = Cuby4EnergyJobConfig(geometry=geometry, charge=charge_whole)
        config.config["qmmm_core"] = core

        # calculation_qm: charge
        charge = round(sum([atom.charge for atom in struct_qm.atoms]))
        config.config["calculation_qm"] = {}
        config.config["calculation_qm"]["charge"] = charge

        # calculation_mm: parm7 mereged
        parm7_whole = generate_path_in_dir(6, self.base_dir, ".parm7")
        struct_whole.save(parm7_whole)

        config.config["calculation_mm"] = {}
        config.config["calculation_mm"]["amber_top_file"] = parm7_whole

        # calculation_qmregion_mm: parm7 lig
        parm7_qm = generate_path_in_dir(6, self.base_dir, ".parm7")
        struct_qm.save(parm7_qm)

        config.config["calculation_qmregion_mm"] = {"amber_top_file": parm7_qm}

        return config