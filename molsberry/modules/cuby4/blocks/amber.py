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
##                          AMBER Blocks                           ##
#####################################################################

class _AMBER_Utils:
    @staticmethod
    def omm_to_parms(omm_rep: OpenMMInputMolRep, base_dir: str) \
        -> Tuple[str, str, str]:
        struct: parmed.Structure = parmed.openmm.load_topology(
            omm_rep.topology,
            omm_rep.system,
            omm_rep.positions
        )
        geometry = generate_path_in_dir(6, base_dir, ".pdb")
        struct.save(geometry)

        prmtop = generate_path_in_dir(6, base_dir, ".parm7")
        struct.save(prmtop)

        return struct, geometry, prmtop
    
    @staticmethod
    def parmed_to_parms(parmed_rep: ParmedMolRep, base_dir: str) \
        -> Tuple[str, str, str]:
        struct = parmed_rep.content

        geometry = generate_path_in_dir(6, base_dir, ".pdb")
        struct.save(geometry)

        prmtop = generate_path_in_dir(6, base_dir, ".parm7")
        struct.save(prmtop)

        return struct, geometry, prmtop

class Cuby4AMBEREnergyCalculator(Cuby4GeneralEnergyCalculator):
    name = "c4amberspc"
    display_name = "Cuby4-AMBER Single Point Calculator"

    def __init__(self, 
                 interface_config = None,
                 debug: bool = False, 
                 save_output: bool = False,
                 work_dir: str = ".",
                 cuby4_exe: str = "auto") -> None:
        
        if interface_config is None:
            interface_config = Cuby4AMBERInterfaceConfig()
        interface_config = deepcopy(interface_config)
        
        super().__init__(
            interface_config=interface_config,
            debug=debug,
            save_output=save_output,
            work_dir=work_dir,
            cuby4_exe=cuby4_exe
        )

    def generate_job_config(self, input_dict):
        omm_rep: ParmedMolRep = input_dict[self.input_keys[0]]
        _, geometry, prmtop = _AMBER_Utils.parmed_to_parms(omm_rep, 
                                                           self.base_dir)

        config = Cuby4EnergyJobConfig(
            geometry=geometry
        )
        
        mol_specifics = {
            "amber_top_file": prmtop
        }
        config.config.update(mol_specifics)

        return config
    
class Cuby4AMBEREnergyOptimizer(Cuby4GeneralEnergyOptimizer):
    name = "c4amberopt"
    display_name = "Cuby4-AMBER Energy Optimizer"

    def __init__(self, 
                 interface_config = None,
                 debug: bool = False, 
                 save_output: bool = False,
                 work_dir: str = ".",
                 cuby4_exe: str = "auto") -> None:
        
        if interface_config is None:
            interface_config = Cuby4AMBERInterfaceConfig()
        interface_config = deepcopy(interface_config)
        
        super().__init__(
            interface_config=interface_config,
            debug=debug,
            save_output=save_output,
            work_dir=work_dir,
            cuby4_exe=cuby4_exe
        )

    def generate_job_config(self, input_dict):
        omm_rep: ParmedMolRep = input_dict[self.input_keys[0]]
        _, geometry, prmtop = _AMBER_Utils.parmed_to_parms(omm_rep, 
                                                           self.base_dir)

        restart_file = generate_path_in_dir(6, self.base_dir, ".pdb")
        config = Cuby4OptimizeJobConfig(
            geometry=geometry,
            restart_file=restart_file
        )
        
        mol_specifics = {
            "amber_top_file": prmtop
        }
        config.config.update(mol_specifics)

        return config