from copy import deepcopy, copy
from typing import Dict, List, Any, Type, Tuple
import os
import pathlib

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

#####################################################################
##                      General (Base) Blocks                      ##
#####################################################################

# TODO: General Optimizer...

class Cuby4GeneralBlock(Cuby4Interface, SimpleBlock):
    def __init__(self, 
                 interface_config: Cuby4InterfaceConfig,
                 debug: bool = False, 
                 save_output: bool = False,
                 work_dir: str = ".",
                 cuby4_exe: str = "auto") -> None:
        
        interface_config = deepcopy(interface_config)

        Cuby4Interface.__init__(self, 
                                interface_config=interface_config,
                                cuby4_exe=cuby4_exe,
                                work_dir=work_dir)
        SimpleBlock.__init__(self, debug=debug, save_output=save_output)

    def operate(self, input_dict: Dict[str, Representation]) \
        -> Dict[str, Representation]:
            
            jc: Cuby4JobConfig = self.generate_job_config(input_dict)
            assert "job" in jc.config

            full_config = Cuby4MergedConfig.from_config_list([
                self.interface_config, jc
            ])

            output: str = self.run(full_config)
            return self.generate_out_dict(input_dict, output, full_config)

    def generate_job_config(self, input_dict: Dict[str, Representation]) \
        -> Cuby4JobConfig:
        raise NotImplementedError()
    
    def generate_out_dict(self,
                          input_dict: Dict[str, Representation],
                          cuby4_out: str, 
                          final_config: Cuby4Config) \
                            -> Dict[str, Representation]:
        raise NotImplementedError()

    @property
    def mol_rep(self) -> Type[Representation]:
        return ParmedMolRep
    
class Cuby4GeneralEnergyCalculator(Cuby4GeneralBlock):
    batch_groups = []

    @property
    def inputs(self):
        return [
            ("molecules", MoleculeData, self.mol_rep, False)
        ]
    
    @property
    def outputs(self):
        return [
            ("energy", NumericData, FloatRep, False)
        ]
    
    def generate_out_dict(self, 
                          input_dict,
                          cuby4_out: str, 
                          final_config) -> Dict[str, Representation]:
        energy = float(cuby4_out.split("Energy:")[-1].split()[0])
        return {"energy": self._get_out_rep("energy")(energy)}
    
class Cuby4GeneralEnergyOptimizer(Cuby4GeneralBlock):
    batch_groups = []

    @property
    def inputs(self):
        return [
            ("molecules", MoleculeData, self.mol_rep, False)
        ]
    
    @property
    def outputs(self):
        return [
            ("molecules", MoleculeData, None, False),
            ("energy", NumericData, FloatRep, False),
        ]
    
    def operate(self, input_dict: Dict[str, Representation]) \
        -> Dict[str, Representation]:
            
        jc: Cuby4JobConfig = self.generate_job_config(input_dict)
        assert "job" in jc.config
        assert jc.config["job"] == "optimize"

        # restart_file = generate_path_in_dir(6, self.base_dir, ".pdb")
        # jc.config["restart_file"] = restart_file
        # moved to generate_job_config method of individual blocks
        
        full_config = Cuby4MergedConfig.from_config_list([
            self.interface_config, jc
        ])

        output: str = self.run(full_config)
        energy = float(output.split("Energy:")[-1].split()[0])

        print(jc.config["restart_file"])
        restart_struct = parmed.load_file(jc.config["restart_file"])
        restart_struct: parmed.Structure
        mol_rep = deepcopy(input_dict["molecules"]) # TODO: Problem cause?
        # assert isinstance(mol_rep, Molecule3DRep)
        mol_rep: Molecule3DRep
        mol_rep.update_coordinates(restart_struct.coordinates)

        return {
            "molecules": mol_rep,
            "energy": self._get_out_rep("energy")(energy)
        }
    
    def generate_out_dict(self, 
                          input_dict: Dict[str, Representation],
                          cuby4_out: str, 
                          full_config: Cuby4Config) -> Dict[str, Representation]:
        energy = float(cuby4_out.split("Energy:")[-1].split()[0])

        parmed_rep: ParmedMolRep = input_dict[self.input_keys[0]]
        structure: parmed.Structure = parmed_rep.content
        structure = structure.copy(parmed.Structure) # TODO: works?
        restart: parmed.Structure = parmed.load_file(
            full_config.config["restart_file"])
        structure.coordinates = restart.coordinates
        parmed_rep = ParmedMolRep(structure)
        return {
            "molecules": parmed_rep,
            "energy": self._get_out_rep("energy")(energy)
        }