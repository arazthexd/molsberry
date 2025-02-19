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