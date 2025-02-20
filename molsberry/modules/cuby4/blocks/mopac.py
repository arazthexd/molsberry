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

from .general import Cuby4GeneralEnergyCalculator, Cuby4GeneralBlock

#####################################################################
##                          MOPAC Blocks                           ##
#####################################################################
    
class Cuby4MOPACEnergyCalculator(Cuby4GeneralEnergyCalculator):
    name = "c4mopacspc"
    display_name = "Cuby4-MOPAC Single Point Calculator"

    def __init__(self, 
                 interface_config: C4MIC = C4MIC(),
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

        self.interface_config: C4MIC

    def generate_job_config(self, 
                            input_dict: Dict[str, ParmedMolRep]) -> C4MIC:
        parmed_rep: ParmedMolRep = input_dict[self.input_keys[0]]
        rdkit_rep: RDKitMolRep = parmed_rep.to_RDKitMolRep()
        mopac_rep = MOPACInputMolRep.from_RDKitMolRep(rdkit_rep) # TODO: Check

        jc = Cuby4EnergyJobConfig(
            geometry=generate_path_in_dir(6, self.base_dir, "_c4mecinp.pdb"),
            charge=mopac_rep.charge,
            mult=1 # TODO: Others..?
        )

        if self.interface_config.mozyme:
            jc = Cuby4MergedConfig.from_config_list([
                jc, MOPACMozymeConfig(
                    setpi=mopac_rep.setpi[:50], # TODO: 50 :(
                    neg_cvb=mopac_rep.neg_cvb
                )
            ])

            if self.interface_config.setpi == False:
                if "mopac_setpi" in jc.config:
                    _ = jc.config.pop("mopac_setpi")

            if self.interface_config.cvb == False:
                if "mopac_keywords" in jc.config:
                    kwds = jc.config["mopac_keywords"].split()
                    kwds = [kwd for kwd in kwds if "CVB" not in kwd]
                    jc.config["mopac_keywords"] = " ".join(kwds)

        with open(jc.config["geometry"], "w") as f:
            f.write(mopac_rep.coordinates)

        if self.interface_config.setcharges:
            jc.config["mopac_setcharge"] = {
                i+1: c for i, c in enumerate(mopac_rep.atom_charges)}

        return jc
    
    def generate_out_dict(self, 
                          input_dict,
                          cuby4_out: str, 
                          final_config) -> Dict[str, Representation]:
        energy = float(cuby4_out.split("Energy:")[-1].split()[0])
        return {"energy": self._get_out_rep("energy")(energy)}

class Cuby4MOPACEnergyOptimizer(Cuby4GeneralBlock): # TODO: general optimizer
    name = "c4mopacopt"
    display_name = "Cuby4-MOPAC Energy Optimizer"
    inputs = [
        ("molecules", MoleculeData, ParmedMolRep, False)
    ]
    outputs = [
        ("molecules", MoleculeData, ParmedMolRep, False),
        ("energy", NumericData, FloatRep, False)
    ]
    batch_groups = []

    def __init__(self, 
                 interface_config: C4MIC = C4MIC(),
                 debug: bool = False, 
                 save_output: bool = False,
                 work_dir: str = ".",
                 cuby4_exe: str = "auto",
                 max_cycles: int = 10) -> None:
        
        interface_config = deepcopy(interface_config)

        super().__init__(
            interface_config=interface_config,
            debug=debug,
            save_output=save_output,
            work_dir=work_dir,
            cuby4_exe=cuby4_exe
        )

        self.interface_config: C4MIC

        self.max_cycles = max_cycles

    def generate_job_config(self, 
                            input_dict: Dict[str, ParmedMolRep]) -> C4MIC:
        # TODO: Mix with other MOPAC job config generation methods

        parmed_rep: ParmedMolRep = input_dict[self.input_keys[0]]
        rdkit_rep: RDKitMolRep = parmed_rep.to_RDKitMolRep()
        mopac_rep = MOPACInputMolRep.from_RDKitMolRep(rdkit_rep) # TODO: Check

        jc = Cuby4OptimizeJobConfig(
                geometry=generate_path_in_dir(6, 
                                              self.base_dir, 
                                              "_c4meoinp.pdb"),
                charge=mopac_rep.charge,
                mult=1, # TODO
                restart_file=generate_path_in_dir(6, 
                                                  self.base_dir, 
                                                  "_c4meoout.pdb"),
                maxcycles=self.max_cycles
            )

        if self.interface_config.mozyme:
            jc = Cuby4MergedConfig.from_config_list([
                jc, MOPACMozymeConfig(
                    setpi=mopac_rep.setpi[:50], # TODO: 50 :(
                    neg_cvb=mopac_rep.neg_cvb
                )
            ])

            if self.interface_config.setpi == False:
                if "mopac_setpi" in jc.config:
                    _ = jc.config.pop("mopac_setpi")

            if self.interface_config.cvb == False:
                    kwds = jc.config["mopac_keywords"].split()
                    kwds = [kwd for kwd in kwds if "CVB" not in kwd]
                    jc.config["mopac_keywords"] = " ".join(kwds)

        with open(jc.config["geometry"], "w") as f:
            f.write(mopac_rep.coordinates)

        if self.interface_config.setcharges:
            jc.config["mopac_setcharge"] = {
                i+1: c for i, c in enumerate(mopac_rep.atom_charges)}

        return jc
    
    def generate_out_dict(self, 
                          input_dict: Dict[str, Representation],
                          cuby4_out: str, 
                          full_config: Cuby4Config) -> Dict[str, Representation]:
        energy = float(cuby4_out.split("Energy:")[-1].split()[0])

        parmed_rep: ParmedMolRep = input_dict[self.input_keys[0]]
        structure: parmed.Structure = parmed_rep.content
        structure = structure.copy(parmed.Structure)
        restart: parmed.Structure = parmed.load_file(
            full_config.config["restart_file"])
        structure.coordinates = restart.coordinates
        parmed_rep = ParmedMolRep(structure)
        return {
            "molecules": parmed_rep,
            "energy": self._get_out_rep("energy")(energy)
        }