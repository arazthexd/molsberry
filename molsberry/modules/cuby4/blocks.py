from copy import deepcopy
from typing import Dict, List

from ...core import *
from ..mopac import MOPACInputMolRep

from .configs import (
    Cuby4MOPACInterfaceConfig as C4MIC, 
    Cuby4EnergyJobConfig, 
    Cuby4MergedConfig,
    MOPACMozymeConfig,
    Cuby4OptimizeJobConfig
)
from .interface import Cuby4Interface

class Cuby4MOPACEnergyCalculator(Cuby4Interface, SimpleBlock):
    name = "c4mopacspc"
    display_name = "Cuby4-MOPAC Single Point Calculator"
    inputs = [
        ("molecules", MoleculeData, MOPACInputMolRep, False)
    ]
    outputs = [
        ("energy", NumericData, FloatRep, False),
    ]
    batch_groups = []
    potential_cuby_input_keys = ["molecules", "molecule", "ligands", "ligand",
                                  "proteins", "protein", "pockets", "pocket"]

    def __init__(self, 
                 interface_config: C4MIC = C4MIC(),
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

        self.interface_config: C4MIC
    
    def operate(self, input_dict: Dict[str, MOPACInputMolRep]) \
        -> Dict[str, Representation]:
        cuby_input_keys = [k for k in self.input_keys
                            if k in self.potential_cuby_input_keys]
        assert len(cuby_input_keys) > 0

        for k in cuby_input_keys:
            frag : MOPACInputMolRep = input_dict[k]
            jc = Cuby4EnergyJobConfig(
                geometry=generate_path_in_dir(6, self.work_dir, ".pdb"),
                charge=frag.charge,
                mult=1 # TODO
            )
            if self.interface_config.mozyme:
                jc = Cuby4MergedConfig.from_config_list([
                    jc, MOPACMozymeConfig(
                        setpi=frag.setpi[:50],
                        neg_cvb=frag.neg_cvb
                    )
                ])

                if self.interface_config.setpi == False:
                    _ = jc.config.pop("mopac_setpi")

                print(jc.config["mopac_keywords"])
                if self.interface_config.cvb == False:
                    kwds = jc.config["mopac_keywords"].split()
                    kwds = [kwd for kwd in kwds if "CVB" not in kwd]
                    jc.config["mopac_keywords"] = " ".join(kwds)

            with open(jc.config["geometry"], "w") as f:
                f.write(frag.coordinates)
            
            full_config = Cuby4MergedConfig.from_config_list([
                self.interface_config, jc
            ])

            if self.interface_config.setcharges:
                full_config.config["mopac_setcharge"] = {
                    i+1: c for i, c in enumerate(frag.atom_charges)}

            output: str = self.run(full_config)
            energy = float(output.split("Energy:")[-1].split()[0])
            return {"energy": self._get_out_rep("energy")(energy)}
        


class Cuby4MOPACEnergyOptimizer(Cuby4Interface, SimpleBlock):
    name = "c4mopacopt"
    display_name = "Cuby4-MOPAC Energy Optimizer"
    inputs = [
        ("molecules", MoleculeData, MOPACInputMolRep, False)
    ]
    outputs = [
        ("molecules", MoleculeData, PDBPathRep, False),
        ("energy", NumericData, FloatRep, False)
    ]
    batch_groups = []
    potential_cuby_input_keys = ["molecules", "molecule", "ligands", "ligand",
                                  "proteins", "protein", "pockets", "pocket"]

    def __init__(self, 
                 interface_config: C4MIC = C4MIC(),
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

        self.interface_config: C4MIC
    
    def operate(self, input_dict: Dict[str, MOPACInputMolRep]) \
        -> Dict[str, Representation]:
        cuby_input_keys = [k for k in self.input_keys
                            if k in self.potential_cuby_input_keys]
        assert len(cuby_input_keys) > 0

        for k in cuby_input_keys:
            frag : MOPACInputMolRep = input_dict[k]
            jc = Cuby4OptimizeJobConfig(
                geometry=generate_path_in_dir(6, self.work_dir, ".pdb"),
                charge=frag.charge,
                mult=1, # TODO
                restart_file=generate_path_in_dir(6, self.work_dir, ".pdb")
            )
            if self.interface_config.mozyme:
                jc = Cuby4MergedConfig.from_config_list([
                    jc, MOPACMozymeConfig(
                        setpi=frag.setpi,
                        neg_cvb=frag.neg_cvb
                    )
                ])

                if self.interface_config.setpi == False:
                    _ = jc.config.pop("mopac_setpi")

            with open(jc.config["geometry"], "w") as f:
                f.write(frag.coordinates)
            
            full_config = Cuby4MergedConfig.from_config_list([
                self.interface_config, jc
            ])

            output: str = self.run(full_config)
            energy = float(output.split("Energy:")[-1].split()[0])
            return {
                "molecules": self._get_out_rep("molecules")(jc.config["restart_file"]),
                "energy": self._get_out_rep("energy")(energy)
            }