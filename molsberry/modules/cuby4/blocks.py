from copy import deepcopy
from typing import Dict, List, Any

from ...core import *
from ..mopac import MOPACInputMolRep
from ..openmm import OpenMMInputMolRep
from ..rdkit import RDKitMolRep

from .configs import (
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
from .interface import Cuby4Interface

class Cuby4GeneralEnergyCalculator(Cuby4Interface, SimpleBlock):
    name = "c4generalspc"
    display_name = "Cuby4-General Single Point Calculator"
    inputs = [
        ("molecules", MoleculeData, None, False)
    ]
    outputs = [
        ("energy", NumericData, FloatRep, False),
    ]
    batch_groups = []

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

    def generate_job_config(input_dict: Dict[str, Representation]) \
        -> Cuby4JobConfig:
        raise NotImplementedError()
    
    def operate(self, input_dict: Dict[str, Representation]) \
        -> Dict[str, Representation]:
            
            jc = self.generate_job_config(input_dict)
            
            full_config = Cuby4MergedConfig.from_config_list([
                self.interface_config, jc
            ])

            output: str = self.run(full_config)
            energy = float(output.split("Energy:")[-1].split()[0])
            return {"energy": self._get_out_rep("energy")(energy)}

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
    
    # TODO: Make it a subclass of General Energy Calculator Block

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

import parmed
import os, shutil

class Cuby4AMBEREnergyCalculator(Cuby4GeneralEnergyCalculator):
    name = "c4amberspc"
    display_name = "Cuby4-AMBER Single Point Calculator"
    inputs = [
        ("molecules", MoleculeData, OpenMMInputMolRep, False),
    ]
    outputs = [
        ("energy", NumericData, FloatRep, False),
    ]
    batch_groups = []

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
        omm_rep: OpenMMInputMolRep = input_dict[self.input_keys[0]]
        struct: parmed.Structure = parmed.openmm.load_topology(
            omm_rep.topology,
            omm_rep.system,
            omm_rep.positions
        )
        geometry = generate_path_in_dir(6, self.base_dir, ".pdb")
        struct.save(geometry)

        prmtop = generate_path_in_dir(6, self.base_dir, ".parm7")
        struct.save(prmtop)

        config = Cuby4EnergyJobConfig(
            geometry=geometry
        )
        
        mol_specifics = {
            "amber_top_file": prmtop
        }

        config.config.update(mol_specifics)
        
        return config

class Cuby4QMMMEnergyCalculator(Cuby4Interface, SimpleBlock):
    name = "c4qmmmspc"
    display_name = "Cuby4-QMMM Single Point Calculator"
    inputs = [
        ("qm_region", MoleculeData, OpenMMInputMolRep, False),
        ("nonqm_region", MoleculeData, OpenMMInputMolRep, False)
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

        Cuby4Interface.__init__(self, 
                                interface_config=interface_config,
                                cuby4_exe=cuby4_exe,
                                work_dir=work_dir)
        SimpleBlock.__init__(self, debug=debug, save_output=save_output)

    def operate(self, input_dict: Dict[str, OpenMMInputMolRep]) -> Dict[str, Representation]:
        config = deepcopy(self.interface_config)

        omm_rep_qm = input_dict["qm_region"]
        omm_rep_nonqm = input_dict["nonqm_region"]
        struct_qm: parmed.Structure = parmed.openmm.load_topology(
            omm_rep_qm.topology,
            omm_rep_qm.system,
            omm_rep_qm.positions
        )
        struct_nonqm: parmed.Structure = parmed.openmm.load_topology(
            omm_rep_nonqm.topology,
            omm_rep_nonqm.system,
            omm_rep_nonqm.positions
        )
        struct_whole: parmed.Structure = struct_qm + struct_nonqm

        # geometry, core
        geometry = generate_path_in_dir(6, self.base_dir, ".pdb")
        struct_whole.save(geometry)

        core = f"1-{len(struct_qm.atoms)}"

        config.config["geometry"] = geometry
        config.config["qmmm_core"] = core

        # calculation_qm: charge
        charge = round(sum([atom.charge for atom in struct_qm.atoms]))
        config.config["calculation_qm"]["charge"] = charge

        # calculation_mm: parm7 mereged
        parm7_whole = generate_path_in_dir(6, self.base_dir, ".parm7")
        struct_whole.save(parm7_whole)

        config.config["calculation_mm"]["amber_top_file"] = parm7_whole

        # calculation_qmregion_mm: parm7 lig
        parm7_qm = generate_path_in_dir(6, self.base_dir, ".parm7")
        struct_qm.save(parm7_qm)

        config.config["calculation_qmregion_mm"] = {"amber_top_file": parm7_qm}

        # job
        config.config["job"] = "energy"
        #config.config["maxcycles"] = 15

        output: str = self.run(config)
        energy = float(output.split("Energy:")[-1].split()[0])
        return {
                "energy": self._get_out_rep("energy")(energy)
            }
