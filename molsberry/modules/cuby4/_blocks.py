from copy import deepcopy, copy
from typing import Dict, List, Any, Type, Tuple

from rdkit import Chem
from openmm import unit
import parmed

from ...core import *
from ..mopac import MOPACInputMolRep
from ..openmm import OpenMMInputMolRep
from ..rdkit import RDKitMolRep
from ..parmed import ParmedMolRep

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

#####################################################################
##                      General (Base) Blocks                      ##
#####################################################################

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
    def generate_out_dict(self, 
                          input_dict,
                          cuby4_out: str, 
                          final_config) -> Dict[str, Representation]:
        energy = float(cuby4_out.split("Energy:")[-1].split()[0])
        return {"energy": self._get_out_rep("energy")(energy)}

#####################################################################
##                          MOPAC Blocks                           ##
#####################################################################
    
class Cuby4MOPACEnergyCalculator(Cuby4GeneralEnergyCalculator):
    name = "c4mopacspc"
    display_name = "Cuby4-MOPAC Single Point Calculator"
    inputs = [
        ("molecules", MoleculeData, MOPACInputMolRep, False)
    ]
    outputs = [
        ("energy", NumericData, FloatRep, False),
    ]
    batch_groups = []

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

    def generate_job_config(self, 
                            input_dict: Dict[str, ParmedMolRep]) -> C4MIC:
        parmed_rep: ParmedMolRep = input_dict[self.input_keys[0]]
        rdkit_rep: RDKitMolRep = parmed_rep.to_RDKitMolRep()
        mopac_rep = MOPACInputMolRep.from_RDKitMolRep(rdkit_rep) # TODO: Check

        jc = Cuby4EnergyJobConfig(
            geometry=generate_path_in_dir(6, self.work_dir, "_c4mecinp.pdb"),
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
                          input_dict,
                          cuby4_out: str, 
                          final_config) -> Dict[str, Representation]:
        energy = float(cuby4_out.split("Energy:")[-1].split()[0])
        return {"energy": self._get_out_rep("energy")(energy)}

class Cuby4MOPACEnergyOptimizer(Cuby4GeneralBlock):
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

        Cuby4Interface.__init__(self, 
                                interface_config=interface_config,
                                cuby4_exe=cuby4_exe,
                                work_dir=work_dir)
        SimpleBlock.__init__(self, debug=debug, save_output=save_output)

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
                                              self.work_dir, 
                                              "_c4meoinp.pdb"),
                charge=mopac_rep.charge,
                mult=1, # TODO
                restart_file=generate_path_in_dir(6, 
                                                  self.work_dir, 
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

#####################################################################
##                          AMBER Blocks                           ##
#####################################################################

class _AMBER_Utils:
    @staticmethod
    def omm_to_parms(omm_rep: OpenMMInputMolRep, base_dir: str) \
        -> Tuple[str, str]:
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

class Cuby4AMBEREnergyCalculator(Cuby4GeneralEnergyCalculator):
    name = "c4amberspc"
    display_name = "Cuby4-AMBER Single Point Calculator"
    inputs = [
        ("molecules", MoleculeData, MOPACInputMolRep, False)
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
        _, geometry, prmtop = _AMBER_Utils.omm_to_parms(omm_rep, self.base_dir)

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
    
    mol_rep = OpenMMInputMolRep

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
        _, geometry, prmtop = _AMBER_Utils.omm_to_parms(omm_rep, self.base_dir)

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
    







class Cuby4GeneralEnergyCalculator(Cuby4GeneralBlock):
    name = "c4generalspc"
    display_name = "Cuby4-General Single Point Calculator"
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
    
    def generate_out_dict(self, cuby4_out: str):
        energy = float(cuby4_out.split("Energy:")[-1].split()[0])
        return {"energy": self._get_out_rep("energy")(energy)}

class Cuby4GeneralEnergyOptimizer(Cuby4GeneralBlock):
    name = "c4generalopt"
    display_name = "Cuby4-General Energy Optimizer"
    inputs = [
        ("molecules", MoleculeData, None, False)
    ]
    outputs = [
        ("molecules", MoleculeData, None, False),
        ("energy", NumericData, FloatRep, False),
    ]
    batch_groups = []
    
    @property
    def inputs(self):
        return [
            ("molecules", MoleculeData, self.mol_rep, False)
        ]
    
    @property
    def outputs(self):
        return [
            ("molecules", MoleculeData, self.mol_rep, False),
            ("energy", NumericData, FloatRep, False)
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
                restart_file=generate_path_in_dir(6, self.work_dir, ".pdb"),
                maxcycles=7
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


class _AMBER_Utils:
    @staticmethod
    def omm_to_parms(omm_rep: OpenMMInputMolRep, base_dir: str) \
        -> Tuple[str, str]:
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

class Cuby4AMBEREnergyCalculator(Cuby4GeneralEnergyCalculator):
    name = "c4amberspc"
    display_name = "Cuby4-AMBER Single Point Calculator"
    
    mol_rep = OpenMMInputMolRep

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
        _, geometry, prmtop = _AMBER_Utils.omm_to_parms(omm_rep, self.base_dir)

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
    
    mol_rep = OpenMMInputMolRep

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
        _, geometry, prmtop = _AMBER_Utils.omm_to_parms(omm_rep, self.base_dir)

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

class Cuby4QMMMEnergyOptimizer(Cuby4Interface, SimpleBlock):
    name = "c4qmmmspcop"
    display_name = "Cuby4-QMMM Optimizer"
    inputs = [
        ("qm_region", MoleculeData, ParmedMolRep, False),
        ("nonqm_region", MoleculeData, ParmedMolRep, False)
    ]
    outputs = [
        ("energy", NumericData, FloatRep, False),
        ("qm_region", MoleculeData, OpenMMInputMolRep, False),
        ("nonqm_region", MoleculeData, OpenMMInputMolRep, False)
    ]
    batch_groups = [("qm_region", "nonqm_region")]

    # TODO : we want metals in qmregion

    def __init__(self, 
                 interface_config,
                 debug: bool = False, 
                 save_output: bool = False,
                 work_dir: str = ".",
                 cuby4_exe: str = "auto",
                 threads: int = 1,
                 maxcycles: int = 10,
                 restart_file: bool = True) -> None:
        
        self.work_dir = work_dir
        self.threads = threads
        self.maxcycles = maxcycles
        self.restart_file = generate_path_in_dir(6, self.work_dir, ".pdb")
        interface_config = deepcopy(interface_config)

        Cuby4Interface.__init__(self, 
                                interface_config=interface_config,
                                cuby4_exe=cuby4_exe,
                                work_dir=work_dir)
        SimpleBlock.__init__(self, debug=debug, save_output=save_output)

    def operate(self, input_dict: Dict[str, ParmedMolRep]) -> Dict[str, Representation]:
        config = deepcopy(self.interface_config)

        omm_rep_qm = input_dict["qm_region"]
        omm_rep_nonqm = input_dict["nonqm_region"]
        qmrep = copy(omm_rep_qm)
        nonqmrep = copy(omm_rep_nonqm) # TODO: make it deep!

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
        config.config["cuby_threads"] = self.threads
        config.config['restart_file'] = self.restart_file

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
        config.config["job"] = "optimize"
        config.config["maxcycles"] = self.maxcycles

        output: str = self.run(config)
        energy = float(output.split("Energy:")[-1].split()[0])

        qmrep, nonqmrep = self.update_pl_coords(qmrep, 
                                                nonqmrep,
                                                config.config["restart_file"])
        
        # qmmol = OpenMMInputMolRep.to_RDKitMolRep(qmrep)
        # nonqmmol = OpenMMInputMolRep.to_RDKitMolRep(nonqmrep)
        
        return {
            "energy": self._get_out_rep("energy")(energy),
            "qm_region": qmrep,
            "nonqm_region": nonqmrep
        }
    
    def update_pl_coords(self, 
                         qmrep: OpenMMInputMolRep,
                         nonqmrep: OpenMMInputMolRep,
                         pdb_path: str) -> Tuple[OpenMMInputMolRep, 
                                                 OpenMMInputMolRep]:
        
        pdb_struct: parmed.Structure = parmed.load_file(pdb_path)
        pdb_coords = pdb_struct.coordinates

        qmcoords = pdb_coords[:qmrep.topology.getNumAtoms()] * unit.angstrom
        nonqmcoords = pdb_coords[qmrep.topology.getNumAtoms():] * unit.angstrom

        qmrep.positions = qmcoords
        nonqmrep.positions = nonqmcoords

        return qmrep, nonqmrep