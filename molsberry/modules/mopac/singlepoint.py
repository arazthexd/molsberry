from typing import Any, Dict, List, Optional, Callable
from abc import ABC, abstractmethod

import os
import glob
import subprocess
from copy import deepcopy

from molsberry.core.data import Representation

from ...core import LigandData, ProteinData, MoleculeData, NumericData, FloatRep, StringData, StringRep
from ...core import SimpleBlock, PLInteractionJob

from .representations import MOPACInputMolRep
from .configs import MOPACConfig, MOPACMozymeConfig
from .interface import MOPACInterface
from .shared import MOPAC_OPTIMIZE_KEYWORDS

class MOPACSinglePointCalculator(MOPACInterface, SimpleBlock):
    name = "mopacspc"
    display_name = "MOPAC Single Point Calculator"
    inputs = [
        ("molecules", MoleculeData, MOPACInputMolRep, False)
    ]
    outputs = [
        ("energy", NumericData, FloatRep, False),
        ("out_path", StringData, StringRep, False),
        ("arc_path", StringData, StringRep, False)
    ]
    batch_groups = []
    potential_mopac_input_keys = ["molecules", "molecule", "ligands", "ligand",
                                  "proteins", "protein", "pockets", "pocket"]

    def __init__(self, config: MOPACConfig = MOPACConfig(),
                 debug: bool = False, save_output: bool = False,
                 num_workers: int = None) -> None:
        config = deepcopy(config)
        config.keywords = [key for key in config.keywords 
                           if key not in MOPAC_OPTIMIZE_KEYWORDS]
        if "NOOPT" not in config.keywords:
            config.keywords.append("NOOPT")

        MOPACInterface.__init__(self, config=config)
        SimpleBlock.__init__(self, debug=debug, save_output=save_output,
                             num_workers=num_workers)

    def operate(self, input_dict: Dict[str, MOPACInputMolRep]) \
        -> Dict[str, Representation]:
        mopac_input_keys = [k for k in self.input_keys
                            if k in self.potential_mopac_input_keys]
        assert len(mopac_input_keys) > 0

        result_dict = self.run_spc([input_dict[k] for k in mopac_input_keys])
        return {k: self._get_out_rep(k)(v) for k, v in result_dict.items()}

    def calc_energy(self, mopac_rep: MOPACInputMolRep | List[MOPACInputMolRep]):
        out = self.run_spc(mopac_rep)
        return out["energy"]

    def run_spc(self, mopac_rep: MOPACInputMolRep | List[MOPACInputMolRep], 
                debug: bool = False) -> Dict[str, Any]:
        self.reset_config()
        self.update_config(mopac_rep=mopac_rep)
        out_path, arc_path = self.run_job(config=self.config, debug=debug,
                                          base_dir=self.base_dir)
        with open(out_path, "r") as f:
            out_str = f.read()

        hf = float(out_str.split("FINAL HEAT OF FORMATION =")[1].split()[0])

        return {
            "out_path": out_path,
            "arc_path": arc_path,
            "energy": hf
        }
        
class MOPACLigandSinglePointCalculator(MOPACSinglePointCalculator):
    name = "mopacspc_lig"
    display_name = "MOPAC Ligand Single Point Calculator"
    inputs = [
        ("ligands", LigandData, MOPACInputMolRep, False)
    ]
    mopac_inputs = ["ligands"]

class MOPACProteinSinglePointCalculator(MOPACSinglePointCalculator):
    name = "mopacspc_prot"
    display_name = "MOPAC Protein Single Point Calculator"
    inputs = [
        ("proteins", ProteinData, MOPACInputMolRep, False)
    ]
    mopac_inputs = ["proteins"]

class MOPACPLInteractionCalculator(PLInteractionJob, 
                                   MOPACSinglePointCalculator):
    name = "mopac_plinter"
    display_name = "MOPAC Protein Ligand Interaction Calculator"
    lig_rep = MOPACInputMolRep
    prot_rep = MOPACInputMolRep

    def __init__(self, config = MOPACConfig(),
                 debug: bool = False, save_output: bool = False, 
                 num_workers: int = None) -> None:
        MOPACSinglePointCalculator.__init__(self,
                                            debug=debug, 
                                            save_output=save_output, 
                                            config=config,
                                            num_workers=num_workers)
        self.calculator = MOPACProteinSinglePointCalculator(config=config)
        energy_fn = self.calculator.calc_energy
        PLInteractionJob.__init__(self, energy_fn=energy_fn)

    def combine_mols(self, ligand: MOPACInputMolRep, protein: MOPACInputMolRep):
        return [ligand, protein]
    
    def run_spc(self, reps: List[MOPACInputMolRep]) -> Dict[str, float]:
        assert len(reps) == 2
        self.calculator._parent = self._parent
        interaction_dict = self.calc_interaction(reps[0], reps[1])
        return interaction_dict

# TODO: Node for multiple molecule single point calculations