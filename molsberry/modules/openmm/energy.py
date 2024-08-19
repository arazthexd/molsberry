from typing import Any, Dict, List, Optional, Callable
from abc import ABC, abstractmethod
from collections.abc import Iterable

import os
import glob
import subprocess
from copy import deepcopy

from openmm.app import ForceField

from ...core import LigandData, ProteinData, MoleculeData
from ...core import SimpleBlock, PLInteractionJob, Representation

from .representations import OpenMMInputMolRep
from .constants import DEFAULT_FORCEFIELDS, DEFAULT_PH
from .interface import OpenMMInterface

class OpenMMEnergyCalculator(OpenMMInterface, SimpleBlock):
    name = "openmm_energy"
    display_name = "OpenMM Energy Calculator Block"
    inputs = [
        ("molecules", MoleculeData, OpenMMInputMolRep, False)
    ]
    outputs = [
        ("energy", None, None, False),
    ]
    batch_groups = []
    potential_openmm_input_keys = ["molecules", "molecule", 
                                    "ligands", "ligand", 
                                    "proteins", "protein", 
                                    "pockets", "pocket"]

    def __init__(self, forcefield: ForceField | None = None, 
                 debug: bool = False, save_output: bool = False) -> None:
        SimpleBlock.__init__(self, debug=debug, save_output=save_output)
        self.forcefield = forcefield

    def operate(self, input_dict: Dict[str, OpenMMInputMolRep]) \
        -> Dict[str, Representation]:
        openmm_input_keys = [k for k in self.input_keys
                            if k in self.potential_openmm_input_keys]
        assert len(openmm_input_keys) > 0

        result_dict = {"energy": self.calc_energy([input_dict[k] 
                                                   for k in openmm_input_keys])}
        return {k: self._get_out_rep(k)(v) for k, v in result_dict.items()}

    def calc_energy(self, omm_rep: OpenMMInputMolRep | List[OpenMMInputMolRep],
                    debug: bool = False) -> Dict[str, Any]:

        if isinstance(omm_rep, Iterable):
            omm_rep = OpenMMInputMolRep.merge_reps(omm_rep, self.forcefield)

        elif self.forcefield is not None:
            omm_rep = OpenMMInputMolRep.merge_reps([omm_rep], self.forcefield)
        
        return self.rep_to_energy(omm_rep)
        