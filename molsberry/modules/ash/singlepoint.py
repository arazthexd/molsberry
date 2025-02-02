from typing import Any, Dict, List, Optional, Callable
from abc import ABC, abstractmethod

import os
import glob
import subprocess
from copy import deepcopy

import ash

from molsberry.core.data import Representation

from ...core import LigandData, ProteinData, MoleculeData
from ...core import SimpleBlock, PLInteractionJob

from .representations import ASHFragmentRep
from .interface import ASHInterface

class ASHSinglePointCalculator(ASHInterface, SimpleBlock):
    name = "ashspc"
    display_name = "ASH Single Point Calculator"
    inputs = [
        ("fragments", MoleculeData, ASHFragmentRep, False)
    ]
    outputs = [
        ("energy", None, None, False)
    ]
    batch_groups = []
    # potential_ash_input_keys = ["molecules", "molecule", "ligands", "ligand",
    #                             "proteins", "protein", "pockets", "pocket",
    #                             "fragments"]
    
    def __init__(self,
                 theory: Any,
                 debug = False, 
                 save_output = False, 
                 num_workers = None):
        self.theory = deepcopy(theory)

        ASHInterface.__init__(self)
        SimpleBlock.__init__(self, debug=debug, 
                             save_output=save_output, 
                             num_workers=num_workers)
    
    def operate(self, input_dict: Dict[str, Any]) -> Dict[str, Representation]:
       
       fragrep: ASHFragmentRep = input_dict["fragments"]
       frag = fragrep.content
       ash_results = ash.Singlepoint(fragment=frag, theory=self.theory, 
                                     printlevel=0)
       out_dict = {"energy": ash_results.energy}
       return {k: self._get_out_rep(k)(v) for k, v in out_dict.items()}