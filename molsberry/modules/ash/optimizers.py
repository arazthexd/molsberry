from typing import Any, Dict, List

import os
import glob
from copy import deepcopy

import numpy as np

import ash
from molsberry.core.data import Representation

from ...core import suppress_prints
from ...core import SimpleBlock, PLJob
from ...core import MoleculeData, ProteinData, PDBPathRep
from ...modules.rdkit import RDKitMolRep
from ...global_conf import DATA_UNIQUE_CODE_LEN

from .representations import ASHFragmentRep
from .interface import ASHInterface

import ash

class ASHFragmentOptimizer(ASHInterface, SimpleBlock):
    name = "ashopt"
    display_name = "ASH Fragment Optimizer"
    inputs = [
        ("fragments", MoleculeData, ASHFragmentRep, False)
    ]
    outputs = [
        ("fragments", MoleculeData, [ASHFragmentRep, RDKitMolRep], False)
    ]
    batch_groups = []

    def __init__(self,
                 theory: Any,
                 maxiter: int = 250,
                 debug = False, 
                 save_output = False, 
                 num_workers = None):
        self.theory = deepcopy(theory)
        self.maxiter = maxiter

        ASHInterface.__init__(self)
        SimpleBlock.__init__(self, debug=debug, 
                             save_output=save_output, 
                             num_workers=num_workers)
    
    def operate(self, input_dict: Dict[str, Any]) -> Dict[str, Representation]:
        fragrep: ASHFragmentRep = input_dict[self.input_keys[0]]
        smi = fragrep.smiles
        frag = deepcopy(fragrep.content)
        with suppress_prints():
            results = ash.Optimizer(theory=self.theory, 
                                    fragment=frag,
                                    convergence_setting="SuperLoose", # temp
                                    maxiter=self.maxiter,
                                    printlevel=0)
        
        frag_key = self.output_keys[0]
        out_reps = self._get_out_rep(frag_key)
        fragrep: ASHFragmentRep = out_reps[0](frag)
        fragrep.smiles = smi
        rdrep = fragrep.to_RDKitMolRep()
        return {frag_key: [fragrep, rdrep]}

