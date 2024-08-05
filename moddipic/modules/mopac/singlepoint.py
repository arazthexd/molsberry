from typing import Any, Dict, List

import os
import glob
from copy import deepcopy

from ...core.data.abstract import SpecialDataClass
from ...core.templates import LigandAnalyzerBlock
from ...core.templates import Contexted
from ...core.data import Ligand, Protein
from ...global_conf import RANDOM_JOB_KEY_LENGTH

from .representations import MOPACInputMolRep
from .configs import MOPACConfig, MOPACMozymeConfig
from .utils import (
    ligand_to_mopacrep, generate_random_input_file, write_and_run_mopac
)
from .shared import MOPAC_OUTPUT_DIR, MOPAC_TMP_DIR, MOPAC_OPTIMIZE_KEYWORDS
    
class MOPACSinglePointCalculator:
    name = "MOPAC Formation Heat Calculator"
    output_keys = ["hformation"]
    output_types = [float]

    def __init__(self, config: MOPACConfig = MOPACConfig(), 
                 debug: bool = False, save_output: bool = False) -> None:
        self.config = deepcopy(config)
        self.config.keywords = [key for key in self.config.keywords 
                                if key not in MOPAC_OPTIMIZE_KEYWORDS]
        if "NOOPT" not in self.config.keywords:
            self.config.keywords.append("NOOPT")
        
        self.init_config = deepcopy(self.config)
    
    def reset_config(self):
        self.config = deepcopy(self.init_config)

    def run_spc(self, mopac_rep: MOPACInputMolRep, 
                debug: bool = False) -> Dict[str, float]:
        self.reset_config()
        self.config.add_fragment(mopac_rep)
        path = generate_random_input_file(base_dir=MOPAC_TMP_DIR, 
                                          key_length=RANDOM_JOB_KEY_LENGTH)
        if debug:
            print("MOPAC input file written:", path)
            # TODO: I want this to be matched to ligand "name" which might need
            # separate implementation.
        
        write_and_run_mopac(path, self.config, debug)

        mopac_out_path = path[:-4] + ".out"
        with open(mopac_out_path, "r") as f:
            out_str = f.read()

        hf = float(out_str.split("FINAL HEAT OF FORMATION =")[1].split()[0])

        return {
            "full_out_path": mopac_out_path,
            "h_formation": hf
        }

    # TODO: save_output implementation with new path system.
        
class MOPACLigandSinglePointCalculator(MOPACSinglePointCalculator,
                                       LigandAnalyzerBlock):
    def __init__(self, config: MOPACConfig = MOPACConfig(), 
                 debug: bool = False, save_output: bool = False) -> None:
        LigandAnalyzerBlock.__init__(self, 
                                     debug=debug, save_output=save_output)
        MOPACSinglePointCalculator.__init__(self, config=config)
    
    def analyze(self, ligand: Ligand):
        mopac_rep = ligand.get_representation(MOPACInputMolRep)
        return self.run_spc(mopac_rep)

# TODO: Node for multiple molecule single point calculations