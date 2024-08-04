from typing import Any, Dict, List

import os
import glob
from copy import deepcopy

from ...core.templates import LigandAnalyzerBlock
from ...core.templates.contexted import Contexted
from ...core.data.special_cls import Ligand, Protein
from .representations import MOPACInputMolRep
from .configs import MOPACConfig, MOPACMozymeConfig
from .utils import (
    ligand_to_mopacrep, generate_random_input_file, write_and_run_mopac
)
from ...global_conf import (
    make_sure_exists, TMP_DIR, MAIN_DIR, RANDOM_JOB_KEY_LENGTH
)

MOPAC_OUTPUT_FOLDERNAME = "mopac_out"
MOPAC_OUTPUT_DIR = os.path.join(MAIN_DIR, MOPAC_OUTPUT_FOLDERNAME)
make_sure_exists(MOPAC_OUTPUT_DIR)

MOPAC_TMP_FOLDERNAME = "mopac"
MOPAC_TMP_DIR = os.path.join(TMP_DIR, MOPAC_TMP_FOLDERNAME)
make_sure_exists(MOPAC_TMP_DIR)

MOPAC_OPTIMIZE_KEYWORDS = ["BFGS", "LBFGS", "DFP", "EF"]
    
class MOPACHFormationCalculator(LigandAnalyzerBlock):
    name = "MOPAC Formation Heat Calculator"
    output_keys = ["hformation"]
    output_types = [float]

    def __init__(self, config: MOPACConfig = MOPACConfig(), 
                 debug: bool = False, save_output: bool = False) -> None:
        super().__init__(debug=debug, save_output=save_output)
        self.config = deepcopy(config)
        self.config.keywords = [key for key in self.config.keywords 
                                if key not in MOPAC_OPTIMIZE_KEYWORDS]
        if "NOOPT" not in self.config.keywords:
            self.config.keywords.append("NOOPT")
        
        self.init_config = deepcopy(self.config)
    
    def reset_config(self):
        self.config = deepcopy(self.init_config)

    def analyze(self, ligand: Ligand) -> Dict[str, float]:
        mopac_rep = ligand_to_mopacrep(ligand)
        self.config.add_fragment(mopac_rep)
        path = generate_random_input_file(base_dir=MOPAC_TMP_DIR, 
                                           key_length=RANDOM_JOB_KEY_LENGTH)
        if self.debug:
            print("MOPAC input file written:", path)
            # TODO: I want this to be matched to ligand "name" which might need
            # separate implementation.
        
        write_and_run_mopac(path, self.config, MOPAC_OUTPUT_DIR, self.debug)

    # TODO: save_output implementation with new path system.
        
        