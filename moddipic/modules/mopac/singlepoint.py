from typing import Any, Dict, List

import os
import glob
import subprocess
from copy import deepcopy

from ...core.data.abstract import SpecialDataClass
from ...core.templates import LigandAnalyzerBlock
from ...core.templates import Contexted
from ...core.data import Ligand, Protein
from ...global_conf import RANDOM_JOB_KEY_LEN

from .representations import MOPACInputMolRep
from .configs import MOPACConfig, MOPACMozymeConfig
from .interface import MOPACInterface
from .shared import MOPAC_OUTPUT_DIR, MOPAC_TMP_DIR, MOPAC_OPTIMIZE_KEYWORDS

class MOPACSinglePointCalculator(MOPACInterface):
    def __init__(self, config: MOPACConfig = MOPACConfig()) -> None:
        config = deepcopy(config)
        config.keywords = [key for key in config.keywords 
                           if key not in MOPAC_OPTIMIZE_KEYWORDS]
        if "NOOPT" not in config.keywords:
            config.keywords.append("NOOPT")

        super().__init__(config=config)

    def run_spc(self, mopac_rep: MOPACInputMolRep | List[MOPACInputMolRep], 
                debug: bool = False) -> Dict[str, Any]:
        self.reset_config()
        self.update_config(mopac_rep=mopac_rep)
        out_path, arc_path = self.run_job(self.config, debug)
        with open(out_path, "r") as f:
            out_str = f.read()

        hf = float(out_str.split("FINAL HEAT OF FORMATION =")[1].split()[0])

        return {
            "out_path": out_path,
            "arc_path": arc_path,
            "energy": hf
        }

    # TODO: save_output implementation with new path system.
        
class MOPACLigandSinglePointCalculator(MOPACSinglePointCalculator,
                                       LigandAnalyzerBlock):
    name = "MOPAC Ligand Single Point Calculator"
    output_keys = ["energy", "out_path", "arc_path"]
    output_types = [float, str]

    def __init__(self, config: MOPACConfig = MOPACConfig(), 
                 debug: bool = False, save_output: bool = False) -> None:
        LigandAnalyzerBlock.__init__(self, 
                                     debug=debug, save_output=save_output)
        MOPACSinglePointCalculator.__init__(self, config=config)
    
    def analyze(self, ligand: Ligand) -> Dict[str, Any]:
        mopac_rep = ligand.get_representation(MOPACInputMolRep)
        return self.run_spc(mopac_rep)

# TODO: Node for multiple molecule single point calculations