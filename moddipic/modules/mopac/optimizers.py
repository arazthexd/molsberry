from typing import Any, Dict, List

import os
import glob
from copy import deepcopy

import numpy as np

from ...core.templates import LigandConverterBlock
from ...core.templates import Contexted
from ...core.data import Ligand, Protein
from ...global_conf import RANDOM_JOB_KEY_LEN

from .representations import MOPACInputMolRep
from .configs import MOPACConfig, MOPACMozymeConfig
from .interface import MOPACInterface
from .shared import (
    MOPAC_OUTPUT_DIR, 
    MOPAC_TMP_DIR, 
    MOPAC_OPTIMIZE_KEYWORDS,
    MOPAC_SINGLEPOINT_KEYWORDS
)

class MOPACOptimizer(MOPACInterface):
    def __init__(self, config: MOPACConfig = MOPACConfig(),
                 opt_algorithm: str | None = None) -> None:
        config = deepcopy(config)
        config.keywords = [key for key in config.keywords
                           if key not in MOPAC_SINGLEPOINT_KEYWORDS]
        n_opt_keys = len(set(
            MOPAC_OPTIMIZE_KEYWORDS).intersection(config.keywords))
        if n_opt_keys == 0:
            if opt_algorithm is None:
                raise ValueError(
                    "Provide opt_algorithm or include keyword in input config")
            if opt_algorithm not in MOPAC_OPTIMIZE_KEYWORDS:
                raise ValueError("Opt Algorithm not known/implemented...")
            config.keywords.append(opt_algorithm)
        if n_opt_keys == 1:
            pass
        if n_opt_keys > 1:
            raise ValueError("Only one keyword for opt should be in config.")
        
        super().__init__(config=config)
    
    def run_opt(self, mopac_rep: MOPACInputMolRep | List[MOPACInputMolRep], 
                debug: bool = False) -> Dict[str, Any]:
        self.reset_config()
        self.update_config(mopac_rep=mopac_rep)
        out_path, arc_path = self.run_job(self.config, debug)
        with open(out_path, "r") as f:
            out_str = f.read()

        print(out_path)
        print(out_str)
        pre_energy = float(out_str.split("CYCLE:     1")[1].split()[8])
        post_energy = float(
            out_str.split("FINAL HEAT OF FORMATION =")[1].split()[0])
        coordinates = np.array([
            [float(coord) for coord in line.split()[2:]] 
            for line in out_str.split(" CARTESIAN COORDINATES\n\n")[1] \
                .split("\n\n           Empirical Formula")[0].splitlines()])

        return {
            "out_path": out_path,
            "pre_energy": pre_energy,
            "post_energy": post_energy,
            "coordinates": coordinates
        }

class MOPACLigandOptimizer(MOPACOptimizer, Contexted, LigandConverterBlock):
    name = "MOPAC Ligand Optimizer"
    output_keys = LigandConverterBlock.output_keys + \
        ["pre_energy", "post_energy"]
    output_types = [LigandConverterBlock.single_data_type] + \
        [float, str]
    output_context_keys = ["pre_energy", "post_energy"]
    output_context_types = [float, float]

    def __init__(self, config: MOPACConfig = MOPACConfig(), 
                 opt_algorithm: str | None = None,
                 debug: bool = False, save_output: bool = False) -> None:
        LigandConverterBlock.__init__(self, 
                                     debug=debug, save_output=save_output)
        MOPACOptimizer.__init__(self, config=config, 
                                opt_algorithm=opt_algorithm)
    
    def convert(self, ligand: Ligand) -> Ligand:
        mopac_rep = ligand.get_representation(MOPACInputMolRep)
        output = self.run_opt(mopac_rep)
        newlig = ligand.return_with_new_coords(output["coordinates"])
        self.output_context = {
            "pre_energy": output["pre_energy"],
            "post_energy": output["post_energy"]
        }
        return newlig
        #     "ligands": newlig, # ligands
        #     "pre_energy": output["pre_energy"],
        #     "post_energy": output["post_energy"]
        # } # TODO: Make it have output context