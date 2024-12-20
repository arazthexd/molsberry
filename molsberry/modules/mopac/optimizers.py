from typing import Any, Dict, List

import os
import glob
from copy import deepcopy

import numpy as np

from molsberry.core.data import Representation

from ...core import SimpleBlock, PLJob
from ...core import LigandData, ProteinData, PDBPathRep
from ...modules.rdkit import RDKitMolRep
from ...global_conf import DATA_UNIQUE_CODE_LEN

from ...core.data.generic import NumericData, FloatRep

from .representations import MOPACInputMolRep
from .configs import MOPACConfig, MOPACMozymeConfig
from .interface import MOPACInterface
from .shared import (
    MOPAC_OPTIMIZE_KEYWORDS,
    MOPAC_SINGLEPOINT_KEYWORDS
)

class MOPACOptimizer(MOPACInterface, SimpleBlock):
    def __init__(self, config: MOPACConfig = MOPACConfig(),
                 opt_algorithm: str | None = None) -> None:
        config = deepcopy(config)
        config.keywords.append("THREADS=1")
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
        out_path, arc_path = self.run_job(self.config, debug=debug,
                                          base_dir=self.base_dir)
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

class MOPACLigandOptimizer(MOPACOptimizer, SimpleBlock):
    name = "mopacopt_lig"
    display_name = "MOPAC Ligand Optimizer"
    inputs = [
        ("ligands", LigandData, [MOPACInputMolRep, RDKitMolRep], False)
    ]
    outputs = [
        ("ligands", LigandData, [MOPACInputMolRep, 
                                 RDKitMolRep, PDBPathRep], False),
        ("e_init", NumericData, FloatRep, False),
        ("e_final", NumericData, FloatRep, False)
    ]
    batch_groups = []

    def __init__(self, config: MOPACConfig = MOPACConfig(), 
                 opt_algorithm: str | None = None,
                 debug: bool = False, save_output: bool = False,
                 num_workers: int = None) -> None:
        SimpleBlock.__init__(self, debug=debug, save_output=save_output,
                              num_workers=num_workers)
        MOPACOptimizer.__init__(self, config=config, 
                                opt_algorithm=opt_algorithm)
    
    def operate(self, input_dict: Dict[str, Representation]) \
        -> Dict[str, Representation]:
        mopacrep, rdrep = input_dict[self.input_keys[0]]
        mopacrep: MOPACInputMolRep
        rdrep: RDKitMolRep
        output = self.run_opt(mopacrep)
        mopacrep.update_coordinates(output["coordinates"])
        rdrep.update_coordinates(output["coordinates"])
        pdbrep = PDBPathRep(output["out_path"][:-4]+".pdb")
        return {
            "ligands": [mopacrep, rdrep, pdbrep],
            "e_init": output["pre_energy"],
            "e_final": output["post_energy"]
        }
