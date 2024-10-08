from typing import List

import os, pathlib
import glob
import subprocess
from copy import deepcopy

from ...core import generate_random_str
from ...global_conf import DATA_UNIQUE_CODE_LEN

from .representations import MOPACInputMolRep
from .configs import MOPACConfig

class MOPACInterface:
    def __init__(self, config: MOPACConfig = MOPACConfig()) -> None:
        self.config = deepcopy(config)
        self.init_config = deepcopy(config)

    def reset_config(self):
        self.config = deepcopy(self.init_config)
    
    def update_config(self, 
                      mopac_rep: MOPACInputMolRep | List[MOPACInputMolRep]):
        if isinstance(mopac_rep, MOPACInputMolRep):
            self.config.add_fragment(mopac_rep)
        elif isinstance(mopac_rep, list):
            [self.config.add_fragment(rep) for rep in mopac_rep]
    
    @staticmethod
    def run_job(config: MOPACConfig, base_dir: str, debug: bool):
        path = MOPACInterface.generate_random_input_file(
            base_dir=base_dir, 
            key_length=DATA_UNIQUE_CODE_LEN
        )

        if debug:
            print("MOPAC input file written:", path)
            # TODO: I want this to be matched to ligand "name" which might need
            # separate implementation.

        MOPACInterface.write_and_run_mopac(path, config, debug)

        out_path = path[:-4] + ".out"
        arc_path = path[:-4] + ".arc"
        
        return out_path, arc_path

    @staticmethod
    def write_and_run_mopac(path: str, mopac_config: MOPACConfig, 
                            debug: bool = False):
        cur_dir = pathlib.Path(os.curdir).absolute()
        with open(path, "w") as f:
            f.write(mopac_config.get_config_str())
        os.chdir(os.path.dirname(path))
        subprocess.run(["mopac", os.path.basename(path)], 
                       capture_output=not debug)
        os.chdir(cur_dir)
        # TODO: What exactly to return and what exactly to move to out folder?
    
    @staticmethod
    def generate_random_input_file(base_dir: str, key_length: int) -> str:
        while True:
            random_key = generate_random_str(key_length)
            matched_to_key = glob.glob(os.path.join(base_dir, 
                                                    f"*{random_key}*"))
            if len(matched_to_key) == 0:
                return os.path.join(base_dir, 
                                    generate_random_str(key_length)+".mop")