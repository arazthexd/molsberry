import os
import glob
import subprocess

from ...core.data.special_cls import Ligand, Protein
from ...utils.iotools import generate_random_str
from .representations import MOPACInputMolRep
from .configs import MOPACConfig

def ligand_to_mopacrep(ligand: Ligand) -> MOPACInputMolRep:
    mopac_mol_rep = ligand.get_representation(MOPACInputMolRep)
    return mopac_mol_rep

def protein_to_mopacrep(protein: Protein) -> MOPACInputMolRep:
    # TODO: implement
    pass

def write_and_run_mopac(path: str, mopac_config: MOPACConfig, 
                        debug: bool = False):
    cur_dir = os.curdir
    with open(path, "w") as f:
        f.write(mopac_config.get_config_str())
    os.chdir(os.path.dirname(path))
    subprocess.run(["mopac", os.path.basename(path)], capture_output=not debug)
    os.chdir(cur_dir)
    # TODO: What exactly to return and what exactly to move to output folder?

def generate_random_input_file(base_dir: str, key_length: int) -> str:
    while True:
        random_key = generate_random_str(key_length)
        matched_to_key = glob.glob(os.path.join(base_dir, f"*{random_key}*"))
        if len(matched_to_key) == 0:
            return os.path.join(base_dir, 
                                generate_random_str(key_length)+".mop")
    