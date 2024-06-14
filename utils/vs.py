from typing import List
from copy import copy
import os

import subprocess

from rdkit import Chem

from .common import GenericPocket

GNINA_CONFIG = {
    "bin_path": "gnina",
    "work_dir": "tmp",
    "other_args": []
}

def perform_docking(mols: List[Chem.Mol], pocket: GenericPocket, docking_engine="gnina", config: dict = dict()):

    if docking_engine == "gnina":
        return _perform_docking_gnina(mols, pocket, config)
    else:
        raise NotImplementedError()

def _perform_docking_gnina(mols: List[Chem.Mol], pocket: GenericPocket, config: dict = dict()):

    # Set up the config...
    c = copy(GNINA_CONFIG)
    c.update(config)
    if isinstance(c["other_args"], str):
        c["other_args"] = c["other_args"].split()
    if not os.path.exists(c["work_dir"]):
        os.mkdir(c["work_dir"])

    # Write the ligands that are gonna be docked to the work directory
    ligs_to_dock_path = os.path.join(c["work_dir"], "ligs_to_dock.sdf")
    writer = Chem.SDWriter(ligs_to_dock_path)
    for mol in mols:
        writer.write(mol)
    writer.close()

    # Specify the pocket coordinates and run the docking
    box = pocket.pocket_location.get_params("boxcxyz")
    docked_ligs_path = os.path.join(c["work_dir"], "docked_ligs.sdf")
    subprocess.run([
        c["bin_path"], "-l", ligs_to_dock_path, "-r", pocket.pdb_path, "--center_x", str(box['c'][0]), 
        "--center_y", str(box['c'][1]), "--center_z", str(box['c'][2]), "--size_x", str(box['size'][0]), 
        "--size_y", str(box['size'][1]), "--size_z", str(box['size'][2]), "-o", docked_ligs_path, 
    ] + c["other_args"], capture_output=True)

    # Load the saved results (.sdf) and return
    # NOTE: Output has no hydrogens...
    try:
        docked_ligs = [mol for mol in Chem.SDMolSupplier(docked_ligs_path)]
        return docked_ligs
    except OSError:
        return []

