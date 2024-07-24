from typing import List
from tqdm import tqdm
from copy import deepcopy
import string, random, glob

from rdkit import Chem
from rdkit.Chem import (
    rdDistGeom, 
    rdForceFieldHelpers, 
)

from .moltools import addhs_based_on_confdim

def write_ligands(ligands: List[Chem.Mol], save_path):
    writer = Chem.SDWriter(save_path)
    for lig in ligands:
        writer.write(lig)
    writer.close()
    return save_path

def load_ligands(ligands_path: str, final_3d: bool = False, addHs: bool = False):

    if ligands_path.endswith(".smi"):
        ligs: List[Chem.Mol] = _load_ligands_smi(ligands_path)
    elif ligands_path.endswith(".sdf"):
        ligs: List[Chem.Mol] = list(Chem.SDMolSupplier(ligands_path))
    else:
        raise NotImplementedError()

    if addHs:
        ligs = [addhs_based_on_confdim(lig) for lig in ligs]
    
    if final_3d:
        [rdDistGeom.EmbedMolecule(lig) for lig in ligs if lig.GetNumConformers() == 0]
        [rdDistGeom.EmbedMolecule(lig) for lig in ligs if not lig.GetConformer().Is3D()]
        [rdForceFieldHelpers.MMFFOptimizeMolecule(lig) for lig in ligs]
    
    return ligs

def _load_ligands_smi(ligands_path: str):

    i = 0
    smis, names = [], []
    with open(ligands_path, "r") as f:
        for line in f.readlines():
            splitted = line.split()
            if len(splitted) == 1:
                smi, name = splitted[0], f"unnamed{i}"
                i += 1
            elif len(splitted) == 2:
                smi, name = splitted[0], splitted[1]
            elif len(splitted) > 2:
                smi = splitted[0]
                name = "".join(splitted[1:])
            smis.append(smi)
            names.append(name)
    
    ligs = []
    for smi, name in zip(smis, names):
        lig = Chem.MolFromSmiles(smi)
        lig.SetProp("_Name", name)
        lig.SetProp("InitialSmiles", smi)
        ligs.append(lig)
    
    return ligs

def save_pl_complex(ligand: Chem.Mol, target_path: str, out_path: str):
    with open(target_path, "r") as f:
        target_pdb_lines = f.readlines()
        target_pdb_lines = [l.split("\n")[0] for l in target_pdb_lines]
        n_target_atoms = 0
        max_res_num = 0
        for line in target_pdb_lines:
            if "HETATM" in line or "ATOM" in line: 
                n_target_atoms += 1
                max_res_num = max(max_res_num, int(line[22:26]))
    
    ligand_pdb_lines = Chem.MolToPDBBlock(ligand).split("\n")
    for i, line in enumerate(deepcopy(ligand_pdb_lines)):
        if "HETATM" in line:
            line = (
                line[:6] + str(int(line[6:11])+n_target_atoms).rjust(5) + 
                line[11:22] + str(int(line[22:26])+max_res_num).rjust(4) +
                line[26:]
            )       
        if "CONECT" in line:
            indices = [int(line[i:i+5])
                        for i in range(6, len(line), 5)]
            line = line[:6] + "".join([str(idx+n_target_atoms).rjust(5) for idx in indices])
        ligand_pdb_lines[i] = line

    all_lines = target_pdb_lines + ligand_pdb_lines
    with open(out_path, "w") as f:
        f.write("\n".join([line for line in all_lines 
                           if "CONECT" not in line and "END" not in line and "TER" not in line]))
        f.write("\n".join([line for line in all_lines if "CONECT" in line]))

def generate_random_str(n: int):
    return "".join(random.choices(string.ascii_uppercase + string.digits, k=n))

def get_unique_full_name(prefix: str, extension: str, n: int):
    sim_files = glob.glob(prefix + "_*")
    while True:
        rnd = generate_random_str(n)
        suggestion = prefix + "_" + rnd + extension
        if suggestion not in sim_files:
            return suggestion