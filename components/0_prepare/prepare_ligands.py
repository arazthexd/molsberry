import argparse
import subprocess
import sys
import os

import pandas as pd

import rdkit
from rdkit import Chem
from rdkit.Chem import rdDistGeom, rdMolDescriptors, PandasTools, rdForceFieldHelpers

ROTOR_QUERY = Chem.MolFromSmarts("[!$(*#*)&!D1]-!@[!$(*#*)&!D1]")

def smi2mols(input_file):
    with open(input_file, "r") as f:
        smi_name_list = f.readlines()
    
    i = 0
    mols = []
    for smi_name in smi_name_list:
        s = smi_name.split()
        smi = s[0]
        if len(s) == 2:
            name = s[1]
        else:
            name = "unnamed" + str(i)
            i += 1
        
        mol = Chem.MolFromSmiles(smi)
        mol.SetProp("_Name", name)
        if rdMolDescriptors.CalcExactMolWt(mol) < 600: # TODO: Set this...
            mols.append(mol) 
    
    return mols

def sdf2mols():
    pass

def clean_mols():
    pass

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(
        prog="Ligand Preparer",
    )
    parser.add_argument("filename", help="path to input file (.smi or .sdf)")
    parser.add_argument("outfile", help="path to output file (.sdf)")
    parser.add_argument("-v", "--verbose", action="store_true")
    parser.add_argument("--work_dir", help="where temporary files are stored", default="tmp")
    parser.add_argument("--gypsum_dl_dir", help="if not specified, a simple preprocessing will be used", default="none")
    parser.add_argument("--num_proc", help="number of processors to use", default=1)
    
    args = parser.parse_args()
    v = args.verbose

    assert args.filename.endswith(".sdf") or args.filename.endswith(".smi")
    assert args.outfile.endswith(".sdf")

    if v: print("RDKit Version Using:", rdkit.__version__)
    
    if not os.path.exists(args.work_dir):
        if v: print("Creating work directory...")
        os.mkdir(args.work_dir)

    # if args.filename.endswith(".smi"):
    #     mols = smi2mols(args.filename)
    # else:
    #     raise NotImplementedError
    # primary_ligs_path = f"{os.path.join(args.work_dir, 'primary_ligands.sdf')}"
    # writer = Chem.SDWriter(primary_ligs_path)
    # [writer.write(mol) for mol in mols]

    gypsum_ligs_path = f"{os.path.join(args.work_dir, 'gypsum_ligands.sdf')}"
    if not args.gypsum_dl_dir == "none":
        sys.path.append(args.gypsum_dl_dir)
        gypsum_path = os.path.join(args.gypsum_dl_dir, 'run_gypsum_dl.py')
        cmd_line = f"python {gypsum_path} --source {args.filename} \
            --use_durrant_lab_filters --job_manager multiprocessing --num_processors {args.num_proc} \
                --min_ph 5 --max_ph 9 --skip_optimize_geometry"
        subprocess.run(cmd_line.split())
        os.rename("gypsum_dl_success.sdf", gypsum_ligs_path)
        if os.path.exists("gypsum_dl_failed.smi"):
            os.rename("gypsum_dl_failed.smi", os.path.join(args.work_dir, "gypsum_fails.smi"))
    else:
        raise NotImplementedError

    suppl = Chem.SDMolSupplier(gypsum_ligs_path, removeHs=True)
    writer = Chem.SDWriter(args.outfile)
    for i, mol in enumerate(suppl):
        if i == 0: continue
        print(Chem.MolToSmiles(mol))
        mol = Chem.AddHs(mol)
        rdForceFieldHelpers.MMFFOptimizeMolecule(mol)
        writer.write(mol)