import os, subprocess
import pathlib, shutil, copy, glob
from typing import Dict

import pandas as pd

from rdkit import Chem
from biopandas.pdb import PandasPdb
from openmm import NonbondedForce, unit
from openmm.app import PDBFile, ForceField
from openff.toolkit import Topology
from openff.units import unit as openff_unit

from .common import calc_protein_fcharge
from .mm import prepare_lig_for_openmm, prepare_target_pdb_for_openmm, STARTING_FORCEFIELDS

class DefaultFormat(dict):
    def __missing__(self, key):
        return "{" + key + "}"
 
WORK_DIR = "tmp"
N_THREADS = 8
MOPAC_BIN_PATH = str(pathlib.Path(shutil.which("mopac")).absolute())
AMBER_HOME = "/home/arazthexd/miniforge3/envs/sqmvscreen2/"
QMMM_OPT_MAX_CYCLES = 10

H2O_G = -50.06
H3O_G = -33.03
HO_G = 150.12
H2O_GS = -6.3
H3O_GS = -110.2
HO_GS = -105.0

### DEFAULT CONFIGS ###
CUBY4_DEFAULT_INTERACTION_CONFIG = """
job: interaction
geometry: {geometry}

molecule_a:
  selection: '{select_a}'
  charge: {charge_a}
  mopac_setcharge:
    {setcharge_lines_a}
molecule_b:
  selection: '{select_b}'
  charge: {charge_b}

interface: mopac
method: pm6
mopac_exe: {mopac_exe}
mopac_mozyme: yes
mopac_corrections: d3h4x
mopac_setcharge:
  {setcharge_lines_all}
solvent_model: {solv_model}

cuby_threads: {n_threads}
""".format_map(DefaultFormat(mopac_exe=MOPAC_BIN_PATH, n_threads=N_THREADS))

CUBY4_DEFAULT_ENERGY_CONFIG = """
job: energy
geometry: {geometry}

interface: mopac
method: pm6
mopac_exe: {mopac_exe}
mopac_mozyme: yes
mopac_corrections: d3h4x
mopac_setcharge:
  {setcharge_lines}
solvent_model: {solv_model}
charge: {charge}

cuby_threads: {n_threads}
""".format_map(DefaultFormat(mopac_exe=MOPAC_BIN_PATH, n_threads=N_THREADS))

CUBY4_DEFAULT_OPTIMIZE_CONFIG = """
job: optimize
geometry: {geometry}

interface: mopac
method: pm6
mopac_exe: {mopac_exe}
mopac_mozyme: yes
mopac_corrections: d3h4x
mopac_setcharge:
  {setcharge_lines}
charge: {charge}

restart_file: {restart_path}
cuby_threads: {n_threads}
""".format_map(DefaultFormat(mopac_exe=MOPAC_BIN_PATH, n_threads=N_THREADS))

CUBY4_DEFAULT_QMMM_CONFIG = """
job: optimize
geometry: {geometry}
maxcycles: {max_cycles}

interface: qmmm
qmmm_core: {qmmm_core}
qmmm_embedding: mechanical
gradient_on_point_charges: 'no'

calculation_qm:
  interface: mopac
  method: pm6
  mopac_exe: {mopac_exe}
  mopac_mozyme: 'yes'
  mopac_corrections: d3h4x
  mopac_setcharge:
    {setcharge_lines_qm}
  gradient_on_point_charges: 'no'
  charge: {charge_qm}

calculation_qmregion_mm:
  amber_top_file: {qmregion_top}

calculation_mm:
  interface: amber
  amber_amberhome: {amber_home}
  amber_top_file: {complex_top}

restart_file: {restart_path}
cuby_threads: {n_threads}
""".format_map(DefaultFormat(mopac_exe=MOPAC_BIN_PATH, n_threads=N_THREADS, 
                             amber_home=AMBER_HOME, max_cycles=QMMM_OPT_MAX_CYCLES))
### END OF DEFAULT CONFIGS ###

def dissociate_complex(complex_path, lig_name, lig_write, poc_write):
    compl = open(complex_path, "r")
    lig = open(lig_write, "w")
    poc = open(poc_write, "w")
    for line in compl.readlines():
        if lig_name in line:
            lig.write(line)
        else:
            poc.write(line)
    lig.write("END\n")
    [x.close() for x in [compl, lig, poc]]

def prepare_mol_and_target_for_cuby(mol: Chem.Mol, target_pdb: str, job_name: str = "job"):

    mol = Chem.AddHs(mol, addCoords=True)
    ppdb_target = PandasPdb().read_pdb(target_pdb)
    ppdb_ligand = PandasPdb().read_pdb_from_list(
        Chem.MolToPDBBlock(mol).split("\n")
    )
    lig_res_name = ppdb_ligand.df["HETATM"]["residue_name"].unique()[0]

    ppdb_target.df["HETATM"] = ppdb_ligand.df["HETATM"]
    target_atom_num = len(ppdb_target.df["ATOM"])
    target_res_max = ppdb_target.df["ATOM"]["residue_number"].max()
    ppdb_target.df["HETATM"]["atom_number"] = \
        ppdb_target.df["HETATM"]["atom_number"].apply(lambda n: n + target_atom_num)
    ppdb_target.df["HETATM"]["line_idx"] = \
        ppdb_target.df["HETATM"]["line_idx"].apply(lambda n: n + target_atom_num)
    ppdb_target.df["HETATM"]["residue_number"] = \
        ppdb_target.df["HETATM"]["residue_number"].apply(lambda n: n + target_res_max)
    out_path = pathlib.Path(os.path.join(WORK_DIR, f"cuby_{job_name}_complex.pdb"))
    ppdb_target.to_pdb(out_path.absolute(), ["ATOM", "HETATM"])

    return {"lig_name": lig_res_name, "geometry_path": str(out_path.absolute())}

def calculate_interaction(mol: Chem.Mol, target_pdb: str, job_name: str = "job"):

    o = prepare_mol_and_target_for_cuby(mol, target_pdb, job_name)
    geometry_path = o["geometry_path"]
    lig_res_name = o["lig_name"]
    target_top = Topology.from_pdb(target_pdb)

    config_str = CUBY4_DEFAULT_INTERACTION_CONFIG
    config_str = config_str.format_map(DefaultFormat(
        geometry=geometry_path,
        select_a=":"+lig_res_name,
        charge_a=Chem.GetFormalCharge(mol),
        select_b=f"%not({':'+lig_res_name})",
        charge_b=calc_protein_fcharge(target_pdb),
        setcharge_lines_a="\n    ".join([f"{a.GetIdx()+1}: {str(a.GetFormalCharge())}" for a in mol.GetAtoms()]),
        setcharge_lines_b="\n    ".join([f"{i+1}: {a.formal_charge.m_as(openff_unit.elementary_charge)}" for i, a in enumerate(target_top.atoms)]),
        setcharge_lines_all="\n  ".join([f"{a.GetIdx()+1}: {str(a.GetFormalCharge())}" for a in mol.GetAtoms()]),
        solv_model="cosmo2"
    ))

    output_str = run_cuby(config_str, full_job_name=f"inter_{job_name}")
    interaction_energy = float(output_str.split("Interaction energy: ")[-1][:7])
    return interaction_energy

def calculate_energy(mol: Chem.Mol, job_name: str = "job"):

    geometry_path = os.path.join(WORK_DIR, f"cuby_energy_{job_name}_ligand.pdb")
    geometry_path = str(pathlib.Path(geometry_path).absolute())
    Chem.MolToPDBFile(mol, geometry_path)

    config_str = CUBY4_DEFAULT_ENERGY_CONFIG
    config_str = config_str.format_map(DefaultFormat(
        geometry=geometry_path,
        setcharge_lines="\n  ".join([f"{a.GetIdx()+1}: {str(a.GetFormalCharge())}" for a in mol.GetAtoms()]),
        charge=Chem.GetFormalCharge(mol),
        solv_model="cosmo2"
    ))

    output_str = run_cuby(config_str, full_job_name=f"energy_{job_name}")
    energy = float(output_str.split("nergy: ")[-1][:7])
    return energy

def optimize_molecule_qm(mol: Chem.Mol, job_name: str = "job"):

    geometry_path = os.path.join(WORK_DIR, f"cuby_optimize_{job_name}_ligand.pdb")
    geometry_path = str(pathlib.Path(geometry_path).absolute())
    Chem.MolToPDBFile(mol, geometry_path)

    config_str = CUBY4_DEFAULT_OPTIMIZE_CONFIG
    restart_file = f"qmmm_{job_name}_restart.pdb"
    config_str = config_str.format_map(DefaultFormat(
        geometry=geometry_path,
        setcharge_lines="\n  ".join([f"{a.GetIdx()+1}: {str(a.GetFormalCharge())}" for a in mol.GetAtoms()]),
        charge=Chem.GetFormalCharge(mol),
        restart_path=restart_file
    ))

    output_str = run_cuby(config_str, full_job_name=f"optimize_{job_name}")
    latest_energy = float(output_str.split("nergy: ")[-1][:7])

    newconf = Chem.MolFromPDBFile(os.path.join(WORK_DIR, restart_file), 
                                         removeHs=False).GetConformer()
    mol_new = Chem.Mol(mol)
    mol_new.RemoveAllConformers()
    mol_new.AddConformer(newconf)

    out_dict = dict()
    out_dict["mol_rdkit"] = mol_new
    out_dict["energy"] = latest_energy
    return out_dict

def optimize_molecule_with_pocket(mol: Chem.Mol, target_pdb: str, job_name: str = "job"):

    o = prepare_mol_and_target_for_cuby(mol, target_pdb, job_name)
    geometry_path = o["geometry_path"]
    lig_res_name = o["lig_name"]

    forcefield = ForceField(*STARTING_FORCEFIELDS)
    lig_mm_related = prepare_lig_for_openmm(mol, forcefield, job_name)
    lig_struct = lig_mm_related["parmed_structure"]
    lig_prm_file = f"qmmm_{job_name}_ligand.parm7"
    lig_struct.save(os.path.join(WORK_DIR, lig_prm_file), overwrite=True)

    targ_mm_related = prepare_target_pdb_for_openmm(target_pdb, forcefield)
    targ_struct = targ_mm_related["parmed_structure"]

    complex_struct = targ_struct + lig_struct
    # complex_pdb_path = os.path.join(WORK_DIR, f"qmmm_{job_name}_complex.pdb")
    complex_prm_file = f"qmmm_{job_name}_complex.parm7"
    complex_struct.save(os.path.join(WORK_DIR, complex_prm_file), overwrite=True)

    config_str = CUBY4_DEFAULT_QMMM_CONFIG
    restart_file = f"qmmm_{job_name}_restart.pdb"
    config_str = config_str.format_map(DefaultFormat(
        geometry=geometry_path,
        qmmm_core=":"+lig_res_name,
        charge_qm=Chem.GetFormalCharge(mol),
        setcharge_lines_qm="\n    ".join([f"{a.GetIdx()+1}: {str(a.GetFormalCharge())}" for a in mol.GetAtoms()]),
        qmregion_top=lig_prm_file,
        complex_top=complex_prm_file,
        restart_path=restart_file
    ))
    
    output_str = run_cuby(config_str, full_job_name=f"qmmm_{job_name}")

    latest_energy = float(output_str.split("Energy: ")[-1][:7])
    ligopted_path = os.path.join(WORK_DIR, f"qmmm_{job_name}_ligopted.pdb")
    pocopted_path = os.path.join(WORK_DIR, f"qmmm_{job_name}_pocopted.pdb")
    dissociate_complex(
        os.path.join(WORK_DIR, restart_file), lig_res_name, ligopted_path, pocopted_path
    )
    
    ligand_newconf = Chem.MolFromPDBFile(ligopted_path, removeHs=False).GetConformer()
    ligand_new = Chem.Mol(mol)
    ligand_new.RemoveAllConformers()
    ligand_new.AddConformer(ligand_newconf)
    
    out_dict = dict()
    out_dict["energy"] = latest_energy
    out_dict["ligand_rdkit"] = ligand_new
    out_dict["ligand_pdb"] = pocopted_path
    out_dict["pocket_pdb"] = pocopted_path
    out_dict["complex_pdb"] = os.path.join(WORK_DIR, restart_file)
    return out_dict
    

def run_cuby(config_str, full_job_name):
    
    cwd = str(pathlib.Path(os.curdir).absolute())
    os.chdir(WORK_DIR)
    config_path = f"cuby_{full_job_name}_config.yml"
    with open(config_path, "w") as f:
        f.write(config_str)

    for filename in glob.glob("job_*"):
        shutil.rmtree(filename) 
    output = subprocess.run(["cuby4", config_path], capture_output=True)

    os.chdir(cwd)
    return output.stdout.decode()

def calc_hydrogen_transfer_energy(num_hs_added: int):
    if num_hs_added > 0:
        return num_hs_added * (HO_G + HO_GS - H2O_G - H2O_GS)
    else:
        return num_hs_added * (H3O_G + H3O_GS - H2O_G - H2O_GS)

# m = next(Chem.SDMolSupplier("/home/arazthexd/projects/002_sqm/output/gbminimized_ligands.sdf", removeHs=False))
# optimize_molecule_with_pocket(m, m.GetProp("target_path"))