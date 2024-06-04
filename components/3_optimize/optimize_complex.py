import argparse
import os, glob, shutil, pathlib
import yaml
import warnings
warnings.simplefilter('ignore')

from collections import OrderedDict
import subprocess
from tqdm import tqdm

from rdkit import Chem

from openff.toolkit import Molecule
from openff.units import unit
from openmm.app import PDBFile, ForceField, Modeller
from openmmforcefields.generators import SMIRNOFFTemplateGenerator
import parmed
from pdbfixer import PDBFixer

from ezdd.utils.protein.pockets import PocketLocation, isolate_pocket

config = {
    "job": "optimize",
    "geometry": None,
    "interface": "qmmm",
    "qmmm_core": ":UNL",
    "qmmm_embedding": "mechanical",
    "gradient_on_point_charges": "no",
    "calculation_mm": {
        "interface": "amber",
        "amber_amberhome": str(pathlib.Path(shutil.which("sander")).parent.parent.absolute()),
        "amber_top_file": None,
    },
    "calculation_qmregion_mm": {
        "amber_top_file": None,
    },
    "calculation_qm": {
        "interface": "mopac",
        "method": "pm6",
        "charge": None,
        "mopac_exe": str(pathlib.Path(shutil.which("mopac")).absolute()),
        "gradient_on_point_charges": "no",
        "mopac_mozyme": "yes"
    },
    "cuby_threads": 8,
    "restart_file": None
}

pdbgen_config = {
    "job": "geometry",

}

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(
        prog="Complex (Pocket) Optimizer",
    )
    parser.add_argument("-l", "--ligands", help="path to ligands file (.sdf)")
    parser.add_argument("-p", "--protein", help="path to protein file (.pdb)")
    parser.add_argument("-o", "--output", help="directory for optimized files (folder -> .pdb)")
    parser.add_argument("--work_dir", help="where temporary files are stored", default="tmp")
    args = parser.parse_args()

    if not os.path.exists(args.work_dir):
        os.mkdir(args.work_dir)
    if not os.path.exists(args.output):
        os.mkdir(args.output)

    print("\nRemoving Potential Previous Jobs...")
    for filename in glob.glob("job_*"):
        shutil.rmtree(filename) 

    print("\nBeginning Iteration On Ligands...")
    ligand_suppl = Chem.SDMolSupplier(args.ligands)
    for i, rdligand in enumerate(ligand_suppl):
        print(f"\n### Ligand No. {i} ###")
        rdligand = Chem.AddHs(rdligand, addCoords=True)

        print("\nPreprocessing Ligand Pocket...")
        pocket_loc = PocketLocation("ligand", ligand=rdligand, radius=10)
        poc_path = os.path.join(args.work_dir, f"pocket{i}.pdb")
        isolate_pocket(args.protein, pocket_loc, poc_path)
        poc_pdb = PDBFixer(poc_path)
        poc_pdb.addMissingHydrogens(7.0)
        poc_ff = ForceField("amber14-all.xml")
        poc = poc_ff.createSystem(poc_pdb.topology)
        rec_struct = parmed.openmm.load_topology(
            topology=poc_pdb.topology,
            system=poc,
            xyz=poc_pdb.positions
        )

        print("\nCreating Ligand Parameters...")
        ligand_molecule = Molecule.from_rdkit(rdligand, allow_undefined_stereo=True)
        ligand_path = os.path.join(args.work_dir, f"ligand{i}.pdb")
        ligand_molecule.to_file(ligand_path, "pdb")

        lig_ff = ForceField()
        smirnoff = SMIRNOFFTemplateGenerator(molecules=ligand_molecule)
        lig_ff.registerTemplateGenerator(smirnoff.generator)

        lig_pdb = PDBFile(ligand_path)
        lig = lig_ff.createSystem(lig_pdb.topology)

        lig_struct = parmed.openmm.load_topology(
            topology=lig_pdb.topology,
            system=lig,
            xyz=lig_pdb.positions
        )

        print("\nGenerating And Saving PL Complex Params...")
        os.chdir(args.work_dir)

        complex_struct = rec_struct + lig_struct
        complex_pdb_path = f"complex{i}.pdb" # os.path.join(args.work_dir, f"complex{i}.pdb")
        complex_prm_path = f"complex{i}.parm7" # os.path.join(args.work_dir, f"complex{i}.parm7")
        ligand_prm_path = f"ligand{i}.parm7" # os.path.join(args.work_dir, f"ligand{i}.parm7")
        output_path = f"optimized{i}.pdb"
        complex_struct.save(complex_pdb_path, overwrite=True)
        complex_struct.save(complex_prm_path, overwrite=True)
        lig_struct.save(ligand_prm_path, overwrite=True)

        print("\nSaving And Executing Optimization Config...")
        config["geometry"] = complex_pdb_path
        config["calculation_mm"]["amber_top_file"] = complex_prm_path
        config["calculation_qmregion_mm"]["amber_top_file"] = ligand_prm_path
        config["calculation_qm"]["charge"] = ligand_molecule.total_charge.m_as(unit.elementary_charge)
        config["restart_file"] = output_path
        # TEMPORARY
        config["maxcycles"] = 2

        config_path = f"complex{i}.yml" # os.path.join(args.work_dir, f"complex{i}.yml")
        with open(config_path, "w") as f:
            yaml.dump(config, f)
        
        o = subprocess.run(["cuby4", config_path])
        os.chdir("..")
        shutil.move(os.path.join(args.work_dir, f"optimized{i}.pdb"), 
                    os.path.join(args.output, f"optimized{i}.pdb"))
        
        



        