import os
from tqdm import tqdm
from typing import List
import warnings
warnings.simplefilter("ignore")

import contextlib, io, joblib
import time

import pandas as pd
import numpy as np

from rdkit import Chem
from rdkit.Chem import PandasTools
from rdkit.Chem.Descriptors import rdMolDescriptors

from ezdd.utils.protein.pockets import PocketLocation, isolate_pocket
from ezdd.utils.protein.prepare import clean_protein

from utils import ligprep, sqm, vs, mm, baseline
from utils.common import GenericPocket
from utils.entropy.featurize import calculate_features

# CONFIG, VARIABLES
TMP_DIR = "tmp"
OUT_DIR = "output"
PRIMARY_LIGS_PATH = "data/ligands/tmp/test_ligands.smi"
PRIMARY_LIGS_MAX_WEIGHT = 600
PROTEIN_PATH = "output/prepared_protein.pdb"
PROTEIN_CHAINS = "all"
GNINA_POCKET_CENTER = (1.16, -0.56,- 2.91)
GNINA_POCKET_SIZE = (20, 20, 20)
GNINA_PATH = "/home/arazthexd/tools/gnina/gnina"
PRIMARY_VS_BATCH_SIZE = 5
PRIMARY_VS_EXHAUSTIVENESS = 10
PRIMARY_VS_NUMMODES = 5
PRIMARY_VS_SELECT_RATIO = 0.1
PRIMARY_VS_SELECT_SORT_BY = "CNN_VS"
PRIMARY_VS_SELECT_VINAAFFINITY_THRESH = -4
PRIMARY_VS_SELECT_CNNPOSESCORE_THRESH = 0.4
PRIMARY_VS_SELECT_MAX_PER_ORIG_LIG = 2
SECONDARY_GBSCREEN_PH = 7.4
SECONDARY_GBSCREEN_PH_THRESH = 3
SECONDARY_GBSCREEN_SELECT_MAX_PER_LIG = 1
SECONDARY_GBSCREEN_SELECT_RATIO = 0.1

ENTROPY_MODEL = "/home/arazthexd/projects/002_sqm/models/mlp8_tanh.pkl"


def prepare_primary_ligands(smi_path):
    mols = []
    with open(smi_path, "r") as f:
        smis, names = zip(*[line.split() for line in f.readlines()])
    for smi, name in tqdm(zip(smis, names), total=len(smis)):
        lig = Chem.MolFromSmiles(smi)
        if rdMolDescriptors.CalcExactMolWt(lig) > PRIMARY_LIGS_MAX_WEIGHT:
            continue
        lig.SetProp("_Name", name)
        lig.SetProp("InitialSmiles", smi)
        tautomers = ligprep.enumerate_tautomers(lig)
        stereoisos = [mol for t in tautomers for mol in ligprep.enumerate_stereoisomers(t)]
        try:
            mols.extend([mol for stereoiso in stereoisos for mol in ligprep.enumerate_ring_confs(stereoiso, minimize=True)])
        except:
            break
    return mols, len(smis)

def select_primary_docked(docked_ligs, n_select):
    # n_select = int(len(ligs_to_dock) * PRIMARY_VS_SELECT_RATIO) + 1
    selected_ligs = []
    selected_names = dict()
    for mol in sorted(docked_ligs, key=lambda lig: float(lig.GetProp(PRIMARY_VS_SELECT_SORT_BY)), reverse=True):
        if len(selected_ligs) >= n_select:
            break
        if float(mol.GetProp("minimizedAffinity")) > PRIMARY_VS_SELECT_VINAAFFINITY_THRESH:
            continue
        if float(mol.GetProp("CNNscore")) < PRIMARY_VS_SELECT_CNNPOSESCORE_THRESH:
            continue
        name = mol.GetProp("_Name")
        if selected_names.get(name) is not None:
            if selected_names.get(name) >= PRIMARY_VS_SELECT_MAX_PER_ORIG_LIG:
                continue
        else:
            selected_names[name] = 0
        selected_ligs.append(mol)
        selected_names[name] += 1
    return selected_ligs

def select_secondary_minimized(minimized_ligands, n_select):
    # n_select = int(len(ligs_to_dock) * SECONDARY_GBSCREEN_SELECT_RATIO) + 1
    selected_ligs = []
    selected_names = dict()
    for mol in sorted(minimized_ligands, key=lambda lig: float(lig.GetProp("GBInteractionEnergy")), reverse=False):
        if len(selected_ligs) >= n_select:
            break
        name = mol.GetProp("_Name")
        if selected_names.get(name) is not None:
            if selected_names.get(name) >= SECONDARY_GBSCREEN_SELECT_MAX_PER_LIG:
                continue
        else:
            selected_names[name] = 0
        selected_ligs.append(mol)
        selected_names[name] += 1
    return selected_ligs

def main():

    third_party_out = io.StringIO()
    tick = time.time()

    # 0) Prepare protein file
    print("### PROTEIN PREPARATION ###")
    prepped_protein_path = os.path.join(OUT_DIR, "prepared_protein.pdb")
    clean_protein(PROTEIN_PATH, PROTEIN_CHAINS, prepped_protein_path)
    tock_0 = time.time()
    print(f"Done! Saved Output: {prepped_protein_path}")
    print(f"Time Spent: {tock_0-tick:.2f}s")

    # 1) Prepare primary ligands by enumemration of tautomers, stereoisomers and ring conformers
    print("\n### PRIMARY LIGAND PREPARATION & ENUMERATION ###")
    print("Preparing Ligands (Note: Primary, No Protomers)...")
    # NOTE: ligs_to_dock contains hydrogen and every other definition needed but no protomers included
    ligs_to_dock, n_original_ligs = prepare_primary_ligands(PRIMARY_LIGS_PATH)
    tock_1 = time.time()
    print(f"Done! Number Of Ligands To Dock: {len(ligs_to_dock)}")
    print(f"Time Spent: {tock_1-tock_0:.2f}s")
    print(f"Time Spent Per Original Lig: {(tock_1-tock_0)/n_original_ligs:.2f}")


    # 2) Perform primary virtual screening on primary ligands
    print("\n### PRIMARY VIRTUAL SCREENING ###")
    lig_batches = [ligs_to_dock[i:i+PRIMARY_VS_BATCH_SIZE] for i in range(0, len(ligs_to_dock), PRIMARY_VS_BATCH_SIZE)]
    print(f"Number Of Batches (Size {PRIMARY_VS_BATCH_SIZE}): {len(lig_batches)}")
    # NOTE: docked_ligs, after being returned from gnina has no Hs (except some)
    # TODO: Make sure what is returned always has Hs
    gnina_pocket = GenericPocket(prepped_protein_path, PocketLocation("boxcxyz", center=GNINA_POCKET_CENTER, xyz_size=GNINA_POCKET_SIZE))
    print("Performing Ligand Docking...")
    docked_ligs: List[Chem.Mol] = []
    for lig_batch in tqdm(lig_batches):
        docked_ligs += vs.perform_docking(lig_batch, gnina_pocket, config={
            "bin_path": GNINA_PATH, 
            "other_args": f"--exhaustiveness {PRIMARY_VS_EXHAUSTIVENESS} --num_modes {PRIMARY_VS_NUMMODES}"
        })
    docked_ligs_path = os.path.join(OUT_DIR, "docked_ligands.sdf")
    writer = Chem.SDWriter(docked_ligs_path)
    [writer.write(l) for l in docked_ligs]
    writer.close()
    tock_2 = time.time()
    print(f"Done! Number Of Resulting Conformations: {len(docked_ligs)}")
    print(f"Time Spent: {tock_2-tock_1:.2f}s")
    print(f"Time Spent Per Original Lig: {(tock_2-tock_1)/n_original_ligs:.2f}")


    # 3) Select primary docked ligands with the provided ratio of selection
    print("\n### PRIMARY DOCKED LIGAND SELECTION ###")
    print(f"Selecting {PRIMARY_VS_SELECT_RATIO*100}% Of Primary Docked Ligands...")
    n_select = int(len(ligs_to_dock) * PRIMARY_VS_SELECT_RATIO) + 1
    selected_ligs = select_primary_docked(docked_ligs, n_select)
    tock_3 = time.time()
    print(f"Done! Number Of Selected Ligands: {len(selected_ligs)}")
    print(f"Time Spent: {tock_3-tock_2:.2f}s")
    

    # 4) Prepare selected ligands for the second run of virtual screening by enumerating their protomers
    print("\n### PROTOMER ENUMERATION OF SELECTED LIGANDS ###")
    print("Generating Protomers Of Selected Ligands...")
    ph_min = SECONDARY_GBSCREEN_PH-SECONDARY_GBSCREEN_PH_THRESH
    ph_max = SECONDARY_GBSCREEN_PH+SECONDARY_GBSCREEN_PH_THRESH
    print(f"pH Range: {ph_min}-{ph_max}")
    with contextlib.redirect_stdout(third_party_out):
        ligs_to_minimize = [protomer for lig in selected_ligs for protomer in ligprep.enumerate_protomers(
            lig, min_ph=ph_min, max_ph=ph_max)]
    tock_4 = time.time()
    print(f"Done! Number Of Protomers Ready To Be Minimized: {len(ligs_to_minimize)}")
    print(f"Time Spent: {tock_4-tock_3:.2f}s")


    # 5) Minimize and Rescore with MM-GBSA
    print("\n### PROTOMERS GB MINIMIZATION ###")
    print("Minimizing Ligands Based On AMBER Forcefield + GB Solvent Model...")
    minimized_ligs: List[Chem.Mol] = []
    minimized_ligs_path = os.path.join(OUT_DIR, f"gbminimized_ligands.sdf")
    writer = Chem.SDWriter(minimized_ligs_path)
    for i, lig in tqdm(enumerate(ligs_to_minimize), total=len(ligs_to_minimize)):
        lig = Chem.AddHs(lig, addCoords=True)
        output = mm.minimize_and_calculate_interaction(lig, prepped_protein_path, f"{i}_{lig.GetProp('_Name')}",
                                                    extra_outputs=["ligand_rdkit", "pocket_pdb", "complex_pdb"])
        minimized_lig: Chem.Mol = output["ligand_rdkit"]
        minimized_lig.SetDoubleProp("GBInteractionEnergy", output["interaction_energy"])
        minimized_lig.SetProp("target_path", output["pocket_pdb"])
        minimized_lig.SetProp("complex_path", output["complex_pdb"])
        minimized_ligs.append(minimized_lig)
        writer.write(minimized_lig)
    writer.close()
    tock_5 = time.time()
    print(f"Done! Saved Results To: {minimized_ligs_path}")
    print(f"Time Spent: {tock_5-tock_4:.2f}s")
    print(f"Time Spent Per Original Lig: {(tock_5-tock_4)/n_original_ligs:.2f}")


    # 6) Select secondary minimized ligands based on interaction energy in solvated system
    print("\n### MINIMIZED PROTOMER SELECTION ###")
    print(f"Selecting {SECONDARY_GBSCREEN_SELECT_RATIO*100}% Of Protomers Based On Minimized Structures...")
    n_select = int(len(ligs_to_dock) * SECONDARY_GBSCREEN_SELECT_RATIO) + 1
    selected_ligs = select_secondary_minimized(minimized_ligs, n_select)
    tock_6 = time.time()
    print(f"Done! Number Of Selected Protomers: {len(selected_ligs)}")
    print(f"Time Spent: {tock_6-tock_5:.2f}s")

    # 7) Optimize the selected ligands using Cuby4's qmmm interface (mopac pm6, amber)
    print("\n### QMMM OPTIMIZATION OF PL COMPLEXES ###")
    print("Optimizing Complex Structures Using QMMM Interface...")
    optimized_ligs = []
    optimized_ligs_path = os.path.join(OUT_DIR, f"qmmmopt_ligands.sdf")
    writer = Chem.SDWriter(optimized_ligs_path)
    for i, lig in enumerate(tqdm(selected_ligs)):
        target_path = lig.GetProp("target_path")
        isopocket_path = os.path.join(TMP_DIR, f"isopocket_{i}_{lig.GetProp('_Name')}.pdb")
        isolate_pocket(target_path, PocketLocation("ligand", ligand=lig, radius=10), 
                    output_path=isopocket_path)
        output = sqm.optimize_molecule_with_pocket(lig, isopocket_path, f"{i}_{lig.GetProp('_Name')}")
        optimized_lig = output["ligand_rdkit"]
        optimized_lig.SetProp("target_path", output["pocket_pdb"])
        optimized_lig.SetDoubleProp("optimized_energy", output["energy"])
        writer.write(optimized_lig)
        optimized_ligs.append(optimized_lig)
    tock_7 = time.time()
    print(f"Done! Results Saved To: {optimized_ligs_path}")
    print(f"Time Spent: {tock_7-tock_6:.2f}s")
    print(f"Time Spent Per Original Lig: {(tock_7-tock_6)/n_original_ligs:.2f}")


    # 8) Calculate interaction energy between ligands and pockets
    print()
    print("#############################################")
    print("########## FINAL SCORE CALCULATION ##########")
    print("#############################################")
    print()

    print("1) INTERACTION ENERGY")
    print("Calculating Interactions...")
    best_interaction = 0
    for i, lig in enumerate(tqdm(optimized_ligs)):
        interaction_energy = sqm.calculate_interaction(
            lig, lig.GetProp("target_path"), f"{i}_{lig.GetProp('_Name')}"
        )
        lig.SetDoubleProp("STERM_INTERACTION", interaction_energy)
        best_interaction = min(best_interaction, interaction_energy)
    print(f"Done! Best Interaction: {best_interaction:.2f}")
    tock_81 = time.time()
    print(f"Time Spent: {tock_81-tock_7:.2f}s")
    print(f"Time Spent Per Original Lig: {(tock_81-tock_7)/n_original_ligs:.2f}")

    ###############################################################################
    
    print("2) BASE STATE ENERGY")
    print("Estimating Baseline Tautomers, Protomers & Their Ratios... (Not Really...)")
    with contextlib.redirect_stdout(third_party_out):
        baselines_list = []
        for i, lig in enumerate(tqdm(optimized_ligs)):
            baselines = baseline.get_molecule_baselines(lig)
            baselines_list.append(baselines)
    print("Optimizing & Calculating Baseline Energies...")
    best_opt_energy = 0
    final_baselines = []
    for i, baselines in enumerate(tqdm(baselines_list)):
        bline = baselines[0] # TODO: ????
        final_baselines.append(bline)
        output = sqm.optimize_molecule_qm(bline, f"{i}_{lig.GetProp('_Name')}")
        opt_energy = sqm.calculate_energy(output["mol_rdkit"], f"{i}_{lig.GetProp('_Name')}")
        optimized_ligs[i].SetDoubleProp("STERM_BASELINE_ENERGY", opt_energy)
        best_opt_energy = min(best_opt_energy, opt_energy)
    print(f"Done! Best Baseline Energy: {best_opt_energy:.2f}")
    tock_82 = time.time()
    print(f"Time Spent: {tock_82-tock_81:.2f}s")
    print(f"Time Spent Per Original Lig: {(tock_82-tock_81)/n_original_ligs:.2f}")

    ###############################################################################

    print("3) BOUND STATE ENERGY")
    print("Calculating Bound State Isolated Energy Of Ligands...")
    best_bound_energy = 0
    for i, lig in enumerate(tqdm(optimized_ligs)):
        energy = sqm.calculate_energy(lig, f"{i}_{lig.GetProp('_Name')}")
        lig.SetDoubleProp("STERM_BOUND_ENERGY", energy)
        best_bound_energy = min(best_bound_energy, energy)
    print(f"Done! Best Bound Energy: {best_bound_energy:.2f}")
    tock_83 = time.time()
    print(f"Time Spent: {tock_83-tock_82:.2f}s")
    print(f"Time Spent Per Original Lig: {(tock_83-tock_82)/n_original_ligs:.2f}")

    ###############################################################################

    print("4) PROTON TRANSFER ENERGIES FOR PROTOMERS")
    hscored_ligs_path = os.path.join(OUT_DIR, f"hscored_ligands.sdf")
    baseline_ligs_path = os.path.join(OUT_DIR, f"baselines.sdf")
    writer_hscored = Chem.SDWriter(hscored_ligs_path)
    writer_baseline = Chem.SDWriter(baseline_ligs_path)
    for lig, bline in zip(optimized_ligs, final_baselines):
        lig_charge = Chem.GetFormalCharge(lig)
        bline_charge = Chem.GetFormalCharge(bline)
        num_hs_added = lig_charge - bline_charge
        htransfer_energy = sqm.calc_hydrogen_transfer_energy(num_hs_added)
        lig.SetDoubleProp("STERM_HTRANSFER_ENERGY", htransfer_energy)
        writer_hscored.write(lig)
        writer_baseline.write(bline)
    writer_baseline.close()
    writer_hscored.close()
    tock_84 = time.time()
    print(f"Done! Baselines Saved To {baseline_ligs_path}, Enthalpy Scored Ligs To {hscored_ligs_path}.")
    
    ###############################################################################

    print("5) ENTROPY & ENTROPIC CONTRIBUTION")
    model = joblib.load(ENTROPY_MODEL)
    df = PandasTools.LoadSDF(hscored_ligs_path, removeHs=False, molColName="Mol")
    df["EntropyFeatures"] = df.apply(calculate_features)
    entropies = np.maximum(model.predict(np.log(df["EntropyFeatures"].values+1)), 0)
    df["PredEntropy"] = pd.Series(entropies)
    df["STERM_ENTROPIC_CONTRIB"] = df["PredEntropy"] * 298
    print("Done!")

    df["BINDING_FREE_ENERGY"] = df["STERM_INTERACTION"] + (df["STERM_BOUND_ENERGY"] - df["STERM_BASELINE_ENERGY"]) + \
        df["STERM_HTRANSFER_ENERGY"] + df["STERM_ENTROPIC_CONTRIB"]
    final_file = os.path.join(OUT_DIR, "final_ligands.sdf")
    PandasTools.WriteSDF(df, final_file, idName="ID", molColName="Mol", properties=list(df.columns))




    # 6) Rescore using MOPAC (with Cuby4) energy job with cosmo2 solvation model
    # scores = []
    # for i, lig in tqdm(enumerate(redocked_ligs)):
    #     lig = Chem.AddHs(lig, addCoords=True)
    #     lig_poc_loc = PocketLocation("ligand", ligand=lig, radius=10)
    #     lig_poc_path = os.path.join(TMP_DIR, f"{i}_{lig.GetProp('_Name')}_pocket.pdb")
    #     isolate_pocket(prepped_protein_path, lig_poc_loc, lig_poc_path)
    #     score = sqm.calculate_interaction(lig, lig_poc_path, job_name=f"rescore_{i}_{lig.GetProp('_Name')}")
    #     scores.append(score)

    # print(scores)

if __name__ == "__main__":
    main()