from time import time
from typing import List, Tuple

from rdkit import Chem

from src.utils.iotools import load_ligands, write_ligands
from src.config import *

protein_path = INPUT_PROTEIN
ligands = load_ligands(INPUT_LIGANDS, final_3d=False, addHs=False)

if PROTEIN_PREPARATION_PERFORM:
    print("#"*30)
    print("##" + str.center("PROTEIN PREPARATION", 25), "##")
    print("#"*30)
    tic = time()
    protein_path = PROTEIN_PREPPER.prepare(protein_path)
    toc = time()
    print(f"Done! Saved Output: {protein_path}")
    print(f"Time Spent: {toc-tic:.2f}s")
    print()

#######################################################################

if LIGANDS_PREPARATION_PERFORM:
    print("#"*30)
    print("##" + str.center("LIGAND PREP & ENUM", 25), "##")
    print("#"*30)
    tic = time()
    n_input_ligs = len(ligands)
    ligands = LIGAND_PREPPER.prepare(ligands)
    out_path = os.path.join(OUT_DIR, f"{PRE_DOCKING_LIGANDS_PREFIX}_ligands.sdf")
    with Chem.SDWriter(out_path) as writer:
        [writer.write(ligand) for ligand in ligands]
    toc = time()
    print(toc-tic)
    print(f"Done! Number Of Prepared Ligands: {len(ligands)}")
    print(f"Time Spent: {toc-tic:.2f}s")
    print(f"Time Spent Per Input Lig: {(toc-tic)/n_input_ligs:.2f}")
    print()

#######################################################################

if DOCKING_PERFORM:
    print("#"*30)
    print("##" + str.center("VIRTUAL SCREENING", 25), "##")
    print("#"*30)
    print("[Performing Virtual Screening...]")
    tic = time()
    n_input_ligs = len(ligands)
    out_path = os.path.join(OUT_DIR, f"{POST_DOCK_LIGANDS_PREFIX}_ligands.sdf")
    ligands = DOCKING_PERFORMER.dock(ligands, protein_path, out_path, PROTEIN_POCKET_LOCATION,
                                    return_mols=True)
    print(f"Done! Number Of Resulting Conformations: {len(ligands)}")
    print("[Selecting Poses...]")
    n_select = int(n_input_ligs * DOCKING_SELECT_RATIO) + 1
    ligands = DOCKING_SELECTOR.select(ligands, n_select)
    print(f"Done! Number Of Selected Conformations: {len(ligands)}")
    toc = time()
    print(f"Time Spent: {toc-tic:.2f}s")
    print(f"Time Spent Per Input Lig: {(toc-tic)/n_input_ligs:.2f}")
    print()

#######################################################################

if PROTOMER_GENERATION_PERFORM:
    print("#"*30)
    print("##" + str.center("PROTOMER ENUMERATION", 25), "##")
    print("#"*30)
    tic = time()
    n_input_ligs = len(ligands)
    min_ph = PROTOMER_GENERATION_PH - PROTOMER_GENERATION_PH_RANGE
    max_ph = PROTOMER_GENERATION_PH + PROTOMER_GENERATION_PH_RANGE
    print(f"pH Range: {min_ph}-{max_ph}")
    ligands = PROTOMER_GENERATION_ENUMERATOR.enumerate_multi(ligands, min_ph=min_ph, max_ph=max_ph)
    out_path = os.path.join(OUT_DIR, f"{PROTOMERS_SELECTED_PREFIX}_ligands.sdf")
    write_ligands(ligands, out_path)
    toc = time()
    print(f"Done! Number Of Resulting Protomers: {len(ligands)}")
    print(f"Time Spent: {toc-tic:.2f}s")
    print(f"Time Spent Per Input Lig: {(toc-tic)/n_input_ligs:.2f}")
    print()

#######################################################################

# TODO: This part needs some reordering
if MM_RESCORING_PERFORM:
    print("#"*30)
    print("##" + str.center("MM MINIMIZATION", 25), "##")
    print("#"*30)
    print("[Rescoring Ligands-Protein Complexes...]")
    tic = time()
    n_input_ligs = len(ligands)
    scores, cs_list = [], []
    for ligand in ligands: 
        score, cs = MM_RESCORER.score(ligand, protein_path, return_used_components=True)
        scores.append(score)
        cs_list.append(cs)
    if DEBUG_MODE:
        print("rescores", scores)
    # scores, cs_list = zip(*[MM_RESCORER.score(ligand, protein_path, return_used_components=True)
    #                         for ligand in ligands])
    print("[Selecting From Scores Complexes...]")
    sel_scores, sel_ligs, sel_targs = [], [], []
    n_select = int(n_input_ligs * MM_SELECT_RATIO) + 1
    sel_counter = 0
    name_counter = dict()
    for score, (lig_component, targ_component), ligand in sorted(
        zip(scores, cs_list, ligands), key=lambda x: x[0]):

        if sel_counter > n_select:
            break

        if score > MM_SELECT_MAX_ENERGY:
            continue
        
        if ligand.GetProp("_Name") not in name_counter.keys():
            name_counter[ligand.GetProp("_Name")] = 0
        if name_counter[ligand.GetProp("_Name")] >= MM_SELECT_MAX_PER_LIG:
            continue
        
        sel_counter += 1
        sel_scores.append(score)
        newlig = lig_component.to_rdkit()
        newlig.SetProp("_Name", ligand.GetProp("_Name"))
        targ_component.name = f"{MM_RESCORE_PREFIX}_prot{sel_counter}_{ligand.GetProp('_Name')}"
        targ_path = targ_component.to_pdb()
        sel_targs.append(targ_path)
        newlig.SetProp("CorrespondingTargPath", targ_path)
        newlig.SetDoubleProp("MMRescore", score)
        sel_ligs.append(newlig)
        
    ligands = sel_ligs
    protein_paths = sel_targs
    out_path = os.path.join(OUT_DIR, f"{MM_RESCORE_PREFIX}_ligands.sdf")
    write_ligands(ligands, out_path)
    toc = time()
    print(f"Done! Number Of Selected Rescored Ligands: {len(ligands)}")
    print(f"Time Spent: {toc-tic:.2f}s")
    print(f"Time Spent Per Input Lig: {(toc-tic)/n_input_ligs:.2f}")
    print()

#######################################################################

print("#"*30)
print("##" + str.center("SQM OPTIMIZATION", 25), "##")
print("#"*30)


    