import os, subprocess, glob, shutil, pathlib
from multiprocessing import Pool
import math
from typing import List, Dict, Tuple

from rdkit import Chem
from rdkit.Chem import rdDistGeom, rdForceFieldHelpers, rdMolAlign
from rdkit.Chem.MolStandardize import rdMolStandardize

from utils.ligprep import enumerate_stereoisomers, sync_mol_flexible_rotors

# CONFIG & VARIABLES
WORK_DIR = str(pathlib.Path("tmp").absolute())
print(WORK_DIR)
BASELINE_PROTOMER_PH = 7.4
BASELINE_PROTOMER_PH_THRESHOLD = 5.0
N_THREADS = 4

# CONSTANTS
LN_TO_LOG = 2.303
BETA = 1 / 2.479 # kJ/mol

def get_baseline_tautomers(mol: Chem.Mol, mode: str = "Prop:InitialSmiles"):
    
    if mode.split(":")[0] == "Prop":
        baseline_taus = _get_baseline_tautomer_prop(mol, mode.split(":")[1])
        return {tau: 1/len(baseline_taus) for tau in baseline_taus}
    else:
        raise NotImplementedError()

def _get_baseline_tautomer_prop(mol: Chem.Mol, prop: str): 

    tau = Chem.MolFromSmiles(mol.GetProp(prop))
    tau = Chem.AddHs(tau)
    tau_stereoisos = enumerate_stereoisomers(tau)
    for tau in tau_stereoisos:
        rdDistGeom.EmbedMolecule(tau)
        rdForceFieldHelpers.MMFFOptimizeMolecule(tau)
        sync_mol_flexible_rotors(tau, mol)
    return tau_stereoisos

def calc_boltzmann_factor(energy: float):
    return math.exp(-BETA * energy)

def pka_to_energy(pka, pka_type, initial_energy=0.0, ph=BASELINE_PROTOMER_PH):

    if pka_type == "acidic" :
        # energy = initial_energy - 2.303 * 2.479 * (float(pka) - ph)
        energy = (LN_TO_LOG / BETA) * (float(pka) - ph) + initial_energy
    else:
        # energy = 2.303 * 2.479 * (float(pka) - ph) + initial_energy
        energy = initial_energy - (LN_TO_LOG / BETA) * (float(pka) - ph)
    return energy

def run_qupkake(args, tautomerize=False):

    mol = args[0]
    mol = Chem.RemoveHs(mol)
    prop_dict: dict = args[1]
    name = prop_dict["_Name"]
    tmp_idx = prop_dict["_tmp_idx"]
    [mol.SetProp(k, v) for k, v in prop_dict.items()]

    root = os.path.join(WORK_DIR, f"{name}_{tmp_idx}")
    try:
        shutil.rmtree(root)
    except:
        pass
    if not os.path.exists(root):
        if not os.path.exists(os.path.join(WORK_DIR, name)):
            os.mkdir(os.path.join(WORK_DIR, name))
        os.mkdir(root)

    input_path = os.path.join(root, "qupkake_input.sdf")
    print(input_path)
    writer = Chem.SDWriter(input_path)
    writer.write(mol)
    writer.close()

    # try:
    #     shutil.rmtree(os.path.join(WORK_DIR, "output"))
    #     shutil.rmtree(os.path.join(WORK_DIR, "raw"))
    #     shutil.rmtree(os.path.join(WORK_DIR, "processed"))
    #     shutil.rmtree(os.path.join(WORK_DIR, "logs"))
    # except: 
    #     pass

    tautomer_args = []
    if tautomerize: tautomer_args.append("-t")
    cmd = ['qupkake', 'file', input_path, "-r", root,
                             "-o", "qupkake_out.sdf", "-mp", str(N_THREADS)] \
                                + tautomer_args
    print(" ".join(cmd))
    _ = subprocess.run(['qupkake', 'file', input_path, "-r", root,
                             "-o", "qupkake_out.sdf", "-mp", str(N_THREADS)] \
                                + tautomer_args, capture_output=False)
    try:
        qupkake_results = [mol for mol in Chem.SDMolSupplier(
            os.path.join(root, "output", "qupkake_out.sdf")
        )]
    except FileNotFoundError:
        return []
    return qupkake_results

def runrun_qupkake(pa):
    return run_qupkake(pa)

def process_qupkake(qupkake_results: List[Chem.Mol]):

    atom2pka: Dict[int, float] = dict() # float: pka
    atom2mol: Dict[int, Chem.Mol] = dict()

    for mol in qupkake_results:
        mol = Chem.RemoveHs(mol)

        pka = float(mol.GetProp('pka')) # [7:-1] if: tensor(5.0304) >> 5.0304
        # if abs(BASELINE_PROTOMER_PH - pka) > BASELINE_PROTOMER_PH_THRESHOLD:
        #     continue

        atom_idx = int(mol.GetProp('idx'))
        if atom_idx in atom2pka.keys():
            if abs(BASELINE_PROTOMER_PH - pka) > abs(BASELINE_PROTOMER_PH - atom2pka[atom_idx]):
                continue
        
        atom = mol.GetAtomWithIdx(atom_idx)
        formal_charge = atom.GetFormalCharge()
        num_hydrogens = atom.GetNumExplicitHs()

        pka_type = mol.GetProp('pka_type')
        if pka_type == "acidic":
            atom.SetFormalCharge(formal_charge - 1)
            atom.SetNumExplicitHs(max(num_hydrogens - 1, 0))
            atom.UpdatePropertyCache(False)
        elif pka_type == "basic":
            atom.SetFormalCharge(formal_charge + 1)
            atom.UpdatePropertyCache(False)

        Chem.SanitizeMol(mol)
        mol.UpdatePropertyCache()
        atom2pka[atom_idx] = pka
        atom2mol[atom_idx] = mol

    return list(atom2mol.values())

def get_baseline_protomers(mol: Chem.Mol): # Must have Hs

    # TODO: Figure a way to also use the tautomer selection?
    # init_baseline = run_qupkake(mol, tautomerize=False)[0]
    # init_baseline.SetDoubleProp("pKaEnergy", 0.0)
    # init_baseline_smi = Chem.MolToSmiles(init_baseline)
    
    checked_smi2data = dict()

    mol = Chem.RemoveHs(mol)
    uncharger = rdMolStandardize.Uncharger()
    uncharger.unchargeInPlace(mol)
    prop_dict = mol.GetPropsAsDict(includePrivate=True)
    if "_Name" not in prop_dict.keys():
        prop_dict["_Name"] = "unnamed"

    tmp_idx = 0
    prop_dict["_tmp_idx"] = str(tmp_idx)
    tmp_idx += 1

    qresults = run_qupkake((mol, prop_dict), tautomerize=True)
    mol = Chem.Mol(qresults[0])
    [mol.SetProp(k, v) for k, v in prop_dict.items()]
    mol.SetDoubleProp("pKaEnergy", 0.0)
    print(mol.GetPropsAsDict(includeComputed=True, includePrivate=True))
    smi = Chem.MolToSmiles(mol)
    checked_smi2data[smi] = {
        "num": 1,
        "energy": 0.0,
        "mol": mol
    }
    protomers_tocheck = [mol]
    
    pool = Pool(N_THREADS)
    energy_min = 0
    i = 0
    while True:
        i += 1
        print(f"Creating New Protomers... (Interation {i})")
        if i > 1:
            [p.SetProp("_tmp_idx", str(tmp_idx+i)) for i, p in enumerate(protomers_tocheck)]
            tmp_idx += len(protomers_tocheck)
            prop_dicts = [p.GetPropsAsDict(includePrivate=True) for p in protomers_tocheck]
            qresults_list = pool.map(runrun_qupkake, list(zip(protomers_tocheck, prop_dicts)))
            print(qresults_list[0].GetPropsAsDict())
        else:
            qresults_list = [qresults]
        for protomer, qresults in zip(protomers_tocheck, qresults_list):
            new_protomers = process_qupkake(qresults)
            print("Parent Protomer: ", Chem.MolToSmiles(protomer))
            print("Child Protomers: ")
            print("\t", "\n\t".join([Chem.MolToSmiles(p) for p in new_protomers]), sep="")

            nprotomers_to_add = []
            for nprotomer in new_protomers:
                smi = Chem.MolToSmiles(nprotomer)
                pka = float(nprotomer.GetProp("pka"))
                pka_type = nprotomer.GetProp("pka_type")
                energy = pka_to_energy(pka, pka_type, 
                                       checked_smi2data[Chem.MolToSmiles(protomer)]["energy"])
                
                # if energy - energy_min > 10.0:
                #     print(f"SMILES {smi} Energy Too High, Discarded")
                #     continue
                if energy - checked_smi2data[Chem.MolToSmiles(protomer)]["energy"] > 8.0:
                    print(f"SMILES {smi} Energy Too High, Discarded (E: {energy}, baseline: {checked_smi2data[Chem.MolToSmiles(protomer)]['energy']})")
                    continue
                if smi in checked_smi2data.keys():
                    e_old = checked_smi2data[smi]["energy"]
                    e_new = energy
                    energy = ((checked_smi2data[smi]["num"] * checked_smi2data[smi]["energy"]) + energy) \
                        / (checked_smi2data[smi]["num"] + 1)
                    checked_smi2data[smi]["num"] += 1
                    checked_smi2data[smi]["energy"] = energy
                    print(f"SMILES {smi} Already Checked >> Updated E: {e_old:.2f}(old), {e_new:.2f}(new) to {energy:.2f}")
                    continue
                if abs(Chem.GetFormalCharge(nprotomer)) >= 2:
                    print(f"SMILES {smi} Has Absolute Charge > 1, Discarded")
                    continue
                
                print(f"SMILES {smi} Added To New Protomers To Check")
                
                # sync_mol_flexible_rotors(nprotomer, protomer)

                checked_smi2data[smi] = {
                    "num": 1,
                    "energy": energy,
                    "mol": Chem.Mol(nprotomer)
                }
                energy_min = min([data["energy"] for data in checked_smi2data.values()])
                nprotomers_to_add.append(Chem.Mol(nprotomer))
        
        print(f"Number of New Protomers To Check: {len(nprotomers_to_add)}")
        print()
        if len(nprotomers_to_add) == 0:
            break
        protomers_tocheck = list(nprotomers_to_add)
    
    energy_sum = sum([data["energy"] for data in checked_smi2data.values()])
    energy_mean = energy_sum / len(checked_smi2data)
    protomers = []
    for data in checked_smi2data.values():
        energy = data["energy"]
        bfactor = calc_boltzmann_factor(energy - energy_mean)
        protomer = data["mol"]
        protomer.SetDoubleProp("BoltzmannFactor", bfactor)
        protomers.append(protomer)
    
    partition = sum([p.GetDoubleProp("BoltzmannFactor") for p in protomers])
    protomers = [p for p in protomers if 
                 p.GetDoubleProp("BoltzmannFactor")/partition > 0.05]
    partition = sum([p.GetDoubleProp("BoltzmannFactor") for p in protomers])
    for protomer in protomers:
        bfactor = protomer.GetDoubleProp("BoltzmannFactor")
        protomer.SetDoubleProp("Probability", bfactor/partition)

    print([Chem.MolToSmiles(mm) for mm in protomers])
    return protomers

m = Chem.MolFromSmiles("NC(=N)NC1=NC(CSCCC(=N)NS(N)(=O)=O)=CS1")
m = Chem.AddHs(m)
rdDistGeom.EmbedMolecule(m)
Chem.AssignStereochemistryFrom3D(m)
writer = Chem.SDWriter("lulu.sdf")
m.SetProp("_Name", "thelulu")
for p in get_baseline_protomers(m):
    writer.write(p)