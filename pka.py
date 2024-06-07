from rdkit import Chem
from rdkit.Chem import PandasTools
import re
import subprocess
from collections import defaultdict
import math
from rdkit.Chem import AllChem


pka_list = []
p_list = []
protonation_list = []
checked_smiles = []

def qupkake(smiles, output = "qupkake_output.sdf"):
    subprocess.run(['qupkake', 'smiles', smiles, output ])
    checked_smiles.append(smiles)

def dict_pka(sdf_file = "./data/output/qupkake_output.sdf", energy_first =0.0):
    with open(sdf_file, 'r') as f:
        sdf_data = f.read()
    suppl = Chem.SDMolSupplier()
    suppl.SetData(sdf_data)
    molecule_data = {}
    for mol in suppl:
        atom_idx = int(mol.GetProp('idx'))
        atom = mol.GetAtomWithIdx(atom_idx)
        formal_charge = atom.GetFormalCharge()
        num_hydrogens = atom.GetNumExplicitHs()
        pka_type = mol.GetProp('pka_type')
        if pka_type == "acidic":
            atom.SetFormalCharge(formal_charge - 1)
            atom.SetNumExplicitHs(num_hydrogens - 1)
        else:
            atom.SetFormalCharge(formal_charge + 1)
            atom.SetNumExplicitHs(num_hydrogens + 1)
        smiles = Chem.MolToSmiles(mol)
        pka_number = float(mol.GetProp('pka'))
        energy = dict_energy(pka_type, pka_number ,energy_first= energy_first)
        bol = energy_bol(energy)
        if pka_range_tester(pka_number):
            molecule_data[smiles] = [pka_type, pka_number, energy, bol]
    molecule_data = unique(molecule_data)
    print(molecule_data)
    return dict(molecule_data)


def pka_range_tester(pka, ph_min = 5, ph_max = 9,sdf_file = "./data/output/qupkake_output.sdf"):
    if ph_min < pka < ph_max:
        return True
    else:
        return False

def unique(dict):
    unique_smiles = defaultdict(list)
    for smiles, properties in dict.items():
        unique_smiles[smiles].append(properties)
    return {smiles: properties for smiles, properties in unique_smiles.items()}



def dict_energy(type_pka, pka, energy_first = 0.0, PH = 7.0):
    if type_pka == "acidic" :
        energy = energy_first - 2.303 * 2.479 * (float(pka) - PH)
    else:
        energy = 2.303 * 2.479 * (float(pka) - PH) + energy_first
    return energy

def dict_update(dict1, dict2):
    return dict(dict1).update(dict(dict2))

def energy_bol(energy):
    energy_bol = math.exp((-1) / 2.479 * energy )
    return energy_bol

def probability(dict, smile):
    bols = []
    probability = []
    for smiles in dict:
        bol = dict[smiles][0][-1]
        bols.append(bol)
    sum_bols = sum(bols)
    for smiles in dict :
        bol = dict[smiles][0][-1]
        probability.append(bol / sum_bols)
    bol = dict[smile][0][-1]
    our_probability = bol / sum_bols
    return probability, our_probability

def iterate(my_dict :dict,checked_smiles = checked_smiles):
    so_dict = my_dict
    for smiles in my_dict:
        energy_first = my_dict[smiles][0][2]
        if smiles not in checked_smiles:
            qupkake(smiles, output="so.sdf")
            second_dict = dict_pka(sdf_file= "./data/output/so.sdf", energy_first = energy_first)
            so_dict = dict_update(so_dict, second_dict)
    return so_dict

def probability_finder(smiles):
    your_smile = smiles
    qupkake(your_smile)
    my_dict = dict_pka()
    print(my_dict)
    while len(my_dict) > len(checked_smiles):
        print("-----------------------------------------------------------------------------------------------------------------------------------------------",my_dict)
        my_dict = iterate(my_dict)
    for smiles in my_dict:
        print(smiles ,probability(my_dict, smiles)[1])
        

probability_finder(("C1(=C(CO)N=C(C(=C1)N(C)C)N)C(=O)O"))