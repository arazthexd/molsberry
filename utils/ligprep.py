from typing import List
from itertools import product

from rdkit import Chem
from rdkit.Chem import (
    rdDistGeom, 
    rdMolTransforms, 
    rdMolAlign, 
    AllChem, 
    rdForceFieldHelpers, 
    rdmolops,
    rdRascalMCES,
    rdFMCS
)
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.ML.Cluster import Butina

import dimorphite_dl as dd

ROTOR_QUERY = Chem.MolFromSmarts("[!$(*#*)&!D1]-!@[!$(*#*)&!D1]")

def get_rotatable_dihedrals(mol):
    dihs = []
    for atom2, atom3 in mol.GetSubstructMatches(ROTOR_QUERY):
        neis2 = mol.GetAtomWithIdx(atom2).GetNeighbors()
        if neis2[0].GetIdx() != atom3:
            atom1 = neis2[0].GetIdx()
        else:
            atom1 = neis2[1].GetIdx()
        neis3 = mol.GetAtomWithIdx(atom3).GetNeighbors()
        if neis3[0].GetIdx() != atom2:
            atom4 = neis3[0].GetIdx()
        else:
            atom4 = neis3[1].GetIdx()
        dihs.append((atom1, atom2, atom3, atom4))
    return dihs

def sync_mol_flexible_rotors(mol_prb: Chem.Mol, mol_ref: Chem.Mol):
    
    rot_dihs_prb = get_rotatable_dihedrals(mol_prb)
    rot_dihs_ref = get_rotatable_dihedrals(mol_ref)
    mcs = rdFMCS.FindMCS([mol_prb, mol_ref])
    match_prb: tuple = mol_prb.GetSubstructMatch(mcs.queryMol)
    match_ref: tuple = mol_ref.GetSubstructMatch(mcs.queryMol)
    for prb_dih in rot_dihs_prb:
        if not all(a in match_prb for a in prb_dih):
            continue
        ref_dih = [match_ref[match_prb.index(a)] for a in prb_dih]
        deg = rdMolTransforms.GetDihedralDeg(mol_ref.GetConformer(), *ref_dih)
        rdMolTransforms.SetDihedralDeg(mol_prb.GetConformer(), *prb_dih, deg)
    return rdMolAlign.AlignMol(mol_prb, mol_ref, atomMap=list(zip(match_prb, match_ref)))

def sync_flexible_rotors(mol: Chem.Mol, conf_id: int, conf_ref: int, rot_dihs: list):
    
    for dih in rot_dihs:
        deg = rdMolTransforms.GetDihedralDeg(mol.GetConformer(conf_ref), *dih)
        rdMolTransforms.SetDihedralDeg(mol.GetConformer(conf_id), *dih, deg)
    return rdMolAlign.AlignMolConformers(mol)

def enumerate_stereoisomers(mol: Chem.Mol):
    return list(EnumerateStereoisomers(mol))

def enumerate_tautomers(mol: Chem.Mol, max_tautomers: int = 4):
    te = rdMolStandardize.TautomerEnumerator()
    te.SetRemoveBondStereo(False)
    te.SetRemoveSp3Stereo(False)
    ts: List[Chem.Mol] = list(te.Enumerate(mol))
    # [t.SetProp("TInput", Chem.MolToSmiles(mol)) for t in ts]
    ts = sorted(ts, key=lambda t: te.ScoreTautomer(t), reverse=True)
    return ts[:max_tautomers]

def enumerate_protomers(mol: Chem.Mol, min_ph=7.2, max_ph=7.6, **kwargs):
    
    assert isinstance(mol, Chem.Mol)
    mols = [Chem.AddHs(mol)]

    protomers = dd.run_with_mol_list(mols, silent=True, min_ph=min_ph, max_ph=max_ph, **kwargs)
    protomers = [Chem.AddHs(Chem.RemoveAllHs(protomer)) for protomer in protomers]
    [rdDistGeom.EmbedMolecule(protomer) for protomer in protomers]
    [sync_mol_flexible_rotors(protomer, mol) for protomer in protomers]
    [p.SetProp("_Name", mol.GetProp("_Name")) for p in protomers]
    return protomers

def enumerate_ring_confs(mol: Chem.Mol, max_per_ring: int = 4, num_conf_gen: int = 30, minimize: bool = False,
                         addHs: bool = False, dist_threshold: float = 0.5):

    mol = Chem.AddHs(mol)
    rdDistGeom.EmbedMultipleConfs(mol, numConfs=num_conf_gen)
    if minimize:
        rdForceFieldHelpers.MMFFOptimizeMoleculeConfs(mol)
    
    rdmolops.FindRingFamilies(mol)
    rinfo = mol.GetRingInfo()
    clusts_per_ring = []
    for ring_atoms in rinfo.AtomRings():
        dists = AllChem.GetConformerRMSMatrix(mol, prealigned=False, atomIds=ring_atoms)
        clusts = Butina.ClusterData(dists, mol.GetNumConformers(), dist_threshold, isDistData=True, 
                                    reordering=True)
        clusts = [set(clust) for clust in clusts]
        clusts_per_ring.append(clusts[:max_per_ring])
    
    rings_combos = list(product(*clusts_per_ring))
    final_mols = []
    # print(rings_combos)
    for combo in rings_combos:
        if len(combo) == 0: 
            continue
        common: set = combo[0]
        for s in combo[1:]:
            common.intersection_update(s)
        # print(common)
        if common:
            final_mols.append(Chem.Mol(mol, confId=common.pop()))
    
    # print(final_mols)
    return final_mols
    
    



# import dimorphite_dl as dd
# m = Chem.MolFromSmiles("COC1=C(OC2=C(OCCO)N=C(N=C2NS(=O)(=O)C2=CC=C(C=C2)C(C)(C)C)C2=NC=CC=N2)C=CC=C1")
# ms = enumerate_stereoisomers(m)
# ms = [enumerate_tautomers(m) for m in ms]
# ms = dd.run_with_mol_list([m for mm in ms for m in mm], min_ph=4.4, max_ph=10.4)
# print([Chem.MolToSmiles(m) for m in ms])
# dd.run
# enumerate_protomers(Chem.MolFromSmiles("COC1=C(OC2=C(OCCO)N=C(N=C2NS(=O)(=O)C2=CC=C(C=C2)C(C)(C)C)C2=NC=CC=N2)C=CC=C1"))
# enumerate_ring_confs(Chem.MolFromSmiles("C[C@H]1O[C@H](O[C@@H]2[C@@H](CO)O[C@H](O[C@@H]3[C@@H](CO)OC(O)[C@H](O)[C@H]3O)[C@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@@H]1N[C@H]1C=C(CO)[C@@H](O)[C@H](O)[C@H]1O"), 
#                      num_conf_gen=60, minimize=True)