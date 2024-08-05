from rdkit import Chem
from rdkit.Chem import (
    rdFMCS,
    rdMolAlign,
    rdMolTransforms
)

from ...core.data.abstract import SpecialDataClass
# from .representations import RDKitMolRep

ROTOR_QUERY = Chem.MolFromSmarts("[!$(*#*)&!D1]-!@[!$(*#*)&!D1]")

class RDKitInterface:

    @staticmethod
    def is_mol_3d(mol: Chem.Mol):
        if mol.GetNumConformers() > 0:
            if mol.GetConformer().Is3D():
                return True
        return False
    
    @staticmethod
    def special_cls_to_rdmol(special_cls: SpecialDataClass):
        rdmol = Chem.Mol(special_cls.get_data("rdkit_mol"))
        return rdmol
    
    @staticmethod
    def addhs_based_on_confdim(lig: Chem.Mol):
        if RDKitInterface.is_mol_3d(lig):
            return Chem.AddHs(lig, addCoords=True)
        else:
            return Chem.AddHs(lig)

    @staticmethod
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

    @staticmethod
    def sync_mol_flexible_rotors(mol_prb: Chem.Mol, mol_ref: Chem.Mol):
        
        rot_dihs_prb = RDKitInterface.get_rotatable_dihedrals(mol_prb)
        rot_dihs_ref = RDKitInterface.get_rotatable_dihedrals(mol_ref)
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

    @staticmethod
    def sync_flexible_rotors(mol: Chem.Mol, conf_id: int, conf_ref: int, rot_dihs: list):
        
        for dih in rot_dihs:
            deg = rdMolTransforms.GetDihedralDeg(mol.GetConformer(conf_ref), *dih)
            rdMolTransforms.SetDihedralDeg(mol.GetConformer(conf_id), *dih, deg)
        return rdMolAlign.AlignMolConformers(mol)