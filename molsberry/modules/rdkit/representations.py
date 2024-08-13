from __future__ import annotations

from typing import List
from numpy import ndarray

from molsberry.core.data.abstract import Representation

from ...core import Molecule3DRep, SmallMolRep, ProteinRep
from ...core import SMILESRep, PDBPathRep
from .interface import RDKitInterface

from rdkit import Chem

class RDKitMolRep(Molecule3DRep):
    rep_name = "rdmol"

    def __init__(self, mol: Chem.Mol):
        assert isinstance(mol, Chem.Mol)
        super().__init__(content=mol)
        self.content: Chem.Mol
    
    @classmethod
    def from_SMILESRep(cls, smiles_rep: SMILESRep):
        assert isinstance(smiles_rep, SMILESRep)
        smiles = smiles_rep.content
        mol = Chem.MolFromSmiles(smiles)
        return cls(mol=mol)
    
    @classmethod
    def from_PDBPathRep(cls, pdb_rep: PDBPathRep):
        assert isinstance(pdb_rep, PDBPathRep)
        pdb_path = pdb_rep.content
        mol = Chem.MolFromPDBFile(pdb_path, removeHs=False, sanitize=True)

        # Fix COO groups not being ionized when read from pdb in rdkit
        # TODO: Probably needs more attention
        query_COO = Chem.MolFromSmarts("[$([O]-C(=O)-C)]")
        for atom, in mol.GetSubstructMatches(query_COO):
            atom: Chem.Atom = mol.GetAtomWithIdx(atom)
            if atom.GetPDBResidueInfo().GetIsHeteroAtom() == False:
                atom.SetFormalCharge(-1)
        return cls(mol=mol)
    
    @classmethod
    def from_RDKitMolRep(cls, rdmol_rep: RDKitMolRep):
        return cls(rdmol_rep.content)
    
    def to_RDKitMolRep(self) -> RDKitMolRep:
        return RDKitMolRep(self.content)
    
    def update_coordinates(self, coords: ndarray) -> None:
        rdmol = self.content
        if RDKitInterface.is_mol_3d(rdmol):
            conformer = rdmol.GetConformer()
            [conformer.SetAtomPosition(i, loc) for i, loc in enumerate(coords)]
        else:
            conformer = Chem.Conformer(rdmol.GetNumAtoms())
            conformer.Set3D(True)
            [conformer.SetAtomPosition(i, loc) for i, loc in enumerate(coords)]
            rdmol.AddConformer(conformer)

class RDKitSmallMolRep(RDKitMolRep, SmallMolRep):
    def save_rep(self, exless_filename: str):
        rep_path = exless_filename + ".sdf"
        writer = Chem.SDWriter(rep_path)
        writer.write(self.content)
        writer.close()
    
    @classmethod
    def save_rep_batch(cls, reps: List[Representation], exless_filename: str):
        rep_path = exless_filename + ".sdf"
        writer = Chem.SDWriter(rep_path)
        for rep in reps:
            writer.write(rep.content)
        writer.close()

class RDKitProtRep(RDKitMolRep, ProteinRep):
    rep_name = "rdprot"
    pass

        
        
