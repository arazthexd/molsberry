from numpy import ndarray

from ...core.data.abstract import Representation
from ...core.data.representations import SMILESRep, PDBPathRep
from .interface import RDKitInterface

from rdkit import Chem

class RDKitMolRep(Representation):
    rep_name = "rdkit_mol"

    def __init__(self, mol: Chem.Mol):
        assert isinstance(mol, Chem.Mol)
        super().__init__(data=mol)
        self.data: Chem.Mol
    
    @classmethod
    def from_SMILESRep(cls, smiles_rep: SMILESRep):
        assert isinstance(smiles_rep, SMILESRep)
        smiles = smiles_rep.data
        mol = Chem.MolFromSmiles(smiles)
        return cls(mol=mol)
    
    @classmethod
    def from_PDBPathRep(cls, pdb_rep: PDBPathRep):
        assert isinstance(pdb_rep, PDBPathRep)
        pdb_path = pdb_rep.data
        mol = Chem.MolFromPDBFile(pdb_path, removeHs=False, sanitize=False)

        # Fix COO groups not being ionized when read from pdb in rdkit
        # TODO: Probably needs more attention
        query_COO = Chem.MolFromSmarts("[$([O]-C(=O)-C)]")
        for atom, in mol.GetSubstructMatches(query_COO):
            atom: Chem.Atom = mol.GetAtomWithIdx(atom)
            if atom.GetPDBResidueInfo().GetIsHeteroAtom() == False:
                atom.SetFormalCharge(-1)
        return cls(mol=mol)
    
    def update_coordinates(self, coords: ndarray) -> None:
        rdmol = self.data
        if RDKitInterface.is_mol_3d(rdmol):
            conformer = rdmol.GetConformer()
            [conformer.SetAtomPosition(i, loc) for i, loc in enumerate(coords)]
        else:
            conformer = Chem.Conformer(rdmol.GetNumAtoms())
            conformer.Set3D(True)
            [conformer.SetAtomPosition(i, loc) for i, loc in enumerate(coords)]
            rdmol.AddConformer(conformer)

        
        
