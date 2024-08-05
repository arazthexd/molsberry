from numpy import ndarray

from ...core.data.abstract import Representation
from ...core.data.representations import SMILESRep, PDBPathProteinRep
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
    def from_PDBPathProteinRep(cls, pdb_rep: PDBPathProteinRep):
        assert isinstance(pdb_rep, PDBPathProteinRep)
        pdb_path = pdb_rep.data
        mol = Chem.MolFromPDBFile(pdb_path, removeHs=False)
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

        
        
