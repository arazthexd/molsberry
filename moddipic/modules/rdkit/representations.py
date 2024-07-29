from ...core.data.abstract import Representation
from ...core.data.representations import SMILESRep

from rdkit import Chem

class RDKitMolRep(Representation):
    rep_name = "rdkit_mol"

    def __init__(self, mol: Chem.Mol):
        assert isinstance(mol, Chem.Mol)
        super().__init__(data=mol)
    
    @classmethod
    def from_SMILESRep(cls, smiles_rep: SMILESRep):
        assert isinstance(smiles_rep, SMILESRep)
        smiles = smiles_rep.data
        mol = Chem.MolFromSmiles(smiles)
        return cls(mol=mol)
        
