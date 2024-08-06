import numpy as np

from .abstract import SpecialDataClass, Representation
from .representations import SMILESRep, PDBPathRep
    
class Ligand(SpecialDataClass):
    @classmethod
    def from_smiles(cls, smiles: str):
        ligand = cls()
        ligand.add_representation(SMILESRep(smiles))
        return ligand
    
    def return_with_new_coords(self, coords: np.ndarray):
        newlig = self.copy()
        for rep in newlig._representations.values():
            rep: Representation
            rep.update_coordinates(coords)
        return newlig

class Protein(SpecialDataClass):
    @classmethod
    def from_pdb_path(cls, pdb_path: str):
        protein = cls()
        protein.add_representation(PDBPathRep(pdb_path))
        return protein
    
    def return_with_new_coords(self, coords: np.ndarray):
        newprot = self.copy()
        for rep in newprot._representations.values():
            rep: Representation
            rep.update_coordinates(coords)
        return newprot

class PLComplex(SpecialDataClass):
    pass

class PLComplexMulti(PLComplex):
    pass

# TODO: Do we need PLComplex and PLComplexMulti?