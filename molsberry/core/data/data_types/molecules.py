from typing import Any
import warnings

import numpy as np

from ..abstract import Data, Representation
from ..representations import SMILESRep, PDBPathRep, Molecule3DRep

class MoleculeData(Data):
    def return_with_new_coords(self, coords: np.ndarray):
        newmol = self.copy()
        for rep in newmol._representations.values():
            if isinstance(rep, Molecule3DRep):
                rep.update_coordinates(coords)
        return newmol

class LigandData(MoleculeData):
    @classmethod
    def from_smiles(cls, smiles: str):
        ligand = cls()
        ligand.add_representation(SMILESRep(smiles))
        return ligand

class ProteinData(MoleculeData):
    @classmethod
    def from_pdb_path(cls, pdb_path: str):
        protein = cls()
        protein.add_representation(PDBPathRep(pdb_path))
        return protein
    
# class PLComplex(SpecialDataClass):
#     pass

# class PLComplexMulti(PLComplex):
#     pass

# TODO: Do we need PLComplex and PLComplexMulti?