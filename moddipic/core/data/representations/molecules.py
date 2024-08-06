from typing import Any
from abc import ABC, abstractmethod

from numpy import ndarray
from ..abstract import Representation

import pathlib

class MoleculeRep(Representation, ABC):
    pass

class Molecule3DRep(MoleculeRep, ABC):
    @abstractmethod
    def update_coordinates(self, coords: ndarray):
        # coords must be n x 3
        pass

class Molecule1DRep(MoleculeRep, ABC):
    pass

class SmallMolRep(MoleculeRep, ABC):
    pass

class MacroMolRep(MoleculeRep, ABC):
    pass

class SMILESRep(SmallMolRep, Molecule1DRep):
    rep_name = "smiles"
    def __init__(self, smiles: str):
        assert isinstance(smiles, str)
        super().__init__(data=smiles)

class PDBPathRep(SmallMolRep, MacroMolRep, Molecule3DRep):
    rep_name = "pdb_path"
    def __init__(self, path: str):
        path = str(pathlib.Path(path).absolute())
        assert isinstance(path, str)
        super().__init__(data=path)