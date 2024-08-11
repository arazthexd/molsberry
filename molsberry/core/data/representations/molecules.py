from typing import Any, List
from abc import ABC, abstractmethod
import os, shutil

from numpy import ndarray
from ..abstract import Representation

import pathlib

class MoleculeRep(Representation, ABC):
    pass

class Molecule3DRep(MoleculeRep, ABC):
    @abstractmethod
    def update_coordinates(self, coords: ndarray):
        assert coords.shape[1] == 3
        pass

class SmallMolRep(MoleculeRep, ABC):
    pass

class MacroMolRep(MoleculeRep, ABC):
    pass

class ProteinRep(MacroMolRep, ABC):
    pass

class SMILESRep(SmallMolRep):
    rep_name = "smiles"
    def __init__(self, smiles: str):
        assert isinstance(smiles, str)
        super().__init__(content=smiles)

    def save_rep(self, exless_filename: str):
        rep_path = exless_filename + ".smi"
        with open(rep_path, "w") as f:
            f.write(self.content)

    @classmethod
    def save_rep_batch(cls, reps: List[Representation], exless_filename: str):
        rep_path = exless_filename + ".smi"
        with open(rep_path, "w") as f:
            f.writelines([rep.content for rep in reps])

class PDBPathRep(SmallMolRep, MacroMolRep, Molecule3DRep):
    rep_name = "pdb_path"
    def __init__(self, path: str):
        path = str(pathlib.Path(path).absolute())
        assert isinstance(path, str)
        super().__init__(content=path)

    def save_rep(self, exless_filename: str):
        rep_path = exless_filename + ".pdb"
        shutil.copy(self.content, rep_path)