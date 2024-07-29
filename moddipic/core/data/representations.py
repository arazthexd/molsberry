from typing import Any
from .abstract import Representation

import pathlib

class SMILESRep(Representation):
    rep_name = "smiles"
    def __init__(self, smiles: str):
        assert isinstance(smiles, str)
        super().__init__(data=smiles)

class PDBPathProteinRep(Representation):
    rep_name = "pdb_path"
    def __init__(self, path: str):
        path = str(pathlib.Path(path).absolute())
        assert isinstance(path, str)
        super().__init__(data=path)