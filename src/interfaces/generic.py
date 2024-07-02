
import os
from typing import List, Tuple
from rdkit import Chem
from ..abstract import *
from ..utils.pockets import isolate_pocket, PocketLocation
from ..utils.iotools import get_unique_full_name

class GenericLigandSelector(LigandSelector):
    name = "Generic Ligand Selection"
    def __init__(self, sort_by: str, identifier: str = "_Name", max_per_id: int = 1, 
                 reverse: bool = True):
        super().__init__(identifier)
        self.max_per_id = max_per_id
        self.sort_by = sort_by
        self.reverse = reverse
    
    def select(self, ligands: List[Chem.Mol]) -> List[Chem.Mol]:
        return sorted(ligands, key=lambda lig: float(lig.GetProp(self.sort_by)), 
                      reverse=self.reverse)[:self.max_per_id]

class PocketIsolator(ComplexConverter):
    name = "Pocket Isolation"
    def __init__(self, work_dir: str = "tmp"):
        super().__init__()
        self.work_dir = work_dir

    def convert(self, ligand: Chem.Mol, target_path: str) -> Tuple[Chem.Mol, str]:
        output_prefix = os.path.join(self.work_dir, "isolated_pocket")
        output_path = get_unique_full_name(output_prefix, ".pdb", n=5)
        pocket_loc = PocketLocation(method="ligand", ligand=ligand)
        isolate_pocket(target_path, pocket_loc, output_path)
        return ligand, output_path