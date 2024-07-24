from typing import List

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from ..core.templates import LigandSelectorBlock

class RDKitMWLigSelector(LigandSelectorBlock):
    name = "RDKit Ligand Weight Filtering"
    def __init__(self, identifier: str = "all", 
                 max_wt: float = 600.0, min_wt: float = 30.0, 
                 debug: bool = False):
        super().__init__(debug)
        self.max_wt = max_wt
        self.min_wt = min_wt
    
    def select(self, ligands: List[Chem.Mol]) -> List[Chem.Mol]:
        selected_ligs = []
        for ligand in ligands:
            if rdMolDescriptors.CalcExactMolWt(ligand) > self.max_wt:
                if self.debug:
                    print(f"mol wt > {self.max_wt}: {Chem.MolToSmiles(ligand)}")
                continue
            if rdMolDescriptors.CalcExactMolWt(ligand) < self.min_wt:
                if self.debug:
                    print(f"mol wt < {self.min_wt}: {Chem.MolToSmiles(ligand)}")
                continue
            selected_ligs.append(ligand)
        
        return selected_ligs