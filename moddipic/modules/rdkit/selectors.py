from typing import List

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from ...core.templates import LigandSelectorBlock

from ...core.data.special_cls import Ligand
from .representations import RDKitMolRep
from .interface import RDKitInterface

class RDKitMWLigSelector(RDKitInterface, LigandSelectorBlock):
    name = "RDKit Ligand Weight Filtering"

    # NOTE: Deleted "identifier" from init...
    def __init__(self, max_wt: float = 600.0, min_wt: float = 30.0, 
                 debug: bool = False, save_output: bool = False):
        super().__init__(debug=debug, save_output=save_output)
        self.max_wt = max_wt
        self.min_wt = min_wt
    
    def select(self, ligands: List[Ligand]) -> List[Ligand]:
        rdmols = [self.special_cls_to_rdmol(ligand) for ligand in ligands]
        selected_ligs = []
        for rdmol in rdmols:
            if rdMolDescriptors.CalcExactMolWt(rdmol) > self.max_wt:
                if self.debug:
                    print(f"mol wt > {self.max_wt}: {Chem.MolToSmiles(rdmol)}")
                continue
            if rdMolDescriptors.CalcExactMolWt(rdmol) < self.min_wt:
                if self.debug:
                    print(f"mol wt < {self.min_wt}: {Chem.MolToSmiles(rdmol)}")
                continue
            selected_ligs.append(rdmol)
        
        return selected_ligs