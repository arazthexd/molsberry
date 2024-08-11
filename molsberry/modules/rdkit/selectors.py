from typing import List
from abc import ABC, abstractmethod

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from ...core.data import MoleculeData, LigandData, BatchedData
from ...core.templates import SimpleBlock

from .representations import RDKitMolRep
from .interface import RDKitInterface

class RDKitSelectorBlock(SimpleBlock, ABC):
    block_in_rep = RDKitMolRep
    block_out_rep = RDKitMolRep
    block_out_dtype = MoleculeData

    @abstractmethod
    def select_from_rdmols(self, rdmols: List[Chem.Mol]) -> List[Chem.Mol]:
        pass

    def select(self, batch: BatchedData) -> BatchedData:
        rdmols: List[Chem.Mol] = \
            batch.get_representation_content(self.block_in_rep)
        rdmols = self.select_from_rdmols(rdmols)
        moldata_list = [self.block_out_dtype(rdmol) for rdmol in rdmols]
        return BatchedData(moldata_list)

class RDKitMolWtSelector(RDKitInterface, RDKitSelectorBlock):
    name = "rdmolwtsel"
    display_name = "RDKit Ligand Weight Filtering"
    input_keys = ["molecules"]
    output_keys = ["molecules"]

    def __init__(self, max_wt: float = 600.0, min_wt: float = 30.0, 
                 debug: bool = False, save_output: bool = False):
        super().__init__(debug=debug, save_output=save_output)
        self.max_wt = max_wt
        self.min_wt = min_wt
    
    def select_from_rdmols(self, rdmols: List[Chem.Mol]) -> List[Chem.Mol]:
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