from typing import List, Dict, Union
from abc import ABC, abstractmethod

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from ...core.data import MoleculeData, BatchedData, Representation, BatchedRep
from ...core.templates import SimpleSelectorBlock

from .representations import RDKitMolRep
from .interface import RDKitInterface

class RDKitSelectorBlock(SimpleSelectorBlock, ABC):
    inputs = [
        ("molecules", MoleculeData, RDKitMolRep, False)
    ]
    outputs = [
        ("molecules", MoleculeData, RDKitMolRep, False)
    ]
    batch_groups = []

    @abstractmethod
    def select(self, rdmols: List[Chem.Mol]) -> Union[Chem.Mol, List[Chem.Mol]]:
        pass

    def operate(self, input_dict: Dict[str, BatchedRep]) \
        -> Dict[str, Representation | BatchedRep]:
        batch = input_dict[self.input_keys[0]] # depth = 1
        rdmols: List[Chem.Mol] = batch.content._item_list
        rdmols: Chem.Mol | List[Chem.Mol] = self.select(rdmols)
        
        main_out_key = self.output_keys[0]

        if isinstance(rdmols, Chem.Mol):
            out = self._get_out_rep(main_out_key)(rdmols)
        else:
            out = BatchedRep([self._get_out_rep(main_out_key)(rdmol) 
                              for rdmol in rdmols])
        
        return {main_out_key: out}

class RDKitMolWtSelector(RDKitInterface, RDKitSelectorBlock):
    name = "rdmolwtsel"
    display_name = "RDKit Molecule Weight Filtering"

    def __init__(self, max_wt: float = 600.0, min_wt: float = 30.0, 
                 debug: bool = False, save_output: bool = False):
        super().__init__(debug=debug, save_output=save_output)
        self.max_wt = max_wt
        self.min_wt = min_wt
    
    def select(self, rdmols: List[Chem.Mol]) -> List[Chem.Mol]:
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