from typing import List, Dict, Union, Literal
from abc import ABC, abstractmethod

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, Descriptors
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams

from ...core import MoleculeData, BatchedData, Representation, BatchedRep
from ...core import SimpleSelectorBlock

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

class RDKitMolFlagSelector(RDKitInterface, RDKitSelectorBlock):
    name = "rdmolflagsel"
    display_name = "RDkit Molecule Flag Filtering"

    def __init__(self, debug: bool = False, save_output: bool = False,
                num_workers: int | None = None, 
                filter_flag = Literal['PAINS_A', 'BRENK', 'NIH', 'ALL', 'CHEMBL', 'ZINC'] ) -> None:
        super().__init__(debug = debug, save_output = save_output, 
                        num_workers = num_workers)
        self.filter_flag = filter_flag
    
    def select(self, rdmols: List[Chem.Mol]) -> Chem.Mol | List[Chem.Mol]:
        params = FilterCatalogParams()

        try:  
            catalog_member = getattr(FilterCatalogParams.FilterCatalogs, self.filter_flag)    
            params.AddCatalog(catalog_member)  
        except AttributeError:  
            print(f"Error: '{self.filter_flag}' is not a valid catalog in FilterCatalogParams.FilterCatalogs.") 
        
        catalog = FilterCatalog(params)
        selected_ligs = []
        for mol in rdmols:
            flag = catalog.HasMatch(mol)    
            if flag:
                if self.debug:
                    get_flag = catalog.GetMatches(mol)
                    for entry in get_flag:
                        print(entry.GetDescription())
                continue
            else:
                selected_ligs.append(mol)
        return selected_ligs

class RDKitMolLIPINSKISelector(RDKitInterface, RDKitSelectorBlock):
    name = 'rdmollipinskisel'
    display_anme = 'RDKit Molecule LIPINSKI Filtering'

    def __init__(self, condition_count: int = 3, debug: bool = False, save_output: bool = False,
                num_workers: int | None = None) -> None:
        super().__init__(debug = debug, save_output = save_output, 
                        num_workers = num_workers)
        
        self.condition_count = condition_count

    def select(self, rdmols: List[Chem.Mol]) -> Chem.Mol | List[Chem.Mol]:
        selected_ligs = []
        for mol in rdmols:    
            MW = Descriptors.MolWt(mol)
            HBA = Descriptors.NOCount(mol)
            HBD = Descriptors.NHOHCount(mol)
            LogP = Descriptors.MolLogP(mol)
            conditions = [MW <= 500, HBA <= 10, HBD <= 5, LogP <= 5]
            fail_ro5 = conditions.count(False) >= self.condition_count
            if fail_ro5:
                if self.debug:
                    continue
            else:
                selected_ligs.append(mol)
        return selected_ligs