from abc import ABC, abstractmethod
from typing import List, Dict, Any, Tuple

from rdkit import Chem

class Data(ABC):
    @abstractmethod
    def get_basic_data_type(self):
        pass

    @abstractmethod
    def get_depth(self):
        pass

class Batched(Data):
    def __init__(self, data: list):
        assert isinstance(data, list)
        self.data = data

    def get_basic_data_type(self):
        for i, element in enumerate(self.data):
            if isinstance(element, Data):
                el_type = element.get_basic_data_type()
            else:
                el_type = type(element)
            if i == 0:
                ref_el_type = el_type
                continue
            assert ref_el_type == el_type
        return ref_el_type
    
    def get_depth(self):
        max_depth = 0
        for i, element in enumerate(self.data):
            if isinstance(element, Data):
                depth = 1 + element.get_depth()
            else:
                depth = 0
            max_depth = max(max_depth, depth)
        return max_depth
    
    # def apply(self, func):
    #     output = []
    #     for i, element in enumerate(self.data):
    #         if isinstance(element, Batched):
    #             output.append(element.apply(func))
    #         else:
    #             output.append(func(element))
    #     return Batched(output)
    
    # def flatten(self):
    #     output = []
    #     if any(not isinstance(element, Batched) for element in self.data):
    #         for element in self.data:
    #             output.extend()
    #     for i, element in enumerate(self.data):
    #         if isinstance(element, Batched):
    #             output.append(element.flatten())
    #         else:
    #             output.append(element)
    #     return Batched(output)
      
    def __getitem__(self, index):
        return self.data[index]
    
    def __iter__(self):
        return iter(self.data)

    def __len__(self):
        return len(self.data)
                
class PLComplex():
    def __init__(self, ligand: Chem.Mol, protein: str):
        self.ligand = ligand
        self.protein = protein
