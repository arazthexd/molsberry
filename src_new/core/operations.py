from abc import ABC, abstractmethod
from typing import List, Dict, Any, Tuple

from rdkit import Chem

from .pipeline import PipelineBlock

class Operator(ABC):
    @abstractmethod
    def operate(self, data: Any) -> Any:
        pass

    @abstractmethod
    def in_check(self, data: Any) -> bool:
        pass

    @abstractmethod
    def out_check(self, data: Any) -> bool:
        pass

class Enumerator(Operator, ABC):
    @abstractmethod
    def enumerate(self, data: Any) -> List[Any]:
        pass

    def operate(self, data: Any) -> List[Any]:
        assert self.in_check(data)
        output = self.enumerate(data)
        assert isinstance(output, list)
        assert all(self.in_check(o) for o in output)
        return output
    
class LigandEnumerationBlock(PipelineBlock, Enumerator, ABC):
    name = "Unnamed Ligand Enumeration"
    output_keys = ["enumerated_ligands"]
    required_input_keys = ["ligands"]
    optional_input_keys = []

    def execute(self, input_dict: Dict[str, Any]) -> Dict[str, Any]:
        
    
    def in_check(self, data: Any) -> bool:
        return isinstance(data, Chem.Mol)
    
    def out_check(self, data: Any) -> bool:
        return isinstance(data, Chem.Mol)


        
    

    