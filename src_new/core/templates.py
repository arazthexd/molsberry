from abc import ABC, abstractmethod
from typing import List, Dict, Any, Tuple
from tqdm import tqdm

from rdkit import Chem

from .pipeline import PipelineBlock
from .data import Data, Batched

# TODO: Change Data to Batched

class LigandOperatorBlock(PipelineBlock, ABC):
    name = "Unnamed Ligand Operator"
    required_input_keys = ["ligands"]
    optional_input_keys = []
    output_keys = ["ligands"]

    def check_input(self, input_dict: Dict[str, Any]):
        super().check_input(input_dict)

        if isinstance(input_dict["ligands"], Data):
            print(input_dict["ligands"].get_basic_data_type())
            assert issubclass(
                input_dict["ligands"].get_basic_data_type(), Chem.Mol)
        else:
            assert isinstance(input_dict["ligands"], Chem.Mol)
    
    def check_output(self, output_dict: Dict[str, Any]):
        super().check_output(output_dict)
        
        if isinstance(output_dict["ligands"], Data):
            assert issubclass(
                output_dict["ligands"].get_basic_data_type(), Chem.Mol)
        else:
            assert isinstance(output_dict["ligands"], Chem.Mol)


class LigandEnumeratorBlock(LigandOperatorBlock, ABC):
    name = "Unnamed Ligand Enumerator"
    def __init__(self, flatten: bool = False, debug: bool = False):
        super().__init__(debug)
        self.flatten = flatten

    def check_input(self, input_dict: Dict[str, Any]):
        super().check_input(input_dict)

        self._input_depth = 0
        if isinstance(input_dict["ligands"], Data):
            self._input_depth = input_dict["ligands"].get_depth()
    
    def check_output(self, output_dict: Dict[str, Any]):
        super().check_output(output_dict)

        assert isinstance(output_dict["ligands"], Batched)
        
        self._output_depth = output_dict["ligands"].get_depth()
        if self.flatten:
            assert self._output_depth == self._input_depth
        else:
            assert self._output_depth - 1 == self._input_depth
    
    @abstractmethod
    def enumerate(self, ligand: Chem.Mol) -> List[Chem.Mol]:
        pass

    def execute(self, ligands: Chem.Mol | Batched, 
                no_bar: bool = False) -> Dict[str, Any]:
        output = []
        
        if isinstance(ligands, Chem.Mol):
            ligands = Batched([ligands])

        if not no_bar:
            iterator = tqdm(ligands)
        else:
            iterator = ligands

        for ligand in iterator:
            if isinstance(ligand, Batched):
                # newligs = ligand.apply(helper_func)
                # output.append(newligs)
                output.append(self.execute(ligand, no_bar=True)["ligands"])
            else:
                newligs = self.enumerate(ligand)
                if self.flatten:
                    output.extend(newligs)
                else:
                    output.append(Batched(newligs))
        return {"ligands": Batched(output)}

class LigandConverterBlock(LigandOperatorBlock, ABC):
    pass

    

    