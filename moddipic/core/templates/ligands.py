from abc import ABC, abstractmethod
from typing import List, Dict, Any, Tuple
from tqdm import tqdm

import os

from rdkit import Chem

from ..pipeline import PipelineBlock
from ..data import Data, Batched
from ...utils.iotools import write_ligands

# TODO: Change Data to Batched

class LigandOperatorBlock(PipelineBlock, ABC):
    name = "Unnamed Ligand Operator"
    required_input_keys = ["ligands"]
    optional_input_keys = []
    output_keys = ["ligands"]

    def check_input(self, input_dict: Dict[str, Any]):
        super().check_input(input_dict)

        if isinstance(input_dict["ligands"], Data):
            if self.debug: print("batched data basic type:",
                                 input_dict["ligands"].get_basic_data_type())
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
    
    def pre_execute(self, ligands):
        if isinstance(ligands, Chem.Mol):
            ligands = Batched([ligands])

        if self.debug: print("Input Ligands Batch Depth:", ligands.get_depth())
        
        return ligands
    
    def save(self) -> str:
        sdf_path = self._get_sdf_path()
        if isinstance(self._output["ligands"], Chem.Mol):
            write_ligands([self._output["ligands"]], sdf_path)
        elif isinstance(self._output["ligands"], Batched):
            d = self._output["ligands"].get_depth()
            if d == 0:
                write_ligands(self._output["ligands"].data, sdf_path)
            elif d == 1:
                for i, ligand_set in enumerate(self._output["ligands"]):
                    path = sdf_path[:-4] + f"_{i}.sdf"
                    write_ligands(ligand_set, path)
            else:
                raise ValueError()
    
    def _get_sdf_path(self) -> str:
        sdf_path = os.path.join(self.base_dir, 
                                f"{self.name.replace(' ', '_').lower()}.sdf")
        return sdf_path

class LigandEnumeratorBlock(LigandOperatorBlock, ABC):
    name = "Unnamed Ligand Enumerator"

    def __init__(self, flatten: bool = False, debug: bool = False,
                 save: bool = False):
        super().__init__(debug=debug, save_output=save)
        self.flatten = flatten

    def check_input(self, input_dict: Dict[str, Any]):
        super().check_input(input_dict)

        self._input_depth = 0
        if isinstance(input_dict["ligands"], Batched):
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
        ligands = self.pre_execute(ligands)

        if not no_bar:
            iterator = tqdm(ligands)
        else:
            iterator = ligands
        
        for ligand in iterator:
            if isinstance(ligand, Batched):
                output.append(self.execute(ligand, no_bar=True)["ligands"])
            else:
                newligs = self.enumerate(ligand)
                if self.flatten:
                    output.extend(newligs)
                else:
                    output.append(Batched(newligs))
        return {"ligands": Batched(output)}

class LigandConverterBlock(LigandOperatorBlock, ABC):
    name = "Unnamed Ligand Converter"

    @abstractmethod
    def convert(self, ligand: Chem.Mol) -> Chem.Mol:
        pass

    def execute(self, ligands: Chem.Mol | Batched, 
                no_bar: bool = False) -> Dict[str, Any]:
        output = []
        ligands = self.pre_execute(ligands)

        if not no_bar:
            iterator = tqdm(ligands)
        else:
            iterator = ligands

        for ligand in iterator:
            if isinstance(ligand, Batched):
                output.append(self.execute(ligand, no_bar=True)["ligands"])
            else:
                newlig = self.convert(ligand)
                output.append(newlig)
        
        return {"ligands": Batched(output)}

class LigandSelectorBlock(LigandOperatorBlock, ABC):
    name = "Unnamed Ligand Selector"

    @abstractmethod
    def select(self, ligands: List[Chem.Mol]) -> List[Chem.Mol]:
        pass

    def check_input(self, input_dict: Dict[str, Any]):
        super().check_input(input_dict)

        assert isinstance(input_dict["ligands"], Batched)
        self._input_depth = input_dict["ligands"].get_depth()
         
    def check_output(self, output_dict: Dict[str, Any]):
        super().check_output(output_dict)

        assert isinstance(output_dict["ligands"], Batched)
        self._output_depth = output_dict["ligands"].get_depth()
        assert self._output_depth == self._input_depth
        

    def execute(self, ligands: Batched, 
                no_bar: bool = False) -> Dict[str, Any]:
        output = []
        ligands = self.pre_execute(ligands)

        if ligands.get_depth() == 0:
            output = self.select(ligands.data)
            return {"ligands": Batched(output)}

        if not no_bar:
            iterator = tqdm(ligands)
        else:
            iterator = ligands

        for ligand in iterator:
            output.append(self.execute(ligand, no_bar=True)["ligands"])
            # if isinstance(ligand, Batched):
            #     if ligand.get_depth() > 0:
            #         output.append(self.execute(ligand, no_bar=True)["ligands"])
            #     else:
            #         output.append(self.select(ligand.data))
            # else:
            #     raise TypeError()
        
        return {"ligands": Batched(output)}
        


    

    