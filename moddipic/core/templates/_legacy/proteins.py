from abc import ABC, abstractmethod
from typing import List, Dict, Any, Tuple
from tqdm import tqdm

import os

from ..pipeline import PipelineBlock
from ..data import Batched

class ProteinOperator(PipelineBlock, ABC):
    name = "Unnamed Protein Operator"
    required_input_keys = ["protein"]
    optional_input_keys = []
    output_keys = ["protein"]

    def __init__(self, debug: bool = False, save: bool = False) -> None:
        super().__init__(debug=debug, save_output=save)

        self._prot_folder_init = False
    
    def reset(self, removeExecutionHistory: bool = True, 
              removeConnections: bool = True) -> None:
        super().reset(removeExecutionHistory, removeConnections)

        self._prot_folder_init = False

    def check_input(self, input_dict: Dict[str, Any]):
        super().check_input(input_dict)

        if isinstance(input_dict["protein"], Batched):
            if self.debug: print("batched data basic type:",
                                 input_dict["protein"].get_basic_data_type())
            assert input_dict["protein"].get_basic_data_type() == str
        else:
            assert isinstance(input_dict["protein"], str)

    def check_output(self, output_dict: Dict[str, Any]):
        super().check_output(output_dict)   

        if isinstance(output_dict["protein"], Batched):
            assert issubclass(
                output_dict["protein"].get_basic_data_type(), str)
        else:
            assert isinstance(output_dict["protein"], str) 
            
    def save(self) -> str:
        txt_path = self._get_txt_path()
        raise NotImplementedError()
    
    def _get_txt_path(self) -> str:
        txt_path = os.path.join(self.base_dir, 
                                f"{self.name.replace(' ', '_').lower()}.txt")
        return txt_path
    
    def get_unique_file(self, suffix: str = ".pdb") -> str:
        raise NotImplementedError()

    @property
    def protfiles_dir(self):
        return os.path.join(self.base_dir, "proteins")


class ProteinConverter(ProteinOperator, ABC):
    name = "Unnamed Protein Converter"

    @abstractmethod
    def convert(self, protein: str) -> str:
        pass

