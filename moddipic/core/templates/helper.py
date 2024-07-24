from abc import ABC, abstractmethod
from typing import List, Dict, Any, Tuple
from tqdm import tqdm

from ..pipeline import PipelineBlock
from ..data import Batched

class SingleDataOperator(PipelineBlock, ABC):

    def __init__(self, debug: bool = False, save_output: bool = False) -> None:
        super().__init__(debug=debug, save_output=save_output)

        assert len(self.required_input_keys) == 1
        assert len(self.output_keys) == 1
        assert self.required_input_keys[0] == self.output_keys[0]
    
    @property
    def key(self) -> str:
        return self.required_input_keys[0]

    @property
    @abstractmethod
    def single_data_type(self) -> Any:
        pass

    def check_input(self, input_dict: Dict[str, Any]):
        super().check_input(input_dict)

        if isinstance(input_dict[self.key], Batched):
            input_dict: Dict[str, Batched]
            if self.debug: print("batched data basic type:",
                                 input_dict[self.key].get_basic_data_type())
            assert issubclass(
                input_dict[self.key].get_basic_data_type(), 
                self.single_data_type
            )
        else:
            assert isinstance(
                input_dict[self.key], 
                self.single_data_type
            )
    
    def check_output(self, output_dict: Dict[str, Any]):
        super().check_output(output_dict)

        if isinstance(output_dict[self.key], Batched):
            assert issubclass(
                output_dict[self.key].get_basic_data_type(), 
                self.single_data_type
            )
        else:
            assert isinstance(
                output_dict[self.key], 
                self.single_data_type
            )
    
    def pre_execute(self, data):
        if isinstance(data, self.single_data_type):
            data = Batched([data])

        if self.debug: print("Input Data Batch Depth:", data.get_depth())
        
        return data
    
    