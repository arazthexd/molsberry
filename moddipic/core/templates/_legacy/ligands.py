from abc import ABC, abstractmethod
from typing import List, Dict, Any, Tuple
from tqdm import tqdm

from ...pipeline import PipelineBlock
from ...hi.collections import Batched

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

class SingleDataConverter(SingleDataOperator, ABC):
    
    @abstractmethod
    def convert(self, data: Any) -> Any:
        pass

    def execute(self, data: Any | Batched, 
                no_bar: bool = False) -> Dict[str, Any]:
        output = []
        data = self.pre_execute(data)

        if not no_bar:
            iterator = tqdm(data)
        else:
            iterator = data

        for element in iterator:
            if isinstance(element, Batched):
                output.append(self.execute(element, no_bar=True)[self.key])
            else:
                new_element = self.convert(element)
                output.append(new_element)
        
        return {self.key: Batched(output)}

class SingleDataEnumerator(SingleDataOperator, ABC):

    def __init__(self, flatten: bool = False, debug: bool = False,
                 save_output: bool = False):
        super().__init__(debug=debug, save_output=save_output)
        self.flatten = flatten

    @abstractmethod
    def enumerate(self, data: Any) -> Any:
        pass

    def check_input(self, input_dict: Dict[str, Any]):
        super().check_input(input_dict)

        self._input_depth = 0
        if isinstance(input_dict[self.key], Batched):
            self._input_depth = input_dict[self.key].get_depth()

    def check_output(self, output_dict: Dict[str, Any]):
        super().check_output(output_dict)

        assert isinstance(output_dict[self.key], Batched)
        
        self._output_depth = output_dict[self.key].get_depth()
        if self.flatten:
            assert self._output_depth == self._input_depth
        else:
            assert self._output_depth - 1 == self._input_depth
    
    def execute(self, data: Any | Batched, 
                no_bar: bool = False) -> Dict[str, Any]:
        output = []
        data = self.pre_execute(data)

        if not no_bar:
            iterator = tqdm(data)
        else:
            iterator = data
        
        for element in iterator:
            if isinstance(element, Batched):
                output.append(self.execute(element, no_bar=True)[self.key])
            else:
                new_element = self.enumerate(element)
                if self.flatten:
                    output.extend(new_element)
                else:
                    output.append(Batched(new_element))
        return {self.key: Batched(output)}

class SingleDataSelector(SingleDataOperator, ABC):
    
    @abstractmethod
    def select(self, data: List[Any]) -> List[Any]:
        pass
    
    def check_input(self, input_dict: Dict[str, Any]):
        super().check_input(input_dict)

        assert isinstance(input_dict[self.key], Batched)
        self._input_depth = input_dict[self.key].get_depth()
         
    def check_output(self, output_dict: Dict[str, Any]):
        super().check_output(output_dict)

        assert isinstance(output_dict[self.key], Batched)
        self._output_depth = output_dict[self.key].get_depth()
        assert self._output_depth == self._input_depth
    
    def execute(self, data: Batched, 
                no_bar: bool = False) -> Dict[str, Any]:
        output = []
        data = self.pre_execute(data)

        if data.get_depth() == 0:
            output = self.select(data.data)
            return {self.key: Batched(output)}

        if not no_bar:
            iterator = tqdm(data)
        else:
            iterator = data

        for element in iterator:
            output.append(self.execute(element, no_bar=True)[self.key])
        
        return {self.key: Batched(output)}