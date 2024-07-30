from abc import ABC, abstractmethod
from typing import List, Dict, Any, Tuple, Type

class Contexted(ABC):
    input_context_keys, input_context_types = [], []
    output_context_keys, output_context_types = [], []
    # TODO: Output context should have a setter similar to `output`

    def __init__(self) -> None:
        self.input_context = {k: None for k in self.input_context_keys}
        self.required_input_keys = list(set(self.required_input_keys + 
                                            self.input_context_keys))
        
        self.output_context = {k: None for k in self.output_context_keys}
        self.output_keys = list(set(self.output_keys + 
                                    self.output_context_keys))

    # @property
    # @abstractmethod
    # def input_context_keys(self):
    #     pass

    # @property
    # @abstractmethod
    # def input_context_types(self):
    #     pass

    # @property
    # @abstractmethod
    # def output_context_keys(self):
    #     pass

    # @property
    # @abstractmethod
    # def output_context_types(self):
    #     pass

    def _auto_execute_carry_exe(self, block_input: Dict[str, Any]) -> None:
        self.input_context = {k: block_input[k] 
                              for k in self.input_context_keys}
        assert all(type(v) == self.input_context_types(
            self.input_context_keys.index(k)) for k, v in 
            self.input_context.items())
        
        self._output = self.execute(data=block_input[self.key])
        assert all(type(v) == self.output_context_types(
            self.output_context_keys.index(k)) for k, v in 
            self.output_context.items())
        
        for key in self.output_context_keys:
            self._output[key] = self.output_context[key]
