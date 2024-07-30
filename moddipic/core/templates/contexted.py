from abc import ABC, abstractmethod
from typing import List, Dict, Any, Tuple, Type

class Contexted(ABC):

    def __init__(self) -> None:
        self.context = {k: None for k in self.context_keys}
        self.required_input_keys = list(set(self.required_input_keys + self.context_keys))

    @property
    @abstractmethod
    def context_keys(self):
        pass

    @property
    @abstractmethod
    def context_types(self):
        pass

    def _auto_execute_carry_exe(self, block_input: Dict[str, Any]) -> None:
        self.context = {k: block_input[k] for k in self.context_keys}
        assert all(type(v) == self.context_types(self.context_keys.index(k)) 
                   for k, v in self.context.items())
        self._output = self.execute(data=block_input[self.key])