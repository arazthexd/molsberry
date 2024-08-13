from __future__ import annotations

from abc import ABC, abstractmethod
from typing import List, Dict, Any, Tuple, Type

from ..pipeline import PipelineBlock
from ..data import BatchedData, Data

from .batchop import BatchOperatorBlock

class SimpleBlock(BatchOperatorBlock, ABC):
    
    def meets_criteria(self, pot_input: Dict[str, Data]) -> bool:
        assert all(isinstance(data, Data) for data in pot_input.values())
        return all(not isinstance(data, BatchedData) for key, data 
                   in pot_input.items() if key in self.input_batch_keys)
    
class SimpleEnumeratorBlock(SimpleBlock, ABC):
    
    def __init__(self, debug: bool = False, save_output: bool = False,
                 flatten: bool = False) -> None:
        super().__init__(debug, save_output)
        self.flatten = flatten

    @abstractmethod
    def enumerate(self, some_inp: Any) -> List[Any]:
        pass

    def execute(self, *args: Tuple[Data], **kwargs: Dict[str, Data]) -> \
        Dict[str, Data | BatchedData]:
        out = super().execute(*args, **kwargs)
        for k, v in out.items():
            if k in self.input_batch_keys:
                out[k] = v.flatten()
        return out

class SimpleSelectorBlock(SimpleBlock, ABC):

    @abstractmethod
    def select(self, some_inps: List[Any]) -> List[Any]:
        pass

    def meets_criteria(self, pot_input: Dict[str, Data | BatchedData]) -> bool:
        assert all(isinstance(data, Data) for data in pot_input.values())
        assert isinstance(pot_input[self.input_keys[0]], BatchedData)
        depth = pot_input[self.input_keys[0]].depth
        return depth == 1
    
    def split_input(self, full_input: Dict[str, Data]): 
        super().split_input(full_input)
        # TODO: Should I do this for everything?
        self.batch_input = {k: v if isinstance(v, BatchedData) 
                            else BatchedData([v])
                            for k, v in self.batch_input.items()}
        self.batch_input = {k: BatchedData([v]) 
                            for k, v in self.batch_input.items()}