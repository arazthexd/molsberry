from __future__ import annotations

from abc import ABC, abstractmethod
from typing import List, Dict, Any, Tuple, Type
from tqdm import tqdm

from ..pipeline import PipelineBlock
from ..data.collections import (
    BatchedData, Data, BatchOperator, 
    Representation, BatchedRep
)

from .batch_operator import BatchOperatorBlock

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

    def execute(self, *args: Tuple[Data], **kwargs: Dict[str, Data]) -> \
        Dict[str, Data | BatchedData]:
        out = super().execute(*args, **kwargs)
        for k, v in out.items():
            if k in self.input_batch_keys:
                out[k] = v.flatten()
        return out
    
    # def operate(self, input_dict: Dict[str, Representation | BatchedRep]) \
    #     -> Dict[str, Representation | BatchedRep]:

    #     if self.flatten:
    #         repbatch: BatchedRep = input_dict[self.input_keys[0]]
    #         out_reps: List[Representation] = []
    #         for rep in repbatch:
    #             rep: Representation
    #             out_contents = self.enumerate(rep.content)
    #             out_reps.extend([self.input_reps[0](c) for c in out_contents])
    #         return {self.input_keys[0]: BatchedRep(out_reps)}
        
    #     rep = input_dict[self.input_keys[0]]
    #     out_contents = self.enumerate(rep.content)
    #     return {self.input_keys[0]: BatchedRep([self.input_reps[0](c) 
    #                                             for c in out_contents])}

    # def meets_criteria(self, pot_input: Dict[str, Data | BatchedData]) -> bool:
    #     if self.flatten:
    #         assert all(isinstance(data, Data) for data in pot_input.values())
    #         is_batched = isinstance(pot_input[self.input_keys[0]], BatchedData)
    #         print(is_batched)
    #         if not is_batched: 
    #             return False
    #         depth = pot_input[self.input_keys[0]].depth
    #         return depth == 2
    #     else:
    #         return SimpleBlock.meets_criteria(self, pot_input)

# class SimpleConverterBlock(SimpleBlock, ABC):

#     @abstractmethod
#     def convert(self, data: Data) -> Data:
#         pass

#     def operate(self, input_dict: Dict[str, Data]) -> Dict[str, Data]:
#         return super().operate(input_dict)
    
# class SimpleEnumeratorBlock(SimpleBlock, ABC):
#     base_target_depth = 1

#     def __init__(self, debug: bool = False, save_output: bool = False,
#                  flatten: bool = False):
#         super().__init__(debug, save_output)
#         self.flatten = flatten

#     @abstractmethod
#     def enumerate(self, data: Data) -> BatchedData:
#         pass

#     def operate(self, batch: BatchedData) -> BatchedData:
#         if self.flatten:
#             data_list = []
#             [data_list.extend(self.enumerate(data) for data in batch)]
#             output = BatchedData(data_list)
#             assert len(output) == len(batch)
#             return output

#         batch_list = [self.enumerate(data) for data in batch]
#         output = BatchedData(batch_list)
#         assert len(output) - 1 == len(batch)
#         return output

# class SimpleSelectorBlock(SimpleBlock, ABC):
#     base_target_depth = 2

#     @abstractmethod
#     def select(self, batch: BatchedData) -> BatchedData:
#         pass

#     def operate(self, batch: BatchedData) -> BatchedData:
#         assert batch.depth == 2
#         newbatch = BatchedData([self.select(child_batch) 
#                                 for child_batch in batch])
#         assert newbatch.depth == 2
#         return newbatch
    
# class SimpleAnalyzerBlock(SimpleBlock, ABC):
#     base_target_depth = 1
#     opbatch_output_ids = []

#     @abstractmethod
#     def analyze(self, data: Data) -> Dict[str, Data]:
#         pass

#     def operate(self, batch: BatchedData) -> Dict[str, BatchedData]:
#         results_dict = {key: list() for key in self.output_context_keys}
#         for data in batch:
#             result: Dict[str, Data] = self.analyze(data)
#             [results_dict[k].append(v) for k, v in result.items()]
#         results_dict = {k: BatchedData(v) for k,v in results_dict.items()}
#         return results_dict

#     def apply_on(self, batch: BatchedData) -> Dict[str, BatchedData]:
#         if batch.depth == self.target_depth:
#             return self.operate(batch)
#         elif batch.depth > self.target_depth:
#             child_results_dict = {key: list() 
#                                   for key in self.output_context_keys}
#             for child in batch._item_list:
#                 child_result: Dict[str, BatchedData] = self.apply_on(child)
#                 [child_results_dict[k].append(v) for k, v in child_result.items()]
#             child_results_dict = {k: BatchedData(v) 
#                                   for k,v in child_results_dict.items()}
#             return child_results_dict
#         else:
#             ValueError()
    
#     def decompress_func(self, 
#                         newbatch: Dict[str, BatchedData]) -> Dict[BatchedData]:
#         return newbatch