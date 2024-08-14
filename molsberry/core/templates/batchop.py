from __future__ import annotations
from abc import ABC, abstractmethod
from typing import List, Dict, Any, Tuple, Type, Iterable, Optional
from tqdm import tqdm
from itertools import product, repeat

from mpire import WorkerPool

from ..pipeline import PipelineBlock
from ..data import (
    BatchOperator,
    BatchedData, Data, UnspecifiedData,
    BatchedRep, Representation, UnspecifiedRep
)

class BatchOperatorBlock(PipelineBlock, ABC):

    def __init__(self, debug: bool = False, save_output: bool = False,
                 num_workers: int | None = None) -> None:
        super().__init__(debug=debug, save_output=save_output)
        if num_workers is None:
            self.parallel = False
        else:
            self.parallel = True
            self.n_workers = num_workers

    @property
    @abstractmethod
    def inputs(self) -> \
        List[Tuple[str, 
                   Optional[Type[Data] | List[Type[Data]]], 
                   Optional[Type[Representation] | List[Type[Representation]]], 
                   bool]]:
        pass

    @property
    @abstractmethod
    def outputs(self) -> \
        List[Tuple[str, 
                   Optional[Type[Data]], 
                   Optional[Type[Representation] | List[Type[Representation]]], 
                   bool]]:
        pass

    @property
    @abstractmethod
    def batch_groups(self) -> List[Tuple[str]]:
        pass

    @property
    def input_keys(self) -> Tuple[str]:
        return tuple(zip(*self.inputs))[0]
    
    @property
    def input_dtypes(self) -> Tuple[Optional[Type[Data] | List[Type[Data]]]]:
        return tuple(zip(*self.inputs))[1]
    
    @property
    def input_reps(self) -> Tuple[Optional[Type[Representation] | 
                                  List[Type[Representation]]]]:
        return tuple(zip(*self.inputs))[2]
    
    @property
    def input_batch_keys(self) -> Tuple[str]:
        return (inp[0] for inp in self.inputs if inp[3] == False)
    
    @property
    def input_context_keys(self) -> Tuple[str]:
        return tuple([key for key in self.input_keys 
                      if key not in self.input_batch_keys])
    
    @property # TODO: Delete this...
    def input_batch_paired_keys(self) -> Tuple[str | Tuple[str]]:
        included = (k for g in self.batch_groups for k in g)
        return tuple(self.batch_groups)+tuple(k for k in self.input_batch_keys 
                                                if k not in included)
    
    @property
    def output_keys(self) -> Tuple[str]:
        return tuple(zip(*self.outputs))[0]
    
    @property
    def output_dtypes(self) -> Tuple[Optional[Type[Data]]]:
        return tuple(zip(*self.outputs))[1]
    
    @property
    def output_reps(self) -> Tuple[Optional[Type[Representation] | 
                                   List[Type[Representation]]]]:
        return tuple(zip(*self.outputs))[2]
    
    @property
    def output_batch_keys(self) -> Tuple[str]:
        return (out[0] for out in self.outputs if out[3] == False)
    
    @property
    def output_context_keys(self) -> Tuple[str]:
        return tuple([key for key in self.output_keys 
                      if key not in self.output_batch_keys])
    
    @abstractmethod
    def operate(self, input_dict: Dict[str, Representation]
                ) -> Dict[str, Representation]:
        pass

    @abstractmethod
    def meets_criteria(self, pot_input: Dict[str, Data]) -> bool:
        pass

    def set_num_workers(self, num_workers: int):
        self.parallel = True
        self.n_workers = num_workers
    
    def disable_parallel(self):
        self.parallel = False

    def split_input(self, full_input: Dict[str, Data]):
        self.batch_input = {k: v for k, v in full_input.items() 
                            if k in self.input_batch_keys}
        self.context_input = {k: v for k, v in full_input.items() 
                              if k in self.input_context_keys}

    def execute(self, *args: Tuple[Data], 
                **kwargs: Dict[str, Data]) -> Dict[str, Data | BatchedData]:

        # Combine positional and keyword arguments to make the primary input.
        full_input: Dict[str, Data] = dict()
        for i, arg in enumerate(args):
            full_input[self.input_keys[i]] = arg
        full_input.update(kwargs)

        if self.debug: print(f"Executing with input: {full_input}")

        # Split input into two groups: context and batch (saved as att.)
        self.split_input(full_input)
        
        # Unwrap the context input data into context input representations.
        self.context_input = self.unwrap_input(self.context_input)
        self.context_input: Dict[str, Representation | List[Representation]]
        
        # Create the context_output which can be modified inside `operate`.
        self.context_output = dict()
        self.context_output: Dict[str, Representation | List[Representation]]
        
        # Apply operation on the batch input, with the context input saved as
        # the instance's attribute.
        batch_output: Dict[str, BatchedData] = \
            self.apply_on_batch(self.batch_input)
        
        # Wrap context output representations into context output data.
        self.context_output = self.wrap_output(self.context_output)
        self.context_output: Dict[str, Data]
        
        # Merge batch and context output together and return results.
        final_output = batch_output.copy()
        final_output.update(self.context_output)
        return final_output

    def apply_on_batch(self, batch_input: Dict[str, Data]):
        
        # Make sure all batch inputs are batches. if not, make them.
        batch_input = {k: v if isinstance(v, BatchedData) else BatchedData([v])
                       for k, v in batch_input.items()}
        
        if self.debug: print(f"Apply on batch: {batch_input}")

        # Define independent and paired groups.
        # Independents would be combined using `product`.
        indep_groups: List[List[Dict[str, Data | BatchedData]]] = []
        for paired_keys in self.input_batch_paired_keys:
            paired_keys: str | Tuple[str]
            if isinstance(paired_keys, str):
                paired_keys = (paired_keys, )
            paired_keys: Tuple[str]

            paired_iter: List[Dict[str, Data | BatchedData]] = [
                dict(zip(paired_keys, vs)) 
                for vs in zip(*[batch_input[k] for k in paired_keys])
            ]

            indep_groups.append(paired_iter)
        
        # Add context to the indep_groups as a single-member iterator
        indep_groups.append([self.context_input])
    
        # Create iteratable consisting of all desired combinations of inputs.
        mixed_iter = product(*indep_groups)
        mixed_iter: Iterable[Tuple[Dict[str, Data]]]

        # The following variable will hold the outputs of `operate` method.
        batch_output: Dict[str, List[Data]] = {
            k: list() for k in self.output_batch_keys}

        # If parallel is on, use mpire for multiprocessing, else just do normal
        if self.parallel: # TODO: Customize multi processing
            self.parallel = False
            with WorkerPool(n_jobs=self.n_workers, shared_objects=self, 
                            enable_insights=True, use_dill=True) as pool:
                result_dicts = pool.map(self.potinpmembers_to_result, 
                                        list(zip(mixed_iter,)))
            self.parallel = True
            if self.debug: print(pool.get_insights())
        else:
            result_dicts = [self.potinpmembers_to_result(self, pim) 
                            for pim in mixed_iter]
        
        # Iterate over `result_dicts` + check + update batch_output
        for result_dict, potinpmembers in zip(result_dicts, product(*indep_groups)):
            if result_dict is None:
                # Create a potential input
                pot_input: Dict[str, Data | BatchedData] = dict()
                for member in potinpmembers:
                    pot_input.update(member)
                result_dict = self.apply_on_batch(pot_input, parallel=True)
            assert set(result_dict.keys()) == set(batch_output.keys())
            [batch_output[k].append(v) for k, v in result_dict.items()]
        
        # Convert batch output to be of type Dict[str, BatchedData]
        batch_output = {k: BatchedData(vl) 
                        for k, vl in batch_output.items()}
        return batch_output
    
    @staticmethod
    def potinpmembers_to_result(
        block: BatchOperatorBlock, 
        potinpmembers: Tuple[Dict[str, Data | BatchedData]]) \
            -> Dict[str, Data]:
        
        # Create a potential input
        pot_input: Dict[str, Data | BatchedData] = dict()
        for member in potinpmembers:
            pot_input.update(member)
        
        if block.debug: print(f"pot_input: {pot_input}")

        # When the potential input is created, check if it meets criteria
        # to be an input for `operate` method.
        meets_crit: bool = block.meets_criteria(pot_input)

        # If it does, unwrap the input into representations, run `operate`
        # and wrap the result from representations into data.
        if meets_crit:
            unwrapped_input = block.unwrap_input(pot_input)
            result_dict: Dict[str, Representation] = \
                block.operate(unwrapped_input)
            result_dict: Dict[str, Data] = block.wrap_output(result_dict)

        # If not, repeat the process once again on the input until we reach
        # a potential input that meets the criteria.
        else:
            # result_dict = None
            result_dict: Dict[str, Data] = block.apply_on_batch(pot_input,
                                                               parallel=True)

        return result_dict

    def unwrap_input(self, 
                     input_dict: Dict[str, Data | BatchedData]) -> \
                        Dict[str, Representation | List[Representation]]:
        
        if self.debug: print(f"Unwrapping input: {input_dict}")
        
        unwrapped_dict: Dict[str, Representation | List[Representation]] = {}
        for k, v in input_dict.items():
            v: BatchedData | Data

            print(self._get_inp_dtype(k))

            if isinstance(v, BatchedData):
                if isinstance(self._get_inp_dtype(k), list):
                    assert any(isinstance(v.basic_itype, t) 
                               for t in self._get_inp_dtype(k)) 
                else:
                    assert issubclass(v.basic_itype, self._get_inp_dtype(k))

            elif isinstance(v, Data):
                if isinstance(self._get_inp_dtype(k), list):
                    assert any(isinstance(v, t) 
                               for t in self._get_inp_dtype(k)) 
                else:
                    print(v, self._get_inp_dtype(k))
                    assert isinstance(v, self._get_inp_dtype(k))

            else:
                raise TypeError()  
            
            print(self._get_inp_rep(k))
            try:
                rep = v.get_representation(self._get_inp_rep(k))
            except:
                rep = [v.get_representation(r) for r in self._get_inp_rep(k)]
            rep: Representation | List[Representation]

            unwrapped_dict[k] = rep
        
        return unwrapped_dict

    def wrap_output(
            self, output_dict: Dict[str, Representation | List[Representation]]
            ) -> Dict[str, Data]:
        
        wrapped_dict: Dict[str, Data] = dict()
        for k, v in output_dict.items():
            v: Representation | BatchedRep | \
                List[Representation] | List[BatchedRep] | Any
            
            if isinstance(v, BatchedRep):
                assert issubclass(v.basic_itype, self._get_out_rep(k))
                data_list = [self._get_out_dtype(k)(rep) for rep in v]
                data = BatchedData(data_list)
            
            elif isinstance(v, Representation):
                assert isinstance(v, self._get_out_rep(k))
                data = self._get_out_dtype(k)(init_rep=v)

            elif not (isinstance(v, list) or isinstance(v, tuple)):
                v = UnspecifiedRep(v)
                data = self._get_out_dtype(k)(init_rep=v)

            else:
                data = self._get_out_dtype(k)(init_rep=None)
                for rep in v:
                    data.add_representation(rep)
            
            wrapped_dict[k] = data
        
        return wrapped_dict
    
    def _get_inp_dtype(self, key: str) -> Type[Data] | List[Type[Data]]:
        dtype = self.input_dtypes[self.input_keys.index(key)]
        if dtype is None:
            dtype = UnspecifiedData
        return dtype
    
    def _get_inp_rep(self, key: str) -> (Type[Representation] | 
                                         List[Type[Representation]]):
        rtype = self.input_reps[self.input_keys.index(key)]
        if rtype is None:
            rtype = UnspecifiedRep
        return rtype
    
    def _get_out_dtype(self, key:str) -> Type[Data]:
        dtype = self.output_dtypes[self.output_keys.index(key)]
        if dtype is None:
            dtype = UnspecifiedData
        return dtype
    
    def _get_out_rep(self, key:str) -> (Type[Representation] | 
                                        List[Type[Representation]]):
        rtype = self.output_reps[self.output_keys.index(key)]
        if rtype is None:
            rtype = UnspecifiedRep
        return rtype