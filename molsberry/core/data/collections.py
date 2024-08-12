from __future__ import annotations
from typing import List, Type, Any, Dict, Callable, Iterable
from abc import ABC, abstractmethod

import os

from ...core.data.abstract import Representation

from .abstract import Data

class BatchOperator(ABC):
    @property
    @abstractmethod
    def target_depth(self) -> int:
        pass
    
    @abstractmethod
    def operate(self, batch: Batched) -> Batched:
        pass

class Batched:
    # TODO: Assert if with every branch of the batched tree, type is same.

    def __init__(self, item_list: List[Any | Batched]):
        self._item_list = item_list
        self.obtain_batch_type()
    
    def __len__(self):
        return len(self._item_list)
    
    def __iter__(self):
        return iter(self._item_list)
    
    def __getitem__(self, index):
        return self._item_list[index]
    
    @property
    def batch_type(self):
        return self._batch_type # TODO: return in a format easier to understand
    
    @property
    def shallow_itype(self):
        return self._batch_type[1]
    
    @property
    def basic_itype(self):
        return self._batch_type[-1]
    
    @property
    def depth(self):
        return len(self._batch_type) - 1
    
    @classmethod
    def from_iterable(cls, iterable: Iterable, basic_itype: Type) -> Batched:
        # TODO: Better selection...
        selected_member = iterable[0]
        if isinstance(selected_member, basic_itype):
            return cls(list(iterable))
        
        batches = []
        for member in iterable:
            assert not isinstance(member, basic_itype)
            batch = cls.from_iterable(member, basic_itype)
            batches.append(batch)
        return cls(batches)
    
    def obtain_batch_type(self):
        selected_item = self._item_list[0]
        if isinstance(selected_item, Batched):
            child_type = selected_item.batch_type
        else:
            child_type = [type(selected_item)]

        self._batch_type: List[Type] = [self.__class__] + child_type
    
    def apply_operator(self, operator: BatchOperator) -> Batched: # DO NOT USE
        if self.depth == operator.target_depth:
            return operator.operate(self)
        elif self.depth > operator.target_depth:
            newbatches = []
            for batch in self._item_list:
                newbatch = batch.apply_operator(operator)
                newbatches.append(newbatch)
            return self.__class__(newbatches)
        else:
            ValueError()
    
    def flatten(self) -> Batched:
        if self.depth == 1:
            raise ValueError("Can't flatten anymore...")
        
        if self.depth == 2:
            flattened_list = [item  for child in self._item_list 
                              for item in child._item_list]
            return self.__class__(flattened_list)
        
        return self.__class__([child.flatten() for child in self._item_list])
    
    @classmethod
    def merge(cls, batches: List[Batched]):
        batch_lens = [len(batch) for batch in batches]
        assert len(set(batch_lens)) == 1
        batch_depths = [batch.depth for batch in batches]
        assert len(set(batch_depths)) == 1
        if batch_depths[0] == 1:
            merged_list = []
            for merged_items in zip(*[batch._item_list for batch in batches]):
                merged_list.append(cls(merged_items))
            return cls(merged_list)
        else:
            merged_children = []
            for child_batches in zip(batches):
                merged_child = cls.merge(list(child_batches))
                merged_children.append(merged_child)
            return cls(merged_children)

class BatchedRep(Batched, Representation):

    def __init__(self, rep_list: List[Representation | BatchedRep]):
        Batched.__init__(self, item_list=rep_list)
        self._item_list: List[Representation | BatchedRep]
        self.basic_itype: Representation
        self.shallow_itype: Representation | BatchedRep
    
    @property
    def content(self) -> Batched: # TODO: Or a simple list?
        return Batched([rep.content for rep in self._item_list])

    @property
    def rep_name(self):
        return self.basic_itype.rep_name
    
    def save_rep(self, exless_filename: str): 
        if self.depth == 1:
            self.shallow_itype.save_rep_batch(self._item_list, exless_filename)
        else:
            assert self.depth > 1
            for i, batchrep in enumerate(self._item_list):
                batchrep.save_rep(exless_filename+f"_{i+1}")

class BatchedData(Batched, Data):

    def __init__(self, data_list: List[Data | BatchedData]):
        Batched.__init__(self, item_list=data_list)
        self._item_list: List[Data | BatchedData]
        self.basic_itype: Data
        self.shallow_itype: Data | BatchedData

        Data.__init__(self, init_rep=None)
        self.obtain_batch_rep()
        self._representations: Dict[str, BatchedRep]
    
    def obtain_batch_rep(self) -> None:
        selected_data = self._item_list[0]
        rep_names = list(selected_data._representations.keys())
        for rep_name in rep_names:
            batchrep = self.get_representation(rep_name)
            self._representations[rep_name] = batchrep
        
    def add_representation(self, rep: BatchedRep) -> None:
        self._representations[rep.rep_name] = rep
        # TODO: This needs much more checking...
    
    def create_representation(self, rep_type: Type[Representation]) -> None:
        assert issubclass(rep_type, Representation)
        for data in self._item_list:
            data.create_representation(rep_type=rep_type)
        self._representations[rep_type.rep_name] = \
            self.get_representation(rep_type)
    
    def get_representation(
            self, rep: str | Type[Representation] | None) -> BatchedRep:
        
        if rep is None:
            return list(self._representations.values())[0]
        
        if isinstance(rep, str):
            out = self._get_representation_from_name(rep)
        elif issubclass(rep, Representation):
            out = self._get_representation_from_type(rep)
        else:
            raise TypeError("rep should be of type Representation.")
        if out:
            return out # TODO: What if it's None and rep is str?
            
        rep_list: List[Representation | BatchedRep] = list()
        if self.depth == 1:
            [rep_list.append(data.get_representation(rep)) 
             for data in self._item_list] 
            return BatchedRep(rep_list=rep_list)
        
        for data in self._item_list:
            batchrep: BatchedRep = data.get_representation(rep)
            rep_list.append(batchrep)
        return BatchedRep(rep_list=rep_list) # TODO: This method can be shorter
    
    def get_representation_content(
            self, rep: str | Type[Representation] | None = None) -> List[Any]:
        batchrep = self.get_representation(rep)
        return batchrep.content
    
    def _get_representation_from_type(
            self, rep_type: Type[Representation]) -> Representation | None:
        for name, rep in self._representations.items():
            if issubclass(rep.basic_itype, rep_type):
                return rep
        return None

            
