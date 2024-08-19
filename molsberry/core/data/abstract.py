from __future__ import annotations

from abc import ABC, abstractmethod
from typing import List, Dict, Any, Tuple, Type

import os
import glob
from copy import deepcopy
import pickle

from ...global_conf import DATA_UNIQUE_CODE_LEN
from ..utils import generate_random_str

class Representation(ABC):
    def __init__(self, content: Any):
        self.content = content
    
    @property
    @abstractmethod
    def rep_name(self):
        pass
    
    def save_rep(self, exless_filename: str): 
        # NOTE: This is only the default...
        # Make sure to rewrite it for most new representations...
        rep_path = exless_filename + ".pkl"
        with open(rep_path, "wb") as f:
            pickle.dump(self, f)

    @classmethod
    def save_rep_batch(cls, reps: List[Representation], exless_filename: str):
        for i, rep in enumerate(reps):
            cls.save_rep(rep, exless_filename+f"_{i+1}")

class Data(ABC):
    def __init__(self, init_rep: Representation | None = None):
        self._representations: Dict[str, Representation] = dict()
        self._default_rep = None
        self._latest_rep = None # TODO: Implement these...
        if init_rep:
            if not isinstance(init_rep, Representation):
                raise TypeError("init_rep should be of type Representation.")
            self.add_representation(init_rep)
        self.unique_code = None
    
    def add_representation(self, rep: Representation):
        assert isinstance(rep, Representation)
        self._representations[rep.rep_name] = rep
    
    def create_representation(self, rep_type: Type[Representation]):
        assert issubclass(rep_type, Representation)
        # assert rep_type.rep_name not in list(self._representations.keys())
        # TODO: Above assertion --> Just a warning?
        if self._get_representation_from_type(rep_type) is not None:
            return
        
        for name, rep in self._representations.items():
            rep_to_newrep_name = "to_" + rep_type.__name__
            newrep_from_rep_name = "from_" + type(rep).__name__
            if hasattr(rep, rep_to_newrep_name):
                newrep = getattr(rep, rep_to_newrep_name)()
                self._representations[rep_type.rep_name] = newrep
                return
            if hasattr(rep_type, newrep_from_rep_name):
                newrep = getattr(rep_type, newrep_from_rep_name)(rep)
                self._representations[rep_type.rep_name] = newrep
                return # TODO: Ways to improve on this (more depth)
            else: # TODO: Mix these two
                for parent_rep in type(rep).__mro__:
                    newrep_from_rep_name = "from_" + parent_rep.__name__
                    if hasattr(rep_type, newrep_from_rep_name):
                        newrep = getattr(rep_type, newrep_from_rep_name)(rep)
                        self._representations[rep_type.rep_name] = newrep
                        return
            
        raise NotImplementedError()
    
    def get_representation(self, rep: str | Type[Representation]):
        if isinstance(rep, str):
            out = self._get_representation_from_name(rep)
        elif issubclass(rep, Representation):
            out = self._get_representation_from_type(rep)
        else:
            raise TypeError("rep should be of type Representation.")
        if out:
            return out # TODO: What if it's None and rep is str?
        
        self.create_representation(rep_type=rep)
        return self._get_representation_from_name(rep.rep_name)
    
    def get_representation_content(
            self, rep: str | Type[Representation] | None = None) -> Any:
        if rep is None:
            # TODO: After default and latest rep implementations
            raise NotImplementedError() 
        return self.get_representation(rep).content
    
    def _get_representation_from_name(self, rep_name: str) -> Representation:
        return self._representations.get(rep_name)
    
    def _get_representation_from_type(
            self, rep_type: Type[Representation]) -> Representation:
        
        for name, rep in self._representations.items():
            if issubclass(type(rep), rep_type):
                return rep
        return None
        
    def copy(self) -> Data:
        return deepcopy(self) # TODO: Any better ways to copy more efficiently?
    
    def save(self, data_prefix: str) -> str:
        """Saves every representation of data given a prefix for path"""

        unique_code = self._get_data_unique_code(os.path.dirname(data_prefix))
        for rep_name, rep in self._representations.items():
            exless_filename = data_prefix + f"_{unique_code}_" + rep_name
            rep.save_rep(exless_filename=exless_filename)
    
    def _get_data_unique_code(self, data_dir: str) -> str:
        if not self.unique_code:
            while True:
                unique_code = generate_random_str(DATA_UNIQUE_CODE_LEN)
                matched_to_key = glob.glob(f"{data_dir}*{unique_code}")
                if len(matched_to_key) == 0:
                    self.unique_code = unique_code
                    return unique_code
        return self.unique_code
