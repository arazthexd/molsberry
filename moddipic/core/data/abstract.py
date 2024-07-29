from abc import ABC, abstractmethod
from typing import List, Dict, Any, Tuple, Type

class Representation(ABC):
    def __init__(self, data: Any):
        self.data = data
    
    @property
    @abstractmethod
    def rep_name(self):
        pass

class SpecialDataClass(ABC):
    def __init__(self, init_rep: Representation | None = None):
        self._representations = dict()
        self._default_rep = None
        self._latest_rep = None # TODO: Implement these...
        if init_rep:
            if not isinstance(init_rep, Representation):
                raise TypeError()
            self.add_representation(init_rep)
    
    def add_representation(self, rep: Representation):
        self._representations[rep.rep_name] = rep
    
    def create_representation(self, rep_type: Type[Representation]):
        assert issubclass(rep_type, Representation)
        if not rep_type in [type(r) for r in self._representations.values()]:
            newrep = self._get_representation_from_type(rep_type)
            assert rep_type.rep_name not in list(self._representations.keys())
            # TODO: Above assertion --> Just a warning?
            self._representations[rep_type.rep_name] = newrep
    
    def get_representation(self, rep: str | Type[Representation]):
        
        if isinstance(rep, str):
            return self._get_representation_from_name(rep)
        
        self.create_representation(rep_type=rep)
        return self._get_representation_from_name(rep.rep_name)
    
    def get_data(self, rep: str | Type[Representation] | None = None):
        if rep is None:
            # TODO: After default and latest rep implementations
            raise NotImplementedError() 
        
        return self.get_representation(rep).data
    
    def _get_representation_from_name(self, rep_name: str) -> Representation:
        # TODO: Add a dictionary for type names as well...
        return self._representations[rep_name]
    
    def _get_representation_from_type(
            self, rep_type: Type[Representation]) -> Representation:
        
        for name, rep in self._representations.items():
            if issubclass(type(rep), rep_type): # TODO: Sure?
                return rep
        
        for name, rep in self._representations.items():
            rep_to_newrep_name = "to_" + rep_type.__name__
            newrep_from_rep_name = "from_" + type(rep).__name__
            if hasattr(rep, rep_to_newrep_name):
                return getattr(rep, rep_to_newrep_name)()
            if hasattr(rep_type, newrep_from_rep_name):
                return getattr(rep_type, newrep_from_rep_name)(rep)
        
        raise NotImplementedError()
