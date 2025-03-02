from typing import Any
import numpy as np

from .abstract import Data, Representation
from .unspecified import UnspecifiedRep

class NumericData(Data):
    def __init__(self, init_rep = None):
        super().__init__(init_rep)

class FloatRep(Representation):
    rep_name = 'FloatRepresentation'

    def __init__(self, num):
        super().__init__(num)

    @classmethod
    def from_UnspecifiedRep(cls, num:UnspecifiedRep):
        assert isinstance(num, UnspecifiedRep)
        num = num.content
        num = float(num)
        return cls(num = num)
        
    def __add__(self, value):
        if isinstance(value, FloatRep):
            return FloatRep(self.content + value.content)
        else:
            return FloatRep(self.content + value)
    
    def __mul__(self, value):
        if isinstance(value, FloatRep):
            return FloatRep(self.content * value.content)
        else:
            return FloatRep(self.content * value)
        
    def __add__(self, value):
        if isinstance(value, FloatRep):
            return FloatRep(self.content + value.content)
        else:
            return FloatRep(self.content + value)
    
    def __sub__(self, value):
        if isinstance(value, FloatRep):
            return FloatRep(self.content - value.content)
        else:
            return FloatRep(self.content - value)
    
    def __mul__(self, value):
        if isinstance(value, FloatRep):
            return FloatRep(self.content * value.content)
        else:
            return FloatRep(self.content * value)
    
    def __truediv__(self, value):
        if isinstance(value, FloatRep):
            if value.content == 0:
                raise ZeroDivisionError
            return FloatRep(self.content / value.content)
        else:
            if value == 0:
                raise ZeroDivisionError
            return FloatRep(self.content / value)
    
    def __eq__(self, value):
        if isinstance(value, FloatRep):
            return FloatRep(self.content == value.content)
        else:
            return FloatRep(self.content== value)
    
    def _gt__(self, value):
        if isinstance(value, FloatRep):
            return self.content > value.content
        else:
            return self.content > value
    
    def _lt__(self, value):
        if isinstance(value, FloatRep):
            return self.content < value.content
        else:
            return self.content < value

    def __ge__(self, value):
        if isinstance(value, FloatRep):
            return self.content >= value.content
        else:
            return self.content >= value
    
    def __le__(self, value):
        if isinstance(value, NumericData):
            return self.content <= value.content
        else:
            return self.content >= value.content
    
    def __repr__(self):
        return str(self.content)
    
    def __str__(self):
        return str(self.content)
    
class NpData(Data):  
    def __init__(self, init_rep = None):
        super().__init__(init_rep)

class NpArrayRep(Representation):
    rep_name = 'NumpyArrayRepresentation'

    def __init__(self, value):  
        if isinstance(value, np.ndarray):  
            self.content = value 
        else:  
            self.content = np.array(value)
    
    @classmethod
    def from_UnspecifiedRep(cls, npa:UnspecifiedRep):
        assert isinstance(npa, UnspecifiedRep)
        npa = npa.content
        npa = np.array(npa)
        return cls(npa = npa)

class StringData(Data):
    def __init__(self, init_rep = None):
        super().__init__(init_rep)

class StringRep(Representation):
    rep_name = 'StringRepresentation'
    def __init__(self, string):
        if not string:
            raise ValueError
        if not isinstance(string, str):
            try:
                string = str(string)
            except:
                raise ValueError
        self.content = string
    
    @classmethod
    def from_UnspecifiedRep(cls, string:UnspecifiedRep):
        assert isinstance(string, UnspecifiedRep)
        string = string.content
        string = str(string)
        return cls(string = string)

class BooleanData(Data):
    pass
    
class LocationData(Data):
    pass
