from typing import Any

from .abstract import Data, Representation

class NumericData(Data):
    def __init__(self, init_rep = None):
        super().__init__(init_rep)


class IntRep(Representation):
    rep_name = 'IntegerRepresentation'

    def __init__(self, num):
        if not num:
            raise ValueError
        elif not isinstance(num, (int, float)):
            try:
                num = int(num)
            except:
                raise TypeError
        self.content = num
    
    def __add__(self, value):
        if isinstance(value, IntRep):
            return IntRep(self.content + value.content)
        else:
            return IntRep(self.content + value)
    
    def __mul__(self, value):
        if isinstance(value, IntRep):
            return IntRep(self.content * value.content)
        else:
            return IntRep(self.content * value)
        
    def __add__(self, value):
        if isinstance(value, IntRep):
            return IntRep(self.content + value.content)
        else:
            return IntRep(self.content + value)
    
    def __sub__(self, value):
        if isinstance(value, IntRep):
            return IntRep(self.content - value.content)
        else:
            return FloatRep(self.content - value)
    
    def __mul__(self, value):
        if isinstance(value, IntRep):
            return IntRep(self.content * value.content)
        else:
            return IntRep(self.content * value)
    
    def __truediv__(self, value):
        if isinstance(value, IntRep):
            if value.content == 0:
                raise ZeroDivisionError
            return IntRep(self.content / value.content)
        else:
            if value == 0:
                raise ZeroDivisionError
            return IntRep(self.content / value)
    
    def __eq__(self, value):
        if isinstance(value, IntRep):
            return IntRep(self.content == value.content)
        else:
            return IntRep(self.content== value)

    def _gt__(self, value):
        if isinstance(value, IntRep):
            return self.content > value.content
        else:
            return self.content > value
    
    def _lt__(self, value):
        if isinstance(value, IntRep):
            return self.content < value.content
        else:
            return self.content < value

    def __ge__(self, value):
        if isinstance(value, IntRep):
            return self.content >= value.content
        else:
            return self.content >= value
    
    def __le__(self, value):
        if isinstance(value, IntRep):
            return self.content <= value.content
        else:
            return self.content >= value.content
    
    def __repr__(self):
        return str(self.content)
    
    def __str__(self):
        return str(self.content)
        
class FloatRep(Representation):
    rep_name = 'FloatRepresentation'

    def __init__(self, num):
        if not num:
            raise ValueError
        elif not isinstance(num, (int, float)):
            try:
                num = float(num)
            except:
                raise TypeError
        self.content = num
    
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

class StringData(Data):
    def __init__(self, string, init_rep = None):
        super().__init__(init_rep)
        if not string:
            raise NotImplementedError
        if not isinstance(string, str):
            try:
                string = str(string)
            except:
                raise ValueError

class BooleanData(Data):
    def __init__(self, boolean, init_rep = None):
        super().__init__(init_rep)
        if not isinstance(boolean, bool):
            try:
                boolean = bool(boolean)
            except:
                raise ValueError
    
class LocationData(Data):
    pass
