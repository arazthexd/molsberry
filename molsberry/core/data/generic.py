from typing import Any

from .abstract import Data, Representation

class NumericData(Data):
    def __init__(self, num, init_rep = None):
        super().__init__(init_rep)
        if not num:
            raise NotImplementedError
        if not isinstance(num, (int, float)):
            try:
                num = float(num)
            except:
                raise TypeError
        self.num = num
    
    def __add__(self, value):
        if isinstance(value, NumericData):
            return self.num + value.num
        else:
            return self.num + value
    
    def __sub__(self, value):
        if isinstance(value, NumericData):
            return self.num - value.num
        else:
            return self.num - value
    
    def __mul__(self, value):
        if isinstance(value, NumericData):
            return self.num * value.num
        else:
            return self.num * value
    
    def __truediv__(self, value):
        if isinstance(value, NumericData):
            if value.num == 0:
                raise ZeroDivisionError
            return self.num / value.num
        else:
            if value == 0:
                raise ZeroDivisionError
            return self.num / value
    
    def __eq__(self, value):
        if isinstance(value, NumericData):
            return self.num == value.num
        else:
            return self.num == value
    
    def _gt__(self, value):
        if isinstance(value, NumericData):
            return self.num > value.num
        else:
            return self.num > value
    
    def _lt__(self, value):
        if isinstance(value, NumericData):
            return self.num < value.num
        else:
            return self.num < value

    def __ge__(self, value):
        if isinstance(value, NumericData):
            return self.num >= value.num
        else:
            return self.num >= value
    
    def __le__(self, value):
        if isinstance(value, NumericData):
            return self.num <= value.num
        else:
            return self.num >= value.num
    
    def __repr__(self):
        return str(self.num)
    
    def __str__(self):
        return str(self.num)

class IntData(NumericData):
    def __init__(self, num, init_rep = None):
        super().__init__(num, init_rep)
        if not isinstance(num, int):
            try:
                num = int(num)
            except:
                raise ValueError

class FloatData(NumericData):
    def __init__(self, num, init_rep = None):
        super().__init__(num, init_rep)
        if not isinstance(num, float):
            try:
                num = float(num)
            except:
                raise ValueError

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
