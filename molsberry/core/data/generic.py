from typing import Any

from .abstract import Data, Representation

class NumericData(Data):
    pass

class IntData(NumericData):
    pass

class FloatData(NumericData):
    pass

class StringData(Data):
    pass

class BooleanData(Data):
    pass

class LocationData(Data):
    pass
