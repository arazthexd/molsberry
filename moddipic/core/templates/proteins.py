from abc import ABC, abstractmethod
from typing import List, Dict, Any, Tuple
from tqdm import tqdm

import os

from ..data.collections import Batched
from .helper import (
    SingleDataOperator,
    SingleDataConverter, 
    SingleDataEnumerator, 
    SingleDataSelector
)

class PDBPathProteinSingleDataType:
    single_data_type = str
    required_input_keys = ["protein"]
    optional_input_keys = []
    output_keys = ["protein"]

class ProteinOperatorBlock(PDBPathProteinSingleDataType, 
                           SingleDataOperator, 
                           ABC):
    name = "Unnamed Protein Operator" 

class ProteinConverterBlock(PDBPathProteinSingleDataType, 
                            SingleDataConverter, 
                            ABC):
    name = "Unnamed Protein Converter"

class ProteinEnumeratorBlock(PDBPathProteinSingleDataType, 
                             SingleDataEnumerator, 
                             ABC):
    name = "Unnamed Protein Enumerator"

class ProteinSelectorBlock(PDBPathProteinSingleDataType, 
                           SingleDataSelector, 
                           ABC):
    name = "Unnamed Protein Selector"