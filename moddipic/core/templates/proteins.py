from abc import ABC, abstractmethod
from typing import List, Dict, Any, Tuple
from tqdm import tqdm

import os

from ..data.collections import Batched
from ..data.special_cls import Protein
from .single_data import (
    SingleDataOperator,
    SingleDataConverter, 
    SingleDataEnumerator, 
    SingleDataSelector,
    SingleDataAnalyzer
)

class ProteinSingleDataType:
    single_data_type = Protein
    required_input_keys = ["protein"]
    optional_input_keys = []
    output_keys = ["protein"]
    key = "protein"

class ProteinOperatorBlock(ProteinSingleDataType, 
                           SingleDataOperator, 
                           ABC):
    name = "Unnamed Protein Operator" 

class ProteinConverterBlock(ProteinSingleDataType, 
                            SingleDataConverter, 
                            ABC):
    name = "Unnamed Protein Converter"

class ProteinEnumeratorBlock(ProteinSingleDataType, 
                             SingleDataEnumerator, 
                             ABC):
    name = "Unnamed Protein Enumerator"

class ProteinSelectorBlock(ProteinSingleDataType, 
                           SingleDataSelector, 
                           ABC):
    name = "Unnamed Protein Selector"

class ProteinAnalyzerBlock(ProteinSingleDataType, 
                           SingleDataAnalyzer, 
                           ABC):
    name = "Unnamed Protein Analyzer"