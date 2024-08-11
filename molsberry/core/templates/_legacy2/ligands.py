from abc import ABC, abstractmethod
from typing import List, Dict, Any, Tuple
from tqdm import tqdm

import os

from ..pipeline import PipelineBlock
from ..data.collections import Batched
from ..data.data_types import Ligand
from ...utils.iotools import write_ligands
from .single import (
    SingleDataOperator,
    SingleDataConverter, 
    SingleDataEnumerator, 
    SingleDataSelector,
    SingleDataAnalyzer
)

class LigandSingleDataType:
    single_data_type = Ligand
    required_input_keys = ["ligands"]
    optional_input_keys = []
    output_keys = ["ligands"]
    key = "ligands"

class LigandOperatorBlock(LigandSingleDataType, 
                          SingleDataOperator, 
                          ABC):
    name = "Unnamed Ligand Operator" 

class LigandConverterBlock(LigandSingleDataType, 
                           SingleDataConverter, 
                           ABC):
    name = "Unnamed Ligand Converter"

class LigandEnumeratorBlock(LigandSingleDataType, 
                            SingleDataEnumerator, 
                            ABC):
    name = "Unnamed Ligand Enumerator"

class LigandSelectorBlock(LigandSingleDataType, 
                          SingleDataSelector, 
                          ABC):
    name = "Unnamed Ligand Selector"

class LigandAnalyzerBlock(LigandSingleDataType, 
                          SingleDataAnalyzer, 
                          ABC):
    name = "Unnamed Ligand Analyzer"

    