from abc import ABC, abstractmethod
from typing import List, Dict, Any, Tuple
from tqdm import tqdm

import os

from ..pipeline import PipelineBlock
from ..data.collections import Batched
from ..data.special_cls import Ligand
from ...utils.iotools import write_ligands
from .helper import (
    SingleDataOperator,
    SingleDataConverter, 
    SingleDataEnumerator, 
    SingleDataSelector
)

class LigandSingleDataType:
    single_data_type = Ligand
    required_input_keys = ["ligands"]
    optional_input_keys = []
    output_keys = ["ligands"]

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

    