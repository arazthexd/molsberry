from abc import ABC, abstractmethod
from typing import List, Dict, Any, Tuple
from tqdm import tqdm

import os

from rdkit import Chem

from ..pipeline import PipelineBlock
from ..data.collections import Data, Batched
from ...utils.iotools import write_ligands
from .helper import (
    SingleDataOperator,
    SingleDataConverter, 
    SingleDataEnumerator, 
    SingleDataSelector
)

class RDKitLigandSingleDataType:
    single_data_type = Chem.Mol
    required_input_keys = ["ligands"]
    optional_input_keys = []
    output_keys = ["ligands"]

class LigandOperatorBlock(RDKitLigandSingleDataType, 
                          SingleDataOperator, 
                          ABC):
    name = "Unnamed Ligand Operator" 

class LigandConverterBlock(RDKitLigandSingleDataType, 
                           SingleDataConverter, 
                           ABC):
    name = "Unnamed Ligand Converter"

class LigandEnumeratorBlock(RDKitLigandSingleDataType, 
                            SingleDataEnumerator, 
                            ABC):
    name = "Unnamed Ligand Enumerator"

class LigandSelectorBlock(RDKitLigandSingleDataType, 
                          SingleDataSelector, 
                          ABC):
    name = "Unnamed Ligand Selector"

    