from typing import Dict
from abc import ABC, abstractmethod
from tqdm import tqdm

from rdkit import Chem
from rdkit.Chem import rdDistGeom

from molsberry.core.data import Representation

from ...core.data import MoleculeData, Data, BatchedData, LigandData
from ...core.templates import SimpleBlock

from .representations import RDKitMolRep, RDKitSmallMolRep
from .interface import RDKitInterface

class RDKitConverterBlock(SimpleBlock, ABC):
    inputs = [
        ("molecules", MoleculeData, RDKitMolRep, False)
    ]
    outputs = [
        ("molecules", MoleculeData, RDKitMolRep, False)
    ]
    batch_groups = []

    @abstractmethod
    def convert(self, rdmol: Chem.Mol) -> Chem.Mol:
        pass

    def operate(self, input_dict: Dict[str, Representation]) \
        -> Dict[str, Representation]:

        rdmol = input_dict[self.input_keys[0]].content
        rdmol = self.convert(rdmol)
        return {self.input_keys[0]: self.input_reps[0](rdmol)}
        
class RDKitHydrogenAdder(RDKitInterface, RDKitConverterBlock):
    name = "rdhadder"
    display_name = "RDKit Hydrogen Adder"

    def convert(self, rdmol: Chem.Mol) -> Chem.Mol:
        return self.addhs_based_on_confdim(rdmol)

class RDKitLigandHAdder(RDKitHydrogenAdder):
    name = "rdhadder_lig"
    display_name = "RDKit Ligand Hydrogen Adder"
    inputs = [
        ("ligands", LigandData, RDKitSmallMolRep, False)
    ]
    outputs = [
        ("ligands", LigandData, RDKitSmallMolRep, False)
    ]

class RDKitLigandEmbedder(RDKitInterface, RDKitConverterBlock):
    name = "rdembedder_lig"
    display_name = "RDKit Ligand Embedder"
    inputs = [
        ("ligands", LigandData, RDKitSmallMolRep, False)
    ]
    outputs = [
        ("ligands", LigandData, RDKitSmallMolRep, False)
    ]

    def __init__(self, debug: bool = False, save_output: bool = False,
                 remove_hs: bool = False):
        super().__init__(debug, save_output)
        self.remove_hs = remove_hs

    def convert(self, rdmol: Chem.Mol) -> Chem.Mol:
        rdmol = Chem.AddHs(rdmol)
        rdDistGeom.EmbedMolecule(rdmol)
        if self.remove_hs:
            rdmol = Chem.RemoveHs(rdmol)
        return rdmol