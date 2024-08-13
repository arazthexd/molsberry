from typing import Dict
from abc import ABC, abstractmethod
from tqdm import tqdm

from rdkit import Chem
from rdkit.Chem import rdDistGeom
from rdkit.Chem import rdMolDescriptors, rdForceFieldHelpers
from rdkit.ForceField import rdForceField

from molsberry.core.data import Representation

from ...core import MoleculeData, Data, BatchedData, LigandData
from ...core import SimpleBlock, PLOptimizeJob, OptimizeJob

from .representations import RDKitMolRep, RDKitSmallMolRep
from .interface import RDKitInterface
from .calculators import RDKitMMFFEnergyCalculator

class RDKitModifierBlock(SimpleBlock, ABC):
    inputs = [
        ("molecules", MoleculeData, RDKitMolRep, False)
    ]
    outputs = [
        ("molecules", MoleculeData, RDKitMolRep, False)
    ]
    batch_groups = []

    @abstractmethod
    def modify(self, rdmol: Chem.Mol) -> Chem.Mol:
        pass

    def operate(self, input_dict: Dict[str, Representation]) \
        -> Dict[str, Representation]:

        rdmol = input_dict[self.input_keys[0]].content
        rdmol = self.modify(rdmol)
        main_out_key = self.output_keys[0]
        return {main_out_key: self._get_out_rep(main_out_key)(rdmol)}

class RDKitLigandConverterBlock(RDKitModifierBlock, ABC):
    inputs = [
        ("ligands", LigandData, RDKitSmallMolRep, False)
    ]
    outputs = [
        ("ligands", LigandData, RDKitSmallMolRep, False)
    ]
        
class RDKitHydrogenAdder(RDKitInterface, RDKitModifierBlock):
    name = "rdhadder"
    display_name = "RDKit Hydrogen Adder"

    def modify(self, rdmol: Chem.Mol) -> Chem.Mol:
        return self.addhs_based_on_confdim(rdmol)

class RDKitLigandHAdder(RDKitLigandConverterBlock, RDKitHydrogenAdder):
    name = "rdhadder_lig"
    display_name = "RDKit Ligand Hydrogen Adder"

class RDKitLigandEmbedder(RDKitInterface, RDKitLigandConverterBlock):
    name = "rdembedder_lig"
    display_name = "RDKit Ligand Embedder"

    def __init__(self, debug: bool = False, save_output: bool = False,
                 remove_hs: bool = False):
        super().__init__(debug, save_output)
        self.remove_hs = remove_hs

    def modify(self, rdmol: Chem.Mol) -> Chem.Mol:
        rdmol = Chem.AddHs(rdmol)
        rdDistGeom.EmbedMolecule(rdmol)
        if self.remove_hs:
            rdmol = Chem.RemoveHs(rdmol)
        return rdmol
        


    