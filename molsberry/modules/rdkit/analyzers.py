from typing import Any, Dict, List, Type, Callable, Optional
from abc import ABC, abstractmethod
from tqdm import tqdm

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, rdForceFieldHelpers
from rdkit.ForceField import rdForceField

from ...core.data import (
    MoleculeData, Data, BatchedData, LigandData, UnspecifiedData, 
    Representation, UnspecifiedRep, ProteinData, PDBPathRep
)
from ...core.templates import SimpleBlock, PLInteractionJob

from .representations import RDKitMolRep, RDKitSmallMolRep
from .interface import RDKitInterface

class RDKitAnalyzerBlock(SimpleBlock, ABC):
    inputs = [
        ("molecules", MoleculeData, RDKitMolRep, False)
    ]
    batch_groups = []

    @abstractmethod
    def analyze(self, rdmol: Chem.Mol) -> Dict[str, Any]:
        pass

    def operate(self, input_dict: Dict[str, Representation]) \
        -> Dict[str, Representation]:

        content_dict = {k: v.content for k, v in input_dict.items()}
        result_dict = self.analyze(*[content_dict[k] for k in self.input_keys])
        return {k: self._get_out_rep(k)(v) for k, v in result_dict.items()}

class RDKitMWCalculator(RDKitInterface, RDKitAnalyzerBlock):
    name = "rdmolwtcalc"
    display_name = "RDKit Molecular Weight Calculator"
    outputs = [
        ("molwt", None, None, False)
    ]

    def analyze(self, rdmol: Chem.Mol) -> Dict[str, float]:
        return {"molwt": rdMolDescriptors.CalcExactMolWt(rdmol)}
    
class RDKitMMFFEnergyCalculator(RDKitInterface, RDKitAnalyzerBlock):
    name = "rdmmffenergycalc"
    display_name = "RDKit MMFF Energy Calculator"
    outputs = [
        ("energy", None, None, False) # TODO: Add numeric data types
    ]

    def calc_energy(self, rdmol:Chem.Mol) -> float:
        assert "H" in Chem.MolToSmiles(rdmol)
        assert self.is_mol_3d(rdmol)
        prop = rdForceFieldHelpers.MMFFGetMoleculeProperties(rdmol)
        ff = rdForceFieldHelpers.MMFFGetMoleculeForceField(rdmol, prop)
        ff: rdForceField.ForceField
        return ff.CalcEnergy()

    def analyze(self, rdmol: Chem.Mol) -> Dict[str, float]:
        return {"energy": self.calc_energy(rdmol)}
    
class RDKitPLInteractionCalculator(RDKitInterface, PLInteractionJob,
                                   RDKitAnalyzerBlock):
    # NOTE: Inheritance order is important...
    name = "rdprotliginter"
    display_name = "RDKit Protein Ligand Interaction Calculator"
    lig_rep = RDKitMolRep
    prot_rep = RDKitMolRep

    def __init__(self, debug: bool = False, save_output: bool = False) -> None:
        RDKitAnalyzerBlock.__init__(self, debug=debug, save_output=save_output)
        self.calculator = RDKitMMFFEnergyCalculator()
        energy_fn = self.calculator.calc_energy
        PLInteractionJob.__init__(self, energy_fn)

    def combine_pl(self, ligand, protein):
        return Chem.CombineMols(ligand, protein)
    
    def analyze(self, rdmol: Chem.Mol, rdprot: Chem.Mol) -> Dict[str, float]:
        interaction_dict = self.calc_interaction(rdmol, rdprot)
        return interaction_dict