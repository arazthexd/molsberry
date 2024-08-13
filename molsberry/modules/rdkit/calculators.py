from typing import Any, Dict, List, Type, Callable, Optional
from abc import ABC, abstractmethod
from tqdm import tqdm

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, rdForceFieldHelpers
from rdkit.ForceField import rdForceField

from ...core import Representation, Molecule3DRep, MoleculeData
from ...core import PLInteractionJob, EnergyJob, SimpleBlock

from .representations import RDKitMolRep, RDKitSmallMolRep
from .interface import RDKitInterface

class RDKitCalculatorBlock(SimpleBlock, ABC):
    inputs = [
        ("molecules", MoleculeData, RDKitMolRep, False)
    ]
    batch_groups = []

    @abstractmethod
    def calculate(self, rdmol: Chem.Mol, *args) -> Dict[str, Any]:
        pass

    def operate(self, input_dict: Dict[str, Representation]) \
        -> Dict[str, Representation]:

        content_dict = {k: v.content for k, v in input_dict.items()}
        result_dict = self.calculate(*[content_dict[k] for k in self.input_keys])
        return {k: self._get_out_rep(k)(v) for k, v in result_dict.items()}

class RDKitMWCalculator(RDKitInterface, RDKitCalculatorBlock):
    name = "rdmolwtcalc"
    display_name = "RDKit Molecular Weight Calculator"
    outputs = [
        ("molwt", None, None, False)
    ]

    def calculate(self, rdmol: Chem.Mol) -> Dict[str, float]:
        return {"molwt": rdMolDescriptors.CalcExactMolWt(rdmol)}
    
class RDKitMMFFEnergyCalculator(RDKitInterface, EnergyJob,
                                RDKitCalculatorBlock):
    name = "rdmmffenergycalc"
    display_name = "RDKit MMFF Energy Calculator"
    inputs = [
        ("molecules", MoleculeData, Molecule3DRep, False)
    ]

    def __init__(self, debug: bool = False, save_output: bool = False) -> None:
        RDKitCalculatorBlock.__init__(self, 
                                    debug=debug, save_output=save_output)
        EnergyJob.__init__(self, self.calc_energy)

    @staticmethod
    def calc_energy(rdmol:Chem.Mol) -> float:
        assert "H" in Chem.MolToSmiles(rdmol)
        assert RDKitInterface.is_mol_3d(rdmol)
        prop = rdForceFieldHelpers.MMFFGetMoleculeProperties(rdmol)
        ff = rdForceFieldHelpers.MMFFGetMoleculeForceField(rdmol, prop)
        ff: rdForceField.ForceField
        return ff.CalcEnergy()

    def calculate(self, rdmol: Chem.Mol) -> Dict[str, float]:
        return {"energy": self.calc_energy(rdmol)}
    
class RDKitPLInteractionCalculator(RDKitInterface, PLInteractionJob,
                                   RDKitCalculatorBlock):
    # NOTE: Inheritance order is important...
    name = "rdprotliginter"
    display_name = "RDKit Protein Ligand Interaction Calculator"
    lig_rep = RDKitMolRep
    prot_rep = RDKitMolRep

    def __init__(self, debug: bool = False, save_output: bool = False) -> None:
        RDKitCalculatorBlock.__init__(self, debug=debug, save_output=save_output)
        self.calculator = RDKitMMFFEnergyCalculator()
        energy_fn = self.calculator.calc_energy
        PLInteractionJob.__init__(self, energy_fn)

    def combine_mols(self, *mols):
        merge = mols[0]
        for mol in mols[1:]:
            merge = Chem.CombineMols(merge, mol)
        return merge
    
    def calculate(self, rdmol: Chem.Mol, rdprot: Chem.Mol) -> Dict[str, float]:
        interaction_dict = self.calc_interaction(rdmol, rdprot)
        return interaction_dict