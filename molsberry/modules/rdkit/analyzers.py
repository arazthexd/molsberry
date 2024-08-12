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
from ...core.templates import SimpleBlock

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
    
class RDKitPLInteractionCalculator(RDKitInterface, RDKitAnalyzerBlock):
    name = "rdprotliginter"
    display_name = "RDKit Protein Ligand Interaction Calculator"
    inputs = [
        ("ligands", LigandData, RDKitMolRep, False),
        ("proteins", ProteinData, RDKitMolRep, False)
    ]
    outputs = [
        ("e_interaction", None, None, False),
        ("e_ligand", None, None, False),
        ("e_protein", None, None, False),
        ("e_complex", None, None, False)
    ]
    batch_groups = [("ligands", "proteins")]

    def __init__(self, energy_fn: Optional[Callable] = None,
                 debug: bool = False, save_output: bool = False) -> None:
        super().__init__(debug=debug, save_output=save_output)
        if energy_fn is None:
            calculator = RDKitMMFFEnergyCalculator()
            energy_fn = calculator.calc_energy
        self.energy_fn = energy_fn
    
    def analyze(self, rdmol: Chem.Mol, rdprot: Chem.Mol) -> Dict[str, float]:
        e_ligand = self.energy_fn(rdmol)
        e_protein = self.energy_fn(rdprot)
        
        pl_complex = Chem.CombineMols(rdmol, rdprot)
        e_complex = self.energy_fn(pl_complex)

        return {
            "e_interaction": e_complex - e_ligand - e_protein,
            "e_ligand": e_ligand,
            "e_protein": e_protein,
            "e_complex": e_complex
        }
        
# class RDKitPLInteractionCalculator(RDKitInterface, RDKitAnalyzerBlock):
#     output_keys = ["interaction"]
#     output_types = [float]
#     input_context_keys = ["protein"]
#     input_context_types = [Protein]
#     name = "rdprotliginter"
#     display_name = "RDKit Protein Ligand Interaction Calculator"
#     input_keys = ["molecules"]
#     output_keys = ["molwt"]

#     def __init__(self, debug: bool = False, save_output: bool = False):
#         LigandAnalyzerBlock.__init__(self, debug=debug, save_output=save_output)
#         Contexted.__init__(self)
    
#     def analyze(self, ligand: Ligand) -> Dict[str, float]:
#         lig = special_cls_to_rdmol(ligand)
#         prot = special_cls_to_rdmol(self.input_context["protein"])
#         merg = Chem.CombineMols(lig, prot)

#         ligprops = rdForceFieldHelpers.MMFFGetMoleculeProperties(lig)
#         ligff = rdForceFieldHelpers.MMFFGetMoleculeForceField(lig, ligprops)
#         ligff: rdForceField.ForceField
#         lig_energy = ligff.CalcEnergy()

#         protprops = rdForceFieldHelpers.MMFFGetMoleculeProperties(prot)
#         protff = rdForceFieldHelpers.MMFFGetMoleculeForceField(prot, protprops)
#         protff: rdForceField.ForceField
#         prot_energy = protff.CalcEnergy()

#         mergprops = rdForceFieldHelpers.MMFFGetMoleculeProperties(merg)
#         mergff = rdForceFieldHelpers.MMFFGetMoleculeForceField(merg, mergprops)
#         mergff: rdForceField.ForceField
#         merg_energy = mergff.CalcEnergy()

#         return {"interaction": merg_energy - lig_energy - prot_energy}