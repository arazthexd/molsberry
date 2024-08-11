from typing import Any, Dict, List, Type
from abc import ABC, abstractmethod
from tqdm import tqdm

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, rdForceFieldHelpers
from rdkit.ForceField import rdForceField

from ...core.data import (
    MoleculeData, Data, BatchedData, LigandData, UnspecifiedData, 
    Representation, UnspecifiedRep
)
from ...core.templates import SimpleBlock

from .representations import RDKitMolRep, RDKitSmallMolRep
from .interface import RDKitInterface

class RDKitAnalyzerBlock(SimpleBlock, ABC):
    block_in_rep = RDKitMolRep

    def __init__(self, debug: bool = False, save_output: bool = False):
        super().__init__(debug, save_output)
        assert len(self.block_out_dtype) == \
            len(self.output_keys) - len(self.opbatch_output_ids)
    
    @property
    def block_out_dtype(self) -> List[Type[Data]]:
        return [UnspecifiedData] * (len(self.output_keys) - 
                                    len(self.opbatch_output_ids))
    
    @property
    def block_out_rep(self) -> List[Type[Representation]]:
        return [UnspecifiedRep] * (len(self.output_keys) - 
                                   len(self.opbatch_output_ids))

    @abstractmethod
    def analyze_rdmol(self, rdmol: Chem.Mol) -> Dict[str, Any]:
        pass

    def analyze(self, data: MoleculeData) -> Dict[str, Data]:
        rdmol = data.get_representation_content(self.block_in_rep)
        result_dict = self.analyze_rdmol(rdmol)
        result_dict = {k: self.block_out_dtype[i](self.block_out_rep[i](v)) 
                       for i, (k, v) in enumerate(result_dict.items())}
        return result_dict

class RDKitMWCalculator(RDKitInterface, RDKitAnalyzerBlock):
    name = "rdmolwtcalc"
    display_name = "RDKit Molecular Weight Calculator"
    input_keys = ["molecules"]
    output_keys = ["molwt"]

    def analyze_rdmol(self, rdmol: Chem.Mol) -> Dict[str, float]:
        return {"molwt": rdMolDescriptors.CalcExactMolWt(rdmol)}
 
class RDKitPLInteractionCalculator(RDKitInterface, RDKitAnalyzerBlock):
    output_keys = ["interaction"]
    output_types = [float]
    input_context_keys = ["protein"]
    input_context_types = [Protein]
    name = "rdprotliginter"
    display_name = "RDKit Protein Ligand Interaction Calculator"
    input_keys = ["molecules"]
    output_keys = ["molwt"]

    def __init__(self, debug: bool = False, save_output: bool = False):
        LigandAnalyzerBlock.__init__(self, debug=debug, save_output=save_output)
        Contexted.__init__(self)
    
    def analyze(self, ligand: Ligand) -> Dict[str, float]:
        lig = special_cls_to_rdmol(ligand)
        prot = special_cls_to_rdmol(self.input_context["protein"])
        merg = Chem.CombineMols(lig, prot)

        ligprops = rdForceFieldHelpers.MMFFGetMoleculeProperties(lig)
        ligff = rdForceFieldHelpers.MMFFGetMoleculeForceField(lig, ligprops)
        ligff: rdForceField.ForceField
        lig_energy = ligff.CalcEnergy()

        protprops = rdForceFieldHelpers.MMFFGetMoleculeProperties(prot)
        protff = rdForceFieldHelpers.MMFFGetMoleculeForceField(prot, protprops)
        protff: rdForceField.ForceField
        prot_energy = protff.CalcEnergy()

        mergprops = rdForceFieldHelpers.MMFFGetMoleculeProperties(merg)
        mergff = rdForceFieldHelpers.MMFFGetMoleculeForceField(merg, mergprops)
        mergff: rdForceField.ForceField
        merg_energy = mergff.CalcEnergy()

        return {"interaction": merg_energy - lig_energy - prot_energy}