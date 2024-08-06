from typing import Any, Dict, List

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, rdForceFieldHelpers
from rdkit.ForceField import rdForceField

from ...core.templates import LigandAnalyzerBlock
from ...core.templates.contexted import Contexted
from ...core.data.special_cls import Ligand, Protein
from .representations import RDKitMolRep
from .interface import RDKitInterface
from .utils import special_cls_to_rdmol

class RDKitMWCalculator(RDKitInterface, LigandAnalyzerBlock):
    name = "RDKit Molecular Weight Calculator"
    output_keys = ["molwt"]
    output_types = [float]

    def analyze(self, ligand: Ligand) -> Dict[str, float]:
        rdmol = special_cls_to_rdmol(ligand)
        return {"molwt": rdMolDescriptors.CalcExactMolWt(rdmol)}
 
class RDKitPLInteractionCalculator(RDKitInterface, Contexted, 
                                   LigandAnalyzerBlock):
    output_keys = ["interaction"]
    output_types = [float]
    input_context_keys = ["protein"]
    input_context_types = [Protein]

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