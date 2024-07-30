from typing import List

from rdkit.Chem import rdMolDescriptors

from ...core.templates import LigandAnalyzerBlock
from ...core.data.special_cls import Ligand
from .representations import RDKitMolRep
from .utils import ligand_to_rdmol

class RDKitMWCalculator(LigandAnalyzerBlock):
    name = "RDKit Molecular Weight Calculator"
    output_keys = ["molwt"]
    output_types = [float]

    def analyze(self, ligand: Ligand) -> float:
        rdmol = ligand_to_rdmol(ligand)
        return {"molwt": rdMolDescriptors.CalcExactMolWt(rdmol)}