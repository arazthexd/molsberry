from typing import List

from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem import rdForceFieldHelpers

from ..core.operations import LigandEnumerator
from ..core.pipeline import PipelineBlock

class RDKitTautEmumerator(PipelineBlock, LigandEnumerator):
    name = "RDKIT Tautomer Enumerator"
    def __init__(self, max_tautomers: int = 4, debug: bool = False):
        super().__init__(debug)
        self.enumerator = rdMolStandardize.TautomerEnumerator()
        self.enumerator.SetRemoveBondStereo(False)
        self.enumerator.SetRemoveSp3Stereo(False)
        self.enumerator.SetMaxTautomers(max_tautomers)

    def enumerate(self, ligand: Chem.Mol) -> List[Chem.Mol]:
        ligand = Chem.RemoveHs(ligand, updateExplicitCount=True)
        ts: List[Chem.Mol] = list(self.enumerator.Enumerate(ligand))
        ts = sorted(ts, key=lambda t: self.enumerator.ScoreTautomer(t), reverse=True)
        ts = [Chem.AddHs(t, addCoords=True) for t in ts]
        [rdForceFieldHelpers.MMFFOptimizeMolecule(t) for t in ts]
        return ts
    