from typing import List, Tuple

import ash
from rdkit import Chem

from ...core import Representation
from ..rdkit import RDKitMolRep

class ASHFragmentRep(Representation):
    rep_name = "ash_inpmol"

    def __init__(self, frag: ash.Fragment, smiles: str = None) -> None:
        super().__init__(frag)
        self.content: ash.Fragment
        self.smiles = smiles
    
    @classmethod
    def from_RDKitMolRep(cls, rdkit_rep: RDKitMolRep):
        rdmol = Chem.Mol(rdkit_rep.content)
        if rdmol.GetNumConformers() > 0: # What about 2D conformer?
            pos = rdmol.GetConformer().GetPositions()
            elems = [atom.GetSymbol() for atom in rdmol.GetAtoms()]
            smi = None
        else:
            pos = None
            elems = None
            smi = Chem.MolToSmiles(rdmol)
        
        frag = ash.Fragment(smiles=smi, coords=pos, elems=elems,
                            mult=1, charge=Chem.GetFormalCharge(rdmol))
        return ASHFragmentRep(frag=frag, smiles=smi)

    def to_RDKitMolRep(self) -> RDKitMolRep:
        rdmol = Chem.MolFromSmiles(self.smiles)
        rdmol.AddConformer(Chem.Conformer(rdmol.GetNumAtoms()))
        rdmol.GetConformer().SetPositions(self.content.coords)
        return RDKitMolRep(rdmol)