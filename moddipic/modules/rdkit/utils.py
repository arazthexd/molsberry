from rdkit import Chem

from ...core.data.special_cls import Ligand
from .representations import RDKitMolRep

def ligand_to_rdmol(ligand: Ligand):
    rdmol = Chem.Mol(ligand.get_data(RDKitMolRep))
    return rdmol