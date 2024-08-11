from rdkit import Chem
from ...core.data import Data
from .representations import RDKitMolRep

def data_to_rdmol(data: Data):
    rdmol = Chem.Mol(data.get_representation_content(RDKitMolRep))
    return rdmol