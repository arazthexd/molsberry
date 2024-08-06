from rdkit import Chem
from ...core.data.abstract import SpecialDataClass
from .representations import RDKitMolRep

def special_cls_to_rdmol(special_cls: SpecialDataClass):
    rdmol = Chem.Mol(special_cls.get_data(RDKitMolRep))
    return rdmol