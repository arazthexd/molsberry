from rdkit import Chem

from ...core.data.special_cls import Ligand, Protein
from .representations import RDKitMolRep

def ligand_to_rdmol(ligand: Ligand):
    rdmol = Chem.Mol(ligand.get_data(RDKitMolRep))
    return rdmol

def protein_to_rdmol(protein: Protein):
    protmol = Chem.Mol(protein.get_data(RDKitMolRep))
    return protmol