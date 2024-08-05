from rdkit import Chem
from rdkit.Chem import rdDistGeom

from ...core.templates import LigandConverterBlock
from ...utils.moltools import addhs_based_on_confdim
from ...core.data.special_cls import Ligand
from .representations import RDKitMolRep
from .interface import RDKitInterface

class RDKitLigandHAdder(RDKitInterface, LigandConverterBlock):
    name = "RDKit Hydrogen Adder"

    def convert(self, ligand: Ligand) -> Ligand:
        rdmol = self.special_cls_to_rdmol(ligand)
        rdmol = addhs_based_on_confdim(rdmol)
        return Ligand(RDKitMolRep(rdmol))

class RDKitLigandEmbedder(RDKitInterface, LigandConverterBlock):
    name = "RDKit Ligand Embedder"

    def convert(self, ligand: Chem.Mol) -> Chem.Mol:
        rdmol = self.special_cls_to_rdmol(ligand)
        rdDistGeom.EmbedMolecule(rdmol)
        return Ligand(RDKitMolRep(rdmol))