

from rdkit import Chem
from rdkit.Chem import rdDistGeom

from ..core.templates import LigandConverterBlock
from ..utils.moltools import addhs_based_on_confdim

class RDKitLigandHAdder(LigandConverterBlock):
    name = "RDKit Hydrogen Adder"

    def convert(self, ligand: Chem.Mol) -> Chem.Mol:
        return addhs_based_on_confdim(ligand)

class RDKitLigandEmbedder(LigandConverterBlock):
    name = "RDKit Ligand Embedder"

    def convert(self, ligand: Chem.Mol) -> Chem.Mol:
        ligand = Chem.Mol(ligand)
        rdDistGeom.EmbedMolecule(ligand)
        return ligand