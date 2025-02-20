import pytest
import os, pathlib, shutil

from rdkit import Chem
from rdkit.Chem import rdDistGeom

from openmm.app import PDBFile
from pdbfixer import PDBFixer

from molsberry.core import PDBPathRep
from molsberry.modules.rdkit import RDKitMolRep

# Input Small Molecule Samples

SMILES = [
    ("neutr", "CCCO"),
    ("n+", "CC(=O)OCC(C)C[NH3+]"),
    ("o-", "CCC=CCC(=O)[O-]"),
    ("no2", "c1ccc(O)cc1C[N+](=O)[O-]")
]

@pytest.fixture(params=SMILES)
def sample_sm_rdrep(request) -> RDKitMolRep:
    name, smiles = request.param
    rdmol = Chem.MolFromSmiles(smiles)
    rdmol = Chem.AddHs(rdmol)
    rdDistGeom.EmbedMolecule(rdmol)
    rdmol.SetProp("_Name", name)
    rdrep = RDKitMolRep(rdmol)
    return rdrep

# Input Protein Samples

PROT_PATHS = [
    pathlib.Path("./tests/data/processed/kguD_prot.pdb").absolute()
]

@pytest.fixture(params=PROT_PATHS)
def sample_prot_pdbrep(request):
    return PDBPathRep(request.param)

# Input Pocket Samples

POC_PATHS = [
    pathlib.Path("./tests/data/processed/kguD_poc.pdb").absolute()
]

@pytest.fixture(params=POC_PATHS)
def sample_poc_pdbrep(request):
    return PDBPathRep(request.param)