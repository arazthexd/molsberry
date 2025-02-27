import pytest
import os, pathlib, shutil, warnings
from git import Repo
import wget

from rdkit import Chem
from rdkit.Chem import rdDistGeom

from openmm.app import PDBFile
from pdbfixer import PDBFixer

from molsberry.core import PDBPathRep, SMILESRep
from molsberry.modules.rdkit import RDKitMolRep

### INPUTS
# 1) SMALL MOLS
ZINC_LIG_URLS = [
    "https://zinc.docking.org/substances/ZINC000000000031.sdf"
]
ZINC_LIG_PATHS = [
    "./tests/data/ZINC/" + url.split("/")[-1] for url in ZINC_LIG_URLS
]

SMILES = [
    ("feat_neutr", "CCCO"),
    ("feat_n+", "CC(=O)OCC(C)C[NH3+]"),
    ("feat_o-", "CCC=CCC(=O)[O-]"),
    ("feat_no2", "c1ccc(O)cc1C[N+](=O)[O-]")
]

PLREX_LIGS = [
    ("PR-CA2-5NXG", "./data/PL-REX/001-CA2/structures_pl-rex/5NXG/ligand.sdf"),
    ("PR-BACE1-5QCR", "./data/PL-REX/006-BACE1/structures_pl-rex/5QCR/ligand.sdf")
]

# 2) PROTEINS
PLREX_PROTS = [
    ("PR-CA2-5NXG", "./data/PL-REX/001-CA2/structures_pl-rex/5NXG/protein.pdb"),
    ("PR-BACE1-5QCR", "./data/PL-REX/006-BACE1/structures_pl-rex/5QCR/protein.pdb")
]

# 3) COMPLEX
PLREX_COMPLEXES = list(zip(PLREX_LIGS, PLREX_PROTS))

### DOWNLOAD DATASETS
# 1) PL-REX
if not os.path.exists("./tests/data/PL-REX"):
    Repo.clone_from("https://github.com/Honza-R/PL-REX", "./tests/data/PL-REX")
    with open("./tests/data/PL-REX/.gitignore", "w") as f:
        f.write("*")

# 2) ZINC
if not os.path.exists("./tests/data/ZINC"):
    os.mkdir("./tests/data/ZINC")

for zinc_path, zinc_url in zip(ZINC_LIG_PATHS, ZINC_LIG_URLS):
    if not os.path.exists(zinc_path):
        wget.download(zinc_url, zinc_path)

### SET UP FIXTURES
def param_id_single(param):
    return param[0]

def param_id_pair(param):
    return param[0][0]+param[1][0]

def sdf2rdmol(sdf: str, name: str = "unnamed"):
    rdmol = next(Chem.SDMolSupplier(sdf, removeHs=False))
    rdmol = Chem.AddHs(rdmol, addCoords=True)
    rdmol.SetProp("_Name", name)
    return rdmol

# 1) SMALL MOLS
@pytest.fixture(params=SMILES, ids=param_id_single)
def sample_sm_rdrep_from_smi(request) -> RDKitMolRep:
    name, smiles = request.param
    rdmol = Chem.MolFromSmiles(smiles)
    rdmol = Chem.AddHs(rdmol)
    rdDistGeom.EmbedMolecule(rdmol, randomSeed=2025)
    rdmol.SetProp("_Name", name)
    rdrep = RDKitMolRep(rdmol)
    return rdrep

@pytest.fixture(params=PLREX_LIGS, ids=param_id_single)
def sample_sm_rdrep_from_plrex(request) -> RDKitMolRep:
    name, sdf = request.param
    rdmol = sdf2rdmol(sdf, name)
    rdrep = RDKitMolRep(rdmol)
    return rdrep

@pytest.fixture(params=ZINC_LIG_PATHS, ids=param_id_single)
def sample_sm_rdrep_from_zinc(request) -> RDKitMolRep:
    sdf = request.param
    name = os.path.basename(sdf).split(".")[0]
    rdmol = sdf2rdmol(sdf, name)
    rdrep = RDKitMolRep(rdmol)
    return rdrep

# 2) PROTEINS
@pytest.fixture(params=PLREX_PROTS, ids=param_id_single)
def sample_prot_pdbrep_from_plrex(request):
    return PDBPathRep(request.param[1])

# 3) COMPLEX
@pytest.fixture(params=PLREX_COMPLEXES, ids=param_id_pair)
def sample_complex_rdpdbrep_from_plrex(request):
    (lname, lsdf), (pname, ppdb) = request.param
    return RDKitMolRep(sdf2rdmol(lsdf, lname)), PDBPathRep(ppdb)