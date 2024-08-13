import pytest

from rdkit import Chem

from molsberry.modules.rdkit.selectors import (
    RDKitMolWtSelector
)
from molsberry.core import MoleculeData
from molsberry.core import SMILESRep
from molsberry.modules.rdkit.representations import RDKitMolRep
from molsberry.core.data.collections import BatchedData

@pytest.fixture
def mwselector():
    return RDKitMolWtSelector(debug=True, save_output=False, max_wt=600)

@pytest.fixture
def input_ligands_moltype():
    smis = ["C1CC(=O)CCC1F", "COC(=O)CCCN", 
            "CCc1cccnc1CC(C)(CC(I)(I)CC(I)CCCCC)I", "C1CCCCC1"]
    mols = [Chem.MolFromSmiles(smi) for smi in smis]
    ligs = [MoleculeData(RDKitMolRep(mol)) for mol in mols]
    return BatchedData(ligs)

@pytest.fixture
def input_ligands_smitype():
    smis = ["C1CC(=O)CCC1F", "COC(=O)CCCN", 
            "CCc1cccnc1CC(C)(CC(I)(I)CC(I)CCCCC)I", "C1CCCCC1"]
    ligs = [MoleculeData(SMILESRep(smi)) for smi in smis]
    return BatchedData(ligs)

def test_rdkit_mwselector(mwselector, 
                          input_ligands_moltype,
                          input_ligands_smitype):
    output = mwselector.execute(input_ligands_moltype)
    assert len(output["molecules"]) < len(input_ligands_moltype)
    
    output = mwselector.execute(input_ligands_smitype)
    assert len(output["molecules"]) < len(input_ligands_smitype)