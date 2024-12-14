import pytest

from rdkit import Chem

from molsberry.modules.rdkit.selectors import (
    RDKitMolWtSelector,
    RDKitMolFlagSelector,
    RDKitMolLIPINSKISelector
)
from molsberry.core import MoleculeData
from molsberry.core import SMILESRep
from molsberry.modules.rdkit.representations import RDKitMolRep
from molsberry.core.data.collections import BatchedData

@pytest.fixture
def mwselector():
    return RDKitMolWtSelector(debug=True, save_output=False, max_wt=600)

@pytest.fixture
def flagselector():
    return RDKitMolFlagSelector(debug=True, save_output=False, filter_flag='PAINS_A')

@pytest.fixture
def lipinskiselector():
    return RDKitMolLIPINSKISelector(debug=True, save_output=False, condition_count=3)

@pytest.fixture
def input_ligands_moltype():
    smis = ["C1CC(=O)CCC1F", "COC(=O)CCCN", 
            "CCc1cccnc1CC(C)(CC(I)(I)CC(I)CCCCC)I", "C1CCCCC1", 
            "CC1=C(C=C(C=C1)N2C(=O)C(=C(N2)C)N=NC3=CC=CC(=C3O)C4=CC(=CC=C4)C(=O)O)C"]
    mols = [Chem.MolFromSmiles(smi) for smi in smis]
    ligs = [MoleculeData(RDKitMolRep(mol)) for mol in mols]
    return BatchedData(ligs)

@pytest.fixture
def input_ligands_smitype():
    smis = ["C1CC(=O)CCC1F", "COC(=O)CCCN", 
            "CCc1cccnc1CC(C)(CC(I)(I)CC(I)CCCCC)I", "C1CCCCC1",
            "CC1=C(C=C(C=C1)N2C(=O)C(=C(N2)C)N=NC3=CC=CC(=C3O)C4=CC(=CC=C4)C(=O)O)C"]
    ligs = [MoleculeData(SMILESRep(smi)) for smi in smis]
    return BatchedData(ligs)

def test_rdkit_flagselector(flagselector,
                            input_ligands_moltype,
                            input_ligands_smitype):
    output = flagselector.execute(input_ligands_moltype)
    assert len(output["molecules"]) < len(input_ligands_moltype)

    output = flagselector.execute(input_ligands_smitype)
    assert len(output["molecules"]) < len(input_ligands_smitype)

def test_rdkit_lipinski_selector(lipinskiselector,
                            input_ligands_moltype,
                            input_ligands_smitype):
    output = lipinskiselector.execute(input_ligands_moltype)
    assert len(output["molecules"]) < len(input_ligands_moltype)

    output = lipinskiselector.execute(input_ligands_smitype)
    assert len(output["molecules"]) < len(input_ligands_smitype)

def test_rdkit_mwselector(mwselector, 
                          input_ligands_moltype,
                          input_ligands_smitype):
    output = mwselector.execute(input_ligands_moltype)
    assert len(output["molecules"]) < len(input_ligands_moltype)
    
    output = mwselector.execute(input_ligands_smitype)
    assert len(output["molecules"]) < len(input_ligands_smitype)
