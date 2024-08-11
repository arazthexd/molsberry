import pytest
from molsberry.core.data.data_types import LigandData
from molsberry.core.data.representations import SMILESRep

@pytest.fixture
def input_smiles():
    return "CCOCCC(=O)NC(C)(C)C"

def test_ligand_class_creation_from_smiles(input_smiles):
    ligand = LigandData.from_smiles(input_smiles)
    assert ligand.get_representation_content(SMILESRep) == input_smiles
    assert list(ligand._representations.keys()) == ["smiles"]

