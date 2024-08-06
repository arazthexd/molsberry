import pytest
from moddipic.core.data.data_types import Ligand
from moddipic.core.data.representations import SMILESRep

@pytest.fixture
def input_smiles():
    return "CCOCCC(=O)NC(C)(C)C"

def test_ligand_class_creation_from_smiles(input_smiles):
    ligand = Ligand.from_smiles(input_smiles)
    assert ligand.get_data(SMILESRep) == input_smiles
    assert list(ligand._representations.keys()) == ["smiles"]

