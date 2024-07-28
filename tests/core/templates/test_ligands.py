import pytest
from moddipic.core.templates import ligands
from moddipic.core.data.collections import Batched
from rdkit import Chem

@pytest.fixture
def ligconverter():
    class LigConverterBlock(ligands.LigandConverterBlock):
        def convert(self, ligand):
            return ligand
    return LigConverterBlock

@pytest.fixture
def ligenumerator():
    class LigEnumeratorBlock(ligands.LigandEnumeratorBlock):
        def enumerate(self, ligand):
            return [ligand, ligand]
    return LigEnumeratorBlock

@pytest.fixture
def ligselector():
    class LigSelectorBlock(ligands.LigandSelectorBlock):
        def select(self, ligands):
            return [ligands[0]]
    return LigSelectorBlock

@pytest.fixture
def input_ligands():
    return Batched([
        Chem.MolFromSmiles("CCCOC"),
        Chem.MolFromSmiles("CC"),
        Chem.MolFromSmiles("O")
    ])

def test_ligand_converter_block(ligconverter, input_ligands):
    block = ligconverter()
    raw_newlig = block.convert(input_ligands.data[0])

    output = block.execute(input_ligands)

    with pytest.raises(AssertionError):
        output = block._auto_execute(input_ligands)
    assert len(output["ligands"]) == len(input_ligands)

def test_ligand_enumerator_block(ligenumerator, input_ligands):
    block = ligenumerator()
    raw_newlig = block.enumerate(input_ligands.data[0])

    output = block.execute(input_ligands)
    with pytest.raises(AssertionError):
        output = block._auto_execute(input_ligands)

    block = ligenumerator(flatten=True)
    output = block.execute(input_ligands)
    assert len(output["ligands"]) == len(input_ligands) * 2

def test_ligand_selector_block(ligselector, input_ligands):
    block = ligselector()
    raw_newlig = block.select(input_ligands.data)
    assert len(raw_newlig) == 1

    output = block.execute(input_ligands)
    assert len(output["ligands"]) == 1

# TODO: Create tests for operations of selector, converter and enumerator.


