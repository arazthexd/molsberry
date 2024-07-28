import pytest
from moddipic.core.templates import proteins
from moddipic.core.data.collections import Batched

@pytest.fixture
def protconverter():
    class ProtConverterBlock(proteins.ProteinConverterBlock):
        def convert(self, protein):
            return protein
    return ProtConverterBlock

@pytest.fixture
def protenumerator():
    class ProtEnumeratorBlock(proteins.ProteinEnumeratorBlock):
        def enumerate(self, protein):
            return [protein, protein]
    return ProtEnumeratorBlock

@pytest.fixture
def protselector():
    class ProtSelectorBlock(proteins.ProteinSelectorBlock):
        def select(self, proteins):
            return [proteins[0]]
    return ProtSelectorBlock

@pytest.fixture
def input_proteins():
    return Batched([
        "Path/To/1.pdb",
        "Path/To/2.pdb"
    ])

def test_ligand_converter_block(protconverter, input_proteins):
    block = protconverter()
    output = block.execute(input_proteins)

    with pytest.raises(AssertionError):
        output = block._auto_execute(input_proteins)
    assert len(output["protein"]) == len(input_proteins)

def test_ligand_enumerator_block(protenumerator, input_proteins):
    block = protenumerator()
    raw_newlig = block.enumerate(input_proteins.data[0])

    output = block.execute(input_proteins)
    with pytest.raises(AssertionError):
        output = block._auto_execute(input_proteins)

    block = protenumerator(flatten=True)
    output = block.execute(input_proteins)
    assert len(output["protein"]) == len(input_proteins) * 2

def test_ligand_selector_block(protselector, input_proteins):
    block = protselector()
    raw_newlig = block.select(input_proteins.data)
    assert len(raw_newlig) == 1

    output = block.execute(input_proteins)
    assert len(output["protein"]) == 1