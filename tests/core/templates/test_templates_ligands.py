import pytest, random
from moddipic.core.templates import ligands
from moddipic.core.data.collections import Batched
from moddipic.core.data.special_cls import Ligand
from moddipic.core.templates.contexted import Contexted

@pytest.fixture
def ligconverter():
    class LigConverterBlock(ligands.LigandConverterBlock):
        def convert(self, ligand):
            return ligand
    return LigConverterBlock

@pytest.fixture
def contexedligconverter(): # TODO: Complete
    class LigConverterBlock(Contexted, ligands.LigandConverterBlock):
        input_context_keys = ["context"]
        input_context_types = [str]
        def __init__(self, debug: bool = False, save_output: bool = False):
            ligands.LigandConverterBlock.__init__(self, debug, save_output)
            Contexted.__init__(self)
        def convert(self, ligand):
            assert self.input_context
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
def liganalyzer():
    class LigAnalyzerBlock(ligands.LigandAnalyzerBlock):
        output_keys = ["random_number", "rand2"]
        output_types = [float, float]
        def analyze(self, ligand):
            return {"random_number": random.random(), "rand2": random.random()}
    return LigAnalyzerBlock

@pytest.fixture
def input_ligands():
    smiles_list = ["CCCOC", "CC", "O"]
    return Batched([Ligand.from_smiles(smi) for smi in smiles_list])

def test_ligand_converter_block(ligconverter, input_ligands):
    block = ligconverter()
    raw_newlig = block.convert(input_ligands.data[0])

    output = block.execute(input_ligands)

    with pytest.raises(AssertionError):
        output = block._auto_execute(input_ligands)
    assert len(output["ligands"]) == len(input_ligands)

def test_ligand_contexted_converter_block(contexedligconverter, input_ligands):
    block = contexedligconverter()
    assert set(block.input_context.keys()) == set(["context"])

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

def test_ligand_analyzer_block(liganalyzer, input_ligands):
    block = liganalyzer()
    raw_score = block.analyze(input_ligands.data[0])

    output = block.execute(input_ligands)

    with pytest.raises(AssertionError):
        output = block._auto_execute(input_ligands)
    assert len(output["random_number"]) == len(input_ligands)
    assert len(output["rand2"]) == len(input_ligands)

# TODO: Create tests for operations of selector, converter and enumerator.


