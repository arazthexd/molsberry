import pytest

from moddipic.core.pipeline import PipelineBlock, Pipeline, OutputBlock
from moddipic.core.templates import (
    LigandConverterBlock, 
    LigandEnumeratorBlock
)
from moddipic.core.data.data_types import Ligand
from moddipic.core.data.representations import SMILESRep
from moddipic.core.data.collections import Batched

@pytest.fixture
def input_ligands_smitype():
    smis = ["CCC", "COC(=O)CCCN", "c1cccnc1CC(C)(C)Cl"]
    ligs = [Ligand(SMILESRep(smi)) for smi in smis]
    return Batched(ligs)

@pytest.fixture
def my_pipeline():
    
    class MyDoubler(LigandEnumeratorBlock):
        name = "my_enum"
        def enumerate(self, ligand):
            return [ligand, ligand]
    class MyEyer(LigandConverterBlock):
        name = "my_converter"
        def convert(self, ligand):
            return ligand
        
    class MyPipeline(Pipeline):
        name = "my_pipeline"
        def build(self):
            self.add_blocks({
                "b1": MyDoubler(flatten=True),
                "b2": MyEyer(),
                "output": OutputBlock(["ligands"])
            })
            self.add_connection("b1", "ligands", "b2", "ligands")
            self.add_connection("b2", "ligands", "output", "ligands")
    
    return MyPipeline()

def test_pipeline_exe(input_ligands_smitype, my_pipeline):
    my_pipeline.execute(ligands=input_ligands_smitype)
