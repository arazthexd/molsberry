import os
import pytest
import subprocess

from molsberry.core.pipeline import Pipeline, OutputBlock, InputBlock
from molsberry.core.data import SMILESRep, LigandData

from rdkit import Chem
from molsberry.modules.rdkit.representations import RDKitMolRep
from molsberry.modules.rdkit import (
    RDKitLigandTautEnumerator, RDKitLigandHAdder, RDKitLigandEmbedder
)

base_dir = os.path.dirname(__file__)

@pytest.fixture
def pipeline1_generator():
    class Pipeline1(Pipeline):
        name = "pipeline1"
        def build(self):
            self.add_block(InputBlock(["smiles"]), "input")
            self.add_block(RDKitLigandTautEnumerator(flatten=True), "tau_enum")
            self.add_block(RDKitLigandHAdder(), "h_adder")
            self.add_block(RDKitLigandEmbedder(), "lig_embedder")
            self.add_block(OutputBlock(["ligands"]), "output")
            self.add_connection("input", "smiles", "tau_enum", "ligands")
            self.add_connection("tau_enum", "ligands", "h_adder", "ligands")
            self.add_connection("h_adder", "ligands", "lig_embedder", "ligands")
            self.add_connection("lig_embedder", "ligands", "output", "ligands")
    return Pipeline1

def test_pipeline1(pipeline1_generator):
    pipeline: Pipeline = pipeline1_generator(base_dir=base_dir, 
                                             debug=True, save_output=True)
    smi_rep = SMILESRep("CCC(=O)CCc1ccccc1")
    output = pipeline.execute(smiles=LigandData(smi_rep))
    assert output["ligands"].depth == 1
    assert isinstance(
        output["ligands"][0].get_representation_content(RDKitMolRep), 
        Chem.Mol)
    assert len(output["ligands"]) > 1