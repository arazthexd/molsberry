import os
import pytest
import subprocess

from molsberry.core.pipeline import Pipeline, OutputBlock, InputBlock
from molsberry.core.data import SMILESRep, LigandData

try:
    from rdkit import Chem
    from molsberry.modules.rdkit.representations import RDKitMolRep
    from molsberry.modules.rdkit import (
        RDKitLigandTautEnumerator, RDKitLigandHAdder, RDKitLigandEmbedder
    )
    RDKIT_SUCCESSFUL_IMPORT = True
except: 
    import warnings
    warnings.warn("RDKit not loaded for basic pipelines.")
    RDKIT_SUCCESSFUL_IMPORT = False

try:
    assert subprocess.run(
        ["mopac", "test.mop"], capture_output=True).stderr.decode().strip() == \
            """MOPAC input data-set file: "test" does not exist."""
    from molsberry.modules.mopac.representations import MOPACInputMolRep
    from molsberry.modules.mopac import (
        MOPACLigandSinglePointCalculator, MOPACLigandOptimizer,
        MOPACConfig, MOPACMozymeConfig
    )
    MOPAC_SUCCESSFUL_IMPORT = True
except: 
    import warnings
    warnings.warn("Mopac not loaded for basic pipelines.")
    MOPAC_SUCCESSFUL_IMPORT = False

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

@pytest.fixture
def pipeline2_generator(pipeline1_generator):
    class Pipeline2(Pipeline):
        name = "pipeline2"
        def build(self):
            self.add_block(pipeline1_generator(base_dir=base_dir, debug=True), 
                           "3der")
            config_sp = MOPACConfig()
            self.add_block(
                MOPACLigandSinglePointCalculator(config_sp, debug=True), "sp1")
            config_opt = MOPACMozymeConfig()
            self.add_block(MOPACLigandOptimizer(
                config_opt, opt_algorithm="LBFGS", debug=True), "opt")
            self.add_block(
                MOPACLigandSinglePointCalculator(config_sp, debug=True), "sp2")
            self.add_block(OutputBlock(["ligands", "pre_e", "post_e"]), "output")
            self.add_connection("3der", "ligands", "sp1", "ligands")
            self.add_connection("sp1", "energy", "output", "pre_e")
            self.add_connection("3der", "ligands", "opt", "ligands")
            self.add_connection("opt", "ligands", "sp2", "ligands")
            self.add_connection("opt", "ligands", "output", "ligands")
            self.add_connection("sp2", "energy", "output", "post_e")
    return Pipeline2
            


@pytest.mark.skipif(not RDKIT_SUCCESSFUL_IMPORT, 
                    reason="RDKit Not Available Always")
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

@pytest.mark.skipif(not RDKIT_SUCCESSFUL_IMPORT or not MOPAC_SUCCESSFUL_IMPORT, 
                    reason="RDKit or MOPAC Not Available Always")
def test_pipeline2(pipeline2_generator):
    pipeline: Pipeline = pipeline2_generator(base_dir=base_dir, 
                                             debug=True, save=True)
    smi_rep = SMILESRep("CCC(=O)CCc1ccccc1")
    output = pipeline.execute(ligands=LigandData(smi_rep))
    assert output["ligands"].get_depth() == 0
    assert isinstance(output["ligands"][0].get_data(RDKitMolRep), Chem.Mol)