import pytest
import os, pathlib, shutil

from rdkit import Chem # TODO: Should we delete this?
from molsberry.modules.rdkit import RDKitMolRep

from molsberry.core import Data, ProteinData, LigandData, BatchedData
from molsberry.core import SMILESRep, PDBPathRep
from molsberry.core.pipeline import Pipeline, OutputBlock, InputBlock

from molsberry.modules.mopac.singlepoint import (
    MOPACSinglePointCalculator,
    MOPACPLInteractionCalculator
)
from molsberry.modules.mopac.representations import MOPACInputMolRep
from molsberry.modules.mopac.configs import MOPACConfig, MOPACMozymeConfig

@pytest.fixture(params=[MOPACConfig(), MOPACMozymeConfig()])
def general_spc(request):
    config = MOPACConfig()
    block = MOPACSinglePointCalculator(config=request.param,
                                        debug=True, save_output=False)
    
    base_dir = os.path.join(os.path.dirname(__file__), "out")
    if os.path.exists(base_dir):
        shutil.rmtree(base_dir)
    os.mkdir(base_dir)
    
    class JustAPipeline(Pipeline):
        name = "pipe"
        def build(self):
            self.add_block(OutputBlock(["out"]), "output")
            self.add_block(InputBlock(["inp"]), "input")

    block._parent = JustAPipeline(base_dir=base_dir)

    return block

@pytest.fixture(params=[MOPACMozymeConfig()])
def interaction_calculator(request):
    block = MOPACPLInteractionCalculator(config=request.param,
                                         debug=True, save_output=False)
    
    base_dir = os.path.join(os.path.dirname(__file__), "out")
    if os.path.exists(base_dir):
        shutil.rmtree(base_dir)
    os.mkdir(base_dir)
    
    class JustAPipeline(Pipeline):
        name = "pipe"
        def build(self):
            self.add_block(OutputBlock(["out"]), "output")
            self.add_block(InputBlock(["inp"]), "input")

    block._parent = JustAPipeline(base_dir=base_dir)
    block.calculator._parent = JustAPipeline(base_dir=base_dir)
    return block

@pytest.fixture
def sample_mols():
    sdf_path = os.path.join(os.path.dirname(__file__), 
                            "files", "sample_mols.sdf")
    rdmols = list(Chem.SDMolSupplier(sdf_path))
    rdreps = [RDKitMolRep(rdmol) for rdmol in rdmols]
    # mopacreps = [MOPACInputMolRep.from_RDKitMolRep(rdrep) for rdrep in rdreps]
    return BatchedData([LigandData(rep) for rep in rdreps])

@pytest.fixture
def interaction_lig():
    sdf_path = os.path.join(os.path.dirname(__file__), 
                            "files", "lig.sdf")
    rdmol = next(Chem.SDMolSupplier(sdf_path, removeHs=False))
    return LigandData(RDKitMolRep(rdmol))

@pytest.fixture
def interaction_prot():
    pdb_path = os.path.join(os.path.dirname(__file__), 
                            "files", "poc.pdb")
    return ProteinData(PDBPathRep(pdb_path))

def test_mopac_lig_spc(general_spc, sample_mols):
    output = general_spc.execute(sample_mols)
    assert len(output["energy"]) == len(sample_mols)
    assert issubclass(output["energy"].basic_itype, Data)

@pytest.mark.slow
def test_mopac_interaction(interaction_calculator, 
                           interaction_lig, interaction_prot):
    output = interaction_calculator.execute(interaction_lig, interaction_prot)
    # TODO: Complete this test...



