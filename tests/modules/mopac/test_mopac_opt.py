import pytest
import os, pathlib, shutil

from rdkit import Chem
from molsberry.modules.rdkit import RDKitMolRep

from molsberry.core import Data, ProteinData, LigandData, BatchedData
from molsberry.core import SMILESRep, PDBPathRep
from molsberry.core import Pipeline, OutputBlock, InputBlock

from molsberry.modules.mopac import MOPACLigandOptimizer
from molsberry.modules.mopac.representations import MOPACInputMolRep
from molsberry.modules.mopac.configs import MOPACConfig, MOPACMozymeConfig

ONETHREAD_CONFIG = MOPACMozymeConfig()
ONETHREAD_CONFIG.keywords.append("THREADS=1")
FOURTHREAD_CONFIG = MOPACMozymeConfig()
FOURTHREAD_CONFIG.keywords.append("THREADS=4")
@pytest.fixture(params=[MOPACConfig(), MOPACMozymeConfig(), 
                        ONETHREAD_CONFIG, FOURTHREAD_CONFIG])
def lig_opter(request):
    block = MOPACLigandOptimizer(config=request.param, opt_algorithm="LBFGS",
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

@pytest.fixture
def sample_mols():
    sdf_path = os.path.join(os.path.dirname(__file__), 
                            "files", "sample_mols_mk.sdf")
    rdmols = list(Chem.SDMolSupplier(sdf_path))[:10]
    rdreps = [RDKitMolRep(rdmol) for rdmol in rdmols]
    return BatchedData([LigandData(rep) for rep in rdreps])

def test_mopac_lig_opt(lig_opter, sample_mols):
    output = lig_opter.execute(sample_mols)
    assert len(output["ligands"]) == len(sample_mols)
    assert issubclass(output["ligands"].basic_itype, Data)