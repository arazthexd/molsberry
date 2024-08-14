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
block = MOPACLigandOptimizer(config=ONETHREAD_CONFIG, opt_algorithm="LBFGS",
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

sdf_path = "tests/modules/mopac/files/sample_mols_mk.sdf" 
rdmols = list(Chem.SDMolSupplier(sdf_path))[:4]
rdreps = [RDKitMolRep(rdmol) for rdmol in rdmols]
inputdata = BatchedData([LigandData(rep) for rep in rdreps])

output = block.execute(inputdata)