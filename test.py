import warnings
warnings.simplefilter("ignore")

from rdkit import Chem

from src.pipeline import *
from src.utils.iotools import load_ligands, write_ligands
# pipe = TestPipe()
# pipe = MMPipe()
pipe = QMMMOptimizePipe()
# ligands = load_ligands("/home/arazthexd/projects/002_sqm/data/ligands/tmp/test_ligands.smi")
ligands = list(Chem.SDMolSupplier("test.sdf", removeHs=False))
# target = "/home/arazthexd/projects/002_sqm/data/targets/kguD.pdb"
target = "/home/arazthexd/projects/002_sqm/tmp/mmopt_prot_6E0ZI.pdb"
results = pipe.run(ligands, [target])
print(results[1])
print(pipe.extra_info)
write_ligands(results[0], "test2.sdf")