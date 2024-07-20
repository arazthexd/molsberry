import warnings
warnings.simplefilter("ignore")

from rdkit import Chem
from rdkit.Chem import rdDistGeom

from src.pipeline import *
from src.utils.iotools import load_ligands, write_ligands

pipe = FinalPipe()
# ligands = [next(Chem.SDMolSupplier("/home/arazthexd/projects/002_sqm/output/protenumed_ligands.sdf", removeHs=False))]
ligands = load_ligands("/home/arazthexd/projects/002_sqm/data/ligands/tmp/test_ligands.smi", final_3d=True, addHs=True)
targets = ["/home/arazthexd/projects/002_sqm/data/targets/kguD.pdb"]
results = pipe.run(ligands, targets)
print(pipe.extra_info)
write_ligands(results[0], "test2.sdf")

# pipe = TestPipe2()
# # pipe = MMPipe()
# # pipe = SQMScorePipe()
# # ligands = load_ligands("/home/arazthexd/projects/002_sqm/data/ligands/tmp/test_ligands.smi")
# ligands = list(Chem.SDMolSupplier("/home/arazthexd/projects/002_sqm/output/protenumed_ligands.sdf", 
#                                   removeHs=False))
# # target = "/home/arazthexd/projects/002_sqm/data/targets/kguD.pdb"
# target = "/home/arazthexd/projects/002_sqm/output/prepared_protein.pdb"
# results = pipe.run(ligands, [target])
# print(results[1])
# print(pipe.extra_info)
# write_ligands(results[0], "test2.sdf")