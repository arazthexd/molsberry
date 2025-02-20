import os
from os.path import join, exists, dirname

from rdkit import Chem
from rdkit.Chem import AllChem

from openmm.app import PDBFile
from pdbfixer import PDBFixer

from molsberry.modules.generic import PocketLocation, RDKitPocketIsolator

self_dir = dirname(__file__)
processed_dir = join(self_dir, "processed")

# CREATE DIRECTORY FOR PROCESSED
if not exists(processed_dir):
    os.mkdir(processed_dir)

# KGUD SEPARATE PROTEIN
fixer = PDBFixer(join(self_dir, "kguD.pdb"))
fixer.removeHeterogens()
fixer.addMissingHydrogens()
PDBFile.writeFile(fixer.topology, fixer.positions, join(processed_dir, 
                                                        "kguD_prot.pdb"))

# KGUD SEPARATE LIG
rdcomp = Chem.MolFromPDBFile(join(self_dir, "kguD.pdb"))
smi = ("C1=CC(=C[N+](=C1)[C@H]2[C@@H]([C@@H]([C@H](O2)COP(=O)(O)OP(=O)(O)O"
       "C[C@@H]3[C@H]([C@H]([C@@H](O3)N4C=NC5=C(N=CN=C54)N)O)O)O)O)C(=O)N")
for submol in Chem.GetMolFrags(rdcomp, asMols=True):
    submol: Chem.Mol
    if submol.GetAtomWithIdx(0).GetPDBResidueInfo().GetIsHeteroAtom():
        try:
            submol = AllChem.AssignBondOrdersFromTemplate(
                Chem.MolFromSmiles(smi), submol
            )
        except:
            continue

        lig = Chem.AddHs(submol, addCoords=True, addResidueInfo=True)
        Chem.SDWriter(join(processed_dir, "kguD_lig.sdf")).write(lig)
        break

# KGUD POCKET FROM LIG (NAD)
loc = PocketLocation(method="ligand", ligand=lig, radius=8)
rdprot = Chem.MolFromPDBFile(join(processed_dir, "kguD_prot.pdb"), 
                             removeHs=False)
rdpoc = RDKitPocketIsolator().isolate(rdprot, loc)
Chem.MolToPDBFile(rdpoc, join(processed_dir, "kguD_poc.pdb"))

