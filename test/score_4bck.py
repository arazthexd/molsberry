import os
import wget
import pathlib
base_path = pathlib.Path(".").absolute()
print(base_path)

from rdkit import Chem

# os.chdir("..")
import sys
sys.path.append("..")
from src.abstract import *
from src.interfaces.rdkit import *
from src.interfaces.durrantlab import *
from src.interfaces.openmm import *
from src.interfaces.gnina import *
from src.interfaces.openmm import *
from src.interfaces.generic import *
from src.interfaces.cuby4 import *
from src.interfaces.cuby4c import *
from src.config import *
from src.utils import baselines, pdbtools
# os.chdir(base_path)

SQM_SCORE_MOPAC_SOLV_CONFIG = Cuby4MOPACFullConfig(
    method="pm6",
    mopac_exe="auto",
    mozyme=True,
    corrections="d3h4x",
    solvent="cosmo2",
)

SQM_SCORE_MOPAC_GAS_CONFIG = Cuby4MOPACFullConfig(
    method="pm6",
    mopac_exe="auto",
    mozyme=True,
    corrections="d3h4x",
    solvent="none",
)

class SQMBaselineScorePipe(Pipeline):
    name = "SQM Ligand Baseline Scoring"
    def __init__(self, return_input: bool = False):
        super().__init__(return_input=False, out_dir=OUT_DIR)
        self.rtin = return_input
        interface1 = Cuby4Interface(intconfig=SQM_SCORE_MOPAC_SOLV_CONFIG, cuby4_path="auto", 
                                    work_dir=TMP_DIR, debug=DEBUG_MODE)
        interface2 = Cuby4Interface(intconfig=SQM_SCORE_MOPAC_GAS_CONFIG, cuby4_path="auto", 
                                    work_dir=TMP_DIR, debug=DEBUG_MODE)
        self.blocks = [ # TODO: Other enumerations...
            RDKitTautEmumerator(debug=DEBUG_MODE, max_tautomers=4),
            # DimorphiteProtoEnumerator(debug=DEBUG_MODE, min_ph=7, max_ph=7.8),
            # RDKitRingEnumerator(debug=DEBUG_MODE, minimize=True, num_confs=30, max_per_ring=2, dist_threshold=0.5),
            # OpenMMLigandOptimizer(debug=DEBUG_MODE, forcefields=["amber14-all.xml"], work_dir=TMP_DIR),
            Cuby4LigandOptimizer(debug=DEBUG_MODE, interface=interface2, n_threads=N_AVAILABLE_THREADS, work_dir=TMP_DIR),
            Cuby4LigandEnergyScorer(debug=DEBUG_MODE, interface=interface1, n_threads=N_AVAILABLE_THREADS, work_dir=TMP_DIR),
        ]
    def run(self, ligands: List[Chem.Mol], targets: List[str], extra_info: dict = ...) -> Tuple[List[Chem.Mol] | List[str]]:
        n_primary_ligs = len(ligands)
        orig_ligs = ligands
        orig_targs = targets
        # [ligand.SetDoubleProp("end_state_energy", self.blocks[-1].score(ligand)) for ligand in ligands]
        [ligand.SetIntProp("SpecialIdx", i) for i, ligand in enumerate(ligands)]
        ligands, targets = super().run(ligands, targets, extra_info)
        prev_id = 0
        final_scores = []
        energies = []
        for ligand in ligands:
            # print(prev_id)
            # print(ligand.GetIntProp("SpecialIdx"))
            if ligand.GetIntProp("SpecialIdx") != prev_id:
                prev_id = ligand.GetIntProp("SpecialIdx")
                final_scores.append(baselines.average_energies(energies))
                energies = []
            fcharge = Chem.GetFormalCharge(ligand)
            energy = ligand.GetDoubleProp(Cuby4LigandEnergyScorer.score_name)
            htrans = baselines.calc_hydrogen_transfer_energy(fcharge)
            energies.append(energy+htrans)
        final_scores.append(baselines.average_energies(energies))
        
        # print(n_primary_ligs)
        # print(final_scores)
        assert n_primary_ligs == len(final_scores)
        extra_info["baseline_energy"] = final_scores

        if self.rtin:
            for i, lig in enumerate(orig_ligs):
                lig.SetDoubleProp("baseline_energy", final_scores[i])
            return orig_ligs, orig_targs
        else:
            return ligands, targets
        

class TestPipe(Pipeline):
    name = "Test Pipeline"
    def __init__(self):
        super().__init__()
        self.interface = Cuby4Interface(intconfig=SQM_SCORE_MOPAC_SOLV_CONFIG, cuby4_path="auto", 
                                        work_dir=TMP_DIR, debug=DEBUG_MODE)
        self.blocks = [
            # PDBFixerProteinPrepper(debug=DEBUG_MODE, ph=7.4, chains="all", keep_water=False, 
            #                        add_missing_residues=False, out_dir=TMP_DIR, save_prefix="prepared"),
            # PocketIsolator(debug=DEBUG_MODE, radius=10.0, work_dir=TMP_DIR),
            # Cuby4ComplexInteractScorer(debug=DEBUG_MODE, interface=self.interface, n_threads = N_AVAILABLE_THREADS, work_dir=TMP_DIR),
            SQMBaselineScorePipe(return_input=True),
            # Cuby4LigandEnergyScorer(score_name="hahahatest", debug=DEBUG_MODE, interface=self.interface, n_threads=N_AVAILABLE_THREADS, work_dir=TMP_DIR),
            # Cuby4LigandHTransScorer(debug=DEBUG_MODE, interface=self.interface, n_threads=N_AVAILABLE_THREADS, work_dir=TMP_DIR)
        ]

prot = "/home/arazthexd/projects/002_sqm/test/4bck_prot.pdb"
lig = Chem.MolFromMol2File("/home/arazthexd/projects/002_sqm/test/4bck_lig.mol2", removeHs=False)

pipeline = TestPipe()
output = pipeline.run([lig], [prot])

print(pipeline.extra_info)