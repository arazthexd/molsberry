
from typing import List, Any, Tuple

from rdkit import Chem

from .abstract import *
from .interfaces.rdkit import *
from .interfaces.durrantlab import *
from .interfaces.openmm import *
from .interfaces.gnina import *
from .interfaces.openmm import *
from .interfaces.generic import *
from .interfaces.cuby4 import *
from .interfaces.cuby4c import *
from .config import *
from .utils import baselines


class LigsPrepPipe(Pipeline):
    name = "Ligands Preparation"
    def __init__(self):
        super().__init__(out_dir=OUT_DIR)
        self.blocks = [
            RDKitMWLigSelector(debug=DEBUG_MODE, identifier="all", max_wt=600.0, min_wt=30.0),
            RDKitStereoEnumerator(debug=DEBUG_MODE),
            RDKitTautEmumerator(debug=DEBUG_MODE, max_tautomers=1),
            DimorphiteProtoEnumerator(debug=DEBUG_MODE, min_ph=7.2, max_ph=7.6),
            RDKitRingEnumerator(debug=DEBUG_MODE, minimize=True, num_confs=10, max_per_ring=1, dist_threshold=0.5)
        ]

class ProtPrepPipe(Pipeline):
    name = "Protein Preparation"
    def __init__(self):
        super().__init__(out_dir=OUT_DIR)
        self.blocks = [
            PDBFixerProteinPrepper(debug=DEBUG_MODE, ph=7.4, chains="all", keep_water=False, 
                                   add_missing_residues=False, out_dir=TMP_DIR, save_prefix="prepared")
        ]

class DockingPipe(Pipeline):
    name = "Molecular Docking"
    def __init__(self):
        super().__init__(out_dir=OUT_DIR)
        interface = GninaInterface(work_dir=TMP_DIR) # TODO: Some configurations are not included... 
        pocket_loc = PocketLocation("boxcxyz", center=(1.16, -0.56,- 2.91), xyz_size=(20, 20, 20))
        self.blocks = [
            GninaDocker(debug=DEBUG_MODE, interface=interface, pocket_loc=pocket_loc, flexible_residues=None, 
                        batch_size=None),
            GninaResultsSelector(debug=DEBUG_MODE, sort_by="CNN_VS", max_vina=-4.0, min_posescore=0.4, 
                                 select_ratio=0.1, max_per_lig=1, n_select=5)
        ]

class PostDockEnumPipe(Pipeline):
    name = "Post Docking Enumeration"
    def __init__(self):
        super().__init__(out_dir=OUT_DIR)
        self.blocks = [
            DimorphiteProtoEnumerator(debug=DEBUG_MODE, min_ph=4, max_ph=10),
        ]
    
class MMPipe(Pipeline):
    name = "MM Optimization & Selection"
    def __init__(self):
        super().__init__(out_dir=OUT_DIR)
        # self.forcefield = ForceField("amber14-all.xml", "implicit/gbn2.xml")
        self.blocks = [
            PDBFixerProteinPrepper(debug=DEBUG_MODE, ph=7.4, chains="all", keep_water=False, add_missing_residues=False,
                                   out_dir=TMP_DIR, save_prefix="prepared2"),
            OpenMMComplexOptimizer(debug=DEBUG_MODE, forcefields=["amber14-all.xml", "implicit/gbn2.xml"], work_dir=TMP_DIR, jobname="mmopt"),
            OpenMMComplexInteractScorer(debug=DEBUG_MODE, forcefields=["amber14-all.xml", "implicit/gbn2.xml"], work_dir=TMP_DIR),
            GenericLigandSelector(debug=DEBUG_MODE, sort_by=OpenMMComplexInteractScorer.score_name, identifier="_Name", max_per_id=2, reverse=False),
            GenericLigandSelector(debug=DEBUG_MODE, sort_by=OpenMMComplexInteractScorer.score_name, identifier="all", max_per_id=2, reverse=False)
        ]

SQM_SCORE_MOPAC_GAS_CONFIG = Cuby4MOPACFullConfig(
    method="pm6",
    mopac_exe="auto",
    mozyme=True,
    corrections="d3h4x",
    solvent="none",
)

SQM_SCORE_MOPAC_SOLV_CONFIG = Cuby4MOPACFullConfig(
    method="pm6",
    mopac_exe="auto",
    mozyme=True,
    corrections="d3h4x",
    solvent="cosmo2",
)

class QMMMOptimizePipe(Pipeline):
    name = "QMMM Optimization"
    def __init__(self):
        super().__init__(out_dir=OUT_DIR)
        interface_config = Cuby4QMMMConfig(
            qm_config=SQM_SCORE_MOPAC_GAS_CONFIG,
            mm_config=Cuby4AmberConfig(amber_home="auto", forcefields=["amber14-all.xml"]),
            qmmm_embedding="mechanical",
            grad_on_points=False
        )
        interface = Cuby4Interface(intconfig=interface_config, cuby4_path="auto", 
                                   work_dir=TMP_DIR, debug=DEBUG_MODE)
        self.blocks = [
            PocketIsolator(debug=DEBUG_MODE, radius=10.0, work_dir=TMP_DIR),
            Cuby4ComplexOptimizer(debug=DEBUG_MODE, interface=interface, n_threads=N_AVAILABLE_THREADS, work_dir=TMP_DIR)
        ]

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
            RDKitTautEmumerator(debug=DEBUG_MODE, max_tautomers=100),
            DimorphiteProtoEnumerator(debug=DEBUG_MODE, min_ph=4.4, max_ph=10.4),
            RDKitRingEnumerator(debug=DEBUG_MODE, minimize=True, num_confs=50, max_per_ring=2, dist_threshold=0.5),
            OpenMMLigandOptimizer(debug=DEBUG_MODE, forcefields=["amber14-all.xml"], work_dir=TMP_DIR),
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

class SQMScorePipe(Pipeline):
    name = "SQM Scoring"
    def __init__(self):
        super().__init__()
        interface = Cuby4Interface(intconfig=SQM_SCORE_MOPAC_SOLV_CONFIG, cuby4_path="auto", 
                                   work_dir=TMP_DIR, debug=DEBUG_MODE)
        self.blocks = [
            Cuby4ComplexInteractScorer(debug=DEBUG_MODE, interface=interface, n_threads = N_AVAILABLE_THREADS, work_dir=TMP_DIR),
            SQMBaselineScorePipe(return_input=True),
            Cuby4LigandEnergyScorer(score_name="hahahatest", debug=DEBUG_MODE, interface=interface, n_threads=N_AVAILABLE_THREADS, work_dir=TMP_DIR),
            Cuby4LigandHTransScorer(debug=DEBUG_MODE, interface=interface, n_threads=N_AVAILABLE_THREADS, work_dir=TMP_DIR)
        ]   
    
    def run(self, ligands: List[Chem.Mol], targets: List[str], extra_info: dict = ...) -> Tuple[List[Chem.Mol] | List[str]]:
        ligands, targets = super().run(ligands, targets, extra_info)
        for ligand in ligands:
            final_score = ligand.GetDoubleProp(Cuby4ComplexInteractScorer.score_name) + \
                          ligand.GetDoubleProp(Cuby4LigandEnergyScorer.score_name) + \
                          ligand.GetDoubleProp(Cuby4LigandHTransScorer.score_name) - \
                          ligand.GetDoubleProp("baseline_energy")
            ligand.SetDoubleProp("final_energy", final_score)

        prefix = os.path.join(OUT_DIR, "".join(self.name.split()) + "Output")
        self.save_mid(ligands, targets, prefix+"Ligs.sdf", prefix+"Targs.txt")    
        return ligands, targets    

    #     ligstosave = []
    #     for ligand, baseline_score in zip(ligands, self.extra_info["baseline_energy"]):
    #         ligand.SetDoubleProp("baseline_energy", baseline_score)
        
    #     prefix = os.path.join(self.out_dir, 
    #                             "".join(self.name.split()) + "Output")
    #     self.save_mid(ligands, targets, prefix+"Ligs.sdf", prefix+"Targs.txt")
    #     return ligands, targets

class FinalPipe(Pipeline):
    name = "SQM Virtual Screening"
    def __init__(self):
        super().__init__(out_dir=OUT_DIR)
        self.blocks = [LigsPrepPipe(), ProtPrepPipe(), DockingPipe(), MMPipe(), 
                       QMMMOptimizePipe(), SQMScorePipe()]




class TestPipe(Pipeline):
    name = "test"
    def __init__(self):
        super().__init__()
        self.blocks = [LigsPrepPipe(), ProtPrepPipe(), DockingPipe(), MMPipe(),
                       QMMMOptimizePipe(), SQMBaselineScorePipe(return_input=False)]

class TestPipe2(Pipeline):
    name = "test"
    def __init__(self):
        super().__init__()
        self.blocks = [SQMBaselineScorePipe()]


class PLREXEvalPipe(Pipeline):
    name = "PL-REX Database Scoring Pipeline"
    def __init__(self):
        super().__init__()
        