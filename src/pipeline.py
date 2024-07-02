
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


class LigsPrepPipe(Pipeline):
    name = "Ligands Preparation"
    def __init__(self):
        super().__init__()
        self.blocks = [
            RDKitMWLigSelector(identifier="all", max_wt=600.0, min_wt=30.0),
            RDKitStereoEnumerator(debug=DEBUG_MODE),
            RDKitTautEmumerator(max_tautomers=1),
            DimorphiteProtoEnumerator(debug=DEBUG_MODE, min_ph=7.2, max_ph=7.6),
            RDKitRingEnumerator(minimize=True, num_confs=10, max_per_ring=1, dist_threshold=0.5)
        ]

class ProtPrepPipe(Pipeline):
    name = "Protein Preparation"
    def __init__(self):
        super().__init__()
        self.blocks = [
            PDBFixerProteinPrepper(ph=7.4, chains="all", keep_water=False, add_missing_residues=False,
                                   out_dir=OUT_DIR, save_prefix="prepared")
        ]

class DockingPipe(Pipeline):
    name = "Molecular Docking"
    def __init__(self):
        super().__init__()
        interface = GninaInterface(work_dir=TMP_DIR) # TODO: Some configurations are not included... 
        pocket_loc = PocketLocation("boxcxyz", center=(1.16, -0.56,- 2.91), xyz_size=(20, 20, 20))
        self.blocks = [
            GninaDocker(interface, pocket_loc, flexible_residues=None, batch_size=None),
            GninaResultsSelector(sort_by="CNN_VS", max_vina=-4.0, min_posescore=0.4, select_ratio=0.1,
                                 max_per_lig=1, n_select=5)
        ]

class MMPipe(Pipeline):
    name = "MM Prot Enum, Opt & Select"
    def __init__(self):
        super().__init__()
        # self.forcefield = ForceField("amber14-all.xml", "implicit/gbn2.xml")
        self.blocks = [
            DimorphiteProtoEnumerator(debug=DEBUG_MODE, min_ph=4, max_ph=10),
            PDBFixerProteinPrepper(ph=7.4, chains="all", keep_water=False, add_missing_residues=False,
                                   out_dir=TMP_DIR, save_prefix="prepared2"),
            OpenMMComplexOptimizer(forcefields=["amber14-all.xml", "implicit/gbn2.xml"], work_dir=TMP_DIR, jobname="mmopt"),
            OpenMMComplexInteractScorer(forcefields=["amber14-all.xml", "implicit/gbn2.xml"], work_dir=TMP_DIR, debug=DEBUG_MODE),
            GenericLigandSelector(sort_by=OpenMMComplexInteractScorer.score_name, identifier="_Name", max_per_id=2, reverse=False),
            GenericLigandSelector(sort_by=OpenMMComplexInteractScorer.score_name, identifier="all", max_per_id=2, reverse=False)
        ]

class QMMMOptimizePipe(Pipeline):
    name = "QMMM Optimization"
    def __init__(self):
        super().__init__()
        interface_config = Cuby4QMMMConfig(
            qm_config=Cuby4MOPACFullConfig(
                method="pm6",
                mopac_exe="auto",
                mozyme=True,
                corrections="d3h4x",
                solvent="none",
            ),
            mm_config=Cuby4AmberConfig(amber_home="auto", forcefield=ForceField("amber14-all.xml")),
            qmmm_embedding="mechanical",
            grad_on_points=False
        )
        interface = Cuby4Interface(intconfig=interface_config, cuby4_path="auto", 
                                   work_dir=TMP_DIR)
        self.blocks = [
            PocketIsolator(work_dir=TMP_DIR),
            Cuby4ComplexOptimizer(interface, n_threads=N_AVAILABLE_THREADS, work_dir=TMP_DIR)
        ]


class TestPipe(Pipeline):
    name = "test"
    def __init__(self):
        super().__init__()
        self.blocks = [LigsPrepPipe(), ProtPrepPipe(), DockingPipe(), MMPipe()]