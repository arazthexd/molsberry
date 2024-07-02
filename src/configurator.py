import os
from os import path

from .prepare import DefaultLigSetPrepper, DefaultProteinPrepper
from .interfaces.rdkit import (
    RDKitTautEmumerator, RDKitRingEnumerator, RDKitStereoEnumerator
)
from .interfaces.durrantlab import DimorphiteProtoEnumerator
from .interfaces.durrantlab import DimorphiteProtoEnumerator
from .utils.pockets import PocketLocation

from .config import *

def general_configurer():
    if not path.exists(OUT_DIR):
        os.mkdir(OUT_DIR)
    if not path.exists(TMP_DIR):
        os.mkdir(TMP_DIR)

def create_pocket_location():
    if PROTEIN_POCKET_TYPE == 1:
        PROTEIN_POCKET_LOCATION = PocketLocation(
            "boxcxyz", 
            center=PROTEIN_POCKET_CENTER, 
            xyz_size=PROTEIN_POCKET_SIZE
        )
    else:
        raise NotImplementedError()

    return PROTEIN_POCKET_LOCATION

def create_protein_preparator():
    PROTEIN_PREPPER = DefaultProteinPrepper(
        ph=PROTEIN_PH, chains=PROTEIN_CHAINS, keep_water=PROTEIN_KEEP_WATER,
        add_missing_residues=PROTEIN_ADD_MISSING_RESIDUES, out_dir=OUT_DIR,
        save_prefix=PREPARED_PROTEIN_PREFIX
    )
    return PROTEIN_PREPPER

# Run general config needed for pipeline...
general_configurer()

# Create pocket location
PROTEIN_POCKET_LOCATION = create_pocket_location()

# Create protein preparator
PROTEIN_PREPPER = create_protein_preparator()

#

if PROTEIN_POCKET_TYPE == 1:
    PROTEIN_POCKET_LOCATION = PocketLocation(
        "boxcxyz", 
        center=PROTEIN_POCKET_CENTER, 
        xyz_size=PROTEIN_POCKET_SIZE
    )
else:
    raise NotImplementedError()

# Creating protein preparator instance
PROTEIN_PREPPER = DefaultProteinPrepper(
    ph=PROTEIN_PH, chains=PROTEIN_CHAINS, keep_water=PROTEIN_KEEP_WATER,
    add_missing_residues=PROTEIN_ADD_MISSING_RESIDUES, out_dir=OUT_DIR,
    save_prefix=PREPARED_PROTEIN_PREFIX
)

# Creating ligand preparator instance
tautomer_enumerator = RDKitTautEmumerator() if LIGANDS_ENUMERATE_TAUTOMERS \
    else None
stereo_enumerator = RDKitStereoEnumerator() if LIGANDS_ENUMERATE_STEREOISOMERS \
    else None
protomer_enumerator = DimorphiteProtoEnumerator() if LIGANDS_ENUMERATE_BEST_PROTOMERS \
    else None
ring_enumerator = RDKitRingEnumerator(
    num_confs=LIGANDS_ENUMERATE_RINGS_MAX_CONFS,
    minimize=LIGANDS_ENUMERATE_RINGS_MINIMIZE,
    max_per_ring=LIGANDS_ENUMERATE_RINGS_MAX_PER_RING,
    dist_threshold=LIGANDS_ENUMERATE_RINGS_DIFF_THRESH
) if LIGANDS_ENUMERATE_RINGS else None

LIGAND_PREPPER = DefaultLigSetPrepper(
    tau_enum=tautomer_enumerator,
    stereo_enum=stereo_enumerator,
    prot_enum=protomer_enumerator,
    ring_enum=ring_enumerator,
    debug=DEBUG_MODE,
    ph=LIGANDS_PH,
    max_wt=LIGANDS_MAX_WEIGHT,
    num_jobs=N_AVAILABLE_THREADS
)

# Create docking performer instance
if DOCKING_ENGINE == 1:
    from .interfaces.gnina import GninaDocker, GninaInterface, GninaResultsSelector
    interface = GninaInterface(GNINA_BIN_PATH, work_dir=TMP_DIR)
    interface.set_options(
        exhaustiveness=GNINA_EXHAUSTIVENESS,
        num_modes=GNINA_MODES_PER_LIG,
        cpu=DOCKING_NUM_CPU,
        others=GNINA_OTHER_PARAMS
    )
    DOCKING_PERFORMER = GninaDocker(interface)
    DOCKING_SELECTOR = GninaResultsSelector(
        sort_by=GNINA_SORT_RESULTS_BY,
        max_vina=GNINA_SELECT_VINAAFFINITY_MAX,
        min_posescore=GNINA_SELECT_CNNPOSESCORE_MIN,
        select_ratio=DOCKING_SELECT_RATIO,
        max_per_lig=DOCKING_SELECT_MAX_PER_LIG
    )
else:
    raise NotImplementedError()

# Create protomers
if PROTOMER_GENERATION_ENGINE == 1:
    PROTOMER_GENERATION_ENUMERATOR = DimorphiteProtoEnumerator(debug=DEBUG_MODE)
else:
    raise NotImplementedError()

if MM_RESCORING_PERFORM:
    if MM_ENGINE == 1:
        from .interfaces.openmm import OpenMMComplexInteractScorer
        MM_RESCORER = OpenMMComplexInteractScorer(MM_BASE_FORCEFIELDS, TMP_DIR, 
                                                  optimize_before=MM_OPTIMIZE_PERFORM, debug=DEBUG_MODE)
    else:
        raise NotImplementedError
else:
    raise NotImplementedError()