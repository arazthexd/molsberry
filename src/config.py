##############################################################################
###########                        GENERAL                         ###########
##############################################################################

DEBUG_MODE = True

N_AVAILABLE_THREADS = 4

# Directories & Names
TMP_DIR = "./tmp"
OUT_DIR = "./output"
JOB_PREFIX = "kguD"

# Environment
ENVIRONMENT_PH = 7.4

# Inputs
INPUT_PROTEIN = "data/targets/kguD.pdb"
INPUT_LIGANDS = "data/ligands/tmp/test_ligands.smi"

##############################################################################
###########                     PREPARATION                        ###########
##############################################################################

# Protein Preparation
PROTEIN_PREPARATION_PERFORM = True
PROTEIN_CHAINS = "all"
PROTEIN_PH = ENVIRONMENT_PH
PROTEIN_KEEP_WATER = False
PROTEIN_ADD_MISSING_RESIDUES = False

# Ligand Preparation
LIGANDS_PREPARATION_PERFORM = True
LIGANDS_PH = ENVIRONMENT_PH
LIGANDS_MAX_WEIGHT = 600

# Enumeration
LIGANDS_ENUMERATE_TAUTOMERS = True
LIGANDS_ENUMERATE_TAUTOMERS_MAX = 4

LIGANDS_ENUMERATE_STEREOISOMERS = True

LIGANDS_ENUMERATE_RINGS = True
LIGANDS_ENUMERATE_RINGS_MAX_PER_RING = 2
LIGANDS_ENUMERATE_RINGS_MAX_CONFS = 50
LIGANDS_ENUMERATE_RINGS_MINIMIZE = True
LIGANDS_ENUMERATE_RINGS_DIFF_THRESH = 0.5

LIGANDS_ENUMERATE_BEST_PROTOMERS = True

# Save
PRE_DOCKING_LIGANDS_PREFIX = "prepared" # None = Don't Save
PREPARED_PROTEIN_PREFIX = "prepared"

##############################################################################
###########                      DOCKING (VS)                      ###########
##############################################################################

DOCKING_PERFORM = True
DOCKING_BATCH_SIZE = 5
DOCKING_NUM_CPU = N_AVAILABLE_THREADS

# Engine
DOCKING_ENGINE = 1 # (1) GNINA

# Pocket Definition
PROTEIN_POCKET_TYPE = 1 # (1) Center, XYZ Size # TODO: Others...
PROTEIN_POCKET_CENTER = (1.16, -0.56,- 2.91)
PROTEIN_POCKET_SIZE = (20, 20, 20)
PROTEIN_POCKET_LIGAND = None # TODO: This should be an option as well.

# Selection
DOCKING_SELECT_RATIO = 0.1
DOCKING_SELECT_MAX_PER_LIG = 2

# Engine Specifics
## (1) GNINA
GNINA_BIN_PATH = "/home/arazthexd/tools/gnina/gnina"
GNINA_SORT_RESULTS_BY = "CNN_VS"
GNINA_SELECT_VINAAFFINITY_MAX = -4.0
GNINA_SELECT_CNNPOSESCORE_MIN = 0.4
GNINA_OTHER_PARAMS = ""
GNINA_EXHAUSTIVENESS = 10
GNINA_MODES_PER_LIG = 5

# Save
POST_DOCK_LIGANDS_PREFIX = "docked"

##############################################################################
###########                PROTOMER ENUMERATION                    ###########
##############################################################################

# Protomer Generation (Inplace)
PROTOMER_GENERATION_PERFORM = True
PROTOMER_GENERATION_ENGINE = 1 # (1) Dimorphite_DL
PROTOMER_GENERATION_PH = ENVIRONMENT_PH # Default
PROTOMER_GENERATION_PH_RANGE = 3

# Protomer Selection
PROTOMER_SELECTION_PERFORM = False
PROTOMER_SELECTION_MAX_PER_LIG = 1
PROTOMER_SELECTION_SORT_BY = None # TODO: What?

# Save
PROTOMERS_SELECTED_PREFIX = "protenumed"

##############################################################################
###########                MM FORCEFIELD RESCORING                 ###########
##############################################################################

MM_RESCORING_PERFORM = True

# Engine
MM_ENGINE = 1 # (1) OpenMM 

# Forcefields
MM_BASE_FORCEFIELDS = ["amber14-all.xml", "implicit/gbn2.xml"]
# MM_LIG_FORCEFIELDS = "smirnoff" # TODO: Implement other types...

# Optimization
MM_OPTIMIZE_PERFORM = True
# MM_OPTIMIZE_MODE = 1 # (1) Protein-Ligand (2) Pocket-Ligand (3) ...
# MM_OPTIMIZE_MAX_STEPS = 200

# Selection
MM_SELECTION_PERFORM = True
MM_SELECT_MAX_PER_LIG = 1
MM_SELECT_RATIO = 0.1
MM_SELECT_MAX_ENERGY = 0.0

MM_RESCORE_PREFIX = "mmrescored"

##############################################################################
###########                   QMMM OPTIMIZATION                    ###########
##############################################################################

QMMM_OPTIMIZE_PERFORM = True

# Engine
QMMM_OPTIMIZE_ENGINE = 1 # (1) Cuby4


# Engine Specifics
CUBY4_QMMM_MM_INTERFACE = 1 # (1) Amber
CUBY4_QMMM_QM_INTERFACE = 1 # (1) MOPAC

##############################################################################
###########               AUTOMATIC (DO NOT CHANGE)                ###########
##############################################################################

import os
from os import path




