from typing import Any, Dict, List

import os
import glob
from copy import deepcopy

from ...core.templates import LigandConverterBlock
from ...core.templates import Contexted
from ...core.data import Ligand, Protein
from ...global_conf import RANDOM_JOB_KEY_LENGTH

from .representations import MOPACInputMolRep
from .configs import MOPACConfig, MOPACMozymeConfig
from .utils import (
    ligand_to_mopacrep, generate_random_input_file, write_and_run_mopac
)
from .shared import MOPAC_OUTPUT_DIR, MOPAC_TMP_DIR

class MOPACOptimizer():
    pass

class MOPACLigandOptimizer():
    