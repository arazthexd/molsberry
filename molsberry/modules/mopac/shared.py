# import os
# from ...global_conf import make_sure_exists, MAIN_DIR, TMP_DIR

# MOPAC_OUTPUT_FOLDERNAME = "mopac_out"
# MOPAC_OUTPUT_DIR = os.path.join(MAIN_DIR, MOPAC_OUTPUT_FOLDERNAME)
# make_sure_exists(MOPAC_OUTPUT_DIR)

# MOPAC_TMP_FOLDERNAME = "mopac"
# MOPAC_TMP_DIR = os.path.join(TMP_DIR, MOPAC_TMP_FOLDERNAME)
# make_sure_exists(MOPAC_TMP_DIR)

MOPAC_OPTIMIZE_KEYWORDS = ["BFGS", "LBFGS", "DFP", "EF"]
MOPAC_SINGLEPOINT_KEYWORDS = ["NOOPT"]