import os

def make_sure_exists(path):
    if not os.path.exists(path):
        os.mkdir(path)

MAIN_DIR = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 
    "out"
)
make_sure_exists(MAIN_DIR)

TMP_FOLDERNAME = "tmp"
TMP_DIR = os.path.join(MAIN_DIR, TMP_FOLDERNAME)
make_sure_exists(TMP_DIR)

RANDOM_JOB_KEY_LEN = 6