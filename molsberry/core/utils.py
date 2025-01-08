import random, string, os, sys
from contextlib import contextmanager

def generate_random_str(n: int):
    return "".join(random.choices(string.ascii_uppercase + string.digits, k=n))

@contextmanager
def suppress_prints():
    original_stdout = sys.stdout
    original_stderr = sys.stderr
    sys.stdout = open(os.devnull, 'w')
    sys.stderr = open(os.devnull, 'w')
    try:
        yield
    finally:
        sys.stdout.close()
        sys.stderr.close()
        sys.stdout = original_stdout
        sys.stderr = original_stderr
