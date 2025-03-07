import random, string, os, sys
from contextlib import contextmanager

def generate_random_str(n: int):
    return "".join(random.choices(string.ascii_uppercase + string.digits, k=n))

def generate_name_in_dir(n: int, dir: str, post: str):
    while True:
        s = generate_random_str(n)
        if not os.path.exists(os.path.join(dir, s+post)):
            break

    return s

def generate_path_in_dir(n: int, dir: str, post: str):
    name = generate_name_in_dir(n, dir, post)
    path = os.path.join(dir, f"{name}{post}")

    return path

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
