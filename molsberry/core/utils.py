import random, string

def generate_random_str(n: int):
    return "".join(random.choices(string.ascii_uppercase + string.digits, k=n))