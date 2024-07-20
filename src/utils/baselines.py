from typing import List
import math

from rdkit import Chem

# CONSTANTS
# LN_TO_LOG = 2.303
BETA = 1 / 2.479 / 4.18 # kcal/mol

H2O_G = -50.06
H3O_G = -33.03
HO_G = 150.12 # Top three are for mopac...
H2O_GS = -6.3
H3O_GS = -110.2
HO_GS = -105.0

def calc_hydrogen_transfer_energy(num_hs_added: int):
    if num_hs_added > 0:
        return num_hs_added * (HO_G + HO_GS - H2O_G - H2O_GS)
    else:
        return num_hs_added * (H3O_G + H3O_GS - H2O_G - H2O_GS)
    
def calculate_probabilities(energies: list):
    boltz_facs = []
    partition = 0.0
    for energy in energies:
        boltz_fac = math.exp(-BETA * energy)
        boltz_facs.append(boltz_fac)
        partition += boltz_fac
    probs = [b / partition for b in boltz_facs]
    return probs

def average_energies(energies: list):
    probs = calculate_probabilities(energies)
    print("probs", probs)
    print("energies", energies)
    avgenergy = sum([prob*energy for prob, energy in zip(probs, energies)])
    return avgenergy