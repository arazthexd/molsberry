from openmm import unit, NonbondedForce
from openmm.app import ForceField, PDBFile

from ezdd.utils.protein.pockets import PocketLocation

class GenericPocket:
    def __init__(self, pdb_path: str, pocket_location: PocketLocation):
        self.pdb_path: str = pdb_path
        self.pocket_location: PocketLocation = pocket_location

def calc_protein_fcharge(pdb_path):

    pdb = PDBFile(pdb_path)
    ff = ForceField("amber14-all.xml")
    system = ff.createSystem(pdb.topology)
    nonbonded = [f for f in system.getForces() if isinstance(f, NonbondedForce)][0]

    total_charge = 0
    for i in range(system.getNumParticles()):
        charge, sigma, epsilon = nonbonded.getParticleParameters(i)
        total_charge += charge.value_in_unit(unit.elementary_charge)

    return round(total_charge)