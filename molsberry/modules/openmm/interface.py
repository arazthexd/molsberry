from openmm.app import ForceField
from openmm import VerletIntegrator, Context, State, unit

from ...core import LigandData, ProteinData, MoleculeData
from ...core import SimpleBlock, PLInteractionJob, Representation

from .representations import OpenMMInputMolRep

class OpenMMInterface:

    @staticmethod
    def rep_to_energy(rep: OpenMMInputMolRep) -> float:
        """Calculate the energy of a representation using OpenMM

        Args:
            rep (OpenMMInputMolRep): The representation to calculate energy of

        Returns:
            float: The energy of the representation
        """
        integrator = VerletIntegrator(0.002)
        context = Context(rep.system, integrator)
        context.setPositions(rep.positions)
        state: State = context.getState(getEnergy=True)
        penergy = state.getPotentialEnergy()\
            .value_in_unit(unit.kilocalorie_per_mole)
        return penergy
