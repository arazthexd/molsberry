from __future__ import annotations

from typing import Any
from openmm import (
    System, Context, VerletIntegrator, Vec3, unit, 
    LocalEnergyMinimizer, State, NonbondedForce
)
from openmm.app import ForceField, PDBFile, Modeller, Topology

from ...core.data.abstract import Representation

try: 
    from rdkit import Chem
    from ...modules.rdkit.representations import RDKitMolRep
except: print("Warning: RDKit failed to load in openmm...")

try:
    from openff.toolkit import Molecule
    from openff.units import unit as offunit
    from openff.units.openmm import to_openmm as openff_unit_to_openmm
    try:
        from espaloma_charge.openff_wrapper import EspalomaChargeToolkitWrapper
        ESPALOMA_CHARGE_AVAILABLE = True
    except:
        ESPALOMA_CHARGE_AVAILABLE = False
        print("Warning: Espaloma charge failed to load in openff in openmm...")
    try:
        from openmmforcefields.generators import SMIRNOFFTemplateGenerator
    except: print("Warning: openmmffs failed to load in openff in openmm...")
except: print("Warning: OpenFF failed to load in openmm...")

class OpenMMComponent:
    def __init__(
            self, 
            topology: Topology, 
            system: System, 
            positions: Any, 
            forcefield: ForceField) -> None:
        assert isinstance(topology, Topology)
        assert isinstance(system, System)
        # assert isinstance(positions, ...) # TODO: Typing...
        assert isinstance(forcefield, ForceField)

        self.topology: Topology = topology
        self.system: System = system
        self.positions = positions # TODO: Typing...
        self.forcefield: ForceField = forcefield
    
    @classmethod
    def from_rdmol(cls, rdmol: Chem.Mol, 
                   forcefield: ForceField, 
                   cache_path: str = None) -> OpenMMComponent:
        rdmol = Chem.AddHs(rdmol)
        molecule = Molecule.from_rdkit(rdmol, allow_undefined_stereo=True,
                                       hydrogens_are_explicit=True)
        if ESPALOMA_CHARGE_AVAILABLE:
            molecule.assign_partial_charges(
                'espaloma-am1bcc', 
                toolkit_registry=EspalomaChargeToolkitWrapper(), 
                use_conformers=True)
        else:
            raise NotImplementedError()
        
        topology = molecule.to_topology().to_openmm()
        positions = openff_unit_to_openmm(molecule.conformers[0])

        smirnoff = SMIRNOFFTemplateGenerator(molecules=molecule, 
                                             cache=cache_path)
        forcefield.registerTemplateGenerator(smirnoff.generator)
        system = forcefield.createSystem(topology) # , *args, **kwargs)
                                                   # TODO: get other args

        return cls(topology=topology, 
                   system=system, 
                   positions=positions, 
                   forcefield=forcefield)

class OpenMMComponentRep(Representation):
    rep_name = "openmm_component"

    def __init__(self, component: OpenMMComponent):
        assert isinstance(component, OpenMMComponent)
        super().__init__(data=component)
    
    @classmethod
    def from_RDKitMolRep(cls, rdmol_rep: RDKitMolRep, 
                         forcefield: ForceField = ForceField()):
        rdmol = rdmol_rep.data
        component = OpenMMComponent.from_rdmol(rdmol, forcefield=forcefield)
        return cls(component=component)
    
    


