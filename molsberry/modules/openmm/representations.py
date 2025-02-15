from __future__ import annotations

from typing import List, Any
from numpy import ndarray

from rdkit import Chem
from openmm import System, unit
from openmm.app import Topology, ForceField, Modeller
from openmmforcefields.generators import SMIRNOFFTemplateGenerator
from openff.toolkit import Molecule
from openff.units.openmm import to_openmm as openff_unit_to_openmm
from pdbfixer import PDBFixer
try:
    from espaloma_charge.openff_wrapper import EspalomaChargeToolkitWrapper
    ESPALOMA_CHARGE_AVAILABLE = True
except:
    ESPALOMA_CHARGE_AVAILABLE = False
    print("Warning: Espaloma charge failed to load in openff in openmm...")

from ...core import Molecule3DRep, SmallMolRep, ProteinRep, Representation
from ...core import SMILESRep, PDBPathRep, SDFPathRep
from ...modules.rdkit import RDKitMolRep

from .constants import DEFAULT_FORCEFIELDS, DEFAULT_PH

class OpenMMPreInputMolRep(Molecule3DRep):
    rep_name = "openmm_pre_mol"

    def __init__(
            self,
            topology: Topology,
            positions: Any,
            forcefield: ForceField | None = None) -> None:
        self.topology = topology
        self.positions = positions
        self.forcefield = forcefield

    @classmethod
    def from_RDKitMolRep(cls, rdkit_rep: RDKitMolRep) -> OpenMMPreInputMolRep:
        rdmol = rdkit_rep.content
        rdmol = Chem.AddHs(rdmol)
        molecule = Molecule.from_rdkit(rdmol, allow_undefined_stereo=True,
                                       hydrogens_are_explicit=False)
        if ESPALOMA_CHARGE_AVAILABLE:
            molecule.assign_partial_charges(
                'espaloma-am1bcc', 
                toolkit_registry=EspalomaChargeToolkitWrapper(), 
                use_conformers=True)
        else:
            molecule.assign_partial_charges('mmff94')
        
        topology = molecule.to_topology().to_openmm()
        positions = openff_unit_to_openmm(molecule.conformers[0])

        smirnoff = SMIRNOFFTemplateGenerator(molecules=molecule)
        forcefield = ForceField(*DEFAULT_FORCEFIELDS)
        forcefield.registerTemplateGenerator(smirnoff.generator)

        return cls(topology=topology, 
                   positions=positions, 
                   forcefield=forcefield)
    
    @classmethod
    def from_PDBPathRep(cls, pdb_rep: PDBPathRep) -> OpenMMPreInputMolRep:
        fixer = PDBFixer(pdb_rep.content)
        # fixer.removeHeterogens(keepWater=False) # TODO: Make it customizable.
        fixer.findMissingResidues()
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        fixer.addMissingHydrogens(pH=DEFAULT_PH)
        return cls(topology=fixer.topology, 
                   positions=fixer.positions, 
                   forcefield=None)
    
    def to_OpenMMInputMolRep(self):
        return OpenMMInputMolRep(
            topology=self.topology,
            positions=self.positions,
            forcefield=self.forcefield,
            system=self.forcefield.createSystem(self.topology)
        ) 
    
    def update_coordinates(self, coords: ndarray):
        raise NotImplementedError()
    
    def save_rep(self, exless_filename: str):
        raise NotImplementedError()
    
    @classmethod
    def save_rep_batch(cls, reps: List[Representation], exless_filename: str):
        raise NotImplementedError()
        
class OpenMMInputMolRep(Molecule3DRep):
    rep_name = "openmm_mol"

    def __init__(
            self, 
            topology: Topology, 
            system: System, 
            positions: Any, 
            forcefield: ForceField) -> None:
        
        self.topology: Topology = topology
        self.system: System = system
        self.positions = positions # TODO: Typing...
        self.forcefield: ForceField = forcefield

    @classmethod
    def from_RDKitMolRep(cls, rdkit_rep: RDKitMolRep) -> OpenMMInputMolRep:
        pre_rep = OpenMMPreInputMolRep.from_RDKitMolRep(rdkit_rep)
        system = pre_rep.forcefield.createSystem(pre_rep.topology)

        return cls(topology=pre_rep.topology, 
                   system=system, 
                   positions=pre_rep.positions, 
                   forcefield=pre_rep.forcefield)
    
    @classmethod
    def from_SDFPathRep(cls, sdf_rep: SDFPathRep) -> OpenMMInputMolRep:
        rdkit_rep = RDKitMolRep.from_SDFPathRep(sdf_rep)
        return cls.from_RDKitMolRep(rdkit_rep)
    
    @classmethod
    def from_PDBPathRep(cls, pdb_rep: PDBPathRep) -> OpenMMInputMolRep:
        pre_rep = OpenMMPreInputMolRep.from_PDBPathRep(pdb_rep)
        if pre_rep.forcefield is None:
            pre_rep.forcefield = ForceField(*DEFAULT_FORCEFIELDS)
        system = pre_rep.forcefield.createSystem(pre_rep.topology)
        return cls(topology=pre_rep.topology, 
                   system=system, 
                   positions=pre_rep.positions, 
                   forcefield=pre_rep.forcefield)
    
    @classmethod
    def merge_reps(cls, reps: List[OpenMMInputMolRep], 
                   forcefield: ForceField | None = None, **kwargs):
        
        if forcefield is None:
            forcefield = ForceField(*DEFAULT_FORCEFIELDS)

        modeller = Modeller(reps[0].topology, reps[0].positions)
        for rep in reps[1:]:
            modeller.add(rep.topology, rep.positions) # * unit.nanometer)

        system = forcefield.createSystem(modeller.topology, **kwargs)

        return cls(topology=modeller.topology, 
                   system=system, 
                   positions=modeller.positions, 
                   forcefield=forcefield)
    
    def update_coordinates(self, coords: ndarray):
        # TODO: Check and test
        self.positions = coords * unit.angstrom
    
    def save_rep(self, exless_filename: str):
        raise NotImplementedError()
    
    @classmethod
    def save_rep_batch(cls, reps: List[Representation], exless_filename: str):
        raise NotImplementedError()
    
    
    
    
    