from __future__ import annotations

import os, pathlib
from typing import List, Tuple, Any
from copy import deepcopy

from rdkit import Chem
from rdkit.Chem import Conformer
from rdkit.Geometry import Point3D

from biopandas.pdb import PandasPdb

from openmm import (
    System, Context, VerletIntegrator, Vec3, unit, 
    LocalEnergyMinimizer, State, NonbondedForce
)
from openmm.app import ForceField, PDBFile, Modeller, Topology
from openmmforcefields.generators import SMIRNOFFTemplateGenerator
from pdbfixer import PDBFixer
from openff.toolkit import Molecule
from openff.units import unit as offunit
from openff.units.openmm import to_openmm as openff_unit_to_openmm
import parmed
from espaloma_charge.openff_wrapper import EspalomaChargeToolkitWrapper

from ..abstract import *
from ..utils.iotools import get_unique_full_name

class OpenMMComponent:
    def __init__(self, topology, system, positions, forcefield, 
                 name="component", work_dir="tmp", structure=None, original_intype=None, **kwargs):
        self.topology: Topology = topology
        self.system: System = system
        self.positions = positions
        self.forcefield: ForceField = forcefield
        self.structure: parmed.Structure = structure
        self.orig_intype: str = original_intype
        self.name = name
        self.work_dir = str(pathlib.Path(work_dir).absolute())
        self.others = kwargs
    
    @classmethod
    def create_component(cls, cinput: str | Chem.Mol, forcefield: ForceField, *args, **kwargs):
        if isinstance(cinput, str):
            if cinput.endswith(".pdb"):
                return cls.create_component_from_pdb(cinput, forcefield, *args, **kwargs)
            else:
                raise NotImplementedError()
        elif isinstance(cinput, Chem.Mol):
            return cls.create_component_from_rdmol(cinput, forcefield, *args, **kwargs)
        else:
            raise NotImplementedError()

    @classmethod
    def create_component_from_pdb(cls, pdb_path: str, forcefield: ForceField, 
                                  name: str = "component", *args, **kwargs):

        pdb = PDBFile(pdb_path)
        system = forcefield.createSystem(pdb.topology, *args, **kwargs)

        struct: parmed.Structure = parmed.openmm.load_topology(
            topology=pdb.topology, system=system, xyz=pdb.positions
        )
        return cls(pdb.topology, system, pdb.positions, forcefield, 
                   structure=struct, original_intype="pdb", name=name)

    @classmethod
    def create_component_from_rdmol(cls, rdmol: Chem.Mol, forcefield: ForceField, cache_path: str = None, 
                                    name: str = "component", *args, **kwargs):

        offmol = Molecule.from_rdkit(rdmol, hydrogens_are_explicit=True, allow_undefined_stereo=True)
        offmol.assign_partial_charges('espaloma-am1bcc', toolkit_registry=EspalomaChargeToolkitWrapper(), use_conformers=True)
        topology = offmol.to_topology().to_openmm()
        # positions = [Vec3(x, y, z) * unit.angstrom for x, y, z in offmol.conformers[0].m_as(offunit.angstrom)]
        positions = openff_unit_to_openmm(offmol.conformers[0])

        smirnoff = SMIRNOFFTemplateGenerator(molecules=offmol, cache=cache_path)
        forcefield.registerTemplateGenerator(smirnoff.generator)
        system = forcefield.createSystem(topology, *args, **kwargs)

        struct: parmed.Structure = parmed.openmm.load_topology(
            topology=topology, system=system, xyz=positions
        )

        return cls(topology, system, positions, forcefield, 
                   structure=struct, original_intype="rdmol", 
                   orig_rdmol=Chem.Mol(rdmol), name=name)
    
    @classmethod
    def merge_components(cls, components: List[OpenMMComponent], forcefield: ForceField, 
                         name: str = "mergedc", *args, **kwargs):
        
        modeller = Modeller(components[0].topology, components[0].positions)
        for c in components[1:]:
            modeller.add(c.topology, c.positions) # * unit.nanometer)
        system = forcefield.createSystem(modeller.topology, *args, **kwargs)
        
        struct: parmed.Structure = parmed.openmm.load_topology(
            topology=modeller.topology, system=system, xyz=modeller.positions
        )
        return cls(modeller.topology, system, modeller.positions, forcefield, 
                   structure=struct, original_intype="merged", name=name)
    
    def convert_to_orig_format(self):
        if self.orig_intype == "pdb":
            return self.to_pdb()
        elif self.orig_intype == "rdmol":
            return self.to_rdkit()
        else:
            raise NotImplementedError()

    def to_pdb(self, path=None):
        if path is None:
            path = get_unique_full_name(
                os.path.join(self.work_dir, f"{self.name}"), ".pdb", n=5
            )
        PDBFile.writeFile(self.topology, self.positions, path)
        return path

    def to_rdkit(self):
        new_mol = Chem.Mol(self.others["orig_rdmol"])
        conformer = new_mol.GetConformer()
        for i in range(new_mol.GetNumAtoms()):
            conformer.SetAtomPosition(i, Point3D(*self.positions[i].value_in_unit(unit.angstrom)))
        return new_mol
    
    def to_amber_prmtop(self, path=None):
        if path is None:
            path = get_unique_full_name(
                os.path.join(self.work_dir, f"{self.name}"), ".parm7", n=5
            )
        self.structure.save(path, overwrite=True)
        return path
        

class OpenMMGeneral:
    def __init__(self, forcefields: List[str] | ForceField, work_dir: str,
                 debug: bool = False):
        self.work_dir: str = work_dir
        if isinstance(forcefields, list):
            self.forcefield: ForceField = ForceField(*forcefields)
            self.forcefield_constructor = lambda: ForceField(*forcefields)
        else:
            self.forcefield: ForceField = forcefields
            self.forcefield_constructor = lambda: self.forcefield
        self.debug = debug

    def _optimize(self, components: List[OpenMMComponent], return_components: bool = False, *args, **kwargs
                 ) -> Tuple[OpenMMComponent | str, List[Any]]:
        
        merged = OpenMMComponent.merge_components(components, self.forcefield, *args, **kwargs)
        if self.debug: merged.to_pdb(os.path.join(self.work_dir, f"debug_preopt.pdb"))
        
        integrator = VerletIntegrator(0.002)
        context = Context(merged.system, integrator)
        context.setPositions(merged.positions)
        LocalEnergyMinimizer.minimize(context)

        minimized_state: State = context.getState(getPositions=True)
        minimized_positions = minimized_state.getPositions(asNumpy=True)
        merged.positions = minimized_positions
        if self.debug: merged.to_pdb(os.path.join(self.work_dir, f"debug_postopt.pdb"))

        start_index = 0
        for component in components:
            n_particles = component.topology.getNumAtoms()
            component.positions = minimized_positions[start_index:start_index+n_particles]
            start_index += n_particles
        
        if return_components:
            return merged, components
        
        out = merged.to_pdb(), [component.convert_to_orig_format() for component in components]
        return out # TODO: similar output to cuby4 optimize
    
    def _energy(self, components: List[OpenMMComponent], *args, **kwargs) -> float:

        merged = OpenMMComponent.merge_components(components, self.forcefield, *args, **kwargs)
        integrator = VerletIntegrator(0.002)
        context = Context(merged.system, integrator)
        context.setPositions(merged.positions)
        state: State = context.getState(getEnergy=True)
        penergy = state.getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole)
        return penergy

class OpenMMComplexOptimizer(OpenMMGeneral, ComplexConverter):
    name = "OpenMM Complex Optimization"
    def __init__(self, forcefields: List[str] | ForceField, work_dir: str = ".", jobname: str = "job"):
        super().__init__(forcefields, work_dir)
        self.jobname = jobname
    
    def convert(self, ligand: Chem.Mol, target: str) -> Tuple[Chem.Mol, str]:

        # atom_portion_path = os.path.join(self.work_dir, f"ommcopt_{self.jobname}_atm.pdb")
        # hetatm_portion_path = os.path.join
        # target_ppdb = PandasPdb().read_pdb(target)
        # target_ppdb.to_pdb()
        # TODO: Ability to have another small molecule included in the receptor/target
        self.forcefield = self.forcefield_constructor()
        components = [OpenMMComponent.create_component(ligand, self.forcefield, name=f"{self.jobname}_lig"), 
                      OpenMMComponent.create_component(target, self.forcefield, name=f"{self.jobname}_prot")]
        merged_pdb, (opt_ligand, opt_target) = self._optimize(components)
        return opt_ligand, opt_target

class OpenMMComplexInteractScorer(OpenMMGeneral, ComplexScorer):
    name = "OpenMM Interaction Energy"
    score_name = "openmm_interaction"
    def __init__(self, forcefields: List[str] | ForceField, work_dir: str = ".", 
                 optimize_before: bool = False, debug: bool = False):
        super().__init__(forcefields, work_dir, debug)
        self.optimize_before = optimize_before
    
    def score(self, ligand: Chem.Mol, target: str, return_used_components: bool = False) -> float:
        # TODO: Ability to have another small molecule included in the receptor/target

        self.forcefield = self.forcefield_constructor()
        components = [OpenMMComponent.create_component(ligand, self.forcefield), 
                      OpenMMComponent.create_component(target, self.forcefield)]
        if self.optimize_before:
            merged, components = self._optimize(components)
        
        if self.debug: print(self._energy(components), self._energy([components[0]]), self._energy([components[1]]))
        interaction_energy = self._energy(components) - self._energy([components[0]]) - self._energy([components[1]])
        
        if return_used_components:
            return interaction_energy, components
        else:
            return interaction_energy

class OpenMMUtils:
    @staticmethod
    def calc_protein_fcharge(pdb_path) -> int:

        pdb = PDBFile(pdb_path)
        ff = ForceField("amber14-all.xml")
        system: System = ff.createSystem(pdb.topology)
        nonbonded = [f for f in system.getForces() if isinstance(f, NonbondedForce)][0]

        total_charge = 0
        for i in range(system.getNumParticles()):
            charge, sigma, epsilon = nonbonded.getParticleParameters(i)
            total_charge += charge.value_in_unit(unit.elementary_charge)

        return round(total_charge)
    
class PDBFixerProteinPrepper(TargetConverter):
    name = "PDBFixer Protein Preparation"
    def __init__(self, ph: float = 7.4, chains: str = "all", keep_water: bool = False,
                 add_missing_residues: bool = False, out_dir: str = ".", save_prefix: str = "prepared"):
        self.chains = chains
        self.keep_water = keep_water
        self.ph = ph
        self.add_missing_residues = add_missing_residues

        if not os.path.exists(out_dir):
            os.mkdir(out_dir)
        self.out_dir: str = str(pathlib.Path(out_dir).absolute())
        self.save_prefix: str = save_prefix

    def convert(self, target_path: str) -> str:
        fixer: PDBFixer = PDBFixer(target_path)
        if self.chains != "all":
            rem_chains = [chain.id for chain in fixer.topology.chains() 
                        if chain.id not in self.chains]
            fixer.removeChains(chainIds=rem_chains)
        fixer.removeHeterogens(self.keep_water)
        if self.add_missing_residues:
            fixer.findMissingResidues()
            fixer.findMissingAtoms()
            fixer.addMissingAtoms()
        # fixer.findNonstandardResidues()
        # fixer.replaceNonstandardResidues()
        fixer.addMissingHydrogens(self.ph)  # TODO: Any need to pH value?

        save_path = os.path.join(self.out_dir, f"{self.save_prefix}_protein.pdb")
        PDBFile.writeFile(fixer.topology, fixer.positions, open(save_path, 'w')) 
        return save_path