import os

from rdkit import Chem, Geometry
from rdkit.Chem import rdDistGeom, rdRascalMCES

from openff.toolkit import Molecule, Topology
from espaloma_charge.openff_wrapper import EspalomaChargeToolkitWrapper
from openmm.app import PDBFile, ForceField, Simulation, Modeller
from openmm import VerletIntegrator, Platform, Context, app, unit, Vec3, System, State
from openmmforcefields.generators import SMIRNOFFTemplateGenerator
import parmed

WORK_DIR = "tmp"
STARTING_FORCEFIELDS = ["amber14-all.xml", "implicit/gbn2.xml"]

def prepare_lig_for_openmm(mol: Chem.Mol, forcefield: ForceField, job_name: str = "job"):

    cache_path = os.path.join(WORK_DIR, "template_cache.json")
    lig_pdb_path = os.path.join(WORK_DIR, f"off_{job_name}_tmplig.pdb")
    # lig_sdf_path = os.path.join(WORK_DIR, f"off_{job_name}_tmplig.sdf")

    ligand = Molecule.from_rdkit(mol, hydrogens_are_explicit=True, allow_undefined_stereo=True)
    ligand.assign_partial_charges('espaloma-am1bcc', toolkit_registry=EspalomaChargeToolkitWrapper(), use_conformers=True)
    ligand.to_file(lig_pdb_path, "pdb")
    # ligand.to_file(lig_sdf_path, "sdf")
    ligand_top = ligand.to_topology().to_openmm()
    ligand_pdb = PDBFile(lig_pdb_path)
    # print(ligand.partial_charges)
    
    smirnoff = SMIRNOFFTemplateGenerator(molecules=ligand, cache=cache_path)
    forcefield.registerTemplateGenerator(smirnoff.generator)
    ligand_sys = forcefield.createSystem(ligand_top)

    ligand_struct: parmed.Structure = parmed.openmm.load_topology(
        topology=ligand_top, system=ligand_sys, xyz=ligand_pdb.positions
    )
    return {
        "parmed_structure": ligand_struct,
        "openff_molecule": ligand,
        "openmm_topology": ligand_top,
        "openmm_pdb": ligand_pdb,
        "openmm_system": ligand_sys
    }

def prepare_target_pdb_for_openmm(target_pdb: str, forcefield: ForceField):

    pocket_pdb = PDBFile(target_pdb)
    
    pocket_sys = forcefield.createSystem(pocket_pdb.topology)

    pocket_struct: parmed.Structure = parmed.openmm.load_topology(
        topology=pocket_pdb.topology, system=pocket_sys, xyz=pocket_pdb.positions
    )
    return {
        "openmm_pdb": pocket_pdb,
        "openmm_system": pocket_sys,
        "parmed_structure": pocket_struct
    }

def minimize_and_calculate_interaction(mol: Chem.Mol, target_pdb: str, job_name: str = "job", extra_outputs: list = []):
    
    # Define Forcefield
    forcefield = ForceField(*STARTING_FORCEFIELDS)

    # Parameterize Ligand
    ligand_omm_related = prepare_lig_for_openmm(mol, forcefield, job_name)
    ligand_struct = ligand_omm_related["parmed_structure"]
    ligand_sys = ligand_omm_related["openmm_system"]
    ligand_pdb: PDBFile = ligand_omm_related["openmm_pdb"]
    ligand_context = Context(ligand_sys, VerletIntegrator(0.002))

    # Parameterize Pocket
    pocket_omm_related = prepare_target_pdb_for_openmm(target_pdb, forcefield)
    # pocket_struct = pocket_omm_related["parmed_structure"]
    pocket_sys: System = pocket_omm_related["openmm_system"]
    pocket_pdb: PDBFile = pocket_omm_related["openmm_pdb"]
    pocket_context = Context(pocket_sys, VerletIntegrator(0.002))

    # Combine Lig and Poc and create the complex system
    complex_modeller = Modeller(pocket_pdb.topology, pocket_pdb.positions)
    complex_modeller.add(ligand_pdb.topology, ligand_pdb.positions)
    complex_sys = forcefield.createSystem(complex_modeller.topology)
    # complex_struct: parmed.Structure = pocket_struct + ligand_struct
    # complex_sys = complex_struct.createSystem(implicitSolvent=implicit)
    # print(complex_struct.get_coordinates())

    # Minimize Energy
    simulation = Simulation(complex_modeller.topology, complex_sys, VerletIntegrator(0.002), 
                            platform=Platform.getPlatformByName("CUDA"))
    simulation.context.setPositions(complex_modeller.positions)
    simulation.minimizeEnergy()
    complex_state: State = simulation.context.getState(getEnergy=True, getPositions=True)
    complex_positions = complex_state.getPositions()

    # Get Individual Protein and Ligand Energies
    pocket_context.setPositions(complex_positions[:pocket_sys.getNumParticles()])
    pocket_state: State = pocket_context.getState(getEnergy=True, getPositions=True)
    ligand_context.setPositions(complex_positions[pocket_sys.getNumParticles():])
    ligand_state: State = ligand_context.getState(getEnergy=True, getPositions=True)

    # PDBFile.writeFile(complex_struct.topology, complex_struct.positions, open("before.pdb", "w"))
    # PDBFile.writeFile(complex_struct.topology, complex_positions, open("after.pdb", "w"))
    # PDBFile.writeFile(ligand_struct.topology, ligand_struct.positions, open("ligbefore.pdb", "w"))
    # PDBFile.writeFile(ligand_struct.topology, ligand_state.getPositions(), open("ligafter.pdb", "w"))
    interaction_energy = complex_state.getPotentialEnergy() - \
        pocket_state.getPotentialEnergy() - ligand_state.getPotentialEnergy()
    
    if not extra_outputs:
        return interaction_energy.value_in_unit(unit.kilocalorie_per_mole) 
    
    output = dict()
    output["interaction_energy"] = interaction_energy.value_in_unit(unit.kilocalorie_per_mole) 
    if "ligand_pdb" in extra_outputs:
        final_lig_pdb_path = os.path.join(WORK_DIR, f"gbmin_{job_name}_ligand.pdb")
        PDBFile.writeFile(ligand_struct.topology, ligand_state.getPositions(), open(final_lig_pdb_path, "w"))
        output["ligand_pdb"] = final_lig_pdb_path
    if "ligand_rdkit" in extra_outputs:
        final_lig = Chem.Mol(mol, confId=-1)
        final_lig_conf = final_lig.GetConformer()
        for i, vec3 in enumerate(ligand_state.getPositions()):
            vec3 = vec3.value_in_unit(unit.angstrom)
            vec3: Vec3
            final_lig_conf.SetAtomPosition(i, Geometry.Point3D(vec3.x, vec3.y, vec3.z))
        output["ligand_rdkit"] = final_lig
    if "pocket_pdb" in extra_outputs:
        final_poc_pdb_path = os.path.join(WORK_DIR, f"gbmin_{job_name}_target.pdb")
        PDBFile.writeFile(pocket_pdb.topology, pocket_state.getPositions(), open(final_poc_pdb_path, "w"))
        output["pocket_pdb"] = final_poc_pdb_path
    if "complex_pdb" in extra_outputs:
        final_complex_pdb_path = os.path.join(WORK_DIR, f"gbmin_{job_name}_complex.pdb")
        PDBFile.writeFile(complex_modeller.topology, complex_state.getPositions(), open(final_complex_pdb_path, "w"))
        output["complex_pdb"] = final_poc_pdb_path
    # TODO: Other possible outputs needed of this function
    return output
    
    
    

# m = list(Chem.SDMolSupplier("/home/arazthexd/projects/002_sqm/tmp/docked_ligs.sdf", removeHs=False))[3]
# m = Chem.AddHs(m, addCoords=True)
# minimize_and_calculate_interaction(m, "output/prepared_protein.pdb", None)

