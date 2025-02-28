# TODO: Needs cleaning up and generalising to any combination of QM,MM

import pytest
import os, pathlib, shutil

from rdkit import Chem

from molsberry.core import (
    Pipeline, InputBlock, OutputBlock, 
    MoleculeData, 
    PDBPathRep
)
from molsberry.modules.rdkit import RDKitMolRep
from molsberry.modules.rdkit import (
    RDKitLigandPocketLocator,
    RDKitPocketIsolator
)
from molsberry.modules.cuby4 import (
    Cuby4QMMMInterfaceConfig,
    Cuby4QMMMEnergyCalculator,
    Cuby4QMMMEnergyOptimizer,
    Cuby4MOPACInterfaceConfig,
    Cuby4AMBERInterfaceConfig
)
from molsberry.modules.parmed import (
    ParmedProteinPocketIsolator,
    ParmedMoleculeCombiner,
    OpenFFSmallMoleculeParameterizer,
    OpenMMProteinParameterizer
)

BASE_DIR = "./test_out"

MOPAMBER_INTERFACE_CONFIGS = [
    Cuby4QMMMInterfaceConfig(
        qm_config=Cuby4MOPACInterfaceConfig(method="pm6", 
                                            mozyme=True, 
                                            setcharges=True),
        mm_config=Cuby4AMBERInterfaceConfig(),
        n_threads=1
    )
]
MOPAMBER_INTERFACE_NAMES = [
    "mopac-pm6-mozyme-setcharge+amber"
]

##########################################################
##                       Fixtures                       ##
##########################################################

@pytest.fixture(params=MOPAMBER_INTERFACE_CONFIGS, ids=MOPAMBER_INTERFACE_NAMES)
def qmmm_mopamber_energy_block(request):
    iconf = request.param
    block = Cuby4QMMMEnergyCalculator(interface_config=iconf)
    return block

@pytest.fixture(params=MOPAMBER_INTERFACE_CONFIGS, ids=MOPAMBER_INTERFACE_NAMES)
def qmmm_mopamber_optimize_block(request):
    iconf = request.param
    block = Cuby4QMMMEnergyOptimizer(interface_config=iconf)
    return block

@pytest.fixture()
def pipeline_ligpoc_mopamber_energy(qmmm_mopamber_energy_block) -> Pipeline:

    class EnergyPipeline(Pipeline):
        name = "mopamber_energy"
        display_name = "Pocket Isolation/MOPAC-AMBER(QMMM) Energy Pipeline"
        def build(self):
            self.add_block(InputBlock(["proteins", "ligands"]), "input")

            self.add_block(RDKitLigandPocketLocator(10), "pocloc")
            self.add_connection("input", "ligands", "pocloc", "ligand")

            self.add_block(ParmedProteinPocketIsolator(), "pociso")
            self.add_connection("input", "proteins", "pociso", "protein")
            self.add_connection("pocloc", "location", "pociso", "location")

            self.add_block(OpenMMProteinParameterizer(), "protparam")
            self.add_connection("pociso", "pocket", "protparam", "proteins")
            
            self.add_block(OpenFFSmallMoleculeParameterizer(), "ligparam")
            self.add_connection("input", "ligands", "ligparam", "molecules")

            self.add_block(qmmm_mopamber_energy_block, "mopamberen")
            self.add_connection("ligparam", "molecules", "mopamberen", "qm_region")
            self.add_connection("protparam", "proteins", "mopamberen", "nonqm_region")

            self.add_block(OutputBlock(["energy"]), "output")
            self.add_connection("mopamberen", "energy", "output", "energy")
        
    return EnergyPipeline(base_dir=BASE_DIR)

@pytest.fixture()
def pipeline_ligpoc_mopamber_optimize(qmmm_mopamber_optimize_block) -> Pipeline:

    class OptimizePipeline(Pipeline):
        name = "mopamber_optimize"
        display_name = "Pocket Isolation/MOPAC-AMBER(QMMM) Optimize Pipeline"
        def build(self):
            self.add_block(InputBlock(["proteins", "ligands"]), "input")

            self.add_block(RDKitLigandPocketLocator(5), "pocloc")
            self.add_connection("input", "ligands", "pocloc", "ligand")

            self.add_block(ParmedProteinPocketIsolator(), "pociso")
            self.add_connection("input", "proteins", "pociso", "protein")
            self.add_connection("pocloc", "location", "pociso", "location")

            self.add_block(OpenMMProteinParameterizer(), "protparam")
            self.add_connection("pociso", "pocket", "protparam", "proteins")
            
            self.add_block(OpenFFSmallMoleculeParameterizer(), "ligparam")
            self.add_connection("input", "ligands", "ligparam", "molecules")

            self.add_block(qmmm_mopamber_optimize_block, "mopamberopt")
            self.add_connection("ligparam", "molecules", "mopamberopt", "qm_region")
            self.add_connection("protparam", "proteins", "mopamberopt", "nonqm_region")

            self.add_block(OutputBlock(["energy"]), "output")
            self.add_connection("mopamberopt", "energy", "output", "energy")
        
    return OptimizePipeline(base_dir=BASE_DIR)

##########################################################
##                         Tests                        ##
##########################################################

def test_cuby4_mopamber_comp_energy_from_plrex(
    pipeline_ligpoc_mopamber_energy: Pipeline, 
    sample_complex_rdpdbrep_from_plrex):

    lig, prot = sample_complex_rdpdbrep_from_plrex

    out = pipeline_ligpoc_mopamber_energy.execute(
        ligands=MoleculeData(lig),
        proteins=MoleculeData(prot)
    )

    print("Final Poc Energy:", out["energy"].get_representation_content()[0])

def test_cuby4_mopamber_comp_optimize_from_plrex(
    pipeline_ligpoc_mopamber_optimize: Pipeline, 
    sample_complex_rdpdbrep_from_plrex):

    lig, prot = sample_complex_rdpdbrep_from_plrex

    out = pipeline_ligpoc_mopamber_optimize.execute(
        ligands=MoleculeData(lig),
        proteins=MoleculeData(prot)
    )

    print("Final Poc Energy:", out["energy"].get_representation_content()[0])