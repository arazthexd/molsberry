import pytest
import os, pathlib, shutil

from rdkit import Chem

from molsberry.core import (
    Pipeline, InputBlock, OutputBlock, 
    MoleculeData, 
    PDBPathRep
)
from molsberry.modules.rdkit import RDKitMolRep
from molsberry.modules.cuby4 import (
    Cuby4AMBEREnergyCalculator,
    Cuby4AMBEREnergyOptimizer,
)
from molsberry.modules.parmed import (
    OpenFFSmallMoleculeParameterizer,
    OpenMMProteinParameterizer,
)

BASE_DIR = "./test_out"

@pytest.fixture()
def amber_energy_block():
    block = Cuby4AMBEREnergyCalculator()
    return block

@pytest.fixture()
def amber_optimize_block():
    block = Cuby4AMBEREnergyOptimizer()
    return block

@pytest.fixture()
def openff_parameterizer_block():
    block = OpenFFSmallMoleculeParameterizer()
    return block

@pytest.fixture()
def openmm_parameterizer_block():
    block = OpenMMProteinParameterizer()
    return block

@pytest.fixture()
def energy_pipeline_sm(amber_energy_block, openff_parameterizer_block) \
    -> Pipeline:

    class EnergyPipeline(Pipeline):
        name = "amber_energy"
        display_name = "Amber Energy Pipeline"
        def build(self):
            self.add_block(InputBlock(["molecules"]), "input")

            self.add_block(openff_parameterizer_block, "offparm")
            self.add_connection("input", "molecules", "offparm", "molecules")

            self.add_block(amber_energy_block, "amben")
            self.add_connection("offparm", "molecules", "amben", "molecules")

            self.add_block(OutputBlock(["energy"]), "output")
            self.add_connection("amben", "energy", "output", "energy")
        
    return EnergyPipeline(base_dir=BASE_DIR)

@pytest.fixture()
def energy_pipeline_prot(amber_energy_block, openmm_parameterizer_block) \
    -> Pipeline:

    class EnergyPipeline(Pipeline):
        name = "amber_energy"
        display_name = "Amber Energy Pipeline"
        def build(self):
            self.add_block(InputBlock(["proteins"]), "input")

            self.add_block(openmm_parameterizer_block, "ommparm")
            self.add_connection("input", "proteins", "ommparm", "proteins")

            self.add_block(amber_energy_block, "amben")
            self.add_connection("ommparm", "proteins", "amben", "molecules")

            self.add_block(OutputBlock(["energy"]), "output")
            self.add_connection("amben", "energy", "output", "energy")
        
    return EnergyPipeline(base_dir=BASE_DIR)

@pytest.fixture()
def optimize_pipeline_sm(amber_optimize_block, openff_parameterizer_block) \
    -> Pipeline:

    class OptimizePipeline(Pipeline):
        name = "amber_energy"
        display_name = "Amber Energy Pipeline"
        def build(self):
            self.add_block(InputBlock(["molecules"]), "input")

            self.add_block(openff_parameterizer_block, "offparm")
            self.add_connection("input", "molecules", "offparm", "molecules")

            self.add_block(amber_optimize_block, "ambopt")
            self.add_connection("offparm", "molecules", "ambopt", "molecules")

            self.add_block(OutputBlock(["energy"]), "output")
            self.add_connection("ambopt", "energy", "output", "energy")
        
    return OptimizePipeline(base_dir=BASE_DIR)

def test_cuby4_amber_energy_block_sm(energy_pipeline_sm: Pipeline, 
                                     sample_sm_rdrep: RDKitMolRep):
    out = energy_pipeline_sm.execute(molecules=MoleculeData(sample_sm_rdrep))
    print("Final SM Energy:", out["energy"].get_representation_content()[0])

def test_cuby4_amber_energy_block_prot(energy_pipeline_prot: Pipeline, 
                                       sample_prot_pdbrep: PDBPathRep):
    out = energy_pipeline_prot.execute(proteins=MoleculeData(sample_prot_pdbrep))
    print("Final Prot Energy:", out["energy"].get_representation_content()[0])

def test_cuby4_amber_optimize_block_sm(optimize_pipeline_sm: Pipeline,
                                       sample_sm_rdrep: RDKitMolRep):
    out = optimize_pipeline_sm.execute(molecules=MoleculeData(sample_sm_rdrep))
    print("Final SM Energy:", out["energy"].get_representation_content()[0])

# TODO: Write test_cuby4_amber_optimize_block_prot
