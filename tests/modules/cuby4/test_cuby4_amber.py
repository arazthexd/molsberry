import pytest
import os, pathlib, shutil

from rdkit import Chem

from molsberry.core import (
    Pipeline, InputBlock, OutputBlock, 
    MoleculeData, 
    PDBPathRep
)
from molsberry.modules.rdkit import (
    RDKitMolRep,
    RDKitLigandPocketLocator,
    RDKitPocketIsolator
)
from molsberry.modules.cuby4 import (
    Cuby4AMBEREnergyCalculator,
    Cuby4AMBEREnergyOptimizer,
    Cuby4AMBERInterfaceConfig
)
from molsberry.modules.parmed import (
    OpenFFSmallMoleculeParameterizer,
    OpenMMProteinParameterizer,
    ParmedProteinPocketIsolator
)

BASE_DIR = "./test_out"

##########################################################
##                       Fixtures                       ##
##########################################################

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
def pipeline_sm_param_energy(amber_energy_block, openff_parameterizer_block) \
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
def pipeline_prot_param_energy(amber_energy_block, openmm_parameterizer_block) \
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
def pipeline_prot_isol_param_energy(amber_energy_block, 
                                    openmm_parameterizer_block) -> Pipeline:

    class EnergyPipeline(Pipeline):
        name = "isloate_amber_energy"
        display_name = "Isolate/Amber Energy Pipeline"
        def build(self):
            self.add_block(InputBlock(["proteins", "ligands"]), "input")

            self.add_block(RDKitLigandPocketLocator(10), "pocloc")
            self.add_connection("input", "ligands", "pocloc", "ligand")

            self.add_block(ParmedProteinPocketIsolator(), "pociso")
            self.add_connection("input", "proteins", "pociso", "protein")
            self.add_connection("pocloc", "location", "pociso", "location")

            self.add_block(openmm_parameterizer_block, "ommparm")
            self.add_connection("pociso", "pocket", "ommparm", "proteins")

            self.add_block(amber_energy_block, "amben")
            self.add_connection("ommparm", "proteins", "amben", "molecules")

            self.add_block(OutputBlock(["energy"]), "output")
            self.add_connection("amben", "energy", "output", "energy")
        
    return EnergyPipeline(base_dir=BASE_DIR)

@pytest.fixture()
def pipeline_sm_param_optimize(amber_optimize_block, 
                               openff_parameterizer_block) -> Pipeline:

    class OptimizePipeline(Pipeline):
        name = "amber_optimize"
        display_name = "Amber Optimize Pipeline"
        def build(self):
            self.add_block(InputBlock(["molecules"]), "input")

            self.add_block(openff_parameterizer_block, "offparm")
            self.add_connection("input", "molecules", "offparm", "molecules")

            self.add_block(amber_optimize_block, "ambopt")
            self.add_connection("offparm", "molecules", "ambopt", "molecules")

            self.add_block(OutputBlock(["energy"]), "output")
            self.add_connection("ambopt", "energy", "output", "energy")
        
    return OptimizePipeline(base_dir=BASE_DIR)

@pytest.fixture()
def pipeline_prot_isol_param_optimize(amber_optimize_block, 
                                      openmm_parameterizer_block) -> Pipeline:

    class OptimizePipeline(Pipeline):
        name = "amber_optimize"
        display_name = "Pocket Isolation/Amber Optimize Pipeline"
        def build(self):
            self.add_block(InputBlock(["proteins", "ligands"]), "input")

            self.add_block(RDKitLigandPocketLocator(10), "pocloc")
            self.add_connection("input", "ligands", "pocloc", "ligand")

            self.add_block(ParmedProteinPocketIsolator(), "pociso")
            self.add_connection("input", "proteins", "pociso", "protein")
            self.add_connection("pocloc", "location", "pociso", "location")

            self.add_block(openmm_parameterizer_block, "ommparm")
            self.add_connection("pociso", "pocket", "ommparm", "proteins")

            self.add_block(amber_optimize_block, "ambopt")
            self.add_connection("ommparm", "proteins", "ambopt", "molecules")

            self.add_block(OutputBlock(["energy"]), "output")
            self.add_connection("ambopt", "energy", "output", "energy")
        
    return OptimizePipeline(base_dir=BASE_DIR)


##########################################################
##                         Tests                        ##
##########################################################

### Small Molecule / Energy
def test_cuby4_amber_energy_sm_from_smi(pipeline_sm_param_energy: Pipeline, 
                                        sample_sm_rdrep_from_smi: RDKitMolRep):
    out = pipeline_sm_param_energy.execute(
        molecules=MoleculeData(sample_sm_rdrep_from_smi)
    )
    print("Final SM Energy:", out["energy"].get_representation_content()[0])

def test_cuby4_amber_energy_sm_from_plrex(pipeline_sm_param_energy: Pipeline, 
                                          sample_sm_rdrep_from_plrex: RDKitMolRep):
    out = pipeline_sm_param_energy.execute(
        molecules=MoleculeData(sample_sm_rdrep_from_plrex)
    )
    print("Final SM Energy:", out["energy"].get_representation_content()[0])

def test_cuby4_amber_energy_sm_from_zinc(pipeline_sm_param_energy: Pipeline, 
                                         sample_sm_rdrep_from_zinc: RDKitMolRep):
    out = pipeline_sm_param_energy.execute(
        molecules=MoleculeData(sample_sm_rdrep_from_zinc)
    )
    print("Final SM Energy:", out["energy"].get_representation_content()[0])

### Small Molecule / Optimize
def test_cuby4_amber_energy_sm_from_plrex(
        pipeline_sm_param_optimize: Pipeline,
        sample_sm_rdrep_from_plrex: RDKitMolRep
    ):

    out = pipeline_sm_param_optimize.execute(
        molecules=MoleculeData(sample_sm_rdrep_from_plrex)
    )
    print("Final SM Energy:", out["energy"].get_representation_content()[0])

### Protein / Energy
def test_cuby4_amber_energy_prot_from_plrex(
        pipeline_prot_param_energy: Pipeline, 
        sample_prot_pdbrep_from_plrex: PDBPathRep
    ):

    out = pipeline_prot_param_energy.execute(
        proteins=MoleculeData(sample_prot_pdbrep_from_plrex)
    )
    print("Final Prot Energy:", out["energy"].get_representation_content()[0])

### Pocket Isolation / Energy
def test_cuby4_amber_pociso_energy_from_plrex(
        pipeline_prot_isol_param_energy: Pipeline, 
        sample_complex_rdpdbrep_from_plrex
    ):

    lig, prot = sample_complex_rdpdbrep_from_plrex
    out = pipeline_prot_isol_param_energy.execute(
        ligands=MoleculeData(lig),
        proteins=MoleculeData(prot)
    )
    print("Final Poc Energy:", out["energy"].get_representation_content()[0])


### Pocket Isolation / Optimize
def test_cuby4_amber_pociso_optimize_from_plrex(
        pipeline_prot_isol_param_optimize: Pipeline, 
        sample_complex_rdpdbrep_from_plrex: PDBPathRep
    ):

    lig, prot = sample_complex_rdpdbrep_from_plrex
    out = pipeline_prot_isol_param_optimize.execute(
        ligands=MoleculeData(lig),
        proteins=MoleculeData(prot)
    )
    print("Final Poc Energy:", out["energy"].get_representation_content()[0])