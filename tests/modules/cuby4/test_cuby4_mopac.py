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
    Cuby4MOPACEnergyCalculator,
    Cuby4MOPACEnergyOptimizer,
    Cuby4MOPACInterfaceConfig
)

BASE_DIR = "./test_out"

def save_result_for_later_inspection(name):
    if os.path.exists(os.path.join(BASE_DIR, "job_MOPAC")):
        old_dir = os.path.join(BASE_DIR, "job_MOPAC")
        new_dir = os.path.join(BASE_DIR, f"_job_MOPAC_{name}")
        print(new_dir)
        os.rename(old_dir, new_dir)

@pytest.fixture()
def mopac_energy_block():
    iconf = Cuby4MOPACInterfaceConfig(mozyme=True)
    block = Cuby4MOPACEnergyCalculator(interface_config=iconf)
    return block

@pytest.fixture()
def mopac_optimize_block():
    block = Cuby4MOPACEnergyOptimizer()
    return block

@pytest.fixture()
def energy_pipeline(mopac_energy_block) -> Pipeline:

    class EnergyPipeline(Pipeline):
        name = "mopac_energy"
        display_name = "MOPAC Energy Pipeline"
        def build(self):
            self.add_block(InputBlock(["molecules"]), "input")

            self.add_block(mopac_energy_block, "mopacen")
            self.add_connection("input", "molecules", "mopacen", "molecules")

            self.add_block(OutputBlock(["energy"]), "output")
            self.add_connection("mopacen", "energy", "output", "energy")
        
    return EnergyPipeline(base_dir=BASE_DIR)

# @pytest.fixture()
# def optimize_pipeline_sm(amber_optimize_block, openff_parameterizer_block) \
#     -> Pipeline:

#     class OptimizePipeline(Pipeline):
#         name = "amber_energy"
#         display_name = "Amber Energy Pipeline"
#         def build(self):
#             self.add_block(InputBlock(["molecules"]), "input")

#             self.add_block(openff_parameterizer_block, "offparm")
#             self.add_connection("input", "molecules", "offparm", "molecules")

#             self.add_block(amber_optimize_block, "ambopt")
#             self.add_connection("offparm", "molecules", "ambopt", "molecules")

#             self.add_block(OutputBlock(["energy"]), "output")
#             self.add_connection("ambopt", "energy", "output", "energy")
        
#     return OptimizePipeline(base_dir=BASE_DIR)

def test_cuby4_mopac_energy_block_sm(energy_pipeline: Pipeline, 
                                     sample_sm_rdrep: RDKitMolRep):
    try:
        out = energy_pipeline.execute(molecules=MoleculeData(sample_sm_rdrep))
    
    except ValueError:
        name = sample_sm_rdrep.content.GetProp("_Name")
        save_result_for_later_inspection(name)
        raise ValueError()

    print("Final SM Energy:", out["energy"].get_representation_content()[0])

def test_cuby4_amber_energy_block_poc(energy_pipeline: Pipeline, 
                                      sample_poc_pdbrep: PDBPathRep):
    try:
        out = energy_pipeline.execute(molecules=MoleculeData(sample_poc_pdbrep))
    
    except ValueError:
        name = os.path.basename(sample_poc_pdbrep.content)
        save_result_for_later_inspection(name)
        raise ValueError()
    
    print("Final Poc Energy:", out["energy"].get_representation_content()[0])

# def test_cuby4_amber_optimize_block_sm(optimize_pipeline_sm: Pipeline,
#                                        sample_sm_rdrep: RDKitMolRep):
#     out = optimize_pipeline_sm.execute(molecules=MoleculeData(sample_sm_rdrep))
#     print("Final SM Energy:", out["energy"].get_representation_content()[0])

# TODO: Write test_cuby4_amber_optimize_block_prot