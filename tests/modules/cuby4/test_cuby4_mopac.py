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
    Cuby4MOPACEnergyCalculator,
    Cuby4MOPACEnergyOptimizer,
    Cuby4MOPACInterfaceConfig
)
from molsberry.modules.parmed import (
    ParmedProteinPocketIsolator
)

BASE_DIR = "./test_out"

MOPAC_INTERFACE_CONFIGS = [
    Cuby4MOPACInterfaceConfig(method="pm6"),
    # Cuby4MOPACInterfaceConfig(method="pm6", mozyme=True),
    Cuby4MOPACInterfaceConfig(method="pm7", mozyme=True, setcharges=True)
] # TODO: setpi, cvb tests
MOPAC_INTERFACE_NAMES = [
    "default-pm6",
    # "default-pm6-mozyme",
    "default-pm7-mozyme-setcharge"
]

##########################################################
##                         Utils                        ##
##########################################################

def save_result_for_later_inspection(name):
    old_dir = os.path.join(BASE_DIR, "job_MOPAC")
    new_dir = os.path.join(BASE_DIR, f"_job_MOPAC_{name}")

    if os.path.exists(os.path.join(BASE_DIR, "job_MOPAC")):
        if os.path.exists(new_dir):
            shutil.rmtree(new_dir)
        os.rename(old_dir, new_dir)

##########################################################
##                       Fixtures                       ##
##########################################################

@pytest.fixture(params=MOPAC_INTERFACE_CONFIGS, ids=MOPAC_INTERFACE_NAMES)
def mopac_energy_block(request):
    iconf = request.param
    block = Cuby4MOPACEnergyCalculator(interface_config=iconf)
    return block

@pytest.fixture(params=MOPAC_INTERFACE_CONFIGS, ids=MOPAC_INTERFACE_NAMES)
def mopac_optimize_block(request):
    iconf = request.param
    block = Cuby4MOPACEnergyOptimizer(interface_config=iconf)
    return block

@pytest.fixture()
def pipeline_sm_energy(mopac_energy_block) -> Pipeline:

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

@pytest.fixture()
def pipeline_sm_optimize(mopac_optimize_block) -> Pipeline:

    class OptimizePipeline(Pipeline):
        name = "mopac_optimize"
        display_name = "MOPAC Optimize Pipeline"
        def build(self):
            self.add_block(InputBlock(["molecules"]), "input")

            self.add_block(mopac_optimize_block, "mopacopt")
            self.add_connection("input", "molecules", "mopacopt", "molecules")

            self.add_block(OutputBlock(["energy"]), "output")
            self.add_connection("mopacopt", "energy", "output", "energy")
        
    return OptimizePipeline(base_dir=BASE_DIR)

@pytest.fixture(ids=["default-pm6-mozyme-setcharge"])
def pipeline_prot_isol_energy() -> Pipeline:

    class EnergyPipeline(Pipeline):
        name = "mopac_energy"
        display_name = "Pocket Isolation/MOPAC Energy Pipeline"
        def build(self):
            self.add_block(InputBlock(["proteins", "ligands"]), "input")

            self.add_block(RDKitLigandPocketLocator(10), "pocloc")
            self.add_connection("input", "ligands", "pocloc", "ligand")

            self.add_block(ParmedProteinPocketIsolator(), "pociso")
            self.add_connection("input", "proteins", "pociso", "protein")
            self.add_connection("pocloc", "location", "pociso", "location")

            iconf = Cuby4MOPACInterfaceConfig(method="pm6", 
                                              mozyme=True, 
                                              setcharges=True)
            self.add_block(Cuby4MOPACEnergyCalculator(iconf), "mopacen")
            self.add_connection("pociso", "pocket", "mopacen", "molecules")

            self.add_block(OutputBlock(["energy"]), "output")
            self.add_connection("mopacen", "energy", "output", "energy")
        
    return EnergyPipeline(base_dir=BASE_DIR)


##########################################################
##                         Tests                        ##
##########################################################

def test_cuby4_mopac_energy_sm_from_plrex(
    pipeline_sm_energy: Pipeline, 
    sample_sm_rdrep_from_plrex: RDKitMolRep):

    try:
        out = pipeline_sm_energy.execute(
            molecules=MoleculeData(sample_sm_rdrep_from_plrex)
        )
    
    except ValueError:
        name = sample_sm_rdrep_from_plrex.content.GetProp("_Name")
        save_result_for_later_inspection(name)
        raise ValueError()

    print("Final SM Energy:", out["energy"].get_representation_content()[0])

def test_cuby4_mopac_energy_sm_from_smi(
    pipeline_sm_energy: Pipeline, 
    sample_sm_rdrep_from_smi: RDKitMolRep):

    try:
        out = pipeline_sm_energy.execute(
            molecules=MoleculeData(sample_sm_rdrep_from_smi)
        )
    
    except ValueError:
        name = sample_sm_rdrep_from_smi.content.GetProp("_Name")
        save_result_for_later_inspection(name)
        raise ValueError()

    print("Final SM Energy:", out["energy"].get_representation_content()[0])

def test_cuby4_mopac_optimize_sm_from_plrex(
    pipeline_sm_optimize: Pipeline, 
    sample_sm_rdrep_from_plrex: RDKitMolRep):

    try:
        out = pipeline_sm_optimize.execute(
            molecules=MoleculeData(sample_sm_rdrep_from_plrex)
        )
    
    except ValueError:
        name = sample_sm_rdrep_from_plrex.content.GetProp("_Name")
        save_result_for_later_inspection(name)
        raise ValueError()

    print("Final SM Energy:", out["energy"].get_representation_content()[0])

def test_cuby4_mopac_pociso_energy_from_plrex(
    pipeline_prot_isol_energy: Pipeline, 
    sample_complex_rdpdbrep_from_plrex):

    lig, prot = sample_complex_rdpdbrep_from_plrex

    try:
        out = pipeline_prot_isol_energy.execute(
            ligands=MoleculeData(lig),
            proteins=MoleculeData(prot)
        )
    
    except ValueError:
        name = "last_poc_energy"
        save_result_for_later_inspection(name)
        raise ValueError()

    print("Final Poc Energy:", out["energy"].get_representation_content()[0])


# def test_cuby4_mopac_energy_block_poc(pipeline_sm_energy: Pipeline, 
#                                       sample_poc_pdbrep: PDBPathRep):
#     try:
#         out = pipeline_sm_energy.execute(molecules=MoleculeData(sample_poc_pdbrep))
    
#     except ValueError:
#         name = os.path.basename(sample_poc_pdbrep.content)
#         save_result_for_later_inspection(name)
#         raise ValueError()
    
#     print("Final Poc Energy:", out["energy"].get_representation_content()[0])

# def test_cuby4_amber_optimize_block_sm(optimize_pipeline_sm: Pipeline,
#                                        sample_sm_rdrep: RDKitMolRep):
#     out = optimize_pipeline_sm.execute(molecules=MoleculeData(sample_sm_rdrep))
#     print("Final SM Energy:", out["energy"].get_representation_content()[0])

# TODO: Write test_cuby4_amber_optimize_block_prot