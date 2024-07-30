import pytest

from rdkit import Chem

from moddipic.modules.rdkit.analyzers import (
    RDKitMWCalculator
)
from moddipic.core.data.special_cls import Ligand
from moddipic.core.data.representations import SMILESRep
from moddipic.modules.rdkit.representations import RDKitMolRep
from moddipic.core.data.collections import Batched

@pytest.fixture
def mwcalcer():
    return RDKitMWCalculator(debug=True, save_output=False)

@pytest.fixture
def input_ligands_smitype():
    smis = ["CCC", "COC(=O)CCCN", "c1cccnc1CC(C)(C)Cl"]
    ligs = [Ligand(SMILESRep(smi)) for smi in smis]
    return Batched(ligs)

@pytest.fixture
def input_lig_single():
    return Ligand(SMILESRep("CCC"))

def test_rdkit_mwcalcer(mwcalcer, 
                        input_lig_single,
                        input_ligands_smitype):
    output = mwcalcer.execute(input_ligands_smitype)
    assert len(output["molwt"]) == len(input_ligands_smitype)
    assert output["molwt"].get_basic_data_type() == float
    
    output = mwcalcer.execute(input_lig_single)
    assert len(output["molwt"]) == 1
    assert output["molwt"].get_basic_data_type() == float