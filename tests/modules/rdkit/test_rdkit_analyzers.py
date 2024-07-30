import pytest
import os, pathlib

from rdkit import Chem

from moddipic.modules.rdkit.analyzers import (
    RDKitMWCalculator,
    RDKitPLInteractionCalculator
)
from moddipic.core.data.special_cls import Ligand, Protein
from moddipic.core.data.representations import SMILESRep
from moddipic.modules.rdkit.representations import (
    RDKitMolRep, PDBPathProteinRep
)
from moddipic.core.data.collections import Batched

@pytest.fixture
def mwcalcer():
    return RDKitMWCalculator(debug=True, save_output=False)

@pytest.fixture
def interactioncalcer():
    return RDKitPLInteractionCalculator(debug=True, save_output=False)

@pytest.fixture
def input_ligands_smitype():
    smis = ["CCC", "COC(=O)CCCN", "c1cccnc1CC(C)(C)Cl"]
    ligs = [Ligand(SMILESRep(smi)) for smi in smis]
    return Batched(ligs)

@pytest.fixture
def input_lig_single():
    return Ligand(SMILESRep("CCC"))

@pytest.fixture
def input_lig_from_file():
    sdf_path = os.path.join(str(pathlib.Path(__file__).parent.absolute()), 
                            "files", "lig.sdf")
    rdmol = next(Chem.SDMolSupplier(sdf_path, removeHs=False))
    return Ligand(RDKitMolRep(rdmol))

@pytest.fixture
def input_prot_from_file():
    pdb_path = os.path.join(str(pathlib.Path(__file__).parent.absolute()), 
                            "files", "rec.pdb")
    return Protein(PDBPathProteinRep(pdb_path))

def test_rdkit_mwcalcer(mwcalcer, 
                        input_lig_single,
                        input_ligands_smitype):
    output = mwcalcer.execute(input_ligands_smitype)
    assert len(output["molwt"]) == len(input_ligands_smitype)
    assert output["molwt"].get_basic_data_type() == float
    
    output = mwcalcer.execute(input_lig_single)
    assert len(output["molwt"]) == 1
    assert output["molwt"].get_basic_data_type() == float

def test_rdkit_interactioncalcer(interactioncalcer, 
                                 input_lig_from_file,
                                 input_prot_from_file):
    interactioncalcer.input_context["protein"] = input_prot_from_file
    output = interactioncalcer.execute(input_lig_from_file)
    print("output['interaction'] =", output["interaction"][0])
    assert len(output["interaction"]) == 1
    assert output["interaction"].get_basic_data_type() == float
    assert set(output.keys()) == {"interaction"}
    assert set(interactioncalcer.required_input_keys) == {"ligands", "protein"}
    