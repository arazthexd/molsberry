import pytest
import os, pathlib

from rdkit import Chem

from molsberry.modules.rdkit.calculators import (
    RDKitMWCalculator,
    RDKitPLInteractionCalculator
)
from molsberry.core import Data, ProteinData, LigandData, BatchedData
from molsberry.core import SMILESRep, PDBPathRep
from molsberry.modules.rdkit import RDKitMolRep

@pytest.fixture
def mwcalcer():
    return RDKitMWCalculator(debug=True, save_output=False)

@pytest.fixture
def interactioncalcer():
    return RDKitPLInteractionCalculator(debug=True, save_output=False)

@pytest.fixture
def input_ligands_smitype():
    smis = ["CCC", "COC(=O)CCCN", "c1cccnc1CC(C)(C)Cl"]
    ligs = [LigandData(SMILESRep(smi)) for smi in smis]
    return BatchedData(ligs)

@pytest.fixture
def input_lig_single():
    return LigandData(SMILESRep("CCC"))

@pytest.fixture
def input_lig_from_file():
    sdf_path = os.path.join(str(pathlib.Path(__file__).parent.absolute()), 
                            "files", "lig.sdf")
    rdmol = next(Chem.SDMolSupplier(sdf_path, removeHs=False))
    return LigandData(RDKitMolRep(rdmol))

@pytest.fixture
def input_prot_from_file():
    pdb_path = os.path.join(str(pathlib.Path(__file__).parent.absolute()), 
                            "files", "rec.pdb")
    return ProteinData(PDBPathRep(pdb_path))

def test_rdkit_mwcalcer(mwcalcer, 
                        input_lig_single,
                        input_ligands_smitype):
    output = mwcalcer.execute(input_ligands_smitype)
    assert len(output["molwt"]) == len(input_ligands_smitype)
    assert issubclass(output["molwt"].basic_itype, Data)
    
    output = mwcalcer.execute(input_lig_single)
    assert len(output["molwt"]) == 1
    assert issubclass(output["molwt"].basic_itype, Data)

def test_rdkit_interactioncalcer(interactioncalcer, 
                                 input_lig_from_file,
                                 input_prot_from_file):
    output = interactioncalcer.execute(input_lig_from_file, 
                                       input_prot_from_file)
    print("output['interaction'] =", 
          output["e_interaction"].get_representation_content()[0])
    assert len(output["e_interaction"]) == 1
    assert issubclass(output["e_interaction"].basic_itype, Data)
    assert isinstance(output["e_interaction"].get_representation_content()[0], 
                      float)
    