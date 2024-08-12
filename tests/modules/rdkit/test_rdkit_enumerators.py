import pytest

from rdkit import Chem

from molsberry.modules.rdkit import (
    RDKitLigEnumeratorBlock,
    RDKitLigandTautEnumerator,
    RDKitLigandStereoEnumerator,
    RDKitLigandRingEnumerator
)
from molsberry.core.data.data_types import LigandData
from molsberry.core.data.representations import SMILESRep
from molsberry.modules.rdkit.representations import RDKitMolRep
from molsberry.core.data.collections import BatchedData
from molsberry.utils.moltools import is_mol_3d

@pytest.fixture
def ringenummer():
    return RDKitLigandRingEnumerator(debug=True, save_output=False, 
                                     flatten=True, num_confs=20, 
                                     dist_threshold=0.3)

@pytest.fixture
def stereoenummer():
    return RDKitLigandStereoEnumerator(debug=True, save_output=False, 
                                       flatten=True)

@pytest.fixture
def tautenummer():
    return RDKitLigandTautEnumerator(debug=True, save_output=False, 
                                     flatten=True)

@pytest.fixture
def input_ligands_moltype():
    smis = ["C1CC(=O)CCC1F", "COC(=O)CCCN", 
            "CCc1cccnc1CC(C)(CC(I)(I)CC(I)CCCCC)I", "C1CCCCC1"]
    mols = [Chem.MolFromSmiles(smi) for smi in smis]
    ligs = [LigandData(RDKitMolRep(mol)) for mol in mols]
    return BatchedData(ligs)

@pytest.fixture
def input_ligands_smitype():
    smis = ["C1CC(=O)CCC1F", "COC(=O)CCCN", 
            "CCc1cccnc1CC(C)(CC(I)(I)CC(I)CCCCC)I", "C1CCCCC1"]
    ligs = [LigandData(SMILESRep(smi)) for smi in smis]
    return BatchedData(ligs)

def test_rdkit_tautenummer(tautenummer, 
                           input_ligands_moltype,
                           input_ligands_smitype):
    output = tautenummer.execute(input_ligands_moltype)
    assert len(output["ligands"]) > len(input_ligands_moltype)
    assert output["ligands"].basic_itype == LigandData
    
    output = tautenummer.execute(input_ligands_smitype)
    assert len(output["ligands"]) > len(input_ligands_smitype)
    assert output["ligands"].basic_itype == LigandData
    
def test_rdkit_stereoenummer(stereoenummer, 
                             input_ligands_moltype,
                             input_ligands_smitype):
    output = stereoenummer.execute(input_ligands_moltype)
    assert len(output["ligands"]) > len(input_ligands_moltype)
    assert output["ligands"].basic_itype == LigandData
    
    output = stereoenummer.execute(input_ligands_smitype)
    assert len(output["ligands"]) > len(input_ligands_smitype)
    assert output["ligands"].basic_itype == LigandData

def test_rdkit_ringenummer(ringenummer, 
                           input_ligands_moltype,
                           input_ligands_smitype):
    output = ringenummer.execute(input_ligands_moltype)
    assert len(output["ligands"]) > len(input_ligands_moltype)
    assert output["ligands"].basic_itype == LigandData
    
    output = ringenummer.execute(input_ligands_smitype)
    assert len(output["ligands"]) > len(input_ligands_smitype)
    assert output["ligands"].basic_itype == LigandData