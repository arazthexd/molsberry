import pytest

from rdkit import Chem

from moddipic.modules.rdkit.enumerators import (
    RDKitRingEnumerator,
    RDKitStereoEnumerator,
    RDKitTautEnumerator
)
from moddipic.core.data.special_cls import Ligand
from moddipic.core.data.representations import SMILESRep
from moddipic.modules.rdkit.representations import RDKitMolRep
from moddipic.core.data.collections import Batched
from moddipic.utils.moltools import is_mol_3d

@pytest.fixture
def ringenummer():
    return RDKitRingEnumerator(debug=True, save_output=False, flatten=True,
                               num_confs=100)

@pytest.fixture
def stereoenummer():
    return RDKitStereoEnumerator(debug=True, save_output=False, flatten=True)

@pytest.fixture
def tautenummer():
    return RDKitTautEnumerator(debug=True, save_output=False, flatten=True)

@pytest.fixture
def input_ligands_moltype():
    smis = ["C1CC(=O)CCC1F", "COC(=O)CCCN", 
            "CCc1cccnc1CC(C)(CC(I)(I)CC(I)CCCCC)I", "C1CCCCC1"]
    mols = [Chem.MolFromSmiles(smi) for smi in smis]
    ligs = [Ligand(RDKitMolRep(mol)) for mol in mols]
    return Batched(ligs)

@pytest.fixture
def input_ligands_smitype():
    smis = ["C1CC(=O)CCC1F", "COC(=O)CCCN", 
            "CCc1cccnc1CC(C)(CC(I)(I)CC(I)CCCCC)I", "C1CCCCC1"]
    ligs = [Ligand(SMILESRep(smi)) for smi in smis]
    return Batched(ligs)

def test_rdkit_tautenummer(tautenummer, 
                           input_ligands_moltype,
                           input_ligands_smitype):
    output = tautenummer.execute(input_ligands_moltype)
    assert len(output["ligands"]) > len(input_ligands_moltype)
    
    output = tautenummer.execute(input_ligands_smitype)
    assert len(output["ligands"]) > len(input_ligands_smitype)
    
def test_rdkit_stereoenummer(stereoenummer, 
                             input_ligands_moltype,
                             input_ligands_smitype):
    output = stereoenummer.execute(input_ligands_moltype)
    assert len(output["ligands"]) > len(input_ligands_moltype)
    
    output = stereoenummer.execute(input_ligands_smitype)
    assert len(output["ligands"]) > len(input_ligands_smitype)

def test_rdkit_ringenummer(ringenummer, 
                           input_ligands_moltype,
                           input_ligands_smitype):
    output = ringenummer.execute(input_ligands_moltype)
    assert len(output["ligands"]) > len(input_ligands_moltype)
    
    output = ringenummer.execute(input_ligands_smitype)
    assert len(output["ligands"]) > len(input_ligands_smitype)