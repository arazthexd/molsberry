import pytest

from rdkit import Chem

from moddipic.modules.rdkit.converters import (
    RDKitLigandEmbedder,
    RDKitLigandHAdder
)
from moddipic.core.data.special_cls import Ligand
from moddipic.core.data.representations import SMILESRep
from moddipic.modules.rdkit.representations import RDKitMolRep
from moddipic.core.data.collections import Batched
from moddipic.utils.moltools import is_mol_3d

@pytest.fixture
def ligembedder():
    return RDKitLigandEmbedder(debug=True, save_output=False)

@pytest.fixture
def lighadder():
    return RDKitLigandHAdder(debug=True, save_output=False)

@pytest.fixture
def input_ligands_moltype():
    smis = ["CCC", "COC(=O)CCCN", "c1cccnc1CC(C)(C)Cl"]
    mols = [Chem.MolFromSmiles(smi) for smi in smis]
    ligs = [Ligand(RDKitMolRep(mol)) for mol in mols]
    return Batched(ligs)

@pytest.fixture
def input_ligands_smitype():
    smis = ["CCC", "COC(=O)CCCN", "c1cccnc1CC(C)(C)Cl"]
    ligs = [Ligand(SMILESRep(smi)) for smi in smis]
    return Batched(ligs)

def test_rdkit_ligembedder(ligembedder, 
                           input_ligands_moltype,
                           input_ligands_smitype):
    output = ligembedder.execute(input_ligands_moltype)
    assert len(output["ligands"]) == len(input_ligands_moltype)
    assert all(is_mol_3d(lig.get_data(RDKitMolRep)) 
               for lig in output["ligands"])
    assert all(not is_mol_3d(lig.get_data(RDKitMolRep)) 
               for lig in input_ligands_moltype)
    
    output = ligembedder.execute(input_ligands_smitype)
    assert len(output["ligands"]) == len(input_ligands_smitype)
    assert all(is_mol_3d(lig.get_data(RDKitMolRep)) 
               for lig in output["ligands"])
    
def test_rdkit_lighadder(lighadder, 
                         input_ligands_moltype,
                         input_ligands_smitype):
    output = lighadder.execute(input_ligands_moltype)
    assert len(output["ligands"]) == len(input_ligands_moltype)
    assert all("H" in Chem.MolToSmiles(lig.get_data(RDKitMolRep)) 
               for lig in output["ligands"])
    assert all("H" not in Chem.MolToSmiles(lig.get_data(RDKitMolRep)) 
               for lig in input_ligands_moltype)
    
    output = lighadder.execute(input_ligands_smitype)
    assert len(output["ligands"]) == len(input_ligands_smitype)
    assert all("H" in Chem.MolToSmiles(lig.get_data(RDKitMolRep)) 
               for lig in output["ligands"])