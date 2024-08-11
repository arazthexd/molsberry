import pytest

from rdkit import Chem

from molsberry.modules.rdkit.converters import (
    RDKitLigandEmbedder,
    RDKitLigandHAdder
)
from molsberry.core.data import LigandData, BatchedData
from molsberry.core.data import SMILESRep
from molsberry.modules.rdkit.representations import RDKitSmallMolRep
from molsberry.modules.rdkit.interface import RDKitInterface

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
    ligs = [LigandData(RDKitSmallMolRep(mol)) for mol in mols]
    return BatchedData(ligs)

@pytest.fixture
def input_ligands_smitype():
    smis = ["CCC", "COC(=O)CCCN", "c1cccnc1CC(C)(C)Cl"]
    ligs = [LigandData(SMILESRep(smi)) for smi in smis]
    return BatchedData(ligs)

def test_rdkit_ligembedder(ligembedder, 
                           input_ligands_moltype,
                           input_ligands_smitype):
    output = ligembedder.execute(input_ligands_moltype)
    assert len(output["ligands"]) == len(input_ligands_moltype)
    assert all(RDKitInterface.is_mol_3d(
        lig.get_representation_content(RDKitSmallMolRep)) 
        for lig in output["ligands"])
    assert output["ligands"].basic_itype == LigandData
    assert all(not RDKitInterface.is_mol_3d(
        lig.get_representation_content(RDKitSmallMolRep)) 
        for lig in input_ligands_moltype)
    
    output = ligembedder.execute(input_ligands_smitype)
    assert len(output["ligands"]) == len(input_ligands_smitype)
    assert all(RDKitInterface.is_mol_3d(
        lig.get_representation_content(RDKitSmallMolRep)) 
        for lig in output["ligands"])
    assert output["ligands"].basic_itype == LigandData
    
def test_rdkit_lighadder(lighadder, 
                         input_ligands_moltype,
                         input_ligands_smitype):
    output = lighadder.execute(input_ligands_moltype)
    assert len(output["ligands"]) == len(input_ligands_moltype)
    assert all("H" in Chem.MolToSmiles(
        lig.get_representation_content(RDKitSmallMolRep)) 
        for lig in output["ligands"])
    assert all("H" not in Chem.MolToSmiles(
        lig.get_representation_content(RDKitSmallMolRep)) 
        for lig in input_ligands_moltype)
    assert output["ligands"].basic_itype == LigandData
    
    output = lighadder.execute(input_ligands_smitype)
    assert len(output["ligands"]) == len(input_ligands_smitype)
    assert all("H" in Chem.MolToSmiles(
        lig.get_representation_content(RDKitSmallMolRep)) 
        for lig in output["ligands"])
    assert output["ligands"].basic_itype == LigandData