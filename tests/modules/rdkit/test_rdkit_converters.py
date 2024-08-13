import pytest
import os

from rdkit import Chem

from molsberry.modules.rdkit.converters import (
    RDKitLigandEmbedder,
    RDKitLigandHAdder,
    RDKitMMFFOptimizer,
    RDKitPLComplexOptimizer
)
from molsberry.core.data import (
    LigandData, BatchedData, MoleculeData, ProteinData
)
from molsberry.core.data import SMILESRep
from molsberry.modules.rdkit.representations import (
    RDKitSmallMolRep, PDBPathRep
)
from molsberry.modules.rdkit.interface import RDKitInterface

@pytest.fixture
def ligembedder():
    return RDKitLigandEmbedder(debug=True, save_output=False)

@pytest.fixture
def lighadder():
    return RDKitLigandHAdder(debug=True, save_output=False)

@pytest.fixture
def ligoptimizer():
    return RDKitMMFFOptimizer(debug=True, save_output=False)

@pytest.fixture
def ploptimizer():
    print(RDKitPLComplexOptimizer().input_keys)
    return RDKitPLComplexOptimizer(debug=True, save_output=False)

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

@pytest.fixture
def input_lig_from_file():
    sdf_path = os.path.join(os.path.dirname(__file__), 
                            "files", "lig.sdf")
    rdmol = next(Chem.SDMolSupplier(sdf_path, removeHs=False))
    return LigandData(RDKitSmallMolRep(rdmol))

@pytest.fixture
def input_prot_from_file():
    pdb_path = os.path.join(os.path.dirname(__file__), 
                            "files", "rec.pdb")
    return ProteinData(PDBPathRep(pdb_path))

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

def test_rdkit_ligoptimizer(ligoptimizer, 
                            input_lig_from_file):
    output = ligoptimizer.execute(input_lig_from_file)
    assert len(output["molecules"]) == 1
    assert all("H" in Chem.MolToSmiles(
        lig.get_representation_content(RDKitSmallMolRep)) 
        for lig in output["molecules"])
    assert output["molecules"].basic_itype == MoleculeData

def test_rdkit_ploptimizer(ploptimizer, 
                           input_lig_from_file,
                           input_prot_from_file):
    output = ploptimizer.execute(input_lig_from_file, input_prot_from_file)
    assert len(output["ligand"]) == 1
    assert all("H" in Chem.MolToSmiles(
        lig.get_representation_content(RDKitSmallMolRep)) 
        for lig in output["ligand"])
    assert output["ligand"].basic_itype == LigandData

# TODO: Write tests for RDKit Optimizer.