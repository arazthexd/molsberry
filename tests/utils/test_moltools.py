# import pytest
# from moddipic.utils.moltools import (
#     is_mol_3d,
#     sync_mol_flexible_rotors,
#     sync_flexible_rotors,
#     get_rotatable_dihedrals,
# )

# from rdkit import Chem
# from rdkit.Chem import rdDepictor, rdDistGeom

# testing_mols = [
#     {
#         "smiles": "COCc1cnc(C)cc1Cl",
#         "rotors": [(0, 1, 2, 3), (1, 2, 3, 4)],
#     }, 
#     {
#         "smiles": "CCCC",
#         "rotors": [(0, 1, 2, 3)]
#     }
# ]

# @pytest.fixture(params=testing_mols)
# def base_mol(request):
#     yield Chem.MolFromSmiles(request.param["smiles"])

# @pytest.fixture(params=testing_mols) # type: ignore
# def base_mol_rotors(request):
#     yield (
#         Chem.MolFromSmiles(request.param["smiles"]), request.param["rotors"]
#     )

# def test_is_mol_3d(base_mol):
#     mol = Chem.Mol(base_mol)

#     assert is_mol_3d(mol) == False

#     rdDistGeom.EmbedMolecule(mol)
#     assert is_mol_3d(mol) == True
    
#     rdDepictor.Compute2DCoords(mol)
#     assert is_mol_3d(mol) == False

# def test_get_rotatable_dihedrals(base_mol_rotors):
#     mol = Chem.Mol(base_mol_rotors[0])

#     assert get_rotatable_dihedrals(mol) == base_mol_rotors[1]


# # TODO: Tests for sync mol methods...

