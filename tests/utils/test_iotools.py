# import pytest, os, pathlib
# from moddipic.utils.iotools import (
#     write_ligands,
#     load_ligands
# )
# from moddipic.utils.moltools import is_mol_3d
# from tests.utils.example_mols import testing_mols

# from rdkit import Chem

# sdf_path = os.path.join(pathlib.Path(__file__).parent, "test_mols.sdf")

# def setup_module():
#     print(testing_mols)
#     mols = [Chem.MolFromSmiles(m["smiles"]) for m in testing_mols]
#     writer = Chem.SDWriter(sdf_path)
#     [writer.write(mol) for mol in mols]
#     writer.close()

# def teardown_module():
#     if os.path.exists(sdf_path):
#         os.remove(sdf_path)

# def test_load_ligands():
#     ligands = load_ligands(sdf_path)
#     for lig, ref in zip(ligands, testing_mols):
#         assert Chem.MolToSmiles(lig) == \
#             Chem.MolToSmiles(Chem.MolFromSmiles(ref["smiles"]))
#         assert is_mol_3d(lig) == False
    
#     ligands = load_ligands(sdf_path, final_3d=True)
#     for lig, ref in zip(ligands, testing_mols):
#         assert Chem.MolToSmiles(lig) == \
#             Chem.MolToSmiles(Chem.MolFromSmiles(ref["smiles"]))
#         assert is_mol_3d(lig) == True
    
#     ligands = load_ligands(sdf_path, final_3d=True, addHs=True)
#     for lig, ref in zip(ligands, testing_mols):
#         assert Chem.MolToSmiles(lig) == \
#             Chem.MolToSmiles(Chem.AddHs(Chem.MolFromSmiles(ref["smiles"])))
#         assert is_mol_3d(lig) == True
    
#     ligands = load_ligands(sdf_path, final_3d=False, addHs=True)
#     for lig, ref in zip(ligands, testing_mols):
#         assert Chem.MolToSmiles(lig) == \
#             Chem.MolToSmiles(Chem.AddHs(Chem.MolFromSmiles(ref["smiles"])))
#         assert is_mol_3d(lig) == False

# # TODO: Add iotools module tests...