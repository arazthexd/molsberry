# from moddipic.modules.openmm.representations import OpenMMComponentRep
# from moddipic.modules.rdkit.representations import RDKitMolRep

# import pytest
# from rdkit import Chem
# from rdkit.Chem import rdDistGeom

# @pytest.fixture
# def rdkit_mol_rep():
#     rdmoll = Chem.AddHs(Chem.MolFromSmiles("CCCCOC"))
#     rdDistGeom.EmbedMolecule(rdmoll)
#     return RDKitMolRep(rdmoll)

# def test_openmm_from_rdkit(rdkit_mol_rep):
#     comp_rep = OpenMMComponentRep.from_RDKitMolRep(rdkit_mol_rep)
