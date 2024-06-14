from rdkit import Chem

def CalcNumQueryMatches(mol: Chem.Mol, query: Chem.Mol) -> int:
    mol = Chem.RemoveHs(mol)
    return len(mol.GetSubstructMatches(query, uniquify=True))


