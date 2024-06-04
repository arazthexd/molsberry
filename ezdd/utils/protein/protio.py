from rdkit import Chem
from pdbfixer import PDBFixer
from biopandas.pdb import PandasPdb
from openmm.app import PDBFile

def read_prot_file(file_path: str, out_type: str):
    assert isinstance(file_path, str)
    if file_path.endswith(".pdb"):
        if out_type == "pdbfixer":
            protein = _read_pdb_to_fixer(file_path)
        elif out_type == "ppdb":
            protein = _read_pdb_to_ppdb(file_path)
        elif out_type == "openmm":
            protein = PDBFile(file_path)
        else:
            raise NotImplementedError()
    else:
        raise NotImplementedError()
    
    return protein
    
def _read_pdb_to_fixer(file_path):
    return PDBFixer(file_path)

def _read_pdb_to_ppdb(file_path):
    return PandasPdb().read_pdb(file_path)

