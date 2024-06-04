import warnings

from rdkit import Chem

def read_mol_file(file_path: str, out_type: str):
    assert isinstance(file_path, str)
    if file_path.endswith(".smi"):
        if out_type == "rdmol":
            mols = _read_smi_to_rdmol(file_path)
            return mols
        else:
            raise NotImplementedError()
    elif file_path.endswith(".sdf"):
        if out_type == "rdmol":
            mols = _read_sdf_to_rdmol(file_path)
            return mols
        else:
            raise NotImplementedError()
    else:
        raise NotImplementedError()

def write_mol_file(input, out_path: str):
    if out_path.endswith(".sdf"):
        if isinstance(input, list):
            assert len(input) > 0
            if isinstance(input[0], Chem.Mol):
                _write_rdmols_to_sdf(input, out_path)
            else:
                raise NotImplementedError()
        elif isinstance(input, Chem.Mol):
            _write_rdmols_to_sdf([input], out_path)
        else:
            raise NotImplementedError()
    else:
        raise NotImplementedError()
    
def _read_smi_to_rdmol(file_path): # TODO: multiprocessing implementation for large libraries
    mols = []
    i = 0
    with open(file_path, "r") as f:
        for line in f.readlines():
            split_line = line.split()
            if len(split_line) > 0:
                mols.append(Chem.MolFromSmiles(split_line[0]))
                if len(split_line) > 1:
                    mols[-1].SetProp("_Name", split_line[1])
                else:
                    mols[-1].SetProp("_Name", f"unnamed{i}")
                    i += 1
    return mols

def _read_sdf_to_rdmol(file_path): # TODO: multiprocessing implementation for large libraries
    no_name_warn = True
    mols = []
    suppl = Chem.SDMolSupplier(file_path)
    for mol in suppl:
        mols.append(mol)
        if not no_name_warn:
            continue
        if mol.GetProp("_Name") == None: # TODO: This needs checking...
            warnings.warn("Molecules with no name observed...")
            no_name_warn = False
    return mols

def _write_rdmols_to_sdf(rdmol_list, out_path): # TODO: multiprocessing implementation for large libraries
    sdwriter = Chem.SDWriter(out_path)
    for mol in rdmol_list:
        sdwriter.write(mol)