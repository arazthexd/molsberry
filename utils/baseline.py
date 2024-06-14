from rdkit import Chem
from rdkit.Chem import rdDistGeom, rdForceFieldHelpers

from .ligprep import enumerate_protomers

DIMORPHITE_PH = 7.4
DIMORPHITE_PH_RANGE = 0.2

def get_atom_pkas(mol: Chem.Mol):
    pass

def get_molecule_base_tautomers(mol: Chem.Mol, mode: str = "Prop:InitialSmiles"):
    if mode.split(":")[0] == "Prop":
        return _get_molecule_base_tautomers_prop(mol, mode.split(":")[1])
    else:
        raise NotImplementedError()
    
def _get_molecule_base_tautomers_prop(mol: Chem.Mol, prop_name: str):
    smi = mol.GetProp(prop_name)
    tau = Chem.MolFromSmiles(smi)
    tau = Chem.AddHs(tau)
    rdDistGeom.EmbedMolecule(tau)
    rdForceFieldHelpers.MMFFOptimizeMolecule(tau)
    return [tau]

def get_molecule_base_protomers(mol: Chem.Mol, mode: str = "dimorphite_dl"):
    if mode.split(":")[0] == "dimorphite_dl":
        return _get_molecule_base_protomers_dimorphite(mol)
    else:
        raise NotImplementedError()

def _get_molecule_base_protomers_dimorphite(mol: Chem.Mol):
    if isinstance(mol, list): mol = mol[0] # TODO: THIS IS TEMPORARY
    pH_min = DIMORPHITE_PH - DIMORPHITE_PH_RANGE
    pH_max = DIMORPHITE_PH + DIMORPHITE_PH_RANGE
    if "_Name" not in mol.GetPropNames():
        mol.SetProp("_Name", "unnamed")
    protomers = enumerate_protomers(mol, pH_min, pH_max)
    return protomers

def get_molecule_baselines(mol: Chem.Mol, tmode: str = "Prop:InitialSmiles", pmode: str = "dimorphite_dl"):
    
    base_tautomers = get_molecule_base_tautomers(mol, tmode)
    base_protomers = get_molecule_base_protomers(base_tautomers, pmode)
    return base_protomers
