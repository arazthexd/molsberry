from typing import Any, List
import os, pathlib, glob

from moddipic.core.data.abstract import Representation
from moddipic.core.data.representations import PDBPathProteinRep

from moddipic.modules.rdkit.representations import RDKitMolRep
from rdkit import Chem
# TODO: Find another method to not import everything if not available

from .constants import CUBY4_TMP_PATH

class Cuby4InputMoleculeRep(Representation):
    rep_name = "cuby4_geometry"

    def __init__(self, geometry_path: str, charge: int, setpi: List[str],
                 setcharge: List[str], cvb: List[str]):
        self.geometry_path: str = geometry_path # NOTE: Not saved right away.
        self.charge: int = charge
        self.setpi: List[str] = setpi
        self.setcharge: List[str] = setcharge
        self.cvb: List[str] = cvb
    
    @classmethod 
    def from_RDKitMolRep(cls, rdmol_rep: RDKitMolRep):
        rdmol: Chem.Mol = Chem.Mol(rdmol_rep.data)
        charge = Chem.GetFormalCharge(rdmol)
        setpi = cls.rdmol_to_setpi(rdmol)
        setcharge = cls.rdmol_to_setcharge(rdmol)
        cvb = cls.rdmol_to_cvb(rdmol)
        geometry_path = cls.get_filename() + ".sdf" # TODO: Fix this...
        return cls(geometry_path, charge, setpi, setcharge, cvb)
    
    @classmethod
    def from_PDBPathProteinRep(cls, pdb_rep: PDBPathProteinRep):
        max_attempts = 5
        pdb_path = pdb_rep.data
        setpi = []
        for i in range(max_attempts):
            try:
                protein_top: Topology = Topology.from_pdb(pdb_path)
                setcharge = [
                    f"{i+1}: {a.formal_charge.m_as(unit.elementary_charge)}" 
                    for i, a in enumerate(protein_top.atoms)
                ]
            except: pass
        

    @staticmethod
    def get_filename():
        while True:
            prefix = "CB4_Geom_"
            rand_str = Representation.get_filename()
            full_prefix = os.path.join(
                str(pathlib.Path(CUBY4_TMP_PATH).absolute()),
                prefix+rand_str
            )
            if not glob.glob(full_prefix+"*"):
                return full_prefix

    @staticmethod
    def rdmol_to_setpi(rdmol: Chem.Mol) -> List[str]:
        setpi = []
        no_res_mol = Chem.Mol(rdmol)
        Chem.Kekulize(no_res_mol)
        for bond in no_res_mol.GetBonds():
            if bond.GetBondType() == Chem.BondType.DOUBLE or \
                bond.GetBondType() == Chem.BondType.TRIPLE:
                setpi.append(f"{bond.GetBeginAtom().GetIdx()+1};{bond.GetEndAtom().GetIdx()+1}")
            if bond.GetBondType() == Chem.BondType.TRIPLE:
                setpi.append(f"{bond.GetBeginAtom().GetIdx()+1};{bond.GetEndAtom().GetIdx()+1}")
        return setpi
    
    @staticmethod
    def rdmol_to_setcharge(rdmol: Chem.Mol) -> List[str]:
        def helper(x):
            if x == 0: return str(x)
            if x > 0:
                return str("'+'") * x
            if x < 0:
                return str("'-'") * (-x)
        setcharge = [f"{a.GetIdx()+1}: {helper(a.GetFormalCharge())}" 
                    for a in rdmol.GetAtoms()]
        return setcharge
    
    @staticmethod
    def rdmol_to_cvb(rdmol: Chem.Mol) -> Any:
        cvb_list = []
        for atom1 in range(rdmol.GetNumAtoms()):
            for atom2 in range(atom1+1,rdmol.GetNumAtoms()):
                if not rdmol.GetBondBetweenAtoms(atom1, atom2):
                    if rdmol.GetAtomWithIdx(atom1).GetSymbol() == "O" and \
                        rdmol.GetAtomWithIdx(atom2).GetSymbol() == "O":
                        cvb_list.append(f"{atom1+1}:-{atom2+1}")

        return cvb_list
    
    def __repr__(self) -> str:
        repr_text = f"<Cuby4InputMoleculeRep(path='{self.geometry_path}', " + \
            f"charge={self.charge})>"
        return repr_text