from __future__ import annotations
from typing import Any, List, Tuple

import os
import pathlib
import subprocess

import numpy as np

from ...core.data.abstract import Representation
from ...core.data.representations import SMILESRep, PDBPathRep

try:
    from rdkit import Chem
    from ...modules.rdkit.representations import RDKitMolRep
    RDKIT_SUCCESSFUL_IMPORT = True
except: 
    RDKIT_SUCCESSFUL_IMPORT = False

class MOPACInputMolRep(Representation):
    rep_name = "mopac_input"

    def __init__(self, charge: int, coordinates: str,
                 setpi: List[Tuple[int, int]], neg_cvb: List[Tuple[int, int]],
                 desciption: str = "some description") -> None:
        data = {
            "charge": charge,
            "description": desciption,
            "coordinates": coordinates,
            "setpi": setpi,
            "neg_cvb": neg_cvb
        }
        super().__init__(data)
        self.charge = charge
        self.coordinates = coordinates
        self.description = desciption
        self.setpi = setpi
        self.neg_cvb = neg_cvb
    
    if RDKIT_SUCCESSFUL_IMPORT:
        @classmethod
        def from_RDKitMolRep(cls, rdkit_rep: RDKitMolRep):
            rdmol = Chem.Mol(rdkit_rep.data)
            coordinates = "\n".join(
                [line for line in Chem.MolToPDBBlock(rdmol).split("\n") 
                 if "ATOM" in line or "HETATM" in line])
            charge = Chem.GetFormalCharge(rdmol)
            setpi = cls.rdmol_to_setpi(rdmol)
            neg_cvb = cls.rdmol_to_neg_cvb(rdmol)
            return MOPACInputMolRep(
                charge=charge,
                coordinates=coordinates,
                setpi=setpi,
                neg_cvb=neg_cvb,
                desciption="Converted from rdkit molecule representation."
            )
        
        @classmethod
        def from_SMILESRep(cls, smiles_rep: SMILESRep):
            rdkit_rep = RDKitMolRep.from_SMILESRep(smiles_rep)
            return cls.from_RDKitMolRep(rdkit_rep)
        
        @classmethod
        def from_PDBPathRep(cls, pdb_rep: PDBPathRep):
            rdkit_rep = RDKitMolRep.from_PDBPathRep(pdb_rep)
            return cls.from_RDKitMolRep(rdkit_rep)
    
        @staticmethod
        def rdmol_to_setpi(rdmol: Chem.Mol) -> List[Tuple[int, int]]:
            setpi: List[Tuple[str, str]] = []
            rdmol = Chem.Mol(rdmol)
            Chem.Kekulize(rdmol)
            for bond in rdmol.GetBonds():
                bond: Chem.Bond
                if bond.GetBondType() == Chem.BondType.DOUBLE or \
                    bond.GetBondType() == Chem.BondType.TRIPLE:
                    setpi.append(
                        (
                            bond.GetBeginAtom().GetIdx()+1, 
                            bond.GetEndAtom().GetIdx()+1
                        )
                    )
                if bond.GetBondType() == Chem.BondType.TRIPLE:
                    setpi.append(
                        (
                            bond.GetBeginAtom().GetIdx()+1, 
                            bond.GetEndAtom().GetIdx()+1
                        )
                    )
            return setpi
        
        @staticmethod
        def rdmol_to_neg_cvb(rdmol: Chem.Mol) -> List[Tuple[int, int]]:
            neg_cvb_list = []
            for atom1 in range(rdmol.GetNumAtoms()):
                if not rdmol.GetAtomWithIdx(atom1).GetSymbol() == "O":
                    continue
                pos1 = rdmol.GetConformer().GetAtomPosition(atom1)
                for atom2 in range(atom1+1,rdmol.GetNumAtoms()):
                    if not rdmol.GetAtomWithIdx(atom2).GetSymbol() == "O":
                        continue
                    pos2 = rdmol.GetConformer().GetAtomPosition(atom2)
                    dist = np.linalg.norm(pos2 - pos1)
                    if dist < 3:
                        neg_cvb_list.append((atom1+1, atom2+1))
            return neg_cvb_list
    
    def update_coordinates(self, coords: np.ndarray):
        lines = self.coordinates.splitlines()
        newlines = []
        for loc, line in zip(coords, lines):
            x = "{:.3f}".format(loc[0]).rjust(6)
            y = "{:.3f}".format(loc[1]).rjust(6)
            z = "{:.3f}".format(loc[2]).rjust(6)
            line = line[:32] + x + "  " + y + "  " + z + line[32+22:]
            newlines.append(line)
        self.coordinates = "\n".join(newlines)


    
            
        
