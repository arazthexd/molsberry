from __future__ import annotations

import parmed
from rdkit import Chem
import os

from ...core import MoleculeRep, PDBPathRep, generate_path_in_dir
from ..rdkit import RDKitMolRep

from openbabel import openbabel, pybel

BOND_ORDER_TO_RDBONDTYPE = {
    1.0: Chem.BondType.SINGLE,
    1.5: Chem.BondType.AROMATIC,
    2.0: Chem.BondType.DOUBLE,
    3.0: Chem.BondType.TRIPLE,
}

METAL_2P_ATOMIC_NUMS = [12, 25, 26, 27, 29, 30]

def check_atom_and_res_names(atom1, atom2, name1, name2, resname=None):
    if resname is not None:
        if atom1.residue.name != resname or atom2.residue.name != resname:
            return False
    return set([atom1.name, atom2.name]) == set([name1, name2])

class ParmedMolRep(MoleculeRep): # TODO: Parameterized vs NonParama
    rep_name = "Parm_mol"
    def __init__(self, structure: parmed.Structure):
        super().__init__(content=structure)
        self.content: parmed.Structure
    
    @property
    def parameterized(self) -> bool:
        return self.content.atoms[0].atom_type != \
            parmed.topologyobjects.UnassignedAtomType

    @staticmethod
    def rdatom2parmedatom(rdatom: Chem.Atom) -> parmed.Atom:
        return parmed.Atom(
            name=rdatom.GetSymbol()+str(rdatom.GetIdx()),
            formal_charge=rdatom.GetFormalCharge(),
            mass=rdatom.GetMass(),
            atomic_number=rdatom.GetAtomicNum()
        )
    
    @classmethod
    def from_RDKitMolRep(cls, rdrep: RDKitMolRep) -> ParmedMolRep:
        rdmol = rdrep.content
        stmol = cls.rdkit2parmed(rdmol)
        return cls(stmol)
    
    @classmethod
    def from_PDBPathRep(cls, pdbrep: PDBPathRep) -> ParmedMolRep:
        pdb = pdbrep.content
        st = cls.pdb2parmed(pdb)
        return cls(st)
    
    @staticmethod
    def pdb2parmed(pdb: str) -> parmed.Structure:
        st: parmed.Structure = parmed.load_file(pdb)
        
        for atom in st.atoms:
            atom: parmed.Atom

            if atom.atomic_number in METAL_2P_ATOMIC_NUMS:
                if atom.formal_charge is None or atom.formal_charge == 0:
                    atom.formal_charge = 2
                    print(f"WARNING: A metal atom from pdb file had unspecified charge. 2+ will be used, unless charge specified in the file.")
                    print("file:", pdb)
                    print("atom idx:", atom.idx)

            if atom.residue.name == "LYS" and atom.name == "NZ":
                atom.formal_charge = 1
            
            if atom.residue.name == "ARG" and atom.name == "NH1":
                atom.formal_charge = 1

            if atom.residue.name == "HIP" and atom.name == "ND1":
                atom.formal_charge = 1

            if atom.residue.name == "GLU" and atom.name == "OE2":
                atom.formal_charge = -1
            
            if atom.residue.name == "ASP" and atom.name == "OD2":
                atom.formal_charge = -1

            if atom.residue == st.residues[0]:
                if atom.residue.name not in ["NME", "ACE"]:
                    if atom.name == "N":
                        atom.formal_charge = 1
            
            if atom.residue == st.residues[-1]:
                if atom.residue.name not in ["NME", "ACE"]:
                    if atom.name == "OXT":
                        atom.formal_charge = -1
            
        
        for bond in st.bonds:
            bond: parmed.Bond
            atom1: parmed.Atom = bond.atom1
            atom2: parmed.Atom = bond.atom2
            res1: parmed.Residue = atom1.residue
            res2: parmed.Residue = atom2.residue

            # C=O for all AAs
            if check_atom_and_res_names(atom1, atom2, "C", "O"):
                bond.order = 2.0
            
            # NH1=CZ for ARG
            if check_atom_and_res_names(atom1, atom2, "NH1", "CZ", "ARG"):
                bond.order = 2.0

            # CG=OD1 for ASN
            if check_atom_and_res_names(atom1, atom2, "CG", "OD1", "ASN"):
                bond.order = 2.0

            # CG=OD1 for ASP
            if check_atom_and_res_names(atom1, atom2, "CG", "OD1", "ASP"):
                bond.order = 2.0

            # CD=OE1 for GLN
            if check_atom_and_res_names(atom1, atom2, "CD", "OE1", "GLN"):
                bond.order = 2.0

            # CD=OE1 for GLU
            if check_atom_and_res_names(atom1, atom2, "CD", "OE1", "GLU"):
                bond.order = 2.0

            # ND1=CE1 and CD2=CG for HIE, HIS, HIP
            if res1.name in ["HIE", "HIS", "HIP"]:
                if check_atom_and_res_names(atom1, atom2, "ND1", "CE1"):
                    bond.order = 2.0
                if check_atom_and_res_names(atom1, atom2, "CD2", "CG"):
                    bond.order = 2.0
            
            # NE2=CE1 and CD2=CG for HID (Histidine with delta nitrogen protonated)
            if check_atom_and_res_names(atom1, atom2, "NE2", "CE1", "HID"):
                bond.order = 2.0
            if check_atom_and_res_names(atom1, atom2, "CD2", "CG", "HID"):
                bond.order = 2.0

            # ND1=CE1 and CD2=CG for HIP (Histidine with both nitrogens protonated)
            if check_atom_and_res_names(atom1, atom2, "ND1", "CE1", "HIP"):
                bond.order = 2.0
            if check_atom_and_res_names(atom1, atom2, "CD2", "CG", "HIP"):
                bond.order = 2.0

            # Aromatic bonds for PHE (Phenylalanine) and TYR (Tyrosine)
            if res1.name in ["PHE", "TYR"]:
                aromatic_bonds = [
                    ("CG", "CD1"), ("CD1", "CE1"), ("CE1", "CZ"),
                    ("CZ", "CE2"), ("CE2", "CD2"), ("CD2", "CG")
                ]
                for atom_pair in aromatic_bonds:
                    if check_atom_and_res_names(atom1, atom2, atom_pair[0], atom_pair[1]):
                        bond.order = 1.5  # Aromatic bond order

            # Aromatic bonds for TRP (Tryptophan) + CG=CD1
            if res1.name == "TRP":
                aromatic_bonds = [
                    ("CD2", "CE2"), ("CE2", "CZ2"), ("CZ2", "CH2"),
                    ("CH2", "CZ3"), ("CZ3", "CE3"), ("CE3", "CD2")
                ]
                for atom_pair in aromatic_bonds:
                    if check_atom_and_res_names(atom1, atom2, atom_pair[0], atom_pair[1]):
                        bond.order = 1.5  # Aromatic bond order
                if check_atom_and_res_names(atom1, atom2, "CG", "CD1"):
                    bond.order = 2
        return st
            

    def to_RDKitMolRep(self) -> RDKitMolRep:
        stmol: parmed.Structure = self.content
        rdmoln = self.parmed2rdkit(stmol)
        return RDKitMolRep(rdmoln)
    
    @staticmethod
    def rdkit2parmed(rdmol: Chem.Mol) -> parmed.Structure:
        stmol = parmed.Structure()
        coords = rdmol.GetConformer().GetPositions()

        for i, rdatom in enumerate(rdmol.GetAtoms()):
            rdatom: Chem.Atom
            statom = ParmedMolRep.rdatom2parmedatom(rdatom)
            pdbinfo = rdatom.GetPDBResidueInfo()
            if pdbinfo is None:
                stmol.add_atom(statom, resname="UNK", resnum=1)
            else:
                stmol.add_atom(statom, 
                               pdbinfo.GetResidueName(), 
                               pdbinfo.GetResidueNumber())

        stmol.coordinates = coords

        # TODO: Kekulize?
        for i, rdbond in enumerate(rdmol.GetBonds()):
            rdbond: Chem.Bond
            at1, at2 = rdbond.GetBeginAtomIdx(), rdbond.GetEndAtomIdx()
            stbond = parmed.Bond(stmol.atoms[at1], stmol.atoms[at2], 
                                order=rdbond.GetBondTypeAsDouble())
            stmol.bonds.append(stbond)

        res_number_list = []
        for res in stmol.residues:
            res: parmed.Residue
            res_num = res.number
            if res_num not in res_number_list:
                res_number_list.append(res.number)
            else:
                print(stmol.residues[res_num])
                atoms = res.atoms
                [stmol.atoms.remove(atom) for atom in atoms]
                [stmol.residues[res_num-1].add_atom(atom) for atom in atoms]
                
        return stmol
    
    @staticmethod
    def parmed2rdkit(stmol: parmed.Structure) -> Chem.Mol:
        rdmoln = Chem.RWMol()
        for i, atom in enumerate(stmol.atoms):
            atom: parmed.Atom
            rdatom = Chem.Atom(atom.atomic_number)
            rdatom.SetFormalCharge(atom.formal_charge)
            pdbinfo = Chem.AtomPDBResidueInfo(
                atomName=atom.name,
                serialNumber=i,
                residueName=atom.residue.name,
                residueNumber=atom.residue.number,
                chainId=atom.residue.chain
            )
            rdatom.SetPDBResidueInfo(pdbinfo)
            rdmoln.AddAtom(rdatom)
        [rdmoln.AddBond(bond.atom1.idx, bond.atom2.idx,
                        order=BOND_ORDER_TO_RDBONDTYPE.get(bond.order)) for bond in stmol.bonds]
        conformer = Chem.Conformer(rdmoln.GetNumAtoms())
        conformer.SetPositions(stmol.coordinates)
        rdmoln.AddConformer(conformer)
        rdmoln = rdmoln.GetMol()
        Chem.SanitizeMol(rdmoln)
        rdmoln.UpdatePropertyCache()
        return rdmoln
    