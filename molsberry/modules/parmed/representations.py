from __future__ import annotations

import parmed
from rdkit import Chem

from ...core import MoleculeRep, PDBPathRep
from ..rdkit import RDKitMolRep

BOND_ORDER_TO_RDBONDTYPE = {
    1.0: Chem.BondType.SINGLE,
    1.5: Chem.BondType.AROMATIC,
    2.0: Chem.BondType.DOUBLE,
    3.0: Chem.BondType.TRIPLE,
}

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
            charge=rdatom.GetFormalCharge(),
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
        stmol = parmed.load_file(pdb)
        return cls(stmol)

    def to_RDKitMolRep(self) -> RDKitMolRep:
        stmol: parmed.Structure = self.content
        rdmoln = self.parmed2rdkit(stmol)
        return RDKitMolRep(rdmoln)
    
    @staticmethod
    def rdkit2parmed(rdmol: Chem.Mol) -> parmed.Structure:
        stmol = parmed.Structure()
        coords = rdmol.GetConformer().GetPositions()

        for i, rdatom in enumerate(rdmol.GetAtoms()):
            statom = ParmedMolRep.rdatom2parmedatom(rdatom)
            stmol.add_atom(statom, resname="UNK", resnum=1)

        stmol.coordinates = coords

        for i, rdbond in enumerate(rdmol.GetBonds()):
            rdbond: Chem.Bond
            at1, at2 = rdbond.GetBeginAtomIdx(), rdbond.GetEndAtomIdx()
            stbond = parmed.Bond(stmol.atoms[at1], stmol.atoms[at2], 
                                order=rdbond.GetBondTypeAsDouble())
            stmol.bonds.append(stbond)

        return stmol
    
    @staticmethod
    def parmed2rdkit(stmol: parmed.Structure) -> Chem.Mol:
        rdmoln = Chem.RWMol()
        for atom in stmol.atoms:
            atom: parmed.Atom
            rdatom = Chem.Atom(atom.atomic_number)
            rdatom.SetFormalCharge(int(atom.charge))
            rdmoln.AddAtom(rdatom)
        [rdmoln.AddBond(bond.atom1.idx, bond.atom2.idx,
                        order=BOND_ORDER_TO_RDBONDTYPE.get(bond.order)) for bond in stmol.bonds]
        conformer = Chem.Conformer(rdmoln.GetNumAtoms())
        conformer.SetPositions(stmol.coordinates)
        rdmoln.AddConformer(conformer)
        rdmoln = rdmoln.GetMol()
        return rdmoln
    