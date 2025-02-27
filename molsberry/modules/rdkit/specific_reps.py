from .representations import RDKitMolRep, SmallMolRep, ProteinRep, PDBPathRep, List, Representation, Chem # :D

class RDKitSmallMolRep(RDKitMolRep, SmallMolRep):
    def save_rep(self, exless_filename: str):
        rep_path = exless_filename + ".sdf"
        writer = Chem.SDWriter(rep_path)
        writer.write(self.content)
        writer.close()
    
    @classmethod
    def save_rep_batch(cls, reps: List[Representation], exless_filename: str):
        rep_path = exless_filename + ".sdf"
        writer = Chem.SDWriter(rep_path)
        for rep in reps:
            writer.write(rep.content)
        writer.close()

from ..parmed import ParmedMolRep

class RDKitProtRep(RDKitMolRep, ProteinRep):
    rep_name = "rdprot"
    
    @classmethod
    def from_PDBPathRep(cls, pdb_rep: PDBPathRep):
        assert isinstance(pdb_rep, PDBPathRep)
        pdb_path = pdb_rep.content
        mol = ParmedMolRep.parmed2rdkit(ParmedMolRep.pdb2parmed(pdb_path))

        # Fix COO groups not being ionized when read from pdb in rdkit
        # TODO: Probably needs more attention
        query_COO = Chem.MolFromSmarts("[$([O]-C(=O)-C)]")
        for atom, in mol.GetSubstructMatches(query_COO):
            atom: Chem.Atom = mol.GetAtomWithIdx(atom)
            if atom.GetPDBResidueInfo().GetIsHeteroAtom() == False:
                atom.SetFormalCharge(-1)

        return cls(mol=mol)
