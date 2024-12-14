from typing import List, Dict
import os

from rdkit import Chem

from ....core import (
    SimpleBlock, ProteinData, PDBPathRep, LocationData, Representation,
    generate_random_str, LigandData
)
from ...rdkit import RDKitProtRep, RDKitMolRep
from .locationrep import PocketLocationRep, PocketLocation

def nterm(emol: Chem.EditableMol, atom: Chem.Atom):
    aname = atom.GetPDBResidueInfo().GetName().strip()

    if aname in ["O", "C"]:
        pdbinfo = atom.GetPDBResidueInfo()
        pdbinfo.SetResidueName("ACE")
        emol.ReplaceAtom(atom.GetIdx(), atom)
        return
    
    emol.RemoveAtom(atom.GetIdx())
    return

def cterm(emol: Chem.EditableMol, atom: Chem.Atom):
    aname = atom.GetPDBResidueInfo().GetName().strip()

    if aname in ["N", "CA"]:
        pdbinfo = atom.GetPDBResidueInfo()
        pdbinfo.SetResidueName("NME")
        emol.ReplaceAtom(atom.GetIdx(), atom)
        return
    
    emol.RemoveAtom(atom.GetIdx())
    return

class RDKitPocketIsolator(SimpleBlock):
    name = "rdpocketisolator"
    display_name = "(RDKit) Pocket Isolator"
    inputs = [
        ("protein", ProteinData, RDKitProtRep, False),
        ("location", LocationData, PocketLocationRep, False)
    ]
    outputs = [
        ("pocket", ProteinData, PDBPathRep, False)
    ]
    batch_groups = [("protein", "location")]

    def operate(self, input_dict: Dict[str, Representation]) \
        -> Dict[str, Representation]:

        rdprot = input_dict[self.input_keys[0]].content
        loc = input_dict[self.input_keys[1]].content
        rdpoc = self.isolate(rdprot, loc)
        main_out_key = self.output_keys[0]
        rdpoc_pdb = generate_random_str(6)
        rdpoc_pdb = os.path.join(self.base_dir, f"{rdpoc_pdb}.pdb")
        Chem.MolToPDBFile(rdpoc, rdpoc_pdb)
        return {main_out_key: self._get_out_rep(main_out_key)(rdpoc_pdb)}

    def isolate(self, rdprot: Chem.Mol, loc: PocketLocation) -> Chem.Mol:
        remres = self.get_remaining_residues(rdprot, loc)
        terminal = self.identify_terminal_cappings(rdprot, remres)

        rdprote = Chem.EditableMol(rdprot)
        rdprote.BeginBatchEdit()
        for atom in rdprot.GetAtoms():
            atom: Chem.Atom
            rn = atom.GetPDBResidueInfo().GetResidueNumber()
            if rn in remres:
                if terminal[remres.index(rn)] == "N":
                    nterm(rdprote, atom)
                if terminal[remres.index(rn)] == "C":
                    cterm(rdprote, atom)
            else:
                rdprote.RemoveAtom(atom.GetIdx())
        rdprote.CommitBatchEdit()

        rdprotee = rdprote.GetMol()
        Chem.SanitizeMol(rdprotee)
        return rdprotee
    
    @staticmethod
    def get_remaining_residues(rdprot: Chem.Mol, loc: PocketLocation):
        remres = set()
        pos = rdprot.GetConformer().GetPositions()
        allres = set(atom.GetPDBResidueInfo().GetResidueNumber() 
                    for atom in rdprot.GetAtoms())
        for atom in rdprot.GetAtoms():
            atom: Chem.Atom
            if atom.GetPDBResidueInfo().GetIsHeteroAtom():
                continue
            if loc.point_is_included(pos[atom.GetIdx()]):
                rn = atom.GetPDBResidueInfo().GetResidueNumber()
                remres.add(rn)
                if rn-1 in allres:
                    remres.add(rn-1)
                if rn+1 in allres:
                    remres.add(rn+1)
        return sorted(remres)
    
    @staticmethod
    def identify_terminal_cappings(rdprot: Chem.Mol, remres: List[int]):
        if remres[0] == 1:
            terminal = [None]
        else:
            terminal = ["N"]

        for i, resnum in enumerate(remres[1:-1], start=1):
            prev_gap = (resnum - remres[i-1])
            next_gap = (remres[i+1] - resnum)
            
            if prev_gap == 1 and next_gap == 1:
                terminal.append(None)
                continue

            if prev_gap > 1 and next_gap > 1:
                terminal.append(None)
                continue

            if prev_gap > 1:
                terminal.append("N")
                continue

            if next_gap > 1:
                terminal.append("C")
                continue
        
        allres = set(atom.GetPDBResidueInfo().GetResidueNumber() 
                    for atom in rdprot.GetAtoms())
        if remres[-1] == max(allres):
            terminal.append(None)
        else:
            terminal.append("C")
        
        return terminal
    
class RDKitLigandPocketLocator(SimpleBlock):
    name = "rdligpoclocator"
    display_name = "(RDKit) Ligand Pocket Locator"
    inputs = [
        ("ligand", LigandData, RDKitMolRep, False)
    ]
    outputs = [
        ("location", LocationData, PocketLocationRep, False)
    ]
    batch_groups = []

    def __init__(self, radius: float, 
                 debug = False, save_output = False, num_workers = None):
        super().__init__(debug, save_output, num_workers)
        self.radius = radius

    def operate(self, input_dict: Dict[str, Representation]) \
        -> Dict[str, Representation]:

        rdmol = input_dict[self.input_keys[0]].content
        main_out_key = self.output_keys[0]
        loc = PocketLocation("ligand", ligand=rdmol, radius=self.radius)
        return {main_out_key: self._get_out_rep(main_out_key)(loc)}