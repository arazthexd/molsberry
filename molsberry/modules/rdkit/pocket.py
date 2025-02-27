from typing import List, Dict
import os

import numpy as np

from rdkit import Chem
from pdbfixer import pdbfixer

from ...core import (
    SimpleBlock, ProteinData, PDBPathRep, LocationData, Representation,
    generate_random_str, MoleculeData, generate_path_in_dir
)
from . import RDKitProtRep, RDKitMolRep
from ...core import PocketLocationRep, PocketLocation

PROT_RESNAMES = (pdbfixer.proteinResidues +
                 list(pdbfixer.substitutions.keys()) +
                 ["HIE", "HID", "GLH", "CYX"]) # TODO: what else

class RDKitLigandPocketLocator(SimpleBlock):
    name = "rdligpoclocator"
    display_name = "(RDKit) Ligand Pocket Locator"
    inputs = [
        ("ligand", MoleculeData, RDKitMolRep, False)
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

def nterm(emol: Chem.RWMol, atom: Chem.Atom):
    aname = atom.GetPDBResidueInfo().GetName().strip()
    pdbinfo = atom.GetPDBResidueInfo()
    if pdbinfo.GetResidueName() == "ACE":
        return

    pdbinfo.SetResidueName("ACE")
    pdbinfo.SetIsHeteroAtom(False)

    coords = emol.GetConformer().GetPositions()
    
    if aname in ["O", "C"]:
        emol.ReplaceAtom(atom.GetIdx(), atom)
        return
    
    if aname == "CA":
        newname = "CH3"
        pdbinfo.SetName(newname)
        emol.ReplaceAtom(atom.GetIdx(), atom)
        return
    
    if aname in ["HA", "HA2", "HA3"]:
        newname = {"HA": "HH31", "HA2": "HH31", "HA3": "HH32"}.get(aname)
        pdbinfo.SetName(newname)
        emol.ReplaceAtom(atom.GetIdx(), atom)
        return

    if aname in ["CB", "N"]:
        newname = {"CB": "HH32", "N": "HH33"}.get(aname)
        newatom = Chem.Atom(1)
        pdbinfo.SetName(newname)
        newatom.SetPDBResidueInfo(pdbinfo)

        coords_h: np.ndarray = coords[atom.GetIdx()]
        atom_ca = [atom for atom in atom.GetNeighbors() 
                   if atom.GetPDBResidueInfo().GetName().strip() in ["CA", 
                                                                     "CH3"]][0]
        coords_ca: np.ndarray = coords[atom_ca.GetIdx()]
        vec = (coords_h - coords_ca) * 1.08 / np.linalg.norm(coords_h - coords_ca) 
        coords_h_new = coords_ca + vec
        coords[atom.GetIdx()] = coords_h_new
        emol.GetConformer().SetPositions(coords)

        emol.ReplaceAtom(atom.GetIdx(), newatom)
        return

    emol.RemoveAtom(atom.GetIdx())
    return

def cterm(emol: Chem.RWMol, atom: Chem.Atom):
    aname = atom.GetPDBResidueInfo().GetName().strip()
    pdbinfo = atom.GetPDBResidueInfo()
    if pdbinfo.GetResidueName() == "NME":
        return
    
    resname = pdbinfo.GetResidueName()
    pdbinfo.SetResidueName("NME")
    pdbinfo.SetIsHeteroAtom(False)

    coords = emol.GetConformer().GetPositions()

    if aname in ["N", "H"]:
        emol.ReplaceAtom(atom.GetIdx(), atom)
        return
    
    if aname == "CA":
        newname = "CH3"
        pdbinfo.SetName(newname)
        emol.ReplaceAtom(atom.GetIdx(), atom)
        return
    
    if aname in ["HA", "HA2", "HA3"]:
        newname = {"HA": "HH31", "HA2": "HH31", "HA3": "HH32"}.get(aname)
        pdbinfo.SetName(newname)
        emol.ReplaceAtom(atom.GetIdx(), atom)
        return

    if aname in ["CB", "C"]:
        newname = {"CB": "HH32", "C": "HH33"}.get(aname)
        newatom = Chem.Atom(1)
        pdbinfo.SetName(newname)
        newatom.SetPDBResidueInfo(pdbinfo)

        coords_h: np.ndarray = coords[atom.GetIdx()]
        atom_ca = [atom for atom in atom.GetNeighbors() 
                   if atom.GetPDBResidueInfo().GetName().strip() in ["CA", 
                                                                     "CH3"]][0]
        coords_ca: np.ndarray = coords[atom_ca.GetIdx()]
        vec = (coords_h - coords_ca) * 1.08 / np.linalg.norm(coords_h - coords_ca) 
        coords_h_new = coords_ca + vec
        coords[atom.GetIdx()] = coords_h_new
        emol.GetConformer().SetPositions(coords)

        emol.ReplaceAtom(atom.GetIdx(), newatom)
        return
    
    if aname == "CD" and resname == "PRO":
        newname = "H"
        newatom = Chem.Atom(1)
        pdbinfo.SetName(newname)
        newatom.SetPDBResidueInfo(pdbinfo)

        coords_h: np.ndarray = coords[atom.GetIdx()]
        atom_n = [atom for atom in atom.GetNeighbors() 
                   if atom.GetPDBResidueInfo().GetName().strip() in ["N"]][0]
        coords_n: np.ndarray = coords[atom_n.GetIdx()]
        vec = (coords_h - coords_n) * 1.08 / np.linalg.norm(coords_h - coords_n) 
        coords_h_new = coords_n + vec
        coords[atom.GetIdx()] = coords_h_new
        emol.GetConformer().SetPositions(coords)

        emol.ReplaceAtom(atom.GetIdx(), newatom)
        return

    
    emol.RemoveAtom(atom.GetIdx())
    return

def cyx_fix(emol: Chem.RWMol, atom: Chem.Atom, loc: PocketLocation):
    if atom.GetSymbol() != "S":
        return
    
    if atom.GetPDBResidueInfo().GetResidueName() != "CYX":
        return
    
    ref_s = atom
    ref_s_rn = atom.GetPDBResidueInfo().GetResidueNumber()
    coords = emol.GetConformer().GetPositions()

    if not any(loc.point_is_included(coords[at.GetIdx()]) 
               for at in emol.GetAtoms()
               if at.GetPDBResidueInfo().GetResidueNumber() == ref_s_rn):
        return
    
    nei_s = [nei for nei in ref_s.GetNeighbors() if nei.GetSymbol() == "S"]
    if nei_s.__len__() == 0:
        print("WARNING: Potential Problem in ss_removal in pocket isolation")
        return
    nei_s: Chem.Atom = nei_s[0]

    nei_s_rn = nei_s.GetPDBResidueInfo().GetResidueNumber()
    if any(loc.point_is_included(coords[at.GetIdx()]) 
           for at in emol.GetAtoms()
           if at.GetPDBResidueInfo().GetResidueNumber() == nei_s_rn):
        return
    
    for at in emol.GetAtoms():
        at: Chem.Atom
        pi = at.GetPDBResidueInfo()
        if pi.GetResidueNumber() == ref_s_rn:
            pi.SetResidueName("CYS")
            at.SetPDBResidueInfo(pi)
            emol.ReplaceAtom(at.GetIdx(), at)
    
    ref_s_pdbinfo = ref_s.GetPDBResidueInfo()
    nei_h = Chem.Atom(1)
    nei_pi = Chem.AtomPDBResidueInfo(
        atomName="HG",
        residueName="CYS",
        residueNumber=ref_s_pdbinfo.GetResidueNumber(),
        chainId=ref_s_pdbinfo.GetChainId(),
    )
    nei_h.SetPDBResidueInfo(nei_pi)
    nei_h_idx = emol.ReplaceAtom(nei_s.GetIdx(), nei_h)
    
    atom_coords: np.ndarray = coords[atom.GetIdx()]
    nei_s_coords = emol.GetConformer().GetPositions()[nei_s.GetIdx()]
    vec = (nei_s_coords - atom_coords) * 1.37 / np.linalg.norm(nei_s_coords - atom_coords) 
    coords_h_new = atom_coords + vec
    # coords = np.concatenate([coords, coords_h_new.reshape(1, 3)])
    coords[nei_s.GetIdx()] = coords_h_new
    emol.GetConformer().SetPositions(coords)
    return


class RDKitPocketIsolator(SimpleBlock):
    name = "rdpocketisolator"
    display_name = "(RDKit) Pocket Isolator"
    inputs = [
        ("protein", MoleculeData, RDKitProtRep, False),
        ("location", LocationData, PocketLocationRep, False)
    ]
    outputs = [
        ("pocket", MoleculeData, RDKitMolRep, False) # Before: PDBPathRep
    ]
    batch_groups = [("protein", "location")]

    def operate(self, input_dict: Dict[str, Representation]) \
        -> Dict[str, Representation]:

        rdprot = input_dict[self.input_keys[0]].content
        for b in rdprot.GetBonds():
            if b.GetBondType() == 2: print("double bond found!"); break
        loc = input_dict[self.input_keys[1]].content
        rdpoc = self.isolate(rdprot, loc)
        main_out_key = self.output_keys[0]

        rdpoc_pdb = generate_path_in_dir(6, self.base_dir, "_pocisoout.pdb")
        Chem.MolToPDBFile(rdpoc, rdpoc_pdb)
        
        return {main_out_key: self._get_out_rep(main_out_key)(rdpoc)}
        # Before: rdpoc_pdb

    def isolate(self, rdprot: Chem.Mol, loc: PocketLocation) -> Chem.Mol:
        remres, nonprotres = self.get_remaining_residues(rdprot, loc)
        terminal = self.identify_terminal_cappings(rdprot, remres)

        # rdprote = Chem.EditableMol(rdprot)
        rdprote = Chem.RWMol(rdprot)

        rdprote.BeginBatchEdit()
        for atom in rdprot.GetAtoms():
            atom: Chem.Atom
            cyx_fix(rdprote, atom, loc)
        rdprote.CommitBatchEdit()

        rdprote.BeginBatchEdit()
        for atom in rdprot.GetAtoms():
            atom: Chem.Atom
            rn = atom.GetPDBResidueInfo().GetResidueNumber()
            if rn in remres:
                if terminal[remres.index(rn)] == "N":
                    nterm(rdprote, atom)
                if terminal[remres.index(rn)] == "C":
                    cterm(rdprote, atom)
            elif rn not in nonprotres:
                rdprote.RemoveAtom(atom.GetIdx())
        rdprote.CommitBatchEdit()

        rdprotee = rdprote.GetMol()
        Chem.SanitizeMol(rdprotee)

        return rdprotee
    
    @staticmethod
    def get_remaining_residues(rdprot: Chem.Mol, loc: PocketLocation):
        remres = set()
        nonprotres = set()
        pos = rdprot.GetConformer().GetPositions()
        allres = set(atom.GetPDBResidueInfo().GetResidueNumber() 
                    for atom in rdprot.GetAtoms())
        for atom in rdprot.GetAtoms():
            atom: Chem.Atom
            # if atom.GetPDBResidueInfo().GetIsHeteroAtom(): # NOTE: What if a metal is hetero?
            #     continue

            if atom.GetPDBResidueInfo().GetResidueName() not in PROT_RESNAMES:
                if atom.GetPDBResidueInfo().GetResidueName() not in ["ACE", "NME"]:
                    nonprotres.add(atom.GetPDBResidueInfo().GetResidueNumber())
                continue

            if loc.point_is_included(pos[atom.GetIdx()]):
                rn = atom.GetPDBResidueInfo().GetResidueNumber()
                remres.add(rn)
                if rn-1 in allres:
                    remres.add(rn-1)
                if rn+1 in allres:
                    remres.add(rn+1)
        return sorted(remres), sorted(nonprotres)
    
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

            if prev_gap > 1 and next_gap > 1: # TODO: Can it ever happen?
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