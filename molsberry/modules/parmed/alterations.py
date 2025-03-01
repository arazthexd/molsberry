from typing import Dict, List, Tuple

import numpy as np

from openmm import app, unit
from pdbfixer import PDBFixer
import parmed
from openff.toolkit import Molecule
from openff.units.openmm import to_openmm as openff_unit_to_openmm
from openmmforcefields.generators import SMIRNOFFTemplateGenerator

from ...core import SimpleBlock, MoleculeData, generate_path_in_dir
from ..openmm import DEFAULT_FORCEFIELDS
from .representations import ParmedMolRep
from ...core import PDBPathRep, LocationData
from ...core import PocketLocationRep, PocketLocation

def shorten_distance(pos_fixed: np.ndarray, pos_free: np.ndarray, dist: float, debug: bool = False):
    diff = pos_free - pos_fixed
    diff = diff / np.linalg.norm(diff)
    diff = diff * dist
    return pos_fixed + diff

def delete_atom(structure: parmed.Structure, atom: parmed.Atom, debug: bool = False):
    for bond in list(atom.bonds):
        if bond in structure.bonds:
            bond: parmed.Bond
            # bond.delete()
            structure.bonds.remove(bond)  
    
    if atom.residue is not None:
        if debug:
            print("removing atom", atom, "from res", atom.residue)
        atom.residue.atoms.remove(atom)
    
    if atom in list(structure.atoms): 
        if debug:
            print("removing atom", atom, "from struct", structure)
        structure.atoms.remove(atom)

def delete_residue_atoms(structure: parmed.Structure, residue: parmed.Residue, debug: bool = False):
    for atom in list(residue.atoms):
        atom: parmed.Atom
        delete_atom(structure, atom, debug=debug)

def clean_up(structure: parmed.Structure, debug: bool = False):
    for atom in list(structure.atoms):
        atom: parmed.Atom

        for bond in list(atom.bonds):
            bond: parmed.Bond

            if bond.atom1.residue is None or bond.atom2.residue is None:
                atom.bonds.remove(bond)

    res_to_del = [res for res in structure.residues if res.is_empty()]
    [structure.residues.remove(res) for res in res_to_del]
    
    for residue in structure.residues:
        if residue.chain == '':
            residue.chain = "A"
    
    # model = app.Modeller(structure.topology, structure.coordinates*unit.angstrom)
    # model.addHydrogens()
    # return parmed.openmm.load_topology(model.topology, xyz=model.positions)
    return structure

class ParmedProteinPocketIsolator(SimpleBlock):
    name = "parmed_isolator"
    display_name = "ParmEd Pocket Isolator"
    inputs = [
        ("protein", MoleculeData, ParmedMolRep, False),
        ("location", LocationData, PocketLocationRep, False)
    ]
    outputs = [
        ("pocket", MoleculeData, ParmedMolRep, False),
    ]
    batch_groups = [("protein", "location")]

    def __init__(self, debug: bool = False, **kwargs):
        super().__init__(**kwargs)
        self.debug = debug

    def operate(self, input_dict: Dict[str, ParmedMolRep | PocketLocationRep]) \
        -> Dict[str, ParmedMolRep]:
        prot = input_dict["protein"].content
        loc = input_dict["location"].content
        poc = self.isolate(prot, loc)
        poc.save(generate_path_in_dir(6, self.base_dir, "_pocisoout.pdb"))
        return {"pocket": ParmedMolRep(poc)}

    def fix_cyx(self, prot: parmed.Structure, loc: PocketLocation):
        cyx_residues = [res for res in prot.residues 
                        if res.name.strip() == "CYX"]
        for res in prot.residues:
            if res.name.strip() != "CYS":
                continue
            at_s = next(at for at in res.atoms if at.element_name == "S")
            has_nei_s = any(at.element_name=="S" for at in at_s.bond_partners)
            if not has_nei_s:
                continue
            res.name = "CYX"
            cyx_residues.append(res)
        
        for cyx_res in cyx_residues:
            cyx_res: parmed.Residue
            if not self._is_residue_included(cyx_res, loc):
                continue

            for cyx_at in cyx_res.atoms:
                cyx_at: parmed.Atom
                if cyx_at.element_name.strip() != "S":
                    continue

                nei_s = [nei for nei in cyx_at.bond_partners 
                         if nei.element_name.strip() == "S"]
                if nei_s.__len__() == 0:
                    if self.debug:
                        print("WARNING: Potential problem in fix_cys!")
                    break
                nei_s: parmed.Atom = nei_s[0]

                if self._is_residue_included(nei_s.residue, loc):
                    break

                # Change neighbor S to H
                nei_s_coords = np.array([nei_s.xx, nei_s.xy, nei_s.xz])
                cyx_s_coords = np.array([cyx_at.xx, cyx_at.xy, cyx_at.xz])

                new_nei_h_coords = shorten_distance(cyx_s_coords,
                                                    nei_s_coords,
                                                    dist=1.337)
                # nei_s.name = "HG"
                # nei_s.atomic_number = 1
                # nei_s.formal_charge = 0
                # delete_atom(prot, nei_s, self.debug)
                # prot.add_atom(nei_s, "CYX", cyx_res.number)
                # print(cyx_res.atoms)
                # nei_s.xx, nei_s.xy, nei_s.xz = new_nei_h_coords

                # non_s_bonds = [b for b in nei_s.bonds if cyx_at not in b]
                # [b.delete() for b in non_s_bonds]
                # [prot.bonds.remove(b) for b in non_s_bonds if b in prot.bonds]

                nei_h = parmed.Atom(atomic_number=1, formal_charge=0)
                nei_h.xx, nei_h.xy, nei_h.xz = new_nei_h_coords
                nei_h.name = "HG"
                prot.add_atom_to_residue(nei_h, cyx_res)
                hs_bond = parmed.Bond(nei_h, cyx_at)
                prot.bonds.append(hs_bond)
                
                # Change residue
                cyx_res.name = "CYS"
                break

    def _is_residue_included(self, 
                             residue: parmed.Residue, 
                             loc: PocketLocation) -> bool:
        for atom in residue.atoms:
            coords = np.array([atom.xx, atom.xy, atom.xz])
            if self.debug:
                print("checking inclusion", "coords:", coords)
            if loc.point_is_included(coords):
                return True
        return False

    def _get_neighboring_residues(self, 
                                  structure: parmed.Structure,
                                  residue: parmed.Residue):
        nei_n = None
        nei_c = None
        for atom in residue.atoms:
            atom: parmed.Atom
            
            if atom.name == "N":
                neis = [nei for nei in atom.bond_partners if nei.name == "C"]
                
                if len(neis) > 1:
                    raise ValueError("More than 1 neighbors for N atom")
                elif len(neis) == 1:
                    nei_n = neis[0].residue
                continue
            
            if atom.name == "C":
                neis = [nei for nei in atom.bond_partners if nei.name == "N"]
                
                if len(neis) > 1:
                    raise ValueError("More than 1 neighbors for N atom")
                elif len(neis) == 1:
                    nei_c = neis[0].residue
                continue
        
        return nei_n, nei_c

    def create_nterm(self, 
                     structure: parmed.Structure, 
                     residue: parmed.Residue):
        assert residue in structure.residues

        if residue.name == "ACE":
            return

        for atom in list(residue.atoms):
            atom: parmed.Atom

            if atom.name.strip() in ["O", "C"]:
                if self.debug:
                    print("atom", atom, "is kept")
                continue

            if atom.name.strip() == "CA":
                if self.debug:
                    print("atom", atom, "name is changed to CH3")
                atom.name = "CH3"
                continue

            if atom.name.strip() in ["HA", "HA2", "HA3"]:
                if self.debug:
                    print("atom", atom, "name is changed to HH31/2")
                atom.name = {"HA": "HH31", 
                             "HA2": "HH31", 
                             "HA3": "HH32"}.get(atom.name)
                continue

            if atom.name.strip() in ["CB", "N"]:
                if self.debug:
                    print("atom", atom, "will turn to H")
                atom.name = {"CB": "HH32", 
                             "N": "HH33"}.get(atom.name)
                atom.atomic_number = 1
                # TODO: anything else?

                if self.debug:
                    print("new atom:", atom)

                atom_ca = next(at for at in residue.atoms 
                               if at.name.strip() in ["CA", "CH3"])
                if self.debug:
                    print("CH3/CA nei", atom_ca, "found")

                non_ca_bonds = [b for b in atom.bonds if atom_ca not in b]
                [b.delete() for b in non_ca_bonds]
                [structure.bonds.remove(b) for b in non_ca_bonds if b in structure.bonds]

                atom.xx, atom.xy, atom.xz = shorten_distance(
                    pos_fixed=np.array([atom_ca.xx, atom_ca.xy, atom_ca.xz]),
                    pos_free=np.array([atom.xx, atom.xy, atom.xz]),
                    dist=1.08
                )
                continue
            
            if self.debug:
                print("deleting atom", atom)
            delete_atom(structure, atom, debug=self.debug)

        residue.name = "ACE"

    def create_cterm(self, 
                     structure: parmed.Structure, 
                     residue: parmed.Residue):
        assert residue in structure.residues

        if residue.name == "NME":
            return

        for atom in list(residue.atoms):
            atom: parmed.Atom

            if atom.name.strip() in ["N", "H"]:
                if self.debug:
                    print("atom", atom, "is kept")
                continue

            if atom.name.strip() == "CA":
                if self.debug:
                    print("atom", atom, "name is changed to CH3")
                atom.name = "CH3"
                continue

            if atom.name.strip() in ["HA", "HA2", "HA3"]:
                if self.debug:
                    print("atom", atom, "name is changed to HH31/2")
                atom.name = {"HA": "HH31", 
                             "HA2": "HH31", 
                             "HA3": "HH32"}.get(atom.name)
                continue

            if atom.name.strip() in ["CB", "C"]:
                if self.debug:
                    print("atom", atom, "will turn to H")
                atom.name = {"CB": "HH32", 
                             "C": "HH33"}.get(atom.name)
                atom.atomic_number = 1
                # TODO: anything else?

                if self.debug:
                    print("new atom:", atom)

                atom_ca = next(at for at in residue.atoms 
                               if at.name.strip() in ["CA", "CH3"])
                if self.debug:
                    print("CH3/CA nei", atom_ca, "found")
                
                non_ca_bonds = [b for b in atom.bonds if atom_ca not in b]
                [b.delete() for b in non_ca_bonds]
                [structure.bonds.remove(b) for b in non_ca_bonds if b in structure.bonds]
                
                atom.xx, atom.xy, atom.xz = shorten_distance(
                    pos_fixed=np.array([atom_ca.xx, atom_ca.xy, atom_ca.xz]),
                    pos_free=np.array([atom.xx, atom.xy, atom.xz]),
                    dist=1.08
                )
                continue
            
            if atom.name.strip() == "CD" and residue.name.strip() == "PRO":
                atom.name = "H"
                atom.atomic_number = 1

                atom_n = next(at for at in residue.atoms 
                              if at.name.strip() == "N")
                
                non_n_bonds = [b for b in atom.bonds if atom_n not in b]
                [b.delete() for b in non_n_bonds]
                [structure.bonds.remove(b) for b in non_n_bonds if b in structure.bonds]
                
                atom.xx, atom.xy, atom.xz = shorten_distance(
                    pos_fixed=np.array([atom_n.xx, atom_n.xy, atom_n.xz]),
                    pos_free=np.array([atom.xx, atom.xy, atom.xz]),
                    dist=1.08
                )
                continue
            
            if self.debug:
                print("deleting atom", atom)
            delete_atom(structure, atom, debug=self.debug)

        residue.name = "NME"

    def isolate(self, prot: parmed.Structure, loc: PocketLocation):
        self.fix_cyx(prot, loc)
        
        for res in prot.residues:
            res: parmed.Residue
            if self.debug:
                print()

            # If current residue is included, keep.
            if self._is_residue_included(res, loc):
                if self.debug:
                    print("keeping", res)
                continue
            
            nei_n, nei_c = self._get_neighboring_residues(prot, res)
            if nei_n is not None:
                is_inc_n = self._is_residue_included(nei_n, loc)
            else:
                is_inc_n = False
            if nei_c is not None:
                is_inc_c = self._is_residue_included(nei_c, loc)
            else:
                is_inc_c = False
            if self.debug:
                print("neis of", res,"are", nei_n, nei_c)
            
            # If both neighbors are included, keep it.
            if is_inc_n and is_inc_c:
                if self.debug:
                    print("both neis included, keeping", res)
                continue
            
            # If neighbors are not included as well as itself, delete
            if not (is_inc_n or is_inc_c):
                if self.debug:
                    print("no neis included, deleting")
                delete_residue_atoms(prot, res, debug=self.debug)
                continue

            if is_inc_n:
                if self.debug:
                    print("nei connecting to n,", nei_n, "is included, turning to NME")
                self.create_cterm(prot, res)
                if self.debug:
                    print("new res:", res)
            
            if is_inc_c:
                if self.debug:
                    print("nei connecting to c,", nei_c, "is included, turning to ACE")
                self.create_nterm(prot, res)
                if self.debug:
                    print("new res:", res)
                    
        return clean_up(prot, debug=self.debug)
    

class ParmedMoleculeCombiner(SimpleBlock):
    name = "parmedmoleculecombiner"

    outputs = [("combined", MoleculeData, ParmedMolRep, False)]

    def __init__(self, num_inputs: int, debug=False, 
                 save_output=False, num_workers=None):  
        super().__init__(debug, save_output, num_workers)  
        self.num_inputs = num_inputs  
        self._inputs = self._generate_inputs()  

    @property  
    def inputs(self) -> List[tuple]:  
        return self._inputs  
    
    @property
    def batch_groups(self) -> List[Tuple[str]]:
        return [tuple(f"molecule_{i+1}" for i in range(self.num_inputs))]

    def _generate_inputs(self) -> List[tuple]:  
        return [(f"molecule_{i+1}", MoleculeData, ParmedMolRep, False) 
                for i in range(self.num_inputs)] 
    
    def operate(self, input_dict: Dict[str, ParmedMolRep]):
        molecules = [input_dict[self.input_keys[i]].content 
            for i in range(self.num_inputs)]  
        merged_structures = sum(molecules[1:], start=molecules[0])
        merged_structures: parmed.Structure
        output = {"combined": ParmedMolRep(merged_structures)}
        return output
        