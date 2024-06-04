import os
from copy import copy, deepcopy
from abc import ABC, abstractmethod
from typing import Callable, Dict

import numpy as np 
import pandas as pd

from rdkit import Chem
from biopandas.pdb import PandasPdb

from .protio import read_prot_file

class PocketLocation: # TODO: Write docstrings for PocketLocation methods...
    """Handles pocket location specifications for different interfaces and converts \
        between them.

    This class can convert between different types of pocket location formats:
    - point:    Point and Radius, 
    - box2p:    3D Box With Two Opposite Corner Points,
    - boxcxyz:  3D Box Given Center Point and Sizes in X, Y, Z Dimensions,
    - ligand:   Ligand and Radius,
    - residues: Protein Residues
    """

    def __init__(self, method, **kwargs): 
        
        self.params = {}
        if method == "point":
            self._init_point(**kwargs)
        elif method == "box2p":
            self._init_box2p(**kwargs)
        elif method == "boxcxyz":
            self._init_boxcxyz(**kwargs)
        elif method == "ligand":
            self._init_ligand(**kwargs)
        elif method == "residues":
            self._init_residues(**kwargs)
        else:
            raise ValueError("The method should be in ('point', 'box', 'ligand', 'residues')")
        self.loc_type = method
    
    def point_is_included(self, point):
        if self.loc_type == "point":
            return self._point_is_included_point(point)
        elif self.loc_type == "ligand":
            return self._point_is_included_ligand(point)
        else:
            raise NotImplementedError() # TODO: Write other types...
    
    def get_params(self, param_type: str):
        if param_type not in self.params.keys():
            for format in self.params.keys():
                try:
                    self.loctype_conversion(format, param_type)
                    return self.params[param_type]
                except NotImplementedError:
                    continue
            raise ValueError(f"Could not convert existing formats to {param_type}")
        else:
            return self.params[param_type]
        
    def loctype_conversion(self, source_format, target_format):
        """
        Converts the pocket location to the target format.

        Args:
            target_format (str): The target format ('point', 'box', or 'ligand').
        """
        if source_format not in self.params.keys():
            raise ValueError("Source format is not defined!")
        if source_format == "ligand" and target_format == "boxcxyz":
            coords = self.params["ligand"]["c"]
            radius = self.params["ligand"]["r"]
            lig_center = (np.max(coords, axis=0) + np.min(coords, axis=0)) / 2
            size = np.max(coords, axis=0) - np.min(coords, axis=0) + radius
            self.params["boxcxyz"] = {"c": lig_center, "size": size}
        elif source_format == "point" and target_format == "boxcxyz":
            point = self.params["point"]["p"]
            radius = self.params["point"]["r"]
            self.params["boxcxyz"] = {"c": point, "size": np.array([radius]*3)}
        else:
            raise NotImplementedError()
        # TODO: Implement other conversions between pocket location formats

    def _init_point(self, point=None, radius=10):
        point = np.array(point)
        self.params["point"] = {"p": point, "r": radius}

    def _init_box2p(self, point1, point2):
        pass

    def _init_boxcxyz(self, center, xyz_size):
        self.params["boxcxyz"] = {"c": np.array(center), "size": np.array(xyz_size)}

    def _init_ligand(self, ligand, radius=10):
        if isinstance(ligand, str):
            if ligand.endswith(".pdb"):
                ligand = Chem.MolFromPDBFile(ligand)
            elif ligand.endswith(".mol2"):
                ligand = Chem.MolFromMol2File(ligand)
            elif ligand.endswith(".sdf"):
                ligand = next(Chem.SDMolSupplier(ligand))
            else:
                raise ValueError("ligand has to be .pdb, .mol2 or .sdf...")
        elif isinstance(ligand, Chem.Mol):
            pass
        else:
            raise ValueError("ligand has to be .pdb, .mol2 or .sdf or an rdmol instance...")
        
        coords = ligand.GetConformer().GetPositions()
        self.params["ligand"] = {
            "c": coords, "r": radius, "l": ligand
        }

    def _init_residues(self, protein_pdb, res_ids):
        pass

    def _point_is_included_point(self, point):
        point = np.array(point)
        if np.linalg.norm(point - self.params["point"]["p"]) <= self.params["point"]["r"]:
            return True
        return False
    
    def _point_is_included_ligand(self, point):
        point = np.array(point)
        diff = self.params["ligand"]["c"] - point
        if np.min(np.linalg.norm(diff, axis=1)) <= self.params["ligand"]["r"]:
            return True
        return False

class TerminalCapping(ABC):

    # TODO: Avoid too many warnings... manage it.

    def __init__(self, name: str = "CAP", addH_final: bool = True):
        self.addH_final = addH_final
        self.work_dir = "."
        self.cleanup = True
        self.capname = name
    
    def res_to_cap(self, residue: PandasPdb, ref: PandasPdb) -> PandasPdb:

        assert len(residue.df["ATOM"]["residue_number"].unique()) == 1
        assert len(ref.df["ATOM"]["residue_number"].unique()) == 1
        resid = residue.df["ATOM"].iloc[0]["residue_number"]

        tmp_path_init = os.path.join(self.work_dir, f"resref{resid}.pdb")
        tmp_path_final = os.path.join(self.work_dir, f"resref{resid}_f.pdb")

        resref = PandasPdb()
        resref.df["ATOM"] = pd.concat([residue.df["ATOM"], ref.df["ATOM"]], axis=0)
        resref.df["ATOM"].sort_index(ascending=True)
        resref.to_pdb(tmp_path_init)
        rdresref = Chem.MolFromPDBFile(tmp_path_init, removeHs=False)
        rdresref_edit = Chem.RWMol(rdresref)

        if self.cleanup:
            os.remove(tmp_path_init)

        rdresref_edit.BeginBatchEdit()

        for atom in rdresref_edit.GetAtoms():
            if int(atom.GetPDBResidueInfo().GetResidueNumber()) != int(resid):
                continue

            self.modify(rdresref_edit, atom)

        rdresref_edit.CommitBatchEdit()
        rdresref = rdresref_edit.GetMol()
        Chem.SanitizeMol(rdresref)

        if self.addH_final:
            rdresref: Chem.Mol = Chem.AddHs(rdresref, addCoords=True, addResidueInfo=False) # I DON'T KNOW WHAT'S WRONG
            hname_counter = 1
        
        for atom in rdresref.GetAtoms():
            atom: Chem.Atom
            
            info = atom.GetPDBResidueInfo()
            if info:
                if info.GetResidueNumber() == resid:
                    info.SetResidueName(self.capname)
                    info.SetIsHeteroAtom(False)
                continue

            if not atom.GetAtomicNum() == 1:
                continue
            
            nei: Chem.Atom = atom.GetNeighbors()[0]
            nei_info = nei.GetPDBResidueInfo()
            
            pdb_info = Chem.AtomPDBResidueInfo()
            hname = f"H{hname_counter}"
            if nei_info.GetResidueNumber() == resid:
                hname_counter += 1
                pdb_info.SetResidueName(self.capname)
            else:
                pdb_info.SetResidueName(nei_info.GetResidueName())
            pdb_info.SetIsHeteroAtom(False)
            pdb_info.SetChainId(nei_info.GetChainId())
            pdb_info.SetResidueNumber(int(nei_info.GetResidueNumber()))
            pdb_info.SetName(hname.center(4))
            atom.SetPDBResidueInfo(pdb_info)

        Chem.MolToPDBFile(rdresref, tmp_path_final)
        cap = PandasPdb().read_pdb(tmp_path_final)
        cap.df["ATOM"] = cap.df["ATOM"][cap.df["ATOM"]["residue_number"] == resid]
        if self.cleanup:
            os.remove(tmp_path_final)

        # TODO: Check that ref and res are not disconnected because of the modifications
        return cap

    @abstractmethod
    def modify(elf, edit_mol: Chem.RWMol, atom: Chem.Atom) -> None:
        pass

class GenericCapping(TerminalCapping):
        
    def __init__(self, name: str, keeping: list, modifications: dict, addH_final: bool = True, 
                 work_dir: str = "."):
        super().__init__(name, addH_final)
        self.keeping = keeping
        self.modifications = modifications
        self.work_dir = work_dir

        assert all([k in self.keeping for k in self.modifications.keys()])
    
    def modify(self, edit_mol: Chem.RWMol, atom: Chem.Atom) -> None:

        atom_name = atom.GetPDBResidueInfo().GetName().strip()
        atom_idx = atom.GetIdx()

        if not atom_name in self.keeping:
            edit_mol.RemoveAtom(atom_idx)
            return
        
        if not atom_name in self.modifications.keys():
            return
        
        mods = self.modifications[atom_name]

        # Cap Name
        pdb_info = atom.GetPDBResidueInfo()
        pdb_info.SetResidueName(self.capname)
        atom.SetPDBResidueInfo(pdb_info) # TODO: Why is it not taking place?

        # Modify
        for mod in mods.items():
            continue # TODO: Implement modifications if needed...
        
        return

# TODO: Delete all hydrogens and add them back at the end...
  
DEFAULT_NTER_CAP = GenericCapping(
    name="ACE",
    keeping=["C", "O", "CA"],
    modifications=dict(),
    addH_final=False
)

DEFAULT_CTER_CAP = GenericCapping(
    name="NME",
    keeping=["N", "CA"],
    modifications=dict(),
    addH_final=False
)

def isolate_pocket(protein_input, pocket_location: PocketLocation, output_path = None, 
                   capping: Dict[str, TerminalCapping] = {"nter": DEFAULT_NTER_CAP, "cter": DEFAULT_CTER_CAP}, 
                   work_dir="tmp"):
    
    if not os.path.exists(work_dir):
        os.mkdir(work_dir)

    # Load the protein file as a PandasPdb object and get all the residues numbers, 
    # atom coords and see which atoms are included in the specified pocket location.
    ppdb = read_prot_file(protein_input, "ppdb")
    resnum_set = set(ppdb.df['ATOM']['residue_number'])
    atom_coords = ppdb.df['ATOM'][['x_coord', 'y_coord', 'z_coord']].values
    ppdb.df["ATOM"]["included"] = [
        pocket_location.point_is_included(coords) for coords in atom_coords
    ]

    resid = -1
    atoms_df: pd.DataFrame = deepcopy(ppdb.df["ATOM"])

    # Iterate over atoms not included in pocket and if they are next to a border 
    # residue (included in the pocket location), make them into a capping.
    for _, row in ppdb.df["ATOM"][ppdb.df["ATOM"]["included"] == False].iterrows():
        
        # If the residue has been processed before (one of the atoms processed), skip...
        if resid == row["residue_number"]:
            continue
        resid = row["residue_number"]

        # If any atoms of current residue are in the pocket location, the residue is
        # included and will not be converted to a capping.
        if ppdb.df["ATOM"][ppdb.df["ATOM"]["residue_number"] == resid]["included"].sum() > 0:
            continue

        # Check whether residue should be converted into an N- or C-terminal
        # If a residue has a larger ID, it is closer to the C-terminal and should be 
        # considered a C-terminal capping and vice versa
        upid, downid = resid + 1, resid - 1
        cap_as_nter = upid in resnum_set and ppdb.df["ATOM"][
            ppdb.df["ATOM"]["residue_number"] == upid
        ]["included"].sum() > 0
        cap_as_cter = downid in resnum_set and ppdb.df["ATOM"][
            ppdb.df["ATOM"]["residue_number"] == downid
        ]["included"].sum() > 0

        # In case a residue has to be converted to both C- and N-TER cappings, just keep it.
        if cap_as_cter and cap_as_nter:
            continue
        
        # In case a residue is not in the desired pocket and not converting to capping, it will
        # be removed from dataframe.
        if not (cap_as_cter or cap_as_nter):
            atoms_df = atoms_df[atoms_df["residue_number"] != resid]
            continue
        
        residue = PandasPdb()
        residue.df["ATOM"] = ppdb.df["ATOM"][ppdb.df["ATOM"]["residue_number"] == resid]
        ref = PandasPdb()
        if cap_as_cter:
            ref.df["ATOM"] = ppdb.df["ATOM"][ppdb.df["ATOM"]["residue_number"] == upid]
            cap: PandasPdb = capping["cter"].res_to_cap(residue, ref)
        if cap_as_nter:
            ref.df["ATOM"] = ppdb.df["ATOM"][ppdb.df["ATOM"]["residue_number"] == downid]
            cap: PandasPdb = capping["nter"].res_to_cap(residue, ref)
        
        res_start = atoms_df.index.get_loc(atoms_df[atoms_df["residue_number"] == resid].index.min())
        res_end = atoms_df.index.get_loc(atoms_df[atoms_df["residue_number"] == resid].index.max())
        atoms_df = pd.concat([atoms_df.iloc[:res_start], cap.df["ATOM"], atoms_df.iloc[res_end+1:]])

    ppdb.df["ATOM"] = atoms_df.sort_values(["residue_number", "atom_number"], ascending=True)
    ppdb.df["ATOM"]["atom_number"] = [i+1 for i in range(len(ppdb.df["ATOM"]))]
    ppdb.df["ATOM"]["line_idx"] = [i+1 for i in range(len(ppdb.df["ATOM"]))] # TODO: Clean this mess...

    if output_path:
        ppdb.to_pdb(output_path, records=["ATOM"])
    
    return ppdb


###########################
##         TESTS         ##
###########################

######## Capping ##########
# CLASS
# print("don't forget to del prints")
# pdb = PandasPdb().read_pdb("protein_clean.pdb")
# res, ref = PandasPdb(), PandasPdb()
# res.df["ATOM"] = pdb.df["ATOM"][pdb.df["ATOM"]["residue_number"] == 21]
# ref.df["ATOM"] = pdb.df["ATOM"][pdb.df["ATOM"]["residue_number"] == 20]
# c = DEFAULT_CTER_CAP.res_to_cap(res, ref)
# print(c.df["ATOM"])

# ISOLATE
# isolate_pocket("protein_clean.pdb", PocketLocation("point", point=(3, 11, 11), radius=20),
#                output_path="pocket.pdb")


###########################
##        DUMPED         ##
###########################

# DEFAULT_NTER = {
#     "keep_atoms": ["C", "O", "CA"],
#     "change_atoms": {
#         "C": "CCP",
#         "O": "OCP",
#         "CA": "HCP"
#     },
#     "atoms_to_hydrogen": ["CA"]
# }

# DEFAULT_NTER = {
#     "capname": "CPN",
#     "modify": {
#         "C": ("C", "C", "")
#     }
# }

# DEFAULT_CTER = {
#     "keep_atoms": ["N", "CA", "C", "CB", "HA", "HA2", "HA3", "H"],
#     "change_atoms": {
#         "N": "NCP",
#         "CA": "CACP",
#         "C": "HCP1",
#         "H": "HCPN",
#         "HA": "HCP2",
#         "CB": "HCP3",
#         "HA2": "HCP2",
#         "HA3": "HCP3"
#     },
#     "atoms_to_hydrogen": ["CB", "C"]
# }