import os, pathlib, subprocess, shutil, warnings
from typing import List, Tuple

import pandas as pd
from rdkit import Chem
from biopandas.pdb import PandasPdb

from ..abstract import *
from ..utils.pockets import PocketLocation
from ..utils.iotools import write_ligands, load_ligands, get_unique_full_name

class GninaInterface(Interface):

    def __init__(self, bin_path: str = "auto", work_dir: str = "tmp") -> None:
        
        if bin_path == "auto":
            try:
                bin_path = shutil.which("gnina")
            except:
                raise FileNotFoundError()
            
        if not os.path.exists(bin_path):
            raise FileNotFoundError()
        
        if not os.path.exists(work_dir):
            os.mkdir(work_dir)
        
        self.bin_path: str = str(pathlib.Path(bin_path).absolute())
        self.work_dir: str = str(pathlib.Path(work_dir).absolute())

        self.exhaustiveness: int = 8
        self.num_modes: int = 9
        self.cpu: int = 8

        self.other_params: list = []
    
    def run(self, cmd: List[str], debug: bool = False):
        output = subprocess.run(cmd, capture_output=not debug)
        if output.returncode != 0:
            print(output.stderr.decode())
            print(output.stdout.decode())
            raise ChildProcessError()

    def set_options(self, exhaustiveness: int = 8, num_modes: int = 9, cpu: int = 1, others: str = "") -> None:
        self.exhaustiveness = exhaustiveness
        self.num_modes = num_modes
        self.cpu = cpu

        self.other_params = others.split()
        if self.other_params is None:
            self.other_params = []
    
    def update_options(self, **kwargs): # TODO: testing...
        for arg, value in kwargs.items():
            if arg in self.__dict__.keys():
                setattr(self, arg, value)
            else:
                warnings.warn("One of the attributes does not exist...")
    
    def _get_params_cmd(self) -> List[str]:
        params = [
            "--exhaustiveness", str(self.exhaustiveness),
            "--num_modes", str(self.num_modes),
            "--cpu", str(self.cpu),
        ] + self.other_params
        return params
    
    @staticmethod
    def _pocloc_to_cmdline(pocket_loc: PocketLocation) -> List[str]:
        box = pocket_loc.get_params("boxcxyz")
        cmdline = ["--center_x", str(box['c'][0]), "--center_y", str(box['c'][1]), "--center_z", str(box['c'][2]), 
                   "--size_x", str(box['size'][0]), "--size_y", str(box['size'][1]), "--size_z", str(box['size'][2])]
        return cmdline
    
    def _get_flexible_models(self, flexibles: str):
        flexibles_prefix = os.path.join(self.work_dir, "gnina_flexible")
        with open(flexibles, "r") as f:
            flexibles_str = f.read()
        flexibles_models = flexibles_str.split("\nENDMDL")
        out_names = []
        for i, model in enumerate(flexibles_models):
            out_file = get_unique_full_name(f"{flexibles_prefix}_{i}", ".pdb", n=5)
            with open(out_file, "w") as f:
                flexibles_str = f.write(model)
            out_names.append(out_file)
        return out_names
    
    # def update_target_with_flexibles(self, target: str, flexibles: str) -> List[str]: # TODO: Needs fixing... DO NOT USE
    #     target_ppdb = PandasPdb().read_pdb(target)
    #     target_df: pd.DataFrame = target_ppdb.df["ATOM"]

    #     flexibles_ppdb = PandasPdb().read_pdb(flexibles)
    #     flexibles_ppdb.label_models()
    #     flex_df: pd.DataFrame = flexibles_ppdb.df["ATOM"]
    #     flex_model_ids = list(flex_df["model_id"].unique())
    #     flex_reses = flex_df[["chain_id", "residue_number"]].drop_duplicates()

    #     target_df.drop(target_df[(target_df["residue_number"].isin(list(flex_reses["residue_number"])))
    #                              & (target_df["chain_id"].isin(list(flex_reses["chain_id"])))].index, inplace=True, axis=0)
        
    #     new_targets = []
    #     for model_id in flex_model_ids:
    #         model_ppdb = flexibles_ppdb.get_model(model_id)
    #         model_df = model_ppdb.df["ATOM"]
    #         target_copy = target_df.copy()
    #         target_copy = pd.concat([target_copy, model_df]).sort_values(["chain_id", "residue_number"])
    #         target_copy.reset_index()
    #         target_copy["line_idx"] = pd.Series(list())
    #         target_ppdb.df["ATOM"] = target_copy
    #         out_path = f"{target[:-4]}_{model_id}.pdb"
    #         target_ppdb.to_pdb(f"{target[:-4]}_{model_id}.pdb")
    #         new_targets.append(out_path)
        
    #     return new_targets

class GninaDocker(ComplexEnumerator):
    name = "Gnina Docking"
    def __init__(self, interface: GninaInterface, pocket_loc: PocketLocation, batch_size: int = None,
                 out_path: str = "auto", in_path_ligand: str = "auto", flexible_residues: str = None, 
                 out_path_flexible: str = "auto", debug: bool = False) -> None:
        super().__init__(debug)
        self.interface = interface
        self.batch_size = batch_size # TODO: Implement batch sizes
        self.pocket_loc = pocket_loc
        if out_path == "auto":
            self.out_path = get_unique_full_name(
                os.path.join(self.interface.work_dir, "gnina_results"), ".sdf", 5
            )
        else:
            self.out_path = out_path
        if in_path_ligand == "auto":
            self.in_path_ligand = get_unique_full_name(
                os.path.join(self.interface.work_dir, "gnina_in_lig"), ".sdf", 5
            )
        else:
            self.in_path_ligand = in_path_ligand

        if flexible_residues:
            if out_path_flexible == "auto":
                out_path_flexible = os.path.join(self.interface.work_dir, "gnina_flex_out.pdb")
            self.interface.other_params.extend(["--flexres", flexible_residues, "--out_flex", out_path_flexible,
                                                "--full_flex_output"])

    def enumerate(self, ligand: Chem.Mol, target_path: str) -> Tuple[List[Chem.Mol], List[str]]:
        pocket_cmd = self.interface._pocloc_to_cmdline(self.pocket_loc)
        params_cmd = self.interface._get_params_cmd()
        write_ligands([ligand], self.in_path_ligand)
        interface_input = [
            self.interface.bin_path, "-l", self.in_path_ligand, "-r", target_path, "-o", self.out_path, 
        ] + pocket_cmd + params_cmd
        self.interface.run(interface_input)
        
        out_ligands = load_ligands(self.out_path, addHs=True)
        if "--flexres" not in self.interface._get_params_cmd():
            return out_ligands, [target_path] * len(out_ligands)
        out_flex = self.interface.other_params[
            self.interface.other_params.index("--out_flex") + 1
        ]
        return out_ligands, self.interface._get_flexible_models(out_flex)
    
    def dock(self, ligand: Chem.Mol, target_path: str) -> Tuple[List[Chem.Mol], List[str]]:
        return self.enumerate(ligand, target_path)

class GninaResultsSelector(LigandSelector):
    name = "Gnina Results Selection"
    def __init__(self, sort_by: str = "CNN_VS", max_vina: float = -4.0, min_posescore: float = 0.4,
                 select_ratio: float = 0.1, max_per_lig: int = 2, n_select: int = 1, debug: bool = False, identifier: str = "_Name"):
        super().__init__(identifier=identifier, debug=debug)
        self.sort_by = sort_by
        self.max_vina = max_vina
        self.min_posescore = min_posescore
        self.select_ratio = select_ratio
        self.max_per_lig = max_per_lig
        self.n_select = n_select
    
    def select(self, ligands: List[Chem.Mol]) -> List[Chem.Mol]:
        selected_ligs = []
        selected_names = dict()
        for mol in sorted(ligands, key=lambda lig: float(lig.GetProp(self.sort_by)), reverse=True):
            if len(selected_ligs) >= self.n_select:
                break
            if float(mol.GetProp("minimizedAffinity")) > self.max_vina:
                continue
            if float(mol.GetProp("CNNscore")) < self.min_posescore:
                continue
            name = mol.GetProp("_Name")
            if selected_names.get(name) is not None:
                if selected_names.get(name) >= self.max_per_lig:
                    continue
            else:
                selected_names[name] = 0
            selected_ligs.append(mol)
            selected_names[name] += 1
        return selected_ligs


