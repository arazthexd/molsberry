from abc import ABC, abstractmethod
import pathlib, shutil, subprocess, os, glob
import string
import random
from typing import List, Tuple

from rdkit import Chem 
from rdkit.Geometry import Point3D

from ..abstract import *
from ..utils.iotools import write_ligands, save_pl_complex, get_unique_full_name
from .cuby4c import (
    Cuby4FullConfig, Cuby4EnergyConfig, Cuby4MOPACConfig, Cuby4OptimizeConfig,
    Cuby4InterfaceConfig, Cuby4JobConfig, Cuby4InteractionConfig
)

class Cuby4Interface(Interface):
    def __init__(self, intconfig: Cuby4InterfaceConfig,
                 cuby4_path: str = "auto", work_dir: str = "."):
        if cuby4_path == "auto":
            self.path = str(pathlib.Path(shutil.which("cuby4")).absolute())
        else:
            self.path = str(pathlib.Path(cuby4_path).absolute())
        self.work_dir = str(pathlib.Path(work_dir).absolute())
        self.intconfig = intconfig
    
    def run(self, config: Cuby4FullConfig, job_name: str = "job"):
        for filename in glob.glob("job_*"):
            shutil.rmtree(filename) 
        config_path = os.path.join(self.work_dir, f"cuby_{job_name}_config.yml")
        config_path = str(pathlib.Path(config_path).absolute())
        config_str = config.get_config_str()
        with open(config_path, "w") as f:
            f.write(config_str)
        output = subprocess.run([self.path, config_path], capture_output=True, check=True)
        return output.stdout.decode()

class Cuby4General:
    def __init__(self, interface: Cuby4Interface, n_threads: int = 1, work_dir: str = "."):
        self.interface = interface
        self.work_dir = str(pathlib.Path(work_dir).absolute())
        self.n_threads = n_threads
    
    def _energy(self, molecules: str | Chem.Mol | Tuple[Chem.Mol | str, str], job_name: str = None, 
                small_mols: List[Chem.Mol] = list()):
        self.jobconfig = Cuby4EnergyConfig(None, self.n_threads)
        self.jobconfig.update_config_for(molecules, unique_mols=small_mols)
        self.interface.intconfig.update_config_for(molecules, unique_mols=small_mols)
        if not job_name:
            job_name = self.jobconfig.geometry.split(".")[0].split("/")[-1]
        output = self.interface.run(Cuby4FullConfig(
            self.jobconfig, self.interface.intconfig), job_name)
        return float(output.split("nergy: ")[-1][:7])
    
    def _optimize(self, molecules: str | Chem.Mol | Tuple[Chem.Mol | str, str], max_cycles: int = 200, 
                  job_name: str = None, small_mols: List[Chem.Mol] = list()):
        restart_path = str(pathlib.Path(os.path.join(self.work_dir, f"{job_name}_restart.pdb")))
        self.jobconfig = Cuby4OptimizeConfig(None, restart_path, max_cycles, self.n_threads)
        self.jobconfig.update_config_for(molecules)
        self.interface.intconfig.update_config_for(molecules, unique_mols=small_mols)
        if not job_name:
            job_name = self.jobconfig.geometry.split(".")[0].split("/")[-1]
        output = self.interface.run(Cuby4FullConfig(
            self.jobconfig, self.interface.intconfig), job_name)
        return float(output.split("nergy: ")[-1][:7]), restart_path

class Cuby4LigandEnergyScorer(Cuby4General, LigandScorer):
    def __init__(self, interface: Cuby4Interface, n_threads: int = 1, work_dir: str = "."):
        super().__init__(interface, n_threads, work_dir)
    
    def score(self, ligand: Chem.Mol):
        return self._energy(ligand)

class Cuby4ComplexEnergyScorer(Cuby4General, ComplexScorer):
    def __init__(self, interface: Cuby4Interface, n_threads: int = 1, work_dir: str = "."):
        super().__init__(interface, n_threads, work_dir)
    
    def score(self, ligand: Chem.Mol, target: str):
        return self._energy((ligand, target), small_mols=[ligand])
    
class Cuby4ComplexInteractScorer(Cuby4General, ComplexScorer):
    def __init__(self, interface: Cuby4Interface, n_threads: int = 1, work_dir: str = "."):
        super().__init__(interface, n_threads, work_dir)
    
    def score(self, ligand: Chem.Mol, target: str):
        # comp_path = os.path.join(self.work_dir, "tmp_complex.pdb")
        # save_pl_complex(ligand, target, comp_path) # TODO: better implementation of job name
        lig_energy = self._energy(ligand)
        targ_energy = self._energy(target)
        comp_energy = self._energy((ligand, target), small_mols=[ligand])
        return comp_energy - lig_energy - targ_energy
    
class Cuby4ComplexOptimizer(Cuby4General, ComplexConverter):
    name = "Cuby4 Complex Optimizer"
    def __init__(self, interface: Cuby4Interface, n_threads: int = 1, work_dir: str = "."):
        super().__init__(interface, n_threads, work_dir)
    
    def convert(self, ligand: Chem.Mol, target: str) -> Tuple[Chem.Mol, str]:
        energy, out_path = self._optimize((ligand, target), job_name="cuby4_comp_optimize_target")
        hetlines: List[str] = []
        atlines: List[str] = []
        with open(out_path, "r") as f:
            for line in f.read().split("\n"):
                if "HETATM" in line:
                    hetlines.append(line)
                if "ATOM" in line:
                    atlines.append(line)

        n_atoms_lig = ligand.GetNumAtoms()
        conformer = ligand.GetConformer()
        for i, line in enumerate(hetlines[-1::-1][:n_atoms_lig]):
            x, y, z = line[30:54].split()
            conformer.SetAtomPosition(n_atoms_lig-i-1, Point3D(float(x), float(y), float(z)))
        
        target_prefix = os.path.join(self.work_dir, f"cuby4_comp_optimize_target")
        target_path = get_unique_full_name(target_prefix, ".pdb", n=5)
        with open(target_path, "w") as f:
            f.writelines(atlines)

        return ligand, target_path
