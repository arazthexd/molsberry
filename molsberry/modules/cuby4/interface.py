from abc import ABC, abstractmethod
import pathlib, shutil, subprocess, os, glob
import string
import random
from typing import List, Tuple

from rdkit import Chem 
from rdkit.Geometry import Point3D

from ...core import generate_name_in_dir
from .configs import Cuby4InterfaceConfig, Cuby4Config, Cuby4MergedConfig


class Cuby4Interface:
    def __init__(self, interface_config: Cuby4InterfaceConfig,
                 cuby4_exe: str = "auto", work_dir: str = "."):
        if cuby4_exe == "auto":
            self.path = shutil.which("cuby4")
        else:
            self.path = str(pathlib.Path(cuby4_exe).absolute())
        self.work_dir = str(pathlib.Path(work_dir).absolute())
        self.interface_config = interface_config
    
    def run(self, config: Cuby4Config, job_name: str = None):
        
        if job_name is None:
            job_name = generate_name_in_dir(6, self.work_dir, ".yml")
        config_path = os.path.join(self.work_dir, f"{job_name}_conf.yml")
        config_path = str(pathlib.Path(config_path).absolute())
        config_final = Cuby4MergedConfig.from_config_list([self.interface_config,
                                                           config])
        config_str = config_final.get_string()
        with open(config_path, "w") as f:
            f.write(config_str)
        
        prev_path = os.curdir
        os.chdir(self.work_dir)

        for filename in glob.glob("job_*"):
            shutil.rmtree(filename) 
        output = subprocess.run([self.path, config_path], capture_output=True)
        
        os.chdir(prev_path)

        if output.returncode != 0:
            print(output)
            raise ValueError()
        return output.stdout.decode()
