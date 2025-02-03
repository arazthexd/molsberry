from typing import Dict, List, Any, Tuple
from shutil import which
from collections import OrderedDict

from ...core import generate_path_in_dir

#####################################################################
##                          Base Classes                           ##
#####################################################################
class Cuby4Config:
    def __init__(self, config: Dict[str, Dict | Any] = None):
        if config is None:
            config = OrderedDict()
        self.config = config

    def get_string(self):
        return self.config_to_string(self.config, "")
    
    @classmethod
    def config_to_string(cls, config: Dict[str, Any], indent: str = ""):
        cstring = ""

        for k, v in config.items():
            if isinstance(v, bool):
                v = {True: "yes", False: "no"}.get(v)

            if v.__class__ not in [bool, str, dict]:
                v = str(v)

            if isinstance(v, str):
                cstring += f"{indent}{k}: {v}\n"

            elif isinstance(v, dict):
                cstring += f"{indent}{k}:\n"
                indent += "  "
                cstring += cls.config_to_string(v, indent)
                print(cls.config_to_string(v, indent))
                indent = indent[:-2]

        return cstring

class Cuby4InterfaceConfig(Cuby4Config):
    def __init__(self, 
                 interface: str = None):
        super().__init__()
        self.config["interface"] = interface

class Cuby4CompositeInterfaceConfig(Cuby4InterfaceConfig):
    pass

class Cuby4JobConfig(Cuby4Config):
    def __init__(self,
                 job: str = None):
        super().__init__()
        self.config["job"] = job

class Cuby4MergedConfig(Cuby4Config):

    @classmethod
    def from_config_list(cls, configs: List[Cuby4Config]):
        dicts = [c.config for c in configs]
        merged = OrderedDict()
        for d in dicts:
            merged.update(d)

        # TODO for mopac keywords
        mopac_kwds = ""
        for d in dicts:
            if "mopac_keywords" in d:
                mopac_kwds = " ".join([mopac_kwds, d["mopac_keywords"]])
        if len(mopac_kwds) > 0:
            merged["mopac_keywords"] = mopac_kwds

        return cls(merged)

#####################################################################
##                           Job Classes                           ##
#####################################################################

class Cuby4EnergyJobConfig(Cuby4JobConfig):
    def __init__(self,
                 geometry: str,
                 charge: str | int = 0,
                 mult: str | int = 1):
        
        super().__init__(job="energy")
        charge = int(charge)
        mult = int(mult)

        self.geometry = geometry
        self.charge = charge
        self.mult = mult

        self.config["geometry"] = geometry
        self.config["charge"] = charge
        self.config["multiplicity"] = mult

class Cuby4OptimizeJobConfig(Cuby4EnergyJobConfig):
    def __init__(self, 
                 geometry: str, 
                 charge: int = 0, 
                 mult: int = 1,
                 restart_file: str = None):
        super().__init__(geometry, charge, mult)
        self.job = "optimize"
        self.config["job"] = "optimize" # TODO
        self.restart_file = restart_file
        self.config["restart_file"] = restart_file

#####################################################################
##                          MOPAC Classes                          ##
#####################################################################

from ..mopac import MOPACInputMolRep

class MOPACConfigUtils:
    def _update_c4mopac_keywords(self, kws: List[str]):
        if hasattr(self, "keywords"):
            for kw in kws:
                if kw not in self.keywords:
                    self.keywords.append(kw)
        else:
            self.keywords = kws
        
        if not hasattr(self, "config"):
            self.config = {}
        self.config["mopac_keywords"] = " ".join(self.keywords)
    
    # def _mopac_to_confdict(self, mop_mol: MOPACInputMolRep) -> dict:
    #     conf_dict = {
    #         "geometry": None,
    #         "charge": None,
    #         "multiplicity": None
    #     }


class Cuby4MOPACInterfaceConfig(Cuby4InterfaceConfig, MOPACConfigUtils):
    def __init__(self,
                 exe: str = "auto",
                 method = "pm6",
                 mozyme: str | bool = False,
                 keywords: str | List[str] = []):
        
        super().__init__(interface="mopac")
        if exe == "auto":
            exe = which("mopac")
        if isinstance(mozyme, str):
            mozyme = {"yes": True, "no": False}.get(mozyme.lower())
        if isinstance(keywords, str):
            keywords = keywords.split()
            
        self.exe = exe
        self.method = method
        self.mozyme = mozyme

        self._update_c4mopac_keywords(keywords)

        self.config["mopac_exe"] = exe
        self.config["method"] = method
        self.config["mopac_mozyme"] = mozyme

class MOPACMozymeConfig(Cuby4Config, MOPACConfigUtils):
    def __init__(self,
                 setpi: List[Tuple[int, int]] = None, 
                 neg_cvb: List[Tuple[int, int]] = None):
        super().__init__()
        
        if setpi is None:
            setpi = []
        if neg_cvb is None:
            neg_cvb = []
        
        self.setpi = setpi
        self.neg_cvb = neg_cvb

        self.compile_setpi()
        self.compile_cvb()

    def compile_setpi(self):
        if len(self.setpi) == 0:
            self.config.pop("mopac_setpi", None)
            return
        
        setpi_str = ""
        for aid1, aid2 in self.setpi:
            setpi_str += f"\n  - {aid1};{aid2}"
        
        self.config["mopac_setpi"] = setpi_str
        return
    
    def compile_cvb(self):

        if len(self.neg_cvb) == 0:
            return
        
        cvb_str = "CVB("
        for aid1, aid2 in self.neg_cvb:
            if cvb_str != "CVB(":
                cvb_str += ";"
            cvb_str += f"{aid1}:-{aid2}"
        cvb_str += ")"

        self._update_c4mopac_keywords([cvb_str])

# class Cuby4MOPACJobConfig(Cuby4JobConfig, MOPACConfigUtils):
    
#     @classmethod
#     def from_mopac_mol(cls, mopac_mol: MOPACInputMolRep, mozyme: bool = False):
#         if mozyme:
#             conf = Cuby4MergedConfig.from_config_list([
#                 Cuby4()
#             ])

# from rdkit import Chem
# from rdkit.Chem import rdDistGeom, rdForceFieldHelpers
# import os

# from abc import ABC, abstractmethod
# import pathlib, shutil, subprocess, os, glob
# import string
# import random
# from typing import List, Tuple

# from rdkit import Chem 
# from rdkit.Geometry import Point3D

# from molsberry.core import generate_name_in_dir

# m = Chem.MolFromSmiles("CCCOC(=O)[O-]")
# m = Chem.AddHs(m)
# rdDistGeom.EmbedMolecule(m, useRandomCoords=True, randomSeed=-1)
# rdForceFieldHelpers.MMFFOptimizeMolecule(m)
# mp = os.path.join(os.path.dirname(__file__), "test_mol.pdb")
# Chem.MolToPDBFile(m, mp)

# c = Cuby4MergedConfig.from_config_list([
#     Cuby4EnergyJobConfig(mp, charge=-1, mult=1),
#     MOPACMozymeConfig(setpi=[], neg_cvb=[])
# ])

# class Cuby4Interface:
#     def __init__(self, interface_config: Cuby4InterfaceConfig,
#                  cuby4_exe: str = "auto", work_dir: str = "."):
#         if cuby4_exe == "auto":
#             self.path = shutil.which("cuby4")
#         else:
#             self.path = str(pathlib.Path(cuby4_exe).absolute())
#         self.work_dir = str(pathlib.Path(work_dir).absolute())
#         self.interface_config = interface_config
    
#     def run(self, config: Cuby4Config, job_name: str = None):
        
#         if job_name is None:
#             job_name = generate_name_in_dir(6, self.work_dir, ".yml")
#         config_path = os.path.join(self.work_dir, f"{job_name}_conf.yml")
#         config_path = str(pathlib.Path(config_path).absolute())
#         config_final = Cuby4MergedConfig.from_config_list([self.interface_config,
#                                                            config])
#         config_str = config_final.get_string()
#         with open(config_path, "w") as f:
#             f.write(config_str)
        
#         prev_path = os.curdir
#         os.chdir(self.work_dir)

#         for filename in glob.glob("job_*"):
#             shutil.rmtree(filename) 
#         output = subprocess.run([self.path, config_path], capture_output=True)
        
#         os.chdir(prev_path)

#         if output.returncode != 0:
#             print(output)
#             raise ValueError()
#         return output.stdout.decode()
    
# interface = Cuby4Interface(
#     Cuby4MOPACInterfaceConfig(mozyme=True), work_dir=os.path.dirname(__file__)
# )
# print(interface.run(c))

# # with open(os.path.join(os.path.dirname(__file__), "c4c.yml"), "w") as f:
# #     f.write(c.get_string())