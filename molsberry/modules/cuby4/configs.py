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

            if v.__class__ not in [bool, str, dict, OrderedDict]:
                v = str(v)

            if isinstance(v, str):
                cstring += f"{indent}{k}: {v}\n"

            elif isinstance(v, dict):
                cstring += f"{indent}{k}:\n"
                indent += "  "
                cstring += cls.config_to_string(v, indent)
                # print(cls.config_to_string(v, indent))
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
                 setpi: bool = False,
                 setcharges: bool = False,
                 cvb: bool = False,
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
        self.setpi = setpi
        self.setcharges = setcharges
        self.cvb = cvb

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

#####################################################################
##                          AMBER Classes                          ##
#####################################################################

class Cuby4AMBERInterfaceConfig(Cuby4InterfaceConfig, MOPACConfigUtils):
    def __init__(self,
                 home: str):
        
        super().__init__(interface="amber")
        self.home = home

        self.config["amber_amberhome"] = home

#####################################################################
##                          QMMM Classes                           ##
#####################################################################

class Cuby4QMMMInterfaceConfig(Cuby4InterfaceConfig):
    def __init__(self,
                 qm_config: Cuby4Config,
                 mm_config: Cuby4Config, 
                 embedding: str = "mechanical",
                 grad_on_point_charges: bool = False):
        
        super().__init__(interface="qmmm")
        
        self.qm_config = qm_config
        self.mm_config = mm_config
        
        if isinstance(grad_on_point_charges, str):
            grad_on_point_charges = {"no": False, 
                                     "yes": True}.get(grad_on_point_charges)
        self.grad_on_point_charges = grad_on_point_charges
        self.embedding = embedding

        self.config["calculation_qm"] = qm_config.config
        self.config["calculation_mm"] = mm_config.config
        self.config["qmmm_embedding"] = embedding
        self.config["gradient_on_point_charges"] = {True: "yes", 
                                                    False: "no"}.get(grad_on_point_charges)




