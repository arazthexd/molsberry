from __future__ import annotations

import os, pathlib
from abc import ABC, abstractmethod
from typing import List, Tuple, Any, Dict
from itertools import repeat
from collections import OrderedDict
from tqdm import tqdm

from paradag import DAG

from rdkit import Chem

from .utils.iotools import write_ligands

def zip_complex(ligands, targets):
    if len(targets) == 1:
        zipped = zip(ligands, repeat(targets[0]))
    elif len(targets) == len(ligands):
        zipped = zip(ligands, targets)
    else:
        raise ValueError()
    return zipped

class Interface(ABC):
    def __init__(self):
        pass

class PipelineBlock(ABC):
    def __init__(self, debug: bool = False):
        self.debug = debug
        self.executed = False
        self.ligands_out = list()
        self.targets_out = list()
        self.extra_info = dict()
    
    def reset(self):
        self.executed = False
        self.ligands_out = list()
        self.targets_out = list()
        self.extra_info = dict()

    @abstractmethod
    def run(self, ligands: List[Chem.Mol], targets: List[str], 
            extra_info: dict = dict()) -> Tuple[List[Chem.Mol], List[str]]:
        if self.debug: print("target len", len(targets))
        assert len(ligands) == len(targets) or len(targets) == 1
        assert all(ligand.HasProp("_Name") for ligand in ligands)

        print()
        print(f"[Running Pipe Block: {self.name}]")
        
        return ligands, targets
    
    @property
    @abstractmethod
    def name(self) -> str:
        return "pipeblock"

class Pipeline(ABC):
    def __init__(self, return_input: bool = False, out_dir: str = "output", 
                 extra_info: dict = dict()):
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)
        self.blocks: List[PipelineBlock] = []
        self.extra_info: dict = extra_info
        self.return_input = return_input
        self.out_dir = str(pathlib.Path(out_dir).absolute())

    def run(self, ligands: List[Chem.Mol], targets: List[str],
            extra_info: dict = dict(), save_output: bool = True) -> Tuple[List[Chem.Mol], List[str]]:
        
        # if extra_info:
        #     self.extra_info = extra_info
        # else:
        #     extra_info = self.extra_info

        l = 55
        print()
        print(">"*l)
        print(">>" + str.center("STARTED: " + self.name, l-5), ">>")
        print(">"*l)
        origligs, origtargs = ligands, targets
        for block in self.blocks:
            ligands, targets = block.run(ligands, targets, extra_info)
        print()
        print("<"*l)
        print("<<" + str.center("ENDED: " + self.name, l-5), "<<")
        print("<"*l)

        self.extra_info = extra_info

        if save_output:
            prefix = os.path.join(self.out_dir, 
                                "".join(self.name.split()) + "Output")
            self.save_mid(ligands, targets, prefix+"Ligs.sdf", prefix+"Targs.txt")

        if not self.return_input:
            return ligands, targets
        return origligs, origtargs
    
    @staticmethod
    def save_mid(ligands, targets, ligname, targname):
        write_ligands(ligands, ligname)
        with open(targname, "w") as f:
            f.write("\n".join(targets))
    
    @property
    @abstractmethod
    def name(self) -> str:
        return "pipeline"
    
class LigandEnumerator(PipelineBlock, ABC):
    def __init__(self, debug: bool = False):
        super().__init__(debug)

    @abstractmethod
    def enumerate(self, ligand: Chem.Mol) -> List[Chem.Mol]:
        return [ligand]
    
    def run(self, ligands: List[Chem.Mol], targets: List[str], 
            extra_info: dict = dict()) -> Tuple[List[Chem.Mol] | List[str]]:
        ligands, targets = PipelineBlock.run(self, ligands, targets, extra_info)
        zipped = zip_complex(tqdm(ligands), targets)

        out_ligands = []
        out_targets = []
        for i, (ligand, target) in enumerate(zipped):
            ligand.SetIntProp("PrevIdx", i) # TODO: how to utilize this for debug?
            newligs = self.enumerate(ligand)
            newtargs = [target for _ in range(len(newligs))]
            out_ligands.extend(newligs)
            out_targets.extend(newtargs)
        
        if len(targets) == 1:
            return out_ligands, targets
        else:
            return out_ligands, out_targets
    
class LigandSelector(PipelineBlock, ABC):
    def __init__(self, identifier: str = "_Name", debug: bool = False):
        super().__init__(debug)
        self.identifier = identifier

    @abstractmethod
    def select(self, ligands: List[Chem.Mol]) -> List[Chem.Mol]:
        return ligands[:1]

    def run(self, ligands: List[Chem.Mol], targets: List[str], 
            extra_info: dict = dict()) -> Tuple[List[Chem.Mol] | List[str]]:
        ligands, targets = PipelineBlock.run(self, ligands, targets, extra_info)
        
        if self.identifier == "all":
            [ligand.SetIntProp("PrevIdx", i) for i, ligand in enumerate(tqdm(ligands))]
            out_ligands = self.select(ligands)
        else:
            groups: Dict[str, List] = dict()
            for i, ligand in enumerate(tqdm(ligands)):
                ligand.SetIntProp("PrevIdx", i)
                identifier = ligand.GetProp(self.identifier)
                if identifier not in groups.keys():
                    groups[identifier] = []
                groups[identifier].append(ligand)
        
            out_ligands = []
            for identifier, ligs in groups.items():
                newligs = self.select(ligs)
                out_ligands.extend(newligs) 
        
        if len(targets) == 1:
            return out_ligands, targets
        
        out_targets = []
        for lig in out_ligands:
            idx = lig.GetIntProp("PrevIdx")
            out_targets.append(targets[idx])
        return out_ligands, out_targets

class LigandConverter(PipelineBlock, ABC):
    def __init__(self, debug: bool = False):
        super().__init__(debug)
        
    @abstractmethod
    def convert(self, ligand: Chem.Mol) -> Chem.Mol:
        return ligand
    
    def run(self, ligands: List[Chem.Mol], targets: List[str], 
            extra_info: dict = dict()) -> Tuple[List[Chem.Mol] | List[str]]:
        ligands, targets = PipelineBlock.run(self, ligands, targets, extra_info)
        out_ligands = []
        for i, ligand in enumerate(tqdm(ligands)):
            ligand.SetIntProp("PrevIdx", i)
            out_ligands.append(self.convert(ligand))
        return out_ligands, targets

class LigandScorer(PipelineBlock, ABC):
    def __init__(self, score_name: str, debug: bool = False):
        super().__init__(debug)
        self._score_name = score_name

    @property
    @abstractmethod
    def score_name(self) -> str:
        return self._score_name
    
    @abstractmethod
    def score(self, ligand: Chem.Mol) -> float:
        return 0.0
    
    def run(self, ligands: List[Chem.Mol], targets: List[str], 
            extra_info: dict = dict()) -> Tuple[List[Chem.Mol] | List[str]]:
        ligands, targets = PipelineBlock.run(self, ligands, targets, extra_info)
        if self.score_name in extra_info.keys(): # TODO: prevent overwrites or at least manage them?
            print("Warning: Overwriting a score in extra_info...")
        extra_info[self.score_name] = []
        for i, ligand in enumerate(tqdm(ligands)):
            score = self.score(ligand)
            ligand.SetDoubleProp(self.score_name, score)
            extra_info[self.score_name].append(score)
        return ligands, targets

class TargetEnumerator(PipelineBlock, ABC):
    pass

class TargetSelector(PipelineBlock, ABC):
    pass
  
class TargetConverter(PipelineBlock, ABC):
    def __init__(self, debug: bool = False):
        super().__init__(debug)

    @abstractmethod
    def convert(self, target_path: str) -> str:
        return target_path
    
    def run(self, ligands: List[Chem.Mol], targets: List[str], 
            extra_info: dict = dict()) -> Tuple[List[Chem.Mol] | List[str]]:
        ligands, targets = PipelineBlock.run(self, ligands, targets, extra_info) 
        out_targets = []
        for target in tqdm(targets):
            out_targets.append(self.convert(target))
        return ligands, out_targets

class TargetScorer(ABC):
    pass

class ComplexEnumerator(PipelineBlock, ABC):
    def __init__(self, debug: bool = False):
        super().__init__(debug)

    @abstractmethod
    def enumerate(self, ligand: Chem.Mol, target_path: str) -> Tuple[List[Chem.Mol], List[str]]:
        return [ligand], [target_path]
    
    def run(self, ligands: List[Chem.Mol], targets: List[str], 
            extra_info: dict = dict()) -> Tuple[List[Chem.Mol] | List[str]]:
        ligands, targets = PipelineBlock.run(self, ligands, targets, extra_info)
        zipped = zip_complex(tqdm(ligands), targets)

        out_ligands = []
        out_targets = []
        for ligand, target in zipped:
            newligs, newtargs = self.enumerate(ligand, target)
            out_ligands.extend(newligs)
            out_targets.extend(newtargs)
        return out_ligands, out_targets

class ComplexSelector(PipelineBlock, ABC):
    pass

class ComplexConverter(PipelineBlock, ABC):
    def __init__(self, debug: bool = False):
        super().__init__(debug)

    @abstractmethod
    def convert(self, ligand: Chem.Mol, target_path: str) -> Tuple[Chem.Mol, str]:
        return ligand, target_path
    
    def run(self, ligands: List[Chem.Mol], targets: List[str], 
            extra_info: dict = dict()) -> Tuple[List[Chem.Mol] | List[str]]:
        ligands, targets = PipelineBlock.run(self, ligands, targets, extra_info)
        zipped = zip_complex(tqdm(ligands), targets)
        
        out_ligands = []
        out_targets = []
        for i, (ligand, target) in enumerate(zipped):
            newlig, newtarg = self.convert(ligand, target)
            out_ligands.append(newlig)
            out_targets.append(newtarg)
        return out_ligands, out_targets

class ComplexScorer(PipelineBlock, ABC):
    def __init__(self, score_name: str, debug: bool = False):
        super().__init__(debug)
        self._score_name = score_name

    @property
    @abstractmethod
    def score_name(self) -> str:
        return self._score_name
    
    @abstractmethod
    def score(self, ligand: Chem.Mol, target: str):
        pass
    
    def run(self, ligands: List[Chem.Mol], targets: List[str], 
            extra_info: dict = dict()) -> Tuple[List[Chem.Mol] | List[str]]:
        ligands, targets = PipelineBlock.run(self, ligands, targets, extra_info)
        zipped = zip_complex(tqdm(ligands), targets)

        if self.score_name in extra_info.keys(): # TODO: prevent overwrites or at least manage them?
            print("Warning: Overwriting a score in extra_info...")
        extra_info[self.score_name] = []

        for ligand, target in zipped:
            score = self.score(ligand, target)
            ligand.SetDoubleProp(self.score_name, score)
            extra_info[self.score_name].append(score)
        return ligands, targets
    



