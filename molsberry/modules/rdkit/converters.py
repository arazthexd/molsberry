from typing import Dict
from abc import ABC, abstractmethod
from tqdm import tqdm

from rdkit import Chem
from rdkit.Chem import rdDistGeom
from rdkit.Chem import rdMolDescriptors, rdForceFieldHelpers
from rdkit.ForceField import rdForceField

from molsberry.core.data import Representation

from ...core.data import MoleculeData, Data, BatchedData, LigandData
from ...core.templates import (
    SimpleBlock, PLOptimizeJob, OptimizeJob
)

from .representations import RDKitMolRep, RDKitSmallMolRep
from .interface import RDKitInterface
from .analyzers import RDKitMMFFEnergyCalculator

class RDKitConverterBlock(SimpleBlock, ABC):
    inputs = [
        ("molecules", MoleculeData, RDKitMolRep, False)
    ]
    outputs = [
        ("molecules", MoleculeData, RDKitMolRep, False)
    ]
    batch_groups = []

    @abstractmethod
    def convert(self, rdmol: Chem.Mol) -> Chem.Mol:
        pass

    def operate(self, input_dict: Dict[str, Representation]) \
        -> Dict[str, Representation]:

        rdmol = input_dict[self.input_keys[0]].content
        rdmol = self.convert(rdmol)
        main_out_key = self.output_keys[0]
        return {main_out_key: self._get_out_rep(main_out_key)(rdmol)}

class RDKitLigandConverterBlock(RDKitConverterBlock, ABC):
    inputs = [
        ("ligands", LigandData, RDKitSmallMolRep, False)
    ]
    outputs = [
        ("ligands", LigandData, RDKitSmallMolRep, False)
    ]
        
class RDKitHydrogenAdder(RDKitInterface, RDKitConverterBlock):
    name = "rdhadder"
    display_name = "RDKit Hydrogen Adder"

    def convert(self, rdmol: Chem.Mol) -> Chem.Mol:
        return self.addhs_based_on_confdim(rdmol)

class RDKitLigandHAdder(RDKitLigandConverterBlock, RDKitHydrogenAdder):
    name = "rdhadder_lig"
    display_name = "RDKit Ligand Hydrogen Adder"

class RDKitLigandEmbedder(RDKitInterface, RDKitLigandConverterBlock):
    name = "rdembedder_lig"
    display_name = "RDKit Ligand Embedder"

    def __init__(self, debug: bool = False, save_output: bool = False,
                 remove_hs: bool = False):
        super().__init__(debug, save_output)
        self.remove_hs = remove_hs

    def convert(self, rdmol: Chem.Mol) -> Chem.Mol:
        rdmol = Chem.AddHs(rdmol)
        rdDistGeom.EmbedMolecule(rdmol)
        if self.remove_hs:
            rdmol = Chem.RemoveHs(rdmol)
        return rdmol

class RDKitMMFFOptimizer(RDKitInterface, OptimizeJob, SimpleBlock):
    name = "rdmmff_opt"
    display_name = "RDKit MMFF Optimizer"
    inputs = [
        ("molecules", MoleculeData, RDKitMolRep, False)
    ]
    batch_groups = []
    mol_rep = RDKitMolRep
    optimize_keys = ["molecules"]

    def __init__(self, max_cycles: int = 200,
                 debug: bool = False, save_output: bool = False) -> None:
        SimpleBlock.__init__(self, debug, save_output)
        self.max_cycles = max_cycles
        self.en_calculator = RDKitMMFFEnergyCalculator()
        OptimizeJob.__init__(self, energy_fn=self.en_calculator.calc_energy,
                             optimize_fn=self.rdmol_optimize)
    
    def operate(self, input_dict: Dict[str, Representation]) \
        -> Dict[str, Representation]:
        
        rdmol = input_dict[self.input_keys[0]].content
        _, _, _, out_dict = self.optimize(rdmol)
        return {k: self._get_out_rep(k)(v) for k, v in out_dict.items()}

    def rdmol_optimize(self, rdmol: Chem.Mol) -> Chem.Mol:
        assert self.is_mol_3d(rdmol)
        rdmol = Chem.Mol(rdmol)
        props = rdForceFieldHelpers.MMFFGetMoleculeProperties(rdmol)
        ff: rdForceField.ForceField = \
            rdForceFieldHelpers.MMFFGetMoleculeForceField(rdmol, props)
        ff.Minimize(maxIts=self.max_cycles)
        return rdmol

class RDKitPLComplexOptimizer(PLOptimizeJob, RDKitMMFFOptimizer):
    name = "rdpl_opt"
    display_name = "RDKit Protein Ligand Complex Optimizer"
    lig_rep = RDKitMolRep
    prot_rep = RDKitMolRep

    def __init__(self, max_cycles: int = 200,
                 debug: bool = False, save_output: bool = False) -> None:
        SimpleBlock.__init__(self, debug, save_output)
        self.max_cycles = max_cycles
        self.en_calculator = RDKitMMFFEnergyCalculator()
        PLOptimizeJob.__init__(self, energy_fn=self.en_calculator.calc_energy,
                               optimize_fn=self.rdmol_optimize)
        
    def operate(self, input_dict: Dict[str, Representation]) \
        -> Dict[str, Representation]:
        
        rdlig = input_dict[self.input_keys[0]].content
        rdprot = input_dict[self.input_keys[1]].content
        out_dict = self.optimize(rdlig, rdprot)
        return {k: self._get_out_rep(k)(v) for k, v in out_dict.items()}
    
    def combine_mols(self, *mols):
        return RDKitMMFFOptimizer.combine_mols(self, *mols)

    def separate_mols(self, merged):
        mols = Chem.GetMolFrags(merged, asMols=True)
        ligand = mols[0]
        protein = self.combine_mols(*mols[1:])
        return ligand, protein
        
        


    