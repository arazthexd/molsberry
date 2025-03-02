from typing import Dict
from abc import ABC, abstractmethod
from tqdm import tqdm

from rdkit import Chem
from rdkit.Chem import rdDistGeom
from rdkit.Chem import rdMolDescriptors, rdForceFieldHelpers
from rdkit.ForceField import rdForceField

from molsberry.core.data import Representation

from ...core import MoleculeData, Data, BatchedData, LigandData
from ...core import SimpleBlock, PLOptimizeJob, OptimizeJob

from .representations import RDKitMolRep
from .specific_reps import RDKitSmallMolRep
from .interface import RDKitInterface
from .calculators import RDKitMMFFEnergyCalculator

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