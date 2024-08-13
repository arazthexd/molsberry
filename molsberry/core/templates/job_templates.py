from typing import Optional, Callable, Tuple, Any, Type, List, Dict
from abc import ABC, abstractmethod

from ..data import ProteinData, LigandData, Representation
from .simple_blocks import SimpleBlock

class PLInteractionJob(ABC):

    def __init__(self, energy_fn: Optional[Callable] = None):
        self.energy_fn = energy_fn

    @property
    @abstractmethod
    def lig_rep(self) -> Type[Representation]:
        pass

    @property
    @abstractmethod
    def prot_rep(self) -> Type[Representation]:
        pass

    @property
    def inputs(self) -> Tuple[Any]:
        return [
            ("ligands", LigandData, self.lig_rep, False),
            ("proteins", ProteinData, self.prot_rep, False)
        ]
    
    @property
    def outputs(self) -> Tuple[Any]:
        # TODO: Update numeric data after adding it to core.
        return [
            ("e_interaction", None, None, False),
            ("e_ligand", None, None, False),
            ("e_protein", None, None, False),
            ("e_complex", None, None, False)
        ]
    
    @property
    def batch_groups(self) -> List[Tuple[str]]:
        return [("ligands", "proteins")]
    
    @abstractmethod
    def combine_pl(self, ligand, protein):
        pass

    def calc_interaction(self, ligand, protein) -> Dict[str, float]:
        e_ligand = self.energy_fn(ligand)
        e_protein = self.energy_fn(protein)
        
        pl_complex = self.combine_pl(ligand, protein)
        e_complex = self.energy_fn(pl_complex)

        return {
            "e_interaction": e_complex - e_ligand - e_protein,
            "e_ligand": e_ligand,
            "e_protein": e_protein,
            "e_complex": e_complex
        }
    

