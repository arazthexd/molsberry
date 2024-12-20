from typing import Optional, Callable, Tuple, Any, Type, List, Dict
from abc import ABC, abstractmethod

from ..data import ProteinData, LigandData, Representation, MoleculeData, NumericData, FloatRep
from .batchop import BatchOperatorBlock

from ...core import NumericData, FloatRep

class PLJob(ABC):
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
    def batch_groups(self) -> List[Tuple[str]]:
        return [("ligands", "proteins")]

class EnergyJob(ABC):
    def __init__(self, energy_fn: Optional[Callable]) -> None:
        self.energy_fn = energy_fn

    @property
    def outputs(self) -> List[Tuple[Any]]:
        return [
            ("energy", NumericData, FloatRep, False)
        ]
    
    def calc_energy(self, mol, *args, **kwargs):        
        en = self.energy_fn(mol, *args, **kwargs)
        return en
    
class InteractionJob(EnergyJob, ABC):
    def __init__(self, energy_fn: Callable[..., Any] | None) -> None:
        super().__init__(energy_fn)
        assert len(self.optimize_keys) == 2
    
    @property
    @abstractmethod
    def optimize_keys(self):
        pass

    @property
    def outputs(self) -> List[Tuple[Any]]:
        # TODO: Update numeric data after adding it to core.
        return [("e_interaction", NumericData, FloatRep, False),
                ("e_combined", NumericData, FloatRep, False)] + \
            self._get_outputs(self.optimize_keys)
    
    @staticmethod
    def _get_outputs(en_keys: List[str]) -> List[Tuple[Any]]:
        keys = [en_keys] if isinstance(en_keys, str) else en_keys
        outs = []
        for en_key in keys:
            outs.extend([ # TODO: Update numeric data after adding it to core.
                (f"e_{en_key}", None, None, False)
            ])
        return outs
    
    @abstractmethod
    def combine_mols(self, mol1, mol2):
        pass

    def calc_interaction(self, mol1, mol2, *args, **kwargs):
        e1 = self.calc_energy(mol1, *args[2:], **kwargs)
        e2 = self.calc_energy(mol2, *args[2:], **kwargs)
        merged = self.combine_mols(mol1, mol2)
        em = self.calc_energy(merged, *args[2:], **kwargs)

        return {
            self.outputs[0][0]: em - e1 - e2,
            self.outputs[1][0]: em,
            self.outputs[2][0]: e1,
            self.outputs[3][0]: e2
        }

class PLInteractionJob(PLJob, InteractionJob, ABC):
    optimize_keys = ["ligand", "protein"]

class OptimizeJob(ABC):
    def __init__(self, energy_fn: Callable[..., Any] | None,
                 optimize_fn: Callable[..., Any]) -> None:
        self.energy_fn = energy_fn
        self.optimize_fn = optimize_fn
    
    @property
    @abstractmethod
    def mol_rep(self):
        pass

    @property
    @abstractmethod
    def optimize_keys(self):
        pass
    
    @property
    def outputs(self):
        return [
            ("molecules", MoleculeData, self.mol_rep, False),
            ("e_init", NumericData, FloatRep, False),
            ("e_final", NumericData, FloatRep, False),
        ]
    
    @abstractmethod
    def combine_mols(self, *args):
        pass

    @abstractmethod
    def separate_mols(self, merged):
        pass

    def optimize(self, *mols):
        i_energies = []
        for mol in mols:
            i_energies.append(self.energy_fn(mol))
        
        merged = self.combine_mols(*mols)
        mols = self.separate_mols(self.optimize_fn(merged))

        f_energies = []
        for mol in mols:
            f_energies.append(self.energy_fn(mol))
        
        out = {}
        for i, mol in enumerate(mols):
            out[self.optimize_keys[i]] = mol
        out["e_init"] = sum(i_energies)
        out["e_final"] = sum(f_energies)
        return mols, i_energies, f_energies, out

class PLOptimizeJob(PLJob, OptimizeJob, ABC):
    optimize_keys = ["ligand", "protein"]

    @property
    def outputs(self):
        return [
            ("ligand", LigandData, self.lig_rep, False),
            ("protein", ProteinData, self.lig_rep, False),
            ("e_ligand_init", NumericData, FloatRep, False),
            ("e_ligand_final", NumericData, FloatRep, False),
            ("e_protein_init", NumericData, FloatRep, False),
            ("e_protein_final", NumericData, FloatRep, False),
        ]
    
    def optimize(self, ligand, protein):
        (ligand, protein), i_energies, f_energies, _ = \
            super().optimize(ligand, protein)
        return {
            "ligand": ligand,
            "protein": protein,
            "e_ligand_init": i_energies[0],
            "e_ligand_final": f_energies[0],
            "e_protein_init": i_energies[1],
            "e_protein_final": f_energies[1]
        }

# TODO: Create an Optimize and a PLOptimize job.