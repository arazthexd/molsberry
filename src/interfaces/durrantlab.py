from typing import List, Tuple
import contextlib, io

from rdkit import Chem
from rdkit.Chem import (
    rdDistGeom, rdForceFieldHelpers
)

try:
    import dimorphite_dl as dd
except ImportError:
    pass

from ..abstract import LigandEnumerator, PipelineBlock
from ..utils.moltools import sync_mol_flexible_rotors

out_io = io.StringIO()

class DimorphiteProtoEnumerator(LigandEnumerator):
    name = "DimorphiteDL Protomer Enumerator"
    def __init__(self, debug: bool = False, align_mols: bool = True, min_ph: float = 7.2,
                 max_ph: float = 7.6, **kwargs):
        super().__init__(debug) # debug should go here
        self.align_mols: bool = align_mols
        self.min_ph = min_ph
        self.max_ph = max_ph
        self.kwargs = kwargs

    def enumerate(self, ligand: Chem.Mol) -> List[Chem.Mol]:
        
        ligand = Chem.RemoveAllHs(ligand)
        
        if not self.debug:
            with contextlib.redirect_stdout(out_io):
                protomers = dd.run_with_mol_list(
                    [ligand], silent=not self.debug, min_ph=self.min_ph, 
                    max_ph=self.max_ph, **self.kwargs
                )
        else:
            protomers = dd.run_with_mol_list(
                [ligand], silent=not self.debug, min_ph=self.min_ph, 
                max_ph=self.max_ph, **self.kwargs
            )

        if len(protomers) == 0:
            return [ligand]
        
        if self.align_mols:
            protomers = [Chem.RemoveAllHs(protomer) for protomer in protomers]
            [protomer.AddConformer(ligand.GetConformer()) for protomer in protomers]
            protomers = [Chem.AddHs(protomer, addCoords=True) for protomer in protomers]
        else:
            [rdDistGeom.EmbedMolecule(protomer) for protomer in protomers]
        
        [protomer.SetProp("_Name", ligand.GetProp("_Name")) for protomer in protomers]
        return protomers

    # def enumerate_multi(self, ligands: List[Chem.Mol]) -> List[Chem.Mol]:

    #     ligands = [Chem.RemoveAllHs(ligand) for ligand in ligands]
    #     [ligand.SetIntProp("dimorphite_input_idx", i) for i, ligand in enumerate(ligands)]
        
    #     if not self.debug:
    #         with contextlib.redirect_stdout(out_io):
    #             protomers = dd.run_with_mol_list(
    #                 ligands, silent=not self.debug, min_ph=self.min_ph, 
    #                 max_ph=self.max_ph, **self.kwargs
    #             )
    #     else:
    #         protomers = dd.run_with_mol_list(
    #             ligands, silent=not self.debug, min_ph=self.min_ph, 
    #             max_ph=self.max_ph, **self.kwargs
    #         )
        
    #     # protomers = [Chem.RemoveAllHs(protomer) for protomer in protomers]
    #     [rdDistGeom.EmbedMolecule(protomer) for protomer in protomers]
    #     # [protomer.AddConformer(ligands[protomer.GetIntProp("dimorphite_input_idx")].GetConformer()) 
    #     #  for protomer in protomers]
    #     protomers = [Chem.AddHs(protomer, addCoords=True) for protomer in protomers]
    #     if self.align_mols:
    #         for protomer in protomers:
    #             ref = ligands[protomer.GetIntProp("dimorphite_input_idx")]
    #             if self.debug: print(protomer.GetIntProp("dimorphite_input_idx"))
    #             if self.debug: print(f"syncing protomer {Chem.MolToSmiles(protomer)} with ref {Chem.MolToSmiles(ref)}")
    #             sync_mol_flexible_rotors(protomer, ref)
    #         # [sync_mol_flexible_rotors(
    #         #     protomer, ligands[protomer.GetIntProp("dimorphite_input_idx")]) 
    #         #     for protomer in protomers]
    #     # [rdForceFieldHelpers.MMFFOptimizeMolecule(protomer) for protomer in protomers]
    #     [protomer.SetProp(
    #         "_Name", ligands[protomer.GetIntProp("dimorphite_input_idx")].GetProp("_Name"))
    #     for protomer in protomers]
    #     return protomers
