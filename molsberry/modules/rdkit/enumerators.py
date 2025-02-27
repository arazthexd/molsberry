from typing import Any, Dict, List
from abc import ABC, abstractmethod

from itertools import product

from rdkit import Chem
from rdkit.Chem import (
    rdDistGeom, rdForceFieldHelpers, rdmolops, AllChem, rdMolAlign, 
    rdMolDescriptors
)
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.EnumerateStereoisomers import (
    EnumerateStereoisomers, StereoEnumerationOptions, 
)
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.ML.Cluster import Butina # type: ignore

from ...core import (
    MoleculeData, Data, BatchedData, LigandData, 
    BatchedRep, Representation
)
from ...core import SimpleEnumeratorBlock

from .specific_reps import RDKitSmallMolRep
from .interface import RDKitInterface

class RDKitLigEnumeratorBlock(SimpleEnumeratorBlock, ABC):
    inputs = [
        ("ligands", LigandData, RDKitSmallMolRep, False)
    ]
    outputs = [
        ("ligands", LigandData, RDKitSmallMolRep, False)
    ]
    batch_groups = []

    @abstractmethod
    def enumerate(self, rdmol: Chem.Mol) -> List[Chem.Mol]:
        pass
    
    def operate(self, input_dict: Dict[str, Representation]) \
        -> Dict[str, Representation | BatchedRep]:

        rep: RDKitSmallMolRep = input_dict[self.input_keys[0]]
        rdmol = rep.content
        rdmols = self.enumerate(rdmol)
        main_out_key = self.output_keys[0]
        return {main_out_key: BatchedRep([self._get_out_rep(main_out_key)(mol) 
                                          for mol in rdmols])}

class RDKitLigandTautEnumerator(RDKitInterface, RDKitLigEnumeratorBlock):
    name = "rdtautenum_lig"
    display_name = "RDKIT Ligand Tautomer Enumerator"

    def __init__(self, max_tautomers: int = 4, minimize: bool = False,
                 flatten: bool = False, debug: bool = False, 
                 save_output: bool = False):
        super().__init__(flatten=flatten, debug=debug, save_output=save_output)

        self.enumerator = rdMolStandardize.TautomerEnumerator()
        self.enumerator.SetRemoveBondStereo(False)
        self.enumerator.SetRemoveSp3Stereo(False)
        self.enumerator.SetMaxTautomers(max_tautomers)

        self.minimize = minimize
    
    def enumerate(self, rdmol: Chem.Mol) -> List[Chem.Mol]:
        rdmol_no_h = Chem.RemoveHs(rdmol, updateExplicitCount=True)
        ts: List[Chem.Mol] = list(self.enumerator.Enumerate(rdmol_no_h))
        ts = sorted(ts, 
                    key=lambda t: self.enumerator.ScoreTautomer(t), 
                    reverse=True)
        
        if self.is_mol_3d(rdmol):
            ts = [Chem.AddHs(t) for t in ts]
            [rdDistGeom.EmbedMolecule(t) for t in ts]
            [self.sync_mol_flexible_rotors(t, rdmol) for t in ts]
            if self.minimize:
                [rdForceFieldHelpers.MMFFOptimizeMolecule(t) for t in ts]
        else:
            ts = [Chem.AddHs(t, addCoords=False) for t in ts]
        
        return ts

class RDKitLigandStereoEnumerator(RDKitInterface, RDKitLigEnumeratorBlock):
    name = "rdstereoenum_lig"
    display_name = "RDKIT Stereoisomer Enumerator"

    def __init__(self, 
                 tryEmbedding: bool = True, 
                 onlyUnassigned: bool = True,
                 flatten: bool = False, debug: bool = False,
                 save_output: bool = False):
        super().__init__(flatten=flatten, debug=debug, save_output=save_output)
        self.options = StereoEnumerationOptions(tryEmbedding=tryEmbedding, 
                                                onlyUnassigned=onlyUnassigned)
    
    def enumerate(self, rdmol: Chem.Mol) -> List[Chem.Mol]:
        rdmol_h = Chem.AddHs(rdmol)
        rdmol_h.RemoveAllConformers()
        newmols = list(EnumerateStereoisomers(rdmol_h, self.options, 
                                              not self.debug))
        
        newmols = [Chem.AddHs(newmol) for newmol in newmols]
        if self.is_mol_3d(rdmol):
            [rdDistGeom.EmbedMolecule(newmol) for newmol in newmols]
            [self.sync_mol_flexible_rotors(newmol, Chem.AddHs(rdmol)) 
             for newmol in newmols]

        return newmols

class RDKitLigandRingEnumerator(RDKitInterface, RDKitLigEnumeratorBlock):
    name = "rdringenum_lig"
    display_name = "RDKIT Ring Conformation Enumerator"
    # TODO: This class needs serious changes...
    # TODO: Implement the method in GypsumDL?

    def __init__(self, 
                 num_confs: int = 30, 
                 minimize: bool = False,
                 max_per_ring: int = 2, 
                 dist_threshold: float = 0.5, 
                 flatten: bool = False, debug: bool = False, 
                 save_output: bool = False):
        super().__init__(flatten=flatten, debug=debug, save_output=save_output)

        self.num_confs: int = num_confs
        self.minimize: bool = minimize
        self.max_per_ring: int = max_per_ring
        self.dist_threshold: float = dist_threshold

    def enumerate(self, rdmol: Chem.Mol) -> List[Chem.Mol]:
        rdmol = Chem.Mol(rdmol)
        rdDistGeom.EmbedMultipleConfs(rdmol, 
                                      numConfs=self.num_confs, 
                                      clearConfs=False)
        if self.minimize:
            rdmol = Chem.AddHs(rdmol, addCoords=True)
            rdForceFieldHelpers.MMFFOptimizeMoleculeConfs(rdmol)

        rdmolops.FindRingFamilies(rdmol)
        rinfo = rdmol.GetRingInfo()
        if rinfo.NumRings() == 0:
            return [rdmol]
        
        clusts_per_ring = []
        for i, ring_atoms in enumerate(rinfo.AtomRings()):
            if all(rdmol.GetAtomWithIdx(a).GetIsAromatic() 
                   for a in ring_atoms):
                continue
            dists = AllChem.GetConformerRMSMatrix(rdmol, prealigned=False, 
                                                  atomIds=ring_atoms)
            clusts = Butina.ClusterData(dists, rdmol.GetNumConformers(), 
                                        self.dist_threshold, 
                                        isDistData=True, reordering=True)
            clusts = [set(clust) for clust in clusts]
            clusts_per_ring.append(clusts[:self.max_per_ring])
        
        # print(clusts_per_ring)
        rings_combos = list(product(*clusts_per_ring))
        final_mols = []
        # print(rings_combos)
        for combo in rings_combos:
            if len(combo) == 0: 
                continue
            common: set = combo[0]
            for s in combo[1:]:
                common.intersection_update(s)
            # print(common)
            if common:
                final_mols.append(Chem.Mol(rdmol, confId=common.pop()))
        
        if final_mols:
            return final_mols
        else:
            return [rdmol]