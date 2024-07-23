from typing import Any, Dict, List

from itertools import product

from rdkit import Chem
from rdkit.Chem import (
    rdDistGeom, rdForceFieldHelpers, rdmolops, AllChem, rdMolAlign, rdMolDescriptors
)
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.EnumerateStereoisomers import (
    EnumerateStereoisomers, StereoEnumerationOptions, 
)
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.ML.Cluster import Butina

from ..core.templates import LigandEnumeratorBlock
from ..utils.moltools import is_mol_3d, sync_mol_flexible_rotors

class RDKitTautEnumerator(LigandEnumeratorBlock):
    name = "RDKIT Tautomer Enumerator"

    def __init__(self, max_tautomers: int = 4, minimize: bool = False,
                 flatten: bool = False, debug: bool = False):
        super().__init__(flatten=flatten, debug=debug)

        self.enumerator = rdMolStandardize.TautomerEnumerator()
        self.enumerator.SetRemoveBondStereo(False)
        self.enumerator.SetRemoveSp3Stereo(False)
        self.enumerator.SetMaxTautomers(max_tautomers)

        self.minimize = minimize
    
    def enumerate(self, ligand: Chem.Mol) -> List[Chem.Mol]:
        ligand_no_h = Chem.RemoveHs(ligand, updateExplicitCount=True)
        ts: List[Chem.Mol] = list(self.enumerator.Enumerate(ligand_no_h))
        ts = sorted(ts, 
                    key=lambda t: self.enumerator.ScoreTautomer(t), 
                    reverse=True)
        
        if is_mol_3d(ligand):
            [rdDistGeom.EmbedMolecule(t) for t in ts]
            ts = [Chem.AddHs(t, addCoords=True) for t in ts]
            [sync_mol_flexible_rotors(t, ligand) for t in ts]
        else:
            ts = [Chem.AddHs(t, addCoords=False) for t in ts]
        
        if self.minimize:
            [rdForceFieldHelpers.MMFFOptimizeMolecule(t) for t in ts]
        
        return ts

class RDKitStereoEnumerator(LigandEnumeratorBlock):
    name = "RDKIT Stereoisomer Enumerator"

    def __init__(self, 
                 tryEmbedding: bool = True, 
                 onlyUnassigned: bool = True,
                 flatten: bool = False, debug: bool = False):
        super().__init__(flatten=flatten, debug=debug)
        self.options = StereoEnumerationOptions(tryEmbedding=tryEmbedding, 
                                                onlyUnassigned=onlyUnassigned)
    
    def enumerate(self, ligand: Chem.Mol) -> List[Chem.Mol]:
        ligand_h = Chem.AddHs(ligand)
        ligand_h.RemoveAllConformers()
        newligs = list(EnumerateStereoisomers(ligand_h, self.options, 
                                              not self.debug))
        
        newligs = [Chem.AddHs(newlig) for newlig in newligs]
        if is_mol_3d(ligand):
            [rdDistGeom.EmbedMolecule(newlig) for newlig in newligs]
            [sync_mol_flexible_rotors(newlig, Chem.AddHs(ligand)) 
             for newlig in newligs]

        return newligs

class RDKitRingEnumerator(LigandEnumeratorBlock):
    name = "RDKIT Ring Conformation Enumerator"
    # TODO: This class needs serious changes...

    def __init__(self, 
                 num_confs: int = 30, 
                 minimize: bool = False,
                 max_per_ring: int = 2, 
                 dist_threshold: float = 0.5, 
                 flatten: bool = False, debug: bool = False):
        super().__init__(flatten=flatten, debug=debug) 

        self.num_confs: int = num_confs
        self.minimize: bool = minimize
        self.max_per_ring: int = max_per_ring
        self.dist_threshold: float = dist_threshold

    def enumerate(self, ligand: Chem.Mol) -> List[Chem.Mol]:
        ligand = Chem.Mol(ligand)
        rdDistGeom.EmbedMultipleConfs(ligand, 
                                      numConfs=self.num_confs, 
                                      clearConfs=False)
        if self.minimize:
            ligand = Chem.AddHs(ligand, addCoords=True)
            rdForceFieldHelpers.MMFFOptimizeMoleculeConfs(ligand)

        rdmolops.FindRingFamilies(ligand)
        rinfo = ligand.GetRingInfo()
        if rinfo.NumRings() == 0:
            return [Chem.Mol(ligand, confId=0)]
        
        clusts_per_ring = []
        for i, ring_atoms in enumerate(rinfo.AtomRings()):
            if all(ligand.GetAtomWithIdx(a).GetIsAromatic() 
                   for a in ring_atoms):
                continue
            dists = AllChem.GetConformerRMSMatrix(ligand, prealigned=False, 
                                                  atomIds=ring_atoms)
            clusts = Butina.ClusterData(dists, ligand.GetNumConformers(), 
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
                final_mols.append(Chem.Mol(ligand, confId=common.pop()))
        
        # if ligand_ref.GetNumConformers() > 0:
        #     if ligand_ref.GetConformer().Is3D():
        #         [rdMolAlign.AlignMol(m, ligand_ref) for m in final_mols]
        
        if final_mols:
            return final_mols
        else:
            return [Chem.Mol(ligand, confId=0)]