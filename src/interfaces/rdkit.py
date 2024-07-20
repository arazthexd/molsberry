from itertools import product
from typing import List, Tuple

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

from ..abstract import (
    LigandSelector, 
    LigandEnumerator, LigandConverter,
    PipelineBlock
)

class RDKitLigandHydrogener(LigandConverter):
    name = "RDKIT Hydrogen Adder"
    def convert(self, ligand: Chem.Mol) -> Chem.Mol:
        return Chem.AddHs(ligand)

class RDKitTautEmumerator(LigandEnumerator):
    name = "RDKIT Tautomer Enumerator"
    def __init__(self, max_tautomers: int = 4, debug: bool = False):
        super().__init__(debug)
        self.enumerator = rdMolStandardize.TautomerEnumerator()
        self.enumerator.SetRemoveBondStereo(False)
        self.enumerator.SetRemoveSp3Stereo(False)
        self.enumerator.SetMaxTautomers(max_tautomers)

    def enumerate(self, ligand: Chem.Mol) -> List[Chem.Mol]:
        ligand = Chem.RemoveHs(ligand, updateExplicitCount=True)
        ts: List[Chem.Mol] = list(self.enumerator.Enumerate(ligand))
        ts = sorted(ts, key=lambda t: self.enumerator.ScoreTautomer(t), reverse=True)
        ts = [Chem.AddHs(t, addCoords=True) for t in ts]
        [rdForceFieldHelpers.MMFFOptimizeMolecule(t) for t in ts]
        return ts

class RDKitRingEnumerator(LigandEnumerator):
    name = "RDKIT Ring Conformation Enumerator"
    def __init__(self, num_confs: int = 30, minimize: bool = False,
                 max_per_ring: int = 2, dist_threshold: float = 0.5, debug: bool = False):
        super().__init__(debug) 
        self.num_confs: int = num_confs
        self.minimize: bool = minimize
        self.max_per_ring: int = max_per_ring
        self.dist_threshold: float = dist_threshold

    def enumerate(self, ligand: Chem.Mol) -> List[Chem.Mol]:
        ligand_ref = Chem.Mol(ligand)
        ligand = Chem.Mol(ligand)
        rdDistGeom.EmbedMultipleConfs(ligand, numConfs=self.num_confs, clearConfs=False)
        if self.minimize:
            ligand = Chem.AddHs(ligand, addCoords=True)
            rdForceFieldHelpers.MMFFOptimizeMoleculeConfs(ligand)

        rdmolops.FindRingFamilies(ligand)
        rinfo = ligand.GetRingInfo()
        if rinfo.NumRings() == 0:
            return [Chem.Mol(ligand, confId=0)]
        
        clusts_per_ring = []
        for i, ring_atoms in enumerate(rinfo.AtomRings()):
            if all(ligand.GetAtomWithIdx(a).GetIsAromatic() for a in ring_atoms):
                continue
            dists = AllChem.GetConformerRMSMatrix(ligand, prealigned=False, atomIds=ring_atoms)
            clusts = Butina.ClusterData(dists, ligand.GetNumConformers(), self.dist_threshold, 
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
        

class RDKitStereoEnumerator(LigandEnumerator):
    name = "RDKIT Stereoisomer Enumerator"
    def __init__(self, debug: bool = False, tryembed: bool = True, **kwargs):
        super().__init__(debug)
        self.options = StereoEnumerationOptions(tryEmbedding=tryembed, **kwargs)
    
    def enumerate(self, ligand: Chem.Mol) -> List[Chem.Mol]:
        ligand = Chem.Mol(ligand)
        ligand = Chem.AddHs(ligand)
        ligand.RemoveAllConformers()
        newligs = list(EnumerateStereoisomers(ligand, self.options, not self.debug))
        # newligs = [Chem.AddHs(newlig) for newlig in newligs]
        # [rdDistGeom.EmbedMolecule(newlig) for newlig in newligs] # TODO: Which one is better? First embed or last embed?
        return newligs

class RDKitMWLigSelector(LigandSelector):
    name = "RDKit Ligand Weight Filtering"
    def __init__(self, identifier: str = "all", max_wt: float = 600.0, min_wt: float = 30.0, 
                 debug: bool = False):
        super().__init__(identifier, debug)
        self.max_wt = max_wt
        self.min_wt = min_wt
    
    def select(self, ligands: List[Chem.Mol]) -> List[Chem.Mol]:
        selected_ligs = []
        for ligand in ligands:
            if rdMolDescriptors.CalcExactMolWt(ligand) > self.max_wt:
                if self.debug:
                    print(f"mol wt > {self.max_wt}: {Chem.MolToSmiles(ligand)}")
                continue
            if rdMolDescriptors.CalcExactMolWt(ligand) < self.min_wt:
                if self.debug:
                    print(f"mol wt < {self.min_wt}: {Chem.MolToSmiles(ligand)}")
                continue
            selected_ligs.append(ligand)
        
        return selected_ligs
    