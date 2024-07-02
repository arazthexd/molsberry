from __future__ import annotations

import os, pathlib
from os import path
from typing import List
from itertools import repeat
from tqdm import tqdm
from joblib import Parallel, delayed, parallel_backend
from mpire import WorkerPool

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.PropertyMol import PropertyMol
from pdbfixer import PDBFixer
from openmm.app import PDBFile

from .abstract import LigSetPreparator, TargetPreparator, LigandEnumerator
from .interfaces.rdkit import (
    RDKitTautEmumerator, RDKitStereoEnumerator, RDKitRingEnumerator
)
from .interfaces.durrantlab import DimorphiteProtoEnumerator
from .utils.moltools import sync_mol_flexible_rotors

class DefaultProteinPrepper(TargetPreparator):
    def __init__(self, ph: float = 7.4, chains: str = "all", keep_water: bool = False,
                 add_missing_residues: bool = False, out_dir: str = ".", save_prefix: str = "prepared"):
        self.chains = chains
        self.keep_water = keep_water
        self.ph = ph
        self.add_missing_residues = add_missing_residues

        if not os.path.exists(out_dir):
            os.mkdir(out_dir)
        self.out_dir: str = str(pathlib.Path(out_dir).absolute())
        self.save_prefix: str = save_prefix

    def prepare(self, target_path: str) -> str:
        fixer: PDBFixer = PDBFixer(target_path)
        if self.chains != "all":
            rem_chains = [chain.id for chain in fixer.topology.chains() 
                        if chain.id not in self.chains]
            fixer.removeChains(chainIds=rem_chains)
        fixer.removeHeterogens(self.keep_water)
        if self.add_missing_residues:
            fixer.findMissingResidues()
            fixer.findMissingAtoms()
            fixer.addMissingAtoms()
        # fixer.findNonstandardResidues()
        # fixer.replaceNonstandardResidues()
        fixer.addMissingHydrogens(self.ph)  # TODO: Any need to pH value?

        save_path = path.join(self.out_dir, f"{self.save_prefix}_protein.pdb")
        PDBFile.writeFile(fixer.topology, fixer.positions, open(save_path, 'w')) 
        return save_path

class DefaultLigSetPrepper(LigSetPreparator):
    def __init__(
            self, 
            tau_enum: LigandEnumerator = RDKitTautEmumerator(), 
            stereo_enum: LigandEnumerator = RDKitStereoEnumerator(), 
            ring_enum: LigandEnumerator = RDKitRingEnumerator(minimize=True),
            prot_enum: LigandEnumerator = DimorphiteProtoEnumerator(align_mols=True),
            debug: bool = False, ph: float = 7.4, max_wt: float = 600.0, num_jobs: int = 4):
        self.te = tau_enum 
        self.se = stereo_enum
        self.re = ring_enum
        self.pe = prot_enum
        # if "debug" in self.te.__dict__.keys():
        #     self.te.debug = debug
        # if "debug" in self.se.__dict__.keys():
        #     self.se.debug = debug
        # if "debug" in self.re.__dict__.keys():
        #     self.re.debug = debug
        # if "debug" in self.pe.__dict__.keys():
        #     self.pe.debug = debug

        self.ph = ph
        self.debug = debug
        self.max_wt = max_wt
        self.num_jobs = num_jobs

    def prepare(self, ligands: List[Chem.Mol]) -> List[Chem.Mol]:
        # prepped_ligs = []
        # for lig in tqdm(ligands):
        #     newligs = self.prepare_single(lig)
        #     if newligs:
        #         prepped_ligs.extend(newligs)
        # return prepped_ligs
        # with parallel_backend('threading', n_jobs=n_threads):
        #     prepped_lig_sets = Parallel()(
        #         delayed(self.prepare_single)(ligand) for ligand in ligands)
        ligands = [PropertyMol(lig) for lig in ligands]
        with WorkerPool(n_jobs=self.num_jobs) as pool:
            results = pool.map(self.prepare_single, ligands, progress_bar=True)
        
        prepped_ligs = []
        for lig_set in results:
            if lig_set is None:
                continue
            for lig in lig_set:
                prepped_ligs.append(lig)
        
        return prepped_ligs
    
    def prepare_single(self, ligand: Chem.Mol) -> List[Chem.Mol]:
        
        if rdMolDescriptors.CalcExactMolWt(ligand) > self.max_wt:
            if self.debug: print(f"ligand {Chem.MolToSmiles(ligand)} wt > {self.max_wt}")
            return None
        
        if ligand.GetNumConformers() > 0:
            if ligand.GetConformer().Is3D():
                Chem.AssignStereochemistryFrom3D(ligand)
        
        initial_ligand = Chem.Mol(ligand)
        initial_smiles = Chem.MolToSmiles(ligand)
        ligands = [ligand]
        
        if self.se:
            ligands = [mol for ligand in ligands 
                       for mol in self.se.enumerate(ligand)]
            if self.debug: print(initial_smiles, "after stereo", len(ligands))

        if self.pe:
            ligands = [mol for ligand in ligands 
                       for mol in self.pe.enumerate(ligand, min_ph=self.ph-0.2, max_ph=self.ph+0.2)]
            if self.debug: print(initial_smiles, "after proto", len(ligands))
        
        if self.te:
            ligands = [mol for ligand in ligands 
                       for mol in self.te.enumerate(ligand)]
            if self.debug: print(initial_smiles, "after tauto", len(ligands))

        if self.re:
            ligands = [mol for ligand in ligands 
                       for mol in self.re.enumerate(ligand)]
            if self.debug: print(initial_smiles, "after ring", len(ligands))

        # ref_3D = Chem.Mol(ligands[0])
        # if initial_ligand.GetNumConformers() > 0:
        #     if initial_ligand.GetConformer().Is3D():
        #         ref_3D = initial_ligand
            
        # [sync_mol_flexible_rotors(ligand, ref_3D) for ligand in ligands]
        if self.debug: print(f"enumerated {initial_smiles} -> {len(ligands)}")
        return ligands
    
    @staticmethod
    def _prepare_multiprocess_helper(self: DefaultLigSetPrepper, ligand: Chem.Mol):
        return self.prepare_single(ligand)
        
