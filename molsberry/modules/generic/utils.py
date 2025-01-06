from typing import Dict, List

from rdkit import Chem

from ...core import (
    SimpleBlock, MoleculeData, LigandData, ProteinData,
    PDBPathRep,
    Representation
)
from ..rdkit import RDKitMolRep, RDKitSmallMolRep

class RDKitProteinLigandSplitter(SimpleBlock):
    name = "rdprotligsplitter"
    display_name = "(RDKit) Protein Ligand Splitter"
    inputs = [
        ("complex", MoleculeData, RDKitMolRep, False)
    ]
    outputs = [
        ("ligands", LigandData, RDKitSmallMolRep, False),
        ("proteins", ProteinData, RDKitMolRep, False)
    ]
    batch_groups = []

    def operate(self, input_dict: Dict[str, Representation]) \
        -> Dict[str, Representation]:

        rdcomp = input_dict[self.input_keys[0]].content
        rdfrags: List[Chem.Mol] = Chem.GetMolFrags(rdcomp, asMols=True)
        rdlig = Chem.Mol()
        rdprot = Chem.Mol()
        for rdfrag in rdfrags:
            if rdfrag.GetAtomWithIdx(0).GetPDBResidueInfo().GetIsHeteroAtom():
                rdlig = Chem.CombineMols(rdlig, rdfrag)
            else:
                rdprot = Chem.CombineMols(rdprot, rdfrag)

        lig_key = self.output_keys[0]
        prot_key = self.output_keys[1]

        return {
            lig_key: self._get_out_rep(lig_key)(rdlig),
            prot_key: self._get_out_rep(prot_key)(rdprot)
        }

class RDKitBondOrderAssigner(SimpleBlock):
    name = "rdbondorder"
    display_name = "(RDKit) Bond Order Assigner From SMILES"
    inputs = [
        ("molecules", MoleculeData, RDKitMolRep, False),
        ("smiles", )
    ]
    outputs = [
        ("ligands", LigandData, RDKitSmallMolRep, False),
        ("proteins", ProteinData, RDKitMolRep, False)
    ]
    batch_groups = []