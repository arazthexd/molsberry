from typing import Dict, List, Tuple

from rdkit import Chem
from rdkit.Chem import AllChem

from ...core import (
    SimpleBlock, 
    MoleculeData, LigandData, ProteinData, StringData,
    PDBPathRep, StringRep, SDFPathRep,
    Representation, generate_path_in_dir,
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
        ("molecules", LigandData, RDKitSmallMolRep, False),
        ("smiles", StringData, StringRep, False)
    ]
    outputs = [
        ("molecules", LigandData, RDKitSmallMolRep, False)
    ]
    batch_groups = [("molecules", "smiles")]

    def operate(self, input_dict: Dict[str, Representation]):
        
        rdmol = input_dict[self.input_keys[0]].content
        smi = input_dict[self.input_keys[1]].content
        smimol = Chem.MolFromSmiles(smi)
        rdmol = AllChem.AssignBondOrdersFromTemplate(refmol=smimol, mol=rdmol)
        
        out_key = self.output_keys[0]
        return {out_key: self._get_out_rep(out_key)(rdmol)}
    
class RDKitProteinLigandCombiner(SimpleBlock):
    name = "rdprotligcombiner"
    display_name = "(RDKit) Protein Ligand Combiner"
    inputs = [
        ("ligands", LigandData, SDFPathRep, False),
        ("proteins", ProteinData, PDBPathRep, False)
    ]
    outputs = [
        ("complex", MoleculeData, PDBPathRep, False)
    ]
    batch_groups = []

    def operate(self, input_dict: Dict[str, Representation]):
        ligandrep = input_dict[self.input_keys[0]].content
        proteinrep = input_dict[self.input_keys[1]].content
        ligandmol = next(Chem.SDMolSupplier(ligandrep))
        proteinmol = Chem.MolFromPDBFile(proteinrep)
        combined_mol = Chem.CombineMols(ligandmol, proteinmol)
        combined_path =  generate_path_in_dir(5, self.base_dir,'.pdb')
        pdbwriter = Chem.PDBWriter(combined_path)
        pdbwriter.write(combined_mol)
        pdbwriter.close()
        output = {"complex": PDBPathRep(combined_path)}
        return output