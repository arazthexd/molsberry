from typing import Dict, List, Tuple

import parmed
from rdkit import Chem
from rdkit.Chem import AllChem

from ...core import (
    SimpleBlock, 
    MoleculeData, LigandData, ProteinData, StringData,
    PDBPathRep, StringRep, SDFPathRep,
    Representation, generate_path_in_dir,
)
from ...modules.parmed import ParmedMolRep

class ParmedMoleculeCombiner(SimpleBlock):
    name = "parmedmoleculecombiner"

    outputs = [("complex", MoleculeData, ParmedMolRep, False)]

    def __init__(self, num_inputs: int, debug=False, 
                 save_output=False, num_workers=None):  
        super().__init__(debug, save_output, num_workers)  
        self.num_inputs = num_inputs  
        self._inputs = self._generate_inputs()  

    @property  
    def inputs(self) -> List[tuple]:  
        return self._inputs  
    
    @property
    def batch_groups(self) -> List[Tuple[str]]:
        return [tuple(f"molecule_{i+1}" for i in range(self.num_inputs))]

    def _generate_inputs(self) -> List[tuple]:  
        return [(f"molecule_{i+1}", MoleculeData, ParmedMolRep, False) 
                for i in range(self.num_inputs)] 
    
    def operate(self, input_dict: Dict[str, Representation]):
        molecules = [input_dict[self.input_keys[i]].content 
            for i in range(self.num_inputs)]  
        merged_structures = sum(molecules[1:], start=molecules[0])
        merged_structures: parmed.Structure
        output = {"complex": ParmedMolRep(merged_structures)}
        return output
        