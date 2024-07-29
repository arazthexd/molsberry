from typing import Any, Dict, List

import os

from rdkit import Chem

from ...core.templates import LigandOperatorBlock

class LigandSaver(LigandOperatorBlock):
    name = "Ligand Saver Block"
    output_keys = []
    force_run: bool = True
    # TODO: Why does the saver cause every node to save outputs?

    def __init__(self, filename: str | None = None, 
                 debug: bool = False) -> None:
        super().__init__(debug=debug, save_output=True)
        if "." not in list(filename):
            self.filename = filename
        else:
            assert len(filename.split(".")) == 2
            self.filename = filename.split(".")[0]
    
    def execute(self, ligands: List[Chem.Mol]) -> Dict[str, Any]:
        return {"ligands": ligands}
    
    def _get_sdf_path(self) -> str:
        sdf_path = os.path.join(self.base_dir, 
                                f"{self.filename.replace(' ', '_').lower()}.sdf")
        return sdf_path