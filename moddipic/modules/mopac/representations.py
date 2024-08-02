from __future__ import annotations
from typing import Any, List

import os
import pathlib

from ...core.data.abstract import Representation
try:
    from rdkit import Chem
    from ...modules.rdkit.representations import RDKitMolRep
    RDKIT_SUCCESS = True
except: 
    RDKIT_SUCCESS = False

class MOPACInputRep(Representation):
    rep_name = "mopac_input"

    def __init__(self, keywords: List[str], coordinates: str,
                  desciption: str = "some description") -> None:
        data = {
            "keywords": keywords,
            "description": desciption,
            "coordinates": coordinates
        }
        super().__init__(data)
        self.keywords = keywords
        self.coordinates = coordinates
        self.description = desciption
    
    if RDKIT_SUCCESS:
        @classmethod
        def from_RDKitMolRep(cls, rdkit_rep: RDKitMolRep):
            rdmol = Chem.Mol(rdkit_rep.data)
            coordinates = Chem.MolToPDBBlock(rdmol)
            charge = Chem.GetFormalCharge(rdmol)
            return MOPACInputRep(
                keywords=[f"CHARGE={charge}"],
                coordinates=coordinates,
                desciption="Converted from rdkit molecule representation."
            )
            
        
