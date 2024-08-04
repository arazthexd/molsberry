from typing import List, Tuple

from .representations import MOPACInputMolRep

MOPAC_TEMPLATE = """{keywords}
{desc1}
{desc2}
{coordinates}
"""

class MOPACConfig:
    def __init__(self):
        self.keywords: List[str] = ["PDBOUT"]
        self.desc1: str = "MOPACConfig "
        self.desc2: str = ""
        self.coordinate_lines: List[str] = list()
        self.charge = 0
        # TODO and NOTE: Currently being written to be compatible with only PDB
    

    def add_fragment(self, input_rep: MOPACInputMolRep):
        # self.keywords.extend(input_rep.keywords)
        # TODO: Any keywords coming from single fragments?
        self.charge += input_rep.charge
        self.desc2 += "| " + input_rep.description + " "
        num_current_atoms = len(self.coordinate_lines)
        self.coordinate_lines.extend([
            self._update_atom_numbers(line, num_current_atoms) for line in 
            input_rep.coordinates.split("\n")
        ])
    
    def get_config_str(self):
        return MOPAC_TEMPLATE.format(
            keywords = " ".join(self.keywords) + f" CHARGE={self.charge}",
            desc1 = self.desc1,
            desc2 = self.desc2,
            coordinates = "\n".join(self.coordinate_lines)
        )
    
    @staticmethod
    def _update_atom_numbers(line, num_current_atoms):
        atom_number = int(line[6:11])
        atom_number += num_current_atoms
        line = line[:6] + str(atom_number).rjust(5) + line[11:]
        return line
    
class MOPACMozymeConfig(MOPACConfig):
    def __init__(self):
        super().__init__()
        self.keywords.append("MOZYME")
        self.keywords.append("GEO-OK")
        self.desc1: str = "MOPACMozymeConfig "
        self.setpi: str = ""
        self.neg_cvb: List[Tuple[int, int]] = list()
    
    def add_fragment(self, input_rep: MOPACInputMolRep):
        num_current_atoms = len(self.coordinate_lines)
        self.setpi = "\n".join(
            [" ".join(
                [str(pi_pair[0]+num_current_atoms), 
                 str(pi_pair[1]+num_current_atoms)]
            ) for pi_pair in input_rep.setpi])
        self.neg_cvb.extend([
            (pair[0]+num_current_atoms, pair[1]+num_current_atoms) for pair
            in input_rep.neg_cvb
        ])
        super().add_fragment(input_rep)
    
    def get_config_str(self):
        old_keywords = self.keywords.copy()
        cvb_text = self._cvblist_to_cvbtxt(self.neg_cvb, mode="negative")
        # TODO: Also consider pos cvb
        self.keywords = self.keywords + [f"CVB({cvb_text})"]
        config = super().get_config_str()
        self.keywords = old_keywords.copy()
        return config + "\n" + self.setpi
    
    @staticmethod
    def _cvblist_to_cvbtxt(cvblist: List[Tuple[int, int]], 
                           mode: str = "negative", atnum_shift: int = 0):
        if mode == "negative":
            cvblist = [":-".join((str(p1+atnum_shift),str(p2+atnum_shift)))
                       for p1, p2 in cvblist]
            return ";".join(cvblist)
        else:
            raise NotImplementedError()

    
    