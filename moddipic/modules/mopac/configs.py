from typing import List

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
        # TODO and NOTE: Currently being written to be compatible with only PDB
    

    def add_fragment(self, input_rep: MOPACInputMolRep):
        self.keywords.extend(input_rep.keywords)
        self.desc2 += "| " + input_rep.description + " "
        num_current_atoms = len(self.coordinate_lines)
        self.coordinate_lines.extend([
            self._update_atom_numbers(line, num_current_atoms) for line in 
            input_rep.coordinates.split("\n")
        ])
    
    def get_config_str(self):
        return MOPAC_TEMPLATE.format(
            keywords = " ".join(self.keywords),
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
        self.desc1: str = "MOPACMozymeConfig "
        self.setpi: str = ""
    
    def add_fragment(self, input_rep: MOPACInputMolRep):
        num_current_atoms = len(self.coordinate_lines)
        self.setpi = "\n".join(
            [" ".join(
                [str(pi_pair[0]+num_current_atoms), 
                 str(pi_pair[1]+num_current_atoms)]
            ) for pi_pair in input_rep.setpi])
        # negcvb_list = []
        # for negcvb_a1, negcvb_a2 in input_rep.neg_cvb:
        #     negcvb_list.append(f"{negcvb_a1}:-{negcvb_a2}")

        # self.keywords.append(f"CVB()")
        super().add_fragment(input_rep)
    
    def get_config_str(self):
        config = super().get_config_str()
        return config + "\n" + self.setpi
    
    