from typing import Dict

from openmm import app
from pdbfixer import PDBFixer
import parmed
from openff.toolkit import Molecule
from openff.units.openmm import to_openmm as openff_unit_to_openmm
from openmmforcefields.generators import SMIRNOFFTemplateGenerator

from ...core import SimpleBlock, MoleculeData, generate_path_in_dir
from ..openmm import DEFAULT_FORCEFIELDS
from .representations import ParmedMolRep
from ...core import PDBPathRep

class OpenMMProteinParameterizer(SimpleBlock):
    name = "openmm_prot_parammer"
    display_name = "OpenMM Protein Parameterizer"
    inputs = [
        ("proteins", MoleculeData, ParmedMolRep, False)
    ]
    outputs = [
        ("proteins", MoleculeData, ParmedMolRep, False),
    ]
    batch_groups = []

    def __init__(self, forcefield: app.ForceField | None = None, 
                 debug: bool = False, save_output: bool = False) -> None:
        SimpleBlock.__init__(self, debug=debug, save_output=save_output)
        if forcefield is None:
            forcefield = app.ForceField(*DEFAULT_FORCEFIELDS)
        self.forcefield = forcefield

    def operate(self, input_dict: Dict[str, ParmedMolRep]) \
        -> Dict[str, ParmedMolRep]:

        stprot = input_dict[self.input_keys[0]].content
        topprot = stprot.topology
        posprot = stprot.coordinates
        pdbfixer_input_path = generate_path_in_dir(6, self.base_dir, "_ommpp.pdb")
        app.PDBFile.writeFile(topprot, posprot, pdbfixer_input_path)
        protmodel = PDBFixer(pdbfixer_input_path)
        protmodel.addMissingHydrogens()

        sysprot = self.forcefield.createSystem(protmodel.topology)
        parammed_prot = parmed.openmm.load_topology(protmodel.topology, 
                                                    sysprot, 
                                                    protmodel.positions)
        
        key = self.output_keys[0]
        return {key: self._get_out_rep(key)(parammed_prot)}

class OpenFFSmallMoleculeParameterizer(SimpleBlock):
    name = "openff_small_mol_parameterizer"
    display_name = "OpenFF Small Molecule Parameterizer"
    inputs = [
        ("molecules", MoleculeData, ParmedMolRep, False)
    ]
    outputs = [
        ("molecules", MoleculeData, ParmedMolRep, False),
    ]
    batch_groups = []

    def __init__(self, forcefield: app.ForceField | None = None, 
                 debug: bool = False, save_output: bool = False) -> None:
        SimpleBlock.__init__(self, debug=debug, save_output=save_output)
        if forcefield is None:
            forcefield = app.ForceField(*DEFAULT_FORCEFIELDS)
        self.forcefield = forcefield
    
    def operate(self, input_dict: Dict[str, ParmedMolRep]) \
        -> Dict[str, ParmedMolRep]:

        stmol = input_dict[self.input_keys[0]].content
        rdmol = ParmedMolRep.parmed2rdkit(stmol)
        offmol = Molecule.from_rdkit(rdmol, hydrogens_are_explicit=True, allow_undefined_stereo=True)
        offmol.assign_partial_charges('mmff94')
        offtop = offmol.to_topology()
        offpos = offmol.conformers[0]
        ommtop = offtop.to_openmm()
        ommpos = openff_unit_to_openmm(offpos)

        smirnoff = SMIRNOFFTemplateGenerator(molecules=offmol)
        self.forcefield.registerTemplateGenerator(smirnoff.generator)

        ommsys = self.forcefield.createSystem(ommtop)
        parammed_struct = parmed.openmm.load_topology(ommtop, ommsys, ommpos)
        
        key = self.output_keys[0]
        return {key: self._get_out_rep(key)(parammed_struct)}
    
class PrmedUpdateCoordinates(SimpleBlock): # TODO: write!!!!
    name = "parmed_update_coordinates"
    display_name = "Parmed Coordinate Updater"
    inputs = [
        ("molecules", MoleculeData, PDBPathRep, False),
        ("molecules", MoleculeData,)
    ]
    outputs = [
        ("molecules", MoleculeData, ParmedMolRep, False),
    ]
    batch_groups = []

    def __init__(self, debug: bool = False, save_output: bool = False) -> None:
        SimpleBlock.__init__(self, debug=debug, save_output=save_output)
        raise NotImplementedError()
    
    def operate(self, input_dict: Dict[str, ParmedMolRep]) \
        -> Dict[str, ParmedMolRep]:
        pass
