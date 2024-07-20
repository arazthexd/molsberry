from abc import ABC, abstractmethod
from typing import List, Tuple
import pathlib, shutil
from rdkit import Chem
import os, subprocess
import string
import random
from copy import deepcopy

# TODO: Should we make the following requisites?
from openmm.app import ForceField
from openmmforcefields.generators import SMIRNOFFTemplateGenerator
from openff.toolkit import Topology, Molecule
from openff.units import unit
from .openmm import OpenMMComponent

from ..abstract import Interface
from ..utils.iotools import write_ligands, save_pl_complex

CUBY4_MAIN_TEMPLATE = """
geometry: {geometry}
cuby_threads: {cuby_threads}"""

CUBY4_INTERACTION_ADDON = """
job: interaction
molecule_a:
  selection: {selection_a}
  {optional_a}
molecule_b:
  selection: {selection_b}
  {optional_b}"""

CUBY4_OPTIMIZE_ADDON = """
job: optimize
maxcycles: {max_cycles}
restart_file: {restart_file}"""

CUBY4_MOPAC_MAIN_ADDON = """
interface: mopac
method: {method}
mopac_exe: {mopac_exe}
mopac_mozyme: {mopac_mozyme}
mopac_corrections: {mopac_corrections}
solvent_model: {solvent_model}"""

CUBY4_MOPAC_MOLSPEC_ADDON = """
mopac_setcharge:
  {mopac_setcharge}
mopac_setpi: [{mopac_setpi}]
charge: {charge}"""

CUBY4_QMMM_MAIN_ADDON = """
interface: qmmm
qmmm_core: {qmmm_core}
qmmm_embedding: {qmmm_embedding}
gradient_on_point_charges: {grad_on_points}
calculation_qm:
  gradient_on_point_charges: {grad_on_points}
  {qm_config}
calculation_mm:
  {mm_config}
calculation_qmregion_mm:
  {qmregion_mm_config}"""

CUBY4_AMBER_MAIN_ADDON = """
interface: amber
amber_amberhome: {amber_home}
{topology_config}"""


class Cuby4Config(ABC):
    def __init__(self):
        pass

    @abstractmethod
    def get_config_str(self) -> str:
        pass

    @abstractmethod
    def update_config_for(self, some_input: Chem.Mol | str | Tuple[Chem.Mol | str, str], 
                          **kwargs):
        pass

class Cuby4JobConfig(Cuby4Config, ABC):
    def __init__(self, geometry: str, cuby_threads: int = 1):
        self.geometry = geometry
        self.cuby_threads = cuby_threads
    
    def get_main_str(self):
        return CUBY4_MAIN_TEMPLATE.format(
            geometry=self.geometry,
            cuby_threads=self.cuby_threads
        )

    @abstractmethod
    def get_job_str(self):
        pass

    def get_config_str(self):
        return self.get_job_str()
    
    # def update_config_for(self, some_input: Chem.Mol | str, work_dir: str = "tmp", **kwargs):
    #     if isinstance(some_input, Chem.Mol):
    #         work_dir = str(pathlib.Path(work_dir).absolute())
    #         tmp_name = "".join(random.choices(string.ascii_uppercase + string.digits, k=10))
    #         tmp_save_path = os.path.join(work_dir, f"{tmp_name}.sdf")
    #         write_ligands([some_input], save_path=tmp_save_path)
    #         self.geometry = tmp_save_path
    #     elif isinstance(some_input, str) and some_input.endswith(".pdb"):
    #         self.geometry = str(pathlib.Path(some_input).absolute())
    #     else:
    #         raise NotImplementedError()

    def update_config_for(self, some_input: Chem.Mol | str | Tuple[Chem.Mol | str, str], 
                          work_dir: str = "tmp", **kwargs):
        if isinstance(some_input, Chem.Mol):
            work_dir = str(pathlib.Path(work_dir).absolute())
            tmp_name = "smallmol_" + "".join(random.choices(string.ascii_uppercase + string.digits, k=10))
            tmp_save_path = os.path.join(work_dir, f"{tmp_name}.sdf")
            write_ligands([some_input], save_path=tmp_save_path)
            self.geometry = tmp_save_path
        elif isinstance(some_input, str) and some_input.endswith(".pdb"):
            self.geometry = str(pathlib.Path(some_input).absolute())
        elif isinstance(some_input, tuple):
            work_dir = str(pathlib.Path(work_dir).absolute())
            ligand = some_input[0]
            target = some_input[1]
            comp_path = os.path.join(work_dir, "cuby_complex_geometry_tmp.pdb")
            save_pl_complex(ligand, target, comp_path) # TODO: better implementation of job name
            self.geometry = comp_path
        else:
            raise NotImplementedError()
        
        return self.geometry

class Cuby4InterfaceConfig(Cuby4Config, ABC):
    def __init__(self):
        pass
    
    @abstractmethod
    def get_interface_str(self) -> str:
        pass

    def get_config_str(self):
        return self.get_interface_str()


class Cuby4InteractionConfig(Cuby4JobConfig):
    def __init__(self, geometry: str, selection_a: str, selection_b: str, 
                 cuby_threads: int = 1):
        super().__init__(geometry, cuby_threads)
        self.selection_a = selection_a
        self.selection_b = selection_b
    
    def get_interact_str(self, opt_a: Cuby4Config = None, opt_b: Cuby4Config = None):
        if opt_a == None:
            opt_a_str = ""
        else:
            opt_a_str = opt_a.get_config_str()
            opt_a_str = opt_a_str.replace("\n", "\n  ")
        if opt_b == None:
            opt_b_str = ""
        else:
            opt_b_str = opt_b.get_config_str()
            opt_b_str = opt_b_str.replace("\n", "\n  ")
        
        return CUBY4_INTERACTION_ADDON.format(
            selection_a=self.selection_a,
            selection_b=self.selection_b,
            optional_a=opt_a_str,
            optional_b=opt_b_str
        )
    
    def get_job_str(self, opt_a: Cuby4Config = None, opt_b: Cuby4Config = None):
        return self.get_main_str() + self.get_interact_str(opt_a, opt_b)
    
    # def get_config_str(self):
    #     return self.get_job_str()

class Cuby4EnergyConfig(Cuby4JobConfig):
    def __init__(self, geometry: str, cuby_threads: int = 1):
        super().__init__(geometry, cuby_threads)
    
    def get_energy_str(self):
        return "\njob: energy"
    
    def get_job_str(self):
        return self.get_main_str() + self.get_energy_str()

class Cuby4OptimizeConfig(Cuby4JobConfig):
    def __init__(self, geometry: str, restart_file: str = "", max_cycles: int = 200,
                 cuby_threads: int = 1):
        super().__init__(geometry, cuby_threads)
        self.restart_file = restart_file
        self.max_cycles = max_cycles
    
    def get_optimize_str(self):
        return CUBY4_OPTIMIZE_ADDON.format(
            restart_file=self.restart_file,
            max_cycles=self.max_cycles
        )
    
    def get_job_str(self):
        return self.get_main_str() + self.get_optimize_str()

class Cuby4MOPACConfig(Cuby4InterfaceConfig):
    def __init__(self, method: str = "pm6", mopac_exe: str = "auto", mozyme: bool = True,
                 corrections: str = "d3h4x", solvent: str = "none"):
        # super().__init__(job, geometry, "mopac", cuby_threads)
        self.method = method
        if mopac_exe == "auto":
            self.mopac_exe = str(pathlib.Path(shutil.which("mopac")).absolute())
        else:
            self.mopac_exe = mopac_exe
        self.mozyme = "yes" if mozyme else "no"
        self.corrections = corrections
        self.solvent = solvent
    
    def get_mopac_main_str(self):
        return CUBY4_MOPAC_MAIN_ADDON.format(
            method=self.method,
            mopac_exe=self.mopac_exe,
            mopac_mozyme=self.mozyme,
            mopac_corrections=self.corrections,
            solvent_model=self.solvent
        )

    def get_interface_str(self) -> str:
        return self.get_mopac_main_str()
    
    def update_config_for(self, some_input: Chem.Mol | str | Tuple[Chem.Mol | str, str], 
                          **kwargs):
        return

class Cuby4MOPACMolSpecConfig(Cuby4Config):
    def __init__(self, setpi: str, setcharge: str, charge: int = 0):
        self.setpi = setpi
        self.setcharge = setcharge
        self.charge = str(charge)
        self.keywords = "\nmopac_keywords: "
    
    def get_mopac_mol_str(self):
        return CUBY4_MOPAC_MOLSPEC_ADDON.format(
            charge=self.charge,
            mopac_setcharge=self.setcharge,
            mopac_setpi=self.setpi
        ) + self.keywords
    
    @classmethod
    def from_rdmol(cls, mol: Chem.Mol):
        setcharge_list = cls.get_mol_setcharge(mol)
        setcharge = "\n  ".join(setcharge_list)
        setpi = ",".join(cls.get_mol_setpi(mol))
        charge = Chem.GetFormalCharge(mol)
        out = cls(setpi, setcharge, charge)
        out.keywords += f"GEO-OK CVB({';'.join(cls.get_mol_cvb(mol))})"
        # print(out.keywords)
        return out

    @classmethod
    def from_pdb(cls, pdb_path: str, unique_molecules: List[Chem.Mol] = list()):
        setcharge_list = cls.get_pdb_setcharge(pdb_path, unique_molecules=unique_molecules)
        setcharge = "\n  ".join(setcharge_list)
        charge = sum([int(sline.split(": ")[1]) for sline in setcharge_list])
        return cls("", setcharge, charge) # TODO

    @staticmethod
    def get_mol_setcharge(mol: Chem.Mol) -> list:
        def helper(x):
            if x == 0: return str(x)
            if x > 0:
                return str("'+'") * x
            if x < 0:
                return str("'-'") * (-x)
        setcharge = [f"{a.GetIdx()+1}: {helper(a.GetFormalCharge())}" 
                    for a in mol.GetAtoms()]
        return setcharge
    
    @staticmethod
    def get_mol_cvb(mol: Chem.Mol) -> list:
        cvb_list = []
        for atom1 in range(mol.GetNumAtoms()):
            for atom2 in range(atom1+1,mol.GetNumAtoms()):
                # if mol.GetBondBetweenAtoms(atom1, atom2):
                #     cvb_list.append(f"{atom1+1}:{atom2+1}")
                # else:
                #     cvb_list.append(f"{atom1+1}:-{atom2+1}")
                if not mol.GetBondBetweenAtoms(atom1, atom2):
                    if mol.GetAtomWithIdx(atom1).GetSymbol() == "O" and \
                        mol.GetAtomWithIdx(atom2).GetSymbol() == "O":
                        cvb_list.append(f"{atom1+1}:-{atom2+1}")

        return cvb_list

    
    @staticmethod
    def get_pdb_setcharge(pdb: str, unique_molecules: List[Chem.Mol] = list()) -> list:
        offmols = [Molecule.from_rdkit(mol, allow_undefined_stereo=True) for mol in unique_molecules]
        target_top = Topology.from_pdb(pdb, unique_molecules=offmols)
        setcharge = [f"{i+1}: {a.formal_charge.m_as(unit.elementary_charge)}" 
                    for i, a in enumerate(target_top.atoms)]
        return setcharge
    
    @staticmethod
    def get_mol_setpi(mol: Chem.Mol):
        setpi = []
        no_res_mol = Chem.Mol(mol)
        Chem.Kekulize(no_res_mol)
        for bond in no_res_mol.GetBonds():
            if bond.GetBondType() == Chem.BondType.DOUBLE or \
                bond.GetBondType() == Chem.BondType.TRIPLE:
                setpi.append(f"{bond.GetBeginAtom().GetIdx()+1};{bond.GetEndAtom().GetIdx()+1}")
            if bond.GetBondType() == Chem.BondType.TRIPLE:
                setpi.append(f"{bond.GetBeginAtom().GetIdx()+1};{bond.GetEndAtom().GetIdx()+1}")
        return setpi
    
    def get_config_str(self) -> str:
        return self.get_mopac_mol_str()
    
    def update_config_for(self, some_input: Chem.Mol | str | Tuple[Chem.Mol | str, str], 
                          **kwargs):
        raise NotImplementedError()

class Cuby4MOPACFullConfig(Cuby4MOPACConfig, Cuby4MOPACMolSpecConfig):
    def __init__(self, method: str = "pm6", mopac_exe: str = "auto", 
                 mozyme: bool = True, corrections: str = "d3h4x", 
                 solvent: str = "none", molspec: Cuby4MOPACMolSpecConfig = None):
        Cuby4MOPACConfig.__init__(self, method, mopac_exe, mozyme, corrections, solvent)
        self.molspec = molspec
    
    def get_config_str(self) -> str:
        return self.get_interface_str() + self.molspec.get_mopac_mol_str()
    
    def update_config_for(self, some_input: Chem.Mol | str | Tuple[Chem.Mol | str, str], 
                          unique_mols: List[Chem.Mol] = None, work_dir: str = "tmp", **kwargs):
        if unique_mols is None:
            unique_mols = []

        if isinstance(some_input, Chem.Mol):
            self.molspec = Cuby4MOPACMolSpecConfig.from_rdmol(some_input)
        elif isinstance(some_input, str) and some_input.endswith(".pdb"):
            self.molspec = Cuby4MOPACMolSpecConfig.from_pdb(some_input, unique_mols)
        elif isinstance(some_input, tuple):
            work_dir = str(pathlib.Path(work_dir).absolute())
            ligand = some_input[0]
            target = some_input[1]
            comp_path = os.path.join(work_dir, "cuby_complex_tomolspec_tmp.pdb")
            save_pl_complex(ligand, target, comp_path) # TODO: better implementation of job name
            self.molspec = Cuby4MOPACMolSpecConfig.from_pdb(comp_path, unique_mols + [ligand])
        else:
            raise NotImplementedError()
    
class Cuby4AmberConfig(Cuby4InterfaceConfig):
    def __init__(self, amber_home: str = "auto", topology_file: str = None, 
                 forcefields: List[str] | ForceField = ForceField()):
        super().__init__()
        if amber_home == "auto":
            self.amber_home = str(pathlib.Path(shutil.which("sander")).parent.parent.absolute())
        else:
            self.amber_home = amber_home
        
        if topology_file:
            self.topology_str = "amber_top_file: " + str(pathlib.Path(topology_file).absolute()) 
        else:
            self.topology_str = ""
        
        if isinstance(forcefields, list):
            self.forcefield: ForceField = ForceField(*forcefields)
            self.forcefield_constructor = lambda: ForceField(*forcefields)
        else:
            self.forcefield: ForceField = forcefields
            self.forcefield_constructor = lambda: self.forcefield
    
    def get_amber_config(self):
        return CUBY4_AMBER_MAIN_ADDON.format(
            amber_home=self.amber_home,
            topology_config=self.topology_str
        )
    
    def get_interface_str(self) -> str:
        return self.get_amber_config()
    
    def update_config_for(self, some_input: Chem.Mol | str | Tuple[Chem.Mol | str, str], **kwargs):

        self.forcefield = self.forcefield_constructor()
        if isinstance(some_input, Chem.Mol):
            omm_component = OpenMMComponent.create_component_from_rdmol(some_input, self.forcefield)
            self.topology_str = "amber_top_file: " + omm_component.to_amber_prmtop()
        elif isinstance(some_input, str) and some_input.endswith(".pdb"):
            omm_component = OpenMMComponent.create_component_from_pdb(some_input, self.forcefield)
            self.topology_str = "amber_top_file: " + omm_component.to_amber_prmtop()
        elif isinstance(some_input, tuple):
            ligand = some_input[0]
            target = some_input[1]
            lig_component = OpenMMComponent.create_component_from_rdmol(ligand, self.forcefield, name="ligand")
            tar_component = OpenMMComponent.create_component_from_pdb(target, self.forcefield, name="target")
            complex_component = OpenMMComponent.merge_components([tar_component, lig_component], self.forcefield)
            self.topology_str = "amber_top_file: " + complex_component.to_amber_prmtop()
        else:
            raise NotImplementedError()

class Cuby4QMMMConfig(Cuby4InterfaceConfig):
    def __init__(self, qm_config: Cuby4InterfaceConfig, mm_config: Cuby4InterfaceConfig, qmmm_core: str = ":UNL", 
                 qmmm_embedding: str = "mechanical", grad_on_points: bool = False):
        self.qmmm_core = qmmm_core
        self.qmmm_embedding = qmmm_embedding
        self.grad_on_points = "yes" if grad_on_points else "no"
        self.qm_config = qm_config
        self.mm_config = mm_config
        self.qm_reg_mm_config = deepcopy(mm_config)
    
    def _get_qmmm_qm_str(self):
        return self.qm_config.get_config_str().replace("\n", "\n  ")
    
    def _get_qmmm_mm_str(self):
        return self.mm_config.get_config_str().replace("\n", "\n  ")
    
    def _get_qmmm_qmreg_str(self):
        return self.qm_reg_mm_config.get_config_str().replace("\n", "\n  ")
    
    def get_qmmm_str(self):
        return CUBY4_QMMM_MAIN_ADDON.format(
            qmmm_core=self.qmmm_core,
            qmmm_embedding=self.qmmm_embedding,
            grad_on_points=self.grad_on_points,
            qm_config=self._get_qmmm_qm_str(),
            mm_config=self._get_qmmm_mm_str(),
            qmregion_mm_config=self._get_qmmm_qmreg_str()
        ) 
    
    def get_interface_str(self) -> str:
        return self.get_qmmm_str()
    
    def update_config_for(self, some_input: Tuple[Chem.Mol | str, str], **kwargs):

        if isinstance(some_input, tuple):
            ligand = some_input[0]
            target = some_input[1]

            self.mm_config.update_config_for((ligand, target))
            self.qm_reg_mm_config.update_config_for(ligand)
            self.qm_config.update_config_for(ligand)
        else:
            raise NotImplementedError()

class Cuby4FullConfig(Cuby4Config):
    def __init__(self, job_config: Cuby4JobConfig, interface_config: Cuby4InterfaceConfig):
        self.job_config = job_config
        self.interface_config = interface_config
    
    def get_config_str(self):
        job_str = self.job_config.get_config_str()
        interface_str = self.interface_config.get_config_str()
        return job_str + interface_str
    
    def update_config_for(self, some_input: Chem.Mol | str, **kwargs):
        raise NotImplementedError()
    