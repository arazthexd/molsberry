from pdbfixer import PDBFixer
from openmm.app import PDBFile

from ezdd.utils.protein.protio import read_prot_file

def clean_protein(protein_input, sel_chains: str = "all", out_path=None):
    """
    Cleans the protein structure using PDBFixer.

    Adds missing atoms or residues, removes heterogens and solvent,
    adds hydrogens, and keeps only the specified chains.
    """
    fixer: PDBFixer = read_prot_file(protein_input, "pdbfixer")
    if sel_chains != "all":
        rem_chains = [chain.id for chain in fixer.topology.chains() if chain.id not in sel_chains]
        fixer.removeChains(chainIds=rem_chains)
    fixer.removeHeterogens(False)
    fixer.findMissingResidues()
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.0)  # TODO: Any need to pH value?
    if out_path:
        PDBFile.writeFile(fixer.topology, fixer.positions, open(out_path, 'w')) 
        # TODO: This overwrites the previous file, right? right?!
    return fixer


