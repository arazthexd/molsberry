import argparse

from rdkit import Chem

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        prog="Docked Ligand Selector",
    )
    parser.add_argument("filename", help="path to docked ligands file (.sdf)")
    parser.add_argument("outfile", help="path to selected docked ligands file (.sdf)")
    args = parser.parse_args()

    suppl = Chem.SDMolSupplier(args.filename, removeHs=False)
    for i, mol in enumerate(suppl):
        