import argparse
import subprocess
import sys
import os

from ezdd.utils.protein.prepare import clean_protein

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(
        prog="Protein Preparer",
    )
    parser.add_argument("filename", help="path to input file (.pdb)")
    parser.add_argument("outfile", help="path to output file (.pdb)")
    parser.add_argument("-v", "--verbose", action="store_true")
    parser.add_argument("--work_dir", help="where temporary files are stored", default="tmp")
    parser.add_argument("--chains", help="chains of protein to keep", default="all")
    
    args = parser.parse_args()

    clean_protein(args.filename, args.chains, args.outfile)