import argparse
import subprocess
import sys
import os

import rdkit
from rdkit import Chem
from rdkit.Chem import rdDistGeom, rdMolDescriptors

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(
        prog="Docking Performer",
    )
    parser.add_argument("-l", "--ligands", help="path to ligands file (.sdf)", required=True)
    parser.add_argument("-p", "--protein", help="path to protein file (.pdb)", required=True)
    parser.add_argument("-o", "--output", help="path to output file (.sdf)", required=True)
    parser.add_argument("--work_dir", help="where temporary files are stored", default="tmp")
    parser.add_argument("--gnina_path", help="path to gnina, in the future might add others")
    parser.add_argument("--center")
    parser.add_argument("--bsize")
    parser.add_argument("--flexres")

    args = parser.parse_args()

    center = args.center.split(",")
    bsize = args.bsize.split(",")

    cmd_line = f"{args.gnina_path} -r {args.protein} -l {args.ligands} --center_x {center[0]} \
        --center_y {center[1]} --center_z {center[2]} --size_x {bsize[0]} --size_y {bsize[1]} \
            --size_z {bsize[2]} -o {args.output} --log gnina.log"
    if args.flexres:
        cmd_line += f" --flexres {args.flexres}"
    print(cmd_line)
    subprocess.run(cmd_line.split())
    