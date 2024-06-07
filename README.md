# Introduction
...

# Installation & Setup
1) Clone this repository
```terminal
git clone https://github.com/arazthexd/sqm-vscreen
cd sqm-vscreen
```
2) Create a new conda/mamba environment using the `environment.yml` file.
```terminal
mamba env create -f environment.yml
mamba activate sqmvscreen
```
3) You will also need a GNINA binary. You can get it from [here](https://github.com/gnina/gnina/releases)
```terminal
wget https://github.com/gnina/gnina/releases/download/v1.1/gnina -P .. # Or any place you would like to store the binary
```
Please note that the above link may not be the latest release of gnina in the future. Always check from their repo for the latest release.

4) Next you need to install [QupKake](https://github.com/Shualdon/QupKake) for pKa predictions used in our code.
```terminal
cd .. # or any other folder you want to put QupKake in
git clone https://github.com/Shualdon/QupKake
cd QupKake
pip install .
```


# Usage
Here is an example case for the sake of reproducibility (should be updated for sure)
This starts when in the main folder of this repo.
1) Preparation of Ligands...
```terminal
python components/0_prepare/prepare_ligands.py data/ligands/tmp/test_ligands.smi output/prepared_ligands.sdf --gypsum_dl_dir ~/tools/gypsum_dl/ --num_proc 8
```
2) Preparation of Protein...
```terminal
export PYTHONPATH=.
python components/0_prepare/prepare_protein.py data/targets/kguD.pdb output/prepared_protein.pdb --chains all --work_dir tmp
```
3) Docking...
```terminal
python components/1_vscreen/dock_ligands.py -l output/prepared_ligands.sdf -p output/prepared_protein.pdb -o output/docked_ligands.sdf --gnina_path ~/tools/gnina/gnina --center 1.16,-0.56,-2.91 --bsize 20,20,20
```
4) Candidate Selection
```terminal
python components/1_vscreen/select_candidates.py output/docked_ligands.sdf output/selected_ligands.sdf
```

# Reference
If you use the code in this repo, please consider citing the below.
- SQMVScreen (This Repo):
```
...
```
- QupKake:
```
@article{qupkake, 
    title={QupKake: Integrating Machine Learning and Quantum Chemistry for micro-pKa Predictions}, 
    DOI={10.26434/chemrxiv-2023-gxplb}, 
    journal={ChemRxiv}, 
    publisher={Cambridge Open Engage}, 
    author={Abarbanel, Omri and Hutchison, Geoffrey}, 
    year={2023}
}
```
