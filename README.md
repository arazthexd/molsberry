# Introduction
...

# Installation & Setup
1) Create a new conda/mamba environment using the `environment.yml` file.
```terminal
mamba env create -f environment.yml
mamba activate sqmvscreen
```
2) Next you need to install [QupKake](https://github.com/Shualdon/QupKake) for pKa predictions used in our code.
```terminal
cd .. # or any other folder you want to put QupKake in
git clone https://github.com/Shualdon/QupKake
cd QupKake
pip install .
```

# Usage
...

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
