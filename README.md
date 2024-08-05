# Introduction
**Mo**dular **L**igand/**S**tructure-**b**ased Drug Discov**er(r)y** Framework in Python

# To-Do
## Core
- [x] Add Analyzer template to the core.
- [x] Add some ligand/protein operators that take a context such as pocket.
- [ ] Create a way of generating list of used project citations (Important, Not Urgent)
- [ ] Add input block for inputs of different names...

## New Modules
- [ ] Cuby4 (On Hold...)
- [ ] OpenMM (Low Priority)
- [x] Independent MOPAC (Urgent)
- [ ] Test + Decide on Adding ASH For QMMM
- [ ] Check out QMMM on MOPAC

## Module Updates
- [x] Add example rdkit analyzer to rdkit module.
- [ ] (Urgent) Add all needed nodes for MOPAC (from different core templates)
- [ ] Write tests for mopac.

## Testing
- [x] Add testing for entire pipelines (just core for now)
- [ ] Start testing for simple example pipelines using rdkit and mopac nodes as well as core.

## Docs
- [ ] Update installation documentation as well as usage.


NOTE: CUBY4 HAS LOTS OF PROBLEMS. BEWARE!

# Installation & Setup
1) Create an environment to install required and optional packages in. 

Using `conda` (preferred):
```terminal
TODO
```
Using `pip`:
```terminal
TODO
```

2) Install requirements for the core package:
TODO

3) Install requirements for each module you're planning to use:
TODO

4) Change global configuration to your liking:
TODO

5) Perform testing to make sure everything is right. 
```terminal
pytest
```

# Use
The framework is designed to be used as a molecular or drug discovery pipeline creator. All is needed to use it is to create a couple of nodes (blocks) in your pipeline and connect their inputs and outputs. You will also need to add an `OutputBlock` if you want your pipeline to output some info generated during the run. At last, you can simply create an instance of your framework, run it with your inputs and the output will be returned. 

Checkout the following example to get a feel of how it all works:
```python
TODO
```

**NOTE**: The framework is designed to be modular. You can create your own modules and blocks or use the modules written by others. If you have ideas for useful modules for different related interfaces, we appreciate and encourage you to contribute and submit a pull request for it. 

## Creating New Blocks
TODO

## Creating New Modules
TODO