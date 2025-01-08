# Introduction
NOTE: This repository is a work in progress and it's far from being usable. At the moment, it's a mini project that hopefully will turn into something useful for our team and anyone else who could benefit from it.
Any feedbacks, collaborations and contributions are absolutely welcome!
:D

**Mo**dular **L**igand and **S**tructure-**b**ased Drug Discov**er(r)y** Framework in Python

# To-Do
## Critical
- [ ] ASH Optimizer --> geomeTRIC optimization converge error

## Short-term
- [ ] Add generic blocks such as mathematics, etc.
- [ ] Add support for MM scoring / optimization. (OpenMM / Ambertools)
- [ ] Add support for QMMM scoring / optimization. (ASH / Cuby4 / Ambertools)
- [ ] Add testing for entire pipelines.


## Long-term
- [ ] Create a way of generating list of used project citations
- [ ] Update installation documentation as well as usage.
- [ ] Improve parallel computing.
- [ ] Create extensive introductory notebooks.


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
