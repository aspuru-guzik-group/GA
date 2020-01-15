# Augmenting genetic algorithms with deep neural networks for exploring the chemical space
This repository contains code for the paper: [Augmenting genetic algorithms with deep neural networks for exploring the chemical space](https://arxiv.org/abs/1909.11655).
Here is a visualization of molecular progress: 

<img align="center" src="./readme_docs/mol_view.gif"/>


# Directory Navigator: 
The code for this repository is arranged based on the experiments of the paper. Particularly: 
- [Experiment 4.1: ](https://github.com/akshat998/GA/tree/master/4.1) Unconstrained optimization and comparison with other generative models
- [Experiment 4.2: ](https://github.com/akshat998/GA/tree/master/4.2) Long term experiment with a time-dependent adaptive penalty
- [Experiment 4.3: ](https://github.com/akshat998/GA/tree/master/4.3) Analysis of molecule classes explored by the GA
- [Experiment 4.4: ](https://github.com/akshat998/GA/tree/master/4.4) Constrained optimization
- [Experiment 4.5: ](https://github.com/akshat998/GA/tree/master/4.5) Simultaneous logP and QED optimization
- [Experiment 4.6: ](https://github.com/akshat998/GA/tree/master/4.6) Modification of the hyperparameter beta

Instructions on running specific experiments are provided in the above links. Please note that the code has been parallelized based on the number of CPU cores for quick property evaluations.

Due to the large size of the repository, we have created a seperate branch that contains outputs from all the eperiment. For this option, please run (note: this is a 4GB branch, and needs 20mins of cloning time): 
```
git clone --single-branch --branch paper_results https://github.com/akshat998/GA.git
```

# Prerequisites

Before running the code, please ensure you have the following:


- [SELFIES (any version)](https://github.com/aspuru-guzik-group/selfies) - 
  The code was run with v0.1.1 (which is the fastest), however, the code is compatible with any version. 
- [RDKit](https://www.rdkit.org/docs/Install.html)
- [tensorboardX](https://pypi.org/project/tensorboardX/)
- [Pytorch v0.4.1](https://pytorch.org/)
- [Python 3.0 or up](https://www.python.org/download/releases/3.0/)
- [numpy](https://pypi.org/project/numpy/)



# Questions, problems?

Make a github issue 😄. Please be as clear and descriptive as possible. Please feel free to reach
out in person: (akshat[DOT]nigam[AT]mail[DOT]utoronto[DOT]ca)

