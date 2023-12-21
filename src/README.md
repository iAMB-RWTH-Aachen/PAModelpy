# PAModelpy - Protein Allocation Model reconstruction in Python

## What is PAModelpy?
Models of metabolism are powerful tools to explore the metabolic potential of microorganism. 
Powerful tools have been created to support Python-based analysis of genome-scale models. 
These models, however, cannot capture all metabolic phenotypes and the simulation results have high flux variability.
Adding protein to each reaction increases the simulation fidelity.
The PAModelpy package is designed to integrate protein constraints and protein sectors as described by
[Alter et al. (2021)](https://journals.asm.org/doi/10.1128/mSystems.00625-20) to metabolic models.
It is the Python implementation of the [PAM MATLAB framework](https://github.com/Spherotob/PAM_public) to create GECKO like ME models.

The PAModelpy package builds upon the community-wide used [COBRApy](https://github.com/opencobra/cobrapy/tree/devel). 
We have extended this package with the following features:
- protein-reaction associations
- infrastructure to include isozymes and promiscuous enzymes
- protein sectors
- specialized objects to build protein allocation models
- the possibility to perform a computational efficient sensitivity analysis

## Installation
[PAModelpy is a PiPy package](https://pypi.org/project/PAModelpy/) which allows for easy installation with pip:

`pip install PAModelpy`

Note that the package has been tested with the [Gurobi](https://www.mathworks.com/products/connections/product_detail/gurobi-optimizer.html) solver.

## Code structure:
- **EnzymeSectors**: The objects which are used to store the data of the different enzyme sectors which are added to the genome-scale model
- **PAModel**: Proteome Allocation (PA) model class. This class builds on to the `cobra.core.Model` class from the COBRApy toolbox with functions to build enzyme sectors, to add enzyme kinetics parameters and in the future to perform a sensitivity analysis on the enzyme variables.
- **Enzyme**: Different classes which relate enzymes to the model with enzyme constraints and variables.
- **CatalyticEvent**: A class which serves as an interface between reactions and enzyme. This allows for easy lookup of Protein-Reaction assocations.
- **PAMValidator**: Functions to validate the model predictions with physiology data and giving a graphical overview. The script uses data for E.coli (found in `./Data/Ecoli_physiology`) by default.

## Dependencies
The dependencies of the PAModelpy package can be found in `src/pyproject.toml`

## License
Copyright institute of Applied Microbiology, RWTH Aachen University, Aachen, Germany (2023)

PAModelpy is free of charge open source software, which can be used and modified for your particular purpose under the [MIT](https://opensource.org/license/mit/)
or [Apache 2.0](https://www.apache.org/licenses/LICENSE-2.0) of the users choice.

Please note that according to these licenses, the software is provided 'as is', WITHOUT WARRANTY OF ANY KIND, without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.