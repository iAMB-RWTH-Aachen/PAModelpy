# PAMpy - Protein Allocation Model reconstruction in Python

## What is PAMpy?
Models of metabolism are powerful tools to explore the metabolic potential of microorganism. 
Powerful tools have been created to support Python-based analysis of genome-scale models. 
These models, however, cannot capture all metabolic phenotypes and the simulation results have high flux variability.
Adding protein to each reaction increases the simulation fidelity.
The PAModelpy package is designed to integrate protein constraints and protein sectors as described by [Alter et al. (2021)](https://journals.asm.org/doi/10.1128/mSystems.00625-20) to metabolic models.
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

## What can you find where in this repository?
This repository contains not only the source code, but also examples and scripts which were used in [publication](linktopublication).
- **Data**
  - *eGFP_expression_Bienick2014*: measured growth rate and eGFP expression by [Bienick et al. (2014)](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0109105)
  - *proteinAllocationModel_iML1515_EnzymaticData_py*: information about the proteinsectors of the PAM for *Escherichia coli* (*E.coli*)
  - *proteome_data_extract_schmidt2016*: quantitative proteomics data from [Schmidt et al. (2016)](https://www.nature.com/articles/nbt.3418) used to parametrize the *E.coli* core PAM
  - *Ecoli_phenotypes/Ecoli_phenotypes_py_rev*: experimental physiology measurements to validate the model simulations
- **Examples**: example notebook on how to build, run and validate a PAM using the PAModelpy package
- **Figures**: scripts used to create Figure 2-3 and supplementary figures
- **Models**: models used (iML1515 and core ecoli model)
- **Results**: results of computational performance analysis
- **Scripts**: scripts used for gathering results
  - computational performance analysis: `compare_computational_efficiency_fac.py` and `numeric_error_estimation_schemes_fac.py`
  - *E.coli* core PAM creation: `analyze_proteome.ipynb` and `create_ecolicore_pam_incl_UE.ipynb`
  - Sensitivity analysis: `toy_ec_pam.py` and `Ecoli_core_sensitivity_analysis.ipynb`
- **src/PAModelpy**: source code for PAModelpy package


## Code structure:
- **EnzymeSectors**: The objects which are used to store the data of the different enzyme sectors which are added to the genome-scale model
- **PAModel**: Proteome Allocation (PA) model class. This class builds on to the `cobra.core.Model` class from the COBRApy toolbox with functions to build enzyme sectors, to add enzyme kinetics parameters and in the future to perform a sensitivity analysis on the enzyme variables.
- **Enzyme**: Different classes which relate enzymes to the model with enzyme constraints and variables.
- **CatalyticEvent**: A class which serves as an interface between reactions and enzyme. This allows for easy lookup of Protein-Reaction assocations.
- **PAMValidator**: Functions to validate the model predictions with physiology data and giving a graphical overview. The script uses data for E.coli (found in `./Data/Ecoli_physiology`) by default.

## Dependencies
Dependencies for the scripts in this repository, not included in the PAModelpy package:
- `PAModelpy`
- `plotly`
- `matplotlib`
- `scipy`
- `time`
- `resource`
- `PIL`

All dependencies can be installed in one go by downloading this repository and running:

`python setup.py install`

from the `PAMpy` directory

The dependencies of the PAModelpy package can be found in `src/pyproject.toml`

## License
Copyright institute of Applied Microbiology, RWTH Aachen University, Aachen, Germany (2023)

PAModelpy is free of charge open source software, which can be used and modified for your particular purpose under the [MIT](https://opensource.org/license/mit/)
or [Apache 2.0](https://www.apache.org/licenses/LICENSE-2.0) of the users choice.

Please note that according to these licenses, the software is provided 'as is', WITHOUT WARRANTY OF ANY KIND, without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.