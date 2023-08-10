# PAMpy
Python implementation of the PAM framework (https://github.com/Spherotob/PAM_public) to create GECKO like ME models.

## Code structure:
- **example_PAM_generation**: (in `./Examples/PAModel_example_script.ipynb`) An example of how to build, run and validate a proteome allocation model for *Escherichia coli*.
- **EnzymeSectors**: The objects which are used to store the data of the different enzyme sectors which are added to the genome-scale model
- **PAModel**: Proteome Allocation (PA) model class. This class builds on to the `cobra.core.Model` class from the COBRApy toolbox with functions to build enzyme sectors, to add enzyme kinetics parameters and in the future to perform a sensitivity analysis on the enzyme variables.
- **Enzyme**: Different classes which relate enzymes to the model with enzyme constraints and variables.
- **CatalyticEvent**: A class which serves as an interface between reactions and enzyme. This allows for easy lookup of Protein-Reaction assocations.
- **PAMValidator**: Functions to validate the model predictions with physiology data and giving a graphical overview. The script uses data for E.coli (found in `./Data/Ecoli_physiology`) by default.

## Dependencies
cobrapy toolbox, Gurobi solver

More dependencies see `pyproject.toml`
