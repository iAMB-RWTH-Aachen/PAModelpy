# Setting up PAMs using PAModelpy

This guide provides step-by-step instructions on how to use the provided Python scripts to build Protein Allocation Models (PAMs) from a genome-scale model and parameter datasets.

---

## Overview

The `pam_generation.py` script allows users to set up a PAM based on a genome-scale metabolic model and an enzyme database. It integrates genome-protein-reaction (GPR) rules with enzyme kinetics to partition cellular proteins into functional sectors.

The accompanying `test_pam_generation.py` provides unit tests to validate the setup process.

---

## Prerequisites

### Required python libraries

Ensure you have `PAModelpy` installed

Install these dependencies via pip:

```bash
pip install cobra PAModelpy
```

### Input files

1. **Genome-scale model**: the path to metabolic model in SBML format (e.g., `iML1515.xml`) or a cobra.Model instance (e.g., build from a json file like `e_coli_core.json`).
2. **Parameter Excel file**: Contains data on enzymatic properties. The file should have at least the following sheets:
   - **ActiveEnzymes**: Enzyme-specific parameters such as reaction IDs, kcat values, and molar masses.
   - **Translational** (optional): Parameters related to the translational protein fraction.
   - **UnusedEnzyme** (optional): Parameters for the unused enzyme sector.

---

## Dataset requirements

### Excel file structure

#### ActiveEnzymes sheet
| **Column**         | **Description**                                         |
|---------------------|---------------------------------------------------------|
| `rxn_id`           | Reaction ID from the genome-scale model.                |
| `enzyme_id`        | Unique enzyme identifier*.                              |
| `gene`             | List of genes associated with the enzyme.               |
| `GPR`              | Gene-protein-reaction association (logical expression). |
| `molMass`          | Molar mass of the enzyme (kDa).                         |
| `kcat_values`      | Turnover number (s⁻¹).                                  |
| `direction`        | Reaction direction (`f` for forward, `b` for backward). |

* A unique enzyme identifier is defined as a single identifier per catalytically active unit. This means that an enzyme complex has a single identifier.
* Enzyme-complex identifiers can be parsed from peptide identifiers and gene-to-protein mapping using the `merge_enzyme_complexes` function.

#### Translational sheet
| **Parameter**      | **Description**                                                               |
|---------------------|-------------------------------------------------------------------------------|
| `id_list`          | Identifier related to protein fraction associated with translational proteins |
| `tps_0`            | Translational protein fraction at zero growth rate. [g_protein/g_CDW]         |
| `tps_mu`           | Change in translational protein fraction per unit change of the associated reaction.  |
| `mol_mass`         | Molar mass of the translational enzymes. [kDa]                                |

#### UnusedEnzyme sheet
| **Parameter**      | **Description**                                                                 |
|---------------------|---------------------------------------------------------------------------------|
| `id_list`          | Identifier related to protein fraction associated with the unused enzyme sector |
| `ups_0`            | Unused enzyme fraction at zero growth rate. [g_protein/g_CDW]                   |
| `ups_mu`           | Change in unused enzyme fraction per unit change of the associated reaction.    |
| `mol_mass`         | Molar mass of unused enzymes.  [kDa]                                            |

---

## Building a PAM

### Steps to create a PAM

1. Prepare the genome-scale model and parameter file.
2. Optional: change defaults in the Config object to match your model identifiers
3. Use the `set_up_pam` function to initialize the model.
4. Optimize the PAM

#### Example usage for *E. coli*
The model defaults are set to generate a PAM for the iML1515 model of *Escherichia coli* K-12.
```python
from Scripts.pam_generation import set_up_pam

#1. Define input paths
model_path = "Models/iML1515.xml"
param_file = "Data/proteinAllocationModel_iML1515_EnzymaticData_new.xlsx"
#2. Config is not required for iML1515
#3. Build the PAM
pam = set_up_pam(pam_info_file=param_file,
                 model=model_path,
                 total_protein=0.258,  # Optional: Total protein concentration (g_prot/g_cdw)
                 active_enzymes=True,
                 translational_enzymes=True,
                 unused_enzymes=True,
                 sensitivity=True,
                 adjust_reaction_ids=False)
#4. Optimize to find the max growth rate at a substrate uptake rate of -10 mmol/gcdw/h
pam.optimize()
print(f"Objective Value: {pam.objective.value}")
```

#### Example usage for other microorganisms
The model defaults have to be adapted to the identifiers for another microorganism. In case of a model in the BiGG namespace,
only the biomass reaction has to be adapted (which is the case for e.g. [iJN1463](https://onlinelibrary.wiley.com/doi/abs/10.1111/1462-2920.14843), 
the genome-scale model for *Pseudomonas putida* KT2440). For some other models, the entire namespace has to be adapted (e.g. for the [Yeast9](https://www.embopress.org/doi/full/10.1038/s44320-024-00060-7) 
model of *Saccharomyces cerevisiae*). This can be accomplished using a [`Config`](api_reference/configuration.md) object. 

##### Example for BiGG namespace: *P. putida* iJN1463
Please note that this is only an example, both the model, as the parameter file do not exist in this repository.

```python
from Scripts.pam_generation import set_up_pam
from PAModelpy import Config

#1. Define input paths
model_path = "Models/iJN1463.xml"
param_file = "Data/proteinAllocationModel_iJN1463_EnzymaticData.xlsx"

#2. Change biomass reaction id in config
config = Config()
config.BIOMASS_REACTION = 'BIOMASS_KT2440_WT3'

#3. Build the PAM
pam = set_up_pam(pam_info_file=param_file,
                 model=model_path,
                 config=config,
                 total_protein=0.258,  # Optional: Total protein concentration (g_prot/g_cdw)
                 active_enzymes=True,
                 translational_enzymes=True,
                 unused_enzymes=True,
                 sensitivity=True,
                 adjust_reaction_ids=False)
#4. Optimize to find the max growth rate at a substrate uptake rate of -10 mmol/gcdw/h
pam.optimize()
print(f"Objective Value: {pam.objective.value}")
```

##### Example for non-BiGG namespace: *S. cerevisia* Yeast9
Please note that this is only an example, both the model, as the parameter file do not exist in this repository.

```python
from Scripts.pam_generation import set_up_pam
from PAModelpy import Config

#1. Define input paths
model_path = "Models/Yeast9.xml"
param_file = "Data/proteinAllocationModel_yeast9_EnzymaticData.xlsx"

#2. Change all the reaction ids in config and the protein regex
config = Config()
    config.TOTAL_PROTEIN_CONSTRAINT_ID = "TotalProteinConstraint"
    config.P_TOT_DEFAULT = 0.388  # g_protein/g_cdw
    config.CO2_EXHANGE_RXNID = "r_1672"
    config.GLUCOSE_EXCHANGE_RXNID = "r_1714"
    config.BIOMASS_REACTION = "r_2111"
    config.OXYGEN_UPTAKE_RXNID = "r_1992"
    config.ACETATE_EXCRETION_RXNID = "r_1634"
    config.PHYS_RXN_IDS = [
    config.BIOMASS_REACTION,
    config.GLUCOSE_EXCHANGE_RXNID,
    config.ACETATE_EXCRETION_RXNID,
    config.CO2_EXHANGE_RXNID,
    config.OXYGEN_UPTAKE_RXNID]
    config.ENZYME_ID_REGEX = r'(Y[A-P][LR][0-9]{3}[CW])'

#3. Build the PAM
pam = set_up_pam(pam_info_file=param_file,
                 model=model_path,
                 config=config,
                 total_protein=0.258,  # Optional: Total protein concentration (g_prot/g_cdw)
                 active_enzymes=True,
                 translational_enzymes=True,
                 unused_enzymes=True,
                 sensitivity=True,
                 adjust_reaction_ids=False)
#4. Optimize to find the max growth rate at a substrate uptake rate of -10 mmol/gcdw/h
pam.optimize()
print(f"Objective Value: {pam.objective.value}")
```

##### A short note on enzyme identifiers
The Config object has an entry which enables the framework to find and extract enzyme identifiers from the PAModel:
```python
Config.ENZYME_ID_REGEX = r'(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})'
```
This is the default [regular expression](https://www.uniprot.org/help/accession_numbers) to find UniProt identifiers, 
as derived from the UniProt database (obtained 2024-08-07). Two other default placeholder identifiers (`E1`, `E10`, `Enzyme_GLC_D`, etx) 
are always included in any regex search with the PAModel:

```python
default_enzyme_regex = r'E[0-9][0-9]*|Enzyme_*'
```
In case you would like to use other placeholders or another form of protein identifiers, please adapt the `Config.ENZYME_ID_REGEX`
attribute with the proper regular expression.

---

## Testing the setup

To verify the correctness of the PAM generation process, run the tests provided in `tests/unit_tests/test_utils/test_pam_generation.py`:

```bash
python -m pytest tests/unit_tests/test_utils/test_pam_generation.py.py
```

---

## Additional features

### Increasing kcat values
The script includes a utility to scale up kcat values in the parameter file:

```python
from PAModelpy.utils.pam_generation import increase_kcats_in_parameter_file

increase_kcats_in_parameter_file(
    kcat_increase_factor=2,
    pam_info_file_path_ori="Data/old_param_file.xlsx",
    pam_info_file_path_out="Data/new_param_file.xlsx"
)
```
### Generate enzyme-complex identifiers from peptide ids (e.g. uniprot annotation)
This helps with setting up a novel parameter file using identifiers obtained from uniprot and a gene-to-protein mapping

```python
from PAModelpy.utils.pam_generation import get_protein_gene_mapping, merge_enzyme_complexes
from cobra.io import read_sbml_model

model = read_sbml_model('Models/iML1515.xml')
enzyme_db = "Data/old_param_file.xlsx"

protein2gene, gene2protein = get_protein_gene_mapping(enzyme_db, model)
# Ensure the enzyme complexes are merged on one row
eco_enzymes_mapped = merge_enzyme_complexes(enzyme_db, gene2protein)
```
---

## Troubleshooting

- **Issue: Missing reaction in enzyme database**  
  Solution: Ensure all reactions in the model are represented in the parameter file or adjust `reaction_ids` using `adjust_reaction_ids=True`.

- **Issue: Objective value is zero after optimization**  
  Solution: Check the input parameter file for consistency and ensure that all reactions are correctly annotated.
