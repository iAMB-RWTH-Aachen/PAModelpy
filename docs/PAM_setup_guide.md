# Setting up PAMs using PAModelpy

This guide provides step-by-step instructions on how to use the provided Python scripts to build Protein Allocation Models (PAMs) from a genome-scale model and parameter datasets.

---

## Overview

The `pam_generation.py` script allows users to set up a PAM based on a genome-scale metabolic model and an enzyme database. It integrates genome-protein-reaction (GPR) rules with enzyme kinetics to partition cellular proteins into functional sectors.

The accompanying `test_pam_generation.py` provides unit tests to validate the setup process.

---

## Prerequisites

### Required Python Libraries

Ensure you have `PAModelpy` installed

Install these dependencies via pip:

```bash
pip install cobra PAModelpy
```

### Input Files

1. **Genome-scale model**: the path to metabolic model in SBML format (e.g., `iML1515.xml`) or a cobra.Model instance (e.g., build from a json file like `e_coli_core.json`).
2. **Parameter Excel file**: Contains data on enzymatic properties. The file should have at least the following sheets:
   - **ActiveEnzymes**: Enzyme-specific parameters such as reaction IDs, kcat values, and molar masses.
   - **Translational** (optional): Parameters related to the translational protein fraction.
   - **UnusedEnzyme** (optional): Parameters for the unused enzyme sector.

---

## Dataset Requirements

### Excel File Structure

#### ActiveEnzymes Sheet
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

#### Translational Sheet
| **Parameter**      | **Description**                                                               |
|---------------------|-------------------------------------------------------------------------------|
| `id_list`          | Identifier related to protein fraction associated with translational proteins |
| `tps_0`            | Translational protein fraction at zero growth rate. [g_protein/g_CDW]         |
| `tps_mu`           | Change in translational protein fraction per unit change of the associated reaction.  |
| `mol_mass`         | Molar mass of the translational enzymes. [kDa]                                |

#### UnusedEnzyme Sheet
| **Parameter**      | **Description**                                                                 |
|---------------------|---------------------------------------------------------------------------------|
| `id_list`          | Identifier related to protein fraction associated with the unused enzyme sector |
| `ups_0`            | Unused enzyme fraction at zero growth rate. [g_protein/g_CDW]                   |
| `ups_mu`           | Change in unused enzyme fraction per unit change of the associated reaction.    |
| `mol_mass`         | Molar mass of unused enzymes.  [kDa]                                            |

---

## Building a PAM

### Steps to Create a PAM

1. Prepare the genome-scale model and parameter file.
2. Use the `set_up_pam` function to initialize the model.

#### Example Usage
```python
from Scripts.pam_generation import set_up_pam

# Define input paths
model_path = "Models/iML1515.xml"
param_file = "Data/proteinAllocationModel_iML1515_EnzymaticData.xlsx"

# Build the PAM
pam = set_up_pam(pam_info_file=param_file,
                 model=model_path,
                 total_protein=0.258,  # Optional: Total protein concentration (g_prot/g_cdw)
                 active_enzymes=True,
                 translational_enzymes=True,
                 unused_enzymes=True,
                 sensitivity=True,
                 adjust_reaction_ids=False)
```

3. Optimize the PAM:
```python
pam.optimize()
print(f"Objective Value: {pam.objective.value}")
```

---

## Testing the Setup

To verify the correctness of the PAM generation process, run the tests provided in `test_pam_generation.py`:

```bash
python -m pytest test_pam_generation.py
```

---

## Additional Features

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
