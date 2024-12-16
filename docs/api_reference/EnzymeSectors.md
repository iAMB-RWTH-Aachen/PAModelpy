# EnzymeSectors

## Sector Objects

```python
class Sector(Object)
```

#### \_\_copy\_\_

```python
def __copy__() -> "Sector"
```

Copy the Sector.

**Returns**:

- `Sector` - A new Sector that is a copy of the original Sector.

#### \_\_deepcopy\_\_

```python
def __deepcopy__(memo: dict) -> "Sector"
```

Copy the Sector with memo.

**Arguments**:

- `memo` _dict_ - Automatically passed parameter.
  

**Returns**:

- `Sector` - A new Sector that is a copy of the original Sector with memo.

## EnzymeSector Objects

```python
class EnzymeSector(Sector)
```

#### DEFAULT\_MOL\_MASS

mean enzymes mass E.coli [g/mol]

## ActiveEnzymeSector Objects

```python
class ActiveEnzymeSector(Sector)
```

#### DEFAULT\_MOL\_MASS

mean enzymes mass E.coli [g/mol]

#### \_\_init\_\_

```python
def __init__(rxn2protein: dict,
             protein2gene: dict = {},
             configuration: Config = Config)
```

_summary_

**Arguments**:

- `rxn2protein` _dict_ - A dictionary with reaction ID, enzymes_dict key, value pairs for each reaction in the active_enzyme sector.
  The enzymes_dict contains the enzyme identifiers of the enzymes related to the specific reaction as keys, and a dict
  with information about the enzyme as values. The information included in this dict includes the turnover number for
  the forward and backward reaction (1/s), molar mass of the enzyme (mol/g), gene identifiers related to the enzyme,
  and with which other enzymes it forms a complex.
- `protein2gene` - dict
  enzyme_id, gene_list key, value pairs for each enzyme.The gene_list value is a list of lists which indicates
  &#x27;and&#x27; or &#x27;or&#x27; relationships between the genes which code for the enzyme(complex).
  

**Example**:

    ```
    For the Parameter rxn2protein a dictionary may look like this:
    {
        'R1':
            {E1:
                        {'f': forward kcat, 'b': backward kcat, 'molmass': molar mass, 'genes': [G1, G2],
                        'protein_reaction_association': [[E1, E2]]
                        },
                    E2:
                        {'f': forward kcat, 'b': backward kcat, 'molmass': molar mass, 'genes': [G3, G4],
                         'protein_reaction_association': [[E1, E2]]
            }
    }

    For the Parameter protein2gene a dictionary may look like this:
    {E1:
            [[gene1], [gene2, gene3]],
    E2:
            [[gene4]]
    }
    where the gene-protein-reaction associations are the following:
    E1: gene1 or (gene2 and gene3)
    E2: gene4
    ```

**Arguments**:

  rxn2protein : nested dict
  reaction id, enzymes_dict key, value pairs for each reaction in the active_enzyme sector.
  The enzymes_dict contains the enzyme identifiers of the enzymes which are related to the specific reaction
  as keys and a dict with information about the enzyme as values. The information included in this dict is:
  turnover number for forward and backward reaction [1/s], molar mass of the enzyme [mol/g], gene identifiers
  related to the enzyme, with which other enzymes it forms a complex.
  
- `protein2gene` - dict
  enzyme_id, gene_list key, value pairs for each enzyme.The gene_list value is a list of lists which indicates
  &#x27;and&#x27; or &#x27;or&#x27; relationships between the genes which code for the enzyme(complex).
  
- `configuration` - Config object, optional
  Information about general configuration of the model including identifier conventions.
  Default as defined in the `PAModelpy.configuration` script for the E.coli iML1515 model.

#### check\_kcat\_values

```python
def check_kcat_values(model, reaction, enzyme_dict)
```

Function to check if the kcat values provided for an enzyme are consistent with the direction of the reaction

**Arguments**:

- `model` _cobra.Model or PAModel_ - Model to which the kcat values should be added.
- `reaction` _cobra.Reaction_ - Reaction that is catalyzed with the enzyme related to the kcats.
- `enzyme_dict` _dict_ - A dictionary with the turnover values for the forward and/or backward reaction for different enzymes [/s].
  

**Example**:

  Example dictionary for the `enzyme_dict` parameter
    ```
    {'E1': {'f': forward kcat, 'b': backward kcat}}
    ```

## TransEnzymeSector Objects

```python
class TransEnzymeSector(EnzymeSector)
```

#### DEFAULT\_MOL\_MASS

default E. coli ribosome molar mass [g/mol]

## UnusedEnzymeSector Objects

```python
class UnusedEnzymeSector(EnzymeSector)
```

#### DEFAULT\_MOL\_MASS

mean enzymes mass E.coli [g/mol]

## CustomSector Objects

```python
class CustomSector(EnzymeSector)
```

#### DEFAULT\_ENZYME\_MOL\_MASS

mean enzymes mass E.coli [g/mol]

