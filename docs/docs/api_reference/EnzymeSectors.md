---
sidebar_label: EnzymeSectors
title: EnzymeSectors
---

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
def __init__(rxn2protein: dict, configuration: Config = Config)
```

_summary_

**Arguments**:

- `rxn2protein` _dict_ - A dictionary with reaction ID, enzymes_dict key, value pairs for each reaction in the active_enzyme sector.
  The enzymes_dict contains the enzyme identifiers of the enzymes related to the specific reaction as keys, and a dict
  with information about the enzyme as values. The information included in this dict includes the turnover number for
  the forward and backward reaction (1/s), molar mass of the enzyme (mol/g), gene identifiers related to the enzyme,
  and with which other enzymes it forms a complex.
- `configuration` _Config object, optional_ - Information about the general configuration of the model, including identifier conventions.
  Default is as defined in the `PAModelpy.configuration` script for the E.coli iML1515 model.
  

**Example**:

    ```
    For the Parameter rxn2protein a dictionary may look like this:
    {
        'R1':
            {
                'E1': {'f': forward kcat, 'b': backward kcat, 'molmass': molar mass, 'genes': ['G1', 'G2'],
                       'complex_with': 'E2'},
                'E2': {'f': forward kcat, 'b': backward kcat, 'molmass': molar mass, 'genes': ['G3', 'G4'],
                       'complex_with': 'E1'}
            }
    }
    ```

#### check\_kcat\_values

```python
def check_kcat_values(model, reaction, kcat)
```

Function to check if the kcat values provided for an enzyme are consistent with the direction of the reaction.

**Arguments**:

- `model` _cobra.Model or PAModel_ - Model to which the kcat values should be added.
- `reaction` _cobra.Reaction_ - Reaction that is catalyzed with the enzyme related to the kcats.
- `kcat` _dict_ - A dictionary with the turnover values for the forward and/or backward reaction for different enzymes [/s].
  

**Example**:

  Example dictionary for the `kcat` parameter
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

