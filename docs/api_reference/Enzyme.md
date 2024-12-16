# Enzyme classes

Classes related to enzymes:
- Enzyme: Constraints relating enzymes to reactions. Including upper and lower bound enzyme constraints
- EnzymeComplex: Constraints relating enzyme complexes to reactions. Including upper and lower bound enzyme constraints
- EnzymeVariable: Variable related to an enzyme. The value of this variable represent the concentration.

## Enzyme Objects

```python
class Enzyme(Object)
```

Upper level Enzyme object containing information about the enzyme and links to the EnzymeVariables for each reaction the enzyme catalyzes.

**Arguments**:

  - id (str): Identifier for the enzyme (e.g., Uniprot ID).
  - rxn2kcat (Dict): Dictionary with reaction ID, kcat value pairs for the forward (f) and backward (b) reaction,
  e.g. `{'PGI': {'f': 30, 'b': 0.1}}`
  - upper_bound (float): Upper bound for the enzyme variable (default 1000.0).
  - lower_bound (float): Lower bound for the enzyme variable (default 0).
  - name (str): Name of the enzyme (default None).
  - molmass (float): Molar mass of the enzyme (default 3.947778784340140e04).
  

**Notes**:

  - This class is used to generate enzyme instances from kcat values and contains information about the forward as well as the backward catalysis.
  - The enzyme is linked to individual cobra.Reaction variables with CatalyticEvent objects.
  - There are two scenarios:
  - Promiscuous enzymes: a single enzyme can catalyze multiple reactions.
  - Other: a single enzyme catalyzes a single reaction.

#### DEFAULT\_ENZYME\_MOL\_MASS

mean enzymes mass E.coli [g/mol]

#### kcat\_values

```python
@property
def kcat_values()
```

Returns a dictionary with kcat values for each associated reaction.

**Returns**:

- `dict` - A dictionary containing kcat values for associated reactions.

#### concentration

```python
@property
def concentration(units: str = "mmol/gDW",
                  return_units: bool = False) -> float
```

Returns the enzyme&#x27;s total concentration considering any associated reactions.

**Arguments**:

- `units` _str, optional_ - Units in which the concentration is calculated (default is &#x27;mmol/gDW&#x27;), other option is &#x27;g/gDW&#x27;.
- `return_units` _bool, optional_ - Determines whether the units should be returned as well.
  

**Returns**:

- `float` - Enzyme concentration as a float.

#### add\_catalytic\_event

```python
def add_catalytic_event(ce: CatalyticEvent, kcats: Dict)
```

Adds a catalytic event associated with a reaction to an enzyme.

**Arguments**:

- `ce` _PAModelpy.Variables.CatalyticEvent_ - A CatalyticEvent object to which the enzyme should be added.
- `kcats` _dict_ - A dictionary containing direction and kcat key-value pairs.
  

**Returns**:

- `NoneType` - None

#### add\_genes

```python
def add_genes(gene_list: list,
              gene_length: list,
              relation: str = 'OR') -> None
```

Add genes to the enzyme and the model related to the enzyme if applicable

**Arguments**:

- `gene_list` _list_ - A list of gene identifiers to be added.
- `gene_length` _list_ - A list of lengths corresponding to each gene.
- `relation` _str, optional_ - The relationship between genes in gene_list.
  Defaults to &#x27;OR&#x27;. Possible values: &#x27;OR&#x27; or &#x27;AND&#x27;.
  

**Raises**:

- `ValueError` - If an invalid relation is provided.
  

**Notes**:

  If relation is &#x27;OR&#x27;, each gene in gene_list will be treated as coding for an individual isozyme
  If relation is &#x27;AND&#x27;, all genes in gene_list will be treated as coding for peptides in an enzyme complex

#### create\_catalytic\_event

```python
def create_catalytic_event(rxn_id: str, kcats: Dict)
```

Creates enzyme variables that link to reactions.

**Arguments**:

- `rxn_id` _str_ - ID of the associated reaction in the model.
- `kcats` _Dict_ - A dictionary containing kcat values for the forward and backward reaction.
  

**Returns**:

- `Variables.CatalyticEvent` - Enzyme variable object of type Variables.CatalyticEvent.

#### create\_enzyme\_variable

```python
def create_enzyme_variable()
```

Creates enzyme variables that link  enzyme to reactions.

#### change\_kcat\_values

```python
def change_kcat_values(rxn2kcat: Dict)
```

Changes the kcat values for the enzyme and updates the enzyme variable (enzymatic reaction) accordingly.

**Arguments**:

- `rxn2kcat` _Dict_ - A dictionary with reaction ID, kcat value pairs for the forward (f) and backward (b) reaction,
  e.g. `{'PGI': {'f': 30, 'b': 0.1}}`

#### get\_kcat\_values

```python
def get_kcat_values(rxn_ids: Union[str, list] = None) -> Dict
```

Returns the kcat values for a specific enzyme and all enzyme-associated reactions.

**Arguments**:

- `rxn_ids` _Union[str, list], optional_ - ID of the reactions for which the kcat values should be returned. It can be a single reaction ID (str) or a list of reaction IDs.
  

**Returns**:

- `Dict` - A dictionary containing kcat values for the forward (f) and backward (b) reactions.

#### remove\_catalytic\_event

```python
def remove_catalytic_event(catalytic_event: Union[CatalyticEvent, str])
```

Function to remove a catalytic event from an enzyme.

**Arguments**:

- `catalytic_event` _Union[CatalyticEvent, str]_ - CatalyticEvent or str, catalytic event or identifier to remove.

#### \_\_copy\_\_

```python
def __copy__() -> "Enzyme"
```

Copy the enzyme variable.

**Returns**:

- `PAModelpy.Enzyme.Enzyme` - A new enzyme that is a copy of the original enzyme.

#### \_\_deepcopy\_\_

```python
def __deepcopy__(memo: dict) -> "Enzyme"
```

Copy the enzyme variable with memo.

**Arguments**:

- `memo` _dict_ - Automatically passed parameter.
  

**Returns**:

- `PAModelpy.Enzyme.Enzyme` - A new enzyme that is a copy of the original enzyme with memo.

## EnzymeComplex Objects

```python
class EnzymeComplex(Enzyme)
```

Upper-level EnzymeComplex object containing information about the enzymes in a complex
and a link to the enzyme variables (CatalyticEvents) for each reaction the enzyme complex catalyzes.

This class is used to generate enzyme instances from kcat values and contains
information about the forward as well as the backward catalysis.

**Arguments**:

- `id` _str_ - Identifier for the enzyme complex (e.g., Uniprot ID).
- `enzymes` _DictList of cobra.core.Enzyme_ - Enzyme objects associated with the enzyme complex.
- `rxn2kcat` _Dict_ - Dictionary with reaction ID, kcat value pairs for the forward (f) and backward (b) reaction,
  e.g. `{'PGI': {'f': 30, 'b': 0.1}}`
- `upper_bound` _float, optional_ - Upper bound for the enzyme variable (default 1000.0).
- `name` _str, optional_ - Name of the enzyme (default None).
- `molmass` _float, optional_ - Molar mass of the enzyme (default 3.947778784340140e04).

#### DEFAULT\_ENZYME\_MOL\_MASS

mean enzymes mass E.coli [g/mol]

## EnzymeVariable Objects

```python
class EnzymeVariable(Reaction)
```

EnzymeVariable is a class for holding information regarding the variable representing an enzyme in the model.
For each reaction, the enzyme variables are summarized in a CatalyticEvent.

There are three different scenarios:
- Enzyme complex: multiple enzymes together are associated with an EnzymeComplex object.
- Isozymes: multiple enzymes independently associated with a single catalytic event.
- Other: a single enzyme is associated with a single catalytic event.

**Arguments**:

- `kcats2rxns` _Dict_ - A dictionary with reaction_id, kcat key, value pairs to connect the enzyme with the associated reaction.
  The kcat is another dictionary with `f` and `b` for the forward and backward reactions, respectively.
- `id` _str, optional_ - The identifier to associate with this enzyme (default None).
- `name` _str, optional_ - A human-readable name for the enzyme (default &quot;&quot;).
- `subsystem` _str, optional_ - Subsystem where the enzyme is meant to function (default &quot;&quot;).
- `lower_bound` _float_ - The lower flux bound (default 0.0).
- `upper_bound` _float, optional_ - The upper flux bound (default None).
- `**kwargs` - Additional keyword arguments are passed on to the parent class.

#### DEFAULT\_ENZYME\_MOL\_MASS

mean enzymes mass E.coli [g/mol]

#### kcat\_values

```python
@property
def kcat_values()
```

Returns a dictionary with kcat values and reactions.

**Returns**:

- `dict` - A dictionary containing kcat values and their associated reactions.

#### flux

```python
@property
def flux() -> float
```

Get the flux value in the most recent solution.

Flux is the primal value of the corresponding variable in the model.

**Returns**:

- `float` - Flux, which is the primal value of the corresponding variable in the model.
  

**Raises**:

- `RuntimeError` - If the underlying model was never optimized beforehand or the reaction is not part of a model.
- `OptimizationError` - If the solver status is anything other than &#x27;optimal&#x27;.
- `AssertionError` - If the flux value is not within the bounds.
  

**Warnings**:

  - Accessing reaction fluxes through a `Solution` object is the safer, preferred, and only guaranteed to be correct way.
  - Reaction flux is retrieved from the currently defined `self._model.solver`. The solver status is checked but there are no guarantees that the current solver state is the one you are looking for.
  - If you modify the underlying model after an optimization, you will retrieve the old optimization values.
  

**Examples**:

    ```
    >>> from cobra.io import load_model
    >>> model = load_model("textbook")
    >>> solution = model.optimize()
    >>> model.variables.PFK.flux
    7.477381962160283
    >>> solution.fluxes.PFK
    7.4773819621602833
    ```

#### concentration

```python
@property
def concentration() -> float
```

Get the enzyme concentration value of the most recent solution.

The enzyme concentration equals the flux value.

**Returns**:

- `float` - Enzyme concentration in mmol/gDW.

#### concentration

```python
@concentration.setter
def concentration(conc: Union[float, int]) -> None
```

Sets the concentration of the enzyme by creating or updating a constraint
that enforces the concentration to be equal to the sum of the forward and reverse
variable primals.

**Arguments**:

  conc : float, int
  The concentration value to be set for the enzyme. This value will be used
  as both the lower and upper bound for the constraint, effectively fixing the
  concentration to this value.
  
  Notes
  -----
  - If a concentration constraint for the enzyme does not already exist in the model,
  this function creates a new constraint named &#x27;&lt;enzyme_id&gt;_conc&#x27;.
  - The concentration constraint is defined as:
  concentration = forward_variable.primal + reverse_variable.primal
  - If the constraint already exists, the linear coefficients for the forward and reverse
  variables are updated to ensure the constraint remains valid.
  
  Raises
  ------
  ValueError
  If `conc` is not a valid numerical value.

#### set\_forward\_concentration

```python
def set_forward_concentration(conc: Union[float, int]) -> None
```

Sets the concentration of the enzyme by creating or updating a constraint
that enforces the concentration to be equal to the sum of only the forward
variable primals. This forces a reaction to be active in the forward direction.
It used the concentration setter functionality and subsequently sets the
coefficient for the reverse variable in the constraint to 0.

**Arguments**:

  conc : float, int
  The concentration value to be set for the enzyme. This value will be used
  as both the lower and upper bound for the constraint, effectively fixing the
  concentration to this value.
  
  Notes
  -----
  - If a concentration constraint for the enzyme does not already exist in the model,
  this function creates a new constraint named &#x27;&lt;enzyme_id&gt;_conc&#x27;.
  - The concentration constraint is defined as:
  concentration = forward_variable.primal
  
  Raises
  ------
  ValueError
  If `conc` is not a valid numerical value.

#### set\_reverse\_concentration

```python
def set_reverse_concentration(conc: Union[float, int]) -> None
```

Sets the concentration of the enzyme by creating or updating a constraint
that enforces the concentration to be equal to the sum of only the reverse
variable primals. This forces a reaction to be active in the reverse direction.
It used the concentration setter functionality and subsequently sets the
coefficient for the forward variable in the constraint to 0.

**Arguments**:

  conc : float, int
  The concentration value to be set for the enzyme. This value will be used
  as both the lower and upper bound for the constraint, effectively fixing the
  concentration to this value.
  
  Notes
  -----
  - If a concentration constraint for the enzyme does not already exist in the model,
  this function creates a new constraint named &#x27;&lt;enzyme_id&gt;_conc&#x27;.
  - The concentration constraint is defined as:
  concentration = reverse_variable.primal
  
  Raises
  ------
  ValueError
  If `conc` is not a valid numerical value.

#### add\_catalytic\_events

```python
def add_catalytic_events(catalytic_events: list, kcats: list)
```

Adding a catalytic event to an enzyme variable

**Arguments**:

- `catalytic_events` _list_ - A list of catalytic events to add.
- `kcats` _list_ - A list with dictionaries containing direction and kcat key-value pairs.

#### add\_reactions

```python
def add_reactions(reaction_kcat_dict: dict)
```

Add reactions to the enzyme variable and create bindings to the related model.
If there are multiple reactions related to a single enzyme, this is an isozyme.

**Arguments**:

- `reaction_kcat_dict` _dict_ - A nested dictionary with the reaction_id, kcat key, value pairs to connect the
  enzyme with the associated reaction. The kcat is another dictionary with `f` and `b` for the forward and
  backward reactions, respectively.

#### remove\_catalytic\_event

```python
def remove_catalytic_event(catalytic_event: Union[CatalyticEvent, str])
```

Remove a catalytic event from an enzyme.

**Arguments**:

- `catalytic_event` _Union[CatalyticEvent, str]_ - CatalyticEvent or str, catalytic event or identifier to remove.

#### remove\_reactions

```python
def remove_reactions(reaction_list: list)
```

Remove reactions from the enzyme variable and remove the reactions from the
constraint expressions related to the enzyme.

**Arguments**:

- `reaction_list` _list_ - A list with Cobra.Reaction objects which should be removed. If a list of identifiers (str)
  is provided, the corresponding enzyme will be obtained from the EnzymeVariables.reaction attribute.

#### change\_kcat\_values

```python
def change_kcat_values(reaction_kcat_dict: dict)
```

Changes kcat values for the enzyme variable.

**Arguments**:

- `reaction_kcat_dict` _dict_ - A nested dictionary with Cobra.Reaction, kcat key, value pairs to connect the
  enzyme with the associated reaction. The kcat is another dictionary with `f` and `b` for the forward and
  backward reactions, respectively.

#### \_\_copy\_\_

```python
def __copy__() -> "PAModelpy.Enzyme.EnzymeVariable"
```

Copy the enzyme variable.

**Returns**:

- `PAModelpy.Enzyme.EnzymeVariable` - A new enzyme variable that is a copy of the original enzyme variable.

#### \_\_deepcopy\_\_

```python
def __deepcopy__(memo: dict) -> "PAModelpy.Enzyme.EnzymeVariable"
```

Copy the enzyme variable with memo.

**Arguments**:

- `memo` _dict_ - Automatically passed parameter.
  

**Returns**:

- `PAModelpy.Enzyme.EnzymeVariable` - A new enzyme variable that is a copy of the original enzyme variable with memo.

