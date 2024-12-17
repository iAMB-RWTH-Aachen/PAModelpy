# CatalyticEvent

CatalyticEvent object which relates Reaction variables to the EnzymeVariable and Enzyme objects.
It contains multiple functions which enable easy mapping and handling of one Event of catalysis
(e.g. one conversion of substrate to product, can be catalyzed by multiple enzymes)

## CatalyticEvent Objects

```python
class CatalyticEvent(Object)
```

CatalyticEvent is a class for holding information regarding the
catalysis of a Reaction in a cobra.Model object. It serves as an interface
between the metabolic reaction and the associated enzyme constraints and variables.

**Notes**:

  There are three different scenarios:
  - Enzyme complex: multiple enzymes together are associated with an EnzymeComplex object
  - Isozymes: multiple enzymes independently associated with a single catalytic event
  - Other: a single enzyme is associated with a single catalytic event
  

**Arguments**:

- `kcats2enzymes` _dict_ - A dictionary with enzyme, kcat key, value pairs to connect the enzyme with the associated reaction.
  The kcat is another dictionary with &#x27;f&#x27; and &#x27;b&#x27; for the forward and backward reactions, respectively.
- `id` _str, optional_ - The identifier to associate with this catalytic event (default None).
- `rxn_id` _str, optional_ - The reaction with which this catalytic event is associated.
- `name` _str, optional_ - A human-readable name for the reaction (default &quot;&quot;).

#### kcat\_values

```python
@property
def kcat_values()
```

returns a dictionary with kcat values and enzymes

#### flux

```python
@property
def flux() -> float
```

Get the flux value in the most recent solution.

Flux is the primal value of the corresponding variable in the model.

**Returns**:

- `flux` _float_ - Flux is the primal value of the corresponding variable in the model.
  

**Warnings**:

  * Accessing reaction fluxes through a `Solution` object is the safer,
  preferred, and only guaranteed to be correct way. You can see how to
  do so easily in the examples.
  * Reaction flux is retrieved from the currently defined
  `self._model.solver`. The solver status is checked but there are no
  guarantees that the current solver state is the one you are looking
  for.
  * If you modify the underlying model after an optimization, you will
  retrieve the old optimization values.
  

**Raises**:

- `RuntimeError` - If the underlying model was never optimized beforehand or the
  reaction is not part of a model.
- `OptimizationError` - If the solver status is anything other than `optimal`.
- `AssertionError` - If the flux value is not within the bounds.
  

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

- `float` - Enzyme concentration in [mmol/gDW].

#### add\_enzymes

```python
def add_enzymes(enzyme_kcat_dict: dict)
```

Add enzymes to the catalytic event and create bindings to the related model.
The enzymes in the enzyme_kcat_dict are individual isozymes. Enzyme complexes
should be added as an EnzymeComplex object with a single kcat value.

**Arguments**:

- `enzyme_kcat_dict` - Dict
  A nested dictionary with enzyme, kcat key, value pairs to connect the
  enzyme with the associated reaction. The kcat is another dictionary with `f` and `b`
  for the forward and backward reactions respectively.

#### remove\_enzymes

```python
def remove_enzymes(enzyme_list: list)
```

Remove enzymes from the catalytic event and remove the catalytic event from the
constraint expressions related to the enzyme.

**Arguments**:

- `enzyme_list` - List[Union[str, PAModelpy.Package.Enzyme]]
  A list with PAModelpy.Package.Enzyme objects to be removed. If a list of identifiers (str)
  is provided, the corresponding enzyme will be obtained from the CatalyticEvent.enzymes attribute.

#### change\_kcat\_values

```python
def change_kcat_values(enzyme_kcat_dict: dict)
```

changes kcat values for the enzyme variable

**Arguments**:

- `enzyme_kcat_dict` - nested Dict
  A Dict with enzyme_id, kcat key, value pairs to connect the
  enzyme with the associated reaction the kcat is another dict with &#x27;f&#x27; and &#x27;b&#x27;
  for the forward and backward reactions respectively.

#### \_\_copy\_\_

```python
def __copy__() -> "CatalyticEvent"
```

Copy the CatalyticEvent.

**Returns**:

  CatalyticEvent:
  A new CatalyticEvent that is a copy of the original CatalyticEvent.

#### \_\_deepcopy\_\_

```python
def __deepcopy__(memo: dict) -> "CatalyticEvent"
```

Copy the CatalyticEvent with memo.

**Arguments**:

- `memo` _dict_ - Automatically passed parameter.
  

**Returns**:

  CatalyticEvent:
  A new CatalyticEvent that is a copy of the original CatalyticEvent with memo.

