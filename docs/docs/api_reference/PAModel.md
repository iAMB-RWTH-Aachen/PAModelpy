---
sidebar_label: PAModel
title: PAModel
---

## PAModel Objects

```python
class PAModel(Model)
```

Class representation for a cobra model extended with enzyme kinetics as published in Alter et al. (2021).

**Arguments**:

- `id_or_model` _str or Model_ - String to use as model id, or actual model to base new model on.
  If a string, it is used as input to load a model from. If a model, a new model object is instantiated with
  the same properties as the original model (default None).
- `name` _str, optional_ - Human-readable string to be model description (default None).
- `p_tot` _float, optional_ - Total protein concentration (condition-dependent) (unit g_prot/g_cdw) (default 0.285).
- `senstitivity` _bool_ - Boolean value whether or not a sensitivity analysis should be performed during each simulation.
  This sensitivity analysis will indicate to which extent individual constraints contribute to the objective value.
  Enzyme sectors (EnzymeSector objects, optional): Information about the different enzyme sectors, including:
  - Active_enzyme: Metabolic active proteins.
  - Transl_enzyme: Enzymes related to translation.
  - Unused_enzymes: Excess enzymes.
  - Custom_enzymes (list): Custom enzyme sectors.
- `configuration` _Config object, optional_ - Information about the general configuration of the model, including
  identifier conventions. Default as defined in the `PAModelpy.configuration` script for the E.coli iML1515 model.
  

**Attributes**:

- `p_tot` _float_ - The fraction of biomass allocated to proteins (units: g_prot/g_cdw).
- `reactions` _DictList_ - A DictList where the key is the reaction identifier and the value is a Reaction.
- `metabolites` _DictList_ - A DictList where the key is the metabolite identifier and the value is a Metabolite.
- `genes` _DictList_ - A DictList where the key is the gene identifier and the value is a Gene.
- `name`0 _DictList_ - A DictList where the key is the group identifier and the value is a Group.
- `name`1 _DictList_ - A DictList where the key is the enzyme identifier and the value is an Enzyme.
- `name`2 _DictList_ - A DictList where the key is the enzyme variable identifier and the value is an EnzymeVariable.
- `name`3 _DictList_ - A DictList where the key is the catalytic event identifier and the value is a CatalyticEvent.
- `name`4 _dict_ - A dictionary containing sector-specific constraints.
- `name`5 _DictList_ - A DictList where the key is the sector identifier and the value is an EnzymeSector.

#### P\_TOT\_DEFAULT

g_protein/g_cdw

#### \_\_init\_\_

```python
def __init__(id_or_model: Union[str, "Model", None] = None,
             name: Optional[str] = None,
             p_tot: Optional[float] = Config.P_TOT_DEFAULT,
             sensitivity: bool = True,
             active_sector: Optional[ActiveEnzymeSector] = None,
             translational_sector: Optional[TransEnzymeSector] = None,
             unused_sector: Optional[UnusedEnzymeSector] = None,
             custom_sectors: Union[List, CustomSector] = None,
             configuration=Config)
```

Constants

#### add\_enzymes

```python
def add_enzymes(enzyme_list: list) -> None
```

Add new enzymes to a model.
Adapted from Cobra.core.model.add_reactions and Cobra.core.model.add_metabolites.

This function will add a DictList of enzymes to the model object and add new variables accordingly.
For each enzyme-associated reaction, a constraint in each direction is added to the model.
The change is reverted upon exit when using the model as a context.

**Arguments**:

- `enzyme_list` _list or Enzyme_ - A list of `Enzyme` objects. If it isn&#x27;t an iterable container, the enzyme will
  be placed into a list.

#### add\_sectors

```python
def add_sectors(sectors: List = None)
```

Adds sector variables to the model and adds these to the total protein constraint.

**Arguments**:

- `sectors` _list_ - A list of PAModelpy.EnzymeSectors to add to the model.

#### add\_sector

```python
def add_sector(sector)
```

Adds the sector variable for a specific sector to the model and adds this to the total protein constraint.
Also stores the sector variables in the model attributes.

**Arguments**:

- `sector` _PAModelpy.EnzymeSector_ - The specific EnzymeSector to add to the model.

#### add\_catalytic\_events

```python
def add_catalytic_events(catalytic_events: Optional[Iterable])
```

Add a new CatalyticEvent to the model.
Will add a list of CatalyticEvent variables to the model object using the function defined in the CatalyticEvent object.

**Arguments**:

- `catalytic_events` _list or variables.CatalyticEvent_ - A list of `variables.CatalyticEvent` objects. If it isn&#x27;t
  an iterable container, the catalytic event will be placed into a list.

#### add\_enzyme\_constraints

```python
def add_enzyme_constraints(constraint_list: Optional[list])
```

Add new enzyme constraints to a model.
Will add a list of constraints to the model object and add new constraints accordingly.
The change is reverted upon exit when using the model as a context.

**Arguments**:

- `constraint_list` _list, str, or constraints.Constraint_ - A list of `constraints.Constraint` objects. If it isn&#x27;t
  an iterable container, the constraint will be placed into a list. Also, a string with the constraint id
  can be provided. A constraint will be created before adding it to the model.

#### add\_sector\_constraints

```python
def add_sector_constraints(constraint_list: Optional[list])
```

Add a new constraint related to a sector to a model.
Will add a list of constraints to the model object and add new constraints accordingly.
The change is reverted upon exit when using the model as a context.

**Arguments**:

- `constraint_list` _list or constraints.Constraint_ - A list of `constraints.Constraint` objects. If it isn&#x27;t an iterable
  container, the constraint will be placed into a list.

#### add\_total\_protein\_constraint

```python
def add_total_protein_constraint(p_tot: Optional[float] = P_TOT_DEFAULT)
```

Function which adds the total protein constraint to the model.
This limits the amount of available enzymes and thus the resulting fluxes.

**Notes**:

  The constraint expression looks like this:
- ``Etot` - sum(E) + E_translprot + E_unusedprot  == p_tot - E_trsn_0 - E_ue_0`
  

**Arguments**:

- `p_tot` _float, optional_ - Fraction of biomass which consists of protein (g_protein/g_cdw).
  Default is 0.258 (E.coli).

#### add\_reactions

```python
def add_reactions(reaction_list: Iterable[Reaction]) -> None
```

Add reactions to the model.
This method is superimposed upon the cobra.Model.add_reactions() function.
As a new feature, it will add constraints to determine the lower and upper bound if a sensitivity analysis should
be performed (which is determined by the model attribute: PAModel.sensitivity).
Reactions with identifiers identical to a reaction already in the model are ignored.
The change is reverted upon exit when using the model as a context.

**Arguments**:

- `reaction_list` _list_ - A list of `cobra.Reaction` objects.

#### add\_lb\_ub\_constraints

```python
def add_lb_ub_constraints()
```

Makes additional constraints for the reaction lower bounds and upperbounds.
By adding these constraints the shadow prices of the reaction bounds can be
calculated and used in sensitivity analysis

#### make\_lb\_ub\_constraint

```python
@staticmethod
def make_lb_ub_constraint(m: Optional[Model], rxn: Reaction,
                          lower_bound: float, upper_bound: float)
```

Adding variables and constraints for the lower and upper bounds of a reaction to a model.
When solving the model, shadow prices for the lower and upper bounds will be calculated.
This allows for the calculation of sensitivity coefficients. The constraints are formulated as follows:

**Notes**:

  Constraints are formulated as follows:
  - `R_ub: R_fwd - R_rev <= UB`
  - `R_lb: -(R_fwd - R_rev) <= -LB`
  

**Arguments**:

- `m` _cobra.Model or PAModelpy.PAModel_ - The model to which the upper and lower bound constraints and variables
  should be added.
- `rxn` _cobra.Reaction_ - The reaction for which upper and lower bound constraints should be generated.
- `lower_bound` _float_ - The value of the lower bound.
- `upper_bound` _float_ - The value of the upper bound.
  

**Returns**:

- `m` _cobra.Model or PAModelpy.PAModel_ - The model with additional constraints and variables for the reactions.

#### make\_enzyme\_min\_max\_constraint

```python
@staticmethod
def make_enzyme_min_max_constraint(m: Optional[Model], enz: Enzyme,
                                   lower_bound: float, upper_bound: float)
```

Adding variables and constraints for the lower and upper bounds of an enzyme&#x27;s concentration to a model.
When solving the model, shadow prices for the lower and upper bounds will be calculated.
This allows for the calculation of sensitivity coefficients.

**Notes**:

  The constraints are formulated as follows:
  - `E_ub: E <= Emax`
  - `E_lb: -E <= -Emin`
  

**Arguments**:

- `m` _cobra.Model or PAModelpy.PAModel_ - The model to which the upper and lower bound constraints and variables
  should be added.
- `rxn` _PAModelpy.Enzyme_ - The enzyme for which minimal and maximal concentration constraints should be generated.
- `lower_bound` _float_ - The value of the lower bound.
- `upper_bound` _float_ - The value of the upper bound.
  

**Returns**:

- `m` _cobra.Model or PAModelpy.PAModel_ - The model with additional constraints and variables for the enzyme&#x27;s
  concentration.

#### parse\_shadow\_prices

```python
@staticmethod
def parse_shadow_prices(shadow_prices)
```

Parse the shadow prices to a DataFrame where each constraint corresponds to a row, and shadow prices and directions are columns.

#### calculate\_csc

```python
def calculate_csc(obj_value, mu, mu_ub, mu_lb, mu_ec_f, mu_ec_b)
```

Calculate the capacity sensitivity coefficient for all inequality constraints in the model.
The sum of all capacity sensitivity coefficients should equal 1 for growth maximization.

Capacity Sensitivity Coefficient Calculation:
Capacity Sensitivity Coefficient = constraint_UB * shadowprice / obj_value

**Arguments**:

- `obj_value` _float_ - The objective value of the model.
- `mu` _DataFrame_ - Shadow prices for all constraints.
- `mu_ub` _DataFrame_ - Shadow prices for the reaction upper bound (UB) constraints.
- `mu_lb` _DataFrame_ - Shadow prices for the reaction lower bound (LB) constraints.
- `mu_ec_f` _DataFrame_ - Shadow prices for the constraints related to enzymatic catalysis of the forward reaction.
- `mu_ec_b` _DataFrame_ - Shadow prices for the constraints related to enzymatic catalysis of the backward reaction.
  

**Returns**:

  None
  
  Results are saved in the `self.capacity_sensitivity_coefficients` attribute as a DataFrame.

#### calculate\_esc

```python
def calculate_esc(obj_value, mu_ec_f, mu_ec_b)
```

Calculate enzyme sensitivity coefficients for the enzyme variables using their primal values,
the objective value, and shadow prices according to the following relations:

Enzyme Sensitivity Coefficient Calculation:
esc = enzyme_variable.primal * constraint.shadowprice / obj_value

**Arguments**:

- `obj_value` _float_ - The objective value from the most recent optimal solution.
- `mu_ec_f` _pd.DataFrame_ - Shadow prices for maximizing enzyme concentrations (forward variables).
- `mu_ec_b` _pd.DataFrame_ - Shadow prices for minimizing enzyme concentrations (reverse variables).
  

**Returns**:

  None
  
  Fills the `PAModel.enzyme_sensitivity_coefficients` dataframe with the calculated enzyme sensitivity coefficients.

#### calculate\_sum\_of\_enzymes

```python
def calculate_sum_of_enzymes()
```

Calculate the sum of all enzyme variables for a feasible solution.

**Returns**:

- `float` - The sum of all enzyme variables in milligrams per gram of cell dry weight per hour (mg/gCDW/h).

#### change\_total\_protein\_constraint

```python
def change_total_protein_constraint(p_tot)
```

Change the fraction of biomass that is allocated to active proteins.

**Arguments**:

- `p_tot` _float_ - The new proteome fraction in grams of protein per gram of cell dry weight (g_protein/g_cdw).

#### change\_reaction\_bounds

```python
def change_reaction_bounds(rxn_id: str,
                           lower_bound: float = None,
                           upper_bound: float = None)
```

Change the reaction bounds. If a sensitivity analysis is required, the bounds of the upper and lower bound
constraints are adjusted.

**Arguments**:

- `rxn_id` _str_ - The string representing the reaction identifier to change.
- `lower_bound` _float, optional_ - The new value for the lower bound of the reaction (default is None).
- `upper_bound` _float, optional_ - The new value for the upper bound of the reaction (default is None).

#### change\_enzyme\_bounds

```python
def change_enzyme_bounds(enzyme_id: str,
                         lower_bound: float = None,
                         upper_bound: float = None)
```

Change the enzyme bounds. If the model should be primed for performing a sensitivity analysis,
the upper bound of the minimum and maximum enzyme concentration constraints are adjusted.

**Arguments**:

- `enzyme_id` _str_ - The string representing the enzyme identifier to change.
- `lower_bound` _float, optional_ - The new value for the minimal enzyme concentration (default is None).
- `upper_bound` _float, optional_ - The new value for the maximal enzyme concentration (default is None).

#### get\_enzymes\_with\_reaction\_id

```python
def get_enzymes_with_reaction_id(rxn_id: str)
```

Return Enzyme objects associated with the reaction identifier through CatalyticEvent objects.

**Arguments**:

- `rxn_id` _str_ - The reaction identifier.
  

**Returns**:

- `DictList` - A DictList of Enzyme objects associated with the reaction.

#### get\_reactions\_with\_enzyme\_id

```python
def get_reactions_with_enzyme_id(enz_id: str)
```

Return a list of reaction identifiers associated with the enzyme identifier (EC number) through CatalyticEvent objects.

**Arguments**:

- `enz_id` _str_ - The enzyme identifier (EC number).
  

**Returns**:

- `List[str]` - A list of reaction identifiers associated with the enzyme.

#### change\_kcat\_value

```python
def change_kcat_value(enzyme_id: str, kcats: dict)
```

Change the turnover number (kcat) of the enzyme for a specific reaction.

**Arguments**:

- `enzyme_id` _str_ - The enzyme identifier.
- `kcats` _dict_ - A dictionary with reaction identifiers as keys and kcat values as values.
  Each kcat value should be a nested dictionary with `f` (forward) and `b` (backward) as keys,
  and the corresponding kcat values as values.
  

**Example**:

  Example dictionary for the `kcat` parameter
    ```
    {'R1': {'f': 10.0, 'b': 5.0}, 'R2': {'f': 7.0, 'b': 3.0}}
    ```

#### remove\_enzymes

```python
def remove_enzymes(
        enzymes: Union[str, Enzyme, List[Union[str, Enzyme]]]) -> None
```

Remove enzymes from the model.

**Arguments**:

- `enzymes` _list, reaction, or str_ - A list with enzymes (`Enzyme`), or their IDs, to remove.
  Enzymes will be placed in a list. Strings will be placed in a list
  and used to find the enzymes in the model.

**Notes**:

  The change is reverted upon exit when using the model as a context.

#### remove\_reactions

```python
def remove_reactions(reactions: Union[str, Reaction, List[Union[str,
                                                                Reaction]]],
                     remove_orphans: bool = False) -> None
```

Remove reactions from the model. Inherited from the cobrapy.core.remove_reactions() function.

**Arguments**:

- `reactions` _list, reaction, or str_ - A list with reactions (`cobra.Reaction`), or their IDs, to remove.
  Reactions will be placed in a list. Strings will be placed in a list
  and used to find the reactions in the model.
- `remove_orphans` _bool, optional_ - Remove orphaned genes and metabolites from the model as well (default False).

**Notes**:

  The change is reverted upon exit when using the model as a context. Also removes associated CatalyticEvents if they exist.

#### remove\_catalytic\_events

```python
def remove_catalytic_events(catalytic_events: Union[
    str, CatalyticEvent, List[Union[str, CatalyticEvent]]],
                            remove_orphans: bool = False) -> None
```

Remove catalytic events from the model.

**Arguments**:

- `reactions` _list, reaction, or str_ - A list with reactions (`cobra.Reaction`), or their IDs, to remove.
  Reactions will be placed in a list. Strings will be placed in a list
  and used to find the reactions in the model.
- `remove_orphans` _bool, optional_ - Remove orphaned genes and metabolites from the model as well (default False).
  

**Notes**:

  The change is reverted upon exit when using the model as a context.

#### remove\_sectors

```python
def remove_sectors(
    sectors: Union[
        str,
        Sector,
        ActiveEnzymeSector,
        List[Union[str, Sector, ActiveEnzymeSector]],
    ]
) -> None
```

Remove sections from the model.

Also removes associated CatalyticEvents if they exist.

**Arguments**:

- `sectors` _list, sector, or str_ - A list with sector (`PAModelpy.Sector` or `PAModelpy.ActiveEnzymeSector`),
  or their IDs, to remove. A single sector will be placed in a list.
  Strings will be placed in a list and used to find the sector in the model.

#### test

```python
def test(glc_flux: Union[int, float] = 10)
```

Test the proteome allocation model.

**Arguments**:

- `glc_flux` _float, optional_ - The glucose flux which limits the growth rate (units: mmol_glc/g_cdw/h, default=10).

#### pfba

```python
def pfba(fraction_of_optimum: float = 1.0,
         proteins: bool = False,
         reactions: bool = True,
         exclude: List["str"] = [],
         objective: Union[Dict, "Objective", None] = None)
```

Perform pFBA (parsimonious Enzyme Usage Flux Balance Analysis) with a custom objective including:
- All reactions
- All proteins
- All proteins and all reactions.

pFBA [1] adds the minimization of all fluxes to the objective of the model. This approach is motivated by the idea that high fluxes have a higher enzyme turnover, and since producing enzymes is costly, the cell will try to minimize overall flux while still maximizing the original objective function, e.g., the growth rate.

**Arguments**:

- `fraction_of_optimum` _float, optional_ - The fraction of optimum which must be maintained. The original objective reaction is constrained to be greater than the maximal value times the `fraction_of_optimum` (default 1.0).
- `objective` _dict or cobra.Model.objective, optional_ - A desired objective to use during optimization in addition to the pFBA objective. Dictionaries (reaction as the key, coefficient as the value) can be used for linear objectives (default None).
- `proteins` _bool, optional_ - Determines whether to include enzyme variables in the pFBA objective.
- `reactions` _bool, optional_ - Determines whether to include reaction variables in the pFBA objective.
- `exclude` _list of reaction ids, optional_ - Reactions to exclude from the minimization objective.
  

**References**:

  - [1] Lewis, N. E., Hixson, K. K., Conrad, T. M., Lerman, J. A., Charusanti, P., Polpitiya, A. D., Palsson, B. O. (2010). Omic data from evolved E. coli are consistent with computed optimal growth from genome-scale models. Molecular Systems Biology, 6, 390. doi:10.1038/msb.2010.47

#### reset\_objective

```python
def reset_objective()
```

Reseting the objective to the standard biomass maximization objective after pFBA

#### optimize

```python
def optimize(objective_sense: Optional[str] = None,
             raise_error: bool = False) -> "Solution"
```

Optimize the model using flux balance analysis. Inherits from the cobra.Model.optimize() function and performs a sensitivity analysis after optimization if this is desired (by setting the PAModel.sensitivity attribute to True).

**Arguments**:

- `objective_sense` _`{None, 'maximize', 'minimize'}`, optional_ - Whether fluxes should be maximized or minimized. In case of None, the previous direction is used (default None).
- `raise_error` _bool_ - If true, raise an OptimizationError if solver status is not optimal (default False).
  

**Returns**:

  Solution
  

**Notes**:

  Only the most commonly used parameters are presented here. Additional parameters for cobra.solver may be available and specified with the appropriate keyword argument.

#### copy

```python
def copy() -> "PAModel"
```

Provide a partial &#x27;deepcopy&#x27; of the Model.

Adjusted from cobra.Model.copy().

All the Metabolite, Gene, Reaction, Enzyme, EnzymeVariable, Sector, and CatalyticEvent objects are created anew but in a faster fashion than deepcopy.

**Returns**:

- `PAModelpy.PAModel` - A new model copy.

