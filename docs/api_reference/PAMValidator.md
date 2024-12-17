# PAMValidator

## PAMValidator Objects

```python
class PAMValidator(object)
```

#### MW\_GLC

g/mol

#### GRADIENT\_MAX

mmol/gdw/h

#### GRADIENT\_STEP

mmol/gdw/h

#### GRADIENT\_MIN

mmol/gdw/h

#### run\_simulations\_glc\_o2\_gradient

```python
def run_simulations_glc_o2_gradient(
        oxygen_gradient: list,
        params_to_save: Union[str, list] = "R_TranslationalProteinSector")
```

Function to run simulations of different oxygen gradients for a range of growth rates.

This will simulate growth for the entire range of glucose concentrations for each oxygen uptake rate as given by the input.

**Arguments**:

- `oxygen_gradient` _list_ - List of upper bounds for the oxygen uptake reaction to loop over.
- `params_to_save` _optional_ - string or list, which parameter(s) to save for further analysis (default: translational protein sector constraint).
  

**Returns**:

- `results` _list of dataframes_ - Saves the growth rate, glucose uptake rate, and the user-defined parameters for each oxygen uptake rate in separate dataframes.

#### run\_simulations\_ups

```python
def run_simulations_ups(
        ups_gradient: list,
        params_to_save: Union[str, list] = "R_TranslationalProteinSector")
```

Function to run simulations with increasing unused enzyme sectors proportions for a range of growth rates.

This will simulate growth for the entire range of glucose concentrations for a range of fractions of ups_0 as given by the input.

**Arguments**:

- `ups_gradient` _list_ - List of upper bounds for the oxygen uptake reaction to loop over.
- `params_to_save` _optional_ - string or list, which parameter(s) to save for further analysis (default: translational protein sector constraint).
  

**Returns**:

- `results` _list of dataframes_ - Saves the growth rate, glucose uptake rate, and the user-defined parameters for each oxygen uptake rate in separate dataframes.

#### custom\_plot

```python
def custom_plot(rxn_ids: list,
                valid_dataframe: pd.DataFrame = None,
                xaxis: str = None,
                c_uptake_rxn: str = GLUCOSE_EXCHANGE_RXNID)
```

Function to plot the results of custom reactions.

**Arguments**:

- `rxn_ids` _list of str_ - Reaction identifiers of the reactions to be plotted.
- `valid_dataframe` _pandas.DataFrame, optional_ - A DataFrame with experimental data to validate the results with.
  The columns should be the same as the rxn_id of the reaction to be plotted and the reaction which should be plotted
  on the x-axis (by default the glucose exchange reaction `EX_glc__D_e_b`). If the DataFrame is not provided,
  only the simulation results will be plotted.
- `xaxis` _str, optional_ - The reaction identifier of the reaction which should be plotted on the x-axis (default: `EX_glc__D_e_b`).
  

**Returns**:

  Prints scatter plots of the model simulations vs. experimental data points (if provided).

