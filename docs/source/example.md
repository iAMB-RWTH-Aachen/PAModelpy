# Example usage of PAModelpy
******
## Example 1: setting up an *Escherichia coli* Protein Allocation model (PAM)
*Escherichia coli* (*E.coli*) is a commonly used model organism in Microbiology. When this microorganism is grown
on increasing glucose concentration, it shifts from a purely respiratory metabolism to a respiro-fermentative metabolic
phenotype. This phenomenon is called 'overflow metabolism'. Interestingly, overflow metabolism cannot be simulated using
normal genome-scale models without additional constraints. With properly parametrized protein-constrained models however,
we are able to simulate this metabolic phenotype. In this example, we'll set-up the *E.coli* PAM, and we'll study the
predicted metabolic phenotypes for a range of glucose uptake rates.

For this entire tutorial, you'll need to load the following packages:

``` python
    #importing the packages
    import os
    from cobra.io import read_sbml_model, load_matlab_model
    import sys
    import pandas as pd
    
    #load PAMpy modules
    from PAModelpy.EnzymeSectors import ActiveEnzymeSector, TransEnzymeSector, UnusedEnzymeSector
    from PAModelpy.PAModel import PAModel
    from PAModelpy.PAMValidator import PAMValidator
    from PAModelpy.configuration import Config
```


### Step 1: Initiate the protein sectors
Each protein-allocation model has three sectors: active enzyme sector (enzymes catalyzing the metabolic reactions), 
translational protein sectors (i.e. ribosomal proteins required for translation) and unused proteins (idle proteins
which help the cell adapt to new conditions). The total of these three sectors is limited by an upperbound. This 
upperbound is determined by the sum of all non-maintenance enzymes, which is assumed to be constant for prokaryotes. 
For examples on how to parametrize these sectors, refer to `Scripts/create_ecolicore_pam_inclUE.ipynb`.

#### 1.1: Active enzyme sector
The active enzyme sector will be build using information about which enzymes catalyzes a specific reaction, 
the turnover rate of the catalysis and the molar mass of the enzyme. In this example we'll use the parameters as 
published by [Alter et al. (2021)](https://journals.asm.org/doi/10.1128/mSystems.00625-20), which can be found in the `Data` folder of the PAModelpy repository

First, we'll define the paths we'll download the data

```python
    protein_sector_info_path = 'Data/proteinAllocationModel_iML1515_EnzymaticData_py.xls'
    active_enzyme_data = pd.read_excel(protein_sector_info_path, sheet_name='ActiveEnzymes'))
```

The data is now in a dataframe with the following columns:
`
rxn_id - rxnName - rxnEquat - EC_nmbr - molMass
`
First, let's add an identifier to the reactions for which the enzyme is unknown, in order to distinguish between the 
enzymes
```python
   #load active enzyme sector information
    active_enzyme_info = pd.read_excel(pam_info_file, sheet_name='ActiveEnzymes')
    
    # replace NaN values with unique identifiers
    #select the NaN values 
    nan_values = active_enzyme_info['EC_nmbr'].isnull()
    #make a list with unique ids
    nan_ids = [f'E{i}' for i in range(nan_values.sum())]
    #replace nan values by unique id
    active_enzyme_info.loc[nan_values, 'EC_nmbr'] = nan_ids
    
    #check if it worked:
    active_enzyme_info[nan_values]
```

We need to collect the data from this table and put it in the correct structure to be parsed into the ActiveEnzymeSector
object. The main input in this object is the rxn2protein dictionary, where all the information about protein-reaction 
associations required to build the protein-reaction relations in the model. It has the following format:
```
        {'R1':
             {'E1':
                  {'f': forward kcat, 'b': backward kcat, 'molmass': molar mass, 'genes': [G1, G2],
                        'complex_with': 'E2'},
                'E2':
                    {'f': forward kcat, 'b': backward kcat, 'molmass': molar mass, 'genes': [G3, G4],
                                 'complex_with': 'E1'}
        }
```

We need to take the following steps to get the right format:

```python
    # parse the enzyme information (kcat values, identifiers and molmasses)
    kcats_dict = active_enzyme_info.set_index(keys='rxnID').loc[:, 'kcat'].to_dict()
    ec_dict = active_enzyme_info.set_index(keys='rxnID').loc[:, 'EC_nmbr'].to_dict()
    molmass_dict = mol_mass=active_enzyme_info.set_index(keys='rxnID').loc[:,'molMass'].to_dict()
    
    
    kcats = {}
    # save fwd and bckw kcats separately in the form of: {rxn_id: {'f': kcat_f, 'b': kcat_b}}
    for rxn, kcat in kcats_dict.items():
        #reversible reaction
        if rxn[-2:] == '_f' or rxn[-2:] == '_b':
            direction = rxn[-1]
            #check if the reaction already exists in the kcat dictionary
            try: 
                kcats[rxn[:-2]][direction] = kcat
            except:
                kcats[rxn[:-2]] = {direction: kcat}
        #irreversible reaction
        else:
            kcats[rxn] = {'f': kcat}
    
    rxn2ec = {}
    #parse the enzyme identifiers for the reactions
    for rxn, ec in ec_dict.items():
        if rxn[-2:] == '_f' or rxn[-2:] == '_b':
            rxn = rxn[:-2]
        for enz in str(ec).split(','):
            rxn2ec[rxn] = enz.strip()
            
    molmass = {}
    #parse the enzyme molmasses for the reactions
    for rxn, mw in molmass_dict.items():
        if rxn[-2:] == '_f' or rxn[-2:] == '_b':
            rxn = rxn[:-2]
        molmass[rxn] = mw
        
    rxn2protein = {}
    for rxn, ec in rxn2ec.items():
        ec_dict = {**kcats[rxn], **{'molmass': molmass[rxn]}}
        #add enzyme to enzymes related to reaction if these are already stored
        if rxn in rxn2protein.keys():
            rxn2protein[rxn] = {**rxn2protein[rxn], **{ec:ec_dict}}
        #if not create new reaction entry
        else:
            rxn2protein[rxn] = {ec:ec_dict}
    
    active_enzyme_sector = ActiveEnzymeSector(rxn2protein=rxn2protein)
```

#### 1.2: Translational protein sector
The translational sector requires less parameters. We only need to define the reaction to which this section is proportional,
(for example the biomass pseudoreaction or substrate uptake rate), defined by `id_list`, and the slope (`tps_mu`) and 
intercept (`tps_0`) of this relation.

```python
translational_info = pd.read_excel(protein_sector_info_path, sheet_name='Translational')
id_list = [translational_info[translational_info.Parameter == 'id_list'].loc[0,'Value']]
translation_enzyme_sector = TransEnzymeSector(id_list=id_list,
                                                tps_0=[translational_info[translational_info.Parameter == 'tps_0'].loc[1,'Value']],
                                                tps_mu=[translational_info[translational_info.Parameter == 'tps_mu'].loc[2,'Value']],
                                                mol_mass=[translational_info[translational_info.Parameter == 'mol_mass'].loc[3,'Value']])
```

#### 1.3 Unused enzyme sector
The unused enzyme sector is defined in a very similar way as the translational protein sector. We assume that this 
section is absent when the microbe is growing at it's highest growth rate. We'll use this assumption to define the slope:
```python
unused_protein_info = pd.read_excel(pam_info_file, sheet_name='ExcessEnzymes')

ups_0 = unused_protein_info[unused_protein_info.Parameter == 'ups_0'].loc[2,'Value']
smax = unused_protein_info[unused_protein_info.Parameter == 's_max_uptake'].loc[1,'Value']
id_list =[unused_protein_info[translational_info.Parameter == 'id_list'].loc[0,'Value']]


unused_protein_sector = UnusedEnzymeSector(id_list=id_list,
                                            ups_mu=[ups_0/smax],
                                            ups_0=[ups_0],
                                            mol_mass=[unused_protein_info[unused_protein_info.Parameter == 'mol_mass'].loc[3,'Value']])
```

### Step 2: Building the model
Now we are ready to build the model! We'll need to determine a maximal protein concentration. Following [Alter et al. (2021)](https://journals.asm.org/doi/10.1128/mSystems.00625-20),
let's take 0.258 mmol/gcdw/h.
As a basis, we'll use the iML1515 model, created by [Monk et al. (2017)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6521705/).

NB: for the more advanced users who want to run the sensitivity analysis, please ensure that the `sensitivity` argument 
is set to `True`. If you are not interested in the sensitivity analysis, you can set `sensitivity` to `False`, which 
will speed up the computation time.

```python
    #load the genome-scale information
    model = read_sbml_model(os.path.join('Models', 'iML1515.xml'))
    
    #load the PAM with the genome-scale information and the information about the enzyme sectors
    pamodel = PAModel(id_or_model=model, 
                      p_tot=0.258,
                       active_sector=active_enzyme_sector, 
                      translational_sector=translation_enzyme_sector,
                      unused_sector=unused_protein_sector,
                            sensitivity =True
                     ) 
```

### Step 3: run and validate the model results
To see if the PAM is working we can run some dummy simulations. Also, the PAMValidator module has functions which 
allow for easy visualization of the model predictions vs measured experimental data.

```python
pamodel.test()
```
This is a simulation with a glucose uptake rate set to 10 mmol/gcdw/h.
We can easily change to a different substrate uptake rate, e.g. 5 mmol/gcdw/h by putting that in as an function argument

```python
pamodel.test(5)
```
In the PAMValidator object, you can find functions to run simulations over a range of substrate uptake rates.
To initiate the PAMValidator, you need to provide the model and a path to an excel file with experimental data. 
In this example,we'll use the experimental data which can be found in 
`Data/Ecoli_phenotypes/Ecoli_phenotypes_py_rev.xls`.

```python
from PAModelpy.PAMValidator import PAMValidator

validator = PAMValidator(pamodel,'Data/Ecoli_phenotypes/Ecoli_phenotypes_py_rev.xls')
#model flux rates of biomass formation, acetate, CO2 and O2 vs glucose uptake rate for a range of growth rates
validator.validate_range()
```

Alternatively, you can run simulations in a good old-fashioned for-loop. For examples on how to do that look at the 
jupyter notebooks in the `Figures` directory.

### Step 4: interpreting the results
What does the result tell you? What is the predicted metabolic phenotype and how does this relate to the experimental 
results. Did the model capture overflow metabolism?

### Outlook
After this tutorial, you know how to apply PAModelpy to this very well-studied *E.coli* example. But how to address you're
specific issue with your specific microbe? In the next example we'll show you how to use the Config() object to give the 
PAModel the right naming conventions for your specific microbe. Are you more interested in performing modifications to 
the model, such as deleting or adding enzymes, changing kcats, changing enzymes upper- and lowerbounds? Then have a look
at the following jupyter notebook: `Examples/PAModel_example_script.ipynb`. Have fun!

### Example 2: Determining the most sensitive enzymes in a toy model

