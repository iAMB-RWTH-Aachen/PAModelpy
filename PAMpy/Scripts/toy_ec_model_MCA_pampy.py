import cobra
import pandas as pd
from cobra import Model, Reaction, Metabolite, core
from cobra.util import Zero
import string
import plotly.express
import sys
from optlang import interface

sys.path.append('../src/')
#importing the tools from the PAModelpy package
from PAMpy import ActiveEnzymeSector, UnusedEnzymeSector, PAModel
#need to have gurobipy installed

#global variables:
global metabolites, n, m, Etot
#global variables:
global metabolites, n, m, Etot
metabolites = ['Glucose', 'ATP', 'CO2', 'Precursor', 'Biomass']
n = 7
m = 5
Etot = 0.0045

#functions:
def build_toy_model():
    '''
    Rebuild the toymodel as in the MATLAB script.
    S = [1,0,0,0,-1,-1,-1,0,0,0;...
    0,1,0,0,1,0,0,-1,-1,0;...
    0,0,0,0,0,1,0,1,0,-1;...
    0,0,0,0,0,0,1,0,0,-1;...
    0,0,0,-1,0,0,0,0,0,1;...
    0,0,-1,0,0,0,0,0,1,1];

    :return: Cobrapy model instance as model
    '''
    #set up model basics
    model = Model('toy_model')
    cobra_config = cobra.Configuration()
    cobra_config.solver = 'gurobi'
    for i in range(1, n + 1):
        rxn = Reaction('R' + str(i))
        lower_bound = 0
        upper_bound = 1e6
        #         if i ==1 :
        #             lower_bound = 10
        #         else:
        #             lower_bound=-10
        if i != 1 and i != 2:
            upper_bound = 100
        else:
            upper_bound = 10

        rxn.lower_bound = lower_bound
        rxn.upper_bound = upper_bound
        model.add_reactions([rxn])


    # add metabolites to the reactions:
    # R1:
    r1 = model.reactions.get_by_id('R1')
    r1.add_metabolites({Metabolite('Glucose'): 1})
    # R2:
    r2 = model.reactions.get_by_id('R2')
    r2.add_metabolites({Metabolite('ATP'): 1, Metabolite('CO2'): 1, Metabolite('Glucose'): -1})
    # R3:
    r3 = model.reactions.get_by_id('R3')
    r3.add_metabolites({Metabolite('Precursor'): 1, Metabolite('Glucose'): -1})
    # R4:
    r4 = model.reactions.get_by_id('R4')
    r4.add_metabolites({Metabolite('ATP'): -1, Metabolite('Biomass'): 1})
    # R5:
    r5 = model.reactions.get_by_id('R5')
    r5.add_metabolites({Metabolite('Precursor'): -1, Metabolite('Biomass'): 1})
    # R6:
    r6 = model.reactions.get_by_id('R6')
    r6.add_metabolites({Metabolite('Biomass'): -1})
    # R7:
    r7 = model.reactions.get_by_id('R7')
    r7.add_metabolites({Metabolite('CO2'): -1})

    return model

def build_active_enzyme_sector():
    kcat_fwd = [1, 0.5, 0.25, 1.5, 0.5]  # High turnover for exchange reactions
    kcat_rev = [kcat / 2 for kcat in kcat_fwd]
    rxn2kcat = {}
    for i in range(n-2):
        rxn_id = f'R{i+1}'
        # 1e-6 to correct for the unit transformation in the model (meant to make the calculations preciser for different unit dimensions)
        #dummy molmass like in MATLAB script
        rxn2kcat = {**rxn2kcat, **{rxn_id: {f'E{i+1}':{'f': kcat_fwd[i]/(3600*1e-6), 'b': kcat_rev[i]/(3600*1e-6), 'molmass': 1e6}}}}

    return ActiveEnzymeSector(rxn2protein = rxn2kcat)

def build_unused_protein_sector():
    return UnusedEnzymeSector(id_list = ['R1'], ups_mu=[0.0001], ups_0=[0.001], mol_mass= [1])

# def build_translational_protein_sector():
#     return TranslationalEnzymeSector(id_list = ['R1'], tps_mu=[0.01*1e-3], tps_0=[0.1*1e-3], mol_mass= [1])

def print_heatmap(xaxis, matrix):
    yaxis = list()
    for i in range(1, n + 1):
        yaxis += [f'R{i}']
    fig = plotly.express.imshow(matrix, aspect="auto",
                                x = xaxis, y = yaxis,
                                labels = dict(x = 'control coefficients', y='maximized reaction'))
    fig.show()

if __name__ == "__main__":
    Cv=list()
    Ca = list()
    x_axis_cc = list()
    x_axis_ac = list()
    model = build_toy_model()
    active_enzyme = build_active_enzyme_sector()
    unused_enzyme = build_unused_protein_sector()
    # translation_enzyme = build_translational_protein_sector()
    pamodel = PAModel(model, name='toy model MCA with enzyme constraints', active_sector=active_enzyme,
                      # translational_sector = translation_enzyme, p_tot = Etot)
                      unused_sector = unused_enzyme, p_tot=Etot)


    for rxn in pamodel.reactions:
        #only get the network reactions, not the upperbound and lowerbound pseudo-reactions
        if 'b' not in rxn.id:
            pamodel.objective = {rxn: 1}
            pamodel.optimize()

            Cv_new = list()
            for cc in ['UB', 'LB', 'EC_f', 'EC_b', 'sector']:
                Cv_new += pamodel.control_coefficients[pamodel.control_coefficients['constraint'] == cc].coefficient.to_list()
            Cv += [Cv_new]

            Ca_new = list()
            for ac in ['rxn', 'enzyme', 'sector']:
                Ca_new += pamodel.allocation_coefficients[pamodel.allocation_coefficients['constraint'] == ac].coefficient.to_list()
            Ca += [Ca_new]

        print('Sum of control coefficients: \t \t \t \t \t \t \t \t', sum(Cv_new))
        print('Sum of allocation coefficients: \t \t \t \t \t \t \t', sum(Ca_new))

    for cc in ['UB', 'LB', 'EC_f', 'EC_b', 'sector']:
        if cc == 'UB' or cc == 'LB':
            x_axis_cc += [rid+'_'+cc for rid in pamodel.control_coefficients[pamodel.control_coefficients['constraint'] == cc].rxn_id.to_list()]
        else:
            x_axis_cc += [rid+'_'+cc[-1] for rid in pamodel.control_coefficients[pamodel.control_coefficients['constraint'] == cc].enzyme_id.to_list()]

    for ac in ['rxn', 'enzyme', 'sector']:
        if ac == 'rxn':
            x_axis_ac += pamodel.allocation_coefficients[pamodel.allocation_coefficients['constraint'] == ac].rxn_id.to_list()
        else:
            x_axis_ac += pamodel.allocation_coefficients[
                pamodel.allocation_coefficients['constraint'] == ac].enzyme_id.to_list()

    print_heatmap(x_axis_cc, Cv)
    print_heatmap(x_axis_ac, Ca)
