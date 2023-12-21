from cobra import Configuration
from cobra import Model, Reaction, Metabolite

import numpy as np

#importing the tools from the PAModelpy package
from PAModelpy.EnzymeSectors import ActiveEnzymeSector, UnusedEnzymeSector, TransEnzymeSector
from PAModelpy.PAModel import PAModel
from PAModelpy.configuration import Config

Config.BIOMASS_REACTION = 'R7'
Config.GLUCOSE_EXCHANGE_RXNID = 'R1'
Config.CO2_EXHANGE_RXNID = 'R8'
Config.ACETATE_EXCRETION_RXNID = 'R9'

#need to have gurobipy installed


#global variables:
global metabolites, n, m, Etot
#global variables:
global metabolites, n, m, Etot
metabolites = ['Substrate', 'ATP', 'CO2', 'Precursor', 'Biomass', 'Byproduct', 'Intermediate']
n = 9
m = 7
Etot = 0.6*1e-3 #will be adjusted in the model with 1e3

#functions:
def build_toy_gem():
    '''
    Rebuild the toymodel as in the MATLAB script.
    sub int byp atp co2 pre bio
R1 = [ 1,  0,  0,  0,  0,  0,  0];
R2 = [-1,  1,  0,  0,  1,  0,  0];
R3 = [ 0, -1,  1,  1,  0,  0,  0];
R3r= -R3;
R4 = [ 0, -1,  0,  2,  1,  0,  0];
R5 = [ 0, -1,  0,  0,  0,  1,  0];
R6 = [ 0,  0,  0, -1,  0, -1,  1];
R7 = [ 0,  0,  0,  0,  0,  0, -1];
R8 = [ 0,  0,  0,  0, -1,  0,  0];
R9 = [ 0,  0, -1,  0,  0,  0,  0];
S  = [R1;R2;R3;R3r;R4;R5;R6;R7;R8;R9]';

    :return: Cobrapy model instance as model
    '''
    #set up model basics
    model = Model('toy_model')
    cobra_config = Configuration()
    cobra_config.solver = 'gurobi'
    for i in range(1, n + 1):
        rxn = Reaction('R' + str(i))
        lower_bound = 0
        upper_bound = 1e6
        #force flux through the system
        if i == 1:
            lower_bound = 1
        #reversible reactions 3, 5 and 9
        if i ==3 or i==5 or 1==9:
            lower_bound =  -1e6
        #constrain nutrient (substrate or byproduct) uptake rate
        if i != 1 or i != 9:
            upper_bound = 100
        else:
            upper_bound = 10

        rxn.lower_bound = lower_bound
        rxn.upper_bound = upper_bound
        model.add_reactions([rxn])


    # add metabolites to the reactions:
    # R1:
    r1 = model.reactions.get_by_id('R1')
    r1.add_metabolites({Metabolite('Substrate'): 1})
    # R2:
    r2 = model.reactions.get_by_id('R2')
    r2.add_metabolites({Metabolite('Substrate'): -1, Metabolite('Intermediate'): 1, Metabolite('CO2'): 1})
    # R3:
    r3 = model.reactions.get_by_id('R3')
    r3.add_metabolites({Metabolite('Intermediate'): -1, Metabolite('Byproduct'):1, Metabolite('ATP'):1})
    # R4:
    r4 = model.reactions.get_by_id('R4')
    r4.add_metabolites({Metabolite('Intermediate'): -1, Metabolite('ATP'): 2, Metabolite('CO2'):1})
    # R5:
    r5 = model.reactions.get_by_id('R5')
    r5.add_metabolites({Metabolite('Intermediate'): -1, Metabolite('Precursor'): 1})
    # R6:
    r6 = model.reactions.get_by_id('R6')
    r6.add_metabolites({Metabolite('ATP'): -1, Metabolite('Precursor'): -1, Metabolite('Biomass'): 1})
    # Exchange reactions
    # R7:
    r7 = model.reactions.get_by_id('R7')
    r7.add_metabolites({Metabolite('Biomass'): -1})
    # R8:
    r8 = model.reactions.get_by_id('R8')
    r8.add_metabolites({Metabolite('CO2'): -1})
    # R9:
    r9 = model.reactions.get_by_id('R9')
    r9.add_metabolites({Metabolite('Byproduct'): -1})

    return model

def build_active_enzyme_sector(Config):
    kcat_fwd = [1, 0.5, 1, 1, 0.5 ,0.45, 1.5]  # High turnover for exchange reactions
    kcat_rev = [kcat for kcat in kcat_fwd]
    rxn2kcat = {}
    for i in range(n-3): # all reactions have an enzyme, except excretion reactions
        rxn_id = f'R{i+1}'
        # 1e-6 to correct for the unit transformation in the model (meant to make the calculations preciser for different unit dimensions)
        #dummy molmass like in MATLAB script
        rxn2kcat = {**rxn2kcat, **{rxn_id: {f'E{i+1}':{'f': kcat_fwd[i]/(3600*1e-6), 'b': kcat_rev[i]/(3600*1e-6), 'molmass': 1e6}}}}

    return ActiveEnzymeSector(rxn2protein = rxn2kcat, configuration=Config)

def build_unused_protein_sector(Config):
    return UnusedEnzymeSector(id_list = ['R1'], ups_mu=[-0.01*1e-3], ups_0=[0.1*1e-3], mol_mass= [1], configuration=Config)

def build_translational_protein_sector(Config):
    return TransEnzymeSector(id_list = ['R7'], tps_mu=[0.01*1e-3], tps_0=[0.01*1e-3], mol_mass= [1], configuration=Config)
def run_simulations(pamodel, substrate_rates):
    substrate_axis = list()
    Ccsc = list()
    Cesc = list()
    x_axis_csc = list()
    mu_list = list()

    for substrate in substrate_rates:
        pamodel.change_reaction_bounds(rxn_id='R1',
                                       lower_bound=0, upper_bound=substrate)
        pamodel.optimize()
        if pamodel.solver.status == 'optimal' and model.objective.value>0:
            print('Running simulations with ', substrate, 'mmol/g_cdw/h of substrate going into the system')
            substrate_axis += [substrate]
            mu_list += [pamodel.objective.value]

            Ccsc_new = list()
            for csc in ['flux_ub', 'flux_lb', 'enzyme_max', 'enzyme_min', 'proteome', 'sector']:
                Ccsc_new += pamodel.capacity_sensitivity_coefficients[pamodel.capacity_sensitivity_coefficients['constraint'] == csc].coefficient.to_list()
            Ccsc += [Ccsc_new]

            Cesc += [pamodel.enzyme_sensitivity_coefficients.coefficient.to_list()]

            print('Sum of capacity sensitivity coefficients: \t \t \t \t \t \t \t ', round(sum(Ccsc_new),6))
            print('Sum of variable sensitivity coefficients: \t \t \t \t \t \t \t ', round(sum(Cesc[-1]), 6), '\n')

    for csc in ['flux_ub', 'flux_lb', 'enzyme_max', 'enzyme_min', 'proteome', 'sector']:
        if csc == 'flux_ub' or csc == 'flux_lb':
            x_axis_csc += [rid +'_' + csc for rid in pamodel.capacity_sensitivity_coefficients[pamodel.capacity_sensitivity_coefficients['constraint'] == csc].rxn_id.to_list()]
        else:
            x_axis_csc += [rid +'_' + csc for rid in pamodel.capacity_sensitivity_coefficients[pamodel.capacity_sensitivity_coefficients['constraint'] == csc].enzyme_id.to_list()]

    x_axis_esc = pamodel.enzyme_sensitivity_coefficients.enzyme_id.to_list()

    return {'substrate_axis': substrate_axis, 'mu_list': mu_list,
            'Ccsc':Ccsc, 'Cesc':Cesc,
            'x_axis_csc': x_axis_csc,'x_axis_esc': x_axis_esc}

def print_heatmap(xaxis, matrix, yaxis = None):
    import plotly.express

    if yaxis is None:
        yaxis = list()
        for i in range(1, n + 1):
            yaxis += [f'R{i}']
    fig = plotly.express.imshow(matrix, aspect="auto",
                                x = xaxis, y = yaxis,
                                labels = dict(x = 'sensitivity coefficients', y='substrate uptake'))
    fig.show()

if __name__ == "__main__":
    model = build_toy_gem()
    active_enzyme = build_active_enzyme_sector(Config)
    unused_enzyme = build_unused_protein_sector(Config)
    translation_enzyme = build_translational_protein_sector(Config)
    pamodel = PAModel(model, name='toy model MCA with enzyme constraints', active_sector=active_enzyme,
                      translational_sector = translation_enzyme,
                      unused_sector = unused_enzyme, p_tot=Etot, configuration=Config)

    #optimize biomass formation
    pamodel.objective={pamodel.reactions.get_by_id('R7') :1}

    substrate_rates = np.arange(1e-3, 1e-1, 1e-3)
    simulation_results = run_simulations(pamodel, substrate_rates)


    print_heatmap(simulation_results['x_axis_csc'], simulation_results['Ccsc'], yaxis=simulation_results['substrate_axis'])
    print_heatmap(simulation_results['x_axis_esc'], simulation_results['Cesc'], yaxis=simulation_results['substrate_axis'])
