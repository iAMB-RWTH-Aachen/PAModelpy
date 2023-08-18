from cobra import Configuration
from cobra import Model, Reaction, Metabolite
import plotly.express

#importing the tools from the PAModelpy package
from PAModelpy.EnzymeSectors import ActiveEnzymeSector, UnusedEnzymeSector, TransEnzymeSector
from PAModelpy.PAModel import PAModel

# from PAModelpy import ActiveEnzymeSector, UnusedEnzymeSector, PAModel
#need to have gurobipy installed

#global variables:
global metabolites, n, m, Etot
#global variables:
global metabolites, n, m, Etot
metabolites = ['Substrate', 'ATP', 'CO2', 'Precursor', 'Biomass', 'Byproduct', 'Intermediate']
n = 9
m = 7
Etot = 0.006

#functions:
def build_toy_gem():
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

def build_active_enzyme_sector():
    kcat_fwd = [1, 0.5, 1, 0.5 ,0.45, 1.5]  # High turnover for exchange reactions
    kcat_rev = [kcat for kcat in kcat_fwd]
    rxn2kcat = {}
    for i in range(n-3): # all reactions have an enzyme, except excretion reactions
        rxn_id = f'R{i+1}'
        # 1e-6 to correct for the unit transformation in the model (meant to make the calculations preciser for different unit dimensions)
        #dummy molmass like in MATLAB script
        rxn2kcat = {**rxn2kcat, **{rxn_id: {f'E{i+1}':{'f': kcat_fwd[i]/(3600*1e-6), 'b': kcat_rev[i]/(3600*1e-6), 'molmass': 1e6}}}}

    return ActiveEnzymeSector(rxn2protein = rxn2kcat)

def build_unused_protein_sector():
    return UnusedEnzymeSector(id_list = ['R1'], ups_mu=[0.01*1e-3], ups_0=[0.1*1e-3], mol_mass= [1])

def build_translational_protein_sector():
    return TransEnzymeSector(id_list = ['R1'], tps_mu=[0.01*1e-3], tps_0=[0.01*1e-3], mol_mass= [1])

def print_heatmap(xaxis, matrix):
    yaxis = list()
    for i in range(1, n + 1):
        yaxis += [f'R{i}']
    fig = plotly.express.imshow(matrix, aspect="auto",
                                x = xaxis, y = yaxis,
                                labels = dict(x = 'control coefficients', y='maximized reaction'))
    fig.show()

if __name__ == "__main__":
    Ccac=list()
    Cfac = list()
    x_axis_cac = list()
    x_axis_fac = list()
    model = build_toy_gem()
    active_enzyme = build_active_enzyme_sector()
    unused_enzyme = build_unused_protein_sector()
    translation_enzyme = build_translational_protein_sector()
    pamodel = PAModel(model, name='toy model MCA with enzyme constraints', active_sector=active_enzyme,
                      translational_sector = translation_enzyme,
                      unused_sector = unused_enzyme, p_tot=Etot)


    for rxn in pamodel.reactions:
        #only get the network reactions, not the upperbound and lowerbound pseudo-reactions
        if 'b' not in rxn.id:
            pamodel.objective = {rxn: 1}
            pamodel.optimize()

            Ccac_new = list()
            for cac in ['UB', 'LB', 'EC_f', 'EC_b', 'sector']:
                Ccac_new += pamodel.capacity_allocation_coefficients[pamodel.capacity_allocation_coefficients['constraint'] == cac].coefficient.to_list()
            Ccac += [Ccac_new]

            Cfac_new = list()
            for fac in ['rxn', 'enzyme', 'sector']:
                Cfac_new += pamodel.flux_allocation_coefficients[pamodel.flux_allocation_coefficients['constraint'] == fac].coefficient.to_list()
            Cfac += [Cfac_new]

        print('Sum of control coefficients: \t \t \t \t \t \t \t \t', round(sum(Ccac_new),6))
        print('Sum of allocation coefficients: \t \t \t \t \t \t \t', round(sum(Cfac_new), 6), '\n')

    for cac in ['UB', 'LB', 'EC_f', 'EC_b', 'sector']:
        if cac == 'UB' or cac == 'LB':
            x_axis_cac += [rid+'_'+cac for rid in pamodel.capacity_allocation_coefficients[pamodel.capacity_allocation_coefficients['constraint'] == cac].rxn_id.to_list()]
        else:
            x_axis_cac += [rid+'_'+cac[-1] for rid in pamodel.capacity_allocation_coefficients[pamodel.capacity_allocation_coefficients['constraint'] == cac].enzyme_id.to_list()]

    for fac in ['rxn', 'enzyme', 'sector']:
        if fac == 'rxn':
            x_axis_fac += pamodel.flux_allocation_coefficients[pamodel.flux_allocation_coefficients['constraint'] == fac].rxn_id.to_list()
        else:
            x_axis_fac += pamodel.flux_allocation_coefficients[
                pamodel.flux_allocation_coefficients['constraint'] == fac].enzyme_id.to_list()

    print_heatmap(x_axis_cac, Ccac)
    print_heatmap(x_axis_fac, Cfac)
