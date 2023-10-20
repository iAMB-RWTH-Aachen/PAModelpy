from cobra import Configuration
from cobra import Model, Reaction, Metabolite
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from PIL import Image
Image.MAX_IMAGE_PIXELS = None #to make sure the image is loaded properly in high quality

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

FONTSIZE = 16

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

def run_simulations(pamodel, substrate_axis):
    substrate_axis = list()
    Ccac = list()
    Cfac = list()
    x_axis_cac = list()
    x_axis_fac = list()
    mu_list = list()

    for substrate in substrate_rates:
        pamodel.change_reaction_bounds(rxn_id='R1',
                                       lower_bound=0, upper_bound=substrate)
        pamodel.optimize()
        if pamodel.solver.status == 'optimal' and model.objective.value>0:
            print('Running simulations with ', substrate, 'mmol/g_cdw/h of substrate going into the system')
            substrate_axis += [substrate]
            mu_list += [pamodel.objective.value]
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

    return {'substrate_axis': substrate_axis, 'mu_list': mu_list,
            'Ccac':Ccac, 'Cfac':Cfac,
            'x_axis_cac': x_axis_cac,'x_axis_fac': x_axis_fac}

def plot_sensitivities(fig, grdspec, glc_rates, mu_list, tot_prot_cac, substrate_cac, e1_fac, substrate_fac):
    gs = gridspec.GridSpecFromSubplotSpec(3, 1,
                                          height_ratios=[1,1,1], hspace=0, subplot_spec=grdspec)
    fac_ax = fig.add_subplot(gs[2, 0])  # FAC linegraph
    mu_ax = fig.add_subplot(gs[0, 0], sharex=fac_ax) # mu vs v_s linegraph
    cac_ax = fig.add_subplot(gs[1, 0], sharex=fac_ax)  # CAC linegraph



    mu_ax.plot(glc_rates, mu_list, color = 'blue', linewidth= 3)#(35/255,158/255,137/255)
    mu_ax.xaxis.set_visible(False)
    mu_ax.set_ylabel('$v_{biomass} (h^{-1})$', fontsize = FONTSIZE)
    # add B panel annotation
    mu_ax.annotate('B', xy=(0.01, 0.01), xycoords='data',
                   xytext=(-0.10, 1.2), textcoords='axes fraction',
                   va='top', ha='left', fontsize=FONTSIZE, weight='bold')

    e_tot = cac_ax.plot(glc_rates, tot_prot_cac, color = 'darkred', linewidth= 3)#(62/255, 64/255, 137/255)
    vs = cac_ax.plot(glc_rates, substrate_cac, color ='gray', linewidth= 3) #(62/255,174/255, 137/255)
    cac_ax.legend([e_tot, vs], labels=['$\phi_{E,0}$','$v_{1}$'], loc='center left')
    cac_ax.xaxis.set_visible(False)
    cac_ax.set_ylim([-0.1,1.3])
    cac_ax.set_ylabel('CAC', fontsize = FONTSIZE)


    e1 = fac_ax.plot(glc_rates, e1_fac, color =(68/255,1/255,84/255), linewidth= 3)
    vs = fac_ax.plot(glc_rates, substrate_fac, color= 'gray', linewidth= 3)#(62/255,174/255, 137/255)
    fac_ax.legend([e1, vs], labels=['$E_{1}$', '$v_{1}$'], loc='center left')
    fac_ax.set_ylim([-0.1,1.3])
    fac_ax.set_ylabel('FAC', fontsize = FONTSIZE)
    fac_ax.set_xlabel('$v_{substrate,max} (mmol_{substrate}/g_{CDW}/h)$', fontsize = FONTSIZE)

    return fig

if __name__ == "__main__":
    width = 13
    height =10
    model = build_toy_gem()
    active_enzyme = build_active_enzyme_sector(Config)
    unused_enzyme = build_unused_protein_sector(Config)
    translation_enzyme = build_translational_protein_sector(Config)
    pamodel = PAModel(model, name='toy model MCA with enzyme constraints', active_sector=active_enzyme,
                      translational_sector = translation_enzyme,
                      unused_sector = unused_enzyme, p_tot=Etot,
                      configuration = Config)
    #optimize biomass formation
    pamodel.objective={pamodel.reactions.get_by_id('R7') :1}
    #over a range of substrate uptake rates
    substrate_rates = np.arange(1e-3, 1e-1, 1e-3)

    simulation_results = run_simulations(pamodel, substrate_rates)

    #plot figure with multiple pannels
    # gridspec inside gridspec
    fig = plt.figure()

    gs0 = gridspec.GridSpec(1, 18, figure=fig, wspace = 25)
    gs_toymodel = gs0[:10]
    gs_sensitivities = gs0[11:]

    image_path = 'Figure1_toy-model.png'
    toy_model = np.asarray(Image.open(image_path))
    ax_fig = fig.add_subplot(gs_toymodel)
    ax_fig.imshow(toy_model)
    ax_fig.annotate('A', xy=(2, 1), xycoords='data',
                    xytext=(-0.1, 1.1), textcoords='axes fraction',
                    va='top', ha='left', fontsize=FONTSIZE, weight='bold')
    ax_fig.axis('off')
    ax_fig.set_xticks([])
    ax_fig.set_yticks([])

    for index, id in enumerate(simulation_results['x_axis_cac']):
        if 'TotalProteinConstraint' in id:
            tot_prot_cac = [row[index] for row in simulation_results['Ccac']]
        if 'R1_UB' in id:
            substrate_cac = [row[index] for row in simulation_results['Ccac']]

    for index, id in enumerate(simulation_results['x_axis_fac']):
        if 'E1' in id:
            e1_fac = [row[index] for row in simulation_results['Cfac']]

        if 'R1' in id:
            substrate_fac = [row[index] for row in simulation_results['Cfac']]


    fig = plot_sensitivities(fig, gs_sensitivities, simulation_results['substrate_axis'], simulation_results['mu_list'],
                             tot_prot_cac, substrate_cac, e1_fac, substrate_fac)
    fig.set_figwidth(width)
    fig.set_figheight(height)

    fig.savefig('Figure1_toy_model-sensitivities.png')
    plt.show()