from cobra import Configuration
from cobra import Model, Reaction, Metabolite
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
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
    r2.add_metabolites({Metabolite('Substrate'): -1, Metabolite('Intermediate'): 1, Metabolite('ATP'): 0.5})
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
    kcat_fwd = [1, 0.5, 5, 0.1, 0.25 ,0.45, 1.5]  # High turnover for exchange reactions [1, 0.5, 1, 1, 0.5 ,0.45, 1.5]
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
    Ccsc = list()
    Cesc = list()
    x_axis_csc = list()
    mu_list = list()
    substrate_consumption = list()
    byproduct_formation = list()
    total_active_protein = list()

    for substrate in substrate_rates:
        pamodel.change_reaction_bounds(rxn_id='R1',
                                       lower_bound=0, upper_bound=substrate)
        pamodel.optimize()
        if pamodel.solver.status == 'optimal' and model.objective.value>0:
            print('Running simulations with ', substrate, 'mmol/g_cdw/h of substrate going into the system')
            substrate_axis += [substrate]
            mu_list += [pamodel.objective.value]
            substrate_consumption += [pamodel.reactions.R1.flux]
            byproduct_formation += [pamodel.reactions.R3.flux]
            #phi_active = total_protein - translational_protein - unused_enzymes
            total_active_protein += [Etot*1e3 - (0.01 + 0.01*pamodel.reactions.R7.flux) - (0.1-0.01*pamodel.reactions.R1.flux)]

            Ccsc_new = list()
            for csc in ['flux_ub', 'flux_lb', 'enzyme_max', 'enzyme_min', 'proteome', 'sector']:
                Ccsc_new += pamodel.capacity_sensitivity_coefficients[pamodel.capacity_sensitivity_coefficients['constraint'] == csc].coefficient.to_list()
            Ccsc += [Ccsc_new]

            Cesc += [pamodel.enzyme_sensitivity_coefficients.coefficient.to_list()]

            print('Sum of capacity sensitivity coefficients: \t \t \t \t \t \t \t ', round(sum(Ccsc_new),6))
            print('Sum of enzyme sensitivity coefficients: \t \t \t \t \t \t \t ', round(sum(Cesc[-1]), 6), '\n')

    for csc in ['flux_ub', 'flux_lb', 'enzyme_max', 'enzyme_min', 'proteome', 'sector']:
        if csc == 'flux_ub' or csc == 'flux_lb':
            x_axis_csc += [rid +'_' + csc for rid in pamodel.capacity_sensitivity_coefficients[pamodel.capacity_sensitivity_coefficients['constraint'] == csc].rxn_id.to_list()]
        else:
            x_axis_csc += [rid +'_' + csc for rid in pamodel.capacity_sensitivity_coefficients[pamodel.capacity_sensitivity_coefficients['constraint'] == csc].enzyme_id.to_list()]

    x_axis_esc = pamodel.enzyme_sensitivity_coefficients.enzyme_id.to_list()

    return {'substrate_axis': substrate_axis, 'mu_list': mu_list,
            'byproduct_formation': byproduct_formation, 'substrate_consumption':substrate_consumption,
            'total_active_protein': total_active_protein,
            'Ccsc':Ccsc, 'Cesc':Cesc,
            'x_axis_csc': x_axis_csc,'x_axis_esc': x_axis_esc}

def plot_sensitivities(fig, grdspec, glc_rates, mu_list, substrate, byproduct, total_active_protein,
                       tot_prot_csc, substrate_csc, e1_esc, pie_data, pie_labels):
    gs = gridspec.GridSpecFromSubplotSpec(2, 1,
                                          height_ratios=[1,1], hspace=0, subplot_spec=grdspec)
    sens_ax = fig.add_subplot(gs[1, 0])  # sensitivity coefficient linegraph
    mu_ax = fig.add_subplot(gs[0, 0], sharex=sens_ax) # mu vs v_s linegraph

    mu = mu_ax.plot(glc_rates, mu_list, color = 'black', linewidth= 3)#(35/255,158/255,137/255)
    sub = mu_ax.plot(glc_rates, substrate, color='orange', linewidth=3)
    byp = mu_ax.plot(glc_rates, byproduct, color='purple', linewidth=3)
    prot_ax = mu_ax.twinx()
    etotact = prot_ax.plot(glc_rates, [p*10 for p in total_active_protein], color = 'darkred', linewidth = 3, linestyle = 'dotted')
    prot_ax.set_ylabel('active enzyme sector (AES) ($g/g_{protein}$)')
    prot_ax.set_ylim([4.9,4.9022])
    prot_ax.annotate('$\cdot 10^{-1}$',xy=(0.98, 0.1), xycoords='figure fraction',
                   xytext=(0.98, 1.1), textcoords='axes fraction',
                   va='top', ha='left', fontsize=FONTSIZE*0.65)
    prot_ax.set_yticks([4.9,4.901,4.902])
    mu_ax.legend([sub, byp, mu, etotact], labels = ['R1', 'R9', 'R7', 'AES'], handles=sub+byp+mu+etotact,loc = 'upper left')
    mu_ax.xaxis.set_visible(False)
    # mu_ax.legend([mu], labels=['growth rate'], loc='center left')
    mu_ax.set_ylabel('flux rate', fontsize = FONTSIZE*3/4)
    # add B panel annotation
    mu_ax.annotate('B', xy=(-0.1, 0.1), xycoords='figure fraction',
                   xytext=(-0.1, 1.1), textcoords='axes fraction',
                   va='top', ha='left', fontsize=FONTSIZE*1.5, weight='bold')
    mu_ax.set_ylim([-0.001,0.07])

    # plot the sensitivity coefficients
    vs = sens_ax.plot(glc_rates, substrate_csc, color ='orange', linewidth= 3) #(62/255,174/255, 137/255)
    e1 = sens_ax.plot(glc_rates, e1_esc, color =(0.00,0.45,0.81), linewidth= 3, linestyle ='dashed')
    e_tot = sens_ax.plot(glc_rates, tot_prot_csc, color ='darkred', linewidth= 3, linestyle ='dotted')#(62/255, 64/255, 137/255)

    sens_ax.legend([vs, e1, e_tot], labels=['$FCSC_{R1}$','$ESC_{E_{1}}$','$PCSC$'], loc='center left')
    sens_ax.set_ylim([-0.01,1.1])
    sens_ax.set_xlim([0,0.1])
    sens_ax.set_ylabel('Sensitivity Coefficients', fontsize = FONTSIZE*3/4)
    sens_ax.set_xlabel('$v_{substrate,max}$ $(mmol_{substrate}/g_{CDW}/h)$', fontsize = FONTSIZE*3/4)

    # Add a pie chart showing the enzyme sensitivities as an inset in the sensitivity linegraph
    inset_ax = inset_axes(sens_ax, width="100%", height="100%",  loc=3, bbox_to_anchor=(0.5, 0.15, 0.7, 0.7),  bbox_transform=sens_ax.transAxes) #TODO
    inset_ax.pie(pie_data, labels=pie_labels,  textprops={'fontweight': 'bold'}, startangle=90, labeldistance = 0.6)

    return fig

if __name__ == "__main__":
    #NOTE:RUN THIS SCRIPT FROM THE MAIN DIRECTORY IN THE COMMAND LINE.
    # otherwise the relative paths set in this script will result in errors
    FONTSIZE = 16
    width = 10
    height = 12
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

    gs0 = gridspec.GridSpec(2, 1, figure=fig, height_ratios=[3,2])
    gs_toymodel = gs0[0]
    gs_sensitivities = gs0[1]

    image_path = 'Figures/Figure1_toy-model_overflow.png'
    toy_model = np.asarray(Image.open(image_path))
    ax_fig = fig.add_subplot(gs_toymodel)
    ax_fig.imshow(toy_model, aspect='equal')
    ax_fig.annotate('A', xy=(2, 1), xycoords='data',
                    xytext=(-0.1, 1), textcoords='axes fraction',
                    va='top', ha='left', fontsize=FONTSIZE*1.5, weight='bold')
    ax_fig.axis('off')
    ax_fig.set_xticks([])
    ax_fig.set_yticks([])

    for index, id in enumerate(simulation_results['x_axis_csc']):
        if 'TotalProteinConstraint_proteome' in id:
            tot_prot_csc = [row[index] for row in simulation_results['Ccsc']]
        if 'R1_flux_ub' in id:
            substrate_csc = [row[index] for row in simulation_results['Ccsc']]

    pielabels = []
    for index, id in enumerate(simulation_results['x_axis_esc']):
        if 'E1' in id:
            e1_esc = [row[index] for row in simulation_results['Cesc']]
        if 'E4' in id:
            pielabels.append('')
        else:
            pielabels.append(id)

    fig = plot_sensitivities(fig, gs_sensitivities, simulation_results['substrate_axis'], simulation_results['mu_list'],
                             simulation_results['substrate_consumption'],simulation_results['byproduct_formation'],
                             simulation_results['total_active_protein'],
                             tot_prot_csc, substrate_csc, e1_esc, simulation_results['Cesc'][-1], pielabels)
    fig.set_figwidth(width)
    fig.set_figheight(height)
    plt.tight_layout()

    fig.savefig('Figures/Figure1_toy_model-sensitivities_overflow.png',bbox_inches='tight', dpi=600)
    plt.show()