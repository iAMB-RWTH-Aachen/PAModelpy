import matplotlib.pyplot as plt
import pandas as pd
import os
import cobra

from Scripts.tam_generation import set_up_toy_tam, set_up_ecolicore_tam
from Scripts.pam_generation import set_up_ecolicore_pam

sinha_ref_conditions = {
    'Holm et al': ['REF', 'NOX', 'ATP'], #WT and NADH  overexpression conditions, mu 0.72, 0.65,0.58 h-1 respectively
    'Ishii et al': ['WT_0.7h-1'],
    'Gerosa et al.': ['Glycerol','Glucose','Acetate', 'Pyruvate','Gluconate','Succinate','Galactose','Fructose'],
}
TRANSCRIPT_FILE_PATH = os.path.join('Data', 'TAModel', 'Sinha-etal_2021_transcript-data.xlsx')
FLUX_FILE_PATH = os.path.join('Data', 'TAModel', 'Sinha-etal_2021_flux-data.xlsx')


mrna_vs_mu_slope = 0.00013049558330984208
mrna_vs_mu_intercept = 1.7750480089801658e-05


def get_transcript_data(transcript_file_path:str = TRANSCRIPT_FILE_PATH, mmol = True,
                        reference: str = 'Holm et al', growth_rates = [0.72, 0.65,0.58]):
    expression_data = pd.read_excel(transcript_file_path, sheet_name=reference, index_col=0)
    # Normalize columns by dividing by the sum of each column
    expression_data_normalized = expression_data.div(expression_data.sum(axis=0), axis=1)
    if mmol:
        mrna_conc = [mrna_vs_mu_intercept + mrna_vs_mu_slope*mu for mu in growth_rates]
        expression_data_mmol = expression_data_normalized.apply(lambda row: row * mrna_conc, axis=1)
        return expression_data_mmol
    else:
        return expression_data_normalized

def get_flux_data(flux_file_path:str = FLUX_FILE_PATH,
                        reference: str = 'Holm et al'):
    expression_data = pd.read_excel(flux_file_path, sheet_name=reference, index_col=0)
    if reference == 'Holm et al':
        #remove R suffix
        expression_data.index = expression_data.index.str.replace('R_', '')
    return expression_data

def get_pam_fluxes(substrate_uptake_rate):
    pam = set_up_ecolicore_pam()
    pam.change_reaction_bounds('EX_glc__D_e', lower_bound=-substrate_uptake_rate, upper_bound=0)
    sol = pam.optimize()
    pam_fluxes = sol.fluxes
    return pam_fluxes,pam

def set_up_tamodel(substrate_uptake_rate, strain ='REF'):
    tam = set_up_ecolicore_tam()
    tam.change_reaction_bounds('EX_glc__D_e', lower_bound=-substrate_uptake_rate, upper_bound=0)
    transcript_data_mmol = get_transcript_data()
    i=0
    infeasible = []
    for transcript in tam.transcripts:
        # print(enzyme, transcript)
        conc = []
        for gene in transcript.genes:
            if gene.id in transcript_data_mmol.index: #and gene.id != 'b2987' and gene.id!= 'b3493': #genes related to phosphate transport seem to mess up the calculations
                conc += [transcript_data_mmol.loc[gene.id, strain]]
        if len(conc) > 0:
            transcript.change_concentration(concentration=min(conc) * 1e-3, error= max(conc)*1e-3* 0.01)

            tam.change_reaction_bounds('EX_glc__D_e', lower_bound=-9, upper_bound=-9)
            tam.optimize()
            if tam.solver.status != 'optimal': infeasible += [transcript]
            else: print(tam.objective.value)

        elif gene.id != 'gene_dummy':
            print(transcript)

    print(infeasible)
            # tam.test()
            # print(transcript.enzymes)
            # print(tam.solver.status)
    # for gene, expression_data in transcript_data_mmol.iterrows():
    #     if i<10000:
    #         transcript_id = 'mRNA_' + gene
    #         if not transcript_id in tam.transcripts: continue
    #
    #         transcript = tam.transcripts.get_by_id('mRNA_' + gene)
    #         # testing wildtype condition
    #         transcript.change_concentration(concentration=expression_data[strain]*1e-3,
    #                                             error=expression_data[strain] *1e-3* 0.01)
    #         i+=1
    #         tam.test()
    #         print(transcript.enzymes)
    #         print(tam.solver.status)
    #         print(transcript)
    #         if tam.solver.status != 'optimal':
    #             print('number of transcripts which were be added to the model: ',i)
    #             for name, cons in transcript._constraints.items():
    #                 print(name, cons, expression_data[strain])
    #             break
    # print('number of transcripts which could be added to the model: ', len(tam.transcripts))
    # enz = tam.enzymes.get_by_id('2.7.3.9')
    # print(enz.name, enz.enzyme_variable.concentration)
    # print(infeasible)
    return tam

def get_tam_fluxes(tam,substrate_uptake_rate):
    tam.change_reaction_bounds('EX_glc__D_e', lower_bound=-substrate_uptake_rate, upper_bound=0)
    sol = tam.optimize()
    if tam.solver.status != 'optimal': return None
    for transcript in tam.transcripts:
        for constr in transcript._constraints.values():
            if constr.dual !=0:
                print(transcript, transcript.mrna_variable.primal, constr.primal)
                print(transcript.f_min*transcript.mrna_variable.primal, transcript.f_max*transcript.mrna_variable.primal)
                for enzyme in transcript.enzymes: print(enzyme.id, enzyme.concentration)
                print(constr)
    tam_fluxes = sol.fluxes
    return tam_fluxes

def compare_flux_data(flux_data, pam_fluxes, tam_fluxes, strain ='REF', abs=True):
    if abs:
        glc_upt_ref = 1
        glc_upt_pam = 1
        glc_upt_tam = 1
    else:
        glc_upt_ref = flux_data[strain]['GLCptspp']
        glc_upt_pam = pam_fluxes['GLCpts']
        glc_upt_tam = tam_fluxes['GLCpts']

    flux_results = flux_data[[strain]]
    flux_results_percentage = flux_results.assign(
        strain=lambda val: val[strain] / glc_upt_ref)
    flux_results_percentage['PAM'] = 0
    flux_results_percentage['TAM'] = 0
    for rxn in flux_data.index:
        ori_rxn = rxn
        if 'pp' in rxn: rxn = rxn.replace('pp', '')
        if 'biomass' in rxn: rxn = 'BIOMASS_Ecoli_core_w_GAM'
        if 'EX_glc' in rxn: rxn = 'EX_glc__D_e'
        if 'EX_ac' in rxn: rxn = 'EX_ac_e'
        flux_results_percentage['PAM'][ori_rxn] = pam_fluxes[rxn] / glc_upt_pam
        flux_results_percentage['TAM'][ori_rxn] = tam_fluxes[rxn] / glc_upt_tam

    print(flux_results_percentage.to_markdown())
    return flux_results_percentage

def compare_fluxes_holm_reference(strain = 'REF', plot = True):
    flux_data =get_flux_data()
    substrate_uptake_rate = flux_data[strain]['GLCptspp']

    pam_fluxes, pam = get_pam_fluxes(substrate_uptake_rate=substrate_uptake_rate)

    tam = set_up_tamodel(substrate_uptake_rate,strain)
    tam_fluxes = get_tam_fluxes(tam, substrate_uptake_rate=substrate_uptake_rate)
    if tam_fluxes is None:
        print('TAModel was infeasible')
        return
    for i,row in tam.capacity_sensitivity_coefficients.iterrows():
        if row.coefficient > 0: print(row)

    # for transcript in tam.transcripts:
    #     if transcript.mrna_variable.lb != 0:
    #         for gene in transcript.genes:
    #             gene_id = gene.id
    #             enzymes = tam.get_enzymes_by_gene_id(gene_id)
    #             for enzyme in enzymes:
    #                 print(enzyme.id, 'TAM',enzyme.concentration,'PAM' ,pam.enzymes.get_by_id(enzyme.id).concentration)
    #                 print()
    # print('RNA constraints primal', tam.constraints[tam.TOTAL_MRNA_CONSTRAINT_ID].primal)
    # print('RNA constraints dual', tam.constraints[tam.TOTAL_MRNA_CONSTRAINT_ID].dual)
    # for enzyme in pam.enzyme_variables[:9]:
    #     tam_enzyme = tam.enzyme_variables.get_by_id(enzyme.id)
    #
    #     print(enzyme.id, enzyme.concentration, tam_enzyme.concentration)
    #     tam_transcripts = tam.get_transcripts_associated_with_enzyme(tam.enzymes.get_by_id(enzyme.id), output='list')
    #     for transcript in tam_transcripts:
    #         print(transcript.id, transcript.mrna_variable.primal)
    #         print('bounds',transcript.mrna_variable.lb, transcript.mrna_variable.ub)
    #         for constraint in transcript._constraints.values():
    #             if constraint.dual>0: print(constraint, constraint.dual)


    flux_relative = compare_flux_data(flux_data, pam_fluxes, tam_fluxes,  strain,abs = False)
    flux_absolute = compare_flux_data(flux_data, pam_fluxes, tam_fluxes,  strain)

    if plot: plot_flux_comparison(flux_absolute, flux_relative, strain)

def plot_flux_comparison(flux_df_abs, flux_df_rel, strain):
    fig, ax = plt.subplots(1,2)

    ax[0].scatter(flux_df_abs['TAM'], flux_df_abs['strain'], color = 'black')
    ax[0].scatter(flux_df_abs['PAM'], flux_df_abs['strain'], color ='red')
    #reference line
    ax[0].plot(flux_df_abs['strain'], flux_df_abs['strain'], linestyle ='dashed')

    ax[0].set_title(strain + ' absolute fluxes')
    ax[0].set_xlabel('simulated flux [$mmol/g_{CDW}/h$]')
    ax[0].set_ylabel('measured flux [$mmol/g_{CDW}/h$]')

    ax[1].scatter(flux_df_rel['TAM'], flux_df_rel['strain'], color='black', label = 'TAM')
    ax[1].scatter(flux_df_rel['PAM'], flux_df_rel['strain'], color='red', label = 'PAM')
    # reference line
    ax[1].plot(flux_df_rel['strain'], flux_df_rel['strain'], linestyle='dashed')

    ax[1].set_title(strain + ' relative fluxes')
    ax[1].set_xlabel('simulated flux [$mmol/g_{CDW}/h$]')
    ax[1].set_ylabel('measured flux [$mmol/g_{CDW}/h$]')

    fig.set_figwidth(20)
    fig.set_figheight(10)
    plt.legend()
    plt.show()


if __name__ == '__main__':

    print('Reference condition')
    compare_fluxes_holm_reference(plot =False)
    print('\n-------------------------------------------------------------------------------------------------')
    # print('mutation 1: NOX strain (overexpression of NADH oxidase)\n')
    # compare_fluxes_holm_reference('NOX', plot=False)
    # TODO print mRNA and protein concentrations to compare with lb
    # TODO print shadowprices of mRNA (are lbs hit? how far can I constrain?)






