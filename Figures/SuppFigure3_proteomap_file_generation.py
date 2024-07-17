import numpy as np
import pandas as pd
from Scripts.pam_generation import set_up_ecoli_pam


# NOTE:RUN THIS SCRIPT FROM THE MAIN DIRECTORY IN THE COMMAND LINE.
# otherwise the relative paths set in this script will result in errors

def calculate_sensitivities(pamodel):
    glc =10

    # disable pyruvate formate lyase (inhibited by oxygen)
    pamodel.change_reaction_bounds(rxn_id='PFL', upper_bound=0)

    print('glucose uptake rate ', glc, ' mmol/gcdw/h')
    with pamodel:
        # change glucose uptake rate
        pamodel.change_reaction_bounds(rxn_id='EX_glc__D_e',
                                           lower_bound=-glc, upper_bound=-glc)
        # solve the model
        pamodel.optimize()
        if pamodel.solver.status == 'optimal':
            genes = []
            enzyme_coeff = pamodel.enzyme_sensitivity_coefficients
            Cesc = enzyme_coeff.coefficient.to_list()

            #get the first gene from the first reaction as this all the reactions and genes in one row relate to the same enzyme
            for i,rxn in enumerate(enzyme_coeff.rxn_id.to_list()):
                rxn_genes = pamodel.reactions.get_by_id(rxn.split(',')[0])._genes
                if all(gene.id in genes for gene in rxn_genes):
                    gene = rxn_genes.pop()
                    genes.append(gene.id)

                else:
                    for gene in rxn_genes:
                        if gene.id not in genes:
                            genes.append(gene.id)
                            break


            print('Sum of enzyme sensitivity coefficients: \t \t \t \t \t \t', round(sum(Cesc), 6), '\n')
    return { 'genes': genes,'Cesc': Cesc}

if __name__ == '__main__':
    result = calculate_sensitivities(set_up_ecoli_pam())
    result_df = pd.DataFrame(result)
    result_df = result_df.groupby('genes', as_index=False)['Cesc'].sum()

    result_df.to_csv('Results/sensitivity_coefficient_proteomap_ecoli_vs10.tsv',sep = '\t',
                     header=False, index = False)
