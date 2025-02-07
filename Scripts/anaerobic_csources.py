import pandas as pd

from Scripts.pam_generation import set_up_ecoli_pam

if __name__ == "__main__":
    condition2uptake = {'Glycerol': 'EX_gly_e', 'Glucose': 'EX_glc__D_e', 'Acetate': 'EX_ac_e', 'Pyruvate': 'EX_pyr_e',
                        'Gluconate': 'EX_glcn_e', 'Succinate': 'EX_succ_e', 'Galactose': 'EX_gal_e',
                        'Fructose': 'EX_fru_e', 'Propionate':"EX_ppa_e"}

    output_df = pd.DataFrame(columns=['carbon_source', 'substrate_uptake', 'growth_rate', 'glycerol'])

    pamodel = set_up_ecoli_pam(sensitivity=False)
    pamodel.change_reaction_bounds("EX_glc__D_e", 0,0)
    pamodel.change_reaction_bounds("EX_o2_e", 0,0)
    pamodel.change_reaction_bounds("EX_gly_e", -1e6, 0)

    for csource, uptake_id in condition2uptake.items():
        pamodel.change_reaction_bounds(uptake_id, -1e6, 0)
        pamodel.optimize()
        if pamodel.solver.status == 'optimal':
            output_df.loc[len(output_df)] = [csource, pamodel.reactions.get_by_id(uptake_id).flux,
                                             pamodel.objective.value, pamodel.reactions.get_by_id('EX_gly_e').flux,]
        else:
            output_df.loc[len(output_df)] = [csource, 'inf', 'inf', 'inf']
        if csource != 'Glycerol':
            pamodel.change_reaction_bounds(uptake_id, 0, 0)

    output_df.to_excel('Results/anaerobic_csources_pam_coglyc.xlsx')

