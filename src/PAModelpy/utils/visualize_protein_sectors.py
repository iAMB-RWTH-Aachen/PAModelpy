import matplotlib.pyplot as plt
from typing import Iterable

from ..PAModel import PAModel

def run_simulations(pamodel:PAModel,
                    substrate_rates:Iterable[float],
                    sub_uptake_id = 'EX_glc__D_e'
                    ) -> list:
    fluxes = []
    total_enzyme_concentration = []
    for substrate in substrate_rates:
        if substrate<0:
            pamodel.change_reaction_bounds(rxn_id=sub_uptake_id,
                                       lower_bound=substrate, upper_bound=0)
        else:
            pamodel.change_reaction_bounds(rxn_id=sub_uptake_id,
                                       lower_bound=0, upper_bound=substrate)

        print('Running simulations with ', substrate, 'mmol/g_cdw/h of substrate going into the system')
        sol_pam =pamodel.optimize()
        tot_prot_conc = 0
        if pamodel.solver.status == 'optimal' and pamodel.objective.value>0:
            fluxes.append(sol_pam.fluxes)
            # save sum of enzymes
            enzyme_variables = [
                var for var in pamodel.constraints[pamodel.TOTAL_PROTEIN_CONSTRAINT_ID].variables
                if not any([
                sector.id_list[0] in var.name for sector in pamodel.sectors if sector.id!='ActiveEnzymeSector'
            ])
                          ]
            for var, coef in pamodel.constraints[pamodel.TOTAL_PROTEIN_CONSTRAINT_ID].get_linear_coefficients(enzyme_variables).items():
                tot_prot_conc += coef*var.primal *1e-3
        total_enzyme_concentration.append(tot_prot_conc)
    return fluxes, total_enzyme_concentration

def visualize_protein_sectors(pamodel:PAModel,
                              substrate_rates:Iterable[float],
                              sub_uptake_id = 'EX_glc__D_e',
                              ax:plt.Axes = None
                              ) -> plt.Axes:
    simulated_fluxes, total_enzyme_concentration = run_simulations(
        pamodel = pamodel,
        substrate_rates = substrate_rates,
        sub_uptake_id = sub_uptake_id
    )
    sector_protein_fraction = {'ActiveEnzymeSector': total_enzyme_concentration}
    for sector in pamodel.sectors:
        if sector.id in sector_protein_fraction: continue
        sector_protein_fraction[sector.id] = [
            fluxes[sector.id_list[0]]*sector.slope*1e-3 + sector.intercept*1e-3 for fluxes in simulated_fluxes
        ]
    sector_protein_fraction['total'] = [sum([k[i] for k in sector_protein_fraction.values()]) for i in range(len(substrate_rates))]

    if ax is None:
        fig, ax = plt.subplot()

    ax.hlines(y=pamodel.p_tot,xmin=0,xmax=max([abs(s) for s in substrate_rates]),
              label = 'metabolic protein fraction',
              color = 'black')
    ax.set_xlabel(sub_uptake_id)
    ax.set_ylabel(r'protein fraction [$\text{g}_{\text{p}}/\text{g}_{\text{CDW}}$]')

    for sector_label, sector_fractions in sector_protein_fraction.items():
        ax.plot([abs(s) for s in substrate_rates], sector_fractions, label = sector_label)

    ax.legend()

    return ax




