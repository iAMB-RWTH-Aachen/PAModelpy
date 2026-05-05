import re
from src.PAModelpy.utils.pam_generation import set_up_pam
from Scripts.create_pamodel_from_diagnostics_file import create_pamodel_from_diagnostics_file

def parse_carbon_count(formula: str) -> int:
    """
    Extract number of carbon atoms from a chemical formula string.
    
    օրինակ:
    C6H12O6 -> 6
    CO2     -> 1
    CH4     -> 1
    H2O     -> 0
    """
    if formula is None:
        return 0

    match = re.search(r'C(?![a-z])(\d*)', formula)
    
    if match:
        # If it's just "C" with no number → count = 1
        return int(match.group(1)) if match.group(1) else 1
    
    return 0

def check_carbon_balance_reaction(rxn):
    carbon_balance = 0
    
    for met, coeff in rxn.metabolites.items():
        nC = parse_carbon_count(met.formula)
        carbon_balance += coeff * nC
    
    return carbon_balance

pam_info_file = "Data/proteinAllocationModel_EnzymaticData_iML1515_10.xlsx"
pamodel = set_up_pam(pam_info_file=pam_info_file, sensitivity=False)
solution = pamodel.optimize()
carbon_balance = 0
imbalanced_reactions = []

for rxn in pamodel.exchanges:
    flux = solution.fluxes[rxn.id]
    for met, coeff in rxn.metabolites.items():
        nC = parse_carbon_count(met.formula)
        print(rxn, met.formula, nC)
        carbon_balance += flux * coeff * nC
    
    formula = met.formula  # e.g. C6H12O6
    nC = parse_carbon_count(formula)
    
    carbon_balance += flux * nC

biomass_rxn = pamodel.reactions.get_by_id("BIOMASS_Ec_iML1515_core_75p37M")  # adjust ID

biomass_flux = solution.fluxes[biomass_rxn.id]

biomass_carbon = 0
for met, coeff in biomass_rxn.metabolites.items():
    nC = parse_carbon_count(met.formula)
    biomass_carbon += flux * coeff * nC

carbon_balance += biomass_carbon

print(carbon_balance)