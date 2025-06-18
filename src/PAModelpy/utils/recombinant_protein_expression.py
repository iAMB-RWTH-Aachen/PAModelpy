from typing import Union, Dict, Optional
from cobra import Model as CobraModel
from cobra import Reaction, Metabolite

from .. import PAModel, Enzyme, EnzymeVariable

DEFAULT_PROTEIN_MOLMASS = 39959.4825 #Da
REFERENCE_GROWTH_RATE = 0.1 #h-1, used to estimate the ribosome utilization per protein being secreted

# aminoacid lookup table
aa_lookup = {'V': 'VAL', 'I': 'ILE', 'L': 'LEU', 'E': 'GLU', 'Q': 'GLN', \
             'D': 'ASP', 'N': 'ASN', 'H': 'HIS', 'W': 'TRP', 'F': 'PHE', 'Y': 'TYR', \
             'R': 'ARG', 'K': 'LYS', 'S': 'SER', 'T': 'THR', 'M': 'MET', 'A': 'ALA', \
             'G': 'GLY', 'P': 'PRO', 'C': 'CYS'}


def check_freq(x:str, total:int) -> Dict[str, Union[float, int]]:
    return {c: x.count(c) / total for c in set(x)}

def read_sequence_from_file(seq_file_path: str)-> str:
    """ Reads the aminoacid composition from a ...file

    Args:
        seq_file_path: path to the file containing the sequence information
    Returns:
        aa_seq: aminoacid sequence (1 Letter code) as string
    """
    # read amino acid sequence
    with open(seq_file_path) as f:
        lines = f.readlines()
    # need to remove document start ('\ufeff') and end ('\n') to only yield the amino acid sequence
    aa_seq = lines[0].strip().replace('\ufeff', '')
    return aa_seq

def match_aminoacid_to_model_identifiers_and_frequency(aa_seq:str
                                                       ) -> Dict[str, float]:
    """ Calculates the frequency of occurance of an aminoacid in a sequence
    Also maps the aminoacids 1 letter codes to the corresponding 3 letter code
    and the BiGG identifier

    Args:
        aa_seq: aminoacid sequence parsed as string
    Returns:
        aa_biggid_freq: dict with bigg_id: a_frequency key: value pairs. The bigg_id represents the
        identifier related to the metabolite related to the specific aminoacid
    """

    # determine amino acid composition
    aa_freq = check_freq(aa_seq, len(aa_seq))

    # match to model identifiers
    aa_biggid_freq = dict()
    for aa, freq in aa_freq.items():
        threeletter = aa_lookup[aa].lower()
        if threeletter != 'gly':
            bigg_id = f'{threeletter}__L_c'
        else:
            bigg_id = f'{threeletter}_c'
        aa_biggid_freq[bigg_id] = freq
    return aa_biggid_freq

def add_aminoacid_sequence_to_intracellular_protein(model: PAModel,
                                                    seq:dict,
                                                    protein:EnzymeVariable) -> PAModel:
    """
    model: COBRA model
    seq: dict with {aminoacid_id: freq} key, value pairs
    protein: enzyme variable
    """
    for aa, freq in seq.items():
        if aa not in model.constraints.keys(): continue
        model.constraints[aa].set_linear_coefficients({
            protein.forward_variable: -freq / protein.molmass,
            protein.reverse_variable: -freq / protein.molmass
        })
    return model

def add_recombinant_intracellular_protein_to_pam(pam:PAModel,
                                   protein:Enzyme,
                                   aa_seq: str) -> PAModel:
    pam.add_enzymes(protein)
    aa_to_freq = match_aminoacid_to_model_identifiers_and_frequency(aa_seq.replace(' ', ''))
    add_aminoacid_sequence_to_intracellular_protein(pam, aa_to_freq, pam.enzyme_variables.get_by_id(protein.id))
    return pam


def add_aminoacid_sequence_to_production_rxns(model: Union[CobraModel, PAModel],
                                              seq:dict,
                                              reaction:Reaction
                                              )-> Union[CobraModel, PAModel]:
    """
    model: COBRA model
    seq: dict with {aminoacid_id: freq} key, value pairs
    reaction: reaction variable to which the aa sequence should be added
    """
    seq_metabolites = {model.metabolites.get_by_id(aa):-coeff for aa, coeff in seq.items()}
    reaction.add_metabolites(seq_metabolites)
    return model
#
def add_protein_export(model: Union[CobraModel, PAModel],
                       protein_name:str
                       )-> Reaction:
    recombinant_protein_c = Metabolite(f"{protein_name}_c",
                                       name= protein_name,
                                       formula='C20H34ClNO',
                                       charge=0,
                                       compartment='c')
    recombinant_protein_e = Metabolite(f"{protein_name}_e",
                                       name= protein_name,
                                       formula='C20H34ClNO',
                                       charge=0,
                                       compartment='e')

    protein_production_rxn = Reaction(f"{protein_name}_production")
    #also need to add a transport reaction
    protein_t_rxn = Reaction(f"{protein_name}t")
    protein_t_rxn.add_metabolites({recombinant_protein_c:-1,recombinant_protein_e:1})
    protein_production_rxn.add_metabolites({recombinant_protein_c:1})
    model.add_reactions([protein_production_rxn, protein_t_rxn])

    #add exchange reaction
    model.add_boundary(recombinant_protein_e,
                       type = 'exchange')


    return protein_production_rxn

def add_ribosome_utilization_for_exported_protein(pam:PAModel,
                                                  protein_production_rxn: Reaction,
                                                  reference_growth_rate: Optional[float] = REFERENCE_GROWTH_RATE
                                                  )-> None:
    transl_sector = pam.sectors.get_by_id('TranslationalProteinSector')

    reference_substrate_rate = _get_subtrate_uptake_rate_for_fixed_growth_rate(pam = pam,
                                                                               substrate_uptake_id=transl_sector.id_list[0],
                                                                               growth_rate=reference_growth_rate
                                                                               )
    total_metabolic_ribosome = transl_sector.tps_mu[0]*reference_substrate_rate+transl_sector.tps_0[0] #g_tps/g_p
    ribosome_per_protein = total_metabolic_ribosome/pam.total_protein_fraction

    pam.constraints[pam.TOTAL_PROTEIN_CONSTRAINT_ID].add_linear_coefficients({
        protein_production_rxn.forward_variable: ribosome_per_protein,
        protein_production_rxn.reverse_variable: -ribosome_per_protein
    })


def _get_subtrate_uptake_rate_for_fixed_growth_rate(pam:PAModel,
                                                    substrate_uptake_id: str,
                                                    growth_rate: float
                                                    )->float:
    if substrate_uptake_id == pam.BIOMASS_REACTION: return growth_rate
    #get max growth in model conditions:
    pam.optimize()
    max_mu = pam.objective.value()
    if growth_rate>max_mu:
        growth_rate=max_mu

    pam.change_reaction_bounds(substrate_uptake_id, lower_bound=-1000)
    pam.change_reaction_bounds(pam.BIOMASS_REACTION,
                               lower_bound=growth_rate, upper_bound=growth_rate)
    pam.objective.direction = 'min'
    pam.objective = {substrate_uptake_id:1}
    pam.optimize()

    return pam.reactions.get_by_id(substrate_uptake_id).flux

def add_protein_export_to_pam(pam:PAModel,
                              protein_name: str,
                              aa_seq: str,
                              reference_growth_rate: Optional[float] = REFERENCE_GROWTH_RATE
                              ) -> Reaction:
    aa_to_freq = match_aminoacid_to_model_identifiers_and_frequency(aa_seq=aa_seq)
    prot_production_rxn = add_protein_export(model=pam,
                                             protein_name=protein_name
                                             )
    add_ribosome_utilization_for_exported_protein(pam = pam,
                                                  protein_production_rxn=prot_production_rxn,
                                                  reference_growth_rate=reference_growth_rate
                                                  )
    add_aminoacid_sequence_to_production_rxns(model=pam,
                                              seq = aa_to_freq,
                                              reaction = prot_production_rxn
                                              )
    return prot_production_rxn


def add_recombinant_protein_production_and_export(aa_txt_file: str,
                                                  pam:PAModel,
                                                  protein_name: str,
                                                  export_efficiency: Union[float, int] = 1, #g_intracellular/g_exported
                                                  molecular_weight: Optional[Union[int,float]] = DEFAULT_PROTEIN_MOLMASS,#Da
                                                  reference_growth_rate: Optional[float] = REFERENCE_GROWTH_RATE
                                                  )-> Reaction:
    aa_seq = read_sequence_from_file(aa_txt_file)

    prot_production_rxn = add_protein_export_to_pam(pam = pam,
                                                    protein_name=protein_name,
                                                    aa_seq = aa_seq,
                                                    reference_growth_rate=reference_growth_rate
                                                    )#TODO add costs of ribosomes to total protein reaction

    #create intracellular enzyme species staying behind in the sell
    enzyme = Enzyme(protein_name,
                    rxn2kcat={prot_production_rxn.id:{'f':export_efficiency/3600}}, #1 protein produced per h if v_production  == 1 mmol/gcdw/h
                    molmass = molecular_weight
                    )
    add_recombinant_intracellular_protein_to_pam(pam=pam,
                                                 protein=enzyme,
                                                 aa_seq=aa_seq
                                                 )
    return prot_production_rxn