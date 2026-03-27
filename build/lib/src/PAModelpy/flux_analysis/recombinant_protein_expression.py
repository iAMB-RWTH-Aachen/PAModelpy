from src.PAModelpy import PAModel, Enzyme, EnzymeVariable

# aminoacid lookup table
aa_lookup = {'V': 'VAL', 'I': 'ILE', 'L': 'LEU', 'E': 'GLU', 'Q': 'GLN', \
             'D': 'ASP', 'N': 'ASN', 'H': 'HIS', 'W': 'TRP', 'F': 'PHE', 'Y': 'TYR', \
             'R': 'ARG', 'K': 'LYS', 'S': 'SER', 'T': 'THR', 'M': 'MET', 'A': 'ALA', \
             'G': 'GLY', 'P': 'PRO', 'C': 'CYS'}


def check_freq(x, total):
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

def match_aminoacid_to_model_identifiers_and_frequency(aa_seq:str) -> dict:
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

def add_aminoacid_sequence(model: PAModel, seq:dict, protein:EnzymeVariable):
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

def add_recombinant_protein_to_pam(pam:PAModel, protein:Enzyme, aa_seq: str) -> PAModel:
    pam.add_enzymes(protein)
    aa_to_freq = match_aminoacid_to_model_identifiers_and_frequency(aa_seq.replace(' ', ''))
    add_aminoacid_sequence(pam, aa_to_freq, pam.enzyme_variables.get_by_id(protein.id))
    return pam