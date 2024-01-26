import pytest
from cobra.io import load_json_model

from src.PAModelpy.PAModel import PAModel
from src.PAModelpy.configuration import Config
from src.PAModelpy.EnzymeSectors import ActiveEnzymeSector, UnusedEnzymeSector, TransEnzymeSector
from src.TAModelpy.TAModel import TAModel

def test_if_genes_are_added_to_pam_along_with_enzymes():
    #Act
    sut = build_toy_pam()

    #Assert
    # By the way the toy model is build, the number of genes should be equal to the nnumber of proteins and they should
    # have the same naming convention (number)
    assert len(sut.genes) == len(sut.enzymes)
    for enzyme in sut.enzymes:
        assert enzyme.id[-1] == enzyme.genes[0][0].id[-1]

def test_if_tam_can_be_build():
    #Act
    pam = build_toy_pam()
    #TODO
    #Act
    sut = build_toy_tam(pam)


############################################################################################################################
#HELPER FUNCTIONS
############################################################################################################################
def build_active_enzyme_sector(Config):
    n = 9
    kcat_fwd = [1, 0.5, 1, 1, 0.5 ,0.45, 1.5]  # High turnover for exchange reactions
    kcat_rev = [kcat for kcat in kcat_fwd]
    protein2gene = {}
    rxn2kcat = {}
    for i in range(n-3): # all reactions have an enzyme, except excretion reactions
        rxn_id = f'R{i+1}'
        # 1e-6 to correct for the unit transformation in the model (meant to make the calculations preciser for different unit dimensions)
        #dummy molmass like in MATLAB script
        rxn2kcat = {**rxn2kcat, **{rxn_id: {f'E{i+1}':{'f': kcat_fwd[i]/(3600*1e-6), 'b': kcat_rev[i]/(3600*1e-6), 'molmass': 1e6}}}}
        protein2gene = {**protein2gene, f'E{i+1}':[[f'g{i+1}']]}
    return ActiveEnzymeSector(rxn2protein = rxn2kcat, protein2gene = protein2gene, configuration=Config)
def build_unused_protein_sector(Config):
    return UnusedEnzymeSector(id_list = ['R1'], ups_mu=[-0.01*1e-3], ups_0=[0.1*1e-3], mol_mass= [1], configuration=Config)

def build_translational_protein_sector(Config):
    return TransEnzymeSector(id_list = ['R7'], tps_mu=[0.01*1e-3], tps_0=[0.01*1e-3], mol_mass= [1], configuration=Config)

def build_toy_pam(sensitivity = True):
    Config.BIOMASS_REACTION = 'R7'
    Config.GLUCOSE_EXCHANGE_RXNID = 'R1'
    Config.CO2_EXHANGE_RXNID = 'R8'
    Config.ACETATE_EXCRETION_RXNID = 'R9'

    model = load_json_model('tests/unit_tests/toy_model.json')
    active_enzyme = build_active_enzyme_sector(Config)
    unused_enzyme = build_unused_protein_sector(Config)
    translation_enzyme = build_translational_protein_sector(Config)
    pamodel = PAModel(model, name='toy model MCA with enzyme constraints', active_sector=active_enzyme,
                      translational_sector = translation_enzyme, sensitivity = sensitivity,
                      unused_sector = unused_enzyme, p_tot=0.6*1e-3, configuration=Config)
    return pamodel

def build_toy_tam(pam:PAModel):
    gene2transcript = {gene: {'id': f'mRNA{gene.id[-1]}', 'length': 100} for gene in pam.genes}
    tamodel = TAModel(pam, gene2transcript)
    return tamodel

