import pytest
import numpy as np
from cobra.io import load_json_model
from cobra import Gene

from src.PAModelpy.PAModel import PAModel
from src.PAModelpy.configuration import Config
from src.PAModelpy.EnzymeSectors import ActiveEnzymeSector, UnusedEnzymeSector, TransEnzymeSector
from src.TAModelpy.TAModel import TAModel
from src.TAModelpy.Transcript import Transcript
from src.TAModelpy.mRNASectors import ActivemRNASector

def test_if_genes_are_added_to_pam_along_with_enzymes():
    #Act
    sut, config = build_toy_pam()

    #Assert
    # By the way the toy model is build, the number of genes should be equal to the number of proteins and they should
    # have the same naming convention (number)
    assert len(sut.genes) == len(sut.enzymes)
    for enzyme in sut.enzymes:
        assert enzyme.id[-1] == enzyme.genes[0][0].id[-1]

def test_if_tam_can_be_build():
    #Act
    sut = build_toy_tam()

def test_if_tam_adds_transcript_to_model():
    # Arrange
    sut = build_toy_tam()
    gene = Gene('test')
    sut.enzymes[0].genes.append([gene])
    enzymes = [sut.enzymes[0]]
    transcript = Transcript(id='test',
                            gene=[gene],
                            enzymes=enzymes,
                            length=10)
    mrna_sector = sut.sectors.get_by_id('ActivemRNASector')
    fmin = mrna_sector.f_min
    fmax = mrna_sector.f_max

    # Act
    sut.add_transcript(transcript, f_min=fmin, f_max=fmax)

    # Assert
    assert_transcript_is_in_model(sut, transcript)

def test_if_tam_makes_mrna_min_max_constraint_correctly():
    # Arrange
    sut = build_toy_tam()
    gene = Gene('test')
    sut.enzymes[0].genes.append(gene)
    enzymes = [sut.enzymes[0]]
    transcript = Transcript(id = 'test',
                            gene = [gene],
                            enzymes = enzymes,
                            length = 10)
    transcript.f_min = sut.f_min
    transcript.f_max = sut.f_max

    transcript._model = sut
    sut.add_cons_vars(transcript.mrna_variable)

    # Act
    sut.make_mrna_max_constraint(enzymes[0], transcript)

    # Assert
    assert_correct_min_max_mrna_constraint(sut, transcript, enzymes[0])

def test_if_enzyme_add_genes_function_add_correct_constraint_to_tamodel():
    # Arrange
    sut = build_toy_tam()
    enzyme_ut = sut.enzymes[0]
    gene_ut = 'test'
    gene_length = 10

    # Act
    enzyme_ut.add_genes(gene_ut, gene_length)
    transcript = sut.transcripts.get_by_id('mRNA_'+gene_ut)

    # Assert
    assert gene_ut in sut.genes
    assert len([gene for gene in enzyme_ut.genes if gene[0].id == gene_ut]) == 1
    assert_correct_min_max_mrna_constraint(sut, transcript, enzyme_ut)

def test_if_tamodel_adds_genes_with_and_relations_to_the_model_correctly():
    # Arrange
    sut = build_toy_tam()
    enzyme_ut = sut.enzymes[0]
    genes_ut = ['test1', 'test2']
    gene_length = [10,11]

    # Act
    enzyme_ut.add_genes(genes_ut, gene_length, relation = 'AND')
    for trans in sut.transcripts:
        for key, value in trans._constraints.items(): print(key, value)
    lumped_transcript = sut.transcripts.get_by_id('mRNA_' + '_'.join(genes_ut))

    # Assert
    assert np.sum(gene_length) == lumped_transcript.length
    assert_correct_min_max_mrna_constraint(sut, lumped_transcript, enzyme_ut)

############################################################################################################################
# HELPER FUNCTIONS
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
    return pamodel, Config

def build_active_mrna_sector(pam, Config):
    gene2transcript = {gene: {'id': f'mRNA{gene.id[-1]}', 'length': 100} for gene in pam.genes}
    return ActivemRNASector(mrnas_0 = 0, mrnas_mu = 1, id_list = ['R7'],
                            gene2transcript = gene2transcript, configuration = Config)

def build_toy_tam():
    pam, config = build_toy_pam()
    active_mrna_sector = build_active_mrna_sector(pam, config)
    tamodel = TAModel(id_or_model = pam,
                      mrna_sector = active_mrna_sector)
    tamodel.objective = 'R7'
    return tamodel

def assert_transcript_is_in_model(tamodel, transcript):
    mrna_sector = tamodel.sectors.get_by_id('ActivemRNASector')
    assert tamodel == transcript.model
    assert transcript in tamodel.transcripts
    # calculating conversion factor, *1e6 to correct for the unit conversion in enzymes (g->mg)
    assert mrna_sector.f_min* (transcript.length**2)/3 *1e6 == transcript.f_min
    assert mrna_sector.f_max* (transcript.length**2)/3 *1e6 == transcript.f_max

def assert_correct_min_max_mrna_constraint(tamodel, transcript, enzyme):

    enz_var_f = enzyme.enzyme_variable.forward_variable
    enz_var_b = enzyme.enzyme_variable.reverse_variable

    ref_min_constraint = {enz_var_f:-1.0, enz_var_b:-1.0, transcript.mrna_variable:transcript.f_min}
    ref_max_constraint = {enz_var_f:1.0, enz_var_b:1.0, transcript.mrna_variable:-transcript.f_max}

    min_constraint_coeff = tamodel.constraints[f'{transcript.id}_min'].get_linear_coefficients(
        [enz_var_f, enz_var_b, transcript.mrna_variable])
    max_constraint_coeff = tamodel.constraints[f'{transcript.id}_max'].get_linear_coefficients(
        [enz_var_f, enz_var_b, transcript.mrna_variable])

    assert all(v == min_constraint_coeff[k] for k,v in ref_min_constraint.items()) and len(ref_min_constraint) == len(min_constraint_coeff)
    assert all(v == max_constraint_coeff[k] for k,v in ref_max_constraint.items()) and len(ref_max_constraint) == len(max_constraint_coeff)