#slope=0.00026415151053080677, intercept=3.5930841561397815e-05
#mrna vs growth

from src.PAModelpy.EnzymeSectors import EnzymeSector
from src.PAModelpy.configuration import Config

class ActivemRNASector(EnzymeSector):


    def __init__(self, mrnas_0:float =3.59*1e-5,
                 mrnas_mu:float = 2.64*1e-4,
                 id_list= [],
                 gene2transcript: dict = None,
                 mol_mass = 310,
                 configuration = Config()):

        """ Setting up the sector containing all the mRNAs related to proteins in the Active Enzymes Sector

        Parameters:
            mrnas_0: intercept of the relation between the mrna concentration in mmol and the reaction in id_list
              (default BIOMASS_REACTION)
            mrnas_mu: slope of the relation between the mrna concentration in mmol_mrna/unit_reaction_speed and the reaction in id_list
              (default BIOMASS_REACTION)
            id_list: id of the reaction to which the mRNA concentrations is related
            gene2transcript (dict(dict)): nested dict containing information about the relation between a gene
                    and a transcript. It has gene_id, transcript_info, key, value pairs were the transcript_info is
                    again a dictionary with the `id` of the trasncript and the `length` in nucleotides as keys.
            mol_mass: molmass of an average nucleotide in g/mol
            configuration: information about the configuration of the model to which the sector will be added

    """
        molmass = 310# g/mol default molmass of mrna nucleotide in ecoli
        if len(id_list) == 0:
            id_list = [configuration.BIOMASS_REACTION]
        super().__init__(id_list, mol_mass, configuration)
        self.mrnas_0 = mrnas_0  # amount of mrna allocated to the active enzyme sector at zero growth (mmol_mrna/g_cdw)
        self.mrnas_mu = mrnas_mu  # amount of mrna allocated to the active enzyme sector at zero growth (mmol_mrna/g_cdw/h)'
        self.id = 'ActivemRNASector'

        self.slope = None # mmol/gcdw
        self.set_slope()
        self.intercept = None #mmol/gcdw
        self.set_intercept()

        self.gene2transcript = gene2transcript
        self.elongation_rates = [12,25] #nt/s
        self.ribosome_spacing = [40,1000]#nt
        #minimal and maximal translation factors per aa in the enzyme
        self.f_min = self.elongation_rates[0]/self.ribosome_spacing[0]
        self.f_max = self.elongation_rates[1]/self.ribosome_spacing[1]

    def set_slope(self):
        self.slope = self.mrnas_mu * 1e3 # unit correction for efficient computing

    def set_intercept(self):
        self.intercept = self.mrnas_0 *1e3 # unit correction for efficient computing

