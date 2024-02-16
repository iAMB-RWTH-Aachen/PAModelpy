import warnings
from typing import List, Optional, Union, Dict, Iterable
from cobra import Object, DictList, Gene
from optlang.symbolics import Zero
import re

from src.TAModelpy.Transcript import Transcript
from src.TAModelpy.mRNASectors import ActivemRNASector
from src.PAModelpy.PAModel import PAModel
from src.PAModelpy.EnzymeSectors import ActiveEnzymeSector, TransEnzymeSector, UnusedEnzymeSector, CustomSector, Sector
from src.PAModelpy.configuration import Config
from src.PAModelpy.Enzyme import Enzyme

class TAModel(PAModel):
    TOTAL_MRNA_CONSTRAINT_ID = 'total_mrna_constraint'
    def __init__(self, id_or_model: Union[str, "Model", None] = None,
                 name: Optional[str] = None,
                 p_tot: Optional[float] = Config.P_TOT_DEFAULT,
                 sensitivity: bool = True,
                 active_sector: Optional[ActiveEnzymeSector] = None,
                 translational_sector: Optional[TransEnzymeSector] = None,
                 unused_sector: Optional[UnusedEnzymeSector] = None,
                 custom_sectors: Union[List, CustomSector] = [None],
                 mrna_sector: ActivemRNASector = None,
                 configuration=Config):
        """ Transcript allocation model

        The transcript allocation model (TAM) is a protein allocation model (PAM) extended with constraints related to
        the translation of mRNA (transcripts) into proteins using a biophysical representation of the translation
        process (based on a simple traffic model). This framework is meant to offer an easy integration of transcriptomics
        data into protein-constrained models.

        Args:
            id_or_model (string, cobra.Model, PAModelpy.PAModel): model to base the reconstruction on or model id
            name: str, optional
                    Human readable string to be model description (default None).
            p_tot : float, optional
                    Total protein concentration (condition dependent) (unit g_p/g_cdw) (default 0.285)
            senstitivity: bool
                    Boolean value wheter or not a sensitivity analysis should be performed during each simulation.
                    This sensitivity analysis will indicate to which extent individual constraints contribute to the
                    objective value.
            Enzyme sectors: EnzymeSector objects, optional
                    Information about the different enzyme sectors, being:
                - Active_enzyme
                        metabolic active proteins
                - Transl_enzyme
                        Enzymes related to translation
                - Unused_enzymes
                        Excess enzymes
                - Custom_enzymes: list
                        custom enzyme sectors
            mrna_sector: relations between mrna and growth rate allocated to the active enzyme sector and the \
                gene2transcript relations
                gene2transcript (dict(dict)): nested dict containing information about the relation between a gene
                    and a transcript. It has gene_id, transcript_info, key, value pairs were the transcript_info is
                    again a dictionary with the `id` of the trasncript and the `length` in nucleotides as keys.
            configuration: Config object, optional
                    Information about general configuration of the model including identifier conventions.
                    Default as defined in the `PAModelpy.configuration` script for the E.coli iML1515 model.

        Example:
            With a premade pamodel
            >>> pamodel = PAModel(model, name='toy model MCA with enzyme constraints', active_sector=active_enzyme,
            >>>          translational_sector = translation_enzyme,sensitivity = sensitivity,
            >>>          unused_sector = unused_enzyme, p_tot=0.6*1e-3, configuration=Config)
            >>> gene2transcript = {gene: {'id': f'mRNA{gene.id[-1]}', 'length': 100} for gene in pamodel.genes}
            >>> tamodel = TAModel(id_or_model = pamodel, mrna_sector = mrna_sector,  gene2transcript= gene2transcript)

            From scratch:
            >>>  tamodel = TAModel(model, name='toy model MCA with enzyme constraints', active_sector=active_enzyme,
            >>>          translational_sector = translation_enzyme, mrna_sector = mrna_sector, sensitivity = sensitivity,
            >>>          unused_sector = unused_enzyme, p_tot=0.6*1e-3, configuration=Config,
            >>>          gene2transcript= gene2transcript)
        """
        #TODO update examples

        self.transcripts = DictList()
        self.f_min = 1
        self.f_max = 2

        if isinstance(id_or_model, PAModel):
            self.init_tamodel_from_pamodel(id_or_model)
        else:
            super().__init__(id_or_model= id_or_model,
                     name= name,
                     p_tot= p_tot,
                     sensitivity =sensitivity,
                     active_sector = active_sector,
                     translational_sector = translational_sector,
                     unused_sector = unused_sector,
                     custom_sectors=custom_sectors,
                     configuration = configuration)

        self.add_mrna_sector(mrna_sector)

    def init_tamodel_from_pamodel(self, pamodel: PAModel) -> None:
        """ Initializes TAModel from a fully parametrized PAModel

        Args:
            pamodel (PAModel): protein allocation model to extend with transcript information
        """
        active_sector = pamodel.sectors.get_by_id('ActiveEnzymeSector')
        unused_sector = pamodel.sectors.get_by_id('UnusedEnzymeSector')
        translational_sector = pamodel.sectors.get_by_id('TranslationalProteinSector')
        custom_sectors = pamodel.sectors.query(re.compile('CustomProteinSector_*'), attribute = 'id')

        super().__init__(id_or_model= pamodel.m_model,
                         name=pamodel.name,
                         p_tot=pamodel.p_tot,
                         sensitivity=pamodel.sensitivity,
                         active_sector=active_sector,
                         translational_sector=translational_sector,
                         unused_sector=unused_sector,
                         custom_sectors=custom_sectors,
                         configuration=pamodel.configuration)

    def add_mrna_sector(self, mrna_sector):
        print('Add the following mRNA sector: ', mrna_sector.id, '\n')
        self.sectors.append(mrna_sector)
        self.add_transcript_information(mrna_sector.gene2transcript,
                                        mrna_sector.f_min, mrna_sector.f_max)
        self.add_total_mrna_constraint(mrna_sector)


    def add_transcript_information(self, gene2transcript: dict,
                                   f_min,f_max) -> None:
        """ Add gene-transcript relationships to the TAModel.

        Args:
            gene2transcript: (dict(dict)): nested dict containing information about the relation between a gene
                and a transcript. It has gene_id, transcript_info, key, value pairs were the transcript_info is
                again a dictionary with the `id` of the trasncript and the `length` in nucleotides as keys.
            f_min: conversion factor between mrna concentration and enzymes. `f_min` determines the smallest number of
                enzymes generated from a single enzyme (in mmol_enzyme/mmol_mrna)
            f_max: conversion factor between mrna concentration and enzymes. `f_max` determines the largest number of
                enzymes generated from a single enzyme (in mmol_enzyme/mmol_mrna)

        Example:
            >>> gene2transcript = {gene: {'id': f'mRNA{gene.id[-1]}', 'length': 100} for gene in pamodel.genes}
            >>> tamodel.add_transcript_information(gene2transcript=gene2transcript)
        """
        for gene, transcript_info in gene2transcript.items():
            if isinstance(gene, str): gene = Gene(gene)
            if not gene in self.genes: self.genes.append(gene)
            gene_object = self.genes.get_by_id(gene.id)
            enzymes = self.get_enzymes_by_gene_id(gene.id)

            transcript_object = Transcript(id = transcript_info['id'],
                                           gene = gene_object,
                                           enzymes = enzymes,
                                           length = transcript_info['length'])
            self.add_transcript(transcript_object, f_min, f_max)

    def add_transcript(self, transcript_object: Transcript,
                       f_min, f_max) -> None:
        """ Add Transcript object to the model

        In the model setter method of Transcript, the genes and enzymes related to the transcript will be added
        to the model if they do not exist already. Also changes the mrna-enzyme relation coefficients to those used
        in the current model

        Args:
            transcript_object: Transcript object to add
            f_min: conversion factor between mrna concentration and enzymes. `f_min` determines the smallest number of
                enzymes generated from a single enzyme (in mmol_enzyme/mmol_mrna)
            f_max: conversion factor between mrna concentration and enzymes. `f_max` determines the largest number of
                enzymes generated from a single enzyme (in mmol_enzyme/mmol_mrna)
        """
        #adjust minimal and maximal protein/transcript ratios based
        transcript_object.f_min = f_min*transcript_object.length**2/3
        transcript_object.f_max = f_max*transcript_object.length**2/3

        print('adding model to transcript ', transcript_object.id)
        transcript_object.model = self
        self.transcripts.append(transcript_object)

    def add_total_mrna_constraint(self, mrna_sector:ActivemRNASector) -> None:

        tot_mrna_constraint = self.problem.Constraint(Zero,name=self.TOTAL_MRNA_CONSTRAINT_ID,
                                                      lb=mrna_sector.intercept, ub = mrna_sector.intercept)
        self.add_cons_vars(tot_mrna_constraint)

        constraint_coefficients = {self.reactions.get_by_id(self.BIOMASS_REACTION).forward_variable:-mrna_sector.slope}
        for transcript in self.transcripts:
            constraint_coefficients = {**constraint_coefficients,
                                       **{transcript.mrna_variable:1}}

        self.constraints[self.TOTAL_MRNA_CONSTRAINT_ID].set_linear_coefficients(constraint_coefficients)


    def get_enzymes_by_gene_id(self, gene_id: str) -> DictList:
        return DictList(enzyme for enzyme in self.enzymes if self._check_if_gene_in_enzyme_genes(gene_id, enzyme))

    def _check_if_gene_in_enzyme_genes(self, gene_id: str,
                                       enzyme: Union[Enzyme, str]) -> bool:
        # check if input is in correct form and if the model exists in the model
        if isinstance(enzyme, str):
            if enzyme in self.enzymes:
                enzyme = self.enzymes.get_by_id(enzyme)
            else:
                warnings.warn('Enzyme is not in the model. Cannot find any genes')
                return False
        #find the genes in the and relationships of the enzyme [[gene1, gene2]] means gene1 and gene2
        for and_relation in enzyme.genes:
            for gene in and_relation:
                if gene.id == gene_id:
                    return True

        return False

    def make_mrna_min_max_constraint(self, enz: Enzyme, transcript:Transcript) -> Transcript:
        """
        Adding variables and constraints for the lower and upperbounds of an Enzyme to a model. The bounds are
        related to the mRNA content present in the cell.
        When solving the model, shadow prices for the lower and upperbounds will be calculated.
        This allows for the calculation of sensitivity coefficients. The constraints are formulated as follows:
        mrna_emzyme_ub : fmax * mRNA - Enzyme >= 0
        mrna_emzyme_lb : fmin * mRNA - Enzyme <= 0

        Parameters
        ----------
        m: cobra.Model or PAModelpy.PAModel
            model to which the upper and lowerbound constraints and variables should be added
        enz: PAModelpy.Enzyme
            Enzyme for which minimal and maximal concentration constraints should be generated

        Returns
        -------
        transcript : Transcript
            updated transcript object
        """
        # make separate constraints for forward and reverse enzyme variables
        fwd_var = self.enzyme_variables.get_by_id(enz.id).forward_variable
        rev_var = self.enzyme_variables.get_by_id(enz.id).reverse_variable

        #
        # if f'{self.id}_{enz.id}_max' in m.constraints.keys():
        #     m.constraints[f'{enz.id}_max'].ub = upper_bound
        # if f'{self.id}_{enz.id}_min' in m.constraints.keys():
        #     m.constraints[f'{enz.id}_min'].ub = -lower_bound
        # else:
        max_constraint = self.problem.Constraint(Zero, name=f'{transcript.id}_{enz.id}_max', lb=0, ub=1e6)
        min_constraint = self.problem.Constraint(Zero, name=f'{transcript.id}_{enz.id}_min', lb=-1e6, ub=0)
        self.add_cons_vars([max_constraint, min_constraint])

        # for key, var in self.variables.iteritems():
        #     print(key,var)
        # self.variables[transcript.id]
        # setting up the constraints
        self.constraints[f'{transcript.id}_{enz.id}_max'].set_linear_coefficients({
            fwd_var: -1,
            rev_var: -1,
            transcript.mrna_variable: transcript.f_max
        })
        self.constraints[f'{transcript.id}_{enz.id}_min'].set_linear_coefficients({
            fwd_var: -1,
            rev_var: -1,
            transcript.mrna_variable: transcript.f_min
        })

        #save the constraints in the enzyme object
        transcript._constraints = {**transcript._constraints,
                             **{f'{transcript.id}_{enz.id}_max':self.constraints[f'{transcript.id}_{enz.id}_max'],
                                f'{transcript.id}_{enz.id}_min':self.constraints[f'{transcript.id}_{enz.id}_min']
                                }}

        return transcript

