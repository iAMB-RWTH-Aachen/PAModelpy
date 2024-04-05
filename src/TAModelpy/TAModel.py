import warnings
from typing import List, Optional, Union, Dict, Iterable

import cobra
import numpy as np
from cobra import Object, DictList, Gene
from optlang.symbolics import Zero
import re

from src.TAModelpy.Transcript import Transcript, LumpedTranscript
from src.TAModelpy.mRNASectors import ActivemRNASector
from src.PAModelpy.PAModel import PAModel
from src.PAModelpy.EnzymeSectors import ActiveEnzymeSector, TransEnzymeSector, UnusedEnzymeSector, CustomSector, Sector
from src.PAModelpy.configuration import Config
from src.PAModelpy.Enzyme import Enzyme, EnzymeComplex

class TAModel(PAModel):
    TOTAL_MRNA_CONSTRAINT_ID = 'total_mrna_constraint'
    MOLMASS_MRNA_NT = 324.3 #g/mol
    def __init__(self, id_or_model: Union[str, "Model", None] = None,
                 name: Optional[str] = None,
                 p_tot: Optional[float] = Config.P_TOT_DEFAULT,
                 sensitivity: bool = True,
                 active_sector: Optional[ActiveEnzymeSector] = None,
                 translational_sector: Optional[TransEnzymeSector] = None,
                 unused_sector: Optional[UnusedEnzymeSector] = None,
                 custom_sectors: Union[List, CustomSector] = [None],
                 mrna_sector: ActivemRNASector = None,
                 total_mrna_constraint: bool = True,
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
            total_mrna_constraint (bool): boolean value to indicate if the total mRNA constraint (upperlimit to the
                    amount of mNRA available for the modeled system) should be added.
            configuration: Config object, optional
                    Information about general configuration of the model including identifier conventions.
                    Default as defined in the `PAModelpy.configuration` script for the E.coli iML1515 model.

        Example:
            With a premade pamodel
            >>> pamodel = PAModel(model, name='toy model MCA with enzyme constraints', active_sector=active_enzyme,
            >>>          translational_sector = translation_enzyme,sensitivity = sensitivity,
            >>>          unused_sector = unused_enzyme, p_tot=0.6*1e-3, configuration=Config)
            >>> gene2transcript = {gene: {'id': f'mRNA{gene.id[-1]}', 'length': 100} for gene in pamodel.genes}
            >>> active_mrna_sector = ActivemRNASector(mrnas_0 = MRNA_0,
            >>>                                  mrnas_mu = MRNA_MU,
            >>>                                  id_list = [BIOMASS_REACTION],
            >>>                                  gene2transcript = gene2transcript,
            >>>                                  configuration = config)
            >>> tamodel = TAModel(id_or_model = pamodel, mrna_sector = active_mrna_sector)

            From scratch:
            >>>  tamodel = TAModel(model, name='toy model MCA with enzyme constraints', active_sector=active_enzyme,
            >>>          translational_sector = translation_enzyme, mrna_sector = active_mrna_sector, sensitivity = sensitivity,
            >>>          unused_sector = unused_enzyme, p_tot=0.6*1e-3, configuration=Config)
        """

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

        print('Setting up the transcript allocation model, ', name)

        self.add_mrna_sector(mrna_sector, total_mrna_constraint)

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

    def add_genes(self, genes:list[list[str]], enzymes:list[Enzyme, EnzymeComplex], gene_lengths:list[int]) -> None:
        """
          Add a gene to the model and connect to associated transcripts.

          Args:
              gene (list): List of lists with the gene objects to be added. The lists represent the OR/AND relations
                    between the individual genes and the enzyme. AND: [[gene1, gene2]], OR: [[gene1], [gene2]]
              enzymes (list): A list of enzyme instances of class PAModelpy.Enzyme associated with the gene.
              gene_lengths (list): A list of lengths corresponding to each gene.

          Returns:
              None

          Notes:
              - The function creates a transcript for each gene in the `genes` list and adds it to the model.
              - If the gene is not already associated with a transcript, a new transcript is created.
              - If the gene is not already in the list of genes associated with the enzyme model,
                it is added to the list.
              - The function also establishes relationships between transcript and enzymes
                associated with the genes.

          """
        mrna_sector = self.sectors.get_by_id('ActivemRNASector')
        transcript_objects = []
        length_index = 0
        for gene_list in genes:
            for gene in gene_list:
                # Add the gene to the model if not already present
                length = gene_lengths[length_index]
                length_index+=1
                if gene not in self.genes:
                    self.genes.append(gene)

                # Create the transcript object
                gene_object = self.genes.get_by_id(gene.id)
                transcript_objects += [Transcript(id='mRNA_' + gene.id,
                                               gene=[gene_object],
                                               enzymes=enzymes,
                                               length=length)]

        # Add the transcripts to the model
        for transcript in transcript_objects:
            self.add_transcript(transcript, mrna_sector.f_min, mrna_sector.f_max)

        # Establish relationships between the transcript and enzymes
        for enzyme in enzymes:
            for gene_relation in genes:
                if gene_relation not in enzyme.genes: enzyme.genes.append(genes)
            self.add_transcript_enzyme_relations(enzyme)

    def add_mrna_sector(self, mrna_sector: ActivemRNASector, total_mrna_constraints: bool = True):
        print('Add the following mRNA sector: ', mrna_sector.id, '\n')
        self.sectors.append(mrna_sector)
        if total_mrna_constraints:
            self.add_total_mrna_constraint(mrna_sector)
        self.add_transcript_information(mrna_sector.gene2transcript,
                                        mrna_sector.f_min, mrna_sector.f_max)



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
        gene_transcript_info = gene2transcript.copy()
        for gene, transcript_info in gene_transcript_info.items():
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
                       f_min:float = None, f_max:float = None) -> None:
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
        mrna_sector = self.sectors.get_by_id('ActivemRNASector')
        if f_min is None: f_min = mrna_sector.f_min
        if f_max is None: f_max = mrna_sector.f_max

        #add to mRNA sector:
        mrna_sector.add_transcript(transcript_object)

        #adjust minimal and maximal protein/transcript ratios based
        transcript_object.f_min = f_min*pow(transcript_object.length,2)/3*1e6
        transcript_object.f_max = f_max*pow(transcript_object.length,2)/3*1e6
        # print('adding model to transcript ', transcript_object.id)
        transcript_object.model = self
        if self.TOTAL_MRNA_CONSTRAINT_ID in self.constraints.keys():
            self.add_transcripts_to_total_mrna_constraint(transcripts = [transcript_object])

    def add_transcript_enzyme_relations(self, enzyme: Union[Enzyme, EnzymeComplex]) -> None:
        """Adds constraints for mRNA abundance and enzyme concentration relations.

            This function checks the gene-protein-reaction (GPR) relations to determine if the type of relation should
            be included in the transcript/enzyme relationship ('and' or 'or').
            Then, it adds constraints based on these relations.

            Args:
                enzyme (Union[Enzyme, EnzymeComplex]): The enzyme or enzyme complex for which constraints are being added.

            Returns:
                None
            """
        # check the gpr relations: is an additional gene required or is there 'competition' with another gene?
        transcript_associated_with_enzyme = self.get_transcripts_associated_with_enzyme(enzyme)
        for relation, transcripts in transcript_associated_with_enzyme.items():
            self.make_mrna_min_max_constraint(enzyme, transcripts, relation)


    def add_total_mrna_constraint(self, mrna_sector:ActivemRNASector) -> None:
        """Adds a constraint to limit the total mRNA abundance.

           This function adds a constraint to limit the total mRNA abundance based on the intercept
           adn slope provided by the the ActivemRNA sector.
           It calculates the coefficients for each transcript and sets the linear coefficients for the constraint.

           Args:
               mrna_sector (ActivemRNASector): The active mRNA sector.

           Returns:
               None
           """
        tot_mrna_constraint = self.problem.Constraint(Zero,name=self.TOTAL_MRNA_CONSTRAINT_ID,
                                                      lb=-1e6, ub = mrna_sector.intercept)
        self.add_cons_vars(tot_mrna_constraint)

        constraint_coefficients = {self.reactions.get_by_id(self.BIOMASS_REACTION).forward_variable:-mrna_sector.slope}
        self.constraints[self.TOTAL_MRNA_CONSTRAINT_ID].set_linear_coefficients(constraint_coefficients)

    def add_transcripts_to_total_mrna_constraint(self, transcripts:list):
        constraint_coefficients = {}
        for transcript in transcripts:
            coeff = transcript.length*self.MOLMASS_MRNA_NT
            # print(transcript.id, coeff)
            constraint_coefficients = {**constraint_coefficients,
                                           **{transcript.mrna_variable:coeff}}
        self.constraints[self.TOTAL_MRNA_CONSTRAINT_ID].set_linear_coefficients(constraint_coefficients)
        self.solver.update()

    def correct_unused_enzymes_for_measured_mrna(self, enzyme: Enzyme):
        enzyme_variable = enzyme.enzyme_variable
        coef = self.constraints[self.TOTAL_PROTEIN_CONSTRAINT_ID].get_linear_coefficients([enzyme_variable.forward_variable])

        self.constraints[self.TOTAL_PROTEIN_CONSTRAINT_ID].set_linear_coefficients({
            enzyme_variable.forward_variable: 0.5*enzyme.molmass * 1e-6,
            enzyme_variable.reverse_variable: 0.5*enzyme.molmass * 1e-6,
        })
        coef2 = self.constraints[self.TOTAL_PROTEIN_CONSTRAINT_ID].get_linear_coefficients([enzyme_variable.forward_variable])
        print(coef, coef2)

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

    def make_mrna_min_max_constraint(self, enz: Enzyme,
                                     transcripts:Union[Transcript, list],
                                     relation = 'or') -> Union[Transcript, list]:
        """
        Adding variables and constraints for the lower and upperbounds of an Enzyme to a model. The bounds are
        related to the mRNA content present in the cell.
        When solving the model, shadow prices for the lower and upperbounds will be calculated.
        This allows for the calculation of sensitivity coefficients. The constraints are formulated as follows:
        mrna_emzyme_ub : fmax * mRNA - Enzyme >= 0
        mrna_emzyme_lb : fmin * mRNA - Enzyme <= 0

        Args:
        ----------
        m: cobra.Model or PAModelpy.PAModel
            model to which the upper and lowerbound constraints and variables should be added
        enz: PAModelpy.Enzyme
            Enzyme for which minimal and maximal concentration constraints should be generated
        transcript: Transcript
            Transcripts or list of Transcripts to relate to the enzyme. If the list of transcripts is longer than 1,
            the relation between the transcripts and enzyme is considered to be an 'or' relation.
        relation: str
            Determines the type of relation with which the transcripts are related to the enzyme. In an `or` relation,
            all transcripts are added to the min/max constraints, in an `and` relation only the lowest available
            transcript is used to constrain the enzyme using a dummy variable:
            and_mrna <= mrna1
            and_mrna <= mrna2
            and_mrna - E <= 0

        Returns:
        -------
        transcript : Transcript
            updated transcript object
        """
        if not isinstance(transcripts, list): transcripts = [transcripts]

        # make separate constraints for forward and reverse enzyme variables
        fwd_var = self.enzyme_variables.get_by_id(enz.id).forward_variable
        rev_var = self.enzyme_variables.get_by_id(enz.id).reverse_variable

        #
        # if f'{self.id}_{enz.id}_max' in m.constraints.keys():
        #     m.constraints[f'{enz.id}_max'].ub = upper_bound
        # if f'{self.id}_{enz.id}_min' in m.constraints.keys():
        #     m.constraints[f'{enz.id}_min'].ub = -lower_bound
        # else:
        transcripts_ids = "_".join([transcript.id for transcript in transcripts])
        # max_constraint = self.problem.Constraint(Zero, name=f'{transcripts_ids}_{enz.id}_max', lb=0, ub=1e6) #TODO check relationships
        # min_constraint = self.problem.Constraint(Zero, name=f'{transcripts_ids}_{enz.id}_min', lb=0, ub=1e6)
        if f'{transcripts_ids}_max' not in self.constraints.keys():
            max_constraint = self.problem.Constraint(Zero, name=f'{transcripts_ids}_max', lb=-1e6,
                                                 ub=0)  # TODO check relationships
            min_constraint = self.problem.Constraint(Zero, name=f'{transcripts_ids}_min', lb=-1e6, ub=0)
            self.add_cons_vars([max_constraint, min_constraint])

        # for key, var in self.variables.iteritems():
        #     print(key,var)
        # self.variables[transcript.id]
        # # setting up the constraints
        self.constraints[f'{transcripts_ids}_max'].set_linear_coefficients({
            fwd_var: 1,
            rev_var: 1
        })
        self.constraints[f'{transcripts_ids}_min'].set_linear_coefficients({
            fwd_var: -1,
            rev_var: -1
        })

        for transcript in transcripts:
            self.set_mrna_enzyme_minmax_constraint_or_relation(enz, transcript, transcripts_ids)
            if transcript not in enz.transcripts: enz.transcripts.append(transcript)

            #save the constraints in the transcript object
            transcript._constraints = {**transcript._constraints,
                                        **{f'{transcripts_ids}_max':self.constraints[f'{transcripts_ids}_max'],
                                    f'{transcripts_ids}_min':self.constraints[f'{transcripts_ids}_min']
                                    }}
            # transcript.mrna_variable.ub =0
        return transcripts

    def _create_mrna_minmax_constraints(self, enzyme: Enzyme, transcript_ids: str,
                                        min_variable, max_variable,
                                        f_min: float, f_max: float) -> None:
        self.constraints[f'{transcript_ids}_max'].set_linear_coefficients({
            max_variable: -f_max
        })
        self.constraints[f'{transcript_ids}_min'].set_linear_coefficients({
            min_variable: f_min
        })

    def set_mrna_enzyme_minmax_constraint_and_relation(self, enzyme: Enzyme,
                                                       transcripts: list, transcript_ids: list) -> None:
        """Sets constraints and relations for mRNA abundance and enzyme concentration.

            The lowest mRNA concentration of the available transcripts is used. This is achieved by adding an
            auxiliary variable for the min and max enzyme concentration in the following format:
            ```
            mrna_aux_min <= mrna1 * fmin_1
            mrna_aux_min <= mrna2 * fmin_2

            mrna_aux_max <= mrna1 * fmax_1
            mrna_aux_max <= mrna2 * fmax_2

            mrna_aux_min - E >= 0
            mrna_aux_max - E <= 0
            ```
           Args:
               enzyme (Enzyme): The enzyme for which constraints and relations are being set.
               transcripts (list): List of Transcript objects.
               transcript_ids (list): List of transcript IDs.

           Returns:
               None
        """
        transcript_aux_min_variable = self.problem.Variable(Zero, name=f'{transcript_ids}_{enzyme.id}_min_aux', lb=0, ub=1e6)
        transcript_aux_max_variable = self.problem.Variable(Zero, name=f'{transcript_ids}_{enzyme.id}_max_aux', lb=0, ub=1e6)

        self.add_cons_vars([transcript_aux_min_variable, transcript_aux_max_variable])

        #Add auxiliary constraints for all transcripts to ensure the minimal mrna abundance in the enzyme complex
        # is used to constrain the enzyme concentration
        for transcript in transcripts:
            aux_min_constraint = self.problem.Constraint(Zero, name=f'{transcript.id}_{enzyme.id}_min_aux', lb=0, ub=1e6)
            aux_max_constraint = self.problem.Constraint(Zero, name=f'{transcript.id}_{enzyme.id}_max_aux', lb=0, ub=1e6)
            self.add_cons_vars([aux_min_constraint, aux_max_constraint])

            # this auxiliary variable is already converted to g_enzyme/g_cdw
            self.constraints[f'{transcript.id}_{enzyme.id}_min_aux'].set_linear_coefficients({
                transcript_aux_min_variable: -1,
                transcripts.mrna_variable: transcript.f_min
            })
            self.constraints[f'{transcript.id}_{enzyme.id}_max_aux'].set_linear_coefficients({
                transcript_aux_max_variable: -1,
                transcripts.mrna_variable: transcript.f_max
            })

            #save the constraints in the transcript object
            transcript._constraints = {**transcript._constraints,
             **{f'{transcript.id}_{enzyme.id}_min_aux': self.constraints[f'{transcript.id}_{enzyme.id}_min_aux'],
                f'{transcript.id}_{enzyme.id}_max_aux': self.constraints[f'{transcript.id}_{enzyme.id}_max_aux']
                }}

        self._create_mrna_minmax_constraints(enzyme = enzyme,
                                             transcript_ids = transcript_ids,
                                             min_variable= transcript_aux_min_variable,
                                             max_variable = transcript_aux_max_variable,
                                             f_min = 1,
                                             f_max = 1
                                             )

    def set_mrna_enzyme_minmax_constraint_or_relation(self, enzyme: Enzyme,
                                                           transcript: Transcript, transcript_id: str):
        """Sets constraints for mRNA abundance and enzyme concentration based on 'or' relationship.
            This is in the following format:
            ```
            mrna_1 * fmin_1 + mrna_2 * fmin_2 - E >= 0
            mrna_1 * fmax_1 + mrna_2 * fmax_2 - E <= 0
            ```

           Args:
               enzyme (Enzyme): The enzyme for which constraints are being set.
               transcript (Transcript): The transcript object.
               transcript_ids (str): the id of the transcript

           Returns:
               None
           """
        self._create_mrna_minmax_constraints(enzyme = enzyme,
                                             transcript_ids = transcript_id,
                                             min_variable= transcript.mrna_variable,
                                             max_variable=transcript.mrna_variable,
                                             f_min = transcript.f_min,
                                             f_max = transcript.f_max
                                             )


    def get_genes_associated_with_enzyme(self, enzyme: Enzyme) -> dict:
        """ Retrieve the model associated genes including their relation from the model

        Goes through all the genes associated with an enzyme and checks if there is an 'and' or 'or' relation
        between the genes. It only saves those genes, which are already associated with the model

        Args:
            enzyme: the enzyme for which the gene list should be evaluated

        Returns:
            gene_relations_dict: a dictionary with all the cobra.Gene objects associated with both the model as
            the enzyme and their relationship

        Example:
            Example of returned object:
            for an enzyme with the following list of genes: `[[gene1, gene2],[gene3, gene4],[gene5]]`
            In this enzyme, there are 2 enzyme complexes with the ids 0 and 1.
            `gene_relation_dict = {'and':[[gene1, gene2], [gene3, gene4]],
                                    'or':[gene5]}`
        """
        gene_relations_dict = {'and':[], 'or':[]}
        # enzyme_complex_id = 0
        for gene_list in enzyme.genes:
            if len(gene_list)>1: # and relationship
                gene_relations_dict['and'] += [[gene for gene in gene_list if gene in self.genes]]
            elif gene_list[0] in self.genes:
                    gene_relations_dict['or'] += [gene_list[0]]
        return gene_relations_dict

    def get_transcripts_associated_with_enzyme(self, enzyme: Enzyme, output: str = 'dict') -> dict:
        """ Retrieve the model associated transcripts including their relation from the model

        Goes through all the genes associated with an enzyme and checks if there is an 'and' or 'or' relation
        between the transcripts. It only saves those transcripts, which are already associated with the model

        Args:
            enzyme: the enzyme for which the gene list should be evaluated

        Returns:
            mrna_relations_dict: a dictionary with all the TAModelpy.Transcript objects associated with both the model as
            the enzyme and their relationship

        Example:
            Example of returned object:
            for an enzyme with the following list of genes: `[[gene1, gene2],[gene3, gene4],[gene5]]`
            In this enzyme, there are 2 enzyme complexes with the ids 0 and 1.
            `mrna_relation_dict = {'and':[[mRNA_gene1, mRNA_gene2], [mRNA_gene3, mRNA_gene4]],
                                    'or':[mRNA_gene5]}`
        """
        if output == 'dict':
            return self.get_transcript_enzyme_associations(enzyme)
        else:
            transcripts = []
            for gene_list in enzyme.genes:
                transcripts += [self.transcripts.get_by_id('mRNA_' + gene.id) for gene in gene_list if
                                                'mRNA_' + gene.id in self.transcripts]
            return transcripts

    def get_transcripts_associated_with_reaction(self, reaction: Union[cobra.Reaction, str]):
        if not isinstance(reaction, str): reaction = reaction.id
        enzymes = self.get_enzymes_with_reaction_id(reaction)
        transcripts = []
        for enzyme in enzymes:
            transcripts += self.get_transcripts_associated_with_enzyme(enzyme, output='list')
        return transcripts

    def get_transcript_enzyme_associations(self, enzyme: Enzyme):
        mrna_relations_dict = {'and': [], 'or': []}
        for gene_list in enzyme.genes:
            gene_list = [gene.id if isinstance(gene, Gene) else gene for gene in gene_list]
            if len(gene_list) > 1:  # and
                lumped_transcript = self.make_lumped_transcript_for_and_relation(enzyme = enzyme, gene_list=gene_list)
                mrna_relations_dict['and'] += [lumped_transcript]
                # mrna_relations_dict['and'] += [[self.transcripts.get_by_id('mRNA_' + gene) for gene in gene_list if
                #                                 'mRNA_' + gene in self.transcripts]]
            elif 'mRNA_' + gene_list[0] in self.transcripts:
                mrna_relations_dict['or'] += [self.transcripts.get_by_id('mRNA_' + gene_list[0])]

        result_dict = mrna_relations_dict.copy()
        # delete empty dict entries
        for key, value in mrna_relations_dict.items():
            if len(value) == 0:
                del result_dict[key]
        return result_dict

    def make_lumped_transcript_for_and_relation(self, enzyme:Enzyme, gene_list: list):
        transcript_list = [self.transcripts.get_by_id('mRNA_' + gene) for gene in gene_list if
                           'mRNA_' + gene in self.transcripts]
        new_transcript_length = np.nansum([trans.length for trans in transcript_list]) # need to translate all mRNAs to get enzyme complex

        new_transcript_id = 'mRNA_'+ '_'.join(gene_list)
        if new_transcript_id in self.transcripts:
            new_transcript = self.transcripts.get_by_id(new_transcript_id)
            new_transcript.length = new_transcript_length
            for trans in transcript_list:
                if new_transcript.id not in trans._lumped_transcripts:
                    trans._lumped_transcripts.append(new_transcript)
                    new_transcript._transcripts.append(trans)
            if enzyme not in new_transcript.enzymes: new_transcript.enzymes.append(enzyme)
        else:
            new_transcript = LumpedTranscript(id = new_transcript_id,
                                    gene = [self.genes.get_by_id(gene_id) for gene_id in gene_list],
                                    enzymes = [enzyme],
                                    length = new_transcript_length)
            self.add_transcript(new_transcript)
            if not new_transcript in enzyme.transcripts: enzyme.transcripts.append(new_transcript)
            for trans in transcript_list:
                trans._lumped_transcripts.append(new_transcript)
                new_transcript._transcripts.append(trans)
        # self.remove_transcripts(transcript_list)
        return new_transcript

    def remove_transcripts(self, transcripts: list):
        for transcript in transcripts:
            to_delete = [cons for cons in transcript._constraints.values() if cons in self.constraints.values()] + [transcript.mrna_variable]
            self.remove_cons_vars(to_delete)
            self.transcripts.remove(transcript)
            self.solver.update()

    def _gene_is_in_enzyme_complex(self, gene: Union[str, Gene], enzyme: Enzyme) -> bool:
        if not isinstance(gene, str): gene = gene.id
        if isinstance(enzyme, EnzymeComplex) and any([gene in genelist for genelist in enzyme.genes]):
            return True
        else:
            return False

    def calculate_csc(self, obj_value, mu, mu_ub, mu_lb, mu_ec_f, mu_ec_b):
        """
        Calculates the capacity sensitivity coefficient for all inequality constraints in the model.
        The sum of all capacity sensitivity coefficients should equal 1 for growth maximization.

        Definition: constraint_UB*shadowprice/obj_value.

        Inherits from the PAModel calculate_csc function and adds the calculation of the csc for the total rna
        constraint.

        :param obj_value: Float
        :param mu: DataFrame
            Shadowprices for all constraints
        :param mu_ub: DataFrame
            Shadowprices for the reaction UB constraints
        :param mu_lb: DataFrame
            Shadowprices for the reaction LB constraints
        :param mu_ec_f: DataFrame
            Shadowprices for the constraint related to an enzymatic catalysis of the forward reaction
        :param mu_ec_b: DataFrame
            Shadowprices for the constraint related to an enzymatic catalysis of the backward reaction

        Results will be saved in the self.capacity_sensitivity_coefficients attribute as a dataframe
        """
        super().calculate_csc(obj_value, mu, mu_ub, mu_lb, mu_ec_f, mu_ec_b)
        if self.TOTAL_MRNA_CONSTRAINT_ID in self.constraints.keys():
            constraint = 'mRNA'
            rxn_id = self.TOTAL_MRNA_CONSTRAINT_ID
            enzyme_id = self.TOTAL_MRNA_CONSTRAINT_ID
            ca_coefficient = self.constraints[enzyme_id].ub * \
                             mu[mu['index'] == self.TOTAL_MRNA_CONSTRAINT_ID]['shadow_prices'].iloc[0] / obj_value

            new_row = [rxn_id, enzyme_id, constraint, ca_coefficient]
            #  add new_row to dataframe
            self.capacity_sensitivity_coefficients.loc[len(self.capacity_sensitivity_coefficients)] = new_row

        #TODO sensitivity analysis for transcripts