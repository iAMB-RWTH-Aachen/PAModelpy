from warnings import warn
from copy import copy, deepcopy

from cobra import Object, Gene, Model
from typing import Union, Literal, Dict, Optional, List
from collections import defaultdict

from .Enzyme import Enzyme, EnzymeComplex
from .configuration import Config


class Sector(Object):
    TOTAL_PROTEIN_CONSTRAINT_ID = Config.TOTAL_PROTEIN_CONSTRAINT_ID

    def get_tpc_metabolite(self, model):
        tpc = self.TOTAL_PROTEIN_CONSTRAINT_ID
        try:
            tpc_metabolite = model.constraints[tpc]
            totprot = True
        except:
            tpc_metabolite = None
            totprot = False
        return totprot, tpc_metabolite

    def __copy__(self) -> "Sector":
        """Copy the Sector.

        Returns:
            Sector: A new Sector that is a copy of the original Sector.
        """

        cop = copy(super(Sector, self))
        return cop

    def __deepcopy__(self, memo: dict) -> "Sector":
        """Copy the Sector with memo.

        Args:
            memo (dict): Automatically passed parameter.

        Returns:
            Sector: A new Sector that is a copy of the original Sector with memo.
        """

        cop = deepcopy(super(Sector, self), memo)
        return cop


class EnzymeSector(Sector):
    DEFAULT_MOL_MASS = 3.947778784340140e04  # mean enzymes mass E.coli [g/mol]

    def __init__(self, id_list, mol_mass, configuration=Config):
        self.id_list = id_list  # list with reaction ids which are associated with the specific sector
        self.mol_mass = (
            mol_mass  # dict with molar masses for each enzyme in the id_list (g/mol)
        )
        self.variables = []
        self.constraints = []
        self.model = None
        self.id = ""
        self.slope = 0
        self.intercept = 0

        self.TOTAL_PROTEIN_CONSTRAINT_ID = configuration.TOTAL_PROTEIN_CONSTRAINT_ID


class ActiveEnzymeSector(Sector):
    DEFAULT_MOL_MASS = 3.947778784340140e04  # mean enzymes mass E.coli [g/mol]
    #class with all the information on enzymes related to metabolic reactions
    def __init__(self, rxn2protein: dict, protein2gene:dict = {},configuration:Config =Config):
        """
        _summary_

        Parameters:
            rxn2protein (dict): A dictionary with reaction ID, enzymes_dict key, value pairs for each reaction in the active_enzyme sector.
                The enzymes_dict contains the enzyme identifiers of the enzymes related to the specific reaction as keys, and a dict
                with information about the enzyme as values. The information included in this dict includes the turnover number for
                the forward and backward reaction (1/s), molar mass of the enzyme (mol/g), gene identifiers related to the enzyme,
                and with which other enzymes it forms a complex.
            protein2gene: dict
            enzyme_id, gene_list key, value pairs for each enzyme.The gene_list value is a list of lists which indicates
            'and' or 'or' relationships between the genes which code for the enzyme(complex).

        Example:
            ```
            For the Parameter rxn2protein a dictionary may look like this:
            {
                'R1':
                    {E1:
                                {'f': forward kcat, 'b': backward kcat, 'molmass': molar mass, 'genes': [G1, G2],
                                'protein_reaction_association': [[E1, E2]]
                                },
                            E2:
                                {'f': forward kcat, 'b': backward kcat, 'molmass': molar mass, 'genes': [G3, G4],
                                 'protein_reaction_association': [[E1, E2]]
                    }
            }

            For the Parameter protein2gene a dictionary may look like this:
            {E1:
                    [[gene1], [gene2, gene3]],
            E2:
                    [[gene4]]
            }
            where the gene-protein-reaction associations are the following:
            E1: gene1 or (gene2 and gene3)
            E2: gene4
            ```
        Args:
        rxn2protein : nested dict
            reaction id, enzymes_dict key, value pairs for each reaction in the active_enzyme sector.
             The enzymes_dict contains the enzyme identifiers of the enzymes which are related to the specific reaction
             as keys and a dict with information about the enzyme as values. The information included in this dict is:
              turnover number for forward and backward reaction [1/s], molar mass of the enzyme [mol/g], gene identifiers
              related to the enzyme, with which other enzymes it forms a complex.

        protein2gene: dict
            enzyme_id, gene_list key, value pairs for each enzyme.The gene_list value is a list of lists which indicates
            'and' or 'or' relationships between the genes which code for the enzyme(complex).

        configuration: Config object, optional
                    Information about general configuration of the model including identifier conventions.
                    Default as defined in the `PAModelpy.configuration` script for the E.coli iML1515 model.
        """
        self.TOTAL_PROTEIN_CONSTRAINT_ID = configuration.TOTAL_PROTEIN_CONSTRAINT_ID

        self.rxn2protein = rxn2protein
        self.protein2gene = protein2gene
        # ID of the sector
        self.id = "ActiveEnzymeSector"
        self.model = None

        # for easy removal
        self.variables = []
        self.constraints = []

    def add(self, model):
        self.model = model
        # make and enzymatic constraint for an existing stoichiometric reaction in the form of: e_i=v_i*kcat_i
        # c.f. the GECKO formulation of Sanchez et al. 2017
        print("Add active protein sector\n")
        # check if the constraint exists
        totprot, tpc_metabolite = self.get_tpc_metabolite(model)
        if not totprot:
            warn(
                "Total protein constraint is not in the model. Skip the Active Enzyme Sector"
            )
            return model

        rxn2protein = self.rxn2protein.copy()
        self.add_rxn2protein(rxn2protein)

        return model

    def add_rxn2protein(self, rxn2protein: Dict[
        str,Dict[
            str,
            Literal['f', 'b', 'molmass',
            'genes', 'protein_reaction_association']
        ]],
                        protein2gene:Optional[Dict[str,list]]= None,
                        verbose: bool = False) -> None:
        """
        Adds enzyme  based on the gene-protein-reaction association as defined in the rxn2protein dictionary
        to an existing PAM. Several checks are performed before adding it to the model associated with this sector.

        Args:
            rxn2protein (dict): dictionary containing gene-protein-reaction association.
            protein2gene (dict): dictionary containing mapping between peptide ids and genes
            verbose (bool): if true, ensures all added enzymes are printed when adding

        Example:
            ```
            For the Parameter rxn2protein a dictionary may look like this
            (E1 and E2 are peptides in a single enzyme complex):
            {
                'R1':
                    {E1_E2:
                                {'f': forward kcat,
                                'b': backward kcat,
                                'molmass': molar mass,
                                'genes': [G1, G2],
                                'protein_reaction_association': [[E1, E2]]
                                }
            }

        The associated protein2gene dictionary would look like this:
        ```
        {E1: [[G1, G2]]}
        ```

        """
        if protein2gene is not None:
            self._add_protein2gene_information(protein2gene)
        self._add_rxn2protein_information(rxn2protein)
        model = self.model

        for rxn_id, enzymes in rxn2protein.copy().items():
            # extract reaction from model
            if rxn_id not in model.reactions:
                #     reaction = model.reactions.get_by_id(rxn_id)
                # else:
                try:
                    reaction = model.reactions.query(rxn_id)
                    for rxn in reaction:
                        rxn2protein[rxn.id] = enzymes
                        if rxn_id in rxn2protein.keys():
                            del rxn2protein[rxn_id]
                except:
                    warn(
                        f"Reaction {rxn_id} is not in the model, this reaction will be skipped"
                    )
                    if rxn_id in rxn2protein.keys():
                        del rxn2protein[rxn_id]

        for rxn_id, enzymes in rxn2protein.items():
            if verbose: print(f'\nAdding an association between reaction {rxn_id} and the following enzymes {list(enzymes.keys())}')

            reaction = model.reactions.get_by_id(rxn_id)
            # skip the reaction if a problem is encountered in check_kcat_values()
            consistent = self.check_kcat_values(model, reaction, enzymes)
            if not consistent:
                continue

            for enzyme_id, enzyme_dict in enzymes.items():
                kcat = dict(
                    (k, v) for k, v in enzyme_dict.items() if k == "f" or k == "b"
                )
                if 'protein_reaction_association' in enzyme_dict.keys():
                    protein_reaction = enzyme_dict['protein_reaction_association']
                else:
                    protein_reaction = enzyme_id

                # get molar mass of enzyme or replace with default value
                if "molmass" in enzyme_dict.keys():
                    molmass = enzyme_dict["molmass"]
                else:
                    molmass = self.DEFAULT_MOL_MASS

                # check if there already exists an Enzyme object for the enzyme associated to this reaction,
                # otherwise create one and store all information
                if enzyme_id in model.enzyme_variables and not self._enzyme_is_enzyme_complex(protein_reaction, enzyme_id):
                    if verbose: print(f'\tEnzyme {enzyme_id} is already associated with {rxn_id}')
                    enzyme = model.enzymes.get_by_id(enzyme_id)
                    self._add_reaction_to_enzyme(model, enzyme, rxn_id, kcat)
                    self.rxn2protein[rxn_id].update({
                        enzyme_id: {
                            **kcat,
                            'genes': enzyme.genes,
                            'protein_reaction_association': protein_reaction
                        }
                    })

                else:
                    if self.protein2gene != {}:
                        gene_list = self._get_model_genes_from_enzyme(enzyme_id, model)
                    else:
                        gene_list = []

                    enzyme = Enzyme(
                        id=enzyme_id,
                        rxn2kcat={rxn_id: kcat},
                        rxn2pr={rxn_id: protein_reaction},
                        molmass=molmass,
                        genes=gene_list
                    )
                    if self._enzyme_is_enzyme_complex(protein_reaction, enzyme_id):
                        for pr in protein_reaction:
                            if len(pr) > 1:
                                #need to sort so the order of the id does not make the enzyme complex unique
                                #e.g.: B1_B2_B3 should not be added if B3_B1_B2 is already in the model
                                enzyme_complex_id = '_'.join(sorted(pr))

                                if enzyme_complex_id not in model.enzyme_variables:
                                    enzyme = EnzymeComplex(
                                        id=enzyme_complex_id,
                                        rxn2kcat={rxn_id: kcat},
                                        molmass=0,
                                        genes=gene_list,
                                        enzymes=[enzyme]
                                    )
                                    model.add_enzymes([enzyme])
                                    if verbose: print(f'\tEnzyme complex {enzyme_complex_id} is added to the model')
                                    #add relation to rxn2protein dictionary

                                    self.rxn2protein[rxn_id] = {**self.rxn2protein[rxn_id],
                                                                **{enzyme_complex_id: {
                                                                    **kcat,
                                                                    'genes': [g.id for g in gene_list[0]],
                                                                    'protein_reaction_association': pr}}}
                                    self.constraints += [enzyme]
                                    self.variables.append(enzyme.enzyme_variable)

                                else:
                                    enz_complex = model.enzymes.get_by_id(enzyme_complex_id)
                                    self._add_reaction_to_enzyme(model, enz_complex, rxn_id, kcat)
                                    enz_complex.reactions.append(rxn_id)
                                    enz_complex.add_enzymes([enzyme])

                    else:
                        # in the add enzymes function the enzyme with corresponding
                        # catalytic events will be added to the model
                        # and connected to the reaction and total protein constraint
                        model.add_enzymes([enzyme])
                        if verbose: print(f'\tEnzyme {enzyme.id} is added to the model')

                        self.constraints += [enzyme]
                        self.variables.append(enzyme.enzyme_variable)

    def check_kcat_values(self, model, reaction, enzyme_dict):
        """
        Function to check if the kcat values provided for an enzyme are consistent with the direction of the reaction
        Args:
            model (cobra.Model or PAModel): Model to which the kcat values should be added.
            reaction (cobra.Reaction): Reaction that is catalyzed with the enzyme related to the kcats.
            enzyme_dict (dict): A dictionary with the turnover values for the forward and/or backward reaction for different enzymes [/s].

        Example:
            Example dictionary for the `enzyme_dict` parameter
            ```
            {'E1': {'f': forward kcat, 'b': backward kcat}}
            ```
        """
        if reaction.id not in model.reactions:
            warn(
                "Reaction " + reaction.id + " not in the model. Skip enzyme constraint"
            )
            return False
        # check if there is an EC number associated to each reaction, otherwise, use the reaction id
        if reaction.id not in self.rxn2protein.keys():
            self.rxn2protein[reaction.id] = reaction.id

        # check consistency between provided kcat values and reaction direction
        directions = []
        kcats = []
        for enzyme_info in enzyme_dict.values():
            # get all directions from the kcat dict
            for key, value in enzyme_info.items():
                if key == 'f' or key =='b':
                    directions += [key]
                    kcats += [value]

        rxn_dir = "consistent"
        if reaction.lower_bound >= 0 and reaction.upper_bound > 0:
            # reaction is irreversible in the forward direction

            if not "f" in directions:
                rxn_dir = "inconsistent"

        elif reaction.lower_bound < 0 and reaction.upper_bound <= 0:
            # reaction is irreversible in the backward direction
            if not "b" in directions:
                rxn_dir = "inconsistent"
        else:
            # reaction is reversible
            if not "f" in directions and "b" in directions:
                rxn_dir = "inconsistent"

        if rxn_dir == "inconsistent":
            warn(
                reaction.id
                + ": reaction directionality does not match provided kcat values. Skip reaction"
            )
            return False

        # check kcat values input
        for k in kcats:
            if k < 0:
                warn(
                    [
                        'Turnover number for reaction "',
                        reaction.id,
                        '" is invalid. Skip for active enzyme sector',
                    ]
                )
                return False

        return rxn_dir == "consistent"

    def _add_protein2gene_information(self, protein2gene_to_add:Dict[str,List])->None:
        merged_dict =self.protein2gene.copy()
        for key, values in protein2gene_to_add.items():
            if key in merged_dict:
                for sublist in values:
                    if sublist not in merged_dict[key]:  # Avoid duplicate sublists
                        merged_dict[key].append(sublist)
            else:
                merged_dict[key] = values

        self.protein2gene =  {k: list(v) for k, v in merged_dict.items()}

    def _add_rxn2protein_information(self, rxn2protein_to_add:Dict[str, Dict]) -> None:
        merged_dict = self.rxn2protein.copy()  # Start with dict1 to avoid modifying it

        for rxn_id, attr_dict in rxn2protein_to_add.items():
            if rxn_id not in merged_dict:
                merged_dict[rxn_id] = attr_dict.copy()
            else:
                for enz_id, enz_dict in attr_dict.items():
                    # Ensure the enzyme entry exists
                    merged_dict[rxn_id].setdefault(enz_id, {})

                    for attr, value in enz_dict.items():
                        if isinstance(value, list):  # Handle protein_reaction_association values (avoid duplicates)
                            merged_dict[rxn_id][enz_id].setdefault(attr, [])
                            merged_dict[rxn_id][enz_id][attr].extend(
                                v for v in value if v not in merged_dict[rxn_id][enz_id][attr]
                            )
                        else:  # Handle numeric or single-value attributes (overwrite)
                            merged_dict[rxn_id][enz_id][attr] = value

        self.rxn2protein = merged_dict

    def _add_reaction_to_enzyme(self, model, enzyme:Union[Enzyme, EnzymeComplex], rxn_id:str, kcat:dict):
        # check if catalytic event already exists:
        catal_event_id = f'CE_{rxn_id}'
        if catal_event_id in enzyme.catalytic_events:
            return

        if catal_event_id in model.catalytic_events:
            ce = model.catalytic_events.get_by_id(catal_event_id)
            ce.add_enzymes({enzyme:kcat})
            enzyme.add_catalytic_event(ce, kcat)
        else:
            enzyme.create_catalytic_event(rxn_id=rxn_id, kcats=kcat)
            model.add_catalytic_events([enzyme.catalytic_events.get_by_id(f'CE_{rxn_id}')])

    def _enzyme_is_enzyme_complex(self, protein_reaction, enzyme_id):
        # make sure default enzyme id is not recognized as enzyme_complex
        if 'Enzyme_' in enzyme_id:
            sorted_enzyme_id = enzyme_id.replace('Enzyme_', '')
        else:
            #make sure the enzyme complex is in the right order, otherwise match will fail
            sorted_enzyme_id = sorted(enzyme_id.split('_'))

        return any(
            [all(
                [sorted(pr) == sorted_enzyme_id and  #enzyme should take part in enzyme complex
                 len(pr)>1] # enzyme complex needs to have at least 2 proteins
            ) for pr in protein_reaction]
        )

    def _get_model_genes_from_enzyme(self, enzyme_id: str, model: Model) -> list:
        """
           Retrieves genes associated with the specified enzyme from the model.

           Args:
               enzyme_id (str): The identifier of the enzyme.
               model (Model): The model containing gene information.

           Returns:
               list: A nested list of gene objects associated with the enzyme.
           """
        if enzyme_id not in self.protein2gene.keys(): return []
        gene_id_list = self.protein2gene[enzyme_id]
        gene_list = []
        for genes_or in gene_id_list:
            genes_and_list = []
            # print(genes_or, enzyme_id)
            for gene_and in genes_or:
                #check if there is an and relation (then gene and should be a list and not a string)
                if isinstance(gene_and, list):
                    for gene in gene_and:
                        if gene not in model.genes:
                            model.genes.append(Gene(gene))
                else:
                    gene = gene_and
                    if gene not in model.genes:
                        model.genes.append(Gene(gene))
                genes_and_list.append(model.genes.get_by_id(gene))
            gene_list.append(genes_and_list)
        return gene_list

    # def __getstate__(self):
    #     # Return the state to be pickled
    #     state = self.__dict__.copy()
    #     # Handle any non-serializable attributes here
    #     return state
    def change_kcat_values(self, rxn_id:str, enzyme_id:str,
                           kcat_f_b:dict[Literal['f', 'b'], float]
                           ) -> None:
        """
            Updates the kcat values for a specific enzyme associated with a given reaction.

            This method modifies the `rxn2protein` dictionary by either adding or updating
            the kcat values (forward and backward) for the specified enzyme in the context
            of the provided reaction.

            Args:
                rxn_id (str): The identifier for the reaction.
                enzyme_id (str): The identifier for the enzyme.
                kcat_f_b (dict[Literal['f', 'b'], float]): A dictionary containing the
                    forward ('f') and backward ('b') kcat values for the enzyme.
                    For example: {'f': 1.5, 'b': 0.8}.

            Returns:
                None: This method updates the instance's `rxn2protein` attribute in-place.

            Example:
                ```python
                model.change_kcat_values(
                    rxn_id="RXN001",
                    enzyme_id="ENZYME123",
                    kcat_f_b={'f': 2.0, 'b': 1.0}
                )
                ```
            """
        self.rxn2protein.setdefault(rxn_id, {}).update({enzyme_id:kcat_f_b})

    def __setstate__(self, state):
        # Restore state from the unpickled state
        self.__dict__.update(state)
        # Handle any attributes that require initialization or special handling here

class TransEnzymeSector(EnzymeSector):
    DEFAULT_MOL_MASS = 4.0590394e05  # default E. coli ribosome molar mass [g/mol]
    BIOMASS_RXNID = Config.BIOMASS_REACTION

    # class with all the information on enzymes related to translation/ribosomes (which is linearly related)
    def __init__(
        self,
        tps_0,
        tps_mu,
        id_list=[BIOMASS_RXNID],
        mol_mass=None,
        configuration=Config,
    ):
        super().__init__(id_list, mol_mass, configuration)
        BIOMASS_RXNID = configuration.BIOMASS_REACTION
        if self.mol_mass is None:
            self.mol_mass = [self.DEFAULT_MOL_MASS]
        # id_list only contains the model identifier of the biomass formation reaction
        # mol_mass is the molar mass of a fictional ribosome, if assigned, the computed enzyme concentration resamble
        # ribosome concentrations
        self.tps_0 = tps_0  # amount of protein allocated to the translational sector at zero growth (g_p/g_cdw)
        self.tps_mu = tps_mu  # amount of protein allocated to the translational sector at zero growth (g_p/g_cdw)'
        self.id = "TranslationalProteinSector"

        # *1000 to convert units from g/g_cdw to mg/g_cdw
        self.tps_0_coeff = self.tps_0[0] * 1e3
        self.slope = None
        self.set_slope()
        self.intercept = None
        self.set_intercept()

        # for easy removal
        self.variables = []
        self.constraints = []

    def add(self, model):
        print("Add translational protein sector\n")
        if self.mol_mass is None:
            self.mol_mass = [self.DEFAULT_MOL_MASS]
        self.set_tps_0_coeff()
        return self.add_sector(
            model=model, slope=self.tps_mu[0] * 1e3, intersect=self.tps_0_coeff
        )

    def set_tps_0_coeff(self):
        # *1000 to convert units from g/g_cdw to mg/g_cdw
        self.tps_0_coeff = self.tps_0[0] * 1e3

    def set_slope(self):
        # *1000 to convert units from g/g_cdw to mg/g_cdw
        self.slope = self.tps_mu[0] * 1e3  # / (self.mol_mass[0] * 1e-3)

    def set_intercept(self):
        self.set_tps_0_coeff()
        self.intercept = self.tps_0_coeff


class UnusedEnzymeSector(EnzymeSector):
    DEFAULT_MOL_MASS = 3.947778784340140e04  # mean enzymes mass E.coli [g/mol]

    # class with all the information on 'excess' proteins (linear dependent on substrate uptake rate)
    def __init__(self, ups_0, ups_mu, id_list, mol_mass=None, configuration=Config):
        super().__init__(id_list, mol_mass, configuration)
        if self.mol_mass is None:
            self.mol_mass = [self.DEFAULT_MOL_MASS]
        # id_list only contains the model identifier of the substrate uptake
        # mol_mass is the molar mass of a concentration unit of the unused protein sector
        self.ups_0 = ups_0  # amount of protein allocated to the excess enzyme sector at zero substrate uptake
        if isinstance(ups_mu, list):
            self.ups_mu = ups_mu[0]  # slope of linear relation with growth
        else:
            self.ups_mu = ups_mu
        # *1000 to convert units from g/g_cdw to mg/g_cdw
        self.ups_0_coeff = self.ups_0[0] * 1e3
        self.id = 'UnusedEnzymeSector'

        self.slope = self.ups_mu * 1e3
        self.intercept = self.ups_0_coeff

        # for easy removal
        self.variables = []
        self.constraints = []

    def add_to_model(self, model):
        print("Add unused protein sector\n")
        if self.mol_mass is None:
            self.mol_mass = [self.DEFAULT_MOL_MASS]
        return self.add_sector(
            model=model, slope=self.ups_mu * 1e3, intersect=self.ups_0_coeff
        )


class CustomSector(EnzymeSector):
    DEFAULT_ENZYME_MOL_MASS = 3.947778784340140e04  # mean enzymes mass E.coli [g/mol]

    # class with all the information on a custom protein allocation secotr (linear dependent on a user defined model variable)
    def __init__(
        self, id_list, name, cps_0, cps_s, mol_mass=None, configuration=Config
    ):
        super().__init__(id_list, mol_mass, configuration)
        # id_list containts the name of the reaction on which the custom protein sector is linearly dependent
        self.name = name  # name of the protein sector
        self.lin_rxn_id = id_list[
            0
        ]  # model identifier of the reaction from which the custom protein sector is linearly dependent
        self.cps_0 = cps_0  # intercept of the linear function describing protein allocation of the custom protein sector
        self.cps_s = cps_s  # slope of the linear function describing protein allocation of the custom protein sector
        self.cps_0_coeff = self.cps_0[0] * 1e3
        self.id = f"CustomProteinSector_{self.name}"

        # *1000 to convert units from g/g_cdw to mg/g_cdw
        self.slope = self.cps_s[0] * 1e3
        self.intercept = self.cps_0_coeff

        # for easy removal
        self.variables = []
        self.constraints = []

    def add_to_model(self, model):
        print(f"Add custom protein sector {self.name}\n")
        return self.add_sector(
            model=model, slope=self.cps_mu[0] * 1e3, intersect=self.cps_0_coeff
        )
