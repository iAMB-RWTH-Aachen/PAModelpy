from warnings import warn
from copy import copy, deepcopy
from cobra import Object

from .Enzyme import Enzyme
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

    def __copy__(self) -> 'Sector':
        """ Copy the Sector
        :return: Sector:
        A new Sector that is a copy of the original Sector
        """

        cop = copy(super(Sector, self))
        return cop

    def __deepcopy__(self, memo: dict) -> 'Sector':
        """ Copy the Sector with memo

        :param: memo:dict:
        Automatically passed parameter

        :return: Sector:
        A new Sector that is a copy of the original Sector with memo
        """

        cop = deepcopy(super(Sector, self), memo)
        return cop

class EnzymeSector(Sector):
    DEFAULT_MOL_MASS = 3.947778784340140e04  # mean enzymes mass E.coli [g/mol]
    def __init__(self, id_list, mol_mass, configuration = Config):
        self.id_list = id_list #list with reaction ids which are associated with the specific sector
        self.mol_mass = mol_mass #dict with molar masses for each enzyme in the id_list (g/mol)
        self.variables = []
        self.constraints = []
        self.model = None
        self.id = ''
        self.slope = 0
        self.intercept = 0

        self.TOTAL_PROTEIN_CONSTRAINT_ID = configuration.TOTAL_PROTEIN_CONSTRAINT_ID


class ActiveEnzymeSector(Sector):
    DEFAULT_MOL_MASS = 3.947778784340140e04  # mean enzymes mass E.coli [g/mol]
    #class with all the information on enzymes related to metabolic reactions
    def __init__(self, rxn2protein: dict, configuration:Config =Config):
        """_summary_

        Parameters
        ----------
        rxn2protein : nested dict
            reaction id, enzymes_dict key, value pairs for each reaction in the active_enzyme sector.
             The enzymes_dict contains the enzyme identifiers of the enzymes which are related to the specific reaction
             as keys and a dict with information about the enzyme as values. The information included in this dict is:
              turnover number for forward and backward reaction [1/s], molar mass of the enzyme [mol/g], gene identifiers
              related to the enzyme, with which other enzymes it forms a complex.
              example: {R1:
                            {E1:
                                {'f': forward kcat, 'b': backward kcat, 'molmass': molar mass, 'genes': [G1, G2],
                                 'complex_with': 'E2'},
                            E2:
                                {'f': forward kcat, 'b': backward kcat, 'molmass': molar mass, 'genes': [G3, G4],
                                 'complex_with': 'E1'}}

        configuration: Config object, optional
                    Information about general configuration of the model including identifier conventions.
                    Default as defined in the `PAModelpy.configuration` script for the E.coli iML1515 model.
        """
        self.TOTAL_PROTEIN_CONSTRAINT_ID = configuration.TOTAL_PROTEIN_CONSTRAINT_ID

        self.rxn2protein = rxn2protein
        # ID of the sector
        self.id = 'ActiveEnzymeSector'
        self.model = None

        # for easy removal
        self.variables = []
        self.constraints = []

    def add(self, model):
        self.model = model
        #make and enzymatic constraint for an existing stoichiometric reaction in the form of: e_i=v_i*kcat_i
        # c.f. the GECKO formulation of Sanchez et al. 2017
        print('Add active protein sector\n')
        #check if the constraint exists
        totprot, tpc_metabolite = self.get_tpc_metabolite(model)
        if not totprot:
            warn('Total protein constraint is not in the model. Skip the Active Enzyme Sector')
            return model

        rxn2protein = self.rxn2protein.copy()

        for rxn_id, enzymes in self.rxn2protein.items():
            # extract reaction from model
            if rxn_id not in model.reactions:
            #     reaction = model.reactions.get_by_id(rxn_id)
            # else:
                try:
                    reaction = model.reactions.query(rxn_id)
                    for rxn in reaction:
                        rxn2protein[rxn.id] = enzymes

                        del rxn2protein[rxn_id]
                except:
                    warn(f'Reaction {rxn_id} is not in the model, this reaction will be skipped')
                    if rxn_id in rxn2protein.keys():
                        del rxn2protein[rxn_id]

        for rxn_id, enzymes in rxn2protein.items():
            reaction = model.reactions.get_by_id(rxn_id)
            # skip the reaction if a problem is encountered in check_kcat_values()
            consistent = self.check_kcat_values(model, reaction, enzymes)
            if not consistent: continue

            for enzyme_id, enzyme_dict in enzymes.items():
                kcat = dict((k, v) for k, v in enzyme_dict.items() if k == 'f' or k == 'b')

                # get molar mass of enzyme or replace with default value
                if 'molmass' in enzyme_dict.keys():
                    molmass = enzyme_dict['molmass']
                else:
                    molmass = self.DEFAULT_MOL_MASS

                # check if there already exists an Enzyme object for the EC numbers associated to this reaction,
                # otherwise create one and store all information
                if enzyme_id in model.enzymes:
                    enzyme = model.enzymes.get_by_id(enzyme_id)

                    #check if catalytic event already exists:
                    catal_event_id = f'CE_{rxn_id}'

                    if catal_event_id in model.catalytic_events:
                        ce = model.catalytic_events.get_by_id(catal_event_id)
                        enzyme.add_catalytic_event(ce, kcat)
                    else:
                        enzyme.create_catalytic_event(rxn_id=rxn_id, kcats=kcat)
                        model.add_catalytic_events([enzyme.catalytic_events.get_by_id(f'CE_{rxn_id}')])
                else:
                    enzyme = Enzyme(
                        id = enzyme_id,
                        rxn2kcat= {rxn_id: kcat},
                        molmass = molmass
                    )
                    #in the add enzymes function the enzyme with corresponding catalytic events will be added to the model
                    #and connected to the reaction and total protein constraint

                    model.add_enzymes([enzyme])

                    #adding to the enzyme sector object for easy removal
                    self.constraints += [enzyme]
                    model.tpc += 1
                self.variables.append(enzyme.enzyme_variable)
        return model

    def check_kcat_values(self, model, reaction, kcat):
        """
        Function to check if the kcat values provided for an enzyme are consistent with the direction of the reaction
        Parameters
        ----------
        model: cobra.Model or PAModel
            model to which the kcat values should be added
        reaction: cobra.Reaction
            reaction which is catalyzed with the enzyme related to the kcats
        kcat: dict
            a dict with the turnover values for the forward and/or backward reaction for different enzymes [/s]
            {'E1':{'f': forward kcat, 'b': backward kcat}}

        Returns
        -------

        """
        if reaction.id not in model.reactions:
            warn('Reaction ' + reaction.id + ' not in the model. Skip enzyme constraint')
            return False
        # check if there is an EC number associated to each reaction, otherwise, use the reaction id
        if reaction.id not in self.rxn2protein.keys():
            self.rxn2protein[reaction.id] = reaction.id

         # check consistency between provided kcat values and reaction direction
        directions = []
        for val in kcat.values():# get all directions from the kcat dict
            directions += list(val.keys())
        kcats = []
        for val in kcat.values():# get all directions from the kcat dict
            kcats += list(val.values())

        rxn_dir = 'consistent'
        if reaction.lower_bound >= 0 and reaction.upper_bound > 0:
            # reaction is irreversible in the forward direction

            if not 'f' in directions:
                rxn_dir = 'inconsistent'

        elif reaction.lower_bound < 0 and reaction.upper_bound <= 0:
            # reaction is irreversible in the backward direction
            if not 'b' in directions:
                rxn_dir = 'inconsistent'
        else:
            # reaction is reversible
            if not 'f' in directions and 'b' in directions:
                rxn_dir = 'inconsistent'

        if rxn_dir == 'inconsistent':
            warn(reaction.id + ': reaction directionality does not match provided kcat values. Skip reaction')
            return False

        # check kcat values input
        for k in kcats:
            if k < 0:
                warn(['Turnover number for reaction "', reaction.id, '" is invalid. Skip for active enzyme sector'])
                return False

        return rxn_dir == 'consistent'

class TransEnzymeSector(EnzymeSector):
    DEFAULT_MOL_MASS = 4.0590394e05  # default E. coli ribosome molar mass [g/mol]
    BIOMASS_RXNID = Config.BIOMASS_REACTION
    #class with all the information on enzymes related to translation/ribosomes (which is linearly related)
    def __init__(self, tps_0, tps_mu,id_list= [BIOMASS_RXNID], mol_mass=None, configuration = Config):
        super().__init__(id_list, mol_mass, configuration)
        BIOMASS_RXNID = configuration.BIOMASS_REACTION
        if self.mol_mass is None:
            self.mol_mass = [self.DEFAULT_MOL_MASS]
        #id_list only contains the model identifier of the biomass formation reaction
        #mol_mass is the molar mass of a fictional ribosome, if assigned, the computed enzyme concentration resamble
        # ribosome concentrations
        self.tps_0 =tps_0 #amount of protein allocated to the translational sector at zero growth (g_p/g_cdw)
        self.tps_mu = tps_mu  # amount of protein allocated to the translational sector at zero growth (g_p/g_cdw)'
        self.id = 'TranslationalProteinSector'

        # *1000 to convert units from g/g_cdw to mg/g_cdw
        self.tps_0_coeff = self.tps_0[0] *1e3
        self.slope = None
        self.set_slope()
        self.intercept = None
        self.set_intercept()

        # for easy removal
        self.variables = []
        self.constraints = []

    def add(self, model):
        print('Add translational protein sector\n')
        if self.mol_mass is None:
            self.mol_mass = [self.DEFAULT_MOL_MASS]
        self.set_tps_0_coeff()
        return self.add_sector(model= model, slope = self.tps_mu[0]* 1e3, intersect=self.tps_0_coeff)

    def set_tps_0_coeff(self):
        # *1000 to convert units from g/g_cdw to mg/g_cdw
        self.tps_0_coeff = self.tps_0[0] *1e3

    def set_slope(self):
        # *1000 to convert units from g/g_cdw to mg/g_cdw
        self.slope = self.tps_mu[0] *1e3#/ (self.mol_mass[0] * 1e-3)

    def set_intercept(self):
        self.set_tps_0_coeff()
        self.intercept = self.tps_0_coeff

class UnusedEnzymeSector(EnzymeSector):
    DEFAULT_MOL_MASS = 3.947778784340140e04  # mean enzymes mass E.coli [g/mol]
    #class with all the information on 'excess' proteins (linear dependent on substrate uptake rate)
    def __init__(self, ups_0, ups_mu ,id_list, mol_mass=None, configuration = Config):
        super().__init__(id_list, mol_mass, configuration)
        if self.mol_mass is None:
            self.mol_mass = [self.DEFAULT_MOL_MASS]
        # id_list only contains the model identifier of the substrate uptake
        # mol_mass is the molar mass of a concentration unit of the unused protein sector
        self.ups_0 = ups_0 #amount of protein allocated to the excess enzyme sector at zero substrate uptake
        if isinstance(ups_mu, list):
            self.ups_mu = ups_mu[0]  # slope of linear relation with growth
        else:
            self.ups_mu = ups_mu
        # *1000 to convert units from g/g_cdw to mg/g_cdw
        self.ups_0_coeff = self.ups_0[0]*1e3

        self.slope = self.ups_mu*1e3
        self.intercept = self.ups_0_coeff

        # for easy removal
        self.variables = []
        self.constraints = []

    def add_to_model(self, model):
        print('Add unused protein sector\n')
        if self.mol_mass is None:
            self.mol_mass = [self.DEFAULT_MOL_MASS]
        return self.add_sector(model= model, slope = self.ups_mu*1e3, intersect=self.ups_0_coeff)

class CustomSector(EnzymeSector):
    DEFAULT_ENZYME_MOL_MASS = 3.947778784340140e04  # mean enzymes mass E.coli [g/mol]

    #class with all the information on a custom protein allocation secotr (linear dependent on a user defined model variable)
    def __init__(self, id_list, name, cps_0,cps_s, mol_mass=None, configuration=Config):
        super().__init__(id_list, mol_mass, configuration)
        #id_list containts the name of the reaction on which the custom protein sector is linearly dependent
        self.name = name #name of the protein sector
        self.lin_rxn_id = id_list[0] #model identifier of the reaction from which the custom protein sector is linearly dependent
        self.cps_0 = cps_0 #intercept of the linear function describing protein allocation of the custom protein sector
        self.cps_s = cps_s #slope of the linear function describing protein allocation of the custom protein sector
        self.cps_0_coeff = self.cps_0[0] *1e3
        self.id = f'CustomProteinSector_{self.name}'

        #*1000 to convert units from g/g_cdw to mg/g_cdw
        self.slope = self.cps_s[0]*1e3
        self.intercept = self.cps_0_coeff

        # for easy removal
        self.variables = []
        self.constraints = []

    def add_to_model(self, model):
        print(f'Add custom protein sector {self.name}\n')
        return self.add_sector(model= model, slope = self.cps_mu[0]*1e3, intersect=self.cps_0_coeff)
