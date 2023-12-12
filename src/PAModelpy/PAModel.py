#cobra tools
import cobra
from cobra import Model, DictList, Reaction, Metabolite, Solution
from cobra.io import load_model
from cobra.util.context import get_context
#type checking
from optlang.symbolics import Zero
from optlang.interface import Objective
from typing import List, Optional, Union, Dict, Iterable
#other
from functools import partial
import warnings
import pandas as pd
from copy import copy, deepcopy
import inspect

# sys.path.append('../')
from .EnzymeSectors import ActiveEnzymeSector, TransEnzymeSector, UnusedEnzymeSector, CustomSector, Sector
from .CatalyticEvent import CatalyticEvent
from .Constraints import Constraint
from .Enzyme import Enzyme
from .configuration import Config


class PAModel(Model):
    """Class representation for a cobra model extended with enyzme kinetics as published in Alter et al. (2021)
                Parameters
                ----------
                id_or_model: str or Model
                    String to use as model id, or actual model to base new model one.
                    If string, it is used as input to load a model from. If model, a new model object is
                    instantiated with the same properties as the original model (default None).
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
                configuration: Config object, optional
                    Information about general configuration of the model including identifier conventions.
                    Default as defined in the `PAModelpy.configuration` script for the E.coli iML1515 model.

                Attributes
                ----------
                p_tot : float
                    The fraction of biomass allocated to proteins (units: g_prot/g_cdw)
                reactions : DictList
                    A DictList where the key is the reaction identifier and the value a
                    Reaction
                metabolites : DictList
                    A DictList where the key is the metabolite identifier and the value a
                    Metabolite
                genes : DictList
                    A DictList where the key is the gene identifier and the value a
                    Gene
                groups : DictList
                    A DictList where the key is the group identifier and the value a
                    Group
                enzymes : DictList
                    A DictList where the key is the enzyme identifier and the value an
                    Enzyme
                enzyme_variables: DictList

                catalytic_events: DictList

                sector_constraints: dict


                sectors : DictList
                    A DictList where the key is the sector identifier and the value an
                    EnyzmeSector
                """
    TOTAL_PROTEIN_CONSTRAINT_ID = Config.TOTAL_PROTEIN_CONSTRAINT_ID
    P_TOT_DEFAULT = Config.P_TOT_DEFAULT  # g_protein/g_cdw
    CO2_EXHANGE_RXNID = Config.CO2_EXHANGE_RXNID
    GLUCOSE_EXCHANGE_RXNID = Config.GLUCOSE_EXCHANGE_RXNID
    BIOMASS_REACTION = Config.BIOMASS_REACTION

    def __init__(self, id_or_model: Union[str, "Model", None] = None,
                 name: Optional[str] = None,
                 p_tot: Optional[float] = Config.P_TOT_DEFAULT,
                 sensitivity: bool = True,
                 active_sector: Optional[ActiveEnzymeSector]=None,
                 translational_sector: Optional[TransEnzymeSector]=None,
                 unused_sector: Optional[UnusedEnzymeSector]=None,
                 custom_sectors: Union[List, CustomSector] =None,
                 configuration = Config):
        """Constants"""
        self.TOTAL_PROTEIN_CONSTRAINT_ID = configuration.TOTAL_PROTEIN_CONSTRAINT_ID
        self.P_TOT_DEFAULT = configuration.P_TOT_DEFAULT  # g_protein/g_cdw
        self.CO2_EXHANGE_RXNID = configuration.CO2_EXHANGE_RXNID
        self.GLUCOSE_EXCHANGE_RXNID = configuration.GLUCOSE_EXCHANGE_RXNID
        self.BIOMASS_REACTION = configuration.BIOMASS_REACTION

        self.configuration = configuration
        """Initialize the Model."""
        if isinstance(id_or_model, str):
            id_or_model = load_model(id_or_model)

        super().__init__(id_or_model=id_or_model, name=name)
        self.m_model = id_or_model.copy() #save a copy of the original m_model

        self.p_tot = p_tot # fraction of biomass allocated to proteins (units: g_prot/g_cdw)
        self.enzymes = DictList() # a list of Enzyme.Enzyme
        self.enzyme_variables = DictList() # a list of variables related to the enzymes (Variables.EnzymeVariable)
        self.catalytic_events = DictList() # a list of objects storing the relation of a single reaction to the enzymes necessary for the catalysis of that reaction (Variables.CatalyticEvent)
        self.enzyme_constraints = {} # a dict with enzyme constraint id (format: 'EC_{ecnmbr}_{direction}'), optlang.Constraint (enzymes) key, value pairs
        self.sectors = DictList() # a list of EnzymeSectors (protein sectors)
        self.sector_constraints = {} # a dict with sector constraint id, optlang.Constraint (enzymes) key, value pairs for constraints related to the protein sectors
        self.tpc = 0 # counter of number of CatalyticEvents which contribute to the total protein constraint
        self.sensitivity = sensitivity
        self.capacity_sensitivity_coefficients = pd.DataFrame() #dataframe to store the result of the sensitivity analysis (capacity sensitivity coefficients for each constraint). The sensitivity coefficients are splitted on LB, UB and the different sectors
        self.enzyme_sensitivity_coefficients = pd.DataFrame() #dataframe to store the result of the sensitivity analysis (sensitivity coefficients for each constraint). The sensitivity coefficients are splitted on LB, UB and the different sectors


        #initialize the model
        print(f'Setting up the proteome allocation model {self.id}\n')
        self.add_total_protein_constraint(p_tot)

        if sensitivity: #perform sensitivity analysis when the model is run
            self.add_lb_ub_constraints()

        for sector in [active_sector, translational_sector, unused_sector, custom_sectors]:
            if sector is not None:
                self.add_sectors([sector])
        print(f'Done with setting up the proteome allocation model {self.id}\n')


    @property
    def total_protein_fraction(self):
        return self.p_tot

    @total_protein_fraction.setter
    def total_protein_fraction(self, p_tot):
        self.change_total_protein_constraint(p_tot)

    @property
    def translational_enzymes(self):
        return self.sectors.get_by_id('TranslationalEnzymeSector')

    @translational_enzymes.setter
    def tranlational_enzymes(self, slope:float, intercept:float, lin_rxn_id: str = 'BIOMASS_Ec_iML1515_WT_75p37M'):
        #input is in g/gDW
        self.change_sector_parameters(self.translational_enzymes, slope, intercept, lin_rxn_id)

    @property
    def unused_enzymes(self):
        return self.sectors.get_by_id('UnusedEnzymeSector')

    @unused_enzymes.setter
    def unused_enzymes(self, slope: float, intercept: float, lin_rxn_id: str = 'EX_glc__D_e_b'):
        self.change_sector_parameters(self.unused_enzymes, slope, intercept, lin_rxn_id)

    @property
    def stoichiometric_matrix(self):
        #TODO check solver, only works for gurobi solver
        matrix = self.problem.solver.getA()
        return matrix.toarray()

    def add_enzymes(self, enzyme_list: list) -> None:
        """Add new enzymes to a model.
            Adapted from Cobra.core.model.add_reactions and Cobra.core.model.add_metabolites()

            Will add a DictList of enzymes to the model object and add new
            variables accordingly.
            For each enzyme-associated reaction a constraint in each direction 
            is added to the model.
            The change is reverted upon exit when using the model as a context.
            :param
            enzyme_list : list or Enzyme.
                   A list of `Enzyme` objects. If it isn't an iterable
                   container, the enzyme will be placed into a list.
               """
        
        def existing_filter(enz: Enzyme) -> bool:
            """Check if the enzyme does not exist in the model.
            Parameters
            ----------
            enz: Variables.Enzyme
            Returns
            -------
            bool
                False if enzyme exists, True if it doesn't.
                If the enzyme exists, will log a warning.
            """
            if enz.id in self.enzymes:
                warnings.warn(f"Ignoring enzyme '{enz.id}' since it already exists.")
                return False
            return True
        
        def parameter_filter(enz: Enzyme) -> bool:
            """Check enzyme parameters are consistent with the model

            Parameters
            ----------
            enz : Enzyme
                Enzyme object inferred from cobra.core.reaction

            Returns
            -------
            bool
                enzyme parameter validity
            """

            # extract parameters from Enzyme
            molmass = enz.molmass
            
            # check molar mass
            if molmass <= 0:
                warnings.warn(f'Molar mass for {enz.id} is invalid: {molmass}')
                return False

            # check if enzyme objects are valid
            for catalytic_event in enz.catalytic_events:
                # extract parameters from enzyme object
                kcats = catalytic_event.kcats
                rxn_id = catalytic_event.rxn_id
                
                # check if reaction associated with the enzyme object exists
                if rxn_id not in self.reactions:
                    warnings.warn('Reaction '+ rxn_id +' not in the model. Skip enzyme constraint')
                    return False

                # check kcat values
                for kcatdict in kcats.values():
                    for kcat in kcatdict.values():
                        if kcat < 0:
                            # invalid kcat value
                            warnings.warn('Turnover number for reaction "'+rxn_id+'" is invalid. Skip for active enzyme sector')
                            return False
                    
                # extract reaction from model
                reaction = self.reactions.get_by_id(rxn_id)    
                for kcats in kcats.values():
                    # check consistency between provided kcat values and reaction direction
                    if self.sensitivity:
                        lower_bound = -self.constraints[f'{rxn_id}_lb'].ub
                        upper_bound = self.constraints[f'{rxn_id}_ub'].ub
                    else:
                        lower_bound = reaction.lower_bound
                        upper_bound = reaction.upper_bound
                    if lower_bound >= 0 and upper_bound > 0:
                        # reaction is irreversible in the forward direction
                        if 'f' not in kcats: # or 'b' in kcats:
                            warnings.warn(rxn_id+': Inconsistencies between the reaction reversibility and the provided kcat values')
                            return False
                    elif lower_bound < 0 and upper_bound <= 0:
                        # reaction is irreversible in the backward direction
                        if 'b' not in kcats or 'f' in kcats:
                            warnings.warn(rxn_id+': Inconsistencies between the reaction reversibility and the provided kcat values')
                            return False
                    else:
                        # reaction is reversible
                        if 'f' not in kcats:# or 'b' not in kcats:
                            warnings.warn(rxn_id+': Inconsistencies between the reaction reversibility and the provided kcat values')
                            return False

            return True

        #check if the input is a list or an enzyme
        if not hasattr(enzyme_list, "__iter__"):
            enzyme_list = [enzyme_list]
        for enz in enzyme_list:
            if not isinstance(enz, Enzyme):
                raise AttributeError('The input was neither an Iterable nor an Enzyme. Please provide only (lists of) Enzyme objects.')

        # First check whether the enzymes exist in the model.
        pruned = filter(existing_filter, enzyme_list)

        # check if the provided parameters (kcats, rxn_id) match with the model reaction
        pruned = DictList(filter(parameter_filter, pruned))

        # add context manager
        context = get_context(self)

        for enzyme in pruned:
            # make enzyme aware of the model
            enzyme._model = self
            if context:
                context(partial(setattr, enzyme, "_model", None))

        # add enzymes to the model
        self.enzymes += pruned

        if context:
            context(partial(self.enzymes.__isub__, pruned))

        # populate solver
        for enzyme in pruned:
            #add constraint related to the enzyme
            self.add_enzyme_constraints([f'EC_{enzyme.id}_f', f'EC_{enzyme.id}_b'])
            #store the constraints connected to the model in the enzyme object
            new_constraints = {}
            new_constraints[f'EC_{enzyme.id}_f'] = self.enzyme_constraints[f'EC_{enzyme.id}_f']
            new_constraints[f'EC_{enzyme.id}_b'] = self.enzyme_constraints[f'EC_{enzyme.id}_b']
            enzyme._constraints = new_constraints

            #add enzyme variable
            enzyme_variable = enzyme.enzyme_variable
            if enzyme_variable.id not in self.variables:
                    # make enzyme variable aware of the model.
                    # The enzyme variable will be added as variable to the model in the model.setter magic function
                    enzyme_variable.model = self
                    # self.enzyme_variables.append(enzyme_variable)
                    if context:
                        context(partial(setattr, enzyme_variable, "_model", None))

            # connect the enzyme to total protein constraint
            enzyme_variable = self.enzyme_variables.get_by_id(enzyme.id)
            self.constraints[self.TOTAL_PROTEIN_CONSTRAINT_ID].set_linear_coefficients(
                { #*1e-6 to make calculations more accurate (solver accuracy)
                enzyme_variable.forward_variable: enzyme.molmass * 1e-6,
                enzyme_variable.reverse_variable: enzyme.molmass * 1e-6,
                }
            )
            self.tpc+=1

            # add enzyme constraints to the model and link to enzyme and reaction variables
            #get the enzyme variable
            enzyme_var_model = self.enzyme_variables.get_by_id(enzyme.id)
            #connect enzyme variable to its upper and lower bound
            if enzyme.upper_bound > self.p_tot*1e3:
                ub = enzyme.upper_bound
            else:
                ub = self.p_tot*1e3
            self, enzyme = self.make_enzyme_min_max_constraint(self, enzyme, lower_bound=enzyme.lower_bound, upper_bound=ub)

            # add the enzyme to the interface between the enzyme constraints, enzyme variables and the reaction
            for catalytic_event in enzyme.catalytic_events:
                if catalytic_event not in self.catalytic_events:
                    # make catalytic event aware of the model.
                    # The catalytic event will be added configured in the model.setter magic function
                    catalytic_event.model = self
                    self.catalytic_events += [catalytic_event]
                    if context:
                        context(partial(setattr, catalytic_event, "_model", None))
                else:
                    # add enzyme to the catalytic event
                    catalytic_event_model = self.catalytic_events.get_by_id(catalytic_event.id)
                    #remove block on the reaction of no enzyme constraint if this is present
                    ce_constraints = catalytic_event_model.constraints.copy()
                    for name, constraint in ce_constraints.items():
                        if 'no_enzyme' in name:
                            del catalytic_event_model.constraints[name]
                            self.remove_cons_vars([constraint])

                    catalytic_event_model.add_enzymes({enzyme: enzyme.rxn2kcat[catalytic_event.rxn_id]})

                    # replace the catalytic event with the already existing event from the model
                    enzyme.catalytic_events._replace_on_id(catalytic_event_model)

            #connect the enzyme to each of the reactions it is associated with
            for rxn, kcatdict in enzyme.rxn2kcat.items():
                # link enzymatic variable to reaction via the turnover number kcat
                for direction, kcat in kcatdict.items():
                    # check direction
                    if direction != 'f' and direction != 'b':
                        warnings.warn(['Invalid kcat direction encountered for ', catalytic_event.id, direction])
                        continue

                    # create enzyme constraint for the reaction if not existent already
                    constraint_id = 'EC_' + enzyme.id + '_' + direction
                    if constraint_id not in self.enzyme_constraints.keys():
                        self.add_enzyme_constraints([constraint_id])

                    # kcat is used as a coefficient for the enzyme concentration
                    # Get existent forward/reverse variables
                    coeff = kcat*3600*1e-6 #3600 to convert to /h to /s *1e-6 to make calculations more accurate
                    if direction == 'f':
                        self.constraints[constraint_id].set_linear_coefficients(
                            {enzyme_var_model.forward_variable: -1,
                              self.reactions.get_by_id(rxn).forward_variable: 1/coeff
                            })
                    elif direction == 'b':
                        self.constraints[constraint_id].set_linear_coefficients(
                            {enzyme_var_model.reverse_variable: -1,
                             self.reactions.get_by_id(rxn).reverse_variable: 1/coeff
                            })

                    # make reaction-enzyme interface and the enzyme variable aware of its participation in the constraint
                    catalytic_event_model.constraints[enzyme.id] = self.constraints[constraint_id]
                    enzyme_var_model.constraints[constraint_id] = self.constraints[constraint_id]
                    self.solver.update()

    def add_sectors(self, sectors: List = None):
        """
        Adds sector variables to the model and add these tot the total protein constraint.

        Parameters
        ----------
        sectors: list
            list of PAModelpy.EnzymeSectors to add to the model

        """
        self.sectors += sectors
        for sector in sectors:
            #different method to add the active enzyme_sector
            if isinstance(sector, ActiveEnzymeSector):
                self = sector.add(self)
            else:
                self.add_sector(sector)

    def add_sector(self, sector):
        """
                Adds the sector variable for a specific sector to the model and add this tot the total protein constraint.
                Also stores the sector variables in the model attributes

                Parameters
                ----------
                sectors: PAModelpy.EnzymeSectors to add to the model

                """
        print(f'Add the following protein sector: {sector.id}\n')
        #make sector aware of the model
        sector.model = self

        lin_rxn = self.reactions.get_by_id(sector.id_list[0])


        #relate the sector to the total protein if it is there
        if self.TOTAL_PROTEIN_CONSTRAINT_ID in self.constraints.keys():
            # add parts of constraint corresponding to the enzyme sector to the total_protein_constraint
            # 1. subtract the intercept value from the sum of protein (Total_protein == Etot - phi_sector
            tpc_ub = self.constraints[self.TOTAL_PROTEIN_CONSTRAINT_ID].ub
            self.constraints[self.TOTAL_PROTEIN_CONSTRAINT_ID].ub = tpc_ub - sector.intercept

            #2. link flux to enzyme concentration
            # link enzyme concentration in the sector to the total enzyme concentration
            self.constraints[self.TOTAL_PROTEIN_CONSTRAINT_ID].set_linear_coefficients({
                lin_rxn.forward_variable: sector.slope,
                lin_rxn.reverse_variable: -sector.slope
            })
        else:
            # constraint corresponding to the enzyme sector
            constraint = self.problem.Constraint(Zero, name=sector.id, lb=0, ub=0)
            self.add_sector_constraints([constraint])

            # add variable which represents the concentration of enzymes in this enzyme sector
            var = self.problem.Variable('R_' + sector.id)
            self.add_cons_vars(var)

            # link flux to enzyme concentration
            self.constraints[sector.id].set_linear_coefficients({
                var: 1,
                lin_rxn.forward_variable: -sector.slope,
                lin_rxn.reverse_variable: sector.slope
            })

            # Make a variable which represents the intercept
            var_0 = self.problem.Variable('R0_' + sector.id)
            var_0.lb, var_0.ub = sector.intercept, sector.intercept
            self.add_cons_vars(var_0)

            #add intercept to the constraint definition
            self.constraints[sector.id].set_linear_coefficients({
                var_0: -1})

            #save the intercept variable for easy removal of the sector
            sector.variables += [var_0]

            # save variable and constraint for easy removal of the sector
            sector.constraints += [self.sector_constraints[sector.id]]
            sector.variables += [var]

        #add sector to sectorlist is it isn't included already
        if not self.sectors.has_id(sector.id):
            self.sectors += [sector]


    def add_catalytic_events(self, catalytic_events: Optional[Iterable]):
        """
        Add a new CatalyticEvent to the model.
        Will add a list of CatalyticEvent variables to the model object using the function defined in the
        CatalyticEvent object.

        Parameters
        ----------
        catalytic_events : list or variables.CatalyticEvent.
            A list of `variables.CatalyticEvent` objects. If it isn't an iterable
            container, the catalytic event will be placed into a list.

        """
        #check if the input is an iterable
        if not hasattr(catalytic_events, "__iter__"):
            catalytic_events = [catalytic_events]
        if len(catalytic_events) == 0:
            return None

        # First check whether the catalytic exist in the model
        catalytic_events = [x for x in catalytic_events if x.id not in self.catalytic_events]

        # if there already exists a catalytic event, the enzymes of the new event should be added to the old one
        present_events = [x for x in catalytic_events if x.id in self.catalytic_events]
        for event in present_events:
            enzkcatdict = {enz: enz.rxn2kcat[event.rxn_id] for enz in event.enzymes}
            model_event = self.catalytic_events.get_by_id(event.id)
            model_event.add_enzymes(enzkcatdict)

        self.catalytic_events += catalytic_events

        bad_ids = [
                m for m in catalytic_events if not isinstance(m.id, str) or len(m.id) < 1
            ]
        if len(bad_ids) != 0:
            raise ValueError(f"invalid identifiers in {repr(bad_ids)}")

        for event in catalytic_events:
            #make catalytic event aware of the model. This will automatically add the catalytic event to the model
            event.model = self

    def add_enzyme_constraints(self, constraint_list:Optional[list]):
       
        """Add new enzyme constraint to a model.
        Will add a list of constraints to the model object and add new
        constraints accordingly.
        The change is reverted upon exit when using the model as a context.
        Parameters
        ----------
        constraint_list : list, str or constraints.Constraint.
            A list of `constraints.Constraint` objects. If it isn't an iterable
            container, the constraint will be placed into a list. Also, a string with
            the constraint id can be provided. A constraint will be created before adding
            it to the model
        """
        if not hasattr(constraint_list, "__iter__"):
            constraint_list = [constraint_list]
        if len(constraint_list) == 0:
            return None

        #check wether the input is a  cobra.Metabolite or string and convert it to constraint
        for i, cons in enumerate(constraint_list):
            if isinstance(cons, Metabolite) or isinstance(cons, Constraint):
                constraint_list[i] = self.problem.Constraint(Zero, name=cons.id, lb=0, ub=0)
            if isinstance(cons, str):
                constraint = self.problem.Constraint(Zero, name=cons, lb=0, ub=0)
                constraint_list[i] = constraint

        # First check whether the metabolites exist in the model
        constraint_list = [x for x in constraint_list if x.name not in self.enzyme_constraints.keys()]
        bad_ids = [
            m for m in constraint_list if not isinstance(m.name, str) or len(m.name) < 1
        ]
        if len(bad_ids) != 0:
            raise ValueError(f"invalid identifiers in {repr(bad_ids)}")

        for x in constraint_list:
            x._model = self

        # from cameo ...
        to_add = []
        for constraint in constraint_list:
            if constraint not in self.enzyme_constraints.values():
                to_add += [constraint]
                self.enzyme_constraints[constraint.name] = constraint
        self.add_cons_vars(to_add)

        context = get_context(self)
        if context:
            context(partial(self.enzyme_constraints.__isub__, constraint_list))
            for x in constraint_list:
                # Do we care?
                context(partial(setattr, x, "_model", None))
        self.solver.update()


    def add_sector_constraints(self, constraint_list:Optional[list]):
        """Add new constraint related to an sector to a model.
               Will add a list of constraint to the model object and add new
               constraints accordingly.
               The change is reverted upon exit when using the model as a context.
               Parameters
               ----------
               constraint_list : list or constraints.Constraint.
                   A list of `constraints.Constraint` objects. If it isn't an iterable
                   container, the constraint will be placed into a list.
               """

        if not hasattr(constraint_list, "__iter__"):
            constraint_list = [constraint_list]
        if len(constraint_list) == 0:
            return None

        # First check whether the metabolites exist in the model
        constraint_list = [x for x in constraint_list if x.name not in self.sector_constraints.keys()]

        bad_ids = [
            m for m in constraint_list if not isinstance(m.name, str) or len(m.name) < 1
        ]
        if len(bad_ids) != 0:
            raise ValueError(f"invalid identifiers in {repr(bad_ids)}")

        for x in constraint_list:
            x._model = self

        # from cameo ...
        to_add = []
        for sect in constraint_list:
            if sect.name not in self.sector_constraints.keys():
                to_add += [sect]
                self.sector_constraints[sect.name] = sect

        self.add_cons_vars(to_add)

        context = get_context(self)
        if context:
            context(partial(self.sector_constraints.__isub__, constraint_list))
            for x in constraint_list:
                # Do we care?
                context(partial(setattr, x, "_model", None))
        self.solver.update()

    def add_total_protein_constraint(self, p_tot: Optional[float]=P_TOT_DEFAULT):
        """
        Function which adds the total protein constraint to the model.
        This limits the amount of available enzymes and thus the resulting fluxes
        constraint expression:
        Etot: sum(E) + E_translprot + E_unusedprot  == ptot - E_trsn_0 - E_ue_0

        :param
            p_tot: float, optional
                Fraction of biomass which consists of protein (g_protein/g_cdw)
                default: 0.258 (E.coli)
        """
        #check if we should add total protein constraint
        if isinstance(self.p_tot, bool) and isinstance(p_tot, bool):
            print('Total condition-dependent protein constraint is not added \n')
            return

        print('Add total condition-dependent protein constraint')
        if isinstance(p_tot, float):
            self.p_tot = p_tot

        print(f'\tTotal protein concentration: {p_tot} g/gDW\n')

        # check if there already is a total protein constraint in the solver
        if self.TOTAL_PROTEIN_CONSTRAINT_ID in self.constraints.keys():
            self.change_total_protein_constraint(p_tot)
        else:
            # Create the pseudometabolite associated with the constraint, this metabolite will be 'produced' in the enzyme 'reactions'
            # *1e3 to convert gp/gcdw to mg/gcdw for easier computing
            tpc_constraint = self.problem.Constraint(Zero, name=self.TOTAL_PROTEIN_CONSTRAINT_ID, lb=0,
                                                     ub=self.p_tot * 1e3)
            self.sector_constraints[self.TOTAL_PROTEIN_CONSTRAINT_ID] = tpc_constraint
            self.add_cons_vars([tpc_constraint])

    def add_reactions(self, reaction_list: Iterable[Reaction]) -> None:
        """Add reactions to the model.
        This method is superimposed upon the cobra.Model.add_reactions() function.
        As a new feature, it will add constraints to determine the lower and upper bound if a sensitivity analysis should
        be performed (which is determined by the model attribute: PAModel.sensitivity).
        Reactions with identifiers identical to a reaction already in the
        model are ignored.
        The change is reverted upon exit when using the model as a context.
        Parameters
        ----------
        reaction_list : list
            A list of `cobra.Reaction` objects
        """
        super().add_reactions(reaction_list=reaction_list)
        #add upper and lower bound constraints if you want to perform a sensitivity analysis
        if self.sensitivity:
            for rxn in reaction_list:
                self = self.make_lb_ub_constraint(self, rxn, rxn.lower_bound, rxn.upper_bound)
                #reset the reaction bounds
                rxn.lower_bound, rxn.upper_bound = -1e6, 1e6

    def add_lb_ub_constraints(self):
        """
        Makes additional constraints for the reaction lower bounds and upperbounds.
        By adding these constraints the shadow prices of the reaction bounds can be
        calculated and used in sensitivity analysis
        """
        for rxn in self.reactions:
            self  = self.make_lb_ub_constraint(self, rxn, rxn.lower_bound, rxn.upper_bound)
            rxn.lower_bound, rxn.upper_bound = -1e6, 1e6

    @staticmethod
    def make_lb_ub_constraint(m: Optional[Model], rxn: Reaction, lower_bound: float, upper_bound:float):
        """
        Adding variables and constraints for the lower and upperbounds of a reaction to a model.
        When solving the model, shadow prices for the lower and upperbounds will be calculated.
        This allows for the calculation of sensitivity coefficients. The constraints are formulated as follows:
        R_ub : R_fwd-R_rev <= UB
        R_lb : -(R_fwd-R_rev) <= -LB

        Parameters
        ----------
        m: cobra.Model or PAModelpy.PAModel
            model to which the upper and lowerbound constraints and variables should be added
        rxn: cobra.Reaction
            reaction for which an upper and lowerbound constraints should be generated
        lower_bound: float
            Value of the lowerbound
        upper_bound: float
            Value of the upperbound

        Returns
        -------
        m : cobra.Model or PAModelpy.PAModel
            model with additional constraints and variables for the reactions
        """
        #check if the constraints already exists in the model. If it does, change the bounds
        if f'{rxn.id}_ub' in m.constraints.keys():
            m.constraints[f'{rxn.id}_ub'].ub = upper_bound
        if f'{rxn.id}_lb' in m.constraints.keys():
            m.constraints[f'{rxn.id}_lb'].ub = -lower_bound
        else:
            ub_constraint = m.problem.Constraint(Zero, name=f'{rxn.id}_ub', lb=-1e6, ub=upper_bound)
            lb_constraint = m.problem.Constraint(Zero, name=f'{rxn.id}_lb', lb=-1e6, ub=-lower_bound)
            m.add_cons_vars([ub_constraint, lb_constraint])

            # setting up the constraints
            m.constraints[f'{rxn.id}_ub'].set_linear_coefficients({
                rxn.forward_variable: 1,
                rxn.reverse_variable: -1
            })
            m.constraints[f'{rxn.id}_lb'].set_linear_coefficients({
                rxn.forward_variable: -1,
                rxn.reverse_variable: 1
            })

        return m

    @staticmethod
    def make_enzyme_min_max_constraint(m: Optional[Model], enz: Enzyme, lower_bound: float, upper_bound:float):
        """
        Adding variables and constraints for the lower and upperbounds of a reaction to a model.
        When solving the model, shadow prices for the lower and upperbounds will be calculated.
        This allows for the calculation of sensitivity coefficients. The constraints are formulated as follows:
        R_ub : E <= Emax
        R_lb : -E <= -Emin

        Parameters
        ----------
        m: cobra.Model or PAModelpy.PAModel
            model to which the upper and lowerbound constraints and variables should be added
        rxn: PAModelpy.Enzyme
            Enzyme for which minimal and maximal concentration constraints should be generated
        lower_bound: float
            Value of the lowerbound
        upper_bound: float
            Value of the upperbound

        Returns
        -------
        m : cobra.Model or PAModelpy.PAModel
            model with additional constraints and variables for the reactions
        """
        # make separate constraints for forward and reverse enzyme variables
        fwd_var = m.enzyme_variables.get_by_id(enz.id).forward_variable
        rev_var = m.enzyme_variables.get_by_id(enz.id).reverse_variable


        if f'{enz.id}_max' in m.constraints.keys():
            m.constraints[f'{enz.id}_max'].ub = upper_bound
        if f'{enz.id}_min' in m.constraints.keys():
            m.constraints[f'{enz.id}_min'].ub = -lower_bound
        else:
            max_constraint = m.problem.Constraint(Zero, name=f'{enz.id}_max', lb=-1e6, ub=upper_bound)
            min_constraint = m.problem.Constraint(Zero, name=f'{enz.id}_min', lb=-1e6, ub=-lower_bound)
            m.add_cons_vars([max_constraint, min_constraint])

            # setting up the constraints
            m.constraints[f'{enz.id}_max'].set_linear_coefficients({
                fwd_var: 1,
                rev_var:1
            })
            m.constraints[f'{enz.id}_min'].set_linear_coefficients({
                fwd_var: -1,
                rev_var: -1
            })

        #save the constraints in the enzyme object
        enz._constraints = {**enz._constraints, **{f'{enz.id}_max':m.constraints[f'{enz.id}_max'], f'{enz.id}_min':m.constraints[f'{enz.id}_min']}}

        return m,enz

    def determine_sensitivity_coefficients(self):
        obj_value = self.objective.value  # v_z
        if obj_value == 0:
            print('Objective value is 0, thus sensitivity coefficients cannot be calculated')
            return
        mu = self.parse_shadow_prices(self.solver.shadow_prices)
        mu_ub = mu[(mu['direction'] == 'ub')].reset_index()
        mu_lb = mu[(
                mu['direction'] == 'lb')].reset_index()
        mu_ec_max = mu[(mu['direction'] == 'max')].reset_index()
        mu_ec_min = mu[(mu['direction'] == 'min')].reset_index()
        mu_ec_f = mu[(mu['direction'] == 'f')].reset_index()
        mu_ec_b = mu[(mu['direction'] == 'b')].reset_index()

        self.calculate_csc(obj_value, mu, mu_ub, mu_lb, mu_ec_max, mu_ec_min)
        self.calculate_esc(obj_value, mu_ec_f, mu_ec_b)
        # self.validate_sensitivity_coefficients()


    @staticmethod
    def parse_shadow_prices(shadow_prices):
        """
        Parsing the shadow prices to a dataframe which has for each constraint a row and the shadowprices and directions
        in the columns
        :param shadow_prices:
        :return:
        """
        df = pd.DataFrame(pd.Series(shadow_prices), columns=['shadow_prices'])
        # extract only reaction bounds
        df['rxn_id'] = df.index

        splitted_rxn_df = pd.DataFrame(
            [x.rsplit('_', 1) for x in df.rxn_id.tolist()],
            columns=['rxn_id', 'direction']
        )
        df = df.drop('rxn_id', axis=1)
        df = df.reset_index()
        df_long = pd.concat([df, splitted_rxn_df], axis=1)
        # df_long[['rxn_id', 'direction']] = df_long['rxn_id'].str.rsplit('_', 1, expand = True).rename(columns=lambda x: 'col{}'.format(x + 1))
        return df_long

    def calculate_csc(self, obj_value, mu, mu_ub, mu_lb, mu_ec_f, mu_ec_b):
        """
        Calculates the capacity sensitivity coefficient for all inequality constraints in the model.
        The sum of all capacity sensitivity coefficients should equal 1 for growth maximization.

        Definition: constraint_UB*shadowprice/obj_value.

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
        self.capacity_sensitivity_coefficients = pd.DataFrame(columns=['rxn_id', 'enzyme_id', 'constraint', 'coefficient'])
        # add capacity sensitivity coefficients for sectors if they are there
        if self.TOTAL_PROTEIN_CONSTRAINT_ID in self.constraints.keys():
            for sector in self.sectors:
                constraint = 'proteome'
                if isinstance(sector, ActiveEnzymeSector):
                    rxn_id = self.TOTAL_PROTEIN_CONSTRAINT_ID
                    enzyme_id = self.TOTAL_PROTEIN_CONSTRAINT_ID
                    ca_coefficient = self.constraints[enzyme_id].ub * mu[mu['rxn_id'] == self.TOTAL_PROTEIN_CONSTRAINT_ID]['shadow_prices'].iloc[0] / obj_value
            new_row = [rxn_id, enzyme_id, constraint, ca_coefficient]
            # add new_row to dataframe
            self.capacity_sensitivity_coefficients.loc[len(self.capacity_sensitivity_coefficients)] = new_row

            #treat sectors separately if there is not a total protein constraint
        else:
            for sector in self.sectors:
                constraint = 'sector'
                rxn_id = 'R_' + sector.id
                enzyme_id = sector.id
                ca_coefficient = self.constraints[enzyme_id].ub * mu[mu['rxn_id'] == sector.id]['shadow_prices'].iloc[0] / obj_value

                new_row = [rxn_id, enzyme_id, constraint, ca_coefficient]
                # add new_row to dataframe
                self.capacity_sensitivity_coefficients.loc[len(self.capacity_sensitivity_coefficients)] = new_row

        for rxn in self.reactions:
            # LB
            sign=1

            if 'EX_' in rxn.id or self.constraints[f'{rxn.id}_lb'].ub<0:
                sign = -1

            ca_coefficient_LB = -sign*self.constraints[f'{rxn.id}_lb'].ub * mu_lb[mu_lb['rxn_id'] == rxn.id]['shadow_prices'].iloc[0] / obj_value
            # UB
            ca_coefficient_UB = self.constraints[f'{rxn.id}_ub'].ub * mu_ub[mu_ub['rxn_id'] == rxn.id]['shadow_prices'].iloc[0] / obj_value

            new_row_UB = [rxn.id,'', 'flux_ub', ca_coefficient_UB]
            new_row_LB = [rxn.id,'', 'flux_lb', ca_coefficient_LB]
            # add new_row to dataframe
            self.capacity_sensitivity_coefficients.loc[len(self.capacity_sensitivity_coefficients)] = new_row_UB
            self.capacity_sensitivity_coefficients.loc[len(self.capacity_sensitivity_coefficients)] = new_row_LB

        for enzyme in self.enzymes:
            # get reactions associated with this enzyme
            reactions = ','.join(self.get_reactions_with_enzyme_id(enzyme.id))
            # get the right row from the shadow price dataframes
            # min constraints can be skipped as they always equal to 0 (the enzymes cannot have a concentration lower than 0)
            mu_ec_max_row = mu_ec_f[mu_ec_f['index'] == f'{enzyme.id}_max']
            mu_ec_min_row = mu_ec_b[mu_ec_b['index'] == f'{enzyme.id}_min']

            # max Enzyme constraint
            ca_coefficient_EC_max= self.constraints[f'{enzyme.id}_max'].ub * mu_ec_max_row['shadow_prices'].iloc[0] / obj_value
            new_enzyme_row_EC_max =[reactions, enzyme.id, 'enzyme_max', ca_coefficient_EC_max]
            # add new_row to dataframe
            self.capacity_sensitivity_coefficients.loc[len(self.capacity_sensitivity_coefficients)] = new_enzyme_row_EC_max

            # min Enzyme constraint
            ca_coefficient_EC_min = self.constraints[f'{enzyme.id}_min'].ub * mu_ec_min_row['shadow_prices'].iloc[0] / obj_value
            new_enzyme_row_EC_min =[reactions, enzyme.id, 'enzyme_min', ca_coefficient_EC_min]
            # add new_row to dataframe
            self.capacity_sensitivity_coefficients.loc[len(self.capacity_sensitivity_coefficients)] = new_enzyme_row_EC_min

    def calculate_esc(self, obj_value, mu_ec_f, mu_ec_b):
        """
        Calculates enzyme sensitivity coefficients for the enzyme variables using their primal values
        the objective value and shadow prices according to the following relations:

        esc = enzyme_variable.primal * constraint.shadowprice/obj.value

        :param obj_value: float: objective value from the most recent optimal solution
        :param mu_ec_f: pd.DataFrame: shadow prices for max enzyme concentrations (forward variables)
        :param mu_ec_b: pd.DataFrame: shadow prices for the min enzyme concentrations (reverse variables)

        :return: fills the PAModel.enzyme_sensitivity_coefficients dataframe
        """
        self.enzyme_sensitivity_coefficients = pd.DataFrame(columns=['rxn_id', 'enzyme_id', 'coefficient'])

        #calculate enzyme sensitivity coefficient
        for enzyme in self.enzymes:
            #get the reactions associated with the enzyme
            reactions = ','.join(self.get_reactions_with_enzyme_id(enzyme.id))

            #get the right row from the shadow price dataframes
            sp_ec_f = mu_ec_f[mu_ec_f['rxn_id']== f'EC_{enzyme.id}']['shadow_prices'].iloc[0]
            sp_ec_b = mu_ec_b[mu_ec_b['rxn_id'] == f'EC_{enzyme.id}']['shadow_prices'].iloc[0]
            e_fwd = self.enzyme_variables.get_by_id(enzyme.id).forward_variable.primal
            e_rev = self.enzyme_variables.get_by_id(enzyme.id).reverse_variable.primal

            #EC: enzyme constraint
            enzyme_sensitivity_coefficient = (e_fwd * sp_ec_f + e_rev * sp_ec_b)/obj_value
            #add new_row to dataframe
            self.enzyme_sensitivity_coefficients.loc[len(self.enzyme_sensitivity_coefficients)] = [reactions, enzyme.id,
                                                                                                   enzyme_sensitivity_coefficient]

    def calculate_sum_of_enzymes(self):
        """
        Calculates the sum of all enzyme variables for a feasible solution

        :return: sum: float: sum of all enzyme variables in mg/gcdw/h
        """
        if self.solver.status != 'optimal':
            warnings.warn('Cannot calculate sum of enzymes: solution is not optimal')
            return
        sum = 0 #mg/gcdw/h
        for enzyme in self.enzyme_variables:
            #convert from mmol to mg using the same formulation as the coefficients in the protein pool constraint
            sum+= enzyme.concentration *enzyme.molmass *1e-6
        return sum


    def change_total_protein_constraint(self, p_tot):
        """
        Changing the fraction of biomass which is allocated to active proteins
        Parameters
        ----------
        p_tot: float
            new proteome fraction in g_protein/g_cdw
        """
        print(f'Change total condition-dependent protein constraint from {self.p_tot} to {p_tot}')
        tot_prot_constraint = self.constraints[self.TOTAL_PROTEIN_CONSTRAINT_ID]
        protein_availability = tot_prot_constraint.ub
        # correct for the difference between old and new total protein to keep the correction for the protein sections (ptot = Etot - phi_t,0 - phi_ue,0)
        new_protein_fraction = p_tot * 1e3
        for sector in self.sectors:
            if hasattr(sector, 'intercept'):
                new_protein_fraction -= sector.intercept
        self.constraints[self.TOTAL_PROTEIN_CONSTRAINT_ID].ub = new_protein_fraction
        self.p_tot = p_tot
        self.solver.update()

    def change_sector_parameters(self, sector, slope:float, intercept:float, lin_rxn_id:str):
        # input in g/gDW
        print(f'Changing the slope and intercept of the {sector.id}')
        print(f'Changing slope from {sector.slope} to {slope*1e3} mg/gcdw/h')
        print(f'Changing intercept from {sector.intercept} to {intercept*1e3} mg/gcdw')

        prev_intercept = sector.intercept
        #*1e3 to convert g to mg
        sector.slope = slope * 1e3
        sector.intercept = intercept * 1e3
        lin_rxn = self.reactions.get_by_id(lin_rxn_id)

        if self.TOTAL_PROTEIN_CONSTRAINT_ID in self.constraints.keys():
            intercept_diff = sector.intercept - prev_intercept
            #set the intercept
            self.constraints[self.TOTAL_PROTEIN_CONSTRAINT_ID].ub = self.constraints[self.TOTAL_PROTEIN_CONSTRAINT_ID].ub - intercept_diff
            #reset the slope
            self.constraints[self.TOTAL_PROTEIN_CONSTRAINT_ID].set_linear_coefficients({
                lin_rxn.forward_variable: -slope,
                lin_rxn.reverse_variable: slope
            })

        else:
            var = self.variables['R_' + sector.id]
            # update the constraint
            self.constraints[sector.id].set_linear_coefficients({
                var: 1,
                lin_rxn.forward_variable: -slope, #/ (sector.mol_mass[0] * 1e-6),
                lin_rxn.reverse_variable: slope #/ (sector.mol_mass[0] * 1e-6)
            })
            # update the sector object
            sector.variables = [var]
            sector.constraints = [self.constraints[sector.id]]

    def change_reaction_bounds(self, rxn_id:str, lower_bound: float = None, upper_bound: float = None):
        """
        Change the reaction bounds. If there should be a sensitivity analysis, the bounds of the upper and lower bound
        constraints are adjusted
        Parameters
        ----------
        rxn_id: str
            string of reaction id to change
        lower_bound: float, optional
            value of the lower bound
        upper_bound: float, optional
            value of the upper bound
        """
        if rxn_id not in self.reactions:
            warnings.warn(f'Reaction {rxn_id} does not exist in the model. Cannot change the upper and lowerbound.')
            return
        rxn = self.reactions.get_by_id(rxn_id)
        #make sure the order of setting is right to prevent errors
        if lower_bound is not None and lower_bound >= rxn.upper_bound:
            if upper_bound is not None: self.change_reaction_ub(rxn_id, upper_bound)
            self.change_reaction_lb(rxn_id,lower_bound)

        elif upper_bound is not None:
            if lower_bound is not None: self.change_reaction_lb(rxn_id, lower_bound)
            self.change_reaction_ub(rxn_id,upper_bound)

    def change_reaction_ub(self, rxn_id:str, upper_bound: float = None):
        if self.sensitivity:
            self.constraints[rxn_id + '_ub'].ub = upper_bound
        else:
            self.reactions.get_by_id(rxn_id).upper_bound = upper_bound

    def change_reaction_lb(self,rxn_id:str, lower_bound: float = None):
        if self.sensitivity:
            self.constraints[rxn_id + '_lb'].ub = -lower_bound
        else:
            self.reactions.get_by_id(rxn_id).lower_bound = lower_bound

    def change_enzyme_bounds(self, enzyme_id: str, lower_bound: float = None, upper_bound: float = None):
        """
                Change the enzyme bounds. If the model should be primed for performing a sensitivity analysis,
                the upper bound of the min and max enzyme concentration constraints are adjusted
                Parameters
                ----------
                enzyme_id: str
                    string of enzyme id to change
                lower_bound: float, optional
                    value of the minimal enzyme concentration
                upper_bound: float, optional
                    value of the maximal enzyme concentration
                """
        if enzyme_id not in self.enzyme_variables:
            warnings.warn(f'Enzyme {enzyme_id} does not exist in the model. Cannot change the minimal and maximal concentrations.')
            return
        enzyme = self.enzyme_variables.get_by_id(enzyme_id)
        # make sure the order of setting is right to prevent errors
        if lower_bound is not None and lower_bound > enzyme.upper_bound:
            if upper_bound is not None: self.change_enzyme_max(enzyme_id, upper_bound)
            self.change_enzyme_min(enzyme_id, lower_bound)

        elif upper_bound is not None:
            if lower_bound is not None: self.change_enzyme_min(enzyme_id, lower_bound)
            self.change_enzyme_max(enzyme_id, upper_bound)

    def change_enzyme_max(self, enzyme_id: str, upper_bound: float = None):
        if self.sensitivity:
            self.constraints[enzyme_id + '_max'].ub = upper_bound
        else:
            self.enzyme_variables.get_by_id(enzyme_id).upper_bound = upper_bound

    def change_enzyme_min(self, enzyme_id: str, lower_bound: float = None):
        if self.sensitivity:
            self.constraints[enzyme_id + '_min'].ub = -lower_bound
        else:
            self.enzyme_variables.get_by_id(enzyme_id).lower_bound = lower_bound

    def get_enzymes_with_reaction_id(self, rxn_id:str):
        """
        Returns Enzyme objects associated with the reaction id through CatalyticEvent objects
        :param rxn_id: str
            reaction identifier

        Returns
        -------
            DictList of Enzyme objects associated with the reaction
        """
        catalytic_event_id = 'CE_' + rxn_id
        catalytic_event = self.catalytic_events.get_by_id(catalytic_event_id)
        enzymes = catalytic_event.enzymes
        return enzymes

    def get_reactions_with_enzyme_id(self, enz_id:str):
        """
        Returns Enzyme objects associated with the reaction id through CatalyticEvent objects
        :param enz_id: str
            enzyme identifier (EC number)

        Returns
        -------
            List of reaction ids associated with the enzyme
        """
        enzyme = self.enzymes.get_by_id(enz_id)
        rxn_ids = list(enzyme.rxn2kcat.keys())
        return rxn_ids

    def change_kcat_value(self, enzyme_id:str, kcats:dict):
        """
        Change the turnover number of the enzyme for a specific reaction
        :param enzyme_id: Enzyme identifier
        :param kcats: dict with reaction id, kcat key, value pairs. kcat is a dict with direction, kcat value
        key, value pairs {'f': fwd_value, 'b': bckwrd_value}
        """
        if self.enzymes.has_id(enzyme_id):
            enzyme = self.enzymes.get_by_id(enzyme_id)
            enzyme.change_kcat_values(kcats)
            #also change the active enzyme sector
            active_enzyme = self.sectors.get_by_id('ActiveEnzymeSector')
            for rxn, kcat_f_b in kcats.items():
                active_enzyme.rxn2protein[rxn][enzyme_id] = kcat_f_b
        else:
            warnings.warn(f'The enzyme {enzyme_id} does not exist in the model. The kcat can thus not be changed.')

    def remove_enzymes(self,
                       enzymes: Union[str, Enzyme, List[Union[str, Enzyme]]]
                       )-> None:
        """Remove enzymes from the model. Adapted from the cobra.core.remove_reactions() function.
        If the reaction is not catalyzed by any enzyme anymore, the reaction ub will become 0

                The change is reverted upon exit when using the model as a context.

                Parameters
                ----------
                enzymes : list or reaction or str
                    A list with reactions (`cobra.Reaction`), or their id's, to remove.
                    Enzyme will be placed in a list. Str will be placed in a list and used to
                    find the reaction in the model.
                """
        if isinstance(enzymes, str) or hasattr(enzymes, "id"):
            warnings.warn("need to pass in a list")
            enzymes = [enzymes]

        context = get_context(self)

        for enzyme in enzymes:
            # Make sure the reaction is in the model
            try:
                enzyme = self.enzymes[self.enzymes.index(enzyme)]
            except ValueError:
                warnings.warn(f"{enzyme.id} not in {self}")
            for constraint in enzyme._constraints.values():
                self.remove_cons_vars([constraint])
                if constraint.name in self.enzyme_constraints.keys():
                    del self.enzyme_constraints[constraint.name]
            self.enzymes.remove(enzyme)
            enzyme._model = None

            #remove variable related to the enzyme
            self.enzyme_variables.remove(enzyme.enzyme_variable)
            forward = enzyme.enzyme_variable.forward_variable
            reverse = enzyme.enzyme_variable.reverse_variable

            self.remove_cons_vars([forward, reverse])

            #remove enzyme from catalytic event
            for ce in enzyme.catalytic_events:
                ce.remove_enzymes([enzyme])

            # remove reference to the enzymes in all groups
            associated_groups = self.get_associated_groups(enzyme)
            for group in associated_groups:
                group.remove_members(enzyme)

    def remove_reactions(
        self,
        reactions: Union[str, Reaction, List[Union[str, Reaction]]],
        remove_orphans: bool = False,
    ) -> None:

        """Remove reactions from the model. Inherited from the cobrapy.core.remove_reactions() function.

                The change is reverted upon exit when using the model as a context. Also removes associated
                CatalyticEvents if they exist.

                Parameters
                ----------
                reactions : list or reaction or str
                    A list with reactions (`cobra.Reaction`), or their id's, to remove.
                    Reaction will be placed in a list. Str will be placed in a list and used to
                    find the reaction in the model.
                remove_orphans : bool, optional
                    Remove orphaned genes and metabolites from the model as
                    well (default False).
                """
        if not hasattr(reactions, "__iter__"):
            reactions = [reactions]
        super().remove_reactions(reactions, remove_orphans)

        #remove the catalytic events if they exist for the reactions
        for rxn in reactions:
            # catalytic_event = self.catalytic_events.get_by_id('CE_' + rxn.id)
            # self.catalytic_events.remove(catalytic_event)
            #
            # # removing catalytic event from the enzymes
            # for enzyme in catalytic_event.enzymes:
            #     enzyme.remove_catalytic_event(catalytic_event)
            #
            # for enzyme_var in catalytic_event.enzyme_variables:
            #     enzyme_var.remove_catalytic_event(catalytic_event)
            #
            #     # removing orphaned enzymes
            #     if remove_orphans and len(enzyme.catalytic_events) == 0:
            #         print(f'Enzyme {enzyme.id} is orphaned. This enzyme will be removed from the model')
            #         self.remove_cons_vars(enzyme._constraints)
            #         enzyme._model = None
            #         self.enzymes.remove(enzyme)
            #         self.remove_cons_vars(enzyme._constraints)

            try:
                catalytic_event = self.catalytic_events.get_by_id('CE_' + rxn.id)
                self.catalytic_events.remove(catalytic_event)

                # removing catalytic event from the enzymes
                for enzyme in catalytic_event.enzymes:
                    enzyme.remove_catalytic_event(catalytic_event)

                for enzyme_var in catalytic_event.enzyme_variables:
                    enzyme_var.remove_catalytic_event(catalytic_event)

                    #removing orphaned enzymes
                    if remove_orphans and len(enzyme.catalytic_events) == 0:
                        print(f'Enzyme {enzyme.id} is orphaned. This enzyme will be removed from the model')
                        self.remove_cons_vars(enzyme._constraints)
                        enzyme._model = None
                        self.enzymes.remove(enzyme)
                        self.remove_cons_vars(enzyme._constraints)

                if self.sensitivity:
                    lb_constraint = self.constraints[rxn.id + '_lb']
                    ub_constraint = self.constraints[rxn.id + '_ub']
                    self.remove_cons_vars([lb_constraint, ub_constraint])

            except:
                continue

    def remove_catalytic_events(self,
                                catalytic_events:Union[str, CatalyticEvent, List[Union[str,CatalyticEvent]]],
                                remove_orphans: bool = False
                                )-> None:
        """Remove catalytic events from the model.


        Parameters
        ----------
        reactions : list or reaction or str
                    A list with reactions (`cobra.Reaction`), or their id's, to remove.
                    Reaction will be placed in a list. Str will be placed in a list and used to
                    find the reaction in the model.
        remove_orphans : bool, optional
                    Remove orphaned genes and metabolites from the model as
                    well (default False).
                """
        for catalytic_event in catalytic_events:
            self.catalytic_events.remove(catalytic_event)

            # removing catalytic event from the enzymes
            for enzyme in catalytic_event.enzymes:
                enzyme.remove_catalytic_event(catalytic_event)

                # removing orphaned enzymes
                if remove_orphans and len(enzyme.catalytic_events) == 0:
                    print(f'Enzyme {enzyme.id} is orphaned. This enzyme will be removed from the model')
                    self.remove_cons_vars(enzyme._constraints)
                    enzyme._model = None
                    self.enzymes.remove(enzyme)
                    self.remove_cons_vars(enzyme._constraints)

    def remove_sectors(self,
                       sectors: Union[str, Sector, ActiveEnzymeSector, List[Union[str, Sector, ActiveEnzymeSector]]]
                       )-> None:
        """Remove sections from the model.

            Also removes associated CatalyticEvents if they exist.

                Parameters
                ----------
                secotrs : list or sector or str
                    A list with sector (`PAModelpy.Sector` or `PAModelpy.ActiveEnzymeSector`), or their id's, to remove.
                    A single sector will be placed in a list. Str will be placed in a list and used to
                    find the secotr in the model.
                """
        if isinstance(sectors, str) or hasattr(sectors, "id"):
            sectors = [sectors]

        for sector in sectors:
            print(f'Removing the following protein sector: {sector.id}\n')
            # remove the connection to the model
            sector._model = None

            # Make sure the sector is in the model
            try:
                sector = self.sectors[self.sectors.index(sector)]
            except ValueError:
                warnings.warn(f"{sector.id} not in {self}")

            # remove the sector from the total protein if it is there
            if self.TOTAL_PROTEIN_CONSTRAINT_ID in self.constraints.keys():
                # remove parts of constraint corresponding to the enzyme sector from the total_protein_constraint
                # 1. add the intercept value from the sum of protein (Total_protein == Etot)
                tpc_ub = self.constraints[self.TOTAL_PROTEIN_CONSTRAINT_ID].ub
                self.constraints[self.TOTAL_PROTEIN_CONSTRAINT_ID].ub = tpc_ub + sector.intercept

                # 2. remove link between flux and enzyme concentration
                # link enzyme concentration in the sector to the total enzyme concentration
                lin_rxn = self.reactions.get_by_id(sector.id_list[0])
                self.constraints[self.TOTAL_PROTEIN_CONSTRAINT_ID].set_linear_coefficients({
                    lin_rxn.forward_variable: 0,
                    lin_rxn.reverse_variable: 0
                })

            else:

                #remove the associated constraints
                for constraint in sector.constraints:
                    if isinstance(constraint, Enzyme):
                        self.remove_enzymes([constraint])
                    #check if constraint is in the solver
                    if constraint in self.constraints.values():
                        self.remove_cons_vars([constraint])
                        self.solver.update()

                #remove the associated variables
                for variable in sector.variables:
                    if isinstance(variable, CatalyticEvent):
                        self.remove_catalytic_events([variable])
                    else:
                        self.remove_cons_vars([variable])
                        self.solver.update()
                        # remove reference to the sector variables in all groups
                        associated_groups = self.get_associated_groups(variable)
                        for group in associated_groups:
                            group.remove_members(variable)

                # remove sector constraint from model easy lookup:
                del self.sector_constraints[sector.id]

            #remove the sector and its connection to the model
            self.sectors.remove(sector)
            sector._model = None

    def test(self, glc_flux: Union[int, float]=10):
        """
        Test the proteome allocation model
        :param glc_flux: glc flux which limits the growth rate (mmol_glc/g_cdw/h, default = 10)
        """
        self.set_glc_uptake_bounds(glc_flux)
        self.optimize()
        if self.solver.status == 'optimal':
            print('Protein allocation model test case was successful.\n')
            print('Optimal objective function value: ',self.objective.value,'\n')
            print('Predicted glucose uptake rate: ',self.reactions.get_by_id(self.GLUCOSE_EXCHANGE_RXNID).flux)

    def set_glc_uptake_bounds(self, bound: Union[int, float]):
        #check if the model is reversible
        if self.GLUCOSE_EXCHANGE_RXNID[-1] == 'b':
            self.change_reaction_bounds(self.GLUCOSE_EXCHANGE_RXNID,bound, bound)
        else:
            self.change_reaction_bounds(self.GLUCOSE_EXCHANGE_RXNID,-bound, -bound)


    def pfba(self,
             fraction_of_optimum: float = 1.0,
             proteins:bool = False,
             reactions:bool = True,
             exclude:List['str'] =[],
             objective: Union[Dict, "Objective", None] = None):
        """Perform pFBA (parsimonious Enzyme Usage Flux Balance Analysis) with a custom objective including:
            all reactions, all proteins, all proteins and all reactions.
            pFBA [1] adds the minimization of all fluxes the objective of the
            model. This approach is motivated by the idea that high fluxes have a
            higher enzyme turn-over and that since producing enzymes is costly,
            the cell will try to minimize overall flux while still maximizing the
            original objective function, e.g. the growth rate.
            Parameters
            ----------
            fraction_of_optimum : float, optional
                The fraction of optimum which must be maintained. The original
                objective reaction is constrained to be greater than maximal value
                times the `fraction_of_optimum` (default 1.0).
            objective : dict or cobra.Model.objective, optional
                A desired objective to use during optimization in addition to the
                pFBA objective. Dictionaries (reaction as key, coefficient as value)
                can be used for linear objectives (default None).
            proteins: bool, optional
                Determines whether to include enzyme variables in the pFBA objective
            reactions: bool, optional
                Determines whether to include reaction variables in the pFBA objective
            exclude: list of reaction ids, optional
                Reactions to exclude fom the minimalization objective


            References
            ----------
            .. [1] Lewis, N. E., Hixson, K. K., Conrad, T. M., Lerman, J. A.,
               Charusanti, P., Polpitiya, A. D., Palsson, B. O. (2010). Omic data
               from evolved E. coli are consistent with computed optimal growth from
               genome-scale models. Molecular Systems Biology, 6,
               390. doi:10.1038/msb.2010.47
            """
        #get the variables which should be included in the objective function
        variables = list()
        if reactions:
            for rxn in self.reactions:
                if rxn not in exclude and rxn not in self.enzymes: variables += [rxn]
        if proteins:
            for enzyme_var in self.enzyme_variables:
                if enzyme_var.id not in exclude: variables += [enzyme_var.forward_variable, enzyme_var.reverse_variable]

        #set custom objective
        if objective is not None:
            self.objective = objective

        #set previous objective as constraint
        if self.solver.objective.name == "_pfba_objective":
            raise ValueError("The model already has a pFBA objective.")
        cobra.util.solver.fix_objective_as_constraint(self, fraction=fraction_of_optimum)

        #add new pFBA objective
        self.objective = self.problem.Objective(
            Zero, direction="min", sloppy=True, name="_pfba_objective"
        )
        self.objective.set_linear_coefficients({v: 1.0 for v in variables})

        #run pFBA
        self.optimize()

    def reset_objective(self):
        """
        Reseting the objective to the standard biomass maximization objective after pFBA
        """
        biomass_rxn = self.reactions.get_by_id(self.BIOMASS_REACTION)
        #reset the biomass reaction to standard UB and LB
        biomass_rxn.lower_bound = 0
        biomass_rxn.upper_bound= 1e3

        #reseting the objective
        self.objective = {biomass_rxn: 1.0}
        self.objective.direction = 'max'

    def optimize(
        self, objective_sense: Optional[str] = None, raise_error: bool = False
    ) -> "Solution":
        """ Optimize the model using flux balance analysis.
            Inherits from the cobra.Model.optimize() function and performs a sensitivity analysis after optimization if
            this is desired (by setting the PAModel.sensitivity attribute to True)
               Parameters
               ----------
               objective_sense : {None, 'maximize' 'minimize'}, optional
                   Whether fluxes should be maximized or minimized. In case of None,
                   the previous direction is used (default None).
               raise_error : bool
                   If true, raise an OptimizationError if solver status is not
                    optimal (default False).
               Returns
               -------
               Solution
               Notes
               -----
               Only the most commonly used parameters are presented here.  Additional
               parameters for cobra.solver may be available and specified with the
               appropriate keyword argument.
               """
        solution = super().optimize(objective_sense,raise_error)
        if self.sensitivity and self.solver.status == 'optimal':
            self.determine_sensitivity_coefficients()
        return solution

    def copy(self) -> 'PAModel':
        """Provide a partial 'deepcopy' of the Model.

                Adjusted from cobra.Model.copy()
               All the Metabolite, Gene, Reaction, Enzyme, EnzymeVariable, Sector and CatalyticEvent
                objects are created anew but in a faster fashion than deepcopy.

               :return: PAModelpy.PAModel: new model copy
               """

        do_not_copy_by_ref = {
            "metabolites",
            "reactions",
            "genes",
            "enzymes",
            "enzyme_variables",
            "enzyme_constraints",
            "sectors",
            "catalytic_events"
            "notes",
            "annotation",
            "groups",
        }

        # setting up the dict for the initialization of the model copy
        model_init_attr = {
            'm_model',
            'name',
            'p_tot',
            'sensitivity',
            'sectors',
            'configuration'
        }

        new_dict = {}
        # for attr in self.__dict__:
        #     if attr not in do_not_copy_by_ref and attr in model_init_attr:
        #         if attr == 'm_model': new_attr = 'id_or_model'
        #         else: new_attr = attr
        #         new_dict[new_attr] = self.__dict__[attr]
        new_dict = self.find_init_args(self)
        new_dict['id_or_model'] = self.m_model
        # initialize new model
        new = self.__class__(**new_dict)
        # also adjust the constants
        for attr in self.__dict__:
            if attr not in do_not_copy_by_ref and not attr in model_init_attr:
                new.__dict__[attr] = self.__dict__[attr]

        new.notes = deepcopy(self.notes)
        new.annotation = deepcopy(self.annotation)

        new.metabolites = DictList()
        do_not_copy_by_ref = {"_reaction", "_model"}
        for metabolite in self.metabolites:
            new_met = metabolite.__class__()
            for attr, value in metabolite.__dict__.items():
                if attr not in do_not_copy_by_ref:
                    new_met.__dict__[attr] = copy(value) if attr == "formula" else value
            new_met._model = new
            new.metabolites.append(new_met)

        new.genes = DictList()
        for gene in self.genes:
            new_gene = gene.__class__(None)
            for attr, value in gene.__dict__.items():
                if attr not in do_not_copy_by_ref:
                    new_gene.__dict__[attr] = (
                        copy(value) if attr == "formula" else value
                    )
            new_gene._model = new
            new.genes.append(new_gene)

        new.reactions = DictList()
        do_not_copy_by_ref = {"_model", "_metabolites", "_genes"}
        for reaction in self.reactions:
            new_reaction = reaction.__class__()
            for attr, value in reaction.__dict__.items():
                if attr not in do_not_copy_by_ref:
                    new_reaction.__dict__[attr] = copy(value)
            new_reaction._model = new
            new.reactions.append(new_reaction)
            # update awareness
            for metabolite, stoic in reaction._metabolites.items():
                new_met = new.metabolites.get_by_id(metabolite.id)
                new_reaction._metabolites[new_met] = stoic
                new_met._reaction.add(new_reaction)
            new_reaction.update_genes_from_gpr()

        #########################################
        # adding protein information
        #########################################
        new.change_total_protein_constraint(self.p_tot)
        new.enzymes = DictList()
        new.catalytic_events = DictList()
        new.enzyme_variables = DictList()
        new.enzyme_constraints = dict()

        # if sectors are added, everything else will be added automatically
        new.sectors = DictList()
        do_not_copy_by_ref = {"model", "variables", "_constraints"}
        for sector in self.sectors:
            d = sector.__dict__
            init_args = self.find_init_args(sector)
            new_sector = sector.__class__(**init_args)
            for attr, value in sector.__dict__.items():
                if attr not in do_not_copy_by_ref:
                    new_sector.__dict__[attr] = copy(value)
            new_sector.model = new
            new.sectors += [sector]

        # if some manual enzymes are added, copy them separately
        do_not_copy_by_ref = {"_model", "enzyme_variables", "catalytic_events", "_constraints"}
        for enzyme in self.enzymes:
            if enzyme not in new.enzymes:  # TODO
                init_args = self.find_init_args(enzyme)
                new_enzyme = enzyme.__class__(**init_args)
                for attr, value in enzyme.__dict__.items():
                    if attr not in do_not_copy_by_ref:
                        new_enzyme.__dict__[attr] = copy(value)
                new_enzyme._model = new
                new.enzymes += [new_enzyme]
                new.enzyme_variables += [new_enzyme.enzyme_variable]
                # new.add_enzymes([enzyme])

        new.groups = DictList()
        do_not_copy_by_ref = {"_model", "_members"}
        # Groups can be members of other groups. We initialize them first and
        # then update their members.
        for group in self.groups:
            new_group: cobra.core.Group = group.__class__(group.id)
            for attr, value in group.__dict__.items():
                if attr not in do_not_copy_by_ref:
                    new_group.__dict__[attr] = copy(value)
            new_group._model = new
            new.groups.append(new_group)
        for group in self.groups:
            new_group = new.groups.get_by_id(group.id)
            # update awareness, as in the reaction copies
            new_objects = []
            for member in group.members:
                if isinstance(member, Metabolite):
                    new_object = new.metabolites.get_by_id(member.id)
                elif isinstance(member, Reaction):
                    new_object = new.reactions.get_by_id(member.id)
                elif isinstance(member, cobra.core.Gene):
                    new_object = new.genes.get_by_id(member.id)
                elif isinstance(member, cobra.core.Group):
                    new_object = new.groups.get_by_id(member.id)
                else:
                    raise TypeError(
                        f"The group member {member!r} is unexpectedly not a "
                        f"metabolite, reaction, gene, nor another group."
                    )
                new_objects.append(new_object)
            new_group.add_members(new_objects)

        try:
            new._solver = deepcopy(self.solver)
            # Cplex has an issue with deep copies
        except Exception:  # pragma: no cover
            new._solver = copy(self.solver)  # pragma: no cover

        # it doesn't make sense to retain the context of a copied model so
        # assign a new empty context
        new._contexts = []

        return new

    def find_init_args(self, object):
        init_args = {}
        for param, default in inspect.signature(object.__init__).parameters.items():
            if param != 'self' and default.default == inspect.Parameter.empty:
                init_args[param] = getattr(object, param)
        return init_args