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
import sys
import numpy as np
import pandas as pd

sys.path.append('../Package/')
import EnzymeSectors
from EnzymeSectors import ActiveEnzymeSector, TransEnzymeSector, UnusedEnzymeSector, CustomSector, Sector
from Variables import Total_protein_variable, Constant, CatalyticEvent, EnzymeVariable
from Constraints import Constraint
from Enzyme import Enzyme
from PAMValidator import PAMValidator
import configuration



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
                    This sensitivity analysis will indicate to which extend individual constraints contribute to the
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
    """Constants"""
    TOTAL_PROTEIN_CONSTRAINT_ID = configuration.TOTAL_PROTEIN_CONSTRAINT_ID
    P_TOT_DEFAULT = configuration.P_TOT_DEFAULT #g_protein/g_cdw
    CO2_EXHANGE_RXNID = configuration.CO2_EXHANGE_RXNID
    GLUCOSE_EXCHANGE_RXNID = configuration.GLUCOSE_EXCHANGE_RXNID
    BIOMASS_REACTION = configuration.BIOMASS_REACTION

    def __init__(self, id_or_model: Union[str, "Model", None] = None,
                 name: Optional[str] = None,
                 p_tot: Optional[float] = P_TOT_DEFAULT,
                 sensitivity: bool = True,
                 active_sector: Optional[ActiveEnzymeSector]=None,
                 translational_sector: Optional[TransEnzymeSector]=None,
                 unused_sector: Optional[UnusedEnzymeSector]=None,
                 custom_sectors: Union[List, CustomSector] =None):

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
        self.control_coefficients = pd.DataFrame() #dataframe to store the result of the sensitivity analysis (control coefficients for each constraint). The sensitivity coefficients are splitted on LB, UB and the different sectors
        self.allocation_coefficients = pd.DataFrame() #dataframe to store the result of the sensitivity analysis (allocation coefficients for each constraint). The sensitivity coefficients are splitted on LB, UB and the different sectors


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

            #add the enzyme to the interface between the enzyme constraints, enzyme variables and the reaction
            for catalytic_event in enzyme.catalytic_events:
                if catalytic_event not in self.catalytic_events:
                    # make catalytic event aware of the model.
                    # The catalytic event will be added configured in the model.setter magic function
                    catalytic_event.model = self
                    self.catalytic_events += [catalytic_event]
                    if context:
                        context(partial(setattr, catalytic_event, "_model", None))
                else:
                    #add enzyme to the catalytic event
                    catalytic_event_model = self.catalytic_events.get_by_id(catalytic_event.id)

                    catalytic_event_model.add_enzymes({enzyme: enzyme.rxn2kcat[catalytic_event.rxn_id]})

                    #replace the catalytic event with the already existing event from the model
                    enzyme.catalytic_events._replace_on_id(catalytic_event_model)
        self.solver.update()

        # connect the enzyme to total protein constraint
        for enzyme in pruned:
            enzyme_variable = self.enzyme_variables.get_by_id(enzyme.id)
            self.constraints[self.TOTAL_PROTEIN_CONSTRAINT_ID].set_linear_coefficients(
                {
                enzyme_variable.forward_variable: enzyme.molmass * 1e-6,
                enzyme_variable.reverse_variable: enzyme.molmass * 1e-6,
                }
            )
            self.tpc+=1

        # add enzyme constraints to the model and link to enzyme and reaction variables
        for enzyme in pruned:
            #get the enzyme variable
            enzyme_var_model = self.enzyme_variables.get_by_id(enzyme.id)
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
                    coeff = kcat*3600*1e-6
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
        self.sectors += sectors
        #making sure the right type is checked if code is run from jupyter notebook or from the terminal
        if not PAMValidator.check_kernel_type():
            # import Package
            # type = Package.EnzymeSectors.ActiveEnzymeSector
            type = ActiveEnzymeSector
        else:
            type = ActiveEnzymeSector
        for sector in sectors:
            #different method to add the active enzyme_sector
            if isinstance(sector, type):
                self = sector.add(self)
            else:
                self.add_sector(sector)

    def add_sector(self, sector):
        print(f'Add the following protein sector: {sector.id}\n')
        #make sector aware of the model
        sector.model = self

        lin_rxn = self.reactions.get_by_id(sector.id_list[0])

        totprot, tpc_metabolite = sector.get_tpc_metabolite(self)

        # constraint corresponding to the enzyme sector
        # constraint is equal to the intercept, like in the mathematical definition
        constraint = self.problem.Constraint(Zero, name = sector.id, lb = sector.intercept, ub = sector.intercept)
        self.add_sector_constraints([constraint])
        sector.constraints += [self.sector_constraints[sector.id]]

        # Protein concentration allocated to enzymes required for translation
        var = self.problem.Variable('R_' + sector.id)
        self.add_cons_vars(var)

        # link flux to enzyme concentration
        self.constraints[sector.id].set_linear_coefficients({
            var: 1,
            lin_rxn.forward_variable: -sector.slope,
            lin_rxn.reverse_variable: sector.slope
        })

        if totprot:
            # linear relation between the amount of translational enzymes and growth rate
            # *1000 to convert units from g/g_cdw to mg/g_cdw
            # /MW to convert to mg to mmol (which is used in the model)
            self.constraints[sector.TOTAL_PROTEIN_CONSTRAINT_ID].set_linear_coefficients({
                var: sector.mol_mass[0] * 1e-6
            })

        # save reaction for easy removal of the sector
        sector.variables += [var]
        sector.constraints.append(self.constraints[sector.id])

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
        Etot: sum(E) + E_translprot + E_unusedprot - ptot == 0

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
        # Create the pseudo 'exchange reaction' which will determine how much of the total protein is consumed
        # tpc_variable = self.problem.Variable(name = 'R_'+self.TOTAL_PROTEIN_CONSTRAINT_ID, lb =0, ub= self.p_tot*1000)

        # Create the pseudometabolite associated with the constraint, this metabolite will be 'produced' in the enzyme 'reactions'
        tpc_constraint= self.problem.Constraint(Zero, name = self.TOTAL_PROTEIN_CONSTRAINT_ID, lb = 0, ub = self.p_tot*1000)
        self.sector_constraints[self.TOTAL_PROTEIN_CONSTRAINT_ID] = tpc_constraint
        # self.add_cons_vars([tpc_constraint, tpc_variable])
        self.add_cons_vars([tpc_constraint])


        #connect the variable and constraint
        # self.constraints[self.TOTAL_PROTEIN_CONSTRAINT_ID].set_linear_coefficients({
        #     tpc_variable: -1
        # })

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
        m: cobra.Model or PAMpy.PAModel
            model to which the upper and lowerbound constraints and variables should be added
        rxn: cobra.Reaction
            reaction for which an upper and lowerbound constraints should be generated
        lower_bound: float
            Value of the lowerbound
        upper_bound: float
            Value of the upperbound

        Returns
        -------
        m : cobra.Model or PAMpy.PAModel
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

    def determine_control_coefficients(self):
        obj_value = self.objective.value  # v_z
        mu = self.parse_shadow_prices(self.solver.shadow_prices)
        mu_ub = mu[(mu['direction'] == 'ub')].reset_index()
        mu_lb = mu[(mu['direction'] == 'lb')].reset_index()
        mu_ec_f = mu[(mu['direction'] == 'f')].reset_index()
        mu_ec_b = mu[(mu['direction'] == 'b')].reset_index()


        self.calculate_control_coefficients(obj_value, mu, mu_ub, mu_lb, mu_ec_f, mu_ec_b)
        self.calculate_allocation_coefficients(obj_value, mu, mu_ub, mu_lb, mu_ec_f, mu_ec_b)


    @staticmethod
    def parse_shadow_prices(shadow_prices):
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

    def calculate_control_coefficients(self, obj_value, mu, mu_ub, mu_lb, mu_ec_f, mu_ec_b):
        self.control_coefficients = pd.DataFrame(columns=['rxn_id', 'enzyme_id', 'constraint', 'coefficient'])
        # add control coefficients for sectors if they are there
        for sector in self.sectors:
            constraint = 'sector'
            if isinstance(sector, ActiveEnzymeSector):
                rxn_id = self.TOTAL_PROTEIN_CONSTRAINT_ID
                enzyme_id = self.TOTAL_PROTEIN_CONSTRAINT_ID
                control_coefficient = self.constraints[enzyme_id].ub * mu[mu['rxn_id'] == self.TOTAL_PROTEIN_CONSTRAINT_ID]['shadow_prices'].iloc[0] / obj_value
            else:
                rxn_id = 'R_' + sector.id
                enzyme_id = sector.id
                control_coefficient = self.constraints[enzyme_id].ub * mu[mu['rxn_id'] == sector.id]['shadow_prices'].iloc[0] / obj_value

            new_row = [rxn_id, enzyme_id, constraint, control_coefficient]
            # add new_row to dataframe
            self.control_coefficients.loc[len(self.control_coefficients)] = new_row

        for rxn in self.reactions:
            new_row = list()
            # LB
            control_coefficient_LB = -self.constraints[f'{rxn.id}_lb'].ub * mu_lb[mu_lb['rxn_id'] == rxn.id]['shadow_prices'].iloc[0] / obj_value
            # UB
            control_coefficient_UB = self.constraints[f'{rxn.id}_ub'].ub * mu_ub[mu_ub['rxn_id'] == rxn.id]['shadow_prices'].iloc[0] / obj_value

            new_row_UB = [rxn.id,'', 'UB', control_coefficient_UB]
            new_row_LB = [rxn.id,'', 'LB', control_coefficient_LB]
            # add new_row to dataframe
            self.control_coefficients.loc[len(self.control_coefficients)] = new_row_UB
            self.control_coefficients.loc[len(self.control_coefficients)] = new_row_LB

            # get all enzymes related to the reaction
            if f'CE_{rxn.id}' in self.catalytic_events:
                for enzyme in self.catalytic_events.get_by_id(f'CE_{rxn.id}').enzymes:

                    # get the right row from the shadow price dataframes
                    mu_ec_f_row = mu_ec_f[mu_ec_f['rxn_id'] == f'EC_{enzyme.id}']
                    mu_ec_b_row = mu_ec_b[mu_ec_b['rxn_id'] == f'EC_{enzyme.id}']

                    # EC_f: forward enzyme constraint
                    control_coefficient_EC_f = self.enzyme_constraints[f'EC_{enzyme.id}_f'].ub * mu_ec_f_row['shadow_prices'].iloc[0] / obj_value

                    new_enzyme_row_EC_f =[rxn.id, enzyme.id, 'EC_f', control_coefficient_EC_f]

                    # EC_b: reverse enzyme constraint
                    control_coefficient_EC_b = self.enzyme_constraints[f'EC_{enzyme.id}_b'].ub * mu_ec_b_row['shadow_prices'].iloc[0] / obj_value

                    new_enzyme_row_EC_b =[rxn.id, enzyme.id, 'EC_b', control_coefficient_EC_b]

                    # add new_row to dataframe
                    self.control_coefficients.loc[len(self.control_coefficients)] = new_enzyme_row_EC_f
                    self.control_coefficients.loc[len(self.control_coefficients)] = new_enzyme_row_EC_b

    def calculate_allocation_coefficients(self,obj_value, mu, mu_ub, mu_lb, mu_ec_f, mu_ec_b):
        self.allocation_coefficients = pd.DataFrame(columns=['rxn_id', 'enzyme_id', 'constraint', 'coefficient'])
        #add allocation coefficients for sectors if they are there
        for sector in self.sectors:
            if not isinstance(sector, ActiveEnzymeSector):
                sector_0 = self.constraints[sector.id].ub
                sector_primal = self.variables['R_' + sector.id].primal
                allocation_coefficient = (sector_0-sector_primal) * mu[mu['rxn_id'] == sector.id]['shadow_prices'].iloc[0] /obj_value
                # add new_row to dataframe
                self.allocation_coefficients.loc[len(self.allocation_coefficients)] = ['R_'+sector.id, sector.id, 'sector', allocation_coefficient]

        for rxn in self.reactions:
            sp_ub = mu_ub[mu_ub['rxn_id'] == rxn.id]['shadow_prices'].iloc[0]
            sp_lb = mu_lb[mu_lb['rxn_id'] == rxn.id]['shadow_prices'].iloc[0]

            rxn_allocation_coefficient = rxn.flux *(sp_ub-sp_lb)/obj_value
            # add new_row to dataframe
            self.allocation_coefficients.loc[len(self.allocation_coefficients)] = [rxn.id, '', 'rxn', rxn_allocation_coefficient]

            #get all enzymes related to the reaction
            if f'CE_{rxn.id}' in self.catalytic_events:
                for enzyme in self.catalytic_events.get_by_id(f'CE_{rxn.id}').enzymes:
                    #get the right row from the shadow price dataframes
                    sp_ec_f = mu_ec_f[mu_ec_f['rxn_id']== f'EC_{enzyme.id}']['shadow_prices'].iloc[0]
                    sp_ec_b = mu_ec_b[mu_ec_b['rxn_id'] == f'EC_{enzyme.id}']['shadow_prices'].iloc[0]
                    e_fwd = self.enzyme_variables.get_by_id(enzyme.id).forward_variable.primal
                    e_rev = self.enzyme_variables.get_by_id(enzyme.id).reverse_variable.primal

                    #EC: enzyme constraint
                    enzyme_allocation_coefficient = (e_fwd * sp_ec_f - e_rev * sp_ec_b)/obj_value
                    #add new_row to dataframe
                    self.allocation_coefficients.loc[len(self.allocation_coefficients)] = [rxn.id, enzyme.id, 'enzyme',
                                                                                           enzyme_allocation_coefficient]

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
        tot_prot_constraint.ub = p_tot
        self.p_tot = p_tot
        self.solver.update()

    def change_sector_parameters(self, sector, slope:float, intercept:float, lin_rxn_id:str):
        # input in g/gDW
        print(f'Changing the slope and intercept of the {sector.id}')
        print(f'Changing slope from {sector.slope} to {slope}')
        print(f'Changing intercept from {sector.intercept} to {intercept}')
        sector.slope = slope *1e3 / (sector.mol_mass[0] * 1e-6)
        sector.intercept = intercept *1e3 / (sector.mol_mass[0] * 1e-6)
        var_0 = self.variables['R_' + sector.id + '_0']
        var_0.set_bounds(intercept, intercept)
        var = self.variables['R_' + sector.id]
        lin_rxn = self.reactions.get_by_id(lin_rxn_id)

        #removing the old constraint
        self.remove(self.constraints[sector.id])

        #add the new updated constraint
        constraint = Constraint(sector.id, compartment='sector')
        self.add_sector_constraints([constraint])

        self.constraints[constraint.id].set_linear_coefficients({
            var: 1, var_0.forward_variable: intercept,
            lin_rxn.forward_variable: -slope,
            lin_rxn.reverse_variable: slope
        })

        #update the sector object
        sector.variables = [var, var_0]
        sector.constraints = [constraint]

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
        if lower_bound is not None:
            if self.sensitivity:
                self.constraints[rxn_id + '_lb'].ub = -lower_bound
            else:
                self.reactions.get_by_id(rxn_id).lower_bound = lower_bound
        if upper_bound is not None:
            if self.sensitivity:
                self.constraints[rxn_id + '_ub'].ub = upper_bound
            else:
                self.reactions.get_by_id(rxn_id).upper_bound = upper_bound

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
        enzymes = catalytic_event.kcats.keys()
        return enzymes

    def get_reactions_with_enzyme_id(self, enz_id:str):
        """
        Returns Enzyme objects associated with the reaction id through CatalyticEvent objects
        :param rxn_id: str
            enzyme identifier (EC number)

        Returns
        -------
            DictList of Enzyme objects associated with the reaction
        """
        enzyme = self.catalytic_events.get_by_id(enz_id)
        enzymes = enzyme.rxn2kcat.keys()
        return enzymes

    def modify_kcat(self, enzyme_id:str, kcats:dict):
        """
        Change the turnover number of the enzyme for a specific reaction
        :param enzyme_id: Enzyme identifier
        :param kcats: dict with reaction id, kcat key, value pairs. kcat is again a dict with direction, kcat value
        key, value pairs
        """
        try:
            enzyme = self.enzymes.get_by_id(enzyme_id)
            enzyme.change_kcat_values(kcats)
        except:
            warnings.warn(f'The enzyme {enzyme_id} does not exist in the model. The kcat can thus not be changed.')

    def remove_enzymes(self,
                       enzymes: Union[str, Enzyme, List[Union[str, Enzyme]]]
                       )-> None:
        """Remove enzymes from the model. Adapted from the cobra.core.remove_reactions() function.

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
                del self.enzyme_constraints[constraint.name]
            self.enzymes.remove(enzyme)
            enzyme._model = None

            #remove variable related to the enzyme
            self.enzyme_variables.remove(enzyme.enzyme_variable)
            forward = enzyme.variables['forward']
            reverse = enzyme.variables['reverse']

            self.remove_cons_vars(forward, reverse)

            #remove enzyme from catalytic event
            for ce in enzymes.catalytic_events:
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
            catalytic_event = self.catalytic_events.get_by_id('CE_' + rxn.id)
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

            try:
                catalytic_event = self.catalytic_events.get_by_id('CE_' + rxn.id)
                self.catalytic_events.remove(catalytic_event)


                #removing catalytic event from the enzymes
                for enzyme in catalytic_event.enzymes:
                    enzyme.remove_catalytic_event(catalytic_event)

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
                reactions : list or reaction or str
                    A list with reactions (`cobra.Reaction`), or their id's, to remove.
                    Reaction will be placed in a list. Str will be placed in a list and used to
                    find the reaction in the model.
                remove_orphans : bool, optional
                    Remove orphaned genes and metabolites from the model as
                    well (default False).
                """
        if isinstance(sectors, str) or hasattr(sectors, "id"):
            warnings.warn("need to pass in a list")
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
        rxn_glc_in = self.reactions.get_by_id(self.GLUCOSE_EXCHANGE_RXNID)
        #check if the model is reversible
        if self.GLUCOSE_EXCHANGE_RXNID[-1] == 'b':
            rxn_glc_in._lower_bound, rxn_glc_in._upper_bound = bound, bound
        else:
            rxn_glc_in._lower_bound, rxn_glc_in._upper_bound = -bound, -bound
        rxn_glc_in.update_variable_bounds()

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
               parameters for cobra.solvers may be available and specified with the
               appropriate keyword argument.
               """
        solution = super().optimize(objective_sense,raise_error)
        if self.sensitivity and self.solver.status == 'optimal':
            self.determine_control_coefficients()
        return solution
