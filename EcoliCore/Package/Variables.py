"""
Classes inherited from cobra.Reaction which enable addition of user-defined constraints
Account for enzyme allocated to forward and reverse reactions
"""
import cobra
from cobra import Reaction, Metabolite, DictList, Model
from cobra.exceptions import OptimizationError
from cobra.util.solver import check_solver_status
from typing import Optional, Union, Dict
from warnings import warn
import sys


class EnzymeVariable(Reaction):
    """
           EnzymeVariable is a class for holding information regarding the
           variable representing an enzyme in the model. For each reaction, the enzyme variables are
           summarized in a CatalyticEvent
           There are three different scenarios:
           - Enzyme complex: multiple enzymes together are associated with an EnzymeComplex object
           - isozymes: multiple enzymes independently associated with a single catalytic event
           - Other: a single enzyme is associated with a single catalytic event

           Parameters
           ----------
           kcats2rxns: Dict
               A Dict with reaction_id, kcat key, value pairs to connect the
               enzyme with the associated reaction the kcat is another dict with 'f' and 'b'
               for the forward and backward reactions respectively.
           id : str, optional
               The identifier to associate with this enzyme (default None)
           name : str, optional
               A human-readable name for the reaction (default "").
           subsystem : str, optional
               Subsystem where the reaction is meant to occur (default "").
           lower_bound : float
               The lower flux bound (default 0.0).
           upper_bound : float, optional
               The upper flux bound (default None).
           **kwargs:
               Further keyword arguments are passed on to the parent class.
           """
    DEFAULT_ENZYME_MOL_MASS = 3.947778784340140e04  # mean enzymes mass E.coli [g/mol]

    def __init__(
            self,
            kcats2rxns: Dict,
            molmass: Union[int, float] = DEFAULT_ENZYME_MOL_MASS,
            id: Optional[str] = None,  # ID of enzymatic reaction,
            name: str = "",
            upper_bound: Optional[float] = None,
            **kwargs
    ):
        super().__init__(
            id=id,
            name=name,
            subsystem='Enzymes',
            lower_bound=-upper_bound,
            upper_bound=upper_bound,
            **kwargs
        )
        self.kcats = kcats2rxns
        self.molmass = molmass
        self.rxn_ids = [rxn for rxn in kcats2rxns.keys()]
        self.enzyme = None #store the enzyme associated to the enzyme variable
        self.catalytic_events = DictList()
        self.reactions = DictList()
        self.constraints = {}  # store IDs of constraint the enzyme is associated with
        self._model = None
        self.annotation = {'type': 'Constraint'}
        self.variables = dict()

    @property
    def kcat_values(self):
        """returns a dictionary with kcat values and reactions
        """
        return self.kcats

    @property
    def flux(self) -> float:
        """
        Get the flux value in the most recent solution.
        Flux is the primal value of the corresponding variable in the model.
        Returns
        -------
        flux: float
            Flux is the primal value of the corresponding variable in the model.
        Warnings
        --------
        * Accessing reaction fluxes through a `Solution` object is the safer,
          preferred, and only guaranteed to be correct way. You can see how to
          do so easily in the examples.
        * Reaction flux is retrieved from the currently defined
          `self._model.solver`. The solver status is checked but there are no
          guarantees that the current solver state is the one you are looking
          for.
        * If you modify the underlying model after an optimization, you will
          retrieve the old optimization values.
        Raises
        ------
        RuntimeError
            If the underlying model was never optimized beforehand or the
            reaction is not part of a model.
        OptimizationError
            If the solver status is anything other than 'optimal'.
        AssertionError
            If the flux value is not within the bounds.
        Examples
        --------
        >>> from cobra.io import load_model
        >>> model = load_model("textbook")
        >>> solution = model.optimize()
        >>> model.variables.PFK.flux
        7.477381962160283
        >>> solution.fluxes.PFK
        7.4773819621602833
        """
        try:
            check_solver_status(self._model.solver.status)
            return self.forward_variable.primal + self.reverse_variable.primal
        except AttributeError:
            raise RuntimeError(f"reaction '{self.id}' is not part of a model")
        # Due to below all-catch, which sucks, need to reraise these.
        except (RuntimeError, OptimizationError) as err:
            raise err
        # Would love to catch CplexSolverError and GurobiError here.
        except Exception as err:
            raise OptimizationError(
                f"Likely no solution exists. Original solver message: {str(err)}."
            ) from err

    @property
    def concentration(self) -> float:
        """
        Get the enzyme concentration value of the most recent solution.
        The enzyme concentration equals the flux value

        Returns
        -------
        float
            enzyme concentration [mmol/gDW]
        """
        return self.flux

    @property
    def model(self):
        return self._model

    @model.setter
    def model(self, model):
        sys.path.append('../Package/')
        from Enzyme import Enzyme
        self._model = model
        # setting up the relations to the model
        # add enzyme instance
        if self.id not in self._model.enzymes:
            enzyme = Enzyme(id =self.id,
                            rxn2kcat = self.kcat,
                            molmass=self.molmass
                            )
            self._model.add_enzymes([enzyme])
        self.enzyme = self._model.enzymes.get_by_id(self.id)
        self.constraints = self.enzyme._constraints
        self._model.enzyme_variables.append(self)

        # create new forward and reverse variables
        forward_variable = self._model.problem.Variable(name=self.id, ub=self.upper_bound, lb=0)
        reverse_variable = self._model.problem.Variable(name=self.reverse_id, ub=-self.lower_bound, lb=0)
        self.variables = {'forward_variable': forward_variable, 'reverse_variable': reverse_variable}
        self._model.add_cons_vars([forward_variable, reverse_variable])


        # add catalytic event and reaction instances
        for rxn_id in self.rxn_ids:
            if not rxn_id in self._model.reactions:
                # create new reaction and add to model
                rxn = cobra.Reaction(id=rxn_id)
                self._model.add_reactions([rxn])
            self.reactions.append(self._model.reactions.get_by_id(rxn_id))
            #create a catalytic event if it doesn't exist already
            if not f'CE_{rxn_id}' in self._model.catalytic_events:
                kcats2enzymes = {self: self.kcats[rxn_id]}
                ce = CatalyticEvent(id = f'CE_{rxn_id}',
                                    kcats2enzymes=kcats2enzymes,
                                    rxn_id=rxn_id)
                self._model.catalytic_events.append(ce)
            # if catalytic event exist, add the enzyme to it
            else:
                self._model.catalytic_events.get_by_id(f'CE_{rxn_id}').add_enzymes({self.enzyme: self.kcats[rxn_id]})
            # if catalytic event is not related to the enzyme variable yet, add it.
            if f'CE_{rxn_id}' not in self.catalytic_events:
                self.catalytic_events.append(self._model.catalytic_events.get_by_id(f'CE_{rxn_id}'))

    @property
    def forward_variable(self):
        if self._model is not None:
            return self._model.variables[self.id]
        else:
            return self.variables['forward_variable']

    @property
    def reverse_variable(self):
        if self._model is not None:
            return self._model.variables[self.reverse_id]
        else:
            return self.variables['reverse_variable']

    def add_catalytic_events(self, catalytic_events:list, kcats:list):
        """
        Adding a catalytic event to an enzyme variable

        Parameters
        ----------
        catalytic_events: list
            Catalytic events to add
        kcats:list
            A list with dicts containing direction, kcat key value pairs
        """

        for i, ce in enumerate(catalytic_events):
            if ce in self.catalytic_events:
                warn(f'Catalytic event {ce.id} is already associated with enzyme variable {self.id}. '
                     f'Continue with other catalytic events')
            else:
                ce.enzyme_variables.append(self)
                self.catalytic_events.append(ce)
                self.add_reactions({ce.rxn_id: kcats[i]})

    def add_reactions(self, reaction_kcat_dict: dict):
        """
        Add reactions to the enzyme variable and create bindings to the related model.
        If there are multiple reactions related to a single enzyme, this is an isozyme.

        Parameters
        ----------
        reaction_kcat_dict: nested dict
            A Dict with the reaction_id, kcat key, value pairs to connect the
            enzyme with the associated reaction the kcat is another dict with 'f' and 'b'
            for the forward and backward reactions respectively.


        """
        for reaction, kcat in reaction_kcat_dict.items():
            # check if the enzyme is already associated to the catalytic event
            try:
                self.reactions.get_by_id(reaction.id)
                warn(
                    f'Reaction {reaction.id} is already associated with the enzyme {self.id}. The enzyme variable will be updated')
                self.change_kcat_values(kcat)
                return
            except:
                pass

            self.kcats[reaction] = kcat

            if self._model is None:
                continue

            # check if enzyme is in the model
            try:
                self._model.reactions.get_by_id(reaction)
            # if not: add the enzyme to the model
            except:
                rxn = Reaction(id = reaction)
                self._model.add_reactions([rxn])

            rxn = self._model.reactions.get_by_id(reaction)
            self.reactions.append(rxn)

            for direction in kcat.keys():
                if direction != 'f' and direction != 'b':
                    warn(f'Invalid direction {direction} for kcat value for enzyme variable {self.id}! Skip!')
                    continue

            # add enzyme to catalytic event and the related variable
            for direction, kcatvalue in kcat.items():
                coeff = kcatvalue * 3600 * 1e-6
                # add enzyme to the associated reaction with kinetic constants
                # and relate enzyme to the catalytic event
                if direction == 'f':
                    self.constraints[f'EC_{self.id}_{direction}'].set_linear_coefficients({
                        rxn.forward_variable: 1 / coeff,
                        self.forward_variable: -1
                    })

                elif direction == 'b':
                    self.constraints[f'EC_{self.id}_{direction}'].set_linear_coefficients({
                        rxn.reverse_variable: 1 / coeff,
                        self.reverse_variable: -1
                    })

    def remove_reactions(self, reaction_list: list):
        """
        Remove reaction from the enzyme variable and remove the reaction from the
        constraint expressions related to the enzyme

        Parameters
        ----------
        reaction_list: list
            A list with Cbra.Reaction objects which should be removed. If a list of identifiers (str)
            is provided, the corresponding enzyme will be obtained from the EnzymeVariables.reaction attribute
        """
        # check the input
        if not hasattr(reaction_list, "__iter__"):
            enzyme_list = [reaction_list]
        if len(reaction_list) == 0:
            return None

        # check whether the input is an PAModelpy.Package.Enzyme or string and find the corresponding enzyme if needed
        for i, rxn in enumerate(reaction_list):
            if isinstance(rxn, str):
                try:
                    reaction_list[i] = self.reactions.get_by_id(rxn)
                except:
                    print(
                        f'Reaction {rxn} is not associated with the enzyme variable {self.id}. This reaction cannot be removed. \n')
                    pass

        for rxn in reaction_list:
            # remove from kcat dict
            del self.kcats[rxn.id]
            # remove from reactions dictlist
            self.reaction.remove(rxn)
            # set coefficient in constraint to 0
            for constraint in self.enzyme._constraints.values():
                self.constraints[constraint.name] = constraint
                coeff = 0
                # set coefficients to 0
                if constraint.name[-1] == 'f':
                    constraint.set_linear_coefficients({
                        rxn.forward_variable: coeff,
                        self.forward_variable: 0
                    })

                elif constraint.name[-1] == 'b':
                    constraint.set_linear_coefficients({
                        rxn.reverse_variable: coeff,
                        self.reverse_variable: 0
                    })
                # remove constraint from list of r=constraints
                del self.constraints[constraint.name]

    def change_kcat_values(self, reaction_kcat_dict: dict):
        """changes kcat values for the enzyme variable
        Parameters
        ----------
        reaction_kcat_dict: nested Dict
            A Dict with Cobra.Reaction, kcat key, value pairs to connect the
            enzyme with the associated reaction the kcat is another dict with 'f' and 'b'
            for the forward and backward reactions respectively.

        """
        # apply changes to internal dicts (one by one to avoid deleting kcat values)
        kcats_change = {}
        for rxn, kcat_dict in reaction_kcat_dict.items():
            # save change in dict
            self.kcats[rxn.id] = kcat_dict
            for direction, kcat in kcat_dict.items():
                if direction != 'f' and direction != 'b':
                    warn(f'Invalid direction {direction} for kcat value for enzyme variable {self.id}! Skip!')
                    continue

                kcats_change[direction] = kcat

            # is enzyme variable already integrated into a model
            if self._model is None:
                warn(f'Catalytic event {self.id} is not integrated into a model!')

            for direction, kcat in kcats_change.items():
                # get constraint
                constraint_id = f'EC_{rxn.id}_{direction}'
                constraint = self.enzyme._constraints[constraint_id]
                # change kcat value in the constraint
                coeff = kcat * 3600 * 1e-6
                if direction == 'f':
                    constraint.set_linear_coefficients({
                        self.reaction.forward_variable: 1 / coeff
                    })
                elif direction == 'b':
                    constraint.set_linear_coefficients({
                        self.reaction.reverse_variable: 1 / coeff
                    })
            self._model.solver.update()


class CatalyticEvent():
    """
        CatalyticEvent is a class for holding information regarding the
        catalysis of a Reaction in a cobra.Model object. It serves as an interface
        between the metabolic reaction and the associated enzyme constraints and variables.
        There are three different scenarios:
        - Enzyme complex: multiple enzymes together are associated with an EnzymeComplex object TODO
        - isozymes: multiple enzymes independently associated with a single catalytic event
        - Other: a single enzyme is associated with a single catalytic event

        Parameters
        ----------
        kcats2enzymes: Dict
            A Dict with enzyme, kcat key, value pairs to connect the
            enzyme with the associated reaction the kcat is another dict with 'f' and 'b'
            for the forward and backward reactions respectively.
        id : str, optional
            The identifier to associate with this catalytic event (default None).
        rxn_id: str, Optional
            The reaction with which this catalytic event is associated
        name : str, optional
            A human-readable name for the reaction (default "").
        """
    def __init__(
            self,
            kcats2enzymes: Dict,
            id: Optional[str] = None,   # ID of enzymatic reaction
            rxn_id: str = '', # ID of reaction
            name: str = ""
    ):
        #identification
        self.id = id
        self.name = name

        #relation to reaction
        self.rxn_id = rxn_id
        self.rxn = None

        #relation to enzymes
        self.kcats = kcats2enzymes
        self.enzymes = DictList()
        for enzyme in kcats2enzymes.keys():
            self.enzymes.append(enzyme)
        self.enzyme_variables = DictList()

        #other attributes
        self.constraints = {} # store IDs of constraint the catalytic event is associated with
        self.variables = dict()
        self._model = None
        self.annotation = {'type':'Constraint'}


    @property
    def kcat_values(self):
        """returns a dictionary with kcat values and enzymes
        """
        return self.kcats

    @property
    def flux(self) -> float:
        """
        Get the flux value in the most recent solution.
        Flux is the primal value of the corresponding variable in the model.
        Returns
        -------
        flux: float
            Flux is the primal value of the corresponding variable in the model.
        Warnings
        --------
        * Accessing reaction fluxes through a `Solution` object is the safer,
          preferred, and only guaranteed to be correct way. You can see how to
          do so easily in the examples.
        * Reaction flux is retrieved from the currently defined
          `self._model.solver`. The solver status is checked but there are no
          guarantees that the current solver state is the one you are looking
          for.
        * If you modify the underlying model after an optimization, you will
          retrieve the old optimization values.
        Raises
        ------
        RuntimeError
            If the underlying model was never optimized beforehand or the
            reaction is not part of a model.
        OptimizationError
            If the solver status is anything other than 'optimal'.
        AssertionError
            If the flux value is not within the bounds.
        Examples
        --------
        >>> from cobra.io import load_model
        >>> model = load_model("textbook")
        >>> solution = model.optimize()
        >>> model.variables.PFK.flux
        7.477381962160283
        >>> solution.fluxes.PFK
        7.4773819621602833
        """
        try:
            check_solver_status(self._model.solver.status)
            total_flux = 0
            for enzyme_variable in self.enzyme_variables:
                fwd_flux = enzyme_variable.variables['forward'].primal
                rev_flux = enzyme_variable.variables['reverse'].primal
                total_flux += fwd_flux + rev_flux
            return total_flux
        except AttributeError:
            raise RuntimeError(f"reaction '{self.id}' is not part of a model")
        # Due to below all-catch, which sucks, need to reraise these.
        except (RuntimeError, OptimizationError) as err:
            raise err
        # Would love to catch CplexSolverError and GurobiError here.
        except Exception as err:
            raise OptimizationError(
                f"Likely no solution exists. Original solver message: {str(err)}."
            ) from err

    @property
    def concentration(self) -> float:
        """
        Get the enzyme concentration value of the most recent solution.
        The enzyme concentration equals the flux value

        Returns
        -------
        float
            enzyme concentration [mmol/gDW]
        """
        return self.flux

    @property
    def model(self):
            return self._model

    @model.setter
    def model(self, model):
        self._model = model
        #add reaction instance
        if self.rxn_id in self._model.reactions:
            self.rxn = self._model.reactions.get_by_id(self.rxn_id)
        else: #create new reaction and add to model
            rxn = cobra.Reaction(id = self.rxn_id)
            self._model.add_reactions([rxn])
            self.rxn = rxn
        #add enzymes to the model if they are not there already
        for enzyme in self.enzymes:
            if enzyme in self._model.enzymes:
                enzyme_model = self._model.enzymes.get_by_id(enzyme.id)
                if self.rxn_id in enzyme_model.rxn2kcat.keys():
                    if self not in enzyme_model.catalytic_events:
                        enzyme_model.add_catalytic_event(self, kcats= {enzyme: enzyme.rxn2kcat[self.rxn_id]})
                else:
                    print(f'Reaction {self.rxn_id} is not related to enzyme {enzyme.id}')
            else:
                self._model.add_enzymes([enzyme])


    def add_enzymes(self,enzyme_kcat_dict: dict):
        """
        Add enzymes to the catalytic event and create bindings to the related model.
        The enzymes in the enzyme_kcat_dict are individual isozymes. Enzyme complexes
        should be added as an EnzymeComplex object with a single kcat value.

        Parameters
        ----------
        enzyme_kcat_dict: nested dict
            A Dict with enzyme, kcat key, value pairs to connect the
            enzyme with the associated reaction the kcat is another dict with 'f' and 'b'
            for the forward and backward reactions respectively.


        """
        for enzyme, kcat in enzyme_kcat_dict.items():
            # check if the enzyme is already associated to the catalytic event
            if enzyme in self.enzymes:
                self.enzymes.get_by_id(enzyme.id)
                # warn(f'Enzyme {enzyme.id} is already associated with catalytic event {self.id}. This enzyme will be updated')
                self.change_kcat_values({enzyme.id: enzyme.get_kcat_values(rxn_ids = [self.rxn_id])[self.rxn_id]})
                continue

            self.enzymes.append(enzyme)
            self.kcats[enzyme] = kcat

            if self._model is None:
                continue
            #check if enzyme is in the model
            try: self._model.enzymes.get_by_id(enzyme.id)
            #if not: add the enzyme to the model
            except: self._model.add_enzymes([enzyme])

            for direction in kcat.keys():
                if direction != 'f' and direction != 'b':
                    warn(f'Invalid direction {direction} for kcat value for enzyme variable {self.id}! Skip!')
                    continue

            #add enzyme to catalytic event and the related variable
            #get enzyme variable
            enzyme_var = self._model.enzyme_variables.get_by_id(enzyme.id)
            self.enzyme_variables.append(enzyme_var)

            #connect the enzyme variable to the enzyme in the model and the reaction
            for direction, kcatvalue in kcat.items():
                self.constraints[f'EC_{enzyme.id}_{direction}'] = enzyme._constraints[f'EC_{enzyme.id}_{direction}']
                coeff = kcatvalue * 3600 * 1e-6
                #add enzyme to the associated reaction with kinetic constants
                #and relate enzyme to the catalytic event
                if direction == 'f':
                    self.constraints[f'EC_{enzyme.id}_{direction}'].set_linear_coefficients({
                        self.rxn.forward_variable: 1/coeff,
                        enzyme_var.forward_variable: -1
                        })

                elif direction == 'b':
                    self.constraints[f'EC_{enzyme.id}_{direction}'].set_linear_coefficients({
                        self.rxn.reverse_variable: 1/coeff,
                        enzyme_var.reverse_variable: -1
                        })

    def remove_enzymes(self, enzyme_list: list):
        """
        Remove enzymes from the catalytic event and remove catalytic event from the
        constraint expressions related to the enzyme

        Parameters
        ----------
        enzyme_list: list
            A list with PAModelpy.Package.Enzyme objects which should be removed. If a list of identifiers (str)
            is provided, the corresponding enzyme will be obtained from the CatalyticEvent.enzymes attribute
        """
        #check the input
        if not hasattr(enzyme_list, "__iter__"):
            enzyme_list = [enzyme_list]
        if len(enzyme_list) == 0:
            return None

        # check wether the input is an PAModelpy.Package.Enzyme or string and find the corresponding enzyme if needed
        for i, enz in enumerate(enzyme_list):
            if isinstance(enz, str):
                try:
                    enzyme_list[i] = self.enzymes.get_by_id(enz)
                except:
                    print(f'Enzyme {enz} is not associated with the catalytic event {self.id}. This enzyme will not be removed. \n')
                    pass

        for enz in enzyme_list:
            #remove from kcat dict
            del self.kcats[enz]
            #remove from enzymes dictlist
            self.enzymes.remove(enz)
            #remove enzyme variable from dictlist
            enzyme_var = self._model.enzyme_variables.get_by_id(enz.id)
            self.enzyme_variables.remoce(enzyme_var)
            #set coefficient in constraint to 0
            for constraint in enz._constraints.values():
                self.constraints[constraint.name] = constraint
                coeff = 0
                #set coefficients to 0
                if constraint.name[-1] == 'f':
                    constraint.set_linear_coefficients({
                        self.rxn.forward_variable: coeff
                        })

                elif constraint.name[-1] == 'b':
                    constraint.set_linear_coefficients({
                        self.rxn.reverse_variable: coeff
                        })
                #remove constraint from list of r=constraints
                del self.constraints[constraint.name]

    def change_kcat_values(self, enzyme_kcat_dict : dict):
        """changes kcat values for the enzyme variable
        Parameters
        ----------
        enzyme_kcat_dict: nested Dict
            A Dict with enzyme, kcat key, value pairs to connect the
            enzyme with the associated reaction the kcat is another dict with 'f' and 'b'
            for the forward and backward reactions respectively.
        
        """
        # apply changes to internal dicts (one by one to avoid deleting kcat values)
        kcats_change = {}
        for enzyme, kcat_dict in enzyme_kcat_dict.items():
            # save change in dict
            self.kcats[enzyme] = kcat_dict
            for direction, kcat in kcat_dict.items():
                if direction != 'f' and direction != 'b':
                    warn(f'Invalid direction {direction} for kcat value for enzyme variable {self.id}! Skip!')
                    continue
            
                kcats_change[direction] = kcat
        
            # is enzyme variable already integrated into a model
            if self._model is None:
                # warn(f'Catalytic event {self.id} is not integrated into a model!')
                return

            enzyme_var = self._model.enzyme_variables.get_by_id(enzyme)
            enzyme_obj = self._model.enzymes.get_by_id(enzyme)

            for direction, kcat in kcats_change.items():
                #change enzyme variable
                enzyme_var.kcats[self.rxn_id][direction] = kcat
                # get constraint
                constraint_id = f'EC_{enzyme}_{direction}'
                constraint = enzyme_obj._constraints[constraint_id]
                # change kcat value in the constraint
                coeff =  kcat * 3600 * 1e-6
                if direction == 'f':
                    constraint.set_linear_coefficients({
                        self.rxn.forward_variable: 1/coeff
                        })   
                elif direction == 'b':    
                    constraint.set_linear_coefficients({
                        self.rxn.reverse_variable: 1/coeff
                        })
            self._model.solver.update()



class Constant(Reaction):

    """
    Constant is a class for holding information regarding
    a constant which can be used to define a constraint in a cobra.Model object.

    Constants are by Reactions which have a constant flux.
    They are related to Constraints which include a constant value.

    Parameters
    ----------
    id : str, optional
        The identifier to associate with this reaction (default None).
    name : str, optional
        A human readable name for the reaction (default "").
    subsystem : str, optional
        Subsystem where the reaction is meant to occur (default "").
    coefficient: float, not optional
        constant value which defines the constant
    direction: str, optional
        Determine in which direction the constraint is related to the constant (addition `+` or substraction `-`)
    constraint: Metabolite or Constraint, not optional
        Constraint to which the constant is related
    **kwargs:
        Further keyword arguments are passed on to the parent class.
    """
    def __init__(
            self,
            id: Optional[str] = None,
            name: str = "",
            coefficient: float = 0.0,
            direction: str = '+',
            constraint: Metabolite = None,
            **kwargs,
    ):
        super().__init__(
            id=id,
            name=name,
            subsystem='Enzymes',
            lower_bound=coefficient,
            upper_bound=coefficient,
            **kwargs
        )
        self._model = None
        self.coeff = 1
        if direction == '-': self.coeff = -1

        self.annotation = {'type':'Constraint'}

    @property
    def model(self):
        return self._model

    @model.setter
    def model(self, model):
        self._model = model
        # create new variables for total_protein objects
        forward_variable = self._model.problem.Variable(self.id, lb=self.lower_bound, ub =self.upper_bound)
        self._model.add_cons_vars([forward_variable])


class Total_protein_variable(Reaction):
    def __init__(
            self,
            id: Optional[str] = None,
            name: str = "",
            coefficient: float = 0.0,
            direction: str = '+',
            constraint: Metabolite = None,
            **kwargs,
    ):
        super().__init__(
            id=id,
            name=name,
            subsystem='Enzymes',
            lower_bound=0,
            upper_bound=coefficient,
            **kwargs
        )
        self._model = None
        self.coeff = 1
        if direction == '-': self.coeff = -1

    @property
    def model(self):
            return self._model

    @model.setter
    def model(self, model):
        self._model = model
        # create new variables for total_protein objects
        forward_variable = self._model.problem.Variable(self.id, lb = self.lower_bound, ub = self.upper_bound)
        self._model.add_cons_vars([forward_variable])

