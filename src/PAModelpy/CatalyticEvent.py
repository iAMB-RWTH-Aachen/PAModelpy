"""
CatalyticEvent object which relates Reaction variables to the EnzymeVariable and Enzyme objects.
It contains multiple functions which enable easy mapping and handling of one Event of catalysis
(e.g. one conversion of substrate to product, can be catalyzed by multiple enzymes)
"""
import cobra
from cobra import DictList, Object
from cobra.exceptions import OptimizationError
from cobra.util.solver import check_solver_status
from optlang.symbolics import Zero
from typing import Optional, Dict
from warnings import warn
from copy import copy, deepcopy

class CatalyticEvent(Object):
    """
        CatalyticEvent is a class for holding information regarding the
        catalysis of a Reaction in a cobra.Model object. It serves as an interface
        between the metabolic reaction and the associated enzyme constraints and variables.
        There are three different scenarios:
        - Enzyme complex: multiple enzymes together are associated with an EnzymeComplex object
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
                self.constraints ={**self.constraints, **enzyme._constraints}
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
                # print(enzyme)
                # warn(f'Enzyme {enzyme.id} is already associated with catalytic event {self.id}. This enzyme will be updated')
                self.change_kcat_values({enzyme.id: enzyme.get_kcat_values(rxn_ids = [self.rxn_id])})
                # print(self.enzymes)
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

            #add constraints to the catalytic event
            self.constraints = {**self.constraints, **enzyme._constraints}

            #connect the enzyme variable to the enzyme in the model and the reaction
            for direction, kcatvalue in kcat.items():
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
            #remove enzyme variable from dictlist if it still exists
            if self.enzyme_variables.has_id(enz.id):
                enzyme_var = self.enzyme_variables.get_by_id(enz.id)
                self.enzyme_variables.remove(enzyme_var)
            #set coefficient in constraint to 0
            for constraint in set([cons for name, cons in self.constraints.items() if enz.id in name]):
                # self.constraints[constraint.name] = constraint
                # coeff = 0
                # #set coefficients to 0
                # if constraint.name[-1] == 'f':
                #     constraint.set_linear_coefficients({
                #         self.rxn.forward_variable: coeff
                #         })
                #
                # elif constraint.name[-1] == 'b':
                #     constraint.set_linear_coefficients({
                #         self.rxn.reverse_variable: coeff
                #         })
                if constraint in self._model.constraints.values():
                    self._model.remove_cons_vars([constraint])
                #remove constraint from list of r=constraints
                del self.constraints[constraint.name]

            #if there are no enzymes associated to the reaction anymore, the reaction flux will be 0
            if len(self.enzymes)==0:
                no_enz_constraint_f = self._model.problem.Constraint(Zero, name = f'{self.rxn_id}_no_enzyme_f', lb=0, ub=0)
                no_enz_constraint_b = self._model.problem.Constraint(Zero, name=f'{self.rxn_id}_no_enzyme_b', lb=0, ub=0)
                self._model.add_cons_vars([no_enz_constraint_f,no_enz_constraint_b])

                self._model.constraints[no_enz_constraint_f.name].set_linear_coefficients({
                    self.rxn.forward_variable:1
                })
                self.constraints[no_enz_constraint_f.name] = no_enz_constraint_f

                self._model.constraints[no_enz_constraint_b.name].set_linear_coefficients({
                    self.rxn.forward_variable: 1
                })
                self.constraints[no_enz_constraint_b.name] = no_enz_constraint_b

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
            if enzyme_var not in self.enzyme_variables: self.enzyme_variables.add(enzyme_var)
            enzyme_obj = self._model.enzymes.get_by_id(enzyme)
            enzyme_var.change_kcat_values({self.rxn: kcat_dict})
            for direction, kcat in kcats_change.items():
                #change enzyme variable
                enzyme_var.kcats[self.rxn_id][direction] = kcat
                # get constraint
                constraint_id = f'EC_{enzyme}_{direction}'
                constraint = enzyme_obj._constraints[constraint_id]
                # change kcat value in the constraint
                coeff =  kcat * 3600 * 1e-6
                if direction == 'f':
                    self._model.constraints[constraint_id].set_linear_coefficients({
                        self.rxn.forward_variable: 1/coeff
                        })
                elif direction == 'b':
                    self._model.constraints[constraint_id].set_linear_coefficients({
                        self.rxn.reverse_variable: 1/coeff
                        })
            self._model.solver.update()

    def __copy__(self) -> 'CatalyticEvent':
        """ Copy the CatalyticEvent
        :return: CatalyticEvent:
        A new CatalyticEvent that is a copy of the original CatalyticEvent
        """

        cop = copy(super(CatalyticEvent, self))
        return cop

    def __deepcopy__(self, memo: dict) -> 'CatalyticEvent':
        """ Copy the CatalyticEvent with memo

        :param: memo:dict:
        Automatically passed parameter

        :return: CatalyticEvent:
        A new CatalyticEvent that is a copy of the original CatalyticEvent with memo
        """

        cop = deepcopy(super(CatalyticEvent, self), memo)
        return cop