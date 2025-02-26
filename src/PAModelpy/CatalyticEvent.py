"""
CatalyticEvent object which relates Reaction variables to the EnzymeVariable and Enzyme objects.
It contains multiple functions which enable easy mapping and handling of one Event of catalysis
(e.g. one conversion of substrate to product, can be catalyzed by multiple enzymes)
"""
import cobra
from cobra import DictList, Object, Reaction
from cobra.exceptions import OptimizationError
from cobra.util.solver import check_solver_status
from optlang.symbolics import Zero
from typing import Optional, Dict
from warnings import warn
from copy import copy, deepcopy
import re
from collections import defaultdict


class CatalyticEvent(Object):
    """
   CatalyticEvent is a class for holding information regarding the
        catalysis of a Reaction in a cobra.Model object. It serves as an interface
        between the metabolic reaction and the associated enzyme constraints and variables.

    Notes:
        There are three different scenarios:
        - Enzyme complex: multiple enzymes together are associated with an EnzymeComplex object
        - Isozymes: multiple enzymes independently associated with a single catalytic event
        - Other: a single enzyme is associated with a single catalytic event

    Parameters:
        kcats2enzymes (dict): A dictionary with enzyme, kcat key, value pairs to connect the enzyme with the associated reaction.
            The kcat is another dictionary with 'f' and 'b' for the forward and backward reactions, respectively.
        id (str, optional): The identifier to associate with this catalytic event (default None).
        rxn_id (str, optional): The reaction with which this catalytic event is associated.
        name (str, optional): A human-readable name for the reaction (default "").
    """

    def __init__(
        self,
        kcats2enzymes: Dict,
        id: Optional[str] = None,  # ID of enzymatic reaction
        rxn_id: str = "",  # ID of reaction
        name: str = "",
    ):
        # identification
        self.id = id
        self.name = name

        # relation to reaction
        self.rxn_id = rxn_id
        self.rxn = None
        self.catalytic_reactions = DictList()

        # relation to enzymes
        self.kcats = kcats2enzymes
        self.enzymes = DictList()
        for enzyme in kcats2enzymes.keys():
            self.enzymes.append(enzyme)
        self.enzyme_variables = DictList()

        # other attributes
        self.constraints = {} # store IDs of constraint the catalytic event is associated with
        self.variables = dict()
        self._model = None
        self.annotation = {"type": "Constraint"}

    @property
    def kcat_values(self):
        """returns a dictionary with kcat values and enzymes"""
        return self.kcats

    @property
    def flux(self) -> float:
        """
        Get the flux value in the most recent solution.

        Flux is the primal value of the corresponding variable in the model.

        Returns:
            flux (float): Flux is the primal value of the corresponding variable in the model.

        Warnings:
            * Accessing reaction fluxes through a `Solution` object is the safer,
              preferred, and only guaranteed to be correct way. You can see how to
              do so easily in the examples.
            * Reaction flux is retrieved from the currently defined
              `self._model.solver`. The solver status is checked but there are no
              guarantees that the current solver state is the one you are looking
              for.
            * If you modify the underlying model after an optimization, you will
              retrieve the old optimization values.

        Raises:
            RuntimeError: If the underlying model was never optimized beforehand or the
                reaction is not part of a model.
            OptimizationError: If the solver status is anything other than `optimal`.
            AssertionError: If the flux value is not within the bounds.

        Examples:
            ```
            >>> from cobra.io import load_model
            >>> model = load_model("textbook")
            >>> solution = model.optimize()
            >>> model.variables.PFK.flux
            7.477381962160283
            >>> solution.fluxes.PFK
            7.4773819621602833
            ```
        """

        try:
            check_solver_status(self._model.solver.status)
            total_flux = 0
            for enzyme_variable in self.enzyme_variables:
                fwd_flux = enzyme_variable.variables["forward"].primal
                rev_flux = enzyme_variable.variables["reverse"].primal
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
        The enzyme concentration equals the flux value.

        Returns:
            float: Enzyme concentration in [mmol/gDW].
        """

        return self.flux

    @property
    def model(self):
        return self._model

    @model.setter
    def model(self, model):
        self._model = model

        # add reaction instance
        if self.rxn_id in self._model.reactions:
            self.rxn = self._model.reactions.get_by_id(self.rxn_id)
        else:  # create new reaction and add to model
            rxn = cobra.Reaction(id=self.rxn_id)
            self._model.add_reactions([rxn])
            self.rxn = rxn
        # add reaction constraint
        ce_constraint = self._model.problem.Constraint(Zero, name= self.id, lb = 0, ub=0)
        self._model.add_cons_vars([ce_constraint])
        self._model.constraints[self.id].set_linear_coefficients({
            self.rxn.forward_variable: 1,
            self.rxn.reverse_variable: -1,
        })
        self.constraints[self.id] = ce_constraint

        for enzyme in self.enzymes:
            if enzyme in self._model.enzymes:
                self.constraints = {**self.constraints, **enzyme._constraints}
                enzyme._constraints[self.id] = ce_constraint
                enzyme_model = self._model.enzymes.get_by_id(enzyme.id)

                if self.rxn_id in enzyme_model.rxn2kcat.keys():
                    if self not in enzyme_model.catalytic_events:
                        enzyme_model.add_catalytic_event(
                            self, kcats={enzyme: enzyme.rxn2kcat[self.rxn_id]}
                        )
                else:
                    print(
                        f"Reaction {self.rxn_id} is not related to enzyme {enzyme.id}"
                    )
            else:
                self._model.add_enzymes([enzyme])

            # add enzyme-reaction association
            self.add_enzyme_reaction_association(enzyme)



    def add_enzymes(self, enzyme_kcat_dict: dict):
        """
        Add enzymes to the catalytic event and create bindings to the related model.
        The enzymes in the enzyme_kcat_dict are individual isozymes. Enzyme complexes
        should be added as an EnzymeComplex object with a single kcat value.

        Parameters:
            enzyme_kcat_dict: Dict
                A nested dictionary with enzyme, kcat key, value pairs to connect the
                enzyme with the associated reaction. The kcat is another dictionary with `f` and `b`
                for the forward and backward reactions respectively.
        """
        # return lists back to dictlist after unpickling
        if isinstance(self.enzymes, list):
            self.enzymes = DictList(self.enzymes)
            self.enzyme_variables = DictList(self.enzyme_variables)

        for enzyme, kcat in enzyme_kcat_dict.items():
            # check if the enzyme is already associated to the catalytic event
            if self.catalytic_reactions.has_id(f'CE_{self.rxn_id}_{enzyme.id}'):
                continue
            if enzyme in self.enzymes:
                self.add_enzyme_reaction_association(enzyme)
                self.change_kcat_values(
                    {enzyme.id: enzyme.get_kcat_values(rxn_ids=[self.rxn_id+"_"+enzyme.id])})
                continue

            self.enzymes.append(enzyme)
            self.kcats[enzyme] = kcat

            if self._model is None:
                continue
            # check if enzyme is in the model
            try:
                self._model.enzymes.get_by_id(enzyme.id)
            # if not: add the enzyme to the model
            except:
                self._model.add_enzymes([enzyme])

            for direction in kcat.keys():
                if direction != "f" and direction != "b":
                    warn(
                        f"Invalid direction {direction} for kcat value for enzyme variable {self.id}! Skip!"
                    )
                    continue

            # add enzyme to catalytic event and the related variable
            # get enzyme variable
            enzyme_var = self._model.enzyme_variables.get_by_id(enzyme.id)
            self.enzyme_variables.append(enzyme_var)

            # add constraints to the catalytic event
            self.constraints = {**self.constraints, **enzyme._constraints}

            # add reaction-protein association
            self.add_enzyme_reaction_association(enzyme)

            #connect the enzyme variable to the enzyme in the model and the reaction
            for direction, kcatvalue in kcat.items():
                #add enzyme to the associated reaction with kinetic constants
                #and relate enzyme to the catalytic event
                if direction == 'f':
                    self.constraints[f'EC_{enzyme.id}_{direction}'].set_linear_coefficients({
                        enzyme_var.forward_variable: -1
                        })

                elif direction == 'b':
                    self.constraints[f'EC_{enzyme.id}_{direction}'].set_linear_coefficients({
                        enzyme_var.reverse_variable: -1
                        })

    def add_enzyme_reaction_association(self, enzyme):
        catalytic_reaction = Reaction(self.id + "_" + enzyme.id, lower_bound=-1e3, upper_bound=1e3)
        self._model.add_reactions([catalytic_reaction])
        self.catalytic_reactions.append(catalytic_reaction)
        self._model.constraints[self.id].set_linear_coefficients({
            catalytic_reaction.forward_variable: -1,
            catalytic_reaction.reverse_variable: 1,
        })
        self.add_catalytic_reaction_to_enzyme_constraint(catalytic_reaction, enzyme)

    def add_catalytic_reaction_to_enzyme_constraint(self, catalytic_reaction:Reaction,
                                                    enzyme):
        #remove relation to metabolic reaction
        for dir in ['f', 'b']:
            self._model._change_kcat_in_enzyme_constraint(self.rxn, enzyme.id,
                                                      dir, 0)
        kcat_dict = self.kcats[enzyme]
        for direction, kcatvalue in kcat_dict.items():
            self._model._change_kcat_in_enzyme_constraint(catalytic_reaction, enzyme.id,
                                                          direction, kcatvalue)

        #change rxn2kcat dict for correct referencing if the reaction was in there
        if self.rxn_id in enzyme.rxn2kcat.keys():
            del enzyme.rxn2kcat[self.rxn_id]
        enzyme.rxn2kcat[catalytic_reaction.id] = kcat_dict
        #also add catalytic reaction to the ActiveEnzymeSector
        self._model.sectors.ActiveEnzymeSector.rxn2protein[catalytic_reaction.id] = {enzyme.id: kcat_dict}


    def remove_enzymes(self, enzyme_list: list):
        """
        Remove enzymes from the catalytic event and remove the catalytic event from the
        constraint expressions related to the enzyme.

        Parameters:
            enzyme_list: List[Union[str, PAModelpy.Package.Enzyme]]
                A list with PAModelpy.Package.Enzyme objects to be removed. If a list of identifiers (str)
                is provided, the corresponding enzyme will be obtained from the CatalyticEvent.enzymes attribute.
        """

        # return lists back to dictlist after unpickling
        if isinstance(self.enzymes, list):
            self.enzymes = DictList(self.enzymes)
            self.enzyme_variables = DictList(self.enzyme_variables)

        # check the input
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
                    print(
                        f"Enzyme {enz} is not associated with the catalytic event {self.id}. This enzyme will not be removed. \n"
                    )
                    pass

            # remove from kcat dict
            del self.kcats[enz]
            # remove from enzymes dictlist
            self.enzymes.remove(enz)
            # remove enzyme variable from dictlist if it still exists
            if self.enzyme_variables.has_id(enz.id):
                enzyme_var = self.enzyme_variables.get_by_id(enz.id)
                self.enzyme_variables.remove(enzyme_var)
            # set coefficient in constraint to 0
            for constraint in set(
                [cons for name, cons in self.constraints.items() if (f'EC_{enz.id}' in name)|(f'{enz.id}_m' in name)]
            ):
                if constraint in self._model.constraints.values():
                    self._model.remove_cons_vars([constraint])
                # remove constraint from list of r=constraints
                del self.constraints[constraint.name]

            # if there are no enzymes associated to the reaction anymore, the reaction flux will be 0
            if len(self.enzymes) == 0 and self._model is not None:
                no_enz_constraint_f = self._model.problem.Constraint(
                    Zero, name=f"{self.rxn_id}_no_enzyme_f", lb=0, ub=0
                )
                no_enz_constraint_b = self._model.problem.Constraint(
                    Zero, name=f"{self.rxn_id}_no_enzyme_b", lb=0, ub=0
                )
                self._model.add_cons_vars([no_enz_constraint_f, no_enz_constraint_b])

                self._model.constraints[
                    no_enz_constraint_f.name
                ].set_linear_coefficients({self.rxn.forward_variable: 1})
                self.constraints[no_enz_constraint_f.name] = no_enz_constraint_f

                self._model.constraints[
                    no_enz_constraint_b.name
                ].set_linear_coefficients({self.rxn.forward_variable: 1})
                self.constraints[no_enz_constraint_b.name] = no_enz_constraint_b

    def change_kcat_values(self, enzyme_kcat_dict : dict):
        """changes kcat values for the enzyme variable

        Args:
        enzyme_kcat_dict: nested Dict
            A Dict with enzyme_id, kcat key, value pairs to connect the
            enzyme with the associated reaction the kcat is another dict with 'f' and 'b'
            for the forward and backward reactions respectively.
        
        """

        # apply changes to internal dicts (one by one to avoid deleting kcat values)
        kcats_change = {}
        kcats = defaultdict(dict, self.kcats)
        for enzyme, kcat_dict in enzyme_kcat_dict.items():
            # save change in dict
            kcats[enzyme] = {**kcats[enzyme],**kcat_dict}
            for direction, kcat in kcat_dict.items():
                if direction != "f" and direction != "b":
                    warn(
                        f"Invalid direction {direction} for kcat value for enzyme variable {self.id}! Skip!"
                    )
                    continue

                kcats_change[direction] = kcat

            # is enzyme variable already integrated into a model
            if self._model is None:
                # warn(f'Catalytic event {self.id} is not integrated into a model!')
                return

            enzyme_var = self._model.enzyme_variables.get_by_id(enzyme)
            if enzyme_var not in self.enzyme_variables:
                self.enzyme_variables.add(enzyme_var)
            ce_reaction = self._model.reactions.get_by_id('CE_'+self.rxn_id+"_"+enzyme)
            enzyme_var.change_kcat_values({ce_reaction: kcat_dict})
            for direction, kcat in kcats_change.items():
                #change enzyme variable
                enzyme_var.kcats[ce_reaction.id][direction] = kcat
                self._model._change_kcat_in_enzyme_constraint(ce_reaction, enzyme, direction, kcat)
        self.kcats = dict(kcats)

    @staticmethod
    def _extract_reaction_id_from_catalytic_reaction_id(input_str: str,
                                                        protein_id_pattern:str = r'(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})',
                                                        default_enzyme_id_pattern:str = r'E[0-9][0-9]*|Enzyme_*'
                                                        ) -> str:

        """
            Extracts the reaction ID from a catalytic reaction string by removing
            any associated protein IDs in various formats.

            The function identifies and strips off protein IDs from the input string.
            Protein IDs can be in one of two formats:
            1. Standard protein identifiers (by default Uniprot ids: e.g., 'P0ABJ1', 'Q1234X').
            2. Protein IDs for unknown enzymes in the format 'E<number>' (e.g., 'E1', 'E10', 'E534').

            If the string starts with 'CE_', the 'CE_' prefix is removed before processing.

            Args:
                input_str (str): The input string representing a catalytic reaction ID,
                    which may contain associated protein IDs.
                protein_id_pattern (str): Regular expression pattern for matching
                    standard protein IDs. Defaults to a UniProt-based regex.
                default_enzyme_id_pattern (str): Regular expression pattern for matching protein
                    IDs in the format 'E<number>' (e.g., 'E1', 'E99'). Defaults to 'E[0-9][0-9]*'.

            Returns:
                str: The extracted reaction ID with any protein IDs removed.

            Examples:
                >>> _extract_reaction_id_from_catalytic_reaction_id('CE_CYTBO3_4pp_P0ABJ1_P0ABJ5_E123')
                'CYTBO3_4pp'

                >>> _extract_reaction_id_from_catalytic_reaction_id('CE_REACTID_E1_E534')
                'REACTID'

                >>> _extract_reaction_id_from_catalytic_reaction_id('LPLIPAL2ATE140_E3')
                'LPLIPAL2ATE140'

            Notes:
                - The function will strip protein IDs from both the standard format
                  (e.g., 'P0ABJ1', 'Q1234X') and the 'E<number>' format.
                - If no valid protein ID is found, the reaction ID is returned as-is.
                - The function assumes that valid reaction IDs are located at the
                  beginning of the string or after the 'CE_' prefix.
            """

        # Remove the 'CE_' prefix if it exists
        if input_str.startswith('CE_'):
            input_str = input_str[3:]

        # Define the regex pattern to match protein IDs
        protein_id_regex = re.compile(r'_'+protein_id_pattern+ r'|' + r'_' + default_enzyme_id_pattern
)

        # split off all protein ids from the reaction
        reaction_id = protein_id_regex.split(input_str)[0]

        # Remove any trailing or leading underscores that might remain
        reaction_id = reaction_id.strip('_')

        return reaction_id

    def __copy__(self) -> "CatalyticEvent":
        """
        Copy the CatalyticEvent.

        Returns:
            CatalyticEvent:
                A new CatalyticEvent that is a copy of the original CatalyticEvent.
        """

        cop = copy(super(CatalyticEvent, self))
        return cop

    def __deepcopy__(self, memo: dict) -> "CatalyticEvent":
        """
        Copy the CatalyticEvent with memo.

        Parameters:
            memo (dict): Automatically passed parameter.

        Returns:
            CatalyticEvent:
                A new CatalyticEvent that is a copy of the original CatalyticEvent with memo.
        """

        cop = deepcopy(super(CatalyticEvent, self), memo)
        return cop

    def __getstate__(self):
        # Return the state to be pickled
        state = self.__dict__.copy()
        # Handling non-serializable attributes
        state['enzyme_variables'] = list(self.enzyme_variables)
        state['enzymes'] = list(self.enzymes)
        return state

    def __setstate__(self, state):
        # Restore state from the unpickled state
        self.__dict__.update(state)



