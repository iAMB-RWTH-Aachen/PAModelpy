import math
from warnings import warn
from copy import copy, deepcopy
from cobra import Object

from .configuration import Config
from .EnzymeSectors import EnzymeSector


class MembraneSector(EnzymeSector):
    def __init__(
            self,
            area_avail_0, #μm2
            area_avail_mu, #μm2/h
            alpha_numbers_dict: {},
            enzyme_location: {},
            cog_class: {} = None,
            max_area: float = 0.033,
            configuration=Config):

        self.id = 'MembraneSector'
        self.area_avail_0 = area_avail_0
        self.area_avail_mu = area_avail_mu
        self.alpha_numbers_dict = alpha_numbers_dict
        self.cog_class = cog_class
        self.enzyme_location = enzyme_location
        self.area_alpha = math.pi * math.pow((0.00023), 2)  # area per alpha helix unit [um]
        self.cdw_per_cell = 0.28 * 1e-12 #0.28 pg
        self.n_a = 6.02214076 * 1e23 #avogadro number
        self.max_membrane_area = max_area #percentage of membrane area that can be covered by proteins
        self.unit_factor = 1e-3 * self.cdw_per_cell * self.n_a

        #Defining the slope and intercept
        self.intercept = self.area_avail_0 #μm2
        self.slope = self.area_avail_mu #μm2/h

    def add(self, model):

        print("Add membrane protein sector \n")
        model.membrane_sector = self
        self._add_membrane_constraint(model)
        pass

    def _add_membrane_constraint(self, model):

        self.membrane_proteins = {}

        self.total_occupied_membrane = 0
        coefficients = {
            model.reactions.get_by_id(model.BIOMASS_REACTION).forward_variable: -self.slope
        }

        for enz_complex in model.enzyme_variables:
            enzymes = enz_complex.id.split("_")
            alpha_numbers_in_complex = [0]  # zero if enzyme is not in membrane

            for enz in enzymes:
                if enz in self.alpha_numbers_dict.keys() and self.enzyme_location[
                    enz] == 'Cell membrane':
                    alpha_numbers_in_complex.append(self.alpha_numbers_dict[enz])
                    self.membrane_proteins[enz_complex.id] = enz_complex.kcats

            alpha_numbers_for_complex = max(alpha_numbers_in_complex)

            coefficients[enz_complex.forward_variable] = (
                        1e-6 * alpha_numbers_for_complex  # correction for the solver issue
                        * self.area_alpha
                        * self.unit_factor
                        / self.max_membrane_area)

            coefficients[enz_complex.reverse_variable] = (
                        1e-6 * alpha_numbers_for_complex  # correction for the solver issue
                        * self.area_alpha
                        * self.unit_factor
                        / self.max_membrane_area)

        occupied_membrane = model.problem.Constraint(0, lb=0, ub=self.intercept, name='membrane')
        model.add_cons_vars(occupied_membrane)
        model.solver.update()
        occupied_membrane.set_linear_coefficients(coefficients=coefficients)

        # # Debugging
        # for rxn, lb in model.rxn_old_bounds_lb.items():
        #     model.reactions.get_by_id(rxn).lower_bound = lb
        #
        # for rxn, ub in model.rxn_old_bounds_ub.items():
        #     model.reactions.get_by_id(rxn).upper_bound = ub

    def calculate_occupied_membrane(self, model):
        occupied_membrane = 0

        for enz_complex in model.enzyme_variables:
            enzymes = enz_complex.id.split("_")
            enz_complex_concentration = enz_complex.forward_variable.primal + enz_complex.reverse_variable.primal
            alpha_numbers_in_complex = [0]

            for enz in enzymes:
                if enz in self.alpha_numbers_dict.keys() and self.enzyme_location[
                    enz] == 'Cell membrane':
                    alpha_numbers_in_complex.append(self.alpha_numbers_dict[enz])

            alpha_numbers_for_complex = max(alpha_numbers_in_complex)
            occupied_membrane += (
                        enz_complex_concentration * 1e-6 * alpha_numbers_for_complex  # correction for the solver issue
                        * self.area_alpha
                        * self.unit_factor)

        available_membrane = self.slope * model.objective.value + self.intercept

        return occupied_membrane, available_membrane

    def change_available_membrane_area(self, new_max_area: float, model):
        self._update_membrane_constraint(new_max_area, model)

    def _update_membrane_constraint(self, new_max_area:float, model):
        self.max_membrane_area = new_max_area
        self.membrane_proteins = {}

        self.total_occupied_membrane = 0
        coefficients = {
            model.reactions.get_by_id(model.BIOMASS_REACTION).forward_variable: -self.slope
        }

        for enz_complex in model.enzyme_variables:
            enzymes = enz_complex.id.split("_")
            alpha_numbers_in_complex = [0]  # zero if enzyme is not in membrane

            for enz in enzymes:
                if enz in self.alpha_numbers_dict.keys() and self.enzyme_location[enz] == 'Cell membrane':
                    alpha_numbers_in_complex.append(self.alpha_numbers_dict[enz])
                    self.membrane_proteins[enz_complex.id] = enz_complex.kcats

            alpha_numbers_for_complex = max(alpha_numbers_in_complex)

            coefficients[enz_complex.forward_variable] = (
                        1e-6 * alpha_numbers_for_complex  # correction for the solver issue
                        * self.area_alpha
                        * self.unit_factor
                        / new_max_area)

            coefficients[enz_complex.reverse_variable] = (
                        1e-6 * alpha_numbers_for_complex  # correction for the solver issue
                        * self.area_alpha
                        * self.unit_factor
                        / new_max_area)

        model.constraints['membrane'].set_linear_coefficients(coefficients=coefficients)
        model.solver.update()

        return self