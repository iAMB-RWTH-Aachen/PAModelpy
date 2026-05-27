import math

class Protein:
    """
    An object which holds functions to calculate the area occupied by this protein species
    by converting the concentration of this protein species into the number of protein molecules
    and then sum the area occupied by the protein molecules

    Args:
        mu(float): the cell's growth rate
        total_protein_per_volume(float) = total protein concentration that the cell can host

    """
    def __init__(self, mu:float=0.67,
                 total_protein_per_volume:float=200):
        self.mu = mu #1/h
        self.total_protein_per_volume = total_protein_per_volume #mg/mL
        self.avogadro = 6.02 * 1e23
        self.alpha_helix_width = 0.46 * 1e-9 #m
        self.alpha_area = (math.pi * math.pow((self.alpha_helix_width / 2), 2))

    def calculate_cell_volume(self):
        self.cell_volume = -0.53*math.pow(self.mu, 2) + 2.41*self.mu + 1.62 #fL or um^3

        return self.cell_volume
    def calculate_protein_disc_area(self, alpha_helix_units):
        self.protein_disc_area = alpha_helix_units * self.alpha_area
        # unit: protein_disc_area [m2], number_of_alpha_helix [-], alpha_helix_width [m]
        return self.protein_disc_area
    def calculate_protein_number_per_volume(self, protein_concentration, p_conc_unit):
        if p_conc_unit == 'fmol/ug':
            protein_concentration = protein_concentration * self.total_protein_per_volume * 1e-9
            #unit: protein_concentration (per volume) [fmol/fL], protein_concentration [fmol/ug], tot_protein_per_volume [ug/ul]  
        elif p_conc_unit == 'fmol/cell':
            protein_concentration = protein_concentration / self.cell_volume
            #unit: protein_concentration (fmol/fL), protein_concentration_per_volume [fmol], cell_volume [fL]
        else:
            ValueError('Please provide provide concentration unit in either fmol/ug or fmol/cell')

        'calculate protein number per cell'
        protein_nr_per_volume = protein_concentration * 1e-15 * self.avogadro
        #unit: protein_nr_per_cell [-/fL], protein_concentration_per_cell [fmol/fL], avogadro [1/mol]'

        return protein_nr_per_volume

    def calculate_protein_area(self, protein_concentration, alpha_helix_units, p_conc_unit:str="fmol/ug"):
        self.calculate_cell_volume()
        self.calculate_protein_disc_area(alpha_helix_units)
        protein_nr_per_volume = self.calculate_protein_number_per_volume(protein_concentration, p_conc_unit) if p_conc_unit != "protein copies/cell" else protein_concentration
        self.S_V_occupied = protein_nr_per_volume * self.protein_disc_area * 1e12
        #unit: S_V_occupied [um2/fL], protein_nr_per_cell [-/fL], protein_disc_area [m2]'
        return self.S_V_occupied

    def is_in_membrane(self, type):
        return type == 'membrane protein'

    def convert_unit(self):
        '''Convert unit from [g protein/L cell volume] to [g protein/g cdw]'''
        total_protein_per_cdw = self.total_protein_per_volume * 0.000855

        return total_protein_per_cdw




