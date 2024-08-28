class Config:
    """
    Object with information about model defaults which are used throughout the package:
    - TOTAL_PROTEIN_CONSTRAINT_ID: str,  `TotalProteinConstraint`
    - P_TOT_DEFAULT: float, 0.258 g_p/g_cdw
    - CO2_EXHANGE_RXNID: str, `EX_co2_e`
    - GLUCOSE_EXCHANGE_RXNID: str, `EX_glc__D_e`
    - BIOMASS_REACTION: str, `BIOMASS_Ec_iML1515_core_75p37M`
    - OXYGEN_UPTAKE_RXNID: str, `EX_o2_e`
    - ACETATE_EXCRETION_RXNID: str, `EX_ac_e`
    - PHYS_RXN_IDS: List of str,  `[BIOMASS_REACTION, GLUCOSE_EXCHANGE_RXNID, ACETATE_EXCRETION_RXNID, CO2_EXHANGE_RXNID, OXYGEN_UPTAKE_RXNID,
                        'PGI', 'G6PDH2r', 'EDA', 'CS', 'ICL', 'PPC', 'ME1', 'ME2']`
    - ENZYME_ID_REGEX: r-str: r'(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})'

    Defaults are configured for the iML1515 E.coli model
    """

    TOTAL_PROTEIN_CONSTRAINT_ID = "TotalProteinConstraint"
    P_TOT_DEFAULT = 0.258  # g_protein/g_cdw
    CO2_EXHANGE_RXNID = "EX_co2_e"
    GLUCOSE_EXCHANGE_RXNID = "EX_glc__D_e"
    BIOMASS_REACTION = "BIOMASS_Ec_iML1515_core_75p37M"
    OXYGEN_UPTAKE_RXNID = "EX_o2_e"
    ACETATE_EXCRETION_RXNID = "EX_ac_e"
    PHYS_RXN_IDS = [
        BIOMASS_REACTION,
        GLUCOSE_EXCHANGE_RXNID,
        ACETATE_EXCRETION_RXNID,
        CO2_EXHANGE_RXNID,
        OXYGEN_UPTAKE_RXNID,
        "PGI",
        "G6PDH2r",
        "EDA",
        "CS",
        "ICL",
        "PPC",
        "ME1",
        "ME2",
    ]

    # Define the regex pattern for protein IDs, obtained from UniProtKB, 2024-08-07
    # https://www.uniprot.org/help/accession_numbers
    ENZYME_ID_REGEX = r'(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})'

    def reset(self):
        """
        Reset the config object to the standard settings for E.coli iML1515.
        """

        self.TOTAL_PROTEIN_CONSTRAINT_ID = "TotalProteinConstraint"
        self.P_TOT_DEFAULT = 0.258  # g_protein/g_cdw
        self.CO2_EXHANGE_RXNID = "EX_co2_e"
        self.GLUCOSE_EXCHANGE_RXNID = "EX_glc__D_e"
        self.BIOMASS_REACTION = "BIOMASS_Ec_iML1515_core_75p37M"
        self.OXYGEN_UPTAKE_RXNID = "EX_o2_e"
        self.ACETATE_EXCRETION_RXNID = "EX_ac_e"
        self.PHYS_RXN_IDS = [
            self.BIOMASS_REACTION,
            self.GLUCOSE_EXCHANGE_RXNID,
            self.ACETATE_EXCRETION_RXNID,
            self.CO2_EXHANGE_RXNID,
            self.OXYGEN_UPTAKE_RXNID,
            "PGI",
            "G6PDH2r",
            "EDA",
            "CS",
            "ICL",
            "PPC",
            "ME1",
            "ME2",
        ]
        self.ENZYME_ID_REGEX = r'(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})'