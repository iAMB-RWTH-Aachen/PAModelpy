---
sidebar_label: configuration
title: configuration
---

## Config Objects

```python
class Config()
```

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

Defaults are configured for the iML1515 E.coli model

#### P\_TOT\_DEFAULT

g_protein/g_cdw

#### reset

```python
def reset()
```

Reset the config object to the standard settings for E.coli iML1515.

