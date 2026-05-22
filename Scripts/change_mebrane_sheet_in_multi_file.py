import pandas as pd

df = pd.DataFrame({
    "Parameter": ["id_list", "sv_0", "sv_slope"],
    "Value": ["membrane", 3.2042, -0.3285],
    "Unit": ["", "µm2/fL", "µm²·h/fL"],
    "Description": [
        "membrane",
        "surface to volume ratio",
        "Increase in ..."
    ]
})


date_list = ['2026_04_01', '2026_05_06', '2026_05_08', '2026_05_09', '2026_05_14']

for date in date_list:
    file_path = f'Results/PAM_parametrizer/Diagnostics_files/{date}/proteinAllocationModel_iML1515_EnzymaticData_multi.xlsx'

    with pd.ExcelWriter(path=file_path,
                        engine='openpyxl',
                        mode='a',
                        if_sheet_exists='replace') as writer:
        df.to_excel(writer, sheet_name='Membrane', index=False)