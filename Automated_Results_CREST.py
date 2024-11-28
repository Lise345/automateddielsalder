import os
import re
import pandas as pd
import numpy as np

def extract_values(file_path):
    pvdz_energy = None
    pvtz_energy = None
    pvqz_energy = None
    gibbs_free_energy = None

    last_scf_done = None

    with open(file_path, 'r') as file:
        lines = file.readlines()
        for i, line in enumerate(lines):
            if 'Error termination via Lnk1e' in line:
                return None, None, None, None  # Discard the file if it contains the error message
            if 'SCF Done:  E(RM062X) =' in line:
                match = re.search(r'SCF Done:  E\(RM062X\) =\s+(-?\d+\.\d+)', line)
                if match:
                    last_scf_done = float(match.group(1))
            if '-------------------------------------------------------' in line and i + 1 < len(lines) and '# m062x cc-pvtz empiricaldispersion=gd3 Geom=Checkpoint' in lines[i + 1]:
                pvdz_energy = last_scf_done
            if '-------------------------------------------------------' in line and i + 1 < len(lines) and '# m062x cc-pvqz empiricaldispersion=gd3 Geom=Checkpoint' in lines[i + 1]:
                pvtz_energy = last_scf_done
            if 'Thermal correction to Gibbs Free Energy=' in line:
                gibbs_free_energy = float(line.split()[-1])

        # Assign the last SCF Done value to pvqz_energy
        pvqz_energy = last_scf_done

    return pvdz_energy, pvtz_energy, pvqz_energy, gibbs_free_energy

directory = './'
data_dict = {}

for filename in os.listdir(directory):
    if filename.startswith('ccpvdz_startgeom-') and (filename.endswith('Complex.log') or filename.endswith('Product.log') or filename.endswith('SP.log')):
        identification_number = filename.split('-')[1].split('_')[0]
        file_path = os.path.join(directory, filename)
        pvdz_energy, pvtz_energy, pvqz_energy, gibbs_free_energy = extract_values(file_path)
        
        if pvdz_energy is None and pvtz_energy is None and pvqz_energy is None and gibbs_free_energy is None:
            continue  # Skip files that contain the error message
        
        if identification_number not in data_dict:
            data_dict[identification_number] = [None] * 13
        
        if filename.endswith('Complex.log'):
            data_dict[identification_number][0] = identification_number
            data_dict[identification_number][1] = pvdz_energy
            data_dict[identification_number][2] = pvtz_energy
            data_dict[identification_number][3] = pvqz_energy
            data_dict[identification_number][4] = gibbs_free_energy
        elif filename.endswith('SP.log'):
            data_dict[identification_number][5] = pvdz_energy
            data_dict[identification_number][6] = pvtz_energy
            data_dict[identification_number][7] = pvqz_energy
            data_dict[identification_number][8] = gibbs_free_energy
        elif filename.endswith('Product.log'):
            data_dict[identification_number][9] = pvdz_energy
            data_dict[identification_number][10] = pvtz_energy
            data_dict[identification_number][11] = pvqz_energy
            data_dict[identification_number][12] = gibbs_free_energy

data_list = [value for key, value in data_dict.items()]

df = pd.DataFrame(data_list, columns=['ID Number', 'Complex PVDZ Energy', 'Complex PVTZ Energy', 'Complex PVQZ Energy', 'Complex Gibbs Correction',
                                      'SP PVDZ Energy', 'SP PVTZ Energy', 'SP PVQZ Energy', 'SP Gibbs Correction',
                                      'Product PVDZ Energy', 'Product PVTZ Energy', 'Product PVQZ Energy', 'Product Gibbs Correction'])

# Add the new columns
energy_of_separate_reagents = -936.284718313083

df['Extrapolated Complex Energy'] = (df['Complex PVDZ Energy'] * df['Complex PVQZ Energy'] - df['Complex PVTZ Energy']**2) / (df['Complex PVDZ Energy'] + df['Complex PVQZ Energy'] - 2 * df['Complex PVTZ Energy'])
df['Extrapolated TS Energy'] = (df['SP PVDZ Energy'] * df['SP PVQZ Energy'] - df['SP PVTZ Energy']**2) / (df['SP PVDZ Energy'] + df['SP PVQZ Energy'] - 2 * df['SP PVTZ Energy'])
df['Extrapolated Product Energy'] = (df['Product PVDZ Energy'] * df['Product PVQZ Energy'] - df['Product PVTZ Energy']**2) / (df['Product PVDZ Energy'] + df['Product PVQZ Energy'] - 2 * df['Product PVTZ Energy'])

df['Complex Energy'] = 627.5*(df['Extrapolated Complex Energy'] + df['Complex Gibbs Correction'] - energy_of_separate_reagents)
df['TS Energy'] = 627.5*(df['Extrapolated TS Energy'] + df['SP Gibbs Correction'] - energy_of_separate_reagents)
df['Product Energy'] = 627.5*(df['Extrapolated Product Energy'] + df['Product Gibbs Correction'] - energy_of_separate_reagents)

df['Pi Value'] = np.exp(-df['Complex Energy'] / (0.001987204259 * 298.15))
df.loc[df['Complex Energy'] > 1, 'Pi Value'] = 0  # Set Pi Value to 0 if Complex Energy is larger than 1

pi_sum = df['Pi Value'].sum()
df['Percentage'] = df['Pi Value'] / pi_sum

df['Rate Constant'] = ((298.15 * 1.380649E-23) / 6.62607015E-34) * np.exp(-(df['TS Energy'] - df['Complex Energy']) * 1000 * 4.184 / (8.314 * 298.15))

df.to_excel('FASTCAR_results.xlsx', index=False)

print("Data extraction complete. The results are saved in 'FASTCAR_results.xlsx'.")
