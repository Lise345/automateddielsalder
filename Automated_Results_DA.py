import os
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

############### Functions to extract values from log files ###############

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
                return None, None, None, None
            if 'SCF Done:  E(RM062X) =' in line:
                match = re.search(r'SCF Done:  E\(RM062X\) =\s+(-?\d+\.\d+)', line)
                if match:
                    last_scf_done = float(match.group(1))
            if '-------------------------------------------------------' in line and i + 1 < len(lines):
                if '# m062x cc-pvtz empiricaldispersion=gd3 Geom=Checkpoint' in lines[i + 1]:
                    pvdz_energy = last_scf_done
                elif '# m062x cc-pvqz empiricaldispersion=gd3 Geom=Checkpoint' in lines[i + 1]:
                    pvtz_energy = last_scf_done
            if 'Thermal correction to Gibbs Free Energy=' in line:
                gibbs_free_energy = float(line.split()[-1])
        pvqz_energy = last_scf_done
    return pvdz_energy, pvtz_energy, pvqz_energy, gibbs_free_energy

def get_separate_reagent_energy(base_id, directory):
    r1_path = next((f for f in os.listdir(directory) if base_id in f and 'R1' in f), None)
    r2_path = next((f for f in os.listdir(directory) if base_id in f and 'R2' in f), None)
    if r1_path and r2_path:
        e1 = extract_values(os.path.join(directory, r1_path))
        e2 = extract_values(os.path.join(directory, r2_path))
        if all(e is not None for e in e1 + e2):
            extrapolated_r1 = (e1[0] * e1[2] - e1[1]**2) / (e1[0] + e1[2] - 2 * e1[1])
            extrapolated_r2 = (e2[0] * e2[2] - e2[1]**2) / (e2[0] + e2[2] - 2 * e2[1])
            return extrapolated_r1 + extrapolated_r2 + e1[3] + e2[3]
    param_path = os.path.join(directory, 'parameters.txt')
    if os.path.exists(param_path):
        with open(param_path, 'r') as f:
            for line in f:
                if line.strip().startswith("Energies of separate reagents"):
                    parts = line.strip().split()
                    if len(parts) >= 5 and parts[4] in base_id:
                        return float(parts[-1])
    return 'Separate reagents not available'

directory = './'
data_dict = {}

for filename in os.listdir(directory):
    if filename.endswith(('Complex.log', 'Product.log', 'SP.log')):
        match = re.search(r'(.*(?:endo|exo))', filename)
        if not match:
            continue
        identification_number = match.group(1)
        file_path = os.path.join(directory, filename)
        pvdz_energy, pvtz_energy, pvqz_energy, gibbs_free_energy = extract_values(file_path)
        if all(v is None for v in [pvdz_energy, pvtz_energy, pvqz_energy, gibbs_free_energy]):
            continue
        if identification_number not in data_dict:
            data_dict[identification_number] = [None] * 14

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

columns = ['ID Number', 'Complex PVDZ Energy', 'Complex PVTZ Energy', 'Complex PVQZ Energy', 'Complex Gibbs Correction',
           'SP PVDZ Energy', 'SP PVTZ Energy', 'SP PVQZ Energy', 'SP Gibbs Correction',
           'Product PVDZ Energy', 'Product PVTZ Energy', 'Product PVQZ Energy', 'Product Gibbs Correction',
           'Energy of Separate Reagents']

df = pd.DataFrame(columns=columns)

for key, val in data_dict.items():
    separate_energy = get_separate_reagent_energy(key, directory)
    val = val[:13] + [separate_energy]
    df.loc[len(df)] = val

mask = df['Energy of Separate Reagents'] != 'Separate reagents not available'

############### Data Cleaning and Calculations ################

# Ensure the original column is numeric (NaN if text like 'Separate reagents not available')
df['Energy of Separate Reagents'] = pd.to_numeric(df['Energy of Separate Reagents'], errors='coerce')

# New column: set to original value if present, or 10 if missing
df['Energy of Separate Reagents'] = df['Energy of Separate Reagents'].fillna(10)

df['Extrapolated Complex Energy'] = (df['Complex PVDZ Energy'] * df['Complex PVQZ Energy'] - df['Complex PVTZ Energy']**2) / (df['Complex PVDZ Energy'] + df['Complex PVQZ Energy'] - 2 * df['Complex PVTZ Energy'])
df['Extrapolated TS Energy'] = (df['SP PVDZ Energy'] * df['SP PVQZ Energy'] - df['SP PVTZ Energy']**2) / (df['SP PVDZ Energy'] + df['SP PVQZ Energy'] - 2 * df['SP PVTZ Energy'])
df['Extrapolated Product Energy'] = (df['Product PVDZ Energy'] * df['Product PVQZ Energy'] - df['Product PVTZ Energy']**2) / (df['Product PVDZ Energy'] + df['Product PVQZ Energy'] - 2 * df['Product PVTZ Energy'])

df['Complex Energy'] = 627.5 * (df['Extrapolated Complex Energy'] + df['Complex Gibbs Correction'] - df['Energy of Separate Reagents'])
df['TS Energy'] = 627.5 * (df['Extrapolated TS Energy'] + df['SP Gibbs Correction'] - df['Energy of Separate Reagents'])
df['Product Energy'] = 627.5 * (df['Extrapolated Product Energy'] + df['Product Gibbs Correction'] - df['Energy of Separate Reagents'])

# Only compute Pi Value where reagent energy is available
df['Pi Value'] = np.where(
    df['Energy of Separate Reagents'].notna(),
    np.exp(-df['Complex Energy'] / (0.001987204259 * 298.15)),
    0
)

df.loc[df['Complex Energy'] > 1, 'Pi Value'] = 0
pi_sum = df['Pi Value'].sum()
df['Percentage'] = df['Pi Value'] / pi_sum

df['TS-C'] = df['TS Energy'] - df['Complex Energy']
df['TS-P'] = df['TS Energy'] - df['Product Energy']
df['P-C'] = df['Product Energy'] - df['Complex Energy']

################ Handling Missing Values ################

df.loc[df['Energy of Separate Reagents'] == 10, 'Energy of Separate Reagents'] = 'Separate energies not available'

################ Plotting the Reaction Pathways ##########

# Extract base ID (e.g., "Rxn1" from "Rxn1_endo")
df['BaseID'] = df['ID Number'].str.extract(r'(.*?)(?:_endo|_exo)', expand=False)
df['Variant'] = df['ID Number'].str.extract(r'_(endo|exo)', expand=False)

# Color mapping
color_map = {'endo': '#24BB7A', 'exo': '#F28500'}

# Plotting
unique_ids = df['BaseID'].dropna().unique()
for base_id in unique_ids:
    subset = df[df['BaseID'] == base_id]
    
    plt.figure(figsize=(8, 6))
    for i, (_, row) in enumerate(subset.iterrows()):
        if isinstance(row['Energy of Separate Reagents'], str):
            continue
        
        variant = row['Variant']
        color = color_map.get(variant, 'gray')
        label = variant if i == 0 or variant not in subset['Variant'][:i].values else None
        
        # X positions: Complex (1), TS (2), Product (3)
        plt.hlines(y=row['Complex Energy'], xmin=0.75, xmax=1.25, colors=color, linewidth=3, label=label)
        plt.hlines(y=row['TS Energy'], xmin=1.75, xmax=2.25, colors=color, linewidth=3)
        plt.hlines(y=row['Product Energy'], xmin=2.75, xmax=3.25, colors=color, linewidth=3)

    plt.hlines(y=0, xmin=0, xmax=4, colors='k', linestyles='dashed', linewidth=1)
    plt.xticks([1, 2, 3], ['Complex', 'TS', 'Product'], fontsize=12)
    plt.title(f"Reaction Path: {base_id}", fontsize=14)
    plt.ylabel('Energy (kcal/mol)', fontsize=12)
    plt.legend(title='Variant', fontsize=11)
    plt.grid(False)
    plt.tight_layout()
    plt.savefig(f'{base_id}_ReactionPath.png')
    plt.close()


################## Making the Bar Plot ##################

# Set the ID column as string if it's not already (for clearer x-axis labels)
df['ID Number'] = df['ID Number'].astype(str)

# Create a new DataFrame for plotting
plot_df = df[['ID Number', 'TS-C', 'TS-P', 'P-C']].set_index('ID Number')

# Transpose for grouped bar plotting (each group per ID)
plot_df = plot_df.T  # Now rows: TS-C, TS-P, P-C ; columns: IDs

# Custom colors for each energy type
colors = ['#065143', "#24BB7A", '#F28500', '#FFBA01',]  # Blue, Orange, Green

# Plot with customizations
ax = plot_df.plot(
    kind='bar',
    figsize=(12, 6),
    color=colors,
    width=0.8
)

# Title and labels with larger font sizes
plt.title('Energy Differences by ID', fontsize=16)
plt.xlabel('Comparison Type', fontsize=14)
plt.ylabel('Energy Difference (kcal/mol)', fontsize=14)

# Tick font sizes
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

# Legend
plt.legend(title='ID', title_fontsize=13, fontsize=11, bbox_to_anchor=(1.05, 1), loc='upper left')

# Grid and layout
plt.grid(True, linestyle='--', alpha=0.7)
plt.tight_layout()
plt.show()

################### Saving the DataFrame to Excel ###################

df.to_excel('AutomatedDA_results.xlsx', index=False)
print("Data extraction complete. The results are saved in 'AutomatedDA_results.xlsx'.")



