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
    enthalpy = None
    last_scf_done = None

    with open(file_path, 'r', encoding='utf-8', errors='replace') as file:
        lines = file.readlines()
        for i, line in enumerate(lines):
            if 'Error termination via Lnk1e' in line:
                return None, None, None, None, None
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
            if 'Thermal correction to Enthalpy=' in line:
                enthalpy = float(line.split()[-1])
        pvqz_energy = last_scf_done
    return pvdz_energy, pvtz_energy, pvqz_energy, gibbs_free_energy, enthalpy

def get_separate_reagent_energy(base_id, directory):
    r1_path = next((f for f in os.listdir(directory) if base_id in f and 'R1' in f), None)
    r2_path = next((f for f in os.listdir(directory) if base_id in f and 'R2' in f), None)
    if r1_path and r2_path:
        try:
            e1 = extract_values(os.path.join(directory, r1_path))
        except UnicodeDecodeError as e:
            print(f"Failed to read file: {r1_path} with error: {e}")
            return 'Separate reagents not available', 'Separate reagents not available'
        e2 = extract_values(os.path.join(directory, r2_path))
        if all(e is not None for e in e1 + e2):
            extrapolated_r1 = (e1[0] * e1[2] - e1[1]**2) / (e1[0] + e1[2] - 2 * e1[1])
            extrapolated_r2 = (e2[0] * e2[2] - e2[1]**2) / (e2[0] + e2[2] - 2 * e2[1])
            gibss_separate_energies = extrapolated_r1 + extrapolated_r2 + e1[3] + e2[3]
            enth_separate_energies = extrapolated_r1 + extrapolated_r2 + e1[4] + e2[4]
            return gibss_separate_energies, enth_separate_energies
    return 'Separate reagents not available', 'Separate reagents not available'

directory = './'
data_dict = {}

matched = 0
parsed = 0
skipped_all_none = 0

for filename in os.listdir(directory):
    if filename.endswith(('Complex.log', 'Product.log', 'SP.log')):
        matched += 1

        name = os.path.splitext(filename)[0]
        identification_number = name.split('-TS', 1)[0]

        file_path = os.path.join(directory, filename)
        vals = extract_values(file_path)
        pvdz_energy, pvtz_energy, pvqz_energy, gibbs_free_energy, enthalpy = vals

        if all(v is None for v in vals):
            skipped_all_none += 1
            continue

        parsed += 1
        
        if identification_number not in data_dict:
            data_dict[identification_number] = [None] * 16

        if filename.endswith('Complex.log'):
            data_dict[identification_number][0:6] = [identification_number, pvdz_energy, pvtz_energy, pvqz_energy, gibbs_free_energy, enthalpy]
        elif filename.endswith('SP.log'):
            data_dict[identification_number][6:11] = [pvdz_energy, pvtz_energy, pvqz_energy, gibbs_free_energy, enthalpy]
        elif filename.endswith('Product.log'):
            data_dict[identification_number][11:16] = [pvdz_energy, pvtz_energy, pvqz_energy, gibbs_free_energy, enthalpy]
            
print("matched log files:", matched)
print("parsed (at least one value found):", parsed)
print("skipped (all None):", skipped_all_none)

columns = ['ID Number', 'Complex PVDZ Energy', 'Complex PVTZ Energy', 'Complex PVQZ Energy', 'Complex Gibbs Correction', 'Complex Enth Correction',
           'SP PVDZ Energy', 'SP PVTZ Energy', 'SP PVQZ Energy', 'SP Gibbs Correction', 'SP Enth Correction',
           'Product PVDZ Energy', 'Product PVTZ Energy', 'Product PVQZ Energy', 'Product Gibbs Correction', 'Product Enth Correction',
           'Gibbs of Separate Reagents', 'Enthalpy of Separate Reagents']

df = pd.DataFrame(columns=columns)

for key, val in data_dict.items():
    separate_energy, enth_separate = get_separate_reagent_energy(key, directory)
    val = val[:16] + [separate_energy] + [enth_separate]
    df.loc[len(df)] = val

df['Gibbs of Separate Reagents'] = pd.to_numeric(df['Gibbs of Separate Reagents'], errors='coerce').fillna(10)
df['Enthalpy of Separate Reagents'] = pd.to_numeric(df['Enthalpy of Separate Reagents'], errors='coerce').fillna(10)


df['Extrapolated Complex Energy'] = (df['Complex PVDZ Energy'] * df['Complex PVQZ Energy'] - df['Complex PVTZ Energy']**2) / (df['Complex PVDZ Energy'] + df['Complex PVQZ Energy'] - 2 * df['Complex PVTZ Energy'])
df['Extrapolated TS Energy'] = (df['SP PVDZ Energy'] * df['SP PVQZ Energy'] - df['SP PVTZ Energy']**2) / (df['SP PVDZ Energy'] + df['SP PVQZ Energy'] - 2 * df['SP PVTZ Energy'])
df['Extrapolated Product Energy'] = (df['Product PVDZ Energy'] * df['Product PVQZ Energy'] - df['Product PVTZ Energy']**2) / (df['Product PVDZ Energy'] + df['Product PVQZ Energy'] - 2 * df['Product PVTZ Energy'])

df['Complex Gibbs'] = 627.5 * (df['Extrapolated Complex Energy'] + df['Complex Gibbs Correction'] - df['Gibbs of Separate Reagents'])
df['TS Gibbs'] = 627.5 * (df['Extrapolated TS Energy'] + df['SP Gibbs Correction'] - df['Gibbs of Separate Reagents'])
df['Product Gibbs'] = 627.5 * (df['Extrapolated Product Energy'] + df['Product Gibbs Correction'] - df['Gibbs of Separate Reagents'])

df['Complex Enthalpy'] = 627.5 * (df['Extrapolated Complex Energy'] + df['Complex Enth Correction'] - df['Enthalpy of Separate Reagents'])
df['TS Enthalpy'] = 627.5 * (df['Extrapolated TS Energy'] + df['SP Enth Correction'] - df['Enthalpy of Separate Reagents'])
df['Product Enthalpy'] = 627.5 * (df['Extrapolated Product Energy'] + df['Product Enth Correction'] - df['Enthalpy of Separate Reagents'])

df['TS-C Enthalpy'] = df['TS Enthalpy'] - df['Complex Enthalpy']
df['TS-P Enthalpy'] = df['TS Enthalpy'] - df['Product Enthalpy']
df['P-C Enthalpy'] = df['Product Enthalpy'] - df['Complex Enthalpy']

df['TS-C'] = df['TS Gibbs'] - df['Complex Gibbs']
df['TS-P'] = df['TS Gibbs'] - df['Product Gibbs']
df['P-C'] = df['Product Gibbs'] - df['Complex Gibbs']

df['Gibbs of Separate Reagents'] = df['Gibbs of Separate Reagents'].astype(object)
df['Enthalpy of Separate Reagents'] = df['Enthalpy of Separate Reagents'].astype(object)


df.loc[df['Gibbs of Separate Reagents'] == 10, 'Gibbs of Separate Reagents'] = 'No separate energies available'
df.loc[df['Enthalpy of Separate Reagents'] == 10, 'Enthalpy of Separate Reagents'] = 'No separate energies available'


# Plotting reaction path and bar plots
df['BaseID'] = df['ID Number'].str.extract(r'(.*?)(?:_endo|_exo)', expand=False)
df['Variant'] = df['ID Number'].str.extract(r'_(endo|exo)', expand=False)
df['ID Number'] = df['ID Number'].astype(str)

color_map = {'endo': '#24BB7A', 'exo': '#F28500'}
os.makedirs("plots", exist_ok=True)

from matplotlib.lines import Line2D

for base_id in df['BaseID'].dropna().unique():
    subset = df[df['BaseID'] == base_id]
    for energy_type in ['Gibbs', 'Enthalpy']:
        plt.figure(figsize=(8, 6))
        used_variants = set()

        for _, row in subset.iterrows():
            if isinstance(row[f'{energy_type} of Separate Reagents'], str):
                continue
            y_vals = [row[f'Complex {energy_type}'], row[f'TS {energy_type}'], row[f'Product {energy_type}']]
            x_vals = [1, 2, 3]
            variant = row['Variant']
            color = color_map.get(variant, 'gray')
            bar_width = 0.2

            # Draw horizontal bars for each point
            for x, y in zip(x_vals, y_vals):
                plt.hlines(y, x - bar_width, x + bar_width, color=color, linewidth=3)

            # Add dummy line to legend only once per variant
            if variant not in used_variants:
                plt.plot([], [], color=color, linewidth=3, label=variant)
                used_variants.add(variant)

        plt.xticks([1, 2, 3], ['Complex', 'TS', 'Product'], fontsize=12)
        plt.title(f"Reaction Path ({energy_type}): {base_id}", fontsize=14)
        plt.ylabel(f'{energy_type} (kcal/mol)', fontsize=12)
        plt.legend(title='Variant', fontsize=11)
        plt.grid(True, linestyle='--', alpha=0.7)
        plt.tight_layout()
        plt.savefig(f'plots/{base_id}_ReactionPath_{energy_type}.png')
        plt.close()


# Bar plots
for energy_type in ['', ' Enthalpy']:
    cols = [f'TS-C{energy_type}', f'TS-P{energy_type}', f'P-C{energy_type}']
    plot_df = df[['ID Number'] + cols].set_index('ID Number').T
    plot_df.plot(kind='bar', figsize=(12, 6), width=0.8)
    plt.title(f'{energy_type.strip()} Differences by ID', fontsize=16)
    plt.xlabel('Comparison Type', fontsize=14)
    plt.ylabel(f'{energy_type.strip()} Difference (kcal/mol)', fontsize=14)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.legend(title='ID', title_fontsize=13, fontsize=11, bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.savefig(f'plots/{energy_type.strip()}_BarPlot.png')

# Create Results subset
results_gibbs = ['ID Number', 'Complex Gibbs', 'TS Gibbs', 'Product Gibbs', 'TS-C', 'TS-P', 'P-C']
resultsgibbs_df = df[results_gibbs]
results_enthalpy = ['ID Number', 'Complex Enthalpy', 'TS Enthalpy', 'Product Enthalpy', 'TS-C Enthalpy', 'TS-P Enthalpy', 'P-C Enthalpy']
resultsenth_df = df[results_enthalpy]


# Save both DataFrames to different sheets
with pd.ExcelWriter('AutomatedDA_results.xlsx', engine='openpyxl') as writer:
    df.to_excel(writer, sheet_name='Data', index=False)
    resultsgibbs_df.to_excel(writer, sheet_name='Results Gibbs', index=False)
    resultsenth_df.to_excel(writer, sheet_name='Results Enthalpy', index=False)
