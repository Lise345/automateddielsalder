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
    files = os.listdir(directory)

    r1_candidates = [f for f in files if base_id in f and 'R1' in f and f.endswith('.log')]
    r2_candidates = [f for f in files if base_id in f and 'R2' in f and f.endswith('.log')]

    print(f"[DEBUG] {base_id} R1 candidates:", r1_candidates)
    print(f"[DEBUG] {base_id} R2 candidates:", r2_candidates)

    r1_path = r1_candidates[0] if r1_candidates else None
    r2_path = r2_candidates[0] if r2_candidates else None

    if not (r1_path and r2_path):
        print(f"[DEBUG] {base_id} missing R1 or R2 selection.")
        return 'Separate reagents not available', 'Separate reagents not available'

    e1 = extract_values(os.path.join(directory, r1_path))
    e2 = extract_values(os.path.join(directory, r2_path))

    print(f"[DEBUG] {base_id} e1:", e1)
    print(f"[DEBUG] {base_id} e2:", e2)

    if not all(e is not None for e in e1 + e2):
        print(f"[DEBUG] {base_id} extraction has None(s).")
        return 'Separate reagents not available', 'Separate reagents not available'

    extrapolated_r1 = (e1[0] * e1[2] - e1[1]**2) / (e1[0] + e1[2] - 2 * e1[1])
    extrapolated_r2 = (e2[0] * e2[2] - e2[1]**2) / (e2[0] + e2[2] - 2 * e2[1])

    gibbs_sep = extrapolated_r1 + extrapolated_r2 + e1[3] + e2[3]
    enth_sep  = extrapolated_r1 + extrapolated_r2 + e1[4] + e2[4]
    return gibbs_sep, enth_sep


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

print("Columns containing Gibbs:", [c for c in df.columns if "Gibbs" in c])
print("Columns containing Enthalpy:", [c for c in df.columns if "Enthalpy" in c])
print("BaseID non-NaN count:", df["BaseID"].notna().sum() if "BaseID" in df.columns else "NO BaseID COL")

# Check if the expected columns exist
for energy_type in ["Gibbs", "Enthalpy"]:
    needed = [f"Complex {energy_type}", f"TS {energy_type}", f"Product {energy_type}"]
    print("Missing for", energy_type, ":", [c for c in needed if c not in df.columns])




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
    
for energy_type in ['Gibbs', 'Enthalpy']:
    plt.figure(figsize=(10, 7))
    used_variants = set()

    # Use unique BaseIDs to define deterministic offsets (so lines don't sit on top of each other)
    base_ids = [b for b in df['ID Number'].dropna().unique()]
    n = len(base_ids)
    if n == 0:
        continue

    # Spread offsets symmetrically around 0; scale keeps things visually readable
    offsets = np.linspace(-0.18, 0.18, n) if n > 1 else [0.0]
    offset_map = {b: offsets[i] for i, b in enumerate(base_ids)}

    bar_width = 0.18  # half-width of the horizontal "state" bars

    # Loop through all rows (all reactions)
    for _, row in df.iterrows():
        base_id = row.get('ID Number', None)
        if pd.isna(base_id):
            continue

        # Skip if separate reagents field is a string (your existing guard)
        if isinstance(row.get(f'{energy_type} of Separate Reagents', None), str):
            continue

        # Make sure the three energies exist and are numeric
        y_vals = [row.get(f'Complex {energy_type}', None),
                  row.get(f'TS {energy_type}', None),
                  row.get(f'Product {energy_type}', None)]
        if any(v is None or (isinstance(v, float) and np.isnan(v)) for v in y_vals):
            continue

        # x positions with BaseID-specific offset
        x_base = np.array([1.0, 2.0, 3.0])
        x_vals = x_base + offset_map.get(base_id, 0.0)

        variant = row.get('Variant', 'Unknown')
        color = color_map.get(variant, 'gray')

        # 1) dotted connector line (requested)
        plt.plot(x_vals, y_vals, linestyle=':', linewidth=1.5, color=color, alpha=0.85)

        # 2) horizontal bars for Complex/TS/Product
        for x, y in zip(x_vals, y_vals):
            plt.hlines(y, x - bar_width, x + bar_width, color=color, linewidth=3, alpha=0.95)

        # Legend entry once per variant
        if variant not in used_variants:
            plt.plot([], [], color=color, linewidth=3, label=variant)
            used_variants.add(variant)

    plt.xticks([1, 2, 3], ['Complex', 'TS', 'Product'], fontsize=12)
    plt.title(f"Reaction Path ({energy_type}) — All Reactions", fontsize=14)
    plt.ylabel(f'{energy_type} (kcal/mol)', fontsize=12)
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.legend(title='Variant', fontsize=10)
    plt.tight_layout()
    plt.savefig(f'plots/ALL_ReactionPath_{energy_type}.png', dpi=300)
    plt.close()

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
