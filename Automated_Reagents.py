import os
import re
import subprocess
import numpy as np

# Atomic symbols mapping
atomic_symbols = {
    1: "H", 6: "C", 7: "N", 8: "O", 9: "F", 16: "S", 17: "Cl"  # Add other atomic symbols if needed
}

# Covalent radii in Å
radii = {
    'H': 1.20, 'He': 1.40, 'Li': 1.82, 'Be': 1.53, 'B': 1.92,
    'C': 1.70, 'N': 1.55, 'O': 1.52, 'F': 1.47, 'Ne': 1.54,
    'Na': 2.27, 'Mg': 1.73, 'Al': 1.84, 'Si': 2.10, 'P': 1.80,
    'S': 2.10, 'Cl': 1.75, 'Ar': 1.88, 'K': 2.75, 'Ca': 2.31,
    'Sc': 2.30, 'Ti': 2.15, 'V': 2.05, 'Cr': 2.05, 'Mn': 2.05,
    'Fe': 2.00, 'Co': 2.00, 'Ni': 1.97, 'Cu': 1.96, 'Zn': 2.01,
    'Ga': 1.87, 'Ge': 2.11, 'As': 1.85, 'Se': 1.90, 'Br': 1.85,
    'Kr': 2.02, 'Rb': 3.03, 'Sr': 2.49, 'Y': 2.40, 'Zr': 2.30,
    'Nb': 2.15, 'Mo': 2.10, 'Tc': 2.05, 'Ru': 2.05, 'Rh': 2.00,
    'Pd': 2.05, 'Ag': 2.03, 'Cd': 2.18, 'In': 1.93, 'Sn': 2.17,
    'Sb': 2.06, 'Te': 2.06, 'I': 1.98, 'Xe': 2.16, 'Cs': 3.43,
    'Ba': 2.68, 'La': 2.50, 'Ce': 2.48, 'Pr': 2.47, 'Nd': 2.45,
    'Pm': 2.43, 'Sm': 2.42, 'Eu': 2.40, 'Gd': 2.38, 'Tb': 2.37,
    'Dy': 2.35, 'Ho': 2.33, 'Er': 2.32, 'Tm': 2.30, 'Yb': 2.28,
    'Lu': 2.27, 'Hf': 2.25, 'Ta': 2.20, 'W': 2.10, 'Re': 2.05,
    'Os': 2.00, 'Ir': 2.00, 'Pt': 2.05, 'Au': 2.10, 'Hg': 2.05,
    'Tl': 1.96, 'Pb': 2.02, 'Bi': 2.07, 'Po': 1.97, 'At': 2.02,
    'Rn': 2.20, 'Fr': 3.48, 'Ra': 2.83, 'Ac': 2.60, 'Th': 2.37,
    'Pa': 2.43, 'U': 2.40, 'Np': 2.39, 'Pu': 2.43, 'Am': 2.44,
    'Cm': 2.45, 'Bk': 2.44, 'Cf': 2.45, 'Es': 2.45, 'Fm': 2.45,
    'Md': 2.45, 'No': 2.45, 'Lr': 2.45
}

def read_parameters(file_path):
    """Reads molecule atom indices from parameters.txt."""
    with open(file_path, 'r') as file:
        file_content = file.read()

    rootdir_match = re.search(r'rootdir\s+(.+)', file_content)
    rootdir = rootdir_match.group(1).strip().strip("'\"")
    
    time_sr = re.search(r'Time for Separate Reagent calcs\s+(\d+)', file_content)
    if time_sr:
        time_for_sr_calcs = int(time_sr.group(1))
        sr_time = f'{time_for_sr_calcs}:00:00'  # Format to HH:MM:SS
    else:
        sr_time = '25:00:00'  # Default value if not found
        print("Time for stationary calculations not found, defaulting to 25:00:00")
        
    functional = re.search(r'Functional\s+(\S+)', file_content, re.IGNORECASE).group(1)
    dispersion = re.search(r'Dispersion\s+(\S+)', file_content, re.IGNORECASE).group(1)
    basis_raw = re.search(r'Basis\s+(\S+)', file_content, re.IGNORECASE).group(1)
    solvent = re.search(r'DFT solvent\s+(\S+)', file_content, re.IGNORECASE).group(1)

    # Translate 'cbs' to cc-pvdz as base level
    basis = "cc-pvdz" if basis_raw.lower() == "cbs" else basis_raw.lower()

        
    # Extract Gaussian module
    gaussian_module_match = re.search(r'Gaussian module (.+)', file_content)
    gaussian_module = gaussian_module_match.group(1).strip()

    return functional.lower(), basis, dispersion.lower(), solvent.lower(), gaussian_module, sr_time, rootdir

# Covalent radii in Å (simplified example)
vdw_radii = {
    'H': 1.20, 'He': 1.40, 'Li': 1.82, 'Be': 1.53, 'B': 1.92,
    'C': 1.70, 'N': 1.55, 'O': 1.52, 'F': 1.47, 'Ne': 1.54,
    'Na': 2.27, 'Mg': 1.73, 'Al': 1.84, 'Si': 2.10, 'P': 1.80,
    'S': 1.80, 'Cl': 1.75, 'Ar': 1.88, 'K': 2.75, 'Ca': 2.31,
    'Sc': 2.30, 'Ti': 2.15, 'V': 2.05, 'Cr': 2.05, 'Mn': 2.05,
    'Fe': 2.00, 'Co': 2.00, 'Ni': 1.97, 'Cu': 1.96, 'Zn': 2.01,
    'Ga': 1.87, 'Ge': 2.11, 'As': 1.85, 'Se': 1.90, 'Br': 1.85,
    'Kr': 2.02, 'Rb': 3.03, 'Sr': 2.49, 'Y': 2.40, 'Zr': 2.30,
    'Nb': 2.15, 'Mo': 2.10, 'Tc': 2.05, 'Ru': 2.05, 'Rh': 2.00,
    'Pd': 2.05, 'Ag': 2.03, 'Cd': 2.18, 'In': 1.93, 'Sn': 2.17,
    'Sb': 2.06, 'Te': 2.06, 'I': 1.98, 'Xe': 2.16, 'Cs': 3.43,
    'Ba': 2.68, 'La': 2.50, 'Ce': 2.48, 'Pr': 2.47, 'Nd': 2.45,
    'Pm': 2.43, 'Sm': 2.42, 'Eu': 2.40, 'Gd': 2.38, 'Tb': 2.37,
    'Dy': 2.35, 'Ho': 2.33, 'Er': 2.32, 'Tm': 2.30, 'Yb': 2.28,
    'Lu': 2.27, 'Hf': 2.25, 'Ta': 2.20, 'W': 2.10, 'Re': 2.05,
    'Os': 2.00, 'Ir': 2.00, 'Pt': 2.05, 'Au': 2.10, 'Hg': 2.05,
    'Tl': 1.96, 'Pb': 2.02, 'Bi': 2.07, 'Po': 1.97, 'At': 2.02,
    'Rn': 2.20, 'Fr': 3.48, 'Ra': 2.83, 'Ac': 2.60, 'Th': 2.37,
    'Pa': 2.43, 'U': 2.40, 'Np': 2.39, 'Pu': 2.43, 'Am': 2.44,
    'Cm': 2.45, 'Bk': 2.44, 'Cf': 2.45, 'Es': 2.45, 'Fm': 2.45,
    'Md': 2.45, 'No': 2.45, 'Lr': 2.45
}




def extract_coordinates_from_log(file_path):
    """Extracts atomic symbols and coordinates from the LAST Standard orientation block."""
    with open(file_path, 'r') as file:
        lines = file.readlines()

    block_starts = [i + 5 for i, line in enumerate(lines) if "Standard orientation" in line]
    if not block_starts:
        print(f"Warning: 'Standard orientation' block not found in {file_path}")
        return []

    start = block_starts[-1]
    for i in range(start, len(lines)):
        if "-----" in lines[i]:
            end = i
            break

    extracted_data = []
    for line in lines[start:end]:
        parts = line.split()
        if len(parts) >= 6:
            atom_number = int(parts[1])
            x, y, z = map(float, parts[3:6])
            atom_symbol = atomic_symbols.get(atom_number, f"Unknown({atom_number})")
            extracted_data.append([atom_symbol, x, y, z])
    return extracted_data



def write_gaussian_input(file_name, molecule, suffix, parameter_file="parameters.txt"):
    """Writes a Gaussian input file based on parameters.txt."""
    functional, basis, dispersion, solvent = read_dft_parameters(parameter_file)
    output_path = f"{file_name}_{suffix}.gjf"

    with open(output_path, 'w') as output_file:
        # Write header
        output_file.write(f"%nprocshared=8\n")
        output_file.write(f"%mem=16GB\n")
        output_file.write(f"%chk={file_name}_{suffix}.chk\n")
        output_file.write(f"# opt=calcfc freq {functional} {basis} {dispersion}\n\n")
        output_file.write(f"{file_name}_{suffix} optfreq\n\n")
        output_file.write("0 1\n")
        
        for atom in molecule:
            output_file.write(f" {atom[0]:<2} {atom[1]:>15.8f} {atom[2]:>15.8f} {atom[3]:>15.8f}\n")
        output_file.write("\n")

        # Only add Link1 if basis was originally "cbs" → i.e. basis was translated to "cc-pvdz"
        if basis == "cc-pvdz":
            link1_text = f"""--Link1--
%nprocshared=8
%mem=16GB
%chk={file_name}_{suffix}.chk
# {functional} cc-pvtz {dispersion} Geom=Checkpoint

{file_name}_{suffix} E_ccpvtz

0 1

--Link1--
%nprocshared=8
%mem=16GB
%chk={file_name}_{suffix}.chk
# {functional} cc-pvqz {dispersion} Geom=Checkpoint

{file_name}_{suffix} E_ccpvqz

0 1

"""
            output_file.write(link1_text)

    return output_path

def create_submission_script(job_name, input_file, output_file, sr_time, gaussian_module):
    """Creates a SLURM submission script."""
    script_name = f"{job_name}.sub"
    with open(script_name, 'w') as script:
        script.write(f"""#!/bin/sh
#SBATCH --job-name={job_name}
#SBATCH --cpus-per-task=12
#SBATCH --output={job_name}.logfile
#SBATCH --time={sr_time}
#SBATCH --partition=zen4
#SBATCH --mem-per-cpu=5GB

module load {gaussian_module}
export GAUSS_SCRDIR=$VSC_SCRATCH_VO_USER/gauss_scrdir$SLURM_JOB_ID
mkdir -p $GAUSS_SCRDIR
g16 -p=$SLURM_CPUS_PER_TASK -m=80GB < {input_file} > {output_file}
""")

    return script_name

def launcher(log_files, parameters_file, dependency_script):
    MAX_JOBS = 100
    """Generates Gaussian input files, submission scripts, and launches jobs."""
    sr_time, rootdir, gaussian_module = read_parameters(parameters_file)
    job_ids = []

    n=0

    for log_file in log_files:
        if n >= MAX_JOBS:
            break

        base_name = os.path.splitext(log_file)[0]
        atoms = extract_coordinates_from_log(log_file)

        if not atoms:
            print(f"Skipping {log_file}: No coordinates extracted.")
            continue
        
        # Build the graph of atoms connected by bonds

        positions = np.array([[x[1], x[2], x[3]] for x in atoms])
        
        # Step 1: Build connectivity graph manually
        def are_bonded(i, j):
            ri = radii.get(atoms[i][0], 1.5)
            rj = radii.get(atoms[j][0], 1.5)
            max_dist = (ri + rj) / 2
            dist = np.linalg.norm(positions[i] - positions[j])
            return dist < max_dist
        
        # Build adjacency list
        adj = {i: [] for i in range(len(atoms))}
        for i in range(len(atoms)):
            for j in range(i + 1, len(atoms)):
                if are_bonded(i, j):
                    adj[i].append(j)
                    adj[j].append(i)
        
        # Step 2: Find connected components (molecules) using BFS
        def find_connected_components(adj):
            visited = set()
            components = []
        
            for start in adj:
                if start not in visited:
                    queue = [start]
                    component = []
                    while queue:
                        node = queue.pop(0)
                        if node not in visited:
                            visited.add(node)
                            component.append(node + 1)  # convert to 1-based index
                            queue.extend(adj[node])
                    components.append(sorted(component))
            return components
        
        # Get molecules
        molecules = find_connected_components(adj)
        
        # Step 3: Assign molecules to local variables
        if len(molecules) < 2:
            print(f"Skipping {log_file}: only {len(molecules)} molecule(s) found.")
            with open("OnlyOneMoleculeFound.txt", "a") as log:
                log.write(f"{log_file} - only {len(molecules)} molecule(s) found\n")
            continue
        
        molecule1_indices = molecules[0]
        molecule2_indices = molecules[1]
        molecule1 = [atoms[i - 1] for i in molecule1_indices]
        molecule2 = [atoms[i - 1] for i in molecule2_indices]
        

        # Write input files
        input_R1 = write_gaussian_input(base_name, molecule1, "R1")
        input_R2 = write_gaussian_input(base_name, molecule2, "R2")
        
        # Create and submit jobs
        for suffix, input_file in zip(["R1", "R2"], [input_R1, input_R2]):
            output_file = input_file.replace(".gjf", ".log")
            job_name = f"{base_name}_{suffix}"
            script_name = create_submission_script(job_name, input_file, output_file, sr_time, gaussian_module)

            result = subprocess.run(f"sbatch {script_name}", shell=True, stdout=subprocess.PIPE, text=True)
            if result.returncode == 0:
                job_id = re.search(r'(\d+)', result.stdout)
                if job_id:
                    job_ids.append(job_id.group(1))
                    print(f"Submitted job {job_name} with ID {job_id.group(1)}")
                    n+=1
            else:
                print(f"Failed to submit job {job_name}: {result.stderr}")

    # Launch dependent script
    if job_ids:
        dependency_str = ":".join(job_ids)

        command = f"sbatch --dependency=afterany:{dependency_str} {dependency_script}"
        
        result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, text=True)
        if result.returncode == 0:
            print(f"Dependent script {dependency_script} submitted successfully.")
        else:
            print(f"Failed to submit dependent script: {result.stderr}")
    else:
        print("No jobs submitted, skipping dependent script.")


# Main workflow
if __name__ == "__main__":
    log_files = [f for f in os.listdir("./") if f.endswith("Complex.log")]
    parameters_file = "parameters.txt"
    functional, basis, dispersion, solvent, gaussian_module, sr_time, rootdir = read_parameters(parameters_file)
    dependency_script = os.path.join(rootdir, '3_Results.sub')
    launcher(log_files, parameters_file, dependency_script)