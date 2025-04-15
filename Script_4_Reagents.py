import os
import re
import subprocess

# Atomic symbols mapping
atomic_symbols = {
    1: "H", 6: "C", 7: "N", 8: "O", 9: "F", 16: "S", 17: "Cl"  # Add other atomic symbols if needed
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

def read_parameters(file_path):
    """Reads molecule atom indices from parameters.txt."""
    with open(file_path, 'r') as file:
        file_content = file.read()

    molecule1_atoms = list(map(int, re.search(r'molecule1_atoms\s*=\s*(.+)', file_content).group(1).split()))
    molecule2_atoms = list(map(int, re.search(r'molecule2_atoms\s*=\s*(.+)', file_content).group(1).split()))

    rootdir_match = re.search(r'rootdir\s+(.+)', file_content)
    rootdir = rootdir_match.group(1).strip().strip("'\"")
    
    time_sr = re.search(r'Time for Separate Reagent calcs\s+(\d+)', file_content)
    if time_sr:
        time_for_sr_calcs = int(time_sr.group(1))
        sr_time = f'{time_for_sr_calcs}:00:00'  # Format to HH:MM:SS
    else:
        sr_time = '25:00:00'  # Default value if not found
        print("Time for stationary calculations not found, defaulting to 25:00:00")

    return molecule1_atoms, molecule2_atoms, sr_time, rootdir

def write_gaussian_input(file_name, molecule, suffix):
    """Writes a Gaussian input file."""
    output_path = f"{file_name}_{suffix}.gjf"
    header = f"""%nprocshared=8
%mem=16GB
%chk={file_name}_{suffix}.chk
# opt=calcfc freq m062x cc-pvdz empiricaldispersion=gd3 

{file_name}_{suffix} optfreq

0 1\n"""
    with open(output_path, 'w') as output_file:
        output_file.write(header)
        for atom in molecule:
            output_file.write(f" {atom[0]:<2} {atom[1]:>15.8f} {atom[2]:>15.8f} {atom[3]:>15.8f}\n")
        output_file.write("\n")

# Append the required "Link1" sections with the correct name
        link1_text = f"""--Link1--
%nprocshared=8
%mem=16GB
%chk={file_name}_{suffix}.chk
# m062x cc-pvtz empiricaldispersion=gd3 Geom=Checkpoint

{file_name}_{suffix} E_ccpvtz

0 1

--Link1--
%nprocshared=8
%mem=16GB
%chk={file_name}_{suffix}.chk
# m062x cc-pvqz empiricaldispersion=gd3 Geom=Checkpoint

{file_name}_{suffix} E_ccpvqz

0 1

"""
        output_file.write(link1_text)

    return output_path

def create_submission_script(job_name, input_file, output_file,sr_time):
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

module load Gaussian/G16.A.03-intel-2022a
export GAUSS_SCRDIR=$VSC_SCRATCH_VO_USER/gauss_scrdir$SLURM_JOB_ID
mkdir -p $GAUSS_SCRDIR
g16 -p=$SLURM_CPUS_PER_TASK -m=80GB < {input_file} > {output_file}
""")
        script.write('rm -r ${VSC_SCRATCH_VO_USER:?}/{gauss_scrdir$SLURM_JOB_ID:?}\n')
        script.write('\n')
    return script_name

def launcher(log_files, parameters_file, dependency_script):
    MAX_JOBS = 100
    """Generates Gaussian input files, submission scripts, and launches jobs."""
    molecule1_indices, molecule2_indices, sr_time, rootdir = read_parameters(parameters_file)
    job_ids = []

    n=0

    while n < MAX_JOBS:
        for log_file in log_files:
            base_name = os.path.splitext(log_file)[0]
            extracted_atoms = extract_coordinates_from_log(log_file)
            
            if not extracted_atoms:
                print(f"Skipping {log_file}: No coordinates extracted.")
                continue

            # Extract molecules
            molecule1 = [extracted_atoms[i - 1] for i in molecule1_indices]
            molecule2 = [extracted_atoms[i - 1] for i in molecule2_indices]

            # Write input files
            input_R1 = write_gaussian_input(base_name, molecule1, "R1")
            input_R2 = write_gaussian_input(base_name, molecule2, "R2")
            
            # Create and submit jobs
            for suffix, input_file in zip(["R1", "R2"], [input_R1, input_R2]):
                output_file = input_file.replace(".gjf", ".log")
                job_name = f"{base_name}_{suffix}"
                script_name = create_submission_script(job_name, input_file, output_file, sr_time)

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
    molecule1_atoms, molecule2_atoms, sr_time, rootdir = read_parameters(parameters_file)
    dependency_script = os.path.join(rootdir, '5_BOLTCAR_Results.sub')
    launcher(log_files, parameters_file, dependency_script)

