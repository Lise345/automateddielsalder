import os
import re
import subprocess

# === 1. Load user parameters ===
with open('./parameters.txt', 'r') as parameters:
    file_content = parameters.read()

    rootdir = re.search(r'rootdir (.+)', file_content).group(1).strip().strip("'\"")
    directory = rootdir

    functional = re.search(r'Functional\s+(\S+)', file_content, re.IGNORECASE).group(1)
    dispersion = re.search(r'Dispersion\s+(\S+)', file_content, re.IGNORECASE).group(1)
    basis_raw = re.search(r'Basis\s+(\S+)', file_content, re.IGNORECASE).group(1).lower()
    solvent = re.search(r'DFT solvent\s+(\S+)', file_content, re.IGNORECASE).group(1).lower()

    basis = "cc-pvdz" if basis_raw == "cbs" else basis_raw

# === 2. Try to get the sgx16 path ===
# Default path from virtualenv bin
venv_path = os.environ.get("VIRTUAL_ENV")
if not venv_path:
    raise EnvironmentError("❌ VIRTUAL_ENV is not set. Please activate your virtual environment.")

sgx16_path = os.path.join("../Scripts", "sgx16")
if not os.path.isfile(sgx16_path):
    raise FileNotFoundError(f"❌ sgx16 not found at: {sgx16_path}")

# === 3. Atomic symbols ===
atomic_symbols = {i: e for i, e in enumerate(
    ["H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
     "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
     "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
     "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr",
     "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
     "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
     "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
     "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
     "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
     "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm",
     "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds",
     "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"], 1)}

# === 4. Geometry extraction ===
def extract_last_geometry(logfile):
    with open(logfile, 'r') as file:
        lines = file.readlines()

    input_orient_indices = [i for i, line in enumerate(lines) if "Input orientation:" in line]
    if not input_orient_indices:
        print(f"❌ No geometry found in {logfile}")
        return []

    geom_start = input_orient_indices[-1] + 5
    geometry = []

    for line in lines[geom_start:]:
        if "---------------------------------------------------------------------" in line:
            break
        geometry.append(line)

    return geometry

def convert_geometry(raw_geometry):
    formatted_geometry = []
    for line in raw_geometry:
        parts = line.split()
        if len(parts) < 6:
            continue
        atom_num = int(parts[1])
        symbol = atomic_symbols.get(atom_num, "X")
        if symbol == "X":
            print(f"⚠ Unknown atom number: {atom_num}")
        formatted_geometry.append([symbol, parts[3], parts[4], parts[5]])
    return formatted_geometry

# === 5. Write .gjf file ===
def write_gjf(geometry, filename):
    name = os.path.splitext(os.path.basename(filename))[0]

    # --- Build route line dynamically ---
    route_parts = ["nmr=giao", "freq", f"{basis}"]

    # Add solvent only if not 'none'
    if solvent.lower() != "none":
        route_parts.append(solvent)

    # Add dispersion only if not 'none'
    if dispersion != "none":
        route_parts.append(dispersion)

    route_line = "# " + " ".join(route_parts) + "\n"

    with open(filename, 'w') as f:
        f.write("%nprocshared=8\n")
        f.write("%mem=32GB\n")
        f.write(f"%chk={filename[:-4]}.chk\n")
        f.write(route_line +  "b3lyp iop(3/76=1000001189,3/77=0961409999,3/78=0000109999) \n\n")
        f.write(f"{filename[:-4]} IR and NMR calculation\n\n")
        f.write("0 1\n")
        for atom in geometry:
            f.write(f"{atom[0]} {atom[1]} {atom[2]} {atom[3]}\n")
        f.write("\n")
    print(f"✅ Input file written: {filename}")


# === 6. Run Gaussian using sgx16 ===
def run_gaussian(gjf_file):
    print(f"Launching Gaussian via sgx16: {gjf_file}", flush=True)
    try:
        result = subprocess.run([sgx16_path, gjf_file, "15"], check=True, capture_output=True, text=True)
        print(f"✅ Gaussian completed for: {gjf_file}", flush=True)
        print(result.stdout, flush=True)
    except subprocess.CalledProcessError as e:
        print(f"❌ Gaussian failed for {gjf_file}", flush=True)
        print("STDOUT:", e.stdout)
        print("STDERR:", e.stderr)
    except FileNotFoundError as e:
        print(f"❌ sgx16 not found: {e}", flush=True)


# === 7. Main loop ===
for file in os.listdir(directory):
    if file.lower().endswith(("product.log", "r1.log", "r2.log")):
        path = os.path.join(directory, file)
        print(f"🔍 Processing {file}", flush=True)

        raw_geom = extract_last_geometry(path)
        if not raw_geom:
            print(f"⚠ Skipping {file} (no geometry found)", flush=True)
            continue

        geom = convert_geometry(raw_geom)
        gjf_path = file.replace(".log", "_IR_NMR.gjf")
        write_gjf(geom, gjf_path)
        run_gaussian(gjf_path)
