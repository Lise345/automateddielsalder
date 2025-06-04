# 🧪 Automated Quantum Chemistry Workflow

This repository automates Gaussian-based DFT analysis of reaction mechanisms, covering TS validation, IRC calculations, extrapolated energy analysis, and optional IR/NMR and reagent calculations. It uses SLURM for job scheduling and supports multi-level theory extrapolation.

---

## 📁 Project Structure

| File / Script                        | Description |
|-------------------------------------|-------------|
| `Automated_IRC_Calculation.py`      | Launches IRC calculations from validated TSs. Converts geometries to `.xyz`, generates `.gjf`, submits jobs, and schedules `2_StationaryPoints_calculator.sub`. |
| `Automated_IRC_Extractor.py`        | Extracts endpoints after IRC. Determines Complex vs Product, submits SP calculations, and triggers `3_Results-2.sub`. |
| `Automated_Results.py`              | Collects all Gaussian output energies, computes extrapolated ΔG and ΔH values, and generates Excel and plots. |
| `Automated_IRandNMR.py`             | Prepares IR/NMR input files for Products and launches Gaussian jobs via `sgx16`. |
| `Automated_ProcessIRandNMR.py`      | Parses `.log` files for IR/NMR data and creates plots. |
| `Automated_Reagents.py`             | *(Optional)* Automatically extracts individual reagent geometries from Complex, submits Gaussian jobs, and schedules `3_Results.sub` once completed. |
| `1_IRC_calculator.sub`              | **Entry point.** Creates a virtual environment, installs dependencies, and launches `Automated_IRC_Calculation.py`. |
| `2_StationaryPoints_calculator.sub` | Triggered **automatically** after IRC jobs. Generates and submits SP jobs. |
| `3_Results-2.sub`                   | Triggered **automatically** after SP jobs. Runs energy analysis. |
| `4_SeparateReagents.sub`            | *(Optional)* Manual launch of `Automated_Reagents.py`. Calculates energies for separate reagents from Complex. Must be followed by manual re-run of `3_Results-2.sub`. |
| `5_IRandNMR.sub`                    | *(Optional)* Submits IR/NMR frequency jobs for each Product. Can be followed by `Automated_ProcessIRandNMR.py` to generate plots. |

---

## ⚙️ How to Use

### 1. 🧾 Edit `parameters.txt`

Example configuration:
```ini
rootdir '/your/project/path'
bin /path/to/bin
Gaussian module gaussian/g16.b01
Functional M06-2X
Dispersion empiricaldispersion=gd3
Basis cc-pvdz
DFT solvent gas
size_molecule = 25
CC1_in = "1 6"
Time for Separate Reagent calcs 25
```

---

### 2. 🚀 Main Workflow

```bash
git clone https://github.com/Lise345/automateddielsalder.git
sbatch 1_IRC_calculator.sub
```

This:
- Creates a virtual environment
- Installs required packages
- Runs `Automated_IRC_Calculation.py`, which:
  - Validates TSs
  - Generates IRC `.gjf` files
  - Submits IRC jobs
  - Schedules `2_StationaryPoints_calculator.sub`

Then:
- `Automated_IRC_Extractor.py` submits Complex/Product calculations
- `3_Results-2.sub` analyzes energy and produces Excel/plots

---

### 3. 🧪 Optional: Calculate Reagent Energies

To calculate energies of separate reagents from Complex:

```bash
sbatch 4_SeparateReagents.sub
```

This launches `Automated_Reagents.py` which:
- Automatically splits reagents using distance-based bonding
- Submits jobs for R1 and R2
- Schedules `3_Results.sub` on completion

> 🧠 If you manually use this script instead of 4_, run:
```bash
python Automated_Reagents.py
```

---

### 4. 🌈 Optional: IR/NMR Spectra

```bash
sbatch 5_IRandNMR.sub
```

Then:
```bash
python Automated_ProcessIRandNMR.py
```

---

## 📦 Dependencies

✅ **No manual installation required**  
`1_IRC_calculator.sub` sets up a virtual environment and installs:
- `numpy`, `pandas`, `matplotlib`, `openpyxl`, etc.

---

## 📊 Output

- Excel: `AutomatedDA_results.xlsx`
- Plots: `plots/` and `spectra_png/`
- SLURM logs: `*.logfile`

---

## 📌 Notes

- TS validation: 1 imaginary frequency between −1000 and −200 cm⁻¹
- Extrapolation: PVDZ–PVTZ–PVQZ 3-point formula
- Job submission handled via SLURM and `sgx16` (if configured)
- Separate reagent analysis is **optional**, but recommended for more accurate ΔG and ΔH values
