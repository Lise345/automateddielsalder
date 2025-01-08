### Sequence of Execution

1. **1_IRC_calculator.sub**
    - This script sets up the environment and launches the `Automated_IRC_Calculation.py` script.
    - It includes loading necessary modules and setting paths.
2. **Automated_IRC_Calculation.py**
    - This script performs several tasks:
        - Reads parameters from `parameters.txt`.
        - Compiles and checks frequencies in transition state log files.
        - Converts Gaussian log files to XYZ format.
        - Generates IRC input files and submits Gaussian jobs.
        - Saves results to an Excel file.
    - The script also sets up job dependencies to launch `2_StationaryPoints_calculator.sub` when all IRC calculations are finished.
3. **2_StationaryPoints_calculator.sub**
    - This script sets up the environment and launches the `Automated_IRC_Extractor.py` script.
    - It includes loading necessary modules and setting paths.
4. **Automated_IRC_Extractor.py**
    - This script processes IRC log files, extracts geometries, and generates input files for stationary point (SP) calculations.
    - Key Functions:
        - `geometryextractor()`: Extracts geometries from log files.
        - `inputgenerator()`: Generates input files for SP calculations.
        - `launcherstatp()`: Submits SP calculation jobs.
        - `launcherTS()`: Submits TS calculation jobs.
        - `launch_dependent_job()`: Submits `3_Results.sub` as a dependent job.
5. **3_Results.sub**
    - This script sets up the environment and launches the `Automated_Results_CREST.py` script.
    - It includes loading necessary modules and setting paths.
6. **Automated_Results_CREST.py**
    - This script processes the results, calculates energies, and generates a summary in an Excel file.
