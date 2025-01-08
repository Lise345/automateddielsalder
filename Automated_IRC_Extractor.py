import os
import sys
import subprocess
import shutil
import re
import numpy as np
import glob
import pandas as pd
import openpyxl
import math

with open('./parameters.txt', 'r') as parameters:
    file_content = parameters.read()

    # Extract the size_molecule value
    size_molecule_match = re.search(r'size_molecule\s*=\s*(\d+)', file_content)
    if size_molecule_match:
        size_molecule = int(size_molecule_match.group(1))

    # Extract the CC1_in value
    CC1_in_match = re.search(r'CC1_in\s*=\s*"(.*?)"', file_content)
    if CC1_in_match:
        CC1_in = list(map(int, CC1_in_match.group(1).split()))

    rootdir = re.search(r'rootdir (.+)', file_content)
    rootdir = rootdir.group(1).strip().strip("'\"")

    binfolder = re.search(r'bin (.+)', file_content)
    binfolder = binfolder.group(1)

CC1_out=CC1_in
CC1=CC1_in

#-------------Prepare TS geometries--------------

workbook = openpyxl.load_workbook("TS_analysis.xlsx")
sheet = workbook.active


data = []
for column in sheet.iter_cols(values_only=True):
    data.append(list(column))

print(data)

listofTS=[]
for caseofTS in data:
    if data.index(caseofTS)==2:
        for i in caseofTS:
            if caseofTS.index(i)>0 and i!=None:
                listofTS.append(i)

def SP_inputgenerator(xyzfile,filename):
    with open(xyzfile, 'r') as file:
        lines = file.readlines()
    with open(filename, 'x') as ip:
        ip.writelines("%nprocshared=4\n")
        ip.writelines("%mem=32GB\n")
        ip.writelines("%chk="+filename[:-4]+".chk"+"\n")
        ip.writelines("# opt=(calcfc,ts,noeigentest) freq cc-pvdz empiricaldispersion=gd3 m062x\n")
        ip.writelines("\n")
        Title=filename[:-4]+" "+"cc-pVTZ_SP"+"\n"
        ip.writelines(Title)
        ip.writelines("\n")
        ip.writelines("0 1\n")
        
        for atom in lines[2:]:
            ip.writelines(atom)
        ip.writelines("\n")
        
        #Writing Link1 part for cc-pVTZ
        ip.writelines("--Link1--\n")
        ip.writelines("%nprocshared=4\n")
        ip.writelines("%mem=32GB\n")
        ip.writelines("%chk="+filename[:-4]+".chk"+"\n")
        ip.writelines("# m062x cc-pvtz empiricaldispersion=gd3 Geom=Checkpoint \n")
        ip.writelines("\n")
        Title=filename[:-4]+" "+"cc-pVQT_SP"+"\n"
        ip.writelines(Title)
        ip.writelines("\n")
        ip.writelines("0 1\n")
        ip.writelines("\n")
        

	    #Writing Link1 part for cc-pVQZ
        ip.writelines("--Link1--\n")
        ip.writelines("%nprocshared=4\n")
        ip.writelines("%mem=32GB\n")
        ip.writelines("%chk="+filename[:-4]+".chk"+"\n")
        ip.writelines("# m062x cc-pvqz empiricaldispersion=gd3 Geom=Checkpoint \n")
        ip.writelines("\n")
        Title=filename[:-4]+" "+"cc-pVQZ_SP"+"\n"
        ip.writelines(Title)
        ip.writelines("\n")
        ip.writelines("0 1\n")
        ip.writelines("\n")
        ip.close()
    return print("Input generated for " + filename[:-4])

#-------------Prepare Complex and Product geometries---------------



atomic_symbols = {
    1: "H", 2: "He", 3: "Li", 4: "Be", 5: "B", 6: "C", 7: "N", 8: "O", 9: "F", 10: "Ne",
    11: "Na", 12: "Mg", 13: "Al", 14: "Si", 15: "P", 16: "S", 17: "Cl", 18: "Ar", 19: "K", 20: "Ca",
    21: "Sc", 22: "Ti", 23: "V", 24: "Cr", 25: "Mn", 26: "Fe", 27: "Co", 28: "Ni", 29: "Cu", 30: "Zn",
    31: "Ga", 32: "Ge", 33: "As", 34: "Se", 35: "Br", 36: "Kr", 37: "Rb", 38: "Sr", 39: "Y", 40: "Zr",
    41: "Nb", 42: "Mo", 43: "Tc", 44: "Ru", 45: "Rh", 46: "Pd", 47: "Ag", 48: "Cd", 49: "In", 50: "Sn",
    51: "Sb", 52: "Te", 53: "I", 54: "Xe", 55: "Cs", 56: "Ba", 57: "La", 58: "Ce", 59: "Pr", 60: "Nd",
    61: "Pm", 62: "Sm", 63: "Eu", 64: "Gd", 65: "Tb", 66: "Dy", 67: "Ho", 68: "Er", 69: "Tm", 70: "Yb",
    71: "Lu", 72: "Hf", 73: "Ta", 74: "W", 75: "Re", 76: "Os", 77: "Ir", 78: "Pt", 79: "Au", 80: "Hg",
    81: "Tl", 82: "Pb", 83: "Bi", 84: "Po", 85: "At", 86: "Rn", 87: "Fr", 88: "Ra", 89: "Ac", 90: "Th",
    91: "Pa", 92: "U", 93: "Np", 94: "Pu", 95: "Am", 96: "Cm", 97: "Bk", 98: "Cf", 99: "Es", 100: "Fm",
    101: "Md", 102: "No", 103: "Lr", 104: "Rf", 105: "Db", 106: "Sg", 107: "Bh", 108: "Hs", 109: "Mt",
    110: "Ds", 111: "Rg", 112: "Cn", 113: "Nh", 114: "Fl", 115: "Mc", 116: "Lv", 117: "Ts", 118: "Og"
}

def geometryextractor(logfile):
    with open(logfile, 'r') as file:
        lines = file.readlines()

    convergence_indices = []
    for i, line in enumerate(lines):
        if "Input orientation" in line:
            convergence_indices.append(i)

    conv_geom=convergence_indices[-1]+5
    
    size_molecule=0
    
    geometry= []
    start=conv_geom
    index=conv_geom

    for index in range(conv_geom, len(lines)):
        if "--------------" in lines[index]:
            break
        geometry.append(lines[index])
        size_molecule =+ 1
    return geometry

def geometryconverter(geometry):
    updated_geometry=[]
    for atom in geometry:
        updated_atom=atom.split()
        atom_number = int(updated_atom[1])
        if atom_number in atomic_symbols:
            updated_atom[1] = atomic_symbols[atom_number]
        else:
            print("Atom does not exist")
        del updated_atom[2]
        updated_geometry.append(updated_atom)
    #print(updated_geometry)
    return updated_geometry

def inputgenerator(geometry, filename):
    with open(filename, 'x') as ip:
        ip.writelines("%nprocshared=8\n")
        ip.writelines("%mem=32GB\n")
        ip.writelines("%chk="+filename[:-4]+".chk"+"\n")
        ip.writelines("# opt=calcfc freq m062x cc-pvdz empiricaldispersion=gd3\n")
        ip.writelines("\n")
        Title=filename[:-4]+" "+"optfreq"+"\n"
        ip.writelines(Title)
        ip.writelines("\n")
        ip.writelines("0 1\n")
        
        for atom in geometry:
            atom_line=atom[1]+" "+atom[2]+" "+atom[3]+" "+atom[4]+"\n"
            ip.writelines(atom_line)
        ip.writelines("\n")
        
        #Writing Link1 part for cc-pVTZ
        ip.writelines("--Link1--\n")
        ip.writelines("%nprocshared=4\n")
        ip.writelines("%mem=32GB\n")
        ip.writelines("%chk="+filename[:-4]+".chk"+"\n")
        ip.writelines("# m062x cc-pvtz empiricaldispersion=gd3 Geom=Checkpoint\n")
        ip.writelines("\n")
        Title=filename[:-4]+" "+"E_ccpvtz"+"\n"
        ip.writelines(Title)
        ip.writelines("\n")
        ip.writelines("0 1\n")
        ip.writelines("\n")
        
        #Writing Link1 part for cc-pVQZ
        ip.writelines("--Link1--\n")
        ip.writelines("%nprocshared=8\n")
        ip.writelines("%mem=32GB\n")
        ip.writelines("%chk="+filename[:-4]+".chk"+"\n")
        ip.writelines("# m062x cc-pvqz empiricaldispersion=gd3 Geom=Checkpoint\n")
        Title=filename[:-4]+" "+"E_ccpvqz"+"\n"
        ip.writelines("\n")
        ip.writelines(Title)
        ip.writelines("\n")
        ip.writelines("0 1\n")
        ip.writelines("\n")
        ip.close()
    return print("Input generated for " + filename[:-4])

def atomwithfloats(atom):
    updated_atom=[]
    for el in atom:
        if atom.index(el)==0 or atom.index(el)==1:
            updated_atom.append(el)
        else:
            up_el=float(el)
            updated_atom.append(up_el)
    return updated_atom

def distance(geometry, coordinates):
    atom1=atomwithfloats(geometry[coordinates[0]])
    atom2=atomwithfloats(geometry[coordinates[1]])
    r_sq=abs(((atom1[2]-atom2[2])**2)+((atom1[3]-atom2[3])**2)+((atom1[4]-atom2[4])**2))
    
    r=math.sqrt(r_sq)
    return r


#-------------Launch calculations-----------------


import os
import re
import subprocess

inp_file_job_ids = []

def launcherstatp(logfilelist):
    print(logfilelist)
    for i, logfile1 in enumerate(logfilelist):
    # Extract the base name of the first logfile
        base_name1 = re.sub(r'(reverse|forward)\.log$', '', logfile1)
    
        for j, logfile2 in enumerate(logfilelist):
            if i < j:  # Ensure each pair is processed only once
            # Extract the base name of the second logfile
                base_name2 = re.sub(r'(reverse|forward)\.log$', '', logfile2)
            
            # Check if both files share the same base name and differ by reverse/forward
                if base_name1 == base_name2 and (
                    logfile1.endswith("reverse.log") and logfile2.endswith("forward.log") or
                    logfile1.endswith("forward.log") and logfile2.endswith("reverse.log")
                ):
                # Process the pair
                    print(f"Matched pair: {logfile1}, {logfile2}")
                
                    geometry1 = geometryextractor(logfile1)
                    updated_geometry1 = geometryconverter(geometry1)
                
                    geometry2 = geometryextractor(logfile2)
                    updated_geometry2 = geometryconverter(geometry2)
                
                    # Determine complex and product filenames based on distance comparison
                    distance1_CC1 = distance(updated_geometry1, CC1)
                    distance2_CC1 = distance(updated_geometry2, CC1)
                
                    if distance1_CC1 > distance2_CC1:
                        filename1 = logfile1[:-4] + "_Complex.gjf"
                        filename2 = logfile2[:-4] + "_Product.gjf"
                    else:
                        filename1 = logfile2[:-4] + "_Complex.gjf"
                        filename2 = logfile1[:-4] + "_Product.gjf"
                
                    inputgenerator(updated_geometry1, filename1)
                    inputgenerator(updated_geometry2, filename2)
                
                # Create submission script
                    reduced_filename = logfile1[:-4] + "_optE"
                    with open(reduced_filename + ".sub", "w") as gsub:
                        gsub.write('#!/bin/sh\n')
                        gsub.write(f'#SBATCH --job-name={reduced_filename}\n')
                        gsub.write('#SBATCH --ntasks=12\n')
                        gsub.write(f'#SBATCH --output={reduced_filename}.logfile\n')
                        gsub.write('#SBATCH --time=40:00:00\n')
                        gsub.write('\n')
                        gsub.write('module load Gaussian/G16.A.03-intel-2022a\n')
                        gsub.write('export GAUSS_SCRDIR=$TMPDIR\n')
                        gsub.write('mkdir -p $GAUSS_SCRDIR\n')
                        gsub.write(f'g16 < {filename1} > {filename1[:-4]}.log\n')
                        gsub.write(f'g16 < {filename2} > {filename2[:-4]}.log\n')

                    
                    # Submit the job
                    result = subprocess.run(
                        f"sbatch {reduced_filename}.sub",
                        shell=True,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE,
                        universal_newlines=True
                    )
                    if result.returncode == 0:
                        job_id_match = re.search(r'(\d+)', result.stdout)
                        if job_id_match:
                            job_id = job_id_match.group(1)
                            inp_file_job_ids.append(job_id)
                            print(f"Job submitted: {job_id}")
                    else:
                         print(f"Failed to submit job: {result.stderr}")
           

def launcherTS(xyzlist):
    for xyzfile in xyzlist:
        filename=xyzfile[:-4]+"_SP.gjf"
        reduced_filename=filename[:-4]
        outputfilename=xyzfile[:-4]+"_SP.log"
        with open(reduced_filename+".sub","w") as gsub:
            gsub.write('#!/bin/sh\n')
            gsub.write(f'#SBATCH --job-name='+filename[:-4]+'\n')
            gsub.write('#SBATCH --ntasks=12\n')
            gsub.write(f'#SBATCH --output='+filename[:-4]+'.logfile\n')
            gsub.write('#SBATCH --time=40:00:00\n')
            gsub.write('\n')
            gsub.write('# Loading modules\n')
            gsub.write('module load Gaussian/G16.A.03-intel-2022a\n')  # Adjust based on the available Gaussian module
            gsub.write('\n')
            gsub.write('# Setting up Gaussian environment\n')
            gsub.write('export GAUSS_SCRDIR=$TMPDIR\n')  # Temporary directory for Gaussian scratch files
            gsub.write('mkdir -p $GAUSS_SCRDIR\n')
            gsub.write('#Launching calculation\n')
            gsub.write('export PATH={binfolder}:$PATH\n')
            
            SP_inputgenerator(xyzfile,filename)
            
            gsub.write(f'g16 < {filename} > {outputfilename}\n')
            gsub.write('\n')
            gsub.close()
        
# Submit the job
        result = subprocess.run(
            f"sbatch {reduced_filename}.sub",
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True
            )
        if result.returncode == 0:
            job_id_match = re.search(r'(\d+)', result.stdout)
            if job_id_match:
                job_id = job_id_match.group(1)
                inp_file_job_ids.append(job_id)
                print(f"Job submitted: {job_id}")
        else:
                print(f"Failed to submit job: {result.stderr}") 



#-----dependent jobs -----
def launch_dependent_job():
    if inp_file_job_ids:
        dependency_str = ":".join(inp_file_job_ids)
        extractor_script = os.path.join(rootdir, '3_Results.sub')
        
        if not os.path.exists(extractor_script):
            raise FileNotFoundError(f"Extractor script not found: {extractor_script}")

        dependency_command = [
            "sbatch",
            f"--dependency=afterany:{dependency_str}",
            extractor_script
        ]
        result = subprocess.run(
            dependency_command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True
        )
        if result.returncode == 0:
            print(f"Extractor job submitted successfully: {result.stdout}")
        else:
            print(f"Failed to submit extractor job: {result.stderr}")
    else:
        print("No jobs were submitted, skipping dependency job submission.")


logfilelist=[]
errorfiles=[]
for file in os.listdir():
    if "IRC" in file and "log" in file and not "logfile" in file:
        with open(file, 'r') as f:
            lines = f.readlines()
            if "Normal termination" in lines[-1]:
                logfilelist.append(file)
            else:
                logfilelist.append(file)
                errorfiles.append(file)

print(logfilelist)

launcherstatp(logfilelist)
launcherTS(listofTS)
launch_dependent_job()
