# -*- coding: utf-8 -*-
"""
Created on Wed Nov 28 15:18:38 2024

@author: Lise
"""



import os
import sys
import subprocess
import shutil
import re
import numpy as np
import glob
import pandas as pd

with open('./parameters.txt', 'r') as parameters:
    file_content = parameters.read()

    rootdir = re.search(r'rootdir (.+)', file_content)
    rootdir = rootdir.group(1).strip().strip("'\"")

    binfolder = re.search(r'bin (.+)', file_content)
    binfolder = binfolder.group(1)




#Part that compiles and checks if the imaginary frequency lies within the #expected range and if we have the number of expected nimag 

def compile_frequencies(lines):
    frequencies=[]
    for line in lines:
        if "Frequencies" in line:
            numbers = re.findall(r'-?\d+\.\d+', line)
            numbers = [float(num) for num in numbers]
            frequencies= frequencies+numbers
    return frequencies

def checkfrequency(filename):
    with open(filename, "r") as readfile:
        lines = readfile.readlines()
        last_line = lines[-1].rstrip()
        last_line_2 = lines[-2].rstrip()
    if "Normal termination" in last_line:
        frequencies=compile_frequencies(lines)
        nimag=0
        for frequency in frequencies:
            if frequency<0:
                nimag=nimag+1
        if nimag != 1:
            print("incorrect number of imag freq for "+filename)
            if nimag>1:
                incorrectTS.append(filename)
        elif frequencies[0]>-1000.0000 and frequencies[0]<-200.0000:
            correctTS.append(filename)
        else:
            incorrectTS.append(filename)
    else:
        errorterm.append(filename)

listfiles=[]
for file in os.listdir():
    if "TS.log" in file:
        listfiles.append(file)
            
incorrectTS= []
correctTS= []
errorterm=[]

for file in listfiles:
    checkfrequency(file)

print("The correct TS are")
print(correctTS)


print("The incorrect TS are")
print(incorrectTS)

print("The error TS are")
print(errorterm)



#Convert files to their xyz

    
def lastgeometry(filename):
    with open(filename, "r") as readfile:
        lines = readfile.readlines()
        indices=[]
        for idx, line in enumerate(lines):  # Track index manually
            if "                         Standard orientation:                        " in line:
                indices.append(idx)
        last_index=indices[-1]
        last_index = indices[-1]
        
        start = last_index + 5
        i = start
        coord = []
        size_molecule = 0
        while i < len(lines):
            if "---------------------------------------------------------------------" in lines[i]:
                break
            strippedline = lines[i].split()
            number_list = [float(num) for num in strippedline]
            coord.append(number_list)
            i += 1
            size_molecule += 1
    return coord, size_molecule

# Dictionary mapping atomic numbers to element symbols 
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

def atomincoord(coord):
    for atom in coord:
        atom_number = atom[1]
        
        if atom_number in atomic_symbols:
            atom[1] = atomic_symbols[atom_number]
        else:
            print("Atom does not exist")
    return coord

IRClist=[]


def convert_gjf_to_xyz(filename):
    coord, size_molecule = lastgeometry(filename)
    newcoord=atomincoord(coord)
    
    reduced_filename=filename[:].strip(".log")
    newfile=reduced_filename + ".xyz"
    if newfile in os.listdir():
        os.remove(newfile)
    xyzfile=open(newfile, "x")
    
    xyzfile.writelines(str(size_molecule) + '\n')
    xyzfile.writelines('\n')
    for atom in newcoord:
        xyzfile.writelines(atom[1]+' '+str(atom[3])+' '+str(atom[4])+' '+str(atom[5])+'\n')
    IRClist.append(newfile)
    return 'Done'
        
for file in correctTS:
    convert_gjf_to_xyz(file)
        

def read_coordinates(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    
    start_index = 0
    for i, line in enumerate(lines):
        if '0 1' in line:
            start_index = i + 1
            break
    
    atoms = []
    for line in lines[start_index:]:
        parts = line.split()
        if len(parts) == 4:
            atoms.append([parts[0], float(parts[1]), float(parts[2]), float(parts[3])])
    return lines[:start_index], atoms

def write_coordinates(file_path, header, atoms):
    with open(file_path, 'w') as file:
        file.writelines(header)
        for atom in atoms:
            file.write(f"{atom[0]:<3} {atom[1]:>15.8f} {atom[2]:>15.8f} {atom[3]:>15.8f}\n")

def get_molecule_atoms(prompt):
    atoms = input(prompt).split()
    return [int(atom) - 1 for atom in atoms]

def get_atom(prompt):
    return int(input(prompt)) - 1

def XYZspliter():
    with open('crest_conformers.xyz', 'r') as rfile:
        lines = rfile.readlines()

    natoms = int(lines[0])
    ngeoms = len(lines) // (natoms + 2)

    for j in range(ngeoms):
        outname = f"xyzfilenum{j+1:04d}.xyz"
        with open(outname, "w") as ow:
            ow.write(str(natoms) + "\n \n")
            ow.writelines(lines[(j * (natoms + 2) + 2):((j + 1) * (natoms + 2))])



def energyfinder(logfile):
    with open(logfile, 'r') as file:
        lines = file.readlines()
    indexlist=[]
    for line in lines:
        if "SCF Done" in line:
            indexlist.append(lines.index(line))
    lastindex=indexlist[len(indexlist)-1]
    linelist=lines[lastindex].split()
    for word in linelist:
        if word=="=":
            Energy=str(float(linelist[linelist.index(word)+1])*627.503)
    return Energy




def save_to_excel(incorrectTS, correctTS, IRClist, filename="TS_analysis.xlsx"):
    # Create a dictionary with lists to save
    data = {
        "Correct TS": pd.Series(correctTS),
        "Incorrect TS": pd.Series(incorrectTS),
        "IRC List": pd.Series(IRClist)
    }

    # Convert the dictionary into a DataFrame
    df = pd.DataFrame(dict([(k, pd.Series(v)) for k,v in data.items()]))

    # Save DataFrame to an Excel file
    df.to_excel(filename, index=False, engine='openpyxl')
    print(f"Data successfully saved to {filename}")

save_to_excel(incorrectTS, correctTS, IRClist)

print("Number of IRC calculations: "+str(len(IRClist)))

#IRC calculation runner

def IRC_inputgenerator(xyzfile, filename, direction):
    with open(xyzfile, 'r') as file:
        lines = file.readlines()
    with open(filename, 'x') as ip:
        ip.writelines("%nprocshared=12\n")
        ip.writelines("%mem=12GB\n")
        ip.writelines("%chk="+filename[:-4]+".chk"+"\n")
        ip.writelines("# irc=("+direction+",calcfc,maxpoints=100,recalc=3,tight) m062x cc-pvdz empiricaldispersion=gd3\n")
        ip.writelines("\n")
        Title=filename+" "+"IRC"+direction+"\n"
        ip.writelines(Title)
        ip.writelines("\n")
        ip.writelines("0 1\n")
        
        for atom in lines[2:]:
            ip.writelines(atom)
        ip.writelines("\n")
        ip.close()
    return print("IRC input generated for " + xyzfile[:-4])
        
def launcher(uplist,rootdir,binfolder):
    inp_file_job_ids = []
    
    n = 1
    
    for xyzfile in uplist:
        reduced_filename=xyzfile[:-4]+"_IRC"
        filename_forward=xyzfile[:-4]+"_IRCforward"+".gjf"
        filename_reverse=xyzfile[:-4]+"_IRCreverse"+".gjf"
        output_forward=xyzfile[:-4]+"_IRCforward"+".log"
        output_reverse=xyzfile[:-4]+"_IRCreverse"+".log"
        IRC_inputgenerator(xyzfile,filename_forward,"forward")
        IRC_inputgenerator(xyzfile,filename_reverse,"reverse")
        
        with open(reduced_filename+".sub","w") as gsub:
            gsub.write('#!/bin/sh\n')
            gsub.write(f'#SBATCH --job-name={reduced_filename}\n')
            gsub.write('#SBATCH --cpus-per-task=12\n')
            gsub.write(f'#SBATCH --output={reduced_filename}.logfile\n')
            gsub.write('#SBATCH --time=15:00:00\n')
            gsub.write('#SBATCH --partition=zen4\n')
            gsub.write('#SBATCH --mem-per-cpu=5GB\n')
            gsub.write('\n')
            gsub.write('# Loading modules\n')
            gsub.write('module load Gaussian/G16.A.03-intel-2022a\n')  # Adjust based on the available Gaussian module
            gsub.write('\n')
            gsub.write('# Setting up Gaussian environment\n')
            gsub.write('export GAUSS_SCRDIR=$TMPDIR\n')  # Temporary directory for Gaussian scratch files
            gsub.write('mkdir -p $GAUSS_SCRDIR\n')
            gsub.write('#Launching calculation\n')
            gsub.write('export PATH={binfolder}:$PATH\n')
            gsub.write('dos2unix {filename_forward}\n')
            gsub.write('dos2unix {filename_reverse}\n')
            gsub.write(f'g16 < {filename_forward} > {output_forward} &\n')
            gsub.write(f'g16 < {filename_reverse} > {output_reverse} &\n')
            gsub.write(f'wait')
            gsub.write('\n')

    
        sbatch_command = f"sbatch {reduced_filename}.sub"
        
        
        result = subprocess.run(
            sbatch_command,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True
        )
        
        if result.returncode == 0:
        # Extract the job ID from the sbatch output
            job_id_match = re.search(r'(\d+)', result.stdout)
            if job_id_match:
                job_id = job_id_match.group(1)
                inp_file_job_ids.append(job_id)  # Collect job IDs for inp_file jobs
                n += 1  # Increment the counter
                print(f"Submitted job {n}")
            else:
                print(f"Failed to extract job ID for {sbatch_command}. Output: {result.stdout}")
        else:
            print(f"Failed to submit job for {sbatch_command}: {result.stderr}")
    
    if inp_file_job_ids:
        dependency_str = ":".join(inp_file_job_ids)
        extractor_script = os.path.join(rootdir,'2_StationaryPoints_calculator.sub')

        if not os.path.exists(extractor_script):
            raise FileNotFoundError(f"Extractor script not found: {extractor_script}")

        dependency_command = [
            "sbatch",
            f"--dependency=afterany:{dependency_str}",
            extractor_script
        ]

        result = subprocess.run(dependency_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        if result.returncode == 0:
            print(f"Extractor job submitted successfully: {result.stdout}")
        else:
            print(f"Failed to submit extractor job: {result.stderr}")
    else:
        print("No jobs were submitted, skipping dependency job submission.")
    

launcher(IRClist,rootdir,binfolder)


    



