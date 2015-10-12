import os, re, subprocess
import numpy as np


def main():

    pdbs = ""
    asa_dir = "asa_files/"
    rsa_dir = "virus_rsa_files/"
    rsa_script = "mod_calc_rsa.py"
    asa_files = []
    if os.path.exists(rsa_dir)  == False:
        os.mkdir(rsa_dir)
    #print os.listdir(asa_dir)
    for asa_file in os.listdir(asa_dir):
        if asa_file.endswith(".txt"):
            asa_files.append(asa_file) 
    for file in asa_files:
        #print asa_files
        name_parts = re.split("_",file)
        print name_parts
        name = name_parts[0] + "_" + name_parts[1] + "_" +  name_parts[2] + "_rsa.txt"
        rsa_command = "python " + rsa_script + " " + asa_dir + file + " " + rsa_dir + name 
        print rsa_command
        subprocess.call(rsa_command, shell = True )
        
        
        #print rsa_command
    #print pdb_files


if __name__ == "__main__":
    main()