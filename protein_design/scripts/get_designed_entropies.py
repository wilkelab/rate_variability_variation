import os, re, subprocess
import numpy as np
'''
Description: This is a script that calculates all of the entropies for the designed proteins. 

'''
def main():
	sequence_files = []
	seq_dir = "../designed_sequences/"
	entropy_dir = "../entropies/"
	entropy_script = "calculate_sequence_entropy.py"
	
	if os.path.exists(entropy_dir) == False: 
		os.mkdir(entropy_dir)
	
	for seq_file in os.listdir(seq_dir): #Get a list of all of the sequenced files 
		if seq_file.endswith(".fasta"):
			sequence_files.append(seq_file)
	print sequence_files
	
	for align_file in sequence_files: #For each sequence file, calculate the entropy at each site
		pdb_name = re.split("_", align_file )[0] + "_" + re.split("_", align_file)[1]	
		entropy_filename = entropy_dir + pdb_name + "_designed_entropies.txt" 
		command = "python " + entropy_script + " " + seq_dir + align_file  + " "  + entropy_filename
		subprocess.call(command, shell = True )

	out = open(entropy_dir + "all_designed_entropies.csv", "w") 
	out.write("pdb_id,designed_entropy\n")#Write a header to the file		

	entropy_files = []
	for file in os.listdir(entropy_dir):
		if file.endswith("entropies.txt"):
			entropy_files.append(file)

	for file in entropy_files: #For each file get the entropy data and write it to the file with all the entropies
		pdb = re.split("_", file )[0]
		chain = re.split("_", file)[1]
	
		#Read this entropy into a numpy array as string
		entropy_data = np.genfromtxt(entropy_dir + file, dtype = str, skip_header = 1, delimiter = "\t")
		rows, cols = entropy_data.shape #Get the shape
		
		for i in xrange(0, rows):
			out.write(pdb + "," + entropy_data[i][cols-1] + "\n") #Write the last column with are the entropies	
	out.close()

if __name__ == "__main__":
	main()