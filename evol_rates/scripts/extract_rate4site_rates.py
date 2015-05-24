import protein_toolbox_helper as ph
import numpy as np
import re, os, sys

'''
Last Edited By: Eleisha Jackson on February 9, 2015
This is a script that extracts the amino acid rates from the individual Rate4Site file for each protein and then places it into a file.
'''

def get_rates(map_file, rate_file):
	'''
	This is a function that extracts and maps the evolutionary rates to the pdb sequence
	
	Args:

		map_file: The file that maps the pdb sequence to the alignment	
		rate_file: The rate4site file
	Returns:
		sites: The sites in the pdb that have been mapped to the alignment	
		aligned_rates: The rates mapped back to sites in the pdb	
	'''
			
	align_sites = []
	rates = []
	amino_acids = []
	ref_seq_array  = []

	input = open(rate_file)
	data = input.readlines()
	data = data[11:]
	data.pop()
	data.pop()

	for line in data: #Get the evolutionary rates from the rate4site file
		site = line[2:8].strip()
		amino = line[9:12].strip()
		rate = line[11:19].strip()
		align_sites.append(site)
		amino_acids.append(amino)
		rates.append(rate)	

	sites = extract_sites(map_file)
	aligned_rates = []
	counter = 0
	ref_seq_array = np.genfromtxt(map_file, delimiter = "\t", dtype = "string", skip_header = 1, usecols = (2)) #Read in the data
	#print ref_seq_array
	#Align the calculated rates to the proper amnio acid in the pdb using the pdb_to_alignment_map created earlier (in match_frequencies.py)
	
	for j in xrange(0, len(ref_seq_array)):		
		if(ref_seq_array[j] == '-'):
			#print str(j+1), ref_seq_array[j], sites[j], "NA"
			aligned_rates.append("NA")
		else:
			#print str(j+1), ref_seq_array[j], sites[j]
			aligned_rates.append(rates[counter])
			counter = counter + 1
	return sites, aligned_rates

def extract_sites(map_file):
	'''
	This is a function extracts the site information from the mapped file
	
	Args:
		pdb: The pdb name
		map_file: The file that maps the pdb sequence to the alignment
		rate_file: the rate4site file
	Returns:
		sites: The sites in the pdb that have been mapped to the alignment		
	'''

	sites = []
	site_data = np.genfromtxt(map_file, delimiter = "\t", dtype = "string", skip_header = 1, usecols = (0)) #Read in the data
	sites = []
	for s in site_data:
		sites.append(s.strip()) #String off the new lines and append the site to the sites list
	return sites

def main():
	if len( sys.argv ) != 4:
		print '''      Usage:      '''
		print "     ", sys.argv[0], "<input rate4site file>",  "<map file for pdb to alignment> " "<output file>"
								
	else:
		rate_file = sys.argv[1] #rate4site file
		map_file = sys.argv[2]  #Name of the file that maps the pdb sequence to the alignment
		outfile = sys.argv[3] #outfile name
		
		out = open(outfile, "w")
		out.write('"site","rate"\n')

		try:
			sites, rates = get_rates(map_file, rate_file) #Extract the rates
		except IOError:
			print "The File: " + pdb + " is not there."

		if (len(sites) != len(rates)):
			print "Sites not mapped correctly to Rates!"
			print "Length of Sites: ", len(sites)
			print "Length of Rates: ", len(rates)

		j = 0
		while ( j < len(sites)): #For each site
			#print j, str(sites[j]), str(rates[j]) + "\n"	
			if (sites[j] != "NA"): #If the site is NOT NA (Meaning the the site is mapped to a site in the pdb...
				out.write(str(sites[j]) + ',' + str(rates[j]) + "\n") #Write out the rate information for that site
			
				j = j + 1
			else: #else do not it
				j = j + 1
		out.close()		

if __name__ == "__main__":
	main()
