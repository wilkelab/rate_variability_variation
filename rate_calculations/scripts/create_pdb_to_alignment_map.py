#!/usr/bin/python
import protein_toolbox_helper as ph
import subprocess, os, sys, re

'''
Last Edited By: Eleisha Jackson (February 2, 2015)
Description: This script create a map. This map maps a position in the alignment to each position in the pdb sequence. Must have MAFFT installed to use this script.
'''

def create_pdb_alignment_map(seqs, pdb_res_nums, map_filename):
	'''
	Creates a file that maps a pdb seqeunce onto a given alignments. 
	Args: 
		seqs: A list of aligned sequences 
		pdb_res_nums: A list of residue numbers from the pdb file
		map_filename: The file name to write the map to
	
	Returns:
		Does not return anything but creates a file with a map
	'''
	
	alignment_length = len(seqs[0]) 
	aligned_pdb_seq = seqs[0]
	#Make map list
	map_outfile = open(map_filename, "w")
	map_outfile.write("align_pos\tpdb_pos\tpdb_aa\n")
	res_counter = 0
	for j in xrange(0, alignment_length):
		map_outfile.write(str(j+1) + "\t")
		if(aligned_pdb_seq[j]!= "-"):
			map_outfile.write(str(pdb_res_nums[res_counter]) + "\t" + aligned_pdb_seq[j])
			res_counter = res_counter + 1
		else:
			map_outfile.write("NA" + "\t" + "-")
		map_outfile.write("\n")

def main():
	if len( sys.argv ) != 6:
		print '''      Usage:      '''
		print "     ", sys.argv[0], "<input alignment>", "<pdb file> " "<output map file> " "<4 Letter pdb> " "<chain>" 
	else:
		align_file = sys.argv[1]
		pdb_file = sys.argv[2]
		map_filename = sys.argv[3]
		pdb = sys.argv[4]
		chain = sys.argv[5]
		
		#fileparts = re.split("/", pdb_file)
		#parts = fileparts[-1]
		#pdb_parts = re.split("_", re.split("\.", parts)[0])
		#pdb = pdb_parts[0]
		#chain = pdb_parts[1]
		#print "chain"
		(pdb_seq, pdb_res_nums) = ph.get_info_from_pdb(pdb_file, chain)	#Get the residue list from pdb
		
		#Get unaligned sequences from the alignment
		(all_seqs, all_headers) = ph.get_unaligned_seqs(align_file)	

		#Then add this to the top of the file
		align_outfile = "temp.txt"
		all_seqs.insert(0, pdb_seq)
		all_headers.insert(0, ">" + pdb + "_" + chain)
		
		#Then write all of the sequences to a file, Align the sequences in a file
		ph.align_seqs_mafft(all_seqs, all_headers, align_outfile)		
		
		#Extract the sequences from the aligned file	
		(seqs, headers) = ph.get_sequences(align_outfile)
		
		#Count the number of non-gap residues and then compare to the residue number length
		if (len(pdb_res_nums) != ( len(seqs[0]) - seqs[0].count("-"))):
			print "Number of residues extracted does not match number of residues in pdb seq"
			print "Number of residues: ",len(pdb_res_nums)
			print "Length of sequence: ", seqs[0].count("-") 
			sys.exit()
	
		create_pdb_alignment_map(seqs, pdb_res_nums, map_filename)
		os.unlink(align_outfile)
if __name__ == "__main__":
	main()