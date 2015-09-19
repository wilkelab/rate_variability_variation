#!/usr/bin/python
import protein_toolbox_helper as ph
import subprocess, os, sys, re

'''
Last Edited By: Eleisha Jackson (May 22, 2015)
Description: This script creates a MAFFT alignment. This alignment contains all sequences from the muscle alignment along with the pdb sequences (the first sequence in the alignment).
Must have MAFFT installed to use this script.
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

def main(argv = sys.argv):
	if (len(argv) != 6):
		print "Usage: " 
		print argv[0], "<pdb file>", "<chain>", "<pdb to alignment map filename>", "<alignment>","<output filename>"
	else:
	

		pdb_file = argv[1]
		chain = argv[2]
		map_file = argv[3]
		align_file = argv[4]	
		align_out_prefix = argv[5]
		align_outfile = align_out_prefix + ".fasta"
		pdb_parts = re.split(".", pdb_file)
		pdb = pdb_parts[0]
		print pdb_parts
		#Open the pdb, Get the sequence from the pdb
		#Get the residue list (make a function)
		(pdb_seq, pdb_res_nums) = ph.get_info_from_pdb(pdb_file, chain)	
		
		#Get unaligned sequences from the alignment
		(all_seqs, all_headers) = ph.get_unaligned_seqs(align_file)
		
		#Then add this to the top of the file
		all_seqs.insert(0, pdb_seq)
		all_headers.insert(0, ">" + align_out_prefix)
		
		#Then write all of the sequences to a file, Align the sequences in a file
		ph.align_seqs_mafft(all_seqs, all_headers, align_outfile)	
	
		#Extract the sequences from the aligned file	
		(seqs, headers) = ph.get_sequences(align_outfile)
		
		#Count the number of non-gap residues and then compare to the residue number length
		if (len(pdb_res_nums) != ( len(seqs[0]) - seqs[0].count("-"))):
			print "Number of residues extracted does not match number of residues in pdb seq"
			print "Number of residues: ",len(pdb_res_nums)
			print "Length of sequence: ", seqs[0].count("-") 
	
		#create_pdb_alignment_map(seqs, pdb_res_nums, map_filename)
		print pdb_parts
if __name__ == "__main__":
	main(sys.argv)