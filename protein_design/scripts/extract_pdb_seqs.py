import os, re
import protein_toolbox_helper as ph
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB.PDBParser import PDBParser
from Bio import AlignIO
from Bio import PDB

def extract_designed_seqs(structure_path, dir, number_designed_proteins, seq_dir):
	seqs = []
	headers = []	
	pdbs = os.listdir(structure_path + dir)
	assert len(pdbs) == number_designed_proteins
	print dir, ",", len(pdbs)
	#print os.listdir(structure_path + dir)	
	for pdb in pdbs:	
		protein_parts = re.split("_", pdb)
		pdb_id = protein_parts[0]
		chain = protein_parts[1]
		seq = ph.get_seq_from_pdb(structure_path + dir + "/" +  pdb, chain)
		seqs.append(seq)
		headers.append(">" + pdb)
	seq_filename = structure_path + seq_dir + dir + "_designed_seqs.fasta" 
	ph.write_seq_to_file(seqs, headers, seq_filename)	

def main():
	seq_dir = "../designed_sequences/" 
	if os.path.exists(seq_dir) == False:
		os.mkdir(seq_dir)
	number_designed_proteins = 500
	structure_path = "../designed_structures/"
	print "something"
	protein_dirs = os.listdir(structure_path)
	print protein_dirs
	for dir in protein_dirs:
		print "In Dir: ", dir
		extract_designed_seqs(structure_path, dir, number_designed_proteins, seq_dir)

if __name__ == "__main__":
	main()
