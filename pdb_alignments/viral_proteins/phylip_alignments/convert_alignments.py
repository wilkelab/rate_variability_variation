import protein_toolbox_helper as ph
import os, re


def main():    
    
    pdb_dir = "../structures/viral_structures/"     
    align_dir = "../old_viral_pdb_alignments/"
    dir_list = os.listdir(align_dir)
    align_list = []
    for a in dir_list:
        if(a.endswith(".fasta")):
            align_list.append(a)  

    print align_list
    for alignment in align_list:
        name_parts = re.split("\.", alignment)
        print name_parts
        
        protein, chain = re.split("_",  name_parts[0])
        pdb_file = pdb_dir + "RepairPDB_" + name_parts[0] + ".pdb"
        if chain == "AB":
            real_chain = "X"
        else:
            real_chain = chain
        pdb_seq = ph.get_seq_from_pdb(pdb_file, real_chain)
        seqs, headers = ph.get_sequences(align_dir + alignment)
        print seqs
        seqs.pop(0)
        headers.pop(0)
        seqs.insert(0, pdb_seq)
        headers.insert(0, ">" + protein + "_" + chain)
        print headers
        out_fasta = "temp_align.txt"
        ph.align_seqs_mafft(seqs, headers, out_fasta)
        
        out_filename = name_parts[0] + "_alignment.phy"
        fasta_filename = name_parts[0] + ".fasta"
        ph.convert_alignment_format(out_fasta, "fasta", "phylip-relaxed", out_filename, all_caps = True)       
        ph.convert_alignment_format(out_fasta, "fasta", "fasta", fasta_filename, all_caps = True)

main()