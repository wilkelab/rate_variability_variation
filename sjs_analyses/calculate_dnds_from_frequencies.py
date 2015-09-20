# SJS
# This script computes a site-wise dN/dS, based on Spielman and Wilke (2015), MBE for each alignment based on amino-acid frequencies. Mutation rates are assumed equal.
# Dependencies: numpy, biopython, pyvolve

from compute_dnds_from_mutsel import *
from pyvolve import Genetics
g = Genetics()
from Bio import AlignIO
import os


def compute_dnds(d):
    c = dNdS_from_MutSel(d)
    dnds = None
    # Assertion error will be raised in dnds calculation if no evolution has occurred. These are uninformative sites.
    try:
        dnds = c.compute_dn()
        if np.isnan(dnds) or np.isinf(dnds):
            dnds = "NA"
    except AssertionError:
        dnds = "NA"
    
    assert(dnds != None), "dnds not computed."
    return dnds


def derive_site_dnds(records):
    
    all_dnds = []
    for x in range(len(records[0])):
        # Collect amino-acid frequencies, ignoring any gaps or noncanonical states.
        column = {}
        column_raw = list(records[:,x])
        for aa in g.amino_acids:
            column[aa] = float(column_raw.count(aa))
        length = np.sum(column.values())
        
        # for totally uninformative columns
        if length == 0.:
            all_dnds.append(0.)
            continue
        
        for key in column:
            column[key] *= 1./length
        assert(abs(1. - np.sum(column.values())) <= 1e-8), "Frequencies don't sum to 1."

        # Normalize dictionary and convert to codon
        codon_dictionary = {}
        for aminoacid in column:
            i = g.amino_acids.index( aminoacid )
            syn_codons = g.genetic_code[i]
            for syn in syn_codons:
                aafreq = column[aminoacid]
                codon_dictionary[ syn ] = aafreq / len(g.genetic_code[i])

        # Calculate the dN/dS  
        all_dnds.append( compute_dnds(codon_dictionary) )       

    return all_dnds
    
    
    
def grab_all_dnds(input_directory, format, extension, output_directory):
    files = os.listdir(input_directory)
    
    for file in files:
        outfile = output_directory + file.split(extension)[0] + "_dnds_from_frequencies.txt"
        if os.path.exists(outfile):
            continue
            
        print file
        if file.endswith(extension):
            with open(input_directory + file, "r") as f:
                aln = AlignIO.read(input_directory + file, format)
        all_dnds = derive_site_dnds(aln)
        assert(len(all_dnds) == len(aln[0])), "\n Not all dN/dS values were caculated."
        
        with open(outfile, "w") as f:
            f.write("site,dnds\n")
            for i in range(len(all_dnds)):
                f.write(str(i+1)+","+str(all_dnds[i])+"\n")
            
            
                

def main():


    output_directory = "dnds_from_frequencies/"
    enzyme_alignments = "../evol_rates/enzyme_proteins/phylip_alignments/" # phylip format!!
    viral_alignments  = "../evol_rates/viral_proteins/pdb_alignments/"     # fasta format!!

    grab_all_dnds(enzyme_alignments, "phylip-relaxed", ".phy", output_directory)
   # grab_all_dnds(viral_alignments, "fasta", ".fasta", output_directory)
    

main()
