# SJS
# This script calculates various tree statistics, including tree length, mean root-to-tip, and mean patristic distance
# Dependencies: dendropy, numpy

from dendropy import Tree
from dendropy.calculate import treemeasure # use this line for PD
import os
import re
import numpy as np



def extract_tree_info(file):
    
    t = Tree.get_from_path(file, 'newick')
    
    # tree length
    tree_length = str(t.length())

    # mean root-to-tip distance
    treetips = t.leaf_nodes()
    rtt = []
    for tip in treetips:
        rtt.append( tip.distance_from_root() )
    mean_rtt = str(np.mean(rtt))
    
    # mean patristic distance
    pd = []
    dist = treemeasure.PatristicDistanceMatrix(tree=t)
    for i, t1 in enumerate(t.taxon_namespace):
        for t2 in t.taxon_namespace[i:]:
            d = dist(t1,t2)
            pd.append( float(d) )
    mean_pairwise = str(np.mean(pd))
        
    return tree_length, mean_rtt, mean_pairwise 



def main():


    output_file = "tree_statistics.txt"
    enzyme_directory = "../virus_trees/"
    viral_directory = "../enzyme_trees/"
    enzyme_trees = os.listdir(enzyme_directory)
    viral_trees  = os.listdir(viral_directory)
    
    with open(output_file,"w") as f:
        f.write("name,tree_length,mean_root_to_tip_distance,mean_pairwise_distance\n")
        
    
    for tree in viral_trees:
        print tree
        values = extract_tree_info(viral_directory + tree)
        getname = re.search(r"tree_(.+)\.tre", tree)
        if getname:
            name = getname.group(1)
        with open(output_file,"a") as f:
            f.write(name + "," + values[0] + "," + values[1] + "," + values[2] + "\n")
        
    for tree in enzyme_trees:
        print tree
        values = extract_tree_info(enzyme_directory + tree)
        getname = re.search(r"RAxML_bestTree\.(\w+_\w)_.+", tree)
        if getname:
            name = getname.group(1)
        with open(output_file,"a") as f:
            f.write(name + "," + values[0] + "," + values[1] + "," + values[2] + "\n")
            

main()
