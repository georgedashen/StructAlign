# save the result from the DALI server and download the newick_unrooted.txt file
# put the downloaded file into the corresponding SwissTree family subfolder (e.g.,ST001)

import os,sys
import toytree

fold = sys.argv[1]

# must be unrooted tree
os.system(f'ls SwissTree/{fold}/structs | sort | uniq > SwissTree/{fold}/identifiers_order.txt')
tre=toytree.tree(f'SwissTree/{fold}/dali_unrooted.nhx',tree_format=0)
with open(f'SwissTree/{fold}/identifiers_order.txt','r') as f:
    ids = [i.strip().replace('.pdb','') for i in f]

name_dict={}
for i,n in enumerate(ids):
    name_dict['s'+str(i+1).zfill(3)+'A']=n

tre=tre.set_node_values(feature='name',values=name_dict)
tre.write(f"SwissTree/{folder}/dali_unrooted_rename.nhx")
