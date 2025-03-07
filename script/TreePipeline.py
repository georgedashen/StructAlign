# conda install -c bioconda fastme
# pip install ete3 colour wget statsmodels matplotlib seaborn

import glob
import ete3
import wget 
import pandas as pd
import os
import time
import json
from src import AFDB_tools, treescore , foldseek2tree
import colour
import numpy as np
import toytree
import pickle
from matplotlib import pyplot as plt
from Bio import SeqIO
import toytree
import toyplot.svg
import csv

cluster = ''
draw = True

# download from existing identifiers.txt
examples = f'../SwissTree{cluster}/*/'
download = False
for folder in glob.glob(examples):
    if not download:
        continue
    idfile = glob.glob(folder + '*identifiers.txt')[0]
    with open(idfile) as idin:
        identifiers = [l for l in idin]
    try:
        os.mkdir(folder + 'structs')
    except:
        for file in glob.glob(folder + 'structs/*.pdb'):
            os.remove(file)
    print('done setup')
    
    if download == True:
        notfound = [ AFDB_tools.grab_struct(uniID, folder + '/structs')  for uniID in identifiers ]
        notfound = [ i for i in notfound if i]
        with open(folder+'unmapped.txt', 'w') as simout:
            simout.write(json.dumps(notfound))

# download from only newick files
download = False
update = 'pruned'
for t in glob.glob(f'../SwissTree{cluster}/*/*updated*tree.nhx'):
    if not download:
        continue
    #print(t)
    try:
        os.mkdir(os.path.dirname(t)+'/structs/')
    except:
        print('already exists')
    outfolder = os.path.dirname(t)+'/structs/'
    
    if download == True:
        if update not in t :
            print(t)
            tre = toytree.tree(t )
            notfound = [ AFDB_tools.grab_struct(uniID, outfolder)  for uniID in tre.get_tip_labels() ]
            notfound = [ i for i in notfound if i]
            with open(os.path.dirname(t)+'/unmapped.txt', 'w') as simout:
                simout.write(json.dumps(notfound))


#prune not found
overwrite = False
for t in glob.glob(f'../SwissTree{cluster}/*/*updated*tree.nhx'):
    if 'pruned' not in t:
        if os.path.exists(t+'pruned.nhx') and not overwrite:
            continue
        with open(os.path.dirname(t)+'/unmapped.txt') as missing:
            missingIDs = json.loads(missing.read())
        #print(missingIDs)
        #print(t)
        tree = ete3.Tree(t)
        #print(len(tree))
        keep = [ n for n in tree.traverse() if n.name not in missingIDs]
        tree.prune(keep, preserve_branch_length=True)
        keep=[ n for n in tree.get_leaves() if n.name!='']
        tree.prune(keep, preserve_branch_length=True)
        #print(len(tree))
        #print(t+'pruned.nhx')
        with open(t+'pruned.nhx' , 'w' ) as treeout:
            treeout.write(tree.write(format=1))


# tree building
overwrite = False
for folder in glob.glob(f'../SwissTree{cluster}/*/structs'):
    #print(folder)
    outfolder = os.path.dirname(folder)
    if os.path.exists(outfolder+'/allvsall.tsv') and not overwrite:
        continue
    cmd = 'foldseek/foldseek' + ' easy-search ' + folder + ' ' + folder +' '+ outfolder+'/allvsall.tsv ' + outfolder+"/tmp --format-output 'query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,lddt,lddtfull,alntmscore' -a -e inf --exhaustive-search"
    os.system(cmd)


def standard_treedraw( tre, sizes= None , colors= None ,fixed_order=None, fixed_position=None , ts = None,  save_file = False  , tiplabels = None):
    if tiplabels is None:
        tiplabels = tre.get_tip_labels()
    canvas, axes, mark = tre.draw(  
        ts = ts,
        node_sizes=sizes,
        node_colors=colors,
        tip_labels_align=True,
        scalebar=True,
        fixed_order=fixed_order, 
        fixed_position=fixed_position,
        tip_labels=tiplabels,
        tip_labels_style={
            "fill": "#262626",
            "font-size": "9px"}
    )
    if save_file:
        toyplot.svg.render(canvas, save_file)


#get trees
swisstreepath = f'../SwissTree{cluster}/*/*nhxpruned.nhx'
swiss_trees = glob.glob(swisstreepath)
#print(swiss_trees)
swiss_toytrees = {}

discretized_levels = 10
maxnodesize = 15
overwrite = False
for t in swiss_trees:
    if os.path.exists(os.path.dirname(t)+'/swiss_toytree.pkl') and not overwrite:
        continue
    print(t)
    treecheck = True
    try:
        tre = toytree.tree(t , format =0 )
    except:
        print('tree err' , t)
        treecheck = False
    if treecheck == True:
        swiss_toytrees[t] = { 'tree' : tre }
        #print(tre)
        swiss_toytrees[t]['labels'] = tre.get_tip_labels()
        swiss_toytrees[t]['coords'] = tre.get_tip_coordinates()
        swiss_toytrees[t]['path'] = t
        labels = swiss_toytrees[t]['labels']
        unidf = AFDB_tools.grab_entries(labels)
        lineages = treescore.make_lineages(unidf)
        #discretized_levels = 10
        #maxnodesize = 15
        taxlabels = dict(zip(unidf['query'] ,unidf.Organism) )
        #species_mapper = 
        tipnames = tre.get_tip_labels()
        taxnames = [ i+ ' '+taxlabels[i] if i in taxlabels else i for i in tipnames]
        #red = colour.Color('blue')
        #blue = colour.Color('red')
        tre = treescore.label_leaves(tre,lineages)
        #color_vals = list(red.range_to(blue, discretized_levels+1))
        overlap = treescore.getTaxOverlap(tre.treenode)
        #treevals = [ node.score for node in tre.treenode.traverse()]
        #print(treevals)
        print('finalscore' , tre.treenode.score)
        swiss_toytrees[t]['score'] = tre.treenode.score
        with open(os.path.dirname(t)+'/swiss_toytree.pkl', 'wb') as f:
            pickle.dump(swiss_toytrees[t], f)
        #treevals = tre.get_node_values('size', True, True)
        #maxval = np.amax(treevals)
        #minval =0
        #bins = [minval + i*(maxval-minval)/discretized_levels for i in range(discretized_levels+1)]
        #inds = np.digitize(treevals, bins)
        #colors = [ color_vals[i-1].hex_l for i in list(inds)]
        #sizes = [ (i / discretized_levels) * maxnodesize for i in list(inds) ]
        #standard_treedraw(tre , sizes= sizes , colors=colors , tiplabels = taxnames , save_file = t + '.svg')
        #standard_treedraw(tre  , sizes= sizes , colors=colors , tiplabels = taxnames )        



## root tree by MAD

import subprocess, shlex
def madroot(t):
    args = 'madroot/mad ' + t
    subprocess.run(shlex.split(args))
    return t+ '.rooted'

## calculate distance
import traceback
#print(swiss_toytrees.keys())

import seaborn
treescores_st = []
treescores_str = []
treescores_struct = []
rfdistances = []

overwrite = False
if not os.path.exists(f'../SwissTree{cluster}.congruence') or overwrite:
    with open(f'../SwissTree{cluster}.congruence', 'w') as f:
        csv.writer(f).writerow(['Tree', 'Model', 'Size', 'Congruence_score'])

if not os.path.exists(f'../SwissTree{cluster}.rf') or overwrite:
    with open(f'../SwissTree{cluster}.rf', 'w') as f:
        csv.writer(f).writerow(['Tree', 'Model', 'Size', 'RF_distance','max_rf'])

#for t in swiss_trees:
for t in swiss_trees:
    folder = os.path.dirname(t)
    swisst = pickle.load(open(folder+'/swiss_toytree.pkl','rb'))
    tmres = folder +'/allvsall.tsv'
    res = pd.read_table(tmres, header = None ,delim_whitespace=True)
    res[0] = res[0].map(lambda x :x.replace('.pdb', ''))
    res[1] = res[1].map(lambda x :x.replace('.pdb', ''))
    ids = list( set(list(res[0].unique()) + list(res[1].unique())))
    ids = list(set(ids).intersection(set(swisst['labels'])))
    #labelorder = [ l for l in swisst['labels'] if l in ids ]
    #positions = [ k for i,k in enumerate(list(swisst['coords'][:,1])) if swisst['labels'][i] in ids ]
    pos = { protid : i for i,protid in enumerate(ids) }
    self_dists = res[res[0] == res[1]]
    self_distmap = dict(zip(self_dists[0] , self_dists[2] ) )    
    kernel_distmat = np.zeros((len(pos), len(pos)))
    distmat = np.zeros((len(pos), len(pos)))
    for idx,row in res.iterrows():
        kernel_distmat[pos[row[0]] , pos[row[1]]] = foldseek2tree.kernelfun(self_distmap[row[0]] , self_distmap[row[1]] , row[2])
        distmat[pos[row[0]] , pos[row[1]]]= row[2] #2:fident
    distmat = 1- distmat
    distmat_txt_kernel = foldseek2tree.distmat_to_txt( ids , distmat , t.replace('nhx','distmat') )
    out_tree_kernel = foldseek2tree.runFastme( 'fastme ' , distmat_txt_kernel )
    out_tree_kernel = foldseek2tree.postprocess(out_tree_kernel, out_tree_kernel.replace('.txt','_PP.txt'))
    out_tree_kernel = madroot(out_tree_kernel)
    tre = toytree.tree(out_tree_kernel)
    swiss_tre = madroot(t)
    swiss_tre = toytree.tree(t)
    
    #tipnames = tre.get_tip_labels()
    unidf = AFDB_tools.grab_entries(ids)
    lineages = treescore.make_lineages(unidf)
    #discretized_levels = 10
    #maxnodesize = 15
    
    #taxlabels = dict(zip(unidf['query'] ,unidf.Organism) )
    #taxnames = [ i+ ' '+taxlabels[i] if i in taxlabels else i for i in tipnames]
    tre = treescore.label_leaves(tre,lineages)
    swiss_tre = treescore.label_leaves(swiss_tre,lineages)
    #color_vals = list(red.range_to(blue, discretized_levels+1))
    overlap = treescore.getTaxOverlap(tre.treenode)
    overlap = treescore.getTaxOverlap(swiss_tre.treenode)
    #treevals = tre.get_node_values('size', True, True)
    
    try:
        rf = tre.treenode.robinson_foulds(swisst['tree'].treenode )
        #rf = tre.treenode.robinson_foulds(swiss_tre.treenode )
        #print( 'RF dist:' , rf )
        rfdistances.append(rf)
    
    except:
        print('rferr')
        rf = tre.treenode.robinson_foulds(swisst['tree'].treenode,unrooted_trees=True )
        rfdistances.append(rf)
        #rf = ['','']

            
    print('struct score: ', tre.treenode.score, 'ST score: ' , swisst['tree'].treenode.score , 'ST_rooted score: ' , swiss_tre.treenode.score)

    treescores_struct.append( tre.treenode.score)
    treescores_st.append( swisst['tree'].treenode.score  )
    treescores_str.append( swiss_tre.treenode.score  )

    
    with open(f'../SwissTree{cluster}.congruence', 'a') as f:
        csv.writer(f).writerow([os.path.basename(t).split('_')[0], 'Swiss', len(swisst['tree']), swisst['tree'].treenode.score])
        csv.writer(f).writerow([os.path.basename(t).split('_')[0], 'foldtree', len(tre), tre.treenode.score])

    with open(f'../SwissTree{cluster}.rf', 'a') as f:
        csv.writer(f).writerow([os.path.basename(t).split('_')[0], 'foldtree', len(tre), rf[0], rf[1] ])
    

if draw:
    plt.scatter( x = treescores_struct , y = treescores_st  , c = 'b', ls = '-')
    plt.scatter( x = treescores_struct , y = treescores_str  , c = 'r', ls = '-')
    plt.plot( [0,2000] , [0,2000])
    plt.xlabel('structure')
    plt.ylabel('swisstree')
    plt.show()

