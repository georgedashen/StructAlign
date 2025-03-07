# conda install -c bioconda fastme
# pip install ete3 colour wget statsmodels matplotlib seaborn

import glob
import ete3
import wget 
import pandas as pd
import os, sys
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

method = sys.argv[1] # deepblast, TMalign, deepalign
cluster = '' # or empty string to disable cluster mode
sim = 'Identity'
subfold = ''
if len(sys.argv)>2:
    sim = sys.argv[2] # 'Identity'(n_ident/n_ali) or 'TM-score' or 'fullIdent'(n_ident/length) or RMSD
if len(sys.argv)>3:
    subfold = sys.argv[3]

#get trees
swisstreepath = f'../SwissTree{cluster}/*/*nhxpruned.nhx'
swiss_trees = glob.glob(swisstreepath)

import subprocess, shlex
def madroot(t):
    args = 'madroot/mad ' + t
    subprocess.run(shlex.split(args))
    return t+ '.rooted'


def standard_treedraw( tre, sizes= None , colors= None ,fixed_order=None, fixed_position=None , ts = None,  save_file = '', tiplabels = None):
    if tiplabels is None:
        tiplabels = tre.get_tip_labels()
    canvas, axes, mark = tre.draw(
        width=2400,
        shrink=1500,
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


def TreeDraw(tre, save_file = '', taxnames = None):
    discretized_levels = 10
    maxnodesize = 15
    red = colour.Color('blue')
    blue = colour.Color('red')
    color_vals = list(red.range_to(blue, discretized_levels+1))
    treevals = tre.get_node_values('size', True, True)
    maxval = np.amax(treevals)
    minval = 0
    bins = [minval + i*(maxval-minval)/discretized_levels for i in range(discretized_levels+1)]
    inds = np.digitize(treevals, bins)
    colors = [ color_vals[i-1].hex_l for i in list(inds)]
    sizes = [ (i / discretized_levels) * maxnodesize for i in list(inds) ]
    standard_treedraw(tre , sizes= sizes , colors=colors , tiplabels = taxnames, save_file = save_file)


## calculate distance
import traceback
import seaborn

treescores_st = []
treescores_str = []
treescores_struct = []
rfdistances = []

overwrite = False
draw = True

if not os.path.exists(f'../SwissTree{cluster}.congruence') or overwrite:
    with open(f'../SwissTree{cluster}.congruence', 'w') as f:
        csv.writer(f).writerow(['Tree', 'Model', 'Size', 'Congruence_score'])

if not os.path.exists(f'../SwissTree{cluster}.rf') or overwrite:
    with open(f'../SwissTree{cluster}.rf', 'w') as f:
        csv.writer(f).writerow(['Tree', 'Model', 'Size', 'RF_distance','max_rf'])

#for t in swiss_trees:
overwrite = True
for t in swiss_trees:
    folder = os.path.dirname(t)
    if subfold:
        if os.path.basename(folder) != subfold:
            continue
    swisst = pickle.load(open(folder+'/swiss_toytree.pkl','rb'))
    if not os.path.exists(f'{folder}/{method}/tmscore.csv'):
        print(f'{folder}/{method}/tmscore.csv not found!')
        continue
    res = pd.read_csv(f'{folder}/{method}/tmscore.csv')
    res = res[['Query','Target', sim]]
    ids = list( set(list(res['Query'].unique()) + list(res['Target'].unique())))
    ids = list(set(ids).intersection(set(swisst['labels'])))
    if len(set(swisst['labels'])):
        if len(set(swisst['labels']))==1:
            ids = list(set(swisst['labels']))
            res.loc[len(res)] = [ids[0],ids[0],1]
    else:
        continue
    if sim == 'RMSD':
        res['RMSD'] = 1 - res['RMSD'] / max(res['RMSD']) # sim = 1 - dist
    pos = { protid : i for i,protid in enumerate(ids) }
    distmat = np.eye(len(pos))
    for i,row in res.iterrows():
        if row['Query'] in pos and row['Target'] in pos:
            distmat[pos[row['Query']] , pos[row['Target']]]= row[sim]
    distmat[distmat>1]=1
    if not bool((distmat.transpose()==distmat).all()):
        if method != 'blast': # if results from combinations: (eg.tmalign,deepalign,plmblast,kpax)
            distmat = distmat + distmat.transpose()
            np.fill_diagonal(distmat,1)
        else: # if results from all permutations: (eg. blast)
            distmat[distmat>distmat.transpose()] = distmat.transpose()[distmat>distmat.transpose()]
            
    distmat = 1- distmat
    distmat_txt_kernel = foldseek2tree.distmat_to_txt( ids , distmat , f'{folder}/{method}/{sim}_distmat' )
    out_tree_kernel = foldseek2tree.runFastme( 'fastme ' , distmat_txt_kernel )
    out_tree_kernel = foldseek2tree.postprocess(out_tree_kernel, out_tree_kernel.replace('.txt','_PP.txt'))
    out_tree_kernel = madroot(out_tree_kernel)
    tre = toytree.tree(out_tree_kernel)
    #swiss_tre = madroot(t)
    #swiss_tre = toytree.tree(t)
    
    unidf = AFDB_tools.grab_entries(ids)
    lineages = treescore.make_lineages(unidf)
    
    tre = treescore.label_leaves(tre,lineages)
    #swiss_tre = treescore.label_leaves(swiss_tre,lineages)
    overlap = treescore.getTaxOverlap(tre.treenode)
    taxlabels = dict(zip(unidf['query'] ,unidf.Organism) )
    tipnames = tre.get_tip_labels()
    taxnames = [ i+ ' '+taxlabels[i] if i in taxlabels else i for i in tipnames]

    if draw:
        tipnames = swisst['tree'].get_tip_labels()
        taxnames = [ i+ ' '+taxlabels[i] if i in taxlabels else i for i in tipnames]
        TreeDraw(swisst['tree'], t+'.svg', taxnames)

        tipnames = tre.get_tip_labels()
        taxnames = [ i+ ' '+taxlabels[i] if i in taxlabels else i for i in tipnames]
        TreeDraw(tre, f'{folder}/{method}/{sim}_tree.svg', taxnames)
    
    try:
        rf = tre.treenode.robinson_foulds(swisst['tree'].treenode )
        rfdistances.append(rf)
    
    except:
        print('rferr')
        rf = tre.treenode.robinson_foulds(swisst['tree'].treenode,unrooted_trees=True )
        rfdistances.append(rf)
            
    print('struct score: ', tre.treenode.score, 'ST score: ' , swisst['tree'].treenode.score )

    treescores_struct.append( tre.treenode.score)
    treescores_st.append( swisst['tree'].treenode.score  )

    with open(f'../SwissTree{cluster}.congruence', 'a') as f:
        csv.writer(f).writerow([os.path.basename(t).split('_')[0], f"{method}_{sim.replace('_','-')}", len(tre), tre.treenode.score])

    with open(f'../SwissTree{cluster}.rf', 'a') as f:
        csv.writer(f).writerow([os.path.basename(t).split('_')[0], f"{method}_{sim.replace('_','-')}", len(tre), rf[0], rf[1] ])
    
if not draw:
    plt.scatter( x = treescores_struct , y = treescores_st  , c = 'b', ls = '-')
    plt.plot( [0,2000] , [0,2000])
    plt.xlabel('structure')
    plt.ylabel('swisstree')
    plt.show()

