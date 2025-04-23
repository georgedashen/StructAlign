import os, sys
import glob
import subprocess
import itertools
import csv

import numpy as np
from deepblast.utils import load_model
from deepblast.dataset.utils import states2alignment

load = True
if load:
    model = load_model('models/deepblast-v3.ckpt', 'models/prot_t5_xl_uniref50', device='cuda')

# for each algorithm to do all vs. all comparison for a set of swiss tree proteins
overwrite = 1
cluster = '_cluster'

# enumerate and then each comparison
fold = sys.argv[1] # ST001 - ST011
os.chdir('../')
folder = f'SwissTree{cluster}/' + fold
struct_list = glob.glob(folder + '/structs/*')
pair_list = list(itertools.combinations(struct_list,2))

os.system(f'mkdir -p {folder}/deepblast')
if not os.path.exists(f'{folder}/deepblast/tmscore.csv') or overwrite:
    with open(f'{folder}/deepblast/tmscore.csv','w') as f:
        csv.writer(f).writerow(['Tree', 'Query', 'Target', 'Model', 'RMSD', 'L_align', 'SO_Identity', 'AA_Identity', 'AA_Identity_full','TM-score'])
for pair in pair_list:
    q_pdb,t_pdb = pair
    q = os.path.basename(q_pdb).replace('.pdb','')
    t = os.path.basename(t_pdb).replace('.pdb','')
    if not os.path.exists(f'{folder}/deepblast/{q}_{t}.deepblast.ali') and not os.path.exists(f'{folder}/deepblast/{q}_{t}.deepblast.tmalign'):
        result = subprocess.run(f"./pdb2fasta {q_pdb}", stdout=subprocess.PIPE, shell=True, text=True)
        x = result.stdout.strip().split('\n')[1]
        result = subprocess.run(f"./pdb2fasta {t_pdb}", stdout=subprocess.PIPE, shell=True, text=True)
        y = result.stdout.strip().split('\n')[1]
        pred_ali = model.align(x,y)
        x_aligned, y_aligned = states2alignment(pred_ali, x, y)
        a1 = '>'+q_pdb.replace('.pdb','')
        a2 = '>'+t_pdb.replace('.pdb','')
        subprocess.run(f"newlines=('{a1}' '{x_aligned}' '{a2}' '{y_aligned}') && printf '%s\n' ${{newlines[@]}} > {folder}/deepblast/{q}_{t}.deepblast.ali.fasta",shell=True,executable='/bin/bash')
        os.system(f"sed '/^>/d' {folder}/deepblast/{q}_{t}.deepblast.ali.fasta > {folder}/deepblast/{q}_{t}.deepblast.ali")


    if not os.path.exists(f'{folder}/deepblast/{q}_{t}.deepblast.tmalign'):
        os.system(f'TMalign {q_pdb} {t_pdb} -I {folder}/deepblast/{q}_{t}.deepblast.ali.fasta > {folder}/deepblast/{q}_{t}.deepblast.tmalign')

    result = subprocess.run(f"grep 'Aligned length' {folder}/deepblast/{q}_{t}.deepblast.tmalign", stdout=subprocess.PIPE, shell=True, text=True)
    stats = result.stdout.strip().split()
    stats = [i.replace(',','') for i in stats]
    lalign, rmsd, aa_identity = stats[2], stats[4], stats[6]

    result = subprocess.run(f"grep 'TM-score=' {folder}/deepblast/{q}_{t}.deepblast.tmalign | cut -d' ' -f2", stdout=subprocess.PIPE, shell=True, text=True)
    tmscore = max(result.stdout.strip().split('\n'))
    
    #result = subprocess.run(f"grep 'Aligned length' {folder}/deepblast/{q}_{t}.deepblast.tmalign", stdout=subprocess.PIPE, shell=True, text=True)
    #identity = result.stdout.strip().split('=')[4].replace(' ','')
    queryLen = subprocess.run(f"grep 'Length of Chain_1' {folder}/deepblast/{q}_{t}.deepblast.tmalign", stdout=subprocess.PIPE, shell=True, text=True).stdout.strip().split(' ')[3]
    targetLen = subprocess.run(f"grep 'Length of Chain_2' {folder}/deepblast/{q}_{t}.deepblast.tmalign", stdout=subprocess.PIPE, shell=True, text=True).stdout.strip().split(' ')[3]
    identity = int(lalign) / sum([int(queryLen),int(targetLen)]) * 2
    aa_identity_full = round(float(aa_identity) * int(lalign)) / sum([int(queryLen),int(targetLen)]) * 2
    
    with open(f'{folder}/deepblast/tmscore.csv','a') as f:
        csv.writer(f).writerow([fold, q, t, 'deepblast', rmsd, lalign, identity, aa_identity, aa_identity_full, tmscore])
