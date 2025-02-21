import os, sys
import glob
import pandas as pd
import subprocess
import itertools
import csv

import torch
import alntools as aln
from Bio import SeqIO
import numpy as np

extr=aln.Extractor(bfactor='global')

# for each algorithm to do all vs. all comparison for a set of swiss tree proteins
overwrite = 1
cluster = '_cluster'

# enumerate and then each comparison
fold = sys.argv[1] # ST001 - ST011
#get embeddings
if not os.path.exists(f'../SwissTree_cluster/{fold}/plmblast_pt'):
    os.system('python embeddings.py start ../SwissTree_cluster/{fold}/seq.fasta ../SwissTree_cluster/{fold}/plmblast_pt -embedder pt --gpu -bs 0 && python scripts/dbtofile.py ../SwissTree_cluster/{fold}/plmblast_pt')
if os.path.exists(f'../SwissTree_cluster/{fold}/plmblast_pt.csv'):
    df = pd.read_csv(f'../SwissTree_cluster/{fold}/plmblast_pt.csv')
    name_dict = {}
    for i in range(len(df)):
        name_dict[str(df['queryid'][i])+'.emb'] = df['id'][i].split()[0].replace(':A','')+'.pt'
    filelist = glob.glob(f'../SwissTree_cluster/{fold}/plmblast_pt/*')
    for f in filelist:
        if os.path.basename(f) in name_dict:
            os.system(f'mv {f} {f.replace(os.path.basename(f),name_dict[os.path.basename(f)])}')
    

os.chdir('../')
folder = f'SwissTree{cluster}/' + fold
struct_list = glob.glob(folder + '/structs/*')
if len(struct_list) < 2:
    print('Not enough sequences!')
    sys.exit(1)
pair_list = list(itertools.combinations(struct_list,2))

os.system(f'mkdir -p {folder}/plmblast')
if not os.path.exists(f'{folder}/plmblast/tmscore.csv') or overwrite:
    with open(f'{folder}/plmblast/tmscore.csv','w') as f:
        csv.writer(f).writerow(['Tree', 'Query', 'Target', 'Model', 'RMSD', 'L_align', 'SO_Identity', 'AA_Identity', 'AA_Identity_full', 'TM-score'])
for pair in pair_list:
    q_pdb,t_pdb = pair
    q = os.path.basename(q_pdb).replace('.pdb','')
    t = os.path.basename(t_pdb).replace('.pdb','')
    if not os.path.exists(f'{folder}/plmblast/{q}_{t}.plmblast.ali'):
        seq1 = subprocess.run(f"./pdb2fasta {q_pdb}", stdout=subprocess.PIPE, shell=True, text=True).stdout.strip().split('\n')[1]
        seq2 = subprocess.run(f"./pdb2fasta {t_pdb}", stdout=subprocess.PIPE, shell=True, text=True).stdout.strip().split('\n')[1]
        # get pt first from embeddings.py
        emb1 = torch.load(f'{folder}/plmblast_pt/{q}.pt').numpy()
        emb2 = torch.load(f'{folder}/plmblast_pt/{t}.pt').numpy()
        
        results = extr.embedding_to_span(emb2,emb1)
        row=results.iloc[results.score.argmax()]
        ali=aln.draw_alignment(row.indices, seq1, seq2, output='str').split('\n')
        
        a1 = '>'+q
        a2 = '>'+t

        subprocess.run(f"newlines=('{a1}' '{ali[0]}' '{a2}' '{ali[2]}') && printf '%s\n' ${{newlines[@]}} > {folder}/plmblast/{q}_{t}.plmblast.ali.fasta",shell=True,executable='/bin/bash')
        os.system(f"sed '/^>/d' {folder}/plmblast/{q}_{t}.plmblast.ali.fasta > {folder}/plmblast/{q}_{t}.plmblast.ali")

    if not os.path.exists(f'{folder}/plmblast/{q}_{t}.plmblast.tmalign'):
        os.system(f'TMalign {q_pdb} {t_pdb} -I {folder}/plmblast/{q}_{t}.plmblast.ali.fasta > {folder}/plmblast/{q}_{t}.plmblast.tmalign')
    stats = subprocess.run(f"grep 'Aligned length' {folder}/plmblast/{q}_{t}.plmblast.tmalign", stdout=subprocess.PIPE, shell=True, text=True).stdout.strip().split()
    stats = [i.replace(',','') for i in stats]
    lalign, rmsd, aa_identity = stats[2], stats[4], stats[6]

    result = subprocess.run(f"cat {folder}/plmblast/{q}_{t}.plmblast.ali", stdout=subprocess.PIPE, shell=True, text=True)
    result = result.stdout.strip().split('\n')
    lalign = sum(1 for i in range(len(result[0])) if result[0][i].isupper() and result[1][i].isupper())

    result = subprocess.run(f"grep 'TM-score=' {folder}/plmblast/{q}_{t}.plmblast.tmalign | cut -d' ' -f2", stdout=subprocess.PIPE, shell=True, text=True)
    tmscore = max(result.stdout.strip().split('\n'))
    
    #result = subprocess.run(f"grep 'Aligned length' {folder}/plmblast/{q}_{t}.plmblast.tmalign", stdout=subprocess.PIPE, shell=True, text=True)
    #identity = result.stdout.strip().split('=')[4].replace(' ','')
    queryLen = subprocess.run(f"grep 'Length of Chain_1' {folder}/plmblast/{q}_{t}.plmblast.tmalign", stdout=subprocess.PIPE, shell=True, text=True).stdout.strip().split(' ')[3]
    targetLen = subprocess.run(f"grep 'Length of Chain_2' {folder}/plmblast/{q}_{t}.plmblast.tmalign", stdout=subprocess.PIPE, shell=True, text=True).stdout.strip().split(' ')[3]
    identity = int(lalign) / sum([int(queryLen),int(targetLen)]) * 2
    aa_identity_full = round(float(aa_identity) * int(lalign)) / sum([int(queryLen),int(targetLen)]) * 2
    
    with open(f'{folder}/plmblast/tmscore.csv','a') as f:
        csv.writer(f).writerow([fold, q, t, 'plmblast', rmsd, lalign, identity, aa_identity, aa_identity_full, tmscore])
