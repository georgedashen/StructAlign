import os, sys
import glob
import subprocess
import itertools
import csv

# for each algorithm to do all vs. all comparison for a set of swiss tree proteins
overwrite = 1
cluster = ''

# enumerate and then each comparison
fold = sys.argv[1] # ST001 - ST011
folder = f'SwissTree{cluster}/' + fold # can switch to SwissTree_cluster
struct_list = glob.glob(folder + '/structs/*')
pair_list = list(itertools.combinations(struct_list,2))

os.system(f'mkdir -p {folder}/TMalign')
if not os.path.exists(f'{folder}/TMalign/tmscore.csv') or overwrite:
    with open(f'{folder}/TMalign/tmscore.csv','w') as f:
        csv.writer(f).writerow(['Tree', 'Query', 'Target', 'Model', 'RMSD', 'L_align', 'SO_Identity', 'AA_Identity', 'AA_Identity_full', 'TM-score'])
for pair in pair_list:
    q_pdb,t_pdb = pair
    q = os.path.basename(q_pdb).replace('.pdb','')
    t = os.path.basename(t_pdb).replace('.pdb','')
    if not os.path.exists(f'{folder}/TMalign/{q}_{t}.tmalign'):
        os.system(f'TMalign {q_pdb} {t_pdb} > {folder}/TMalign/{q}_{t}.tmalign')
    result = subprocess.run(f"grep 'Aligned length' {folder}/TMalign/{q}_{t}.tmalign", stdout=subprocess.PIPE, shell=True, text=True)
    stats = result.stdout.strip().split()
    stats = [i.replace(',','') for i in stats]
    lalign, rmsd, aa_identity = stats[2], stats[4], stats[6] #identity is calculated from matched nongapped positions
    queryLen = subprocess.run(f"grep 'Length of Chain_1' {folder}/TMalign/{q}_{t}.tmalign", stdout=subprocess.PIPE, shell=True, text=True).stdout.strip().split(' ')[3]
    targetLen = subprocess.run(f"grep 'Length of Chain_2' {folder}/TMalign/{q}_{t}.tmalign", stdout=subprocess.PIPE, shell=True, text=True).stdout.strip().split(' ')[3]
    result = subprocess.run(f"grep 'TM-score=' {folder}/TMalign/{q}_{t}.tmalign | cut -d' ' -f2", stdout=subprocess.PIPE, shell=True, text=True)
    aa_identity_full = round(float(aa_identity) * int(lalign)) / sum([int(queryLen),int(targetLen)]) * 2
    identity = int(lalign) / sum([int(queryLen),int(targetLen)]) * 2
    tmscore = max(result.stdout.strip().split('\n'))
    with open(f'{folder}/TMalign/tmscore.csv','a') as f:
        csv.writer(f).writerow([fold, q, t, 'TMalign', rmsd, lalign, identity, aa_identity, aa_identity_full, tmscore])
