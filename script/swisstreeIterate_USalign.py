import os, sys
import glob
import subprocess
import itertools
import csv

# USalign -mm 5 fully-non-sequential
# for each algorithm to do all vs. all comparison for a set of swiss tree proteins
overwrite = 1
cluster = ""

# enumerate and then each comparison
fold = sys.argv[1] # ST001 - ST011
folder = f'SwissTree{cluster}/' + fold

struct_list = glob.glob(folder + '/structs/*')
pair_list = list(itertools.combinations(struct_list,2))

os.system(f'mkdir -p {folder}/USalign')
if not os.path.exists(f'{folder}/USalign/tmscore.csv') or overwrite:
    with open(f'{folder}/USalign/tmscore.csv','w') as f:
        csv.writer(f).writerow(['Tree', 'Query', 'Target', 'Model', 'RMSD', 'L_align', 'SO_Identity', 'AA_Identity', 'AA_Identity_full', 'TM-score'])
for pair in pair_list:
    q_pdb,t_pdb = pair
    q = os.path.basename(q_pdb).replace('.pdb','')
    t = os.path.basename(t_pdb).replace('.pdb','')
    if not os.path.exists(f'{folder}/USalign/{q}_{t}.usalign'):
        os.system(f'USalign {q_pdb} {t_pdb} -mm 5 -outfmt 2 > {folder}/USalign/{q}_{t}.usalign')
    stats = subprocess.run(f"sed -n '2p' {folder}/USalign/{q}_{t}.usalign", stdout=subprocess.PIPE, shell=True, text=True).stdout.strip().split()
    queryLen, targetLen, TM1, TM2, rmsd, lalign, aa_identity = stats[8], stats[9], stats[2], stats[3], stats[4], stats[10], stats[7] #identity is calculated from matched nongapped positions
    aa_identity_full = round(float(aa_identity) * int(lalign)) / sum([int(queryLen),int(targetLen)]) * 2
    identity = int(lalign) / sum([int(queryLen),int(targetLen)]) * 2
    tmscore = max(TM1,TM2)
    with open(f'{folder}/USalign/tmscore.csv','a') as f:
        csv.writer(f).writerow([fold, q, t, 'USalign', rmsd, lalign, identity, aa_identity, aa_identity_full, tmscore])
