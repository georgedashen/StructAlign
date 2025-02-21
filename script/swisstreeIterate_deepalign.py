import os, sys
import glob
import subprocess
import itertools
import csv

# for each algorithm to do all vs. all comparison for a set of swiss tree proteins
overwrite = 1

# enumerate and then each comparison
fold = sys.argv[1] # ST001 - ST011
folder = 'SwissTree_cluster/' + fold
struct_list = glob.glob(folder + '/structs/*')
pair_list = list(itertools.combinations(struct_list,2))

os.system(f'mkdir -p {folder}/deepalign')
if not os.path.exists(f'{folder}/deepalign/tmscore.csv') or overwrite:
    with open(f'{folder}/deepalign/tmscore.csv','w') as f:
        csv.writer(f).writerow(['Tree', 'Query', 'Target', 'Model', 'RMSD', 'L_align', 'SO_Identity', 'AA_Identity', 'AA_Identity_full', 'TM-score'])
for pair in pair_list:
    q_pdb,t_pdb = pair
    q = os.path.basename(q_pdb).replace('.pdb','')
    t = os.path.basename(t_pdb).replace('.pdb','')
    if not os.path.exists(f'{folder}/deepalign/{q}_{t}.deepalign.ali.score'):
        os.system(f'DeepAlign {q_pdb} {t_pdb} -o {folder}/deepalign/{q}_{t}.deepalign.ali -P 0')
    with open(f'{folder}/deepalign/{q}_{t}.deepalign.ali.score','r') as f:
        f.readline()
        scores = f.readline().strip()
    score=scores.split('->')[2].split()
    _, _, qlen, tlen = scores.split('->')[0].split()
    lalign, rmsd, tmscore = score
    identity = int(lalign) / sum([int(qlen),int(tlen)]) * 2
    aa = scores.split('->')[4].split()[0]
    aa_identity = int(aa) / int(lalign)
    aa_identity_full = int(aa) / sum([int(qlen),int(tlen)]) * 2
    with open(f'{folder}/deepalign/tmscore.csv','a') as f:
        csv.writer(f).writerow([fold, q, t, 'deepalign', rmsd, lalign, identity, aa_identity, aa_identity_full, tmscore])

os.system(f'rm {folder}/deepalign/*pdb {folder}/deepalign/*local {folder}/deepalign/*fasta')
