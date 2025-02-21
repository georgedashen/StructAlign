import os, sys
import glob
import subprocess
import itertools
import csv

# for each algorithm to do all vs. all comparison for a set of swiss tree proteins
cluster = '' # _cluster
# enumerate and then each comparison
fold = sys.argv[1] # ST001 - ST011
folder = f'SwissTree{cluster}/' + fold
struct_list = glob.glob(folder + '/structs/*')
pair_list = list(itertools.combinations(struct_list,2))

os.system(f'mkdir -p {folder}/dali')
if not os.path.exists(f'{folder}/dali/tmscore.csv'):
    with open(f'{folder}/dali/tmscore.csv','w') as f:
        csv.writer(f).writerow(['Tree', 'Query', 'Target', 'Model', 'RMSD', 'L_align', 'Identity', 'TM-score'])
for pair in pair_list:
    q_pdb,t_pdb = pair
    q = os.path.basename(q_pdb).replace('.pdb','')
    t = os.path.basename(t_pdb).replace('.pdb','')
    os.system(f"dali.pl --pdbfile1 {q_pdb} --pdbfile2 {t_pdb} --dat1 . --dat2 . --outfmt 'summary,alignments,transrot' --clean >out 2> error")
    os.system(f"python ../extractDALI.py mol1A.txt {entry}.dali.ali")
    os.system(f'TMalign {q_pdb} {t_pdb} > {folder}/TMalign/{q}_{t}.tmalign')
    result = subprocess.run(f"grep 'Aligned length' {folder}/TMalign/{q}_{t}.tmalign", stdout=subprocess.PIPE, shell=True, text=True)
    stats = result.stdout.strip().split(' ')
    stats = [i.replace(',','') for i in stats]
    lalign, rmsd, identity = stats[2], stats[6], stats[8] #identity is calculated from matched nongapped positions
    result = subprocess.run(f"grep 'TM-score=' {folder}/TMalign/{q}_{t}.tmalign | cut -d' ' -f2", stdout=subprocess.PIPE, shell=True, text=True)
    tmscore = max(result.stdout.strip().split('\n'))
    with open(f'{folder}/dali/tmscore.csv','a') as f:
        csv.writer(f).writerow([fold, q, t, 'DALI', rmsd, lalign, identity, tmscore])
