import os, sys
import glob
import subprocess
import itertools
import csv

# for each algorithm to do all vs. all comparison for a set of swiss tree proteins
overwrite = 1

# enumerate and then each comparison
fold = sys.argv[1] # ST001 - ST011
folder = 'SwissTree/' + fold
os.chdir(folder)

struct_list = glob.glob('structs/*pdb')
pair_list = list(itertools.combinations(struct_list,2))

os.system(f'mkdir -p kpax')
os.system(f'mkdir -p kpax_flex')
if not os.path.exists(f'kpax/tmscore.csv') or overwrite:
    with open(f'kpax/tmscore.csv','w') as f:
        csv.writer(f).writerow(['Tree', 'Query', 'Target', 'Model', 'RMSD', 'L_align', 'SO_Identity', 'AA_Identity', 'AA_Identity_full', 'TM-score'])
    with open(f'kpax_flex/tmscore.csv','w') as f:
        csv.writer(f).writerow(['Tree', 'Query', 'Target', 'Model', 'RMSD', 'L_align', 'SO_Identity', 'AA_Identity', 'AA_Identity_full', 'TM-score'])

for pair in pair_list:
    q_pdb,t_pdb = pair
    q = os.path.basename(q_pdb).replace('.pdb','')
    t = os.path.basename(t_pdb).replace('.pdb','')
    if not os.path.exists(f'kpax/{q}_{t}.kpax.ali.fasta'):
        os.system(f'kpax {q_pdb} {t_pdb}')
        os.system(f'cp kpax_results/{q}/{t}_{q}.fasta kpax/{q}_{t}.kpax.ali.fasta')
    if not os.path.exists(f'kpax_flex/{q}_{t}.kpax_flex.ali.fasta'):
        os.system(f'kpax -flex {q_pdb} {t_pdb}')
        os.system(f'cp kpax_results/{q}/{t}_{q}_flex.fasta kpax_flex/{q}_{t}.kpax_flex.ali.fasta')
    for model in ['kpax','kpax_flex']:
        os.system(f"sed '/^>/d' {model}/{q}_{t}.{model}.ali.fasta > {model}/{q}_{t}.{model}.ali")
        with open(f'{model}/{q}_{t}.{model}.ali','r') as f:
            lines = f.read().split('\n\n')
        lines = [i.replace('\n','') for i in lines]
        subprocess.run(f"newlines=('{lines[0]}' '{lines[1]}') && printf '%s\n' ${{newlines[@]}} > {model}/{q}_{t}.{model}.ali",shell=True,executable='/bin/bash')

    for model in ['kpax', 'kpax_flex']:
        if not os.path.exists(f'{model}/{q}_{t}.{model}.tmalign'):
            os.system(f'TMalign {q_pdb} {t_pdb} -I {model}/{q}_{t}.{model}.ali.fasta > {model}/{q}_{t}.{model}.tmalign')

        stats = subprocess.run(f"grep 'Aligned length' {model}/{q}_{t}.{model}.tmalign", stdout=subprocess.PIPE, shell=True, text=True).stdout.strip().split()
        stats = [i.replace(',','') for i in stats]
        lalign, rmsd, aa_identity = stats[2], stats[4], stats[6]
        result = subprocess.run(f"cat {model}/{q}_{t}.{model}.ali", stdout=subprocess.PIPE, shell=True, text=True)
        result = result.stdout.strip().split('\n')
        lalign = sum(1 for i in range(len(result[0])) if result[0][i].isupper() and result[1][i].isupper())

        result = subprocess.run(f"grep 'TM-score=' {model}/{q}_{t}.{model}.tmalign | cut -d' ' -f2", stdout=subprocess.PIPE, shell=True, text=True)
        tmscore = max(result.stdout.strip().split('\n'))

        queryLen = subprocess.run(f"grep 'Length of Chain_1' {model}/{q}_{t}.{model}.tmalign", stdout=subprocess.PIPE, shell=True, text=True).stdout.strip().split(' ')[3]
        targetLen = subprocess.run(f"grep 'Length of Chain_2' {model}/{q}_{t}.{model}.tmalign", stdout=subprocess.PIPE, shell=True, text=True).stdout.strip().split(' ')[3]
        identity = int(lalign) / sum([int(queryLen),int(targetLen)]) * 2
        aa_identity_full = round(float(aa_identity) * int(lalign)) / sum([int(queryLen),int(targetLen)]) * 2

        with open(f'{model}/tmscore.csv','a') as f:
            csv.writer(f).writerow([fold, q, t, model, rmsd, lalign, identity, aa_identity, aa_identity_full, tmscore])
