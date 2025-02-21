import os, sys
import pandas as pd
import subprocess
import itertools
import csv

# Usage: python classification_kpax.py 64 0

batch = int(sys.argv[1]) #64
idx = int(sys.argv[2]) #0-63
overwrite = True

os.system("mkdir -p kpax_results")

if not os.path.exists(f'kpax_results/kpax_scores_b{batch}_{idx}.csv') or overwrite:
    with open(f'kpax_results/kpax_scores_b{batch}_{idx}.csv','w') as f:
        csv.writer(f).writerow(['Query', 'Target', 'TM-score', 'RMSD', 'Identity', 'AA_Identity', 'AA_Identity_full', 'K-score', 'G-score', 'J-score'])
if not os.path.exists(f'kpax_results/kpax_flex_scores_b{batch}_{idx}.csv') or overwrite:
    with open(f'kpax_results/kpax_flex_scores_b{batch}_{idx}.csv','w') as f:
        csv.writer(f).writerow(['Query', 'Target', 'TM-score', 'RMSD', 'Identity', 'AA_Identity', 'AA_Identity_full', 'K-score', 'G-score', 'J-score'])

overwrite = False
df = pd.read_csv('comparison.list')
if len(df) % batch == 0:
    batch_size = len(df) // batch
else:
    batch_size = len(df) // batch + 1

df = df.iloc[batch_size*idx:batch_size*(idx+1),:]

for i in range(len(df)):
    q = df.iloc[i,0]
    t = df.iloc[i,1]
    q_pdb = f'structure_data/scope140_pdb/{q}.ent'
    t_pdb = f'structure_data/pdb70_and_scope/{t}.pdb'
    if not os.path.exists(f'kpax_results/{q}/{t}_{q}.kalign') or overwrite:
        try:
            os.system(f'kpax -nostdout -nowrite -kalign {q_pdb} {t_pdb}')
        except:
            print('{t}_{q} error!')
    if os.path.exists(f'kpax_results/{q}/{t}_{q}.kalign'):
        tmscore = subprocess.run(f"grep 'Tscore' kpax_results/{q}/{t}_{q}.kalign | cut -d':' -f2", stdout=subprocess.PIPE, shell=True, text=True).stdout.strip()
        rmsd = subprocess.run(f"grep 'RMSD-align' kpax_results/{q}/{t}_{q}.kalign | cut -d':' -f2", stdout=subprocess.PIPE, shell=True, text=True).stdout.strip()
        nident = subprocess.run(f"grep 'Nidentity' kpax_results/{q}/{t}_{q}.kalign | cut -d':' -f2 | cut -d' ' -f4", stdout=subprocess.PIPE, shell=True, text=True).stdout.strip()
        nq = subprocess.run(f"grep 'Nquery' kpax_results/{q}/{t}_{q}.kalign | cut -d':' -f2", stdout=subprocess.PIPE, shell=True, text=True).stdout.strip()
        nt = subprocess.run(f"grep 'Ntarget' kpax_results/{q}/{t}_{q}.kalign | cut -d':' -f2", stdout=subprocess.PIPE, shell=True, text=True).stdout.strip()
        na = subprocess.run(f"grep 'Naligned' kpax_results/{q}/{t}_{q}.kalign | cut -d':' -f2", stdout=subprocess.PIPE, shell=True, text=True).stdout.strip()
        identity = int(na) / int(nq)
        if int(na):
            aa_identity = int(nident) / int(na)
        else:
            aa_identity = 0
        aa_identity_full = int(nident) / int(nq)
        kscore = subprocess.run(f"grep 'Kscore' kpax_results/{q}/{t}_{q}.kalign | cut -d':' -f2", stdout=subprocess.PIPE, shell=True, text=True).stdout.strip()
        gscore = subprocess.run(f"grep 'Gscore' kpax_results/{q}/{t}_{q}.kalign | cut -d':' -f2", stdout=subprocess.PIPE, shell=True, text=True).stdout.strip()
        jscore = subprocess.run(f"grep 'Jscore' kpax_results/{q}/{t}_{q}.kalign | cut -d':' -f2", stdout=subprocess.PIPE, shell=True, text=True).stdout.strip()
        with open(f'kpax_results/kpax_scores_b{batch}_{idx}.csv','a') as f:
                csv.writer(f).writerow([q, t, tmscore, rmsd, identity, aa_identity, aa_identity_full, kscore, gscore, jscore])
    if not os.path.exists(f'kpax_results/{q}/{t}_{q}_flex.kalign') or overwrite:
        try:
            os.system(f'kpax -flex -nostdout -nowrite -kalign {q_pdb} {t_pdb}')
        except:
            print('{t}_{q}_flex error!')
    if os.path.exists(f'kpax_results/{q}/{t}_{q}_flex.kalign'):
        tmscore = subprocess.run(f"grep 'Tscore' kpax_results/{q}/{t}_{q}_flex.kalign | cut -d':' -f2", stdout=subprocess.PIPE, shell=True, text=True).stdout.strip()
        rmsd = subprocess.run(f"grep 'RMSD-align' kpax_results/{q}/{t}_{q}_flex.kalign | cut -d':' -f2", stdout=subprocess.PIPE, shell=True, text=True).stdout.strip()
        nident = subprocess.run(f"grep 'Nidentity' kpax_results/{q}/{t}_{q}_flex.kalign | cut -d':' -f2 | cut -d' ' -f4", stdout=subprocess.PIPE, shell=True, text=True).stdout.strip()
        nq = subprocess.run(f"grep 'Nquery' kpax_results/{q}/{t}_{q}_flex.kalign | cut -d':' -f2", stdout=subprocess.PIPE, shell=True, text=True).stdout.strip()
        nt = subprocess.run(f"grep 'Ntarget' kpax_results/{q}/{t}_{q}_flex.kalign | cut -d':' -f2", stdout=subprocess.PIPE, shell=True, text=True).stdout.strip()
        na = subprocess.run(f"grep 'Naligned' kpax_results/{q}/{t}_{q}_flex.kalign | cut -d':' -f2", stdout=subprocess.PIPE, shell=True, text=True).stdout.strip()
        identity = int(na) / int(nq)
        if int(na):
            aa_identity = int(nident) / int(na)
        else:
            aa_identity = 0
        aa_identity_full = int(nident) / int(nq)
        kscore = subprocess.run(f"grep 'Kscore' kpax_results/{q}/{t}_{q}_flex.kalign | cut -d':' -f2", stdout=subprocess.PIPE, shell=True, text=True).stdout.strip()
        gscore = subprocess.run(f"grep 'Gscore' kpax_results/{q}/{t}_{q}_flex.kalign | cut -d':' -f2", stdout=subprocess.PIPE, shell=True, text=True).stdout.strip()
        jscore = subprocess.run(f"grep 'Jscore' kpax_results/{q}/{t}_{q}_flex.kalign | cut -d':' -f2", stdout=subprocess.PIPE, shell=True, text=True).stdout.strip()
        with open(f'kpax_results/kpax_flex_scores_b{batch}_{idx}.csv','a') as f:
            csv.writer(f).writerow([q, t, tmscore, rmsd, identity, aa_identity, aa_identity_full, kscore, gscore, jscore])

