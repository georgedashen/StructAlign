import os, sys
import pandas as pd
import glob
import subprocess
import csv

batch = int(sys.argv[1]) # number of batch
idx = int(sys.argv[2]) # index of batch, start from 0

overwrite = 0
start_idx = 1

os.chdir('CAFA3_MF')
os.system('mkdir -p kpax_results')
if not os.path.exists(f'kpax_results/kpax_scores_b{batch}_{idx}.csv') or overwrite:
    with open(f'kpax_results/kpax_scores_b{batch}_{idx}.csv','w') as f:
        csv.writer(f).writerow(['Query', 'Target', 'TM-score', 'RMSD','SO_Identity','AA_Identity','AA_Identity_full','K-score','G-score','J-score'])
else:
    start_idx = subprocess.run(f"wc -l kpax_results/kpax_scores_b{batch}_{idx}.csv", stdout=subprocess.PIPE, shell=True, text=True).stdout.strip().split(' ')[0]

if not os.path.exists(f'kpax_results/kpax_flex_scores_b{batch}_{idx}.csv') or overwrite:
    with open(f'kpax_results/kpax_flex_scores_b{batch}_{idx}.csv','w') as f:
        csv.writer(f).writerow(['Query', 'Target', 'TM-score', 'RMSD','SO_Identity', 'AA_Identity','AA_Identity_full','K-score', 'G-score', 'J-score'])

df=pd.read_csv('comparison.list')
if len(df) % batch == 0:
    bs = len(df) // batch
else:
    bs = len(df)//batch + 1

df = df.iloc[bs*idx:bs*(idx+1),:]
df = df.reset_index()
df = df.iloc[int(start_idx)-1:,:]
df = df.reset_index()

import time
#total_time=0

for i in range(len(df)):
    q,t = df['Query'][i], df['Target'][i]
    q_pdb, t_pdb = f'test/{q}.pdb', f'train/{t}.pdb'
    os.system(f'mkdir -p kpax_results/{q}')
    #t0=time.time()
    if not os.path.exists(f'kpax_results/{q}/{t}_{q}.kalign'):
        try:
            os.system(f'kpax -nostdout -nowrite -kalign {q_pdb} {t_pdb}')
        #t1=time.time()-t0
        #total_time += t1
        except:
            print(f'{t}_{q} Error!')
    if os.path.exists(f'kpax_results/{q}/{t}_{q}.kalign'):
        tmscore = subprocess.run(f"grep 'Tscore' kpax_results/{q}/{t}_{q}.kalign | cut -d':' -f2", stdout=subprocess.PIPE, shell=True, text=True).stdout.strip()
        rmsd = subprocess.run(f"grep 'RMSD-align' kpax_results/{q}/{t}_{q}.kalign | cut -d':' -f2", stdout=subprocess.PIPE, shell=True, text=True).stdout.strip()
        nq = subprocess.run(f"grep 'Nquery' kpax_results/{q}/{t}_{q}.kalign | cut -d':' -f2", stdout=subprocess.PIPE, shell=True, text=True).stdout.strip()
        nt = subprocess.run(f"grep 'Ntarget' kpax_results/{q}/{t}_{q}.kalign | cut -d':' -f2", stdout=subprocess.PIPE, shell=True, text=True).stdout.strip()
        na = subprocess.run(f"grep 'Naligned' kpax_results/{q}/{t}_{q}.kalign | cut -d':' -f2", stdout=subprocess.PIPE, shell=True, text=True).stdout.strip()
        nident = subprocess.run(f"grep 'Nidentity' kpax_results/{q}/{t}_{q}.kalign | cut -d':' -f2 | cut -d' ' -f4", stdout=subprocess.PIPE, shell=True, text=True).stdout.strip()
        identity = int(na) / max(int(nq),int(nt)) #should use nq?
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
    
    if not os.path.exists(f'kpax_results/{q}/{t}_{q}_flex.kalign'):
        try:
            os.system(f'kpax -flex -nostdout -nowrite -kalign {q_pdb} {t_pdb}')
        #t1=time.time()-t0
        #total_time += t1
        except:
            print(f'{t}_{q} Error!')
    if os.path.exists(f'kpax_results/{q}/{t}_{q}_flex.kalign'):
        tmscore = subprocess.run(f"grep 'Tscore' kpax_results/{q}/{t}_{q}_flex.kalign | cut -d':' -f2", stdout=subprocess.PIPE, shell=True, text=True).stdout.strip()
        rmsd = subprocess.run(f"grep 'RMSD-align' kpax_results/{q}/{t}_{q}_flex.kalign | cut -d':' -f2", stdout=subprocess.PIPE, shell=True, text=True).stdout.strip()
        nq = subprocess.run(f"grep 'Nquery' kpax_results/{q}/{t}_{q}_flex.kalign | cut -d':' -f2", stdout=subprocess.PIPE, shell=True, text=True).stdout.strip()
        nt = subprocess.run(f"grep 'Ntarget' kpax_results/{q}/{t}_{q}_flex.kalign | cut -d':' -f2", stdout=subprocess.PIPE, shell=True, text=True).stdout.strip()
        na = subprocess.run(f"grep 'Naligned' kpax_results/{q}/{t}_{q}_flex.kalign | cut -d':' -f2", stdout=subprocess.PIPE, shell=True, text=True).stdout.strip()
        nident = subprocess.run(f"grep 'Nidentity' kpax_results/{q}/{t}_{q}_flex.kalign | cut -d':' -f2 | cut -d' ' -f4", stdout=subprocess.PIPE, shell=True, text=True).stdout.strip()
        identity = int(na) / max(int(nq),int(nt)) #should use nq?
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
