import os, sys
import pandas as pd
import glob
import subprocess
import csv

batch = int(sys.argv[1]) # number of batch
idx = int(sys.argv[2]) # index of batch, start from 0
#Q = str(sys.argv[3])

overwrite = 1
start_idx = 1

os.chdir('CAFA3_MF')
os.system(f'mkdir -p TMalign')
if not os.path.exists(f'TMalign/tmscore_b{batch}_{idx}.csv') or overwrite:
    with open(f'TMalign/tmscore_b{batch}_{idx}.csv','w') as f:
        csv.writer(f).writerow(['Query', 'Target', 'TM-score', 'RMSD', 'SO-Identity', 'AA_Identity', 'AA_Identity_full'])
else:
    start_idx = subprocess.run(f"wc -l TMalign/tmscore_b{batch}_{idx}.csv", stdout=subprocess.PIPE, shell=True, text=True).stdout.strip().split(' ')[0]


df=pd.read_csv(f'comparison.list')
#df=pd.read_csv(f'comparison_{Q}.list')
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
    #t0=time.time()
    os.system(f'mkdir -p TMalign/{q}')
    if not os.path.exists(f'TMalign/{q}/{q}_{t}.tmalign'):
        try:
            os.system(f'TMalign {q_pdb} {t_pdb} > TMalign/{q}/{q}_{t}.tmalign')
    #t1=time.time()-t0
    #total_time += t1
        except:
            print(f'{q}_{t} Error!')
            continue
    result = subprocess.run(f"grep 'TM-score=' TMalign/{q}/{q}_{t}.tmalign | cut -d' ' -f2", stdout=subprocess.PIPE, shell=True, text=True).stdout.strip()
    if not result:
        os.system(f'TMalign {q_pdb} {t_pdb} > TMalign/{q}/{q}_{t}.tmalign')
        result = subprocess.run(f"grep 'TM-score=' TMalign/{q}/{q}_{t}.tmalign | cut -d' ' -f2", stdout=subprocess.PIPE, shell=True, text=True).stdout.strip()
    tmscore = result.split('\n')[0] # Normalized by query length
    stats = subprocess.run(f"grep 'Aligned length=' TMalign/{q}/{q}_{t}.tmalign", stdout=subprocess.PIPE, shell=True, text=True).stdout.strip().split()
    stats = [i.replace(',','') for i in stats]
    nali, rmsd, aa_identity = stats[2], stats[4], stats[6]
    queryLen = subprocess.run(f"grep 'Length of Chain_1' TMalign/{q}/{q}_{t}.tmalign", stdout=subprocess.PIPE, shell=True, text=True).stdout.strip().split(' ')[3]
    identity = int(nali) / int(queryLen)
    aa_identity_full = round(float(aa_identity) * int(nali)) / int(queryLen)
    with open(f'TMalign/tmscore_b{batch}_{idx}.csv','a') as f:
        csv.writer(f).writerow([q, t, tmscore, rmsd, identity, aa_identity, aa_identity_full])

