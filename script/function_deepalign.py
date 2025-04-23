import os, sys
import pandas as pd
import glob
import subprocess
import csv

batch = int(sys.argv[1]) # number of batch
idx = int(sys.argv[2]) # index of batch, start from 0

check_idx = False
start_idx = 1

os.chdir('CAFA3_MF')
os.system(f'mkdir -p DeepAlign')
overwrite = True
if not os.path.exists(f'DeepAlign/tmscore_b{batch}_{idx}.csv') or overwrite:
    with open(f'DeepAlign/tmscore_b{batch}_{idx}_less.csv','w') as f:
        csv.writer(f).writerow(['Query', 'Target', 'TM-score', 'RMSD', 'SO_Identity', 'AA_Identity', 'AA_Identity_full'])
elif check_idx:
    start_idx = subprocess.run(f"wc -l DeepAlign/tmscore_b{batch}_{idx}.csv", stdout=subprocess.PIPE, shell=True, text=True).stdout.strip().split(' ')[0]


df=pd.read_csv('comparison.list')
if len(df) % batch == 0:
    bs = len(df) // batch
else:
    bs = len(df)//batch + 1

df = df.iloc[bs*idx:bs*(idx+1),:]
df = df.reset_index()
df = df.iloc[int(start_idx)-1:,:]
df = df.reset_index()


import numpy as np
#total_time=0

for i in range(len(df)):
    q,t = df['Query'][i], df['Target'][i]
    q_pdb, t_pdb = f'test/{q}.pdb', f'train/{t}.pdb'
    os.system(f'mkdir -p DeepAlign/{q}')
    if not os.path.exists(f'DeepAlign/{q}/{q}_{t}.deepalign.ali.score'):
        if check_mode:
            with open(f'DeepAlign/check_b{batch}_{idx}.csv','a') as f:
                csv.writer(f).writerow([q, t])
            continue
        if not extract_mode:
            try:
                os.system(f'DeepAlign {q_pdb} {t_pdb} -o DeepAlign/{q}/{q}_{t}.deepalign.ali -P 0')
            except:
                print(f'{q}_{t} Error!')
                continue
    if os.path.exists(f'DeepAlign/{q}/{q}_{t}.deepalign.ali.score') and not check_mode:
        with open(f'DeepAlign/{q}/{q}_{t}.deepalign.ali.score','r') as f:
            f.readline()
            scores = f.readline().strip()
        if not scores:
            os.system(f'DeepAlign {q_pdb} {t_pdb} -o DeepAlign/{q}/{q}_{t}.deepalign.ali -P 0')
            with open(f'DeepAlign/{q}/{q}_{t}.deepalign.ali.score','r') as f:
                f.readline()
                scores = f.readline().strip()
        score=scores.split('->')[2].split()
        _, _, qlen, tlen = scores.split('->')[0].split()
        lalign, rmsd, tmscore = score
        identity = int(lalign) / int(qlen)
        aa = scores.split('->')[4].split()[0]
        if int(lalign):
            aa_identity = int(aa) / int(lalign)
        else:
            aa_identity = 0
        aa_identity_full = int(aa) / int(qlen)
        with open(f'DeepAlign/tmscore_b{batch}_{idx}.csv','a') as f:
            csv.writer(f).writerow([q, t, tmscore, rmsd, identity, aa_identity, aa_identity_full])
        if os.path.exists(f'DeepAlign/{q}/{q}_{t}.deepalign.ali.pdb'):
            os.system(f'rm DeepAlign/{q}/{q}_{t}.deepalign.ali.pdb')
