import os, sys
import pandas as pd
import subprocess
import itertools
import csv

batch = int(sys.argv[1]) #64
idx = int(sys.argv[2]) #0-63

overwrite = False
os.system("mkdir -p USalign")
if not os.path.exists(f'USalign/tmscores_b{batch}_{idx}.csv') or overwrite:
    with open(f'USalign/tmscores_b{batch}_{idx}.csv','w') as f:
        csv.writer(f).writerow(['Query', 'Target', 'TM-score','L_align', 'SO_Identity', 'AA_Identity', 'AA_Identity_full', 'RMSD'])

df = pd.read_csv('comparison.list')
if len(df) % batch == 0:
    batch_size = len(df) // batch
else:
    batch_size = len(df) // batch + 1

df = df.iloc[batch_size*idx:batch_size*(idx+1),:]

overwrite = False
for i in range(len(df)):
    q = df.iloc[i,0]
    t = df.iloc[i,1]
    q_pdb = f'structure_data/scope140_pdb/{q}.ent'
    t_pdb = f'structure_data/pdb70_and_scope/{t}.pdb'
    
    os.system(f'mkdir -p USalign/{q}')
    if not os.path.exists(f'USalign/{q}/{q}_{t}.usalign') or overwrite:
        try:
            os.system(f'USalign {q_pdb} {t_pdb} -mm 5 -outfmt 2 > USalign/{q}/{q}_{t}.usalign')
        except:
            print('{q}_{t} error!')
    if os.path.exists(f'USalign/{q}/{q}_{t}.usalign'):
        stats = subprocess.run(f"sed -n '2p' USalign/{q}/{q}_{t}.usalign", stdout=subprocess.PIPE, shell=True, text=True).stdout.strip().split()
        try:
            queryLen, targetLen, TM1, TM2, rmsd, lalign, aa_identity = stats[8], stats[9], stats[2], stats[3], stats[4], stats[10], stats[7]
        except:
            print(f'{q}_{t} error')
            continue
        aa_identity_full = round(float(aa_identity) * int(lalign)) / sum([int(queryLen),int(targetLen)]) * 2
        identity = int(lalign) / sum([int(queryLen),int(targetLen)]) * 2
        tmscore = max(TM1,TM2)
        with open(f'USalign/tmscores_b{batch}_{idx}.csv','a') as f:
            csv.writer(f).writerow([q, t, tmscore, lalign, identity, aa_identity, aa_identity_full, rmsd])

