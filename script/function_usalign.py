import os, sys
import pandas as pd
import glob
import subprocess
import csv

batch = int(sys.argv[1]) # number of batch
idx = int(sys.argv[2]) # index of batch, start from 0

os.chdir('CAFA3_MF')
os.system(f'mkdir -p USalign')

querylist = glob.glob('test/*pdb')

if len(querylist) % batch == 0:
    bs = len(querylist) // batch
else:
    bs = len(querylist)//batch + 1

querylist = querylist[bs*idx:bs*(idx+1)]
for i in querylist:
    q = os.path.basename(i).replace('.pdb','')
    if not os.path.exists(f'USalign/usalign_{q}.result'):
        os.system(f"USalign {i} -dir2 train/ train.list -suffix .pdb -outfmt 2 -mm 5 -fast > USalign/usalign_{q}.result")
