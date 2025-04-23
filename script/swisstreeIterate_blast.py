import os, sys
import glob
import pandas as pd
import csv

# for each algorithm to do all vs. all comparison for a set of swiss tree proteins
overwrite = 1

# enumerate and then each comparison
fold = sys.argv[1] # ST001 - ST011
folder = 'SwissTree/' + fold
os.chdir(folder)

os.system(f'mkdir -p blast')

df = pd.read_table('blast.result',header=None).iloc[:,:3]
df[0]=[i.replace(':A','') for i in df[0]]
df[1]=[i.replace(':A','') for i in df[1]]
df[2] = df[2]/100
df.columns = ['Query','Target','Identity']
df.drop_duplicates(subset=['Query','Target'],inplace=True)

if not os.path.exists('blast/tmscore.csv') or overwrite:
    df.to_csv('blast/tmscore.csv',index=False)
