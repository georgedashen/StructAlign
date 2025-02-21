# from command line: foldseek easy-search ...
# foldseek createdb SCOP140/structure_data/scope140_pdb/ database/scope140
# foldseek createdb SCOP140/structure_data/pdb70_and_scope/ database/pdb70
# foldseek easy-search database/scope140 database/pdb70 SCOP140_results/allvsall.tsv tmp --format-output 'query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,lddt,alntmscore'

import pandas as pd
import os

info=pd.read_table('SCOP140/scope_140_targets.list',header=None)
df=pd.read_table('SCOP140_results/allvsall.tsv',header=None)

df[0]=df[0].str.replace(r'_MODEL_[0-9]*_A','',regex=True)
df_tm=df[[1,0,13]]
df_fident=df[[1,0,2]]
df_fident.columns=['target','query','fident']
df_tm.columns=['target','query','tmscore']

df_tm = df_tm.sort_values('tmscore',ascending=False)
df_fident = df_fident.sort_values('fident',ascending=False)

df_tm.to_csv('SCOP140/ordered_pooled/foldseek_TMscore',header=None,index=False,sep='\t')
df_fident.to_csv('SCOP140/ordered_pooled/foldseek_Fident',header=None,index=False,sep='\t')
