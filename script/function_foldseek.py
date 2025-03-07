import pandas as pd
import os

df=pd.read_table('CAFA3_MF_results/allvsall.tsv',header=None)

df[0]=df[0].str.replace(r'_MODEL_[0-9]*_A','',regex=True)
df_tm=df[[0,1,13]]
df_fident=df[[0,1,2]]
df_fident.columns=['Query','Target','Fident']
df_tm.columns=['Query','Target','TMscore']

df_tm = df_tm.sort_values('TMscore',ascending=False)
df_fident = df_fident.sort_values('Fident',ascending=False)

df_tm.to_csv('CAFA3_MF/foldseek_TMscore', index=False)
df_fident.to_csv('CAFA3_MF/foldseek_Fident', index=False)
