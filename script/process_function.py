import pandas as pd
import sys

#sim = sys.argv[1] #which score, can be TM-score, RMSD, SO_Identity, AA_Identity, and AA_Identity_full
method = sys.argv[1]

df=pd.read_csv('tmscore_full.csv')
#df=df[['Query','Target',sim]]
dft = df.iloc[:,[0,1,2]]
dfs = df.iloc[:,[0,1,4]]
dfr = df.iloc[:,[0,1,3]]
dfa = df.iloc[:,[0,1,5]]
dfb = df.iloc[:,[0,1,6]]

dfr['RMSD'] = -dfr['RMSD'] / max(dfr['RMSD'])

dft=df.sort_values(by='TM-score',ascending=False)
dfr=dfr.sort_values(by='RMSD',ascending=False)
dfs=dfs.sort_values(by='SO_Identity',ascending=False)
dfa=dfa.sort_values(by='AA_Identity',ascending=False)
dfb=dfb.sort_values(by='AA_Identity_full',ascending=False)

dft.to_csv(f'../{method}_TMscore',index=False)
dfr.to_csv(f'../{method}_RMSD',index=False)
dfs.to_csv(f'../{method}_SO-Identity',index=False)
dfa.to_csv(f'../{method}_AA-Identity',index=False)
dfb.to_csv(f'../{method}_AA-Identity-full',index=False)

