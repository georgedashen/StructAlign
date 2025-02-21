import pandas as pd

df=pd.read_csv('tmscore_full.csv')
dft = df[['Target','Query','TM-score']]
dft.sort_values('TM-score',ascending=False,inplace=True)
dft.columns = ['target','query','tmscore']
dft.to_csv('/data2/zhuoyang/StructAlign/SCOP140/ordered_pooled_new/USalign_TMscore',header=None,index=False,sep='\t')

dft = df[['Target','Query','SO_Identity']]
dft.sort_values('SO_Identity',ascending=False,inplace=True)
dft.columns = ['target','query','SO-Identity']
dft.to_csv('/data2/zhuoyang/StructAlign/SCOP140/ordered_pooled_new/USalign_SO-Identity',header=None,index=False,sep='\t')

dft = df[['Target','Query','AA_Identity']]
dft.sort_values('AA_Identity',ascending=False,inplace=True)
dft.columns = ['target','query','AA-Identity']
dft.to_csv('/data2/zhuoyang/StructAlign/SCOP140/ordered_pooled_new/USalign_AA-Identity',header=None,index=False,sep='\t')

dft = df[['Target','Query','AA_Identity_full']]
dft.sort_values('AA_Identity_full',ascending=False,inplace=True)
dft.columns = ['target','query','AA-Identity-full']
dft.to_csv('/data2/zhuoyang/StructAlign/SCOP140/ordered_pooled_new/USalign_AA-Identity-full',header=None,index=False,sep='\t')

dft = df[['Target','Query','RMSD']]
dft['RMSD'] = -dft['RMSD']
dft.sort_values('RMSD',ascending=False,inplace=True)
dft.columns = ['target','query','RMSD']
dft.to_csv('/data2/zhuoyang/StructAlign/SCOP140/ordered_pooled_new/USalign_RMSD',header=None,index=False,sep='\t')
