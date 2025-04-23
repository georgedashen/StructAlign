import os,sys
import glob
import pandas as pd

filelist = glob.glob('*result')
for f in filelist:
    df = pd.read_csv(f,header=None,skiprows=1,sep='\t')
    df.columns=['Query','Target','TM1','TM2','RMSD','ID1','ID2','IDali','L1','L2','Lali']
    dft = df.iloc[:,:2]
    dft['Query'] = [i.replace('CAFA3_MF/test/','').replace('.pdb:A','') for i in df['Query']]
    dft['Target'] = [i.replace('.pdb:A','') for i in df['Target']]
    dft['TM-score'] = df['TM1']
    dft['RMSD'] = df['RMSD']
    dft['SO_Identity'] = df['Lali'] / df['L1']
    dft['AA_Identity'] = df['IDali']
    dft['AA_Identity_full'] = (df['IDali']*df['Lali']).astype('int')/df['L1']

    query = dft['Query'][0]
    dft.to_csv(f'tmscores_b_{query}.csv',index=False)

