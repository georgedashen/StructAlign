import pandas as pd
from Bio import SeqIO

infile = sys.argv[1]
qfa = sys.argv[2]
tfa = sys.argv[3]
outdir = sys.argv[4]
task = sys.argv[5]

target_dict={}
query_dict={}

for idx,record in enumerate(SeqIO.parse(tfa, "fasta")):
    target_dict[idx] = record.id

for idx,record in enumerate(SeqIO.parse(qfa, "fasta")):
    query_dict[idx] = record.id

df=pd.read_csv(infile,sep=';')
df['qid'] = [query_dict[i] for i in df['qid']]
df['sid'] = [target_dict[i] for i in df['sid']]
df = df[['qid','sid','score','ident','similarity']]

if task == 'classification':
    df[['sid','qid','score']].to_csv(f'{outdir}/plmblast_TMscore',header=None,index=False,sep='\t')
    df[['sid','qid','ident']].to_csv(f'{outdir}/plmblast_Identity',header=None,index=False,sep='\t')
    df[['sid','qid','similarity']].to_csv(f'{outdir}/plmblast_Similarity',header=None,index=False,sep='\t')
else:
    df.columns = ['Query','Target','Score','Identity','Similarity']
    df[['Query','Target','Score']].to_csv(f'{outdir}/plmblast_TMscore',index=False)
    df[['Query','Target','Identity']].to_csv(f'{outdir}/plmblast_Identity',index=False)
    df[['Query','Target','Similarity']].to_csv(f'{outdir}/plmblast_Similarity',index=False)
