import json
import os,sys
import pandas as pd
import numpy as np

injson = sys.argv[1]
qcsv = sys.argv[2]
tcsv = sys.argv[3]
outfile = sys.argv[4]
task = sys.argv[5] #classification or function

with open(injson,'r') as handle:
    df=json.load(handle)

queryid = []
targetid = []
score = []

for i in df.keys():
    for j in df[i].keys():
        queryid.append(i)
        targetid.append(os.path.basename(df[i][j]['file']))
        score.append(df[i][j]['score'])

dft = pd.DataFrame({'Query':queryid, 'Target':targetid, 'Score':score})
query_info = pd.read_csv(qcsv)
query_list = {}
for i in range(len(query_info)):
    query_list[query_info['index'][i]] = query_info['id'][i]
target_info = pd.read_csv(tcsv)
target_info['filename'] = [str(i)+'.emb' for i in target_info.iloc[:,0].values]
target_list = {}
for i in range(len(target_info)):
    target_list[target_info['filename'][i]] = target_info['id'][i]

dft['Query'] = [query_list[int(i)] for i in dft['Query']]
dft['Target'] = [target_list[i] for i in dft['Target']]

if task == 'classification':
    dft = dft[['Target','Query','Score']]
    dft.to_csv(outfile, index=False, header=None, sep='\t')
else:
    dft.to_csv(outfile, index=False)
