import pandas as pd
import os

id_list = pd.read_csv('CAFA3_MF/id_mapping',header=None)
id_dict = {}
for i in range(len(id_list)):
    id_dict[id_list.iloc[i,0]] = id_list.iloc[i,1]

df = pd.read_table('CAFA3_MF/tmvec_results/tmvec_test_1000.txt')
df = df.loc[:,['query_id','database_id','tm-score']]
df.columns = ['Query', 'Target', 'TM-score']
df = df.sort_values('TM-score',ascending=False)
df['Query'] = [id_dict[i] for i in df['Query']]

df.to_csv('CAFA3_MF/TM-Vec_TMscore', index=False)

