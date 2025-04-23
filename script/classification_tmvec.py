# from command line: tm-vec-search ...
# tmvec-build-database --input-fasta SCOP140/pdb70.fasta --tm-vec-model /home/zhuoyang/StructAlign-evaluator/models/tm_vec_cath_model.ckpt --tm-vec-config-path /home/zhuoyang/StructAlign-evaluator/models/tm_vec_cath_model_params.json --device 'gpu' --output SCOP140/pdb70_tmvec_database --batch 32
# tmvec-search --query SCOP140/SCOP140.fasta --tm-vec-model /home/zhuoyang/StructAlign-evaluator/models/tm_vec_cath_model.ckpt --tm-vec-config /home/zhuoyang/StructAlign-evaluator/models/tm_vec_cath_model_params.json --database SCOP140/pdb70_tmvec_database/db --metadata SCOP140/pdb70_tmvec_database/meta.npy --database-fasta SCOP140/pdb70.fasta --device 'gpu' --output-format tabular --output SCOP140/tmvec_results/pdb70_tabular.txt --output-embeddings SCOP140/tmvec_results/query_embeddings.npy

import pandas as pd
import os

df = pd.read_table('SCOP140/tmvec_results/SCOP140_tabular_1000.txt')
df = df.loc[:,['database_id','query_id','tm-score']]
df = df.sort_values('tm-score',ascending=False)

df.to_csv('SCOP140/ordered_pooled/TM-Vec_TMscore',header=None,index=False,sep='\t')

# mkdir SCOP140/ordered_pooled_new
# cp SCOP140/ordered_pooled/TM-Vec_TMscore SCOP140/ordered_pooled_new/TM-Vec_TMscore/
# cd SCOP140
# bin/evaluate_ordered_lists.pl ordered_pooled_new/ combinetable.pdb70 scope_140_targets.list pooled > evaluation_results/pooled_new_pdb70
# cat evaluation_results/pooled_pdb70 evaluation_results/pooled_new_pdb70 > evaluation_results/total_pdb70
