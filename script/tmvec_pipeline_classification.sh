## code has modified in tm_vec_utils.py and tmvec-build-database to enable large batch encoding
tmvec-build-database --input-fasta SCOP140/pdb70.fasta --tm-vec-model models/tm_vec_cath_model.ckpt --tm-vec-config-path models/tm_vec_cath_model_params.json --device 'gpu' --output SCOP140/pdb70_tmvec_database --batch 32

tmvec-search --query SCOP140/SCOP140.fasta --tm-vec-model /home/zhuoyang/StructAlign-evaluator/models/tm_vec_cath_model.ckpt --tm-vec-config models/tm_vec_cath_model_params.json --database SCOP140/pdb70_tmvec_database/db --metadata SCOP140/pdb70_tmvec_database/meta.npy --database-fasta SCOP140/pdb70.fasta --device 'gpu' --output-format tabular --output SCOP140/tmvec_results/pdb70_tabular.txt --output-embeddings SCOP140/tmvec_results/query_embeddings.npy

python script/classification_tmvec.py
