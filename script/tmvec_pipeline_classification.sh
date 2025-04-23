## code has been modified in tm_vec_utils.py and tmvec-build-database to enable large batch encoding, but currently not available
tmvec-build-database --input-fasta pdb70_1000.fasta --tm-vec-model models/tm_vec_cath_model.ckpt --tm-vec-config-path models/tm_vec_cath_model_params.json --device 'gpu' --output SCOP140/pdb70_tmvec_database

tmvec-search --query SCOP140_1000.fasta --tm-vec-model models/tm_vec_cath_model.ckpt --tm-vec-config models/tm_vec_cath_model_params.json --database SCOP140/pdb70_tmvec_database/db --metadata SCOP140/pdb70_tmvec_database/meta.npy --database-fasta pdb70_1000.fasta --device 'gpu' --output-format tabular --output SCOP140/tmvec_results/SCOP140_tabular.txt --output-embeddings SCOP140/tmvec_results/query_embeddings.npy

python script/classification_tmvec.py
