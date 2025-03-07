## code has modified in tm_vec_utils.py and tmvec-build-database to enable large batch encoding
tmvec-build-database --input-fasta train_1000.fasta --tm-vec-model models/tm_vec_cath_model.ckpt --tm-vec-config-path models/tm_vec_cath_model_params.json --device 'gpu' --output CAFA3_MF/train_tmvec_database

mkdir -p CAFA3_MF/tmvec_results

tmvec-search --query test-1000.fasta --tm-vec-model models/tm_vec_cath_model.ckpt --tm-vec-config models/tm_vec_cath_model_params.json --database CAFA3_MF/train_tmvec_database/db.npy --metadata CAFA3_MF/train_tmvec_database/meta.npy --database-fasta train_1000.fasta --device 'gpu' --output-format tabular --output CAFA3_MF/tmvec_results/tmvec_test_1000.txt --output-embeddings CAFA3_MF/tmvec_results/tmvec_test_1000_embeddings.npy

python script/function_tmvec.py
