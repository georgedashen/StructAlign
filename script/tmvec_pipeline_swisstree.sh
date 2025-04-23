# Taking ST001 family as an example

tmvec-build-database --input-fasta SwissTree/ST001/seq.fasta --tm-vec-model models/tm_vec_cath_model.ckpt --tm-vec-config-path models/tm_vec_cath_model_params.json --device 'gpu' --output SwissTree/ST001/tmvec_database

tmvec-search --query SwissTree/ST001/seq.fasta --tm-vec-model models/tm_vec_cath_model.ckpt --tm-vec-config models/tm_vec_cath_model_params.json --database SwissTree/ST001/tmvec_database/db.npy --metadata SwissTree/ST001/tmvec_database/meta.npy --database-fasta SwissTree/ST001/seq.fasta --device 'gpu' --output-format tabular --output SwissTree/ST001/tmvec_database/tabular.txt --output-embeddings SwissTree/ST001/tmvec_database/embedding.npy
