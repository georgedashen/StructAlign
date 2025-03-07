python embeddings.py start ../mf-test.fasta ../CAFA3_MF/test_plmblast_database -embedder pt --gpu -bs 0 --asdir
python embeddings.py start ../mf-train.fasta ../CAFA3_MF/train_plmblast_database -embedder pt --gpu -bs 0 --asdir

python scripts/dbtofile.py ../CAFA3_MF/test_plmblast_database

python scripts/plmblast.py ../CAFA3_MF/train_plmblast_database ../CAFA3_MF/test_plmblast_database ../CAFA3_MF/plmblast_prefilter.json -oc

python ../script/processJson.py ../CAFA3_MF/plmblast_prefilter.json ../mf-test.fasta ../mf-train.fasta ../CAFA3_MF/ function
