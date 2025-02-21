python embeddings.py start ../CAFA3_MF/test.fasta ../CAFA3_MF/test_plmblast_database -embedder pt --gpu -bs 0 --asdir

python script/dbtofile.py

python scripts/plmblast.py ../CAFA3_MF/train_plmblast_database ../CAFA3_MF/test_plmblast_database ../CAFA3_MF/test_plmblast_database/all.oc.json -oc

python processJson.py
