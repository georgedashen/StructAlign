python embeddings.py start ../SCOP140.fasta ../SCOP140/SCOP140_plmblast_database -embedder pt --gpu -bs 0 --asdir

python script/dbtofile.py ../SCOP140/SCOP140_plmblast_database

python embeddings.py start ../pdb70.fasta ../SCOP140/pdb70_plmblast_database -embedder pt --gpu -bs 0 --asdir

# use prefilter to get only similarity score, fast
python script/plmblast.py ../SCOP140/pdb70_plmblast_database ../SCOP140/SCOP140_plmblast_database ../SCOP140/plmblast_prefilter.json -oc
python ../script/processJson.py ../SCOP140/plmblast_prefilter.json ../SCOP140/SCOP140_plmblast_database.csv ../SCOP140/pdb70_plmblast_database.csv ../SCOP140/plmblast_prefilter classification

# or get detailed alignment results, slow
python script/plmblast.py ../SCOP140/pdb70_plmblast_database ../SCOP140/SCOP140_plmblast_database ../SCOP140/plmblast_results.csv
python ../script/processPLM.py ../SCOP140/plmblast_results.csv ../SCOP140.fasta ../pdb70.fasta ../SCOP140/ordered_pooled/ classification
