cat kpax_scores_b64* > kpax_scores_full.csv
sed -i '/Query/d' kpax_scores_full.csv
sed -i '1 s/^/Query,Target,TM-score,RMSD,Identity,AA_Identity,AA_Identity_full,K-score,G-score,J-score\n/' kpax_scores_full.csv

cat kpax_flex_scores_b64* > kpax_flex_scores_full.csv
sed -i '/Query/d' kpax_flex_scores_full.csv
sed -i '1 s/^/Query,Target,TM-score,RMSD,Identity,AA_Identity,AA_Identity_full,K-score,G-score,J-score\n/' kpax_flex_scores_full.csv
