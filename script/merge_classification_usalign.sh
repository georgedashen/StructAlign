cat tmscores_b64_*.csv > tmscore_full.csv
sed -i '/Query/d' tmscore_full.csv
sed -i '1 s/^/Query,Target,TM-score,L_align,SO_Identity,AA_Identity,AA_Identity_full,RMSD\n/' tmscore_full.csv
wc -l tmscore_full.csv

