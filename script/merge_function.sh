cat tmscore_b*.csv > tmscore_full.csv
sed -i '/Query/d' tmscore_full.csv
sed -i '1 s/^/Query,Target,TM-score,RMSD,SO_Identity,AA_Identity,AA_Identity_full\n/' tmscore_full.csv
#wc -l tmscore_full_less.csv
