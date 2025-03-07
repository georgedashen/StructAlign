foldseek createdb CAFA3_MF/test/ database/cafa_mf_test

foldseek createdb CAFA3_MF/train/ database/cafa_mf_train

foldseek easy-search database/cafa_mf_test database/cafa_mf_train CAFA3_MF_result/allvsall.tsv tmp --format-output 'query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,lddt,alntmscore'

python script/function_foldseek.py
