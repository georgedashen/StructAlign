foldseek createdb SCOP140/structure_data/scope140_pdb/ database/scope140

foldseek createdb SCOP140/structure_data/pdb70_and_scope/ database/pdb70

foldseek easy-search database/scope140 database/pdb70 SCOP140_results/allvsall.tsv tmp --format-output 'query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,lddt,alntmscore'

python script/classification_foldseek.py
