# MalisamPDB(MalidupPDB) is a folder that contains all structures in the Malisam(Malidup) dataset
  
mkdir -p MalisamPDB
mkdir -p MalidupPDB

find Malisam -type f -name '*_*.pdb' ! -name '*.*.*' -exec cp {} MalisamPDB/ \;
find Malidup -type f -name '*_*.pdb' ! -name '*.*.*' -exec cp {} MalidupPDB/ \;

foldseek createdb MalisamPDB/ database/Malisam
foldseek createdb MalidupPDB/ database/Malidup

mkdir -p Malidup_result
foldseek easy-search --exhaustive-search 1 -e 1000 database/Malidup database/Malidup Malidup_result/aln.tsv tmpFolder --format-output "query,target,qaln,taln"
cp Malidup_result/aln.tsv Malidup_foldseek.tsv

mkdir -p Malisam_result
foldseek easy-search --exhaustive-search 1 -e 1000 database/Malisam database/Malisam Malisam_result/aln.tsv tmpFolder --format-output "query,target,qaln,taln"
cp Malisam_result/aln.tsv Malisam_foldseek.tsv
