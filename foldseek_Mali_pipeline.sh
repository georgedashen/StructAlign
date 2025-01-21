# MalisamPDB(Malidup) is a folder that contains all structures in the Malisam(Malidup) dataset

mkdir database
mkdir -p database/Malisam
mkdir -p database/Malidup

find Malisam/ -type f -name "*_*.pdb" -exec cp {} database/Malisam +
find Malidup/ -type f -name "*_*.pdb" -exec cp {} database/Malidup +

foldseek createdb MalisamPDB/ database/Malisam
foldseek createdb MalidupPDB/ database/Malidup

mkdir Malidup_result
foldseek easy-search database/Malidup database/Malidup Malidup_result/aln.tsv tmpFolder --format-output "query,target,qaln,taln"
cp Malidup_result/aln.tsv Malidup_foldseek.tsv

mkdir Malisam_result
foldseek easy-search database/Malisam database/Malisam Malisam_result/aln.tsv tmpFolder --format-output "query,target,qaln,taln"
cp Malisam_result/aln.tsv Malisam_foldseek.tsv