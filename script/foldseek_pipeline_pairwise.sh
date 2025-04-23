foldseek createdb MalisamPDB/ database/Malisam
foldseek createdb MalidupPDB/ database/Malidup

mkdir Malidup_result
foldseek easy-search database/Malidup database/Malidup Malidup_result/aln.tsv tmpFolder --format-output "query,target,qaln,taln"

mkdir Malisam_result
foldseek easy-search database/Malisam database/Malisam Malisam_result/aln.tsv tmpFolder --format-output "query,target,qaln,taln"
