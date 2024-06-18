# Process the SCOP dataset so that the PDB files are able to be processed by foldseek

cd SCOP200
pdblist=$(ls | grep pdb)
for i in $pdblist; do
	awk {'printf ("%-6s%5d  %-4s% 3s%6d %11.3f%8.3f%8.3f\n",$1,$2,$3,$4,$5,$6,$7,$8)'} $i > ../SCOP200_truncated/$i
done
cd ../
