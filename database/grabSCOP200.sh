filelist=$(cat SCOP200.txt)
cd SCOP200
for i in $filelist; do
	wget $i
done
