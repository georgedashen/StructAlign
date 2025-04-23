import glob
import os, sys

#usage: python genFasta.py SwissTree/ST001/structs

folder = sys.argv[1] # a folder contain pdb files
outfile = sys.argv[2] #seq.fasta for SwissTree

pdb_list = glob.glob(folder+'/*pdb')

if os.path.exists(f'{os.path.basename(folder)}/{outfile}'):
    os.system(f"rm {os.path.basename(folder)}/{outfile}")

for i in pdb_list:
    prot = os.path.basename(i).replace(".pdb","")
    os.system(f"./pdb2fasta {i} >> {os.path.basename(folder)}/{outfile}")
