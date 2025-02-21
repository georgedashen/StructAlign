import glob
import os, sys

#usage: python createFastaFold.py MalidupPDB MalidupFasta

folder = sys.argv[1]
#folder = 'test'
output = sys.argv[2]
#output = 'testFasta'

pdb_list = glob.glob(folder+'/*pdb')
os.system(f'mkdir -p {output}')

for i in pdb_list:
    prot = os.path.basename(i).replace(".pdb","")
    os.system(f"./pdb2fasta {i} > {output}/{prot}.fasta")
