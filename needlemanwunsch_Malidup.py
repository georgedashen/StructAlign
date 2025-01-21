import os
import subprocess
import numpy as np
import csv

from needleman_wunsch import nw

# enter the folder first
# return q_pdb, t_pdb, entry
# q_pdb: query protein
# t_pdb: target protein
# entry: name for output file

folderlist=os.listdir('.')
for folder in folderlist:
    if os.path.isfile(folder):
        continue
    os.chdir(folder)
    FNL = os.listdir('./')
    FNL = np.array(FNL)

    # Retrieve pdb files
    # make sure the protein will smaller chain index is the query
    pdblist = FNL[[i.endswith('.pdb') and i.count('.')==1 for i in FNL]]
    pdblist.sort()
    q_pdb, t_pdb = pdblist
    entry = FNL[[i.endswith('.pdb') and i.count('.')==2 for i in FNL]][0].split('.')[0]

    # algo q_pdb t_pdb -o entry.ali
    # further processing is required according to the output format
    result = subprocess.run(f"../../pdb2fasta {q_pdb}", stdout=subprocess.PIPE, shell=True, text=True)
    x = result.stdout.strip().split('\n')[1]
    result = subprocess.run(f"../../pdb2fasta {t_pdb}", stdout=subprocess.PIPE, shell=True, text=True)
    y = result.stdout.strip().split('\n')[1]
    
    aln = nw(x,y).split('\n')

    a1 = '>'+q_pdb.replace('.pdb','')
    a2 = '>'+t_pdb.replace('.pdb','')
    subprocess.run(f"newlines=('{a1}' '{aln[0]}' '{a2}' '{aln[1]}') && printf '%s\n' ${{newlines[@]}} > {entry}.needlemanwunsch.ali.fasta",shell=True,executable='/bin/bash')
    os.system(f"sed '/^>/d' {entry}.needlemanwunsch.ali.fasta > {entry}.needlemanwunsch.ali")

    gt = entry+'.manual.ali'
    nwa = entry+'.needlemanwunsch.ali'
    os.system(f'python ../../accuracy.py {gt} {nwa} --verbose 0 --folder {folder} --query {q_pdb} --target {t_pdb} --model needlemanwunsch --outfile {entry}.accuracy')
    
    # run TM-align
    os.system(f'TMalign {q_pdb} {t_pdb} -I {entry}.needlemanwunsch.ali.fasta > {entry}.needlemanwunsch.tmalign')
    result = subprocess.run(f"grep 'User-specified initial alignment' {entry}.needlemanwunsch.tmalign | cut -d'=' -f2", stdout=subprocess.PIPE, shell=True, text=True)
    stats = result.stdout.strip().split(',')
    stats = [i.replace(' ','') for i in stats] # TM-score, L_align, RMSD
    tmscore, lalign, rmsd = stats

    result = subprocess.run(f"cat {entry}.needlemanwunsch.ali", stdout=subprocess.PIPE, shell=True, text=True)
    result = result.stdout.strip().split('\n')
    lalign = sum(1 for i in range(len(result[0])) if result[0][i].isupper() and result[1][i].isupper())

    result = subprocess.run(f"grep 'TM-score=' {entry}.needlemanwunsch.tmalign | cut -d' ' -f2", stdout=subprocess.PIPE, shell=True, text=True)
    tmscore = max(result.stdout.strip().split('\n'))

    if not os.path.exists(f'{entry}.tmscore'):
        with open(f'{entry}.tmscore','w') as f:
            csv.writer(f).writerow(['Folder', 'Query', 'Target', 'Model', 'RMSD', 'L_align', 'TM-score'])
    else:
        with open(f'{entry}.tmscore','a') as f:
            csv.writer(f).writerow([folder, q_pdb, t_pdb, 'needlemanwunsch', rmsd, lalign, tmscore])

    os.chdir('../')
