import os
import subprocess
import numpy as np
import csv

# enter the folder first
# return q_pdb, t_pdb, entry
# q_pdb: query protein
# t_pdb: target protein
# entry: name for output file

os.chdir('Malidup')
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
    os.system(f'DeepAlign {q_pdb} {t_pdb} -o {entry}.deepalign.ali -P 0')
    os.system(f"sed '/^>/d' {entry}.deepalign.ali.fasta > {entry}.deepalign.ali")
    
    gt = entry+'.manual.ali'
    da = entry+'.deepalign.ali'
    os.system(f'python ../../accuracy.py {gt} {da} --verbose 0 --folder {folder} --query {q_pdb} --target {t_pdb} --model deepalign --outfile {entry}.accuracy')
    with open(f'{entry}.deepalign.ali.score','r') as f:
        f.readline()
        scores = f.readline().strip()
    scores=scores.split('->')[2].split()
    lalign, rmsd, tmscore = scores

    if not os.path.exists(f'{entry}.tmscore'):
        with open(f'{entry}.tmscore','w') as f:
            csv.writer(f).writerow(['Folder', 'Query', 'Target', 'Model', 'RMSD', 'L_align', 'TM-score'])
    else:
        with open(f'{entry}.tmscore','a') as f:
            csv.writer(f).writerow([folder, q_pdb, t_pdb, 'deepalign', rmsd, lalign, tmscore])

    os.chdir('../')
