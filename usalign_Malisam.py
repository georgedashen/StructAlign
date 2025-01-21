import os
import subprocess
import numpy as np
import csv

record_accuracy = 1
record_tmscore = 1
overwrite = 1
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
    q = folder[:7]
    t = folder[7:]
    qp = q[1:5]
    tp = t[1:5]
    FNL = os.listdir('./')
    FNL = np.array(FNL)

    # Retrieve pdb files
    q_pdb = max(FNL[[qp in i and i.endswith('.pdb') and i.count('.')==1 for i in FNL]], key=len)
    t_pdb = max(FNL[[tp in i and i.endswith('.pdb') and i.count('.')==1 for i in FNL]], key=len)
    
    entry = folder
    # algo q_pdb t_pdb -o entry.ali
    # further processing is required according to the output format
    if not os.path.exists(f'{entry}.usalign') or overwrite:
        os.system(f'USalign {q_pdb} {t_pdb} -mm 5 > {entry}.usalign')
    
    gt = entry+'.manual.ali'
    ua = entry+'.usalign'
    if record_accuracy:
        os.system(f"sed -i '/USalign/d' {entry}.accuracy")
        os.system(f'python ../../accuracy.py {gt} {ua} --verbose 0 --folder {folder} --query {q_pdb} --target {t_pdb} --model USalign --outfile {entry}.accuracy  --input_type res')

    if not os.path.exists(f'{entry}.tmscore') and record_tmscore:
        with open(f'{entry}.tmscore','w') as f:
            csv.writer(f).writerow(['Folder', 'Query', 'Target', 'Model', 'RMSD', 'L_align', 'TM-score'])
    
    if record_tmscore:
        tmscore = max(subprocess.run(f"grep 'TM-score=' {entry}.usalign | cut -d' ' -f2", stdout=subprocess.PIPE, shell=True, text=True).stdout.strip().split())
        stats = subprocess.run(f"grep 'Aligned length=' {entry}.usalign",stdout=subprocess.PIPE, shell=True, text=True).stdout.strip().split()
        stats = [i.replace(',','') for i in stats]
        lalign, rmsd, aa_identity = stats[2], stats[4], stats[6]
        os.system(f"sed -i '/USalign/d' {entry}.tmscore")
        with open(f'{entry}.tmscore','a') as f:
            csv.writer(f).writerow([folder, q_pdb, t_pdb, 'USalign', rmsd, lalign, tmscore])

    os.chdir('../')
