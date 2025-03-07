import os
import subprocess
import numpy as np
import csv

from deepblast.utils import load_model
from deepblast.dataset.utils import states2alignment

model = load_model('../models/deepblast-v3.ckpt', '../models/prot_t5_xl_uniref50', device='cuda')

# enter the folder first
# return q_pdb, t_pdb, entry
# q_pdb: query protein
# t_pdb: target protein
# entry: name for output file
os.chdir('../Malidup')
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
    pred_ali = model.align(x,y)
    x_aligned, y_aligned = states2alignment(pred_ali, x, y)
    a1 = '>'+q_pdb.replace('.pdb','')
    a2 = '>'+t_pdb.replace('.pdb','')
    subprocess.run(f"newlines=('{a1}' '{x_aligned}' '{a2}' '{y_aligned}') && printf '%s\n' ${{newlines[@]}} > {entry}.deepblast.ali.fasta",shell=True,executable='/bin/bash')
    os.system(f"sed '/^>/d' {entry}.deepblast.ali.fasta > {entry}.deepblast.ali")

    gt = entry+'.manual.ali'
    db = entry+'.deepblast.ali'
    os.system(f'python ../../accuracy.py {gt} {db} --verbose 0 --folder {folder} --query {q_pdb} --target {t_pdb} --model deepblast --outfile {entry}.accuracy')
    
    # run TM-align
    os.system(f'TMalign {q_pdb} {t_pdb} -I {entry}.deepblast.ali.fasta > {entry}.deepblast.tmalign')
    result = subprocess.run(f"grep 'User-specified initial alignment' {entry}.deepblast.tmalign | cut -d'=' -f2", stdout=subprocess.PIPE, shell=True, text=True)
    stats = result.stdout.strip().split(',')
    stats = [i.replace(' ','') for i in stats] # TM-score, L_align, RMSD
    tmscore, lalign, rmsd = stats

    result = subprocess.run(f"cat {entry}.deepblast.ali", stdout=subprocess.PIPE, shell=True, text=True)
    result = result.stdout.strip().split('\n')
    lalign = sum(1 for i in range(len(result[0])) if result[0][i].isupper() and result[1][i].isupper())

    result = subprocess.run(f"grep 'TM-score=' {entry}.deepblast.tmalign | cut -d' ' -f2", stdout=subprocess.PIPE, shell=True, text=True)
    tmscore = max(result.stdout.strip().split('\n'))

    if not os.path.exists(f'{entry}.tmscore'):
        with open(f'{entry}.tmscore','w') as f:
            csv.writer(f).writerow(['Folder', 'Query', 'Target', 'Model', 'RMSD', 'L_align', 'TM-score'])
    else:
        with open(f'{entry}.tmscore','a') as f:
            csv.writer(f).writerow([folder, q_pdb, t_pdb, 'deepblast', rmsd, lalign, tmscore])

    os.chdir('../')
