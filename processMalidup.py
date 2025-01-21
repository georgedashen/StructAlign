# Tm-score (the higher one) is normalized by the length of the smaller (query)

import os
import subprocess
import numpy as np
import csv

#delete existing results
os.system('find -name *.accuracy -type f -delete')
os.system('find -name *.tmscore -type f -delete')

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
    #print(q_pdb)
    #os.system(f"tail -2 {q_pdb} | head -1 | awk {{'print $6'}}")
    #print(t_pdb)
    #os.system(f"tail -2 {t_pdb} | head -1 | awk {{'print $6'}}")
    entry = FNL[[i.endswith('.pdb') and i.count('.')==2 for i in FNL]][0].split('.')[0]
    
    gt = entry+'.manual.ali'
    tm = entry+'.tm.ali'
    dali = entry+'.dali.ali'
    fast = entry+'.fast.ali'
    
    gt_fasta = entry+'.manual.ali.fasta'
    tm_fasta = entry+'.tm.ali.fasta'
    dali_fasta = entry+'.dali.ali.fasta'
    fast_fasta = entry+'.fast.ali.fasta'

    # Create fasta input for TM-align
    a1 = '>'+q_pdb.replace('.pdb','')
    a2 = '>'+t_pdb.replace('.pdb','')
    if not os.path.exists(gt_fasta):
        subprocess.run(f"mapfile -t array < {gt} && newlines=('{a1}' ${{array[0]}} '{a2}' ${{array[1]}}) && printf '%s\n' ${{newlines[@]}} > {gt_fasta}",shell=True,executable='/bin/bash')
    if not os.path.exists(tm_fasta):
        subprocess.run(f"mapfile -t array < {tm} && newlines=('{a1}' ${{array[0]}} '{a2}' ${{array[1]}}) && printf '%s\n' ${{newlines[@]}} > {tm_fasta}",shell=True,executable='/bin/bash')
    if not os.path.exists(dali_fasta):
        subprocess.run(f"mapfile -t array < {dali} && newlines=('{a1}' ${{array[0]}} '{a2}' ${{array[1]}}) && printf '%s\n' ${{newlines[@]}} > {dali_fasta}",shell=True,executable='/bin/bash')
    if not os.path.exists(fast_fasta):
        subprocess.run(f"mapfile -t array < {fast} && newlines=('{a1}' ${{array[0]}} '{a2}' ${{array[1]}}) && printf '%s\n' ${{newlines[@]}} > {fast_fasta}",shell=True,executable='/bin/bash')
    
    # TM-score calculation
    
    #os.system(f'TMalign {q_pdb} {t_pdb} -I {gt_fasta} > {entry}.manual.tmalign')
    #os.system(f'TMalign {q_pdb} {t_pdb} -I {tm_fasta} > {entry}.tm.tmalign')
    #os.system(f'TMalign {q_pdb} {t_pdb} -I {dali_fasta} > {entry}.dali.tmalign')
    #os.system(f'TMalign {q_pdb} {t_pdb} -I {fast_fasta} > {entry}.fast.tmalign')
    os.system(f'TMalign {q_pdb} {t_pdb} -I {fast_fasta} > {entry}.fast.tmalign')
    
    def getTMscore(model):
        if not os.path.exists(f'{entry}.{model}.tmalign'):
            return None
        result = subprocess.run(f"grep 'User-specified initial alignment' {entry}.{model}.tmalign | cut -d'=' -f2", stdout=subprocess.PIPE, shell=True, text=True)
        stats = result.stdout.strip().split(',')
        stats = [i.replace(' ','') for i in stats] # TM-score, L_align, RMSD
        tmscore, lalign, rmsd = stats
        with open(f'{entry}.tmscore','a') as f:
            csv.writer(f).writerow([folder, q_pdb, t_pdb, model, rmsd, lalign, tmscore])
    
    # Header
    with open(f'{entry}.tmscore','w') as f:
        csv.writer(f).writerow(['Folder', 'Query', 'Target', 'Model', 'RMSD', 'L_align', 'TM-score'])
    
    getTMscore('manual')
    getTMscore('tm')
    try:
        getTMscore('dali')
    except:
        pass
    getTMscore('fast')
    
    # Alignment accuracy
    #os.system(f'python ../../accuracy.py {gt} {tm} --verbose 0 --folder {folder} --query {q_pdb} --target {t_pdb} --model tm --outfile {entry}.accuracy')
    #os.system(f'python ../../accuracy.py {gt} {dali} --verbose 0 --folder {folder} --query {q_pdb} --target {t_pdb} --model dali --outfile {entry}.accuracy')
    #os.system(f'python ../../accuracy.py {gt} {fast} --verbose 0 --folder {folder} --query {q_pdb} --target {t_pdb} --model fast --outfile {entry}.accuracy')
    os.chdir('../')
