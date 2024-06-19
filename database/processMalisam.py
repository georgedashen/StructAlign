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
    # Parse protein names
    os.chdir(folder)
    q = folder[:7]
    t = folder[7:]
    qp = q[1:5]
    tp = t[1:5]
    
    FNL = os.listdir('./')
    FNL = np.array(FNL)
    
    # Retrieve pdb files
    q_pdb = max(FNL[[qp in i and i.endswith('.pdb') and i.count('.')==1 for i in FNL]], key=len)
    #print(q_pdb)
    #os.system(f"tail -2 {q_pdb} | head -1 | awk {{'print $6'}}")
    t_pdb = max(FNL[[tp in i and i.endswith('.pdb') and i.count('.')==1 for i in FNL]], key=len)
    #print(t_pdb)
    #os.system(f"tail -2 {t_pdb} | head -1 | awk {{'print $6'}}")
    
    gt_pattern = folder+'.aln'
    gt = folder+'.manual.ali'
    tm = folder+'.tm.ali'
    dali = folder+'.dali.ali'
    fast = folder+'.fast.ali'
    
    gt_fasta = folder+'.manual.ali.fasta'
    tm_fasta = folder+'.tm.ali.fasta'
    dali_fasta = folder+'.dali.ali.fasta'
    fast_fasta = folder+'.fast.ali.fasta'
    '''
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
    os.system(f'TMalign {q_pdb} {t_pdb} -I {gt_fasta} > {folder}.manual.tmalign')
    os.system(f'TMalign {q_pdb} {t_pdb} -I {tm_fasta} > {folder}.tm.tmalign')
    os.system(f'TMalign {q_pdb} {t_pdb} -I {dali_fasta} > {folder}.dali.tmalign')
    os.system(f'TMalign {q_pdb} {t_pdb} -I {fast_fasta} > {folder}.fast.tmalign')
    '''
    def getTMscore(model):
    # The L_aligned and RSMD obtained are different from original results
        result = subprocess.run(f"grep 'User-specified initial alignment' {folder}.{model}.tmalign | cut -d'=' -f2", stdout=subprocess.PIPE, shell=True, text=True)
        stats = result.stdout.strip().split(',')
        stats = [i.replace(' ','') for i in stats] # TM-score, L_align, RMSD
        tmscore, lalign, rmsd = stats

        result = subprocess.run(f"cat {folder}.{model}.ali", stdout=subprocess.PIPE, shell=True, text=True)
        result = result.stdout.strip().split('\n')
        lalign = sum(1 for i in range(len(result[0])) if result[0][i].isupper() and result[1][i].isupper())
            
        with open(f'{folder}.tmscore','a') as f:
            csv.writer(f).writerow([folder, q_pdb, t_pdb, model, rmsd, lalign, tmscore])
    
    # Header
    with open(f'{folder}.tmscore','w') as f:
        csv.writer(f).writerow(['Folder', 'Query', 'Target', 'Model', 'RMSD', 'L_align', 'TM-score'])
    
    getTMscore('manual')
    getTMscore('tm')
    getTMscore('dali')
    getTMscore('fast')
    
    # Alignment accuracy
    os.system(f'python ../../accuracy.py {gt} {tm} --verbose 0 --folder {folder} --query {q_pdb} --target {t_pdb} --model tm --outfile {folder}.accuracy')
    os.system(f'python ../../accuracy.py {gt} {dali} --verbose 0 --folder {folder} --query {q_pdb} --target {t_pdb} --model dali --outfile {folder}.accuracy')
    os.system(f'python ../../accuracy.py {gt} {fast} --verbose 0 --folder {folder} --query {q_pdb} --target {t_pdb} --model fast --outfile {folder}.accuracy')
    
    os.chdir('../')
