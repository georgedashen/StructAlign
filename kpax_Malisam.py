import os
import subprocess
import numpy as np
import csv

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
    q = q_pdb.replace('.pdb','')
    t = t_pdb.replace('.pdb','')
    entry = folder

    # algo q_pdb t_pdb -o entry.ali
    # further processing is required according to the output format
    if not os.path.exists(f'kpax_results/{q}/{t}_{q}.kalign'):
        os.system(f'kpax {q_pdb} {t_pdb}')
        os.system(f'cp kpax_results/{q}/{t}_{q}.fasta {entry}.kpax.ali.fasta')
    if not os.path.exists(f'kpax_results/{q}/{t}_{q}_flex.kalign'):
        os.system(f'kpax -flex {q_pdb} {t_pdb}')
        os.system(f'cp kpax_results/{q}/{t}_{q}_flex.fasta {entry}.kpax_flex.ali.fasta')
    '''
    for model in ['kpax','kpax_flex']:
        os.system(f"sed '/^>/d' {entry}.{model}.ali.fasta > {entry}.{model}.ali")
        with open(f'{entry}.{model}.ali','r') as f:
            lines = f.read().split('\n\n')
        lines = [i.replace('\n','') for i in lines]
        subprocess.run(f"newlines=('{lines[0]}' '{lines[1]}') && printf '%s\n' ${{newlines[@]}} > {entry}.{model}.ali",shell=True,executable='/bin/bash')
    
    gt = entry+'.manual.ali'
    kp = entry+'.kpax.ali'
    kpf = entry+'.kpax_flex.ali'
    os.system(f'python ../../accuracy.py {gt} {kp} --verbose 0 --folder {folder} --query {q_pdb} --target {t_pdb} --model kpax --outfile {entry}.accuracy')
    os.system(f'python ../../accuracy.py {gt} {kpf} --verbose 0 --folder {folder} --query {q_pdb} --target {t_pdb} --model kpax_flex --outfile {entry}.accuracy')
'''
    # run TM-align
    for model in ['', '_flex']:
        '''
        os.system(f'TMalign {q_pdb} {t_pdb} -I {entry}.{model}.ali.fasta > {entry}.{model}.tmalign')
        result = subprocess.run(f"grep 'User-specified initial alignment' {entry}.{model}.tmalign | cut -d'=' -f2", stdout=subprocess.PIPE, shell=True, text=True)
        stats = result.stdout.strip().split(',')
        stats = [i.replace(' ','') for i in stats] # TM-score, L_align, RMSD
        tmscore, lalign, rmsd = stats

        result = subprocess.run(f"cat {entry}.{model}.ali", stdout=subprocess.PIPE, shell=True, text=True)
        result = result.stdout.strip().split('\n')
        lalign = sum(1 for i in range(len(result[0])) if result[0][i].isupper() and result[1][i].isupper())

        result = subprocess.run(f"grep 'TM-score=' {entry}.{model}.tmalign | cut -d' ' -f2", stdout=subprocess.PIPE, shell=True, text=True)
        tmscore = max(result.stdout.strip().split('\n'))
        '''
        tmscore = subprocess.run(f"grep 'Tscore' kpax_results/{q}/{t}_{q}{model}.kalign | cut -d':' -f2", stdout=subprocess.PIPE, shell=True, text=True).stdout.strip()
        rmsd = subprocess.run(f"grep 'RMSD-align' kpax_results/{q}/{t}_{q}{model}.kalign | cut -d':' -f2", stdout=subprocess.PIPE, shell=True, text=True).stdout.strip()
        lalign = subprocess.run(f"grep 'Naligned' kpax_results/{q}/{t}_{q}{model}.kalign | cut -d':' -f2", stdout=subprocess.PIPE, shell=True, text=True).stdout.strip()
        if not os.path.exists(f'{entry}.tmscore'):
            with open(f'{entry}.tmscore','w') as f:
                csv.writer(f).writerow(['Folder', 'Query', 'Target', 'Model', 'RMSD', 'L_align', 'TM-score'])
        with open(f'{entry}.tmscore','a') as f:
            csv.writer(f).writerow([folder, q_pdb, t_pdb, 'kpax'+model, rmsd, lalign, tmscore])

    os.chdir('../')
