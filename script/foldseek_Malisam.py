import os
import pandas as pd
import subprocess
import numpy as np
import csv

record_accuracy = 1
record_tmscore = 1

# enter the folder first Malisam
# return q_pdb, t_pdb, entry
# q_pdb: query protein
# t_pdb: target protein
# entry: name for output file

overwrite = True

os.chdir('Malisam')
df = pd.read_table('../Malisam_foldseek.tsv',header=None)
df.columns = ['Query','Target', 'qaln', 'taln']

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
    try:
        _,_,aln1,aln2=df[(df['Query']==q_pdb.replace('.pdb',''))&(df['Target']==t_pdb.replace('.pdb',''))].values[0]
    except:
        print(f'{entry} result not found!')
        os.chdir('../')
        continue
    
    x = subprocess.run(f"../../pdb2fasta {q_pdb}", stdout=subprocess.PIPE, shell=True, text=True).stdout.strip().split('\n')[1]
    y = subprocess.run(f"../../pdb2fasta {t_pdb}", stdout=subprocess.PIPE, shell=True, text=True).stdout.strip().split('\n')[1]
    qstart = np.char.find(x,aln1.replace('-',''))
    qend = qstart + len(aln1.replace('-',''))
    tstart = np.char.find(y,aln2.replace('-',''))
    tend = tstart + len(aln2.replace('-',''))
    if qstart >= tstart:
        aln1 = x[:qstart].lower() + aln1
        aln2 = (qstart-tstart)*'-' + y[:tstart].lower() + aln2
    else:
        aln1 = (tstart-qstart)*'-' + x[:qstart].lower() + aln1
        aln2 = y[:tstart].lower() + aln2
    if len(x)-qend >= len(y)-tend:
        aln1 = aln1 + x[qend:].lower()
        aln2 = aln2 + y[tend:].lower() + (len(x)-qend-(len(y)-tend))*'-'
    else:
        aln1 = aln1 + x[qend:].lower() + (len(y)-tend-(len(x)-qend))*'-'
        aln2 = aln2 + y[tend:].lower()

    if not os.path.exists(f'{entry}.foldseek.ali'):
        subprocess.run(f"newlines=('>{q}' '{aln1}' '>{t}' '{aln2}') && printf '%s\n' ${{newlines[@]}} > {entry}.foldseek.ali.fasta",shell=True,executable='/bin/bash')
        os.system(f"sed '/^>/d' {entry}.foldseek.ali.fasta > {entry}.foldseek.ali")
    if record_accuracy:
        gt = entry+'.manual.ali'
        fs = entry+'.foldseek.ali'
        os.system(f"sed -i '/foldseek/d' {entry}.accuracy")
        os.system(f'python ../../accuracy.py {gt} {fs} --verbose 0 --folder {folder} --query {q_pdb} --target {t_pdb} --model foldseek --outfile {entry}.accuracy')

    # run TM-align
    try:
        if not os.path.exists(f'{entry}.foldseek.tmalign'):
            os.system(f'TMalign {q_pdb} {t_pdb} -I {entry}.foldseek.ali.fasta > {entry}.foldseek.tmalign')
        result = subprocess.run(f"grep 'User-specified initial alignment' {entry}.foldseek.tmalign | cut -d'=' -f2", stdout=subprocess.PIPE, shell=True, text=True)
        stats = result.stdout.strip().split(',')
        stats = [i.replace(' ','') for i in stats] # TM-score, L_align, RMSD
        tmscore, lalign, rmsd = stats
    except:
        print(f'{entry} TMalign Error!')
        os.chdir('../')
        continue

    result = subprocess.run(f"grep 'TM-score=' {entry}.foldseek.tmalign | cut -d' ' -f2", stdout=subprocess.PIPE, shell=True, text=True)
    tmscore = max(result.stdout.strip().split('\n'))

    if not os.path.exists(f'{entry}.tmscore') and record_tmscore:
        with open(f'{entry}.tmscore','w') as f:
            csv.writer(f).writerow(['Folder', 'Query', 'Target', 'Model', 'RMSD', 'L_align', 'TM-score'])
    if record_tmscore:
        os.system(f"sed -i '/foldseek/d' {entry}.tmscore")
        with open(f'{entry}.tmscore','a') as f:
            csv.writer(f).writerow([folder, q_pdb, t_pdb, 'foldseek', rmsd, lalign, tmscore])

    os.chdir('../')
