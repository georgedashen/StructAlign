import torch
import alntools as aln
import os
from Bio import SeqIO
import subprocess
import numpy as np
import csv

extr=aln.Extractor(bfactor='global')

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
    
    q = q_pdb.replace('.pdb','')
    t = t_pdb.replace('.pdb','')

    seq1 = str(list(SeqIO.parse(f'../../MalidupFasta/{q}.fasta', format='fasta'))[0].seq)
    seq2 = str(list(SeqIO.parse(f'../../MalidupFasta/{t}.fasta', format='fasta'))[0].seq)
    
    emb1 = torch.load(f'../../Malidup_plmblast_pt/{q}.pt').numpy()
    emb2 = torch.load(f'../../Malidup_plmblast_pt/{t}.pt').numpy()

    results = extr.embedding_to_span(emb2,emb1)
    row=results.iloc[results.score.argmax()]
    ali=aln.draw_alignment(row.indices, seq1, seq2, output='str').split('\n')

    a1 = '>'+q
    a2 = '>'+t
    subprocess.run(f"newlines=('{a1}' '{ali[0]}' '{a2}' '{ali[2]}') && printf '%s\n' ${{newlines[@]}} > {entry}.plmblast.ali.fasta",shell=True,executable='/bin/bash')
    os.system(f"sed '/^>/d' {entry}.plmblast.ali.fasta > {entry}.plmblast.ali")

    gt = entry+'.manual.ali'
    pb = entry+'.plmblast.ali'
    os.system(f'python ../../accuracy.py {gt} {pb} --verbose 0 --folder {folder} --query {q_pdb} --target {t_pdb} --model plmblast --outfile {entry}.accuracy')

    # run TM-align
    os.system(f'TMalign {q_pdb} {t_pdb} -I {entry}.plmblast.ali.fasta > {entry}.plmblast.tmalign')
    result = subprocess.run(f"grep 'User-specified initial alignment' {entry}.plmblast.tmalign | cut -d'=' -f2", stdout=subprocess.PIPE, shell=True, text=True)
    stats = result.stdout.strip().split(',')
    stats = [i.replace(' ','') for i in stats] # TM-score, L_align, RMSD
    tmscore, lalign, rmsd = stats

    result = subprocess.run(f"cat {entry}.plmblast.ali", stdout=subprocess.PIPE, shell=True, text=True)
    result = result.stdout.strip().split('\n')
    lalign = sum(1 for i in range(len(result[0])) if result[0][i].isupper() and result[1][i].isupper())

    result = subprocess.run(f"grep 'TM-score=' {entry}.plmblast.tmalign | cut -d' ' -f2", stdout=subprocess.PIPE, shell=True, text=True)
    tmscore = max(result.stdout.strip().split('\n'))

    if not os.path.exists(f'{entry}.tmscore'):
        with open(f'{entry}.tmscore','w') as f:
            csv.writer(f).writerow(['Folder', 'Query', 'Target', 'Model', 'RMSD', 'L_align', 'TM-score'])
    else:
        with open(f'{entry}.tmscore','a') as f:
            csv.writer(f).writerow([folder, q_pdb, t_pdb, 'plmblast', rmsd, lalign, tmscore])

    os.chdir('../')

