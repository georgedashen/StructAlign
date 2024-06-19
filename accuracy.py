# on work:
# the tm, dali, manual alignments in *ali files only include the ground truth region
# while the provided ground truth *aln file has start and end positions

import sys
import numpy as np

gtfile=sys.argv[1] #file ended with manual.ali
infile=sys.argv[2]

b=np.loadtxt(gtfile,dtype=str,delimiter='\n')
gt={}
qc=1 #results not start with '-'
tc=0

# make sure to set the smaller protein as the query protein
swap = False
q_len = len(b[0].replace('-',''))
t_len = len(b[1].replace('-',''))
if q_len > t_len:
    swap = True

if swap:
    b[0],b[1] = b[1],b[0]
for q,t in zip(b[0],b[1]):
    if t!='-':
        tc+=1
    if q!='-':
        qc+=1
        if q.isupper() and t.isupper():
            gt[str(qc)+q]=str(tc)+t

#a=np.loadtxt('d1a05a_d1dgsa3.tm.ali',dtype=str,delimiter='\n')
a=np.loadtxt(infile,dtype=str,delimiter='\n')
pred={}
qc=1 #results could start with '-'
tc=0

if swap:
    a[0],a[1] = a[1],a[0]
for q,t in zip(a[0],a[1]):
    if t!='-':
        tc+=1
    if q!='-':
        qc+=1
        if q.isupper() and t.isupper():
            pred[str(qc)+q]=str(tc)+t
#print(pred)

#accuracy
len_gt=len(gt)
len_pred=len(pred)
TP=0

if len_pred>0:
    common = set(gt.keys()) & set(pred.keys())
    for i in common:
        if pred[i] == gt[i]:
            TP+=1
    precision = TP / len_pred
else:
    precision = 0

recall = accuracy = TP / len_gt

print('\n')
print(f'accuracy:{accuracy}, recall:{recall}, precision:{precision}')

#output header
print('\n')
print('dataset','gt','pred','precision','recall','accuracy')
print(f'Malisam, {gtfile}, {infile}, {precision}, {recall}, {accuracy}')
print('\n')
