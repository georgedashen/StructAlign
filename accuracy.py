# on work:
# the tm, dali, manual alignments in *ali files only include the ground truth region
# while the provided ground truth *aln file has start and end positions

import sys
import numpy as np

gtfile=sys.argv[1] #file ended with manual.ali
infile=sys.argv[2]

#gt='d1a05a_d1dgsa3.aln'
b=np.loadtxt(gtfile,dtype=str,delimiter='\n')
gt={}
qc=1 #results not start with '-'
tc=0
for q,t in zip(b[0],b[1]):
    if t!='-':
        tc+=1
    if q!='-':
        gt[str(qc)+q.upper()]=str(tc)+t.upper()
        qc+=1

#a=np.loadtxt('d1a05a_d1dgsa3.tm.ali',dtype=str,delimiter='\n')
a=np.loadtxt(infile,dtype=str,delimiter='\n')
pred={}
pred_conf={}
qc=1
#results could start with '-'
tc=0

for q,t in zip(a[0],a[1]):
    if t!='-':
        tc+=1
    if q!='-':
        pred[str(qc)+q.upper()]=str(tc)+t.upper()
        if q.isupper():
            pred_conf[str(qc)+q]=str(tc)+t
        qc+=1
#print(pred)
#print(pred_conf)

#accuracy
len_gt=len(gt)
TP=0
TPC=0
for i in gt.keys():
    if gt[i] == pred[i]:
        TP+=1

if len(pred) == len(pred_conf):
    TPC = TP
    len_predconf = len_gt
else:
    #common = intersect(gt.keys(),pred_conf.keys())
    common = set(gt.keys()) & set(pred_conf.keys())
    len_predconf = len(common)
    for i in common:
        if pred_conf[i] == gt[i]:
            TPC+=1

precision = recall = accuracy = TP / len_gt
recall_conf = accuracy_conf = TPC / len_gt
precision_conf = TPC / len_predconf

print('\n')
print('Results without confidence levels:')
print(f'accuracy:{accuracy}, recall:{recall}, precision:{precision}')
if len(pred) != len(pred_conf):
    print('Results with confidence levels:')
    print(f'accuracy:{accuracy_conf}, recall:{recall_conf}, precision:{precision_conf}')

#output header
print('\n')
print('dataset','gt','pred','accuracy','c-accuracy','c-recall','c-precision')
print(f'Malisam, {gtfile}, {infile}, {accuracy}, {accuracy_conf}, {recall_conf}, {precision_conf}')
print('\n')
