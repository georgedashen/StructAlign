# on work:
# the tm, dali, manual alignment patterns only include the ground truth region
# while the provided ground truth has start and end positions

import sys
import numpy as np

gtfile=sys.argv[1]
infile=sys.argv[2]

#gt='d1a05a_d1dgsa3.aln'
with open(gtfile,'r') as file:
    content=file.read().strip()

patterns=content.split('\n\n')
patterns=patterns[1:]

gt={}
qc=tc=1 #results not start with '-'
for pattern in patterns:
    lines=pattern.strip().split('\n')
    first=lines[0].split('|')[1]
    q_offset=np.array(lines[0].split(' '))
    q_offset=q_offset[q_offset!=''][2]
    q_offset=q_offset.split('|')[0]
    q_offset=len(q_offset)
    #qc=int(lines[0].split('  ')[1])+offset
    #for results with only gt aligned regions
    qc+=q_offset

    second=lines[2].split('|')[1]
    t_offset=np.array(lines[2].split(' '))
    t_offset=t_offset[t_offset!=''][2]
    t_offset=t_offset.split('|')[0]
    t_offset=len(t_offset)
    #tc=int(lines[2].split('  ')[1])+offset
    #for results with only gt aligned regions
    tc+=t_offset

    for q,t in zip(first,second):
        if q!='-':
            gt[str(qc)+q]=str(tc)+t
            qc+=1
        if t!='-':
            tc+=1

    q_offset=lines[0].split('|')[2]
    q_offset=q_offset.split(' ')[0]
    q_offset=len(q_offset)
    t_offset=lines[2].split('|')[2]
    t_offset=t_offset.split(' ')[0]
    t_offset=len(t_offset)
    qc+=q_offset
    tc+=t_offset

#print(gt)
#print(gt.keys())

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
