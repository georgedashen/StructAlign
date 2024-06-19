# on work:
# the tm, dali, manual alignments in *ali files only include the ground truth region
# while the provided ground truth *aln file has start and end positions

import os, sys
import argparse
import numpy as np
import csv

example_usage = 'python accuracy.py <*.manual.ali> <*.tm.ali>'
example_usage2 = 'python accuracy.py d1a05a_d1dgsa3.manual.ali d1a05a_d1dgsa3.tm.ali --verbose 0 --folder d1a05a_d1dgsa3 --outfile d1a05a_d1dgsa3.accuracy --model tm --query 1a05_Ai --target 1dgs_Ac'
parser = argparse.ArgumentParser(description='Calculate accuracy', epilog=example_usage)
parser.add_argument('gtfile', type=str, help='should be *.manual.ali')
parser.add_argument('infile', type=str, help='*.model.ali')
parser.add_argument('--verbose', type=int, default=1)
parser.add_argument('--outfile', type=str, default='accuracy_result')
parser.add_argument('--model', type=str, default='')
parser.add_argument('--query', type=str, default='query')
parser.add_argument('--target', type=str, default='target')
parser.add_argument('--folder', type=str, default='')
args = parser.parse_args()

if not args.folder:
    args.folder = os.getcwd()
if not args.model:
    args.model = args.infile.split('.')[-2]


b=np.loadtxt(args.gtfile,dtype=str,delimiter='\n')
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
a=np.loadtxt(args.infile,dtype=str,delimiter='\n')
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

if args.verbose:
    print('')
    print(f'accuracy:{accuracy}, recall:{recall}, precision:{precision}')


if not os.path.exists(args.outfile):
    with open(args.outfile,'w') as f:
        csv.writer(f).writerow(['Folder','Query','Target','Model','Precision','Recall','Accuracy'])

with open(args.outfile,'a') as f:
    csv.writer(f).writerow([args.folder, args.query, args.target, args.model, precision, recall, accuracy])
