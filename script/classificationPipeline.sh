# 0 31 32 63
method=$1
batch=$2

for i in `seq 0 ${batch}`; do
	nohup python classification_${method}.py ${batch} ${i} &
	#echo ${method}_${i}
done

#bin/evaluate_ordered_lists.pl ordered_pooled_new/ combinetable.pdb70 scope_140_targets.list pooled > evaluation_results/pooled_new_pdb70
